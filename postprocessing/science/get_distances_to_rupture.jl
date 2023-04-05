# Script to calculate precise distances from SeisSol's surface/volume cells to the rupturing part of the fault
# Output will be saved as a numpy file (.npy)
# Multithreading in julia can be enabled with the flag -t
# Example start command: mpiexecjl -np 5 julia -t 10 get_distances_to_rupture.jl data-surface.xdmf

using SeisSolXDMF
using ArgParse
using NPZ

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--slip_threshold"
            help = "slip threshold to define the rupturing part of the fault"
            arg_type = Float32
            default = 0.1f0
        "--use_gpu"
            help = "enable GPU acceleration with CUDA.jl"
            action = :store_true
        "--timestep_range"
            help = "parse the timestep range that determines the rupturing part of the fault \
                    (e.g., needed if the output file contains more than one event)"
            nargs = 2
            arg_type = Int
            default = [1,-1]
        "--noMPI"
            help = "disable MPI; this flag will increase the performance of a single rank"
            action = :store_true
        "--joyner_boore_distance"
            help = "compute Joyner-Boore distance instead of rupture distance (ignoring topography)"
            action = :store_true
        "filename"
            help = "path to SeisSol's surface or volume.xdmf"
            required = true
    end
    return parse_args(s)
end


const parsed_args = parse_commandline()
const filename = parsed_args["filename"]
const threshold = parsed_args["slip_threshold"]
const timestep_range = parsed_args["timestep_range"]
const use_gpu = parsed_args["use_gpu"]
const noMPI = parsed_args["noMPI"]
const joyner_boore_distance = parsed_args["joyner_boore_distance"]


function get_cell_barycenters(xdmf)
    topology, geometry = grid_of(xdmf)
    xyz = (1. / 3.) .* (geometry[:,topology[1,:]] .+ geometry[:,topology[2,:]] .+ geometry[:,topology[3,:]])
    return copy(transpose(convert(Matrix{Float32}, xyz)))
end


function get_ASl(xdmf, timestep_range::Vector{Int})
    
    ntimesteps = size(timesteps_of(xdmf))[1]
    timestep_range[2] = (timestep_range[2] == -1 ? ntimesteps : timestep_range[2])
    
    if !(timestep_range[1] > 0 && timestep_range[1] < timestep_range[2] && 
         timestep_range[2] <= ntimesteps)
        error("Parsed timestep range is not compatible with the fault output.")
    end
    
    if timestep_range[1] == 1
        ASl = data_of(xdmf, timestep_range[2], "ASl")
    else
        ASl = data_of(xdmf, timestep_range[2], "ASl") .- data_of(xdmf, timestep_range[1], "ASl")
    end
    return convert(Vector{Float32}, ASl)
end


function get_distances_to_rupture(ASl::Vector{Float32}, xyz_fault::Matrix{Float32}, xyz_receiver::Matrix{Float32}, 
        slip_threshold::Float32, JB_distance::Bool=false)
    
    distances = Vector{Float32}(undef, size(xyz_receiver)[1])
    xyz_rupturing_fault = xyz_fault[ASl .> slip_threshold,:]
    coord_index = (JB_distance ? 2 : 3)
    
    if rank::Int == 0
        println("Using $(Threads.nthreads()) threads per MPI rank...")
    end
    Threads.@threads for i in 1:size(xyz_receiver)[1]
        distances[i] = minimum(sqrt.(sum(((xyz_rupturing_fault[:,1:coord_index] .- reshape(xyz_receiver[i,1:coord_index], 
                                1,coord_index)).^2), dims=2)))
    end
    return distances
end


if use_gpu::Bool
    using CUDA
    function get_distances_to_rupture(ASl::CuArray{Float32, 1, CUDA.Mem.DeviceBuffer}, 
                        xyz_fault::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, 
                        xyz_receiver::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, 
                        slip_threshold::Float32, JB_distance::Bool=false)
        # gpu version of the function (employing multiple dispatch)

        xyz_rupturing_fault = xyz_fault[ASl .> slip_threshold,:]
        distances = Vector{Float32}(undef, size(xyz_receiver)[1])
        coord_index = (JB_distance ? 2 : 3)
        
        for i in 1:size(xyz_receiver)[1]
            distances[i] = minimum(sqrt.(sum(((xyz_rupturing_fault[:,1:coord_index] .- reshape(xyz_receiver[i,1:coord_index], 
                                    1,coord_index)).^2), dims=2)))
        end
        return distances
    end
end


function get_chunk_sizes_and_offsets(N::Integer, n::Integer)
    q,r = divrem(N, n)
    lengths = [i <= r ? q+1 : q for i = 1:n]
    offsets = [i == 1 ? 0 : sum(lengths[1:i-1]) for i = 1:n]
    return lengths, offsets
end


if noMPI::Bool
    const rank = 0
    const nranks = 1
else
    using MPI

    MPI.Init()
    const comm = MPI.COMM_WORLD
    const rank = MPI.Comm_rank(comm)
    const nranks = MPI.Comm_size(comm)
end


function main()

    timeOutput = @timed begin

        xdmf = XDMFFile(filename::String)

        if split(filename::String, '-')[end] == "surface.xdmf" || split(filename::String, '-')[end] == "volume.xdmf"
            filename2 = rsplit(filename::String, "-"; limit = 2)[1]*"-fault.xdmf"
        elseif split(filename::String, '.')[end] == "xdmf"
            filename2 = rsplit(filename::String, "."; limit = 2)[1]*"-fault.xdmf"
        else
            error("The specified filename needs to be a SeisSol .xdmf file")
        end
        
        if !isfile(filename2)
            error("Filename deviates from standard SeisSol output, can not find ...-fault.xdmf")
        end

        xdmf_fault = XDMFFile(filename2::String)

        if use_gpu::Bool
            if !noMPI::Bool
                if length(devices()) < nranks::Int
                    error("Please only use one MPI rank per GPU.")
                end
                comm_l = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank::Int)
                rank_l = MPI.Comm_rank(comm_l)
                gpu_id = CUDA.device!(rank_l)
            end
            xyz_fault = CuArray(get_cell_barycenters(xdmf_fault))
            xyz_receiver = CuArray(get_cell_barycenters(xdmf))
            ASl = CuArray(get_ASl(xdmf_fault, timestep_range::Vector{Int}))
        else
            xyz_fault = get_cell_barycenters(xdmf_fault)
            xyz_receiver = get_cell_barycenters(xdmf)
            ASl = get_ASl(xdmf_fault, timestep_range::Vector{Int})
        end
        
        if rank::Int == 0            
            distance_type = (joyner_boore_distance ? "JB_" : "to_rupture_")
            outputname = "distance_"*distance_type*split(filename::String, '/')[end][1:end-5]
            println("Writing receiver xyz values to $(outputname)_xyz.npy")
            npzwrite(outputname*"_xyz.npy", xyz_receiver)
        end

        nElements = size(xyz_receiver)[1]
        mpi_lengths, mpi_offsets = get_chunk_sizes_and_offsets(nElements, nranks::Int)
        xyz_receiver = xyz_receiver[1+mpi_offsets[(rank::Int)+1]:mpi_lengths[(rank::Int)+1]+mpi_offsets[(rank::Int)+1],:]
    end

    if rank::Int == 0
        println("Time to load and prepare data on rank $(rank): $(timeOutput[2])s")
        println("Number of receivers: $(nElements)")
        println("Fault array size: $(size(xyz_fault))")
        println("Receiver chunk size: $(size(xyz_receiver))")
        println("Using $(nranks) MPI ranks to compute distances...")
    end

    @time begin
        distances = get_distances_to_rupture(ASl, xyz_fault, xyz_receiver, threshold::Float32, joyner_boore_distance::Bool)
        println("Time to compute distances on rank $(rank):")
    end

    if noMPI::Bool
        output = distances
    else
        MPI.Barrier(comm)

        output = (rank::Int == 0 ? Vector{Float32}(undef, nElements) : nothing)

        MPI.Gatherv!(distances, VBuffer(output, mpi_lengths, mpi_offsets, MPI.FLOAT), comm; root=0)
        MPI.Finalize()
    end

    if rank::Int == 0
        println("Writing output to $(outputname).npy")
        npzwrite(outputname*".npy", output)
        println("Done.")
    end
end

main()
