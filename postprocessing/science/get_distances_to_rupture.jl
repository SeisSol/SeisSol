# Script to calculate precise distances from SeisSol's surface/volume cells to the rupturing part of the fault
# Output will be saved as a numpy file (.npy)
# Multithreading in julia can be enabled with the flag -t 10
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
        "--event"
            help = "usually 0; 1 and 2 events are needed for Ridgecrest where one simulation contains two events"
            arg_type = Int
            default = 0
        "--noMPI"
            help = "disable MPI; this flag will increase the performance of a single rank"
            action = :store_true
        "filename"
            help = "path to SeisSol's surface or volume.xdmf"
            required = true
    end
    return parse_args(s)
end


parsed_args = parse_commandline()
filename = parsed_args["filename"]
threshold = parsed_args["slip_threshold"]
event = parsed_args["event"]
use_gpu = parsed_args["use_gpu"]
noMPI = parsed_args["noMPI"]


function get_xyz(xdmf)
    topology, geometry = grid_of(xdmf)
    xyz = (1. / 3.) .* (geometry[:,topology[1,:]] .+ geometry[:,topology[2,:]] .+ geometry[:,topology[3,:]])
    out = copy(transpose(convert(Matrix{Float32}, xyz)))
    return out
end


function get_ASl(xdmf, event)
    timesteps = timesteps_of(xdmf)
    if event == 0
        ASl = data_of(xdmf, size(timesteps)[1], "ASl")
    elseif event == 1
        ASl = data_of(xdmf, Int(round(size(timesteps)[1] / 2) - 1), "ASl")
    elseif event == 2
        ASl = data_of(xdmf, size(timesteps)[1], "ASl") .- data_of(xdmf, Int(round(size(timesteps)[1] / 2) - 1), "ASl")
    else
        println("Only 0, 1, 2 are allowed for the event argument.")
        ccall(:jl_exit, Cvoid, (Int32,), 86) # Exit code 86 to the O/S
    end
    return convert(Vector{Float32}, ASl)
end


function get_distances_to_rupture(ASl::Vector{Float32}, xyz_fault::Matrix{Float32}, xyz_receiver::Matrix{Float32}, 
        slip_threshold::Float32)
    
    distances = Vector{Float32}(undef, size(xyz_receiver)[1])
    rupturing_fault = xyz_fault[ASl .> slip_threshold,:]
    if rank == 0
        println("Using $(Threads.nthreads()) threads per MPI rank...")
    end
    Threads.@threads for i in 1:size(xyz_receiver)[1]
        distances[i] = minimum(sqrt.(sum(((rupturing_fault .- reshape(xyz_receiver[i,:], 1,3)).^2), dims=2)))
    end
    return distances
end


if use_gpu
    using CUDA
    function get_distances_to_rupture(ASl::CuArray{Float32, 1, CUDA.Mem.DeviceBuffer}, 
                        xyz_fault::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, 
                        xyz_receiver::CuArray{Float32, 2, CUDA.Mem.DeviceBuffer}, 
                        slip_threshold::Float32)
        # gpu version of the function (employing multiple dispatch)

        rupturing_fault = xyz_fault[ASl .> slip_threshold,:]
        distances = Vector{Float32}(undef, size(xyz_receiver)[1])
        for i in 1:size(xyz_receiver)[1]
            distances[i] = minimum(sqrt.(sum(((rupturing_fault .- reshape(xyz_receiver[i,:], 1,3)).^2), dims=2)))
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


if noMPI
    rank = 0
    nRanks = 1
else
    using MPI
    
    MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nRanks = MPI.Comm_size(comm)
end


timeOutput = @timed begin
    
    xdmf = XDMFFile(filename)
    
    if split(filename, '-')[end] == "surface.xdmf"
        filename2 = filename[1:end-12]*"fault.xdmf"
    elseif split(filename, '-')[end] == "volume.xdmf"
        filename2 = filename[1:end-11]*"fault.xdmf"
    elseif split(filename, '.')[end] == "xdmf"
        filename2 = filename[1:end-5]*"-fault.xdmf"
    else
        println("Filename deviates from standard SeisSol output, can not find ...-fault.xdmf")
        ccall(:jl_exit, Cvoid, (Int32,), 86) # Exit code 86 to the O/S
    end
    
    xdmf_fault = XDMFFile(filename2)
    
    if use_gpu
        if !noMPI
            if length(devices()) < nRanks
                println("Please only use one MPI rank per GPU.\nExit.")
                ccall(:jl_exit, Cvoid, (Int32,), 86)   
            end
            comm_l = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank)
            rank_l = MPI.Comm_rank(comm_l)
            gpu_id = CUDA.device!(rank_l)
        end
        xyz_fault = CuArray(get_xyz(xdmf_fault))
        xyz_receiver = CuArray(get_xyz(xdmf))
        ASl = CuArray(get_ASl(xdmf_fault, event))
    else
        xyz_fault = get_xyz(xdmf_fault)
        xyz_receiver = get_xyz(xdmf)
        ASl = get_ASl(xdmf_fault, event)
    end
    
    nElements = size(xyz_receiver)[1]
    mpi_lengths, mpi_offsets = get_chunk_sizes_and_offsets(nElements, nRanks)
    xyz_receiver = xyz_receiver[1+mpi_offsets[rank+1]:mpi_lengths[rank+1]+mpi_offsets[rank+1],:]
end

if rank == 0
    println("Time to load and prepare data on rank $(rank): $(timeOutput[2])s")
    println("Number of receivers: $(nElements)")
    println("Fault array size: $(size(xyz_fault))")
    println("Receiver chunk size: $(size(xyz_receiver))")
    println("Using $(nRanks) MPI ranks to compute distances...")
end

@time begin
    distances = get_distances_to_rupture(ASl, xyz_fault, xyz_receiver, threshold)
    println("Time to compute distances on rank $(rank):")
end

if noMPI
    output = distances
else
    MPI.Barrier(comm)
    
    output = (rank == 0 ? Vector{Float32}(undef, nElements) : nothing)

    MPI.Gatherv!(distances, VBuffer(output, mpi_lengths, mpi_offsets, MPI.FLOAT), comm; root=0)
    MPI.Finalize()
end

if rank == 0
    
    outputname = "distance_to_rupture_"*split(filename, '/')[end][1:end-5]
    
    outputname = (event == 0 ? outputname*".npy" : outputname*"_e"*string(event)*".npy")
    
    println("Writing output to $(outputname)")
    npzwrite(outputname, output)
    println("Done.")
end