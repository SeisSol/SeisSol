add_library(SeisSol-lib

src/Initializer/ParameterDB.cpp
src/Initializer/PointMapper.cpp
src/Initializer/GlobalData.cpp
src/Initializer/InternalState.cpp
src/Initializer/MemoryAllocator.cpp
src/Initializer/CellLocalMatrices.cpp

src/Initializer/time_stepping/LtsLayout.cpp
src/Initializer/tree/Lut.cpp
src/Initializer/MemoryManager.cpp
src/Initializer/InitialFieldProjection.cpp
src/Modules/Modules.cpp
src/Modules/ModulesC.cpp
src/Model/common.cpp
src/Numerical_aux/Functions.cpp
src/Numerical_aux/Transformation.cpp
src/Numerical_aux/Statistics.cpp

src/generated_code/subroutine.h
src/generated_code/tensor.cpp
src/generated_code/subroutine.cpp
src/generated_code/tensor.h
src/generated_code/init.cpp
#src/generated_code/KernelTest.t.h
src/generated_code/init.h
src/generated_code/kernel.h
src/generated_code/kernel.cpp

src/Solver/Simulator.cpp
src/Solver/FreeSurfaceIntegrator.cpp
src/Solver/Interoperability.cpp
src/Solver/time_stepping/MiniSeisSol.cpp
src/Solver/time_stepping/TimeCluster.cpp
src/Solver/time_stepping/TimeManager.cpp
src/Solver/Pipeline/DrTuner.cpp
src/Kernels/DynamicRupture.cpp
src/Kernels/Plasticity.cpp
src/Kernels/TimeCommon.cpp
src/Kernels/Receiver.cpp
src/SeisSol.cpp
src/SourceTerm/Manager.cpp

src/SourceTerm/PointSource.cpp
src/Parallel/Pin.cpp
src/Parallel/MPI.cpp
src/Parallel/mpiC.cpp
src/Parallel/FaultMPI.cpp
src/Geometry/GambitReader.cpp

src/Geometry/MeshReaderFBinding.cpp
src/Geometry/MeshTools.cpp
src/Monitoring/FlopCounter.cpp
src/Monitoring/LoopStatistics.cpp
src/Reader/readparC.cpp
#Reader/StressReaderC.cpp
src/Checkpoint/Manager.cpp


# Checkpoint/sionlib/Wavefield.cpp
# Checkpoint/sionlib/Fault.cpp

src/Checkpoint/Backend.cpp
src/Checkpoint/Fault.cpp
src/Checkpoint/posix/Wavefield.cpp
src/Checkpoint/posix/Fault.cpp
src/ResultWriter/AnalysisWriter.cpp
src/ResultWriter/FreeSurfaceWriterExecutor.cpp
src/ResultWriter/PostProcessor.cpp
src/ResultWriter/FaultWriterC.cpp
src/ResultWriter/ReceiverWriter.cpp
src/ResultWriter/FaultWriterExecutor.cpp
src/ResultWriter/FaultWriter.cpp
src/ResultWriter/WaveFieldWriter.cpp
src/ResultWriter/FreeSurfaceWriter.cpp

# Fortran:
src/Geometry/allocate_mesh.f90
src/Geometry/MeshReaderCBinding.f90
src/Solver/close_seissol.f90
src/Solver/calc_deltat.f90
src/Solver/mpiexchangevalues.f90
src/Solver/prak_clif_mod.f90
src/Solver/calc_seissol.f90
src/Solver/f_ctof_bind_interoperability.f90
src/Solver/f_ftoc_bind_interoperability.f90
src/Numerical_aux/quadpoints.f90
src/Numerical_aux/jacobinormal.f90
src/Numerical_aux/convertxieta2xy.f90
src/Numerical_aux/create_fault_rotationmatrix.f90
src/Numerical_aux/trilinearinterpolation.f90
src/Numerical_aux/typesdef.f90
src/Numerical_aux/dgbasis.f90
src/Numerical_aux/gauss.f90
src/Numerical_aux/operators.f90
src/Numerical_aux/ODEInt.cpp
src/Numerical_aux/ODEVector.cpp
src/Modules/ModulesF.f90
src/seissolxx.f90
src/Physics/ini_model.f90
src/Physics/Evaluate_friction_law.f90
src/Physics/ini_model_DR.f90
src/Physics/NucleationFunctions.f90
src/Physics/thermalpressure.f90
src/Physics/InitialField.cpp
src/Reader/readpar.f90
src/Reader/read_backgroundstress.f90
src/ResultWriter/inioutput_seissol.f90
src/ResultWriter/magnitude_output.f90
src/ResultWriter/output_rupturefront.f90
src/ResultWriter/ini_faultoutput.f90
src/ResultWriter/energies.f90
src/ResultWriter/FaultWriterF.f90
src/ResultWriter/faultoutput.f90
src/ResultWriter/common_fault_receiver.f90
src/ResultWriter/receiver.f90
src/Initializer/dg_setup.f90
src/Initializer/ini_optionalfields.f90
src/Initializer/ini_seissol.f90
src/Parallel/mpiF.f90

src/Equations/poroelastic/Model/datastructures.cpp
src/Equations/elastic/Kernels/GravitationalFreeSurfaceBC.cpp

src/DynamicRupture/Factory.cpp
src/DynamicRupture/FrictionLaws/AgingLaw.cpp
src/DynamicRupture/FrictionLaws/BaseFrictionLaw.cpp
src/DynamicRupture/FrictionLaws/ImposedSlipRates.cpp
src/DynamicRupture/FrictionLaws/LinearSlipWeakening.cpp
src/DynamicRupture/FrictionLaws/NoFault.cpp
src/DynamicRupture/FrictionLaws/RateAndState.cpp
src/DynamicRupture/FrictionLaws/SlipLaw.cpp
src/DynamicRupture/FrictionLaws/VelocityWeakening.cpp
src/DynamicRupture/Initializers/BaseDRInitializer.cpp
src/DynamicRupture/Initializers/ImposedSlipInitializer.cpp
src/DynamicRupture/Initializers/LinearSlipWeakeningInitializer.cpp
src/DynamicRupture/Initializers/RateAndStateInitializer.cpp

src/DynamicRupture/Output/Base.cpp
src/DynamicRupture/Output/FaultRefiner/FaultRefiners.cpp
src/DynamicRupture/Output/OutputAux.cpp
)
if (MPI)
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/Wavefield.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/FaultAsync.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/Fault.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/mpio/WavefieldAsync.cpp
)
endif()

target_compile_options(SeisSol-lib PUBLIC ${EXTRA_CXX_FLAGS})

if (HDF5)
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/h5/Wavefield.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Checkpoint/h5/Fault.cpp
    )
endif()

if (HDF5 AND METIS AND MPI)
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry/PUMLReader.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/time_stepping/LtsWeights/LtsWeights.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/time_stepping/LtsWeights/WeightsModels.cpp
    )
endif()


if (NETCDF)
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/SourceTerm/NRFReader.cpp # if netCDF
    )
endif()

if (ASAGI)
  target_sources(SeisSol-lib PUBLIC
    #todo:
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Reader/AsagiModule.cpp
    )
endif()

# Eqations have to be set at compile time currently.
if ("${EQUATIONS}" STREQUAL "elastic")
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/elastic)
  target_compile_definitions(SeisSol-lib PUBLIC USE_ELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic")
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic/Kernels/Time.cpp
    )
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic)
  target_compile_definitions(SeisSol-lib PUBLIC USE_VISCOELASTIC)

elseif ("${EQUATIONS}" STREQUAL "viscoelastic2")
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/viscoelastic2)
  target_compile_definitions(SeisSol-lib PUBLIC USE_VISCOELASTIC2)

elseif ("${EQUATIONS}" STREQUAL "anisotropic")
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/DirichletBoundary.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/anisotropic)
  target_compile_definitions(SeisSol-lib PUBLIC USE_ANISOTROPIC)

elseif ("${EQUATIONS}" STREQUAL "poroelastic")
  target_sources(SeisSol-lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Neighbor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Local.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic/Kernels/Time.cpp
  )
  target_include_directories(SeisSol-lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/Equations/poroelastic)
  target_compile_definitions(SeisSol-lib PUBLIC USE_STP)
  target_compile_definitions(SeisSol-lib PUBLIC USE_POROELASTIC)
endif()

target_include_directories(SeisSol-lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders)

if ("${DEVICE_BACKEND}" STREQUAL "CUDA")

  target_sources(SeisSol-lib PUBLIC
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/LocalIntegrationRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/NeighIntegrationRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/PlasticityRecorder.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/src/Initializer/BatchRecorders/DynamicRuptureRecorder.cpp)

  find_package(CUDA REQUIRED)
  set(CUDA_NVCC_FLAGS -std=c++14;
                      -Xptxas -v;
                      -arch=${DEVICE_SUB_ARCH};
                      -DREAL_SIZE=${REAL_SIZE_IN_BYTES};
                      --compiler-options ${EXTRA_CXX_FLAGS};
                      -O3;)

  set(DEVICE_SRC ${DEVICE_SRC} ${CMAKE_BINARY_DIR}/src/generated_code/gpulike_subroutine.cpp
                               ${CMAKE_CURRENT_SOURCE_DIR}/src/Kernels/DeviceAux/cuda/PlasticityAux.cu)
  set_source_files_properties(${DEVICE_SRC} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)

  cuda_add_library(Seissol-device-lib STATIC ${DEVICE_SRC})
  target_include_directories(Seissol-device-lib PUBLIC ${DEVICE_INCLUDE_DIRS}
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/submodules/yateto/include
                                                       ${CMAKE_BINARY_DIR}/src/generated_code
                                                       ${CMAKE_CURRENT_SOURCE_DIR}/src
                                                       ${CUDA_TOOLKIT_ROOT_DIR})
  target_compile_options(Seissol-device-lib PRIVATE ${EXTRA_CXX_FLAGS})

  target_link_libraries(SeisSol-lib PUBLIC Seissol-device-lib)
  add_dependencies(Seissol-device-lib SeisSol-lib)

endif()
