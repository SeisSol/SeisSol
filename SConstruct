# operation system (required for execution environment)
import os
import sys

# import helpers
import arch
import libs
import utils.gitversion

# print the welcome message
print('********************************************')
print('** Welcome to the build script of SeisSol **')
print('********************************************')
print('Copyright (c) 2012-2018, SeisSol Group')

helpMode = GetOption('help')

def ConfigurationError(msg):
  """Print the error message and exit. Continue only
  if the user wants to show the help message"""

  if not helpMode:
    print(msg)
    Exit(1)

def DeprecatedWarning(option, msg):
  if not helpMode:
    option_str = '** Deprecated option: {}'.format(option)
    msg_str = '** {}'.format(msg)
    star_len = max(len(option_str), len(msg_str))
    stars = '*' * star_len
    print(stars)
    print(option_str)
    print(msg_str)
    print(stars)

#
# set possible variables
#
vars = Variables()

# read parameters from a file if given
vars.AddVariables(
  PathVariable( 'buildVariablesFile', 'location of the python file, which contains the build variables', None, PathVariable.PathIsFile )
)
env = Environment(variables=vars)
if 'buildVariablesFile' in env:
  vars = Variables(env['buildVariablesFile'])

# SeisSol specific variables
vars.AddVariables(
  EnumVariable( 'equations',
                'system of PDEs that will be solved',
                'elastic',
                allowed_values=('elastic', 'viscoelastic', 'viscoelastic2', 'anisotropic')
              ),

  EnumVariable( 'order',
                'convergence order of the ADER-DG method',
                'none',
                allowed_values=('none', '2', '3', '4', '5', '6', '7', '8')
              ),

  ( 'numberOfMechanisms', 'Number of anelastic mechanisms (needs to be set if equations=viscoelastic).', '0' ),
  
  ( 'multipleSimulations', 'Fuse multiple simulations in one run.', '1' ),

  PathVariable( 'memLayout', 'Path to memory layout file.', None, PathVariable.PathIsFile),

  ( 'programName', 'name of the executable', 'none' ),

  PathVariable( 'buildDir', 'where to build the code', 'build', PathVariable.PathIsDirCreate ),

  EnumVariable( 'compileMode', 'mode of the compilation', 'release',
                allowed_values=('debug', 'release', 'relWithDebInfo')
              ),

  EnumVariable( 'parallelization', 'level of parallelization', 'none',
                allowed_values=('none', 'omp', 'mpi', 'hybrid')
              ),

  BoolVariable( 'vecReport', 'print verbose vectorization report when using Intel Compiler suite', False ),

  BoolVariable( 'hdf5', 'use hdf5 library for data output', False ),

  EnumVariable( 'netcdf', 'use netcdf library for mesh input',
                'no',
                allowed_values=('yes', 'no', 'passive')
              ),

  BoolVariable( 'metis', 'use Metis for partitioning', False ),

  BoolVariable( 'sionlib', 'use sion library for checkpointing', False ),

  BoolVariable( 'asagi', 'use asagi for material input', False ),

  BoolVariable( 'memkind', 'use memkind library for hbw memory support', False ),

  EnumVariable( 'unitTests', 'builds additional unit tests',
                'none',
                allowed_values=('none', 'fast', 'all')
              ),

  EnumVariable( 'logLevel',
                'logging level. \'debug\' runs assertations and prints all information available, \'info\' prints information at runtime (time step, plot number), \'warning\' prints warnings during runtime, \'error\' is most basic and prints errors only',
                'info',
                allowed_values=('debug', 'info', 'warning', 'error')
              ),

  EnumVariable( 'logLevel0',
                'logging level for rank 0. Default is same as logLevel',
                'none',
                allowed_values=('none', 'debug', 'info', 'warning', 'error')
              ),

  BoolVariable( 'commThread', 'use communication thread for MPI progression (option has no effect when not compiling hybrid target)', False ),

  BoolVariable( 'plasticity', 'enable plasticity (generated kernels only)', False ),

  BoolVariable( 'integrateQuants', 'enable computation and storage of integrated quantities (generated kernels only)', False ),

  EnumVariable( 'dynamicRuptureMethod',
                'Use quadrature here, cellaverage is EXPERIMENTAL.',
                'quadrature',
                allowed_values=('quadrature', 'cellaverage')
              ),

  EnumVariable( 'PlasticityMethod',
                'choose between two plasticity methods, nodal one in general faster',
                'nb',
                allowed_values=('ip', 'nb')
              ),
  
  EnumVariable( 'GemmTools',
                'choose a gemm tool(s) for the code generator',
                'LIBXSMM,PSpaMM',
                 allowed_values=('LIBXSMM,PSpaMM', 'LIBXSMM', 'MKL', 'OpenBLAS', 'BLIS')
              )
)

# external variables
vars.AddVariables(
  PathVariable( 'memkindDir',
                'memkind installation directory',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'netcdfDir',
                'NetCDF installation directory',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'hdf5Dir',
                'HDF5 installation directory',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'zlibDir',
                'zlib installation directory',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'metisDir',
                '(Par)Metis installation directory',
                None,
                PathVariable.PathAccept ),

  PathVariable( 'sionlibDir',
                'sionlib installation directory',
                None,
                PathVariable.PathAccept ),

  EnumVariable( 'compiler',
                'Select the compiler (default: intel)',
                'intel',
                allowed_values=('intel', 'gcc', 'cray_intel', 'cray_gcc')),

  BoolVariable( 'useExecutionEnvironment',
                'set variables set in the execution environment',
                True ),

  EnumVariable( 'arch',
                'precision -- s for single- and d for double precision -- and architecture used. Warning: \'noarch\' calls the fall-back code and is outperformed by architecture-specific optimizations (if available) greatly.',
                'dnoarch',
                allowed_values=arch.getArchitectures()
              ),

  EnumVariable( 'scalasca', 'instruments code with scalasca. \n \'default\': instruments only outer loops. \n'+\
                                                              ' \'kernels\': additionally instruments inner kernels.\n'+\
                                                              ' \'default_2.x\': outer loops with Scalasca version 2.x\n'+\
                                                              ' \'kernels_2.x\': loops and kernels with Scalasca version 2.x\n',
                'none',
                allowed_values=('none', 'default', 'kernels', 'default_2.x', 'kernels_2.x')
              ),
)

env.Tool('MPITool', vars=vars)

# set environment
env = Environment(variables=vars)

if env['useExecutionEnvironment']:
  env['ENV'] = os.environ

# generate help text
Help(vars.GenerateHelpText(env))

# handle unknown, maybe misspelled variables
unknownVariables = vars.UnknownVariables()

# remove the buildVariablesFile from the list of unknown variables (used before)
if 'buildVariablesFile' in unknownVariables:
  unknownVariables.pop('buildVariablesFile')

# exit in the case of unknown variables
if unknownVariables:
  ConfigurationError("*** The following build variables are unknown: " + str(unknownVariables.keys()))

if env['order'] == 'none':
  ConfigurationError("*** Convergence order not set.")

if env['equations'].startswith('viscoelastic'):
  if env['numberOfMechanisms'] == '0':
    ConfigurationError("*** Number of mechanisms not set.")
elif env['numberOfMechanisms'] != '0':
  ConfigurationError("*** Number of mechanisms must be 0 for elastic equations.")

if env['equations'] == 'viscoelastic':
  DeprecatedWarning('viscoelastic', 'Please use viscoelastic2.')


if int(env['multipleSimulations']) != 1 and int(env['multipleSimulations']) % arch.getAlignedReals(env['arch']) != 0:
  ConfigurationError("*** multipleSimulations must be a multiple of {}.".format(arch.getAlignedReals(env['arch'])))

# check for architecture
if env['arch'] == 'snoarch' or env['arch'] == 'dnoarch':
  print("*** Warning: Using fallback code for unknown architecture. Performance will suffer greatly if used by mistake and an architecture-specific implementation is available.")

if not 'memLayout' in env:
  env['memLayout'] = 'auto'

# Detect SeisSol version
seissol_version = utils.gitversion.get(env)
print('Compiling SeisSol version: {}'.format(seissol_version))

#
# preprocessor, compiler and linker
#

numberOfQuantities = {
  'elastic' : 9,
  'anisotropic' : 9,
  'viscoelastic' : 9 + int(env['numberOfMechanisms']) * 6
}
numberOfQuantities['viscoelastic2'] = numberOfQuantities['viscoelastic']

# Basic compiler setting
if env['compiler'] == 'intel':
    env['CC'] = 'icc'
    env['CXX'] = 'icpc'
    env['F90'] = 'ifort'
elif env['compiler'] == 'gcc':
    env['CC'] = 'gcc'
    env['CXX'] = 'g++'
    env['F90'] = 'gfortran'
elif env['compiler'].startswith('cray_'):
    env['CC'] = 'cc'
    env['CXX'] = 'CC'
    env['F90'] = 'ftn'
else:
    assert(false)

# Parallel compiler required?
if env['parallelization'] in ['mpi', 'hybrid']:
    env.Tool('MPITool')

    # Do not include C++ MPI Bindings
    env.Append(CPPDEFINES=['OMPI_SKIP_MPICXX'])

# Use dynamic linking on Cray and remove any special compiler configuration
# Do this after the MPI tool is called, because the MPI Tool checks for special compilers
if env['compiler'].startswith('cray'):
    env.Append(LINKFLAGS=['-dynamic'])
    env['compiler'] = env['compiler'].replace('cray_', '')

# Include preprocessor in all Fortran builds
env['F90COM'] = env['F90PPCOM']

# Use Fortran for linking
env['LINK'] = env['F90']

# Linker-flags for Fortran linking
if env['compiler'] == 'intel':
    env.Append(LINKFLAGS=['-nofor-main', '-cxxlib']) #Add -ldmalloc for ddt
elif env['compiler'] == 'gcc':
    env.Append(LIBS=['stdc++'])

#
# Scalasca
#

# instrument code with scalasca
if env['scalasca'] in ['default', 'kernels']:
  for mode in ['CC', 'CXX', 'F90', 'LINK']:
    l_scalascaPrelude = 'scalasca -instrument -comp=none -user '

    if env['parallelization'] in ['mpi', 'hybrid']:
      l_scalascaPrelude = l_scalascaPrelude + '-mode=MPI '

    env[mode] = l_scalascaPrelude + env[mode]

if env['scalasca'] in ['default_2.x', 'kernels_2.x']:
  l_scorepArguments = " --noonline-access --nocompiler --user "
  if env['parallelization'] == 'none':
    l_scorepArguments = l_scorepArguments + ' --mpp=none '
  if env['parallelization'] in ['mpi', 'hybrid']:
    l_scorepArguments = l_scorepArguments + ' --mpp=mpi '
    # The following line is required for tests
    env['CONF_PREFIX'] = """
#include <mpi.h>
void init() { MPI_Init(0, 0L); }
"""

  if env['parallelization'] in ['mpi', 'none']:
    l_scorepCxxArguments = l_scorepArguments + ' --thread=none '
  else:
    if env['commThread']:
      # Seems to work with "RF_output_on = 0"
      l_scorepCxxArguments = l_scorepArguments + ' --thread=pthread '
    else:
      if libs.find(env, 'openmp', required=False, version='3.0'):
        l_scorepCxxArguments = l_scorepArguments + ' --thread=omp:ancestry '
      else:
        l_scorepCxxArguments = l_scorepArguments + ' --thread=omp '

  for mode in ['F90']:
    env[mode] = 'scorep' + l_scorepArguments + ' --thread=none ' + env[mode]
  for mode in ['CC', 'CXX', 'LINK']:
    env[mode] = 'scorep' + l_scorepCxxArguments + env[mode]

# kernel instrumentation with scalasca
if env['scalasca'] == 'kernels_2.x':
  env.Append(CPPDEFINES=['INSTRUMENT_KERNELS'])

#
# Common settings
#

# enforce restrictive C/C++-Code
env.Append(CFLAGS   = ['-Wall', '-Werror', '-ansi'],
           CXXFLAGS = ['-Wall', '-Werror', '-ansi'])
if env['compiler'] == 'intel':
    env.Append(CXXFLGAS = ['-wd13379'])
elif env['compiler'] == 'gcc':
    # TODO Fix kernel generation
    env.Append(CXXFLAGS = ['-Wno-error=unknown-pragmas'])

# generate vector report (only if requested)
if env['vecReport']:
  env.Append(CXXFLAGS = ['-vec-report3'])
  env.Append(CFLAGS   = ['-vec-report3'])

# run preprocessor before compiling
env.Append(F90FLAGS=['-cpp'])

# Use complete line
if env['compiler'] == 'gcc':
    env.Append(F90FLAGS=['-ffree-line-length-none'])

# enforce 8 byte precision for reals (required in SeisSol) and compile time boundary check
if env['compiler'] == 'intel':
    env.Append(F90FLAGS=['-r8', '-WB'])
elif env['compiler'] == 'gcc':
    env.Append(F90FLAGS=['-fdefault-real-8'])

# Align structs and arrays
if env['compiler'] == 'intel':
    # TODO Check if Fortran alignment is still necessary in the latest version
    env.Append(F90LFAGS=['-align', '-align', 'array64byte'])

#
# Architecture dependent settings
#
archFlags = arch.getFlags(env['arch'], env['compiler'])
env.Append( CFLAGS    = archFlags,
            CXXFLAGS  = archFlags,
            F90FLAGS  = archFlags,
            LINKFLAGS = archFlags )
env.Append(CPPDEFINES=arch.getDefines(env['arch']))

#
# Compile mode settings
#

# set (pre-)compiler flags for the compile modes
if env['compileMode'] == 'debug':
  env.Append(F90FLAGS = ['-O0'],
             CLFGAGS  = ['-O0'],
             CXXFLAGS = ['-O0'])
  if env['compiler'] == 'intel':
      env.Append(F90FLAGS = ['-shared-intel', '-check'],
                 CLFGAGS  = ['-shared-intel'],
                 CXXFLAGS = ['-shared-intel'])
  else:
      env.Append(F90FLAGS = ['-fcheck=all'])
if env['compileMode'] in ['debug', 'relWithDebInfo']:
  env.Append(F90FLAGS  = ['-g'],
             CFLAGS    = ['-g'],
             CXXFLAGS  = ['-g'],
             LINKFLAGS = ['-g', '-rdynamic'])
  if env['compiler'] == 'intel':
      env.Append(F90FLAGS  = ['-traceback'],
                 CFLAGS    = ['-traceback'],
                 CXXFLAGS  = ['-traceback'])
  else:
      env.Append(F90FLAGS  = ['-fbacktrace'])

if env['compileMode'] in ['relWithDebInfo', 'release']:
    env.Append(CPPDEFINES = ['NDEBUG'])
    env.Append(F90FLAGS = ['-O2'],
               CFLAGS   = ['-O2'],
               CXXFLAGS = ['-O2'])
    if env['compiler'] == 'intel':
        env.Append(F90FLAGS = ['-fno-alias'])

# C++ Standard
env.Append(CXXFLAGS=['-std=c++11'])

#
# Basic preprocessor defines
#
# set precompiler mode for the number of quantities and basis functions
env.Append(CPPDEFINES=['CONVERGENCE_ORDER='+env['order']])
env.Append(CPPDEFINES=['NUMBER_OF_QUANTITIES=' + str(numberOfQuantities[ env['equations'] ]), 'NUMBER_OF_RELAXATION_MECHANISMS=' + str(env['numberOfMechanisms'])])

env.Append(CPPDEFINES=['USE_{}'.format(env['equations'].upper())])

env.Append(CPPDEFINES=['ENABLE_MATRIX_PREFETCH'])

if int(env['multipleSimulations']) > 1:
  env.Append(CPPDEFINES=['MULTIPLE_SIMULATIONS={}'.format(env['multipleSimulations'])])

# add parallel flag for mpi
if env['parallelization'] in ['mpi', 'hybrid']:
    # TODO rename PARALLEL to USE_MPI in the code
    env.Append(CPPDEFINES=['PARALLEL', 'USE_MPI'])

# add OpenMP flags
if env['parallelization'] in ['omp', 'hybrid']:
    env.Append(CPPDEFINES=['OMP'])
    env.Append(CFLAGS    = ['-fopenmp'],
               CXXFLAGS  = ['-fopenmp'],
               F90FLAGS  = ['-fopenmp'],
               LINKFLAGS = ['-fopenmp'])

if( env['plasticity'] ):
  env.Append(CPPDEFINES=['USE_PLASTICITY'])
  if env['PlasticityMethod'] == 'ip':
     env.Append(CPPDEFINES=['USE_PLASTIC_IP'])
  elif env['PlasticityMethod'] == 'nb':
     env.Append(CPPDEFINES=['USE_PLASTIC_NB'])

if( env['integrateQuants'] ):
  env.Append(CPPDEFINES=['INTEGRATE_QUANTITIES'])

# set pre compiler flags for matrix optimizations
env.Append(CPPDEFINES=['GENERATEDKERNELS', 'CLUSTERED_LTS'])

# set pre compiler flags commuincation thread
# pthread is linked after the other libraries
if env['commThread']:
  env.Append(CPPDEFINES=['USE_COMM_THREAD'])

if env['dynamicRuptureMethod'] == 'cellaverage':
  env.Append(CPPDEFINES=['USE_DR_CELLAVERAGE'])

# Default log level for rank 0 is same as logLevel
if env['logLevel0'] == 'none':
  env['logLevel0'] = env['logLevel']

# set level of logger for fortran
if env['logLevel'] == 'debug':
  env.Append(CPPDEFINES=['LOGLEVEL=3'])
elif env['logLevel'] == 'info':
  env.Append(CPPDEFINES=['LOGLEVEL=2'])
elif env['logLevel'] == 'warning':
  env.Append(CPPDEFINES=['LOGLEVEL=1'])
elif env['logLevel'] == 'error':
  env.Append(CPPDEFINES=['LOGLEVEL=0'])
else:
  assert(false)

# set level of logger for rank 0 and C++
if env['logLevel0'] == 'debug':
  env.Append(CPPDEFINES=['LOGLEVEL0=3', 'LOG_LEVEL=3'])
elif env['logLevel0'] == 'info':
  env.Append(CPPDEFINES=['LOGLEVEL0=2', 'LOG_LEVEL=2'])
elif env['logLevel0'] == 'warning':
  env.Append(CPPDEFINES=['LOGLEVEL0=1', 'LOG_LEVEL=1'])
elif env['logLevel0'] == 'error':
  env.Append(CPPDEFINES=['LOGLEVEL0=0', 'LOG_LEVEL=0'])
else:
  assert(false)

# add include path for submodules
env.Append( CPPPATH=['#/submodules', '#/submodules/yateto/include', '#/submodules/eigen3'] )
#
# add libraries
#

env.Tool('cmake')

# Libxsmm
env.Tool('LibxsmmTool', required=True)

# Library pathes
env.Tool('DirTool', fortran=True)

# netCDF
if env['netcdf'] == 'yes':
    libs.find(env, 'netcdf', required=(not helpMode), parallel=(env['parallelization'] in ['hybrid', 'mpi']))
    env.Append(CPPDEFINES=['USE_NETCDF'])
elif env['netcdf'] == 'passive':
    env.Append(CPPDEFINES=['USE_NETCDF', 'NETCDF_PASSIVE'])

# HDF5
if env['hdf5']:
    libs.find(env, 'hdf5', required=(not helpMode), parallel=(env['parallelization'] in ['hybrid', 'mpi']))

    env.Append(CPPDEFINES=['USE_HDF'])

# memkind
if env['memkind']:
  env.Tool('MemkindTool')
  env.Append(CPPDEFINES=['USE_MEMKIND'])

# Metis
if env['metis'] and env['parallelization'] in ['hybrid', 'mpi']:
	libs.find(env, 'metis', required=(not helpMode), parallel=True)
	env.Append(CPPDEFINES=['USE_METIS'])

# sionlib still need to create a Tool for autoconfiguration
if env['sionlib']:
  env.Tool('SionTool', parallel=(env['parallelization'] in ['hybrid', 'mpi']))
else:
  env['sionlib'] = False

# ASAGI
if env['asagi']:
    if env['scalasca'] in ['default', 'kernels']:
        ConfigurationError("*** ASAGI can not run with Scalasca 1.x")

    libs.find(env, 'asagi', parallel=(env['parallelization'] in ['hybrid', 'mpi']), required=(not helpMode))
    env.Append(CPPDEFINES=['USE_ASAGI'])

# yaml-cpp
yaml_cpp = env.CMake( source=[Glob(path + '*/.cpp') for path, dirs, files in os.walk('submodules/yaml-cpp/src')],
                      target=['#/{}/external/yaml-cpp/libyaml-cpp.a'.format(env['buildDir'])],
                      CMakeProject = Dir('submodules/yaml-cpp'),
                      CMakeOpts = ['-DYAML_CPP_BUILD_TOOLS=no', '-DCMAKE_CXX_STANDARD=11', '-DYAML_CPP_BUILD_TESTS=OFF'],
                      cc = env['CC'],
                      cxx = env['CXX'])
env.Append(CPPPATH=['#/submodules/yaml-cpp/include'])
env.Append(LIBS=yaml_cpp)

# impalajit
impalajit = env.CMake( source=[Glob(path + '*/.cc') for path, dirs, files in os.walk('submodules/ImpalaJIT')],
                       target=['#/{}/external/impalajit/libimpalajit.a'.format(env['buildDir'])],
                       CMakeProject = Dir('submodules/ImpalaJIT'),
                       CMakeOpts = ['-DCMAKE_CXX_STANDARD=11'],
                       CMakeBuildDir = Dir('#/{}/external/impalajit/'.format(env['buildDir'])),
                       cc = env['CC'],
                       cxx = env['CXX'])
env.Append(CPPPATH=['#/submodules/ImpalaJIT/include'])
env.Append(LIBS=impalajit)

# easi
env.Append(CPPPATH=['#/submodules/easi/include'])

# ASYNC I/O
env.Append(CPPPATH=['#/submodules/async'])

# Required for communication thread and asynchronous I/O
# pthread has to be appended after other libraries
# this only appears when compiling with scalasca and hdf5/netcdf
env.Append(LIBS=['pthread'])

# add pathname to the list of directories wich are search for include
env.Append(F90FLAGS=['-Isrc'])
env.Append(CPPPATH=['#/src', '#/src/Equations/' + env['equations']])

#
# setup the program name and the build directory
#
if env['programName'] == 'none':
  program_suffix = '%s_%s_%s_%s_%s_%s_%s' %(
    env['compileMode'],
    'generatedKernels',
    env['arch'],
    env['parallelization'],
    'scalasca' if env['scalasca'] != 'none' else 'none',
    numberOfQuantities[ env['equations'] ],
    env['order']
  )
  if int(env['multipleSimulations']) > 1:
    program_suffix = '{}_{}sims'.format(program_suffix, env['multipleSimulations'])
  env['programFile'] = '%s/SeisSol_%s' %(
    env['buildDir'],
    program_suffix
  )
else:
  program_suffix = env['programName']
  env['programFile'] = env['buildDir']+'/'+env['programName']

# build directory

env['buildDir'] = '%s/build_%s' %(env['buildDir'], program_suffix)

# set module path
if env['compiler'] == 'intel':
    env.Append(F90FLAGS='-module ${TARGET.dir}')
elif env['compiler'] == 'gcc':
    env.Append(F90FLAGS='-J ${TARGET.dir}')

# get the source files
env.sourceFiles = []
env.generatedSourceFiles = []
env.generatedTestSourceFiles = []

# Generate the version file
utils.gitversion.generateHeader(env, target=env['buildDir']+'/version.h')

Export('env')
SConscript('generated_code/SConscript', variant_dir=env['buildDir'] + '/generated_code', duplicate=0)
SConscript('src/SConscript', variant_dir=env['buildDir'] + '/src', duplicate=0)
SConscript('submodules/SConscript', variant_dir=env['buildDir']+'/submodules', duplicate=0)
Import('env')

# add generated source files
env.sourceFiles.extend(env.generatedSourceFiles)

# remove .mod entries for the linker
modDirectories = []
sourceFiles = []
for sourceFile in env.sourceFiles:
  sourceFiles.append(sourceFile[0])
  if len(sourceFile) > 1:
    modDir = os.path.dirname(str(sourceFile[1]))
    modDirectories.append(modDir)
for directory in set(modDirectories):
  Execute(Mkdir(directory))
env.AppendUnique(F90PATH=list(map(lambda x: '#/' + x, modDirectories)))

#print(env.Dump())

# build standard version
seissol_build = env.Program('#/'+env['programFile'], sourceFiles)
env.Depends(seissol_build, yaml_cpp)
env.Depends(seissol_build, impalajit)

# build unit tests
if env['unitTests'] != 'none':
  # Anything done here should only affect tests
  env = env.Clone()

  # Compile XDMF cube check
  if env['hdf5']:
      Export('env')
      SConscript('postprocessing/validation/XDMFCubeCheck/SConscript', duplicate=0)

  # define location of cxxtest
  env['CXXTEST'] = 'submodules/cxxtest'

  # Continue testing if tests fail
  env['CXXTEST_SKIP_ERRORS'] = True

  # Replace main with MPI main
  env['CXXTEST_OPTS'] = '--template=' + Dir('.').srcnode().abspath + '/src/tests/mpirunner.tpl'

  # Fail on error (as we can't see OK messages in the output)
  env.Append(CPPDEFINES=['CXXTEST_HAVE_EH', 'CXXTEST_ABORT_TEST_ON_FAIL'])
  env.Append(CPPDEFINES={'SEISSOL_TESTS': '"\\"' + Dir('.').srcnode().abspath + '/src/tests/\\""'})

  # Try to remove weird linker errors
  if env['compiler'] == 'intel':
    env.Append(CXXFLAGS = ['-ffreestanding'])

  # add cxxtest-tool
  env.Tool('cxxtest')

  # Get test source files
  env.sourceFiles = []
  env.testSourceFiles = []
  env.mpiTestSourceFiles = dict()
  for i in range(2, 5):
    env.mpiTestSourceFiles[i] = []

  Export('env')
  SConscript('src/tests/SConscript', variant_dir='#/'+env['buildDir']+'/tests', src_dir='#/')
  Import('env')

  # Remove main() to avoid double definition
  sourceFiles = list(filter(lambda sf: os.path.basename(str(sf)) != 'main.o', sourceFiles))
  # Remove .mod files from additional Fortran files
  for sourceFile in env.sourceFiles:
    sourceFiles.append(sourceFile[0])

  if env.generatedTestSourceFiles:
    if env['parallelization'] in ['mpi', 'hybrid']:
      env['CXXTEST_COMMAND'] = 'mpirun -np 1 %t'
    mpiObject = list(filter(lambda sf: os.path.basename(str(sf)) == 'MPI.o', sourceFiles))
    env.CxxTest(target='#/'+env['buildDir']+'/tests/generated_kernels_test_suite', source=env.generatedSourceFiles+env.generatedTestSourceFiles+mpiObject)

  if env.testSourceFiles:
    if env['parallelization'] in ['mpi', 'hybrid']:
      env['CXXTEST_COMMAND'] = 'mpirun -np 1 %t'
    env.CxxTest(target='#/'+env['buildDir']+'/tests/serial_test_suite', source=sourceFiles+env.testSourceFiles)

  if env['parallelization'] in ['mpi', 'hybrid']:
    for ranks, mpiTestSourceFiles in env.mpiTestSourceFiles.items():
      if mpiTestSourceFiles:
        env['CXXTEST_COMMAND'] = 'mpirun -np {0} %t'.format(ranks)
        env.CxxTest(target='#/'+env['buildDir']+'/tests/parallel_test_suite_{0}'.format(ranks), source=sourceFiles+mpiTestSourceFiles)
