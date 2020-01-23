#  PSpaMM - Code Generator for Sparse Matrix Multiplication
#  source code: https://github.com/peterwauligmann/PSpaMM
#
#  NOTE: it is not an official cmake search file
#  author: Ravil Dorozhinskii
#  email: ravil.dorozhinskii@tum.de 
#
#  PSpaMM_FOUND - system has Libxsmm_executable
#
#  Additional env. variables that may be used by this module. 
#  They can change the default behaviour and
#  could to be set before calling find_package:
#
#  PSpaMM_DIR - the root directory of the BLIS installation


include(FindPackageHandleStandardArgs)

find_program(PSpaMM_PROGRAM pspamm.py
  HINTS ENV PSpaMM_DIR
  DOC "Directory where the PSpaMM python script is located"
)

find_package_handle_standard_args(PSpaMM
                                  REQUIRED_VARS PSpaMM_PROGRAM)