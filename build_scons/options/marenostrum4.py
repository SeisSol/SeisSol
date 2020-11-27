#! /usr/bin/python

# @section DESCRIPTION
# Build parameters for MareNostrum4
#

# build options
compileMode                 = 'release'
parallelization             = 'hybrid'
generatedKernels            = 'yes'
measureNodeLevelPerformance = 'none'
useExecutionEnvironment     = 'yes'
logLevel                    = 'warning'
logLevel0                   = 'info'

# machine dependent options (SNB-EP, Bulldozer)
arch = 'dskx'
