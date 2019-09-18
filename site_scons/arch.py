import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))
codegen_include_path = os.path.join(script_dir, '..', 'generated_code')
sys.path.append(codegen_include_path)
from generated_code.arch import getArchitectures, getRealSize, getCpu

def getAlignment(architecture):
  alignments = {
      'noarch': 16,
      'wsm': 16,
      'snb': 32,
      'hsw': 32,
      'knc': 64,
      'knl': 64,
      'skx': 64
  }
  return alignments[ getCpu(architecture) ]

def getAlignedReals(architecture):
  bytesPerReal = 8 if architecture[0] == 'd' else 4
  return getAlignment(architecture) // bytesPerReal
  
def getFlags(architecture, compiler):
  if architecture not in getArchitectures():
    raise ValueError('Unknown architecture.')
  
  cpu = getCpu(architecture)
  
  if cpu == 'wsm':
    flags = ['-msse3']
  elif cpu == 'snb':
    flags =  ['-mavx']
  elif cpu == 'hsw':
    if compiler == 'intel':
      flags = ['-xCORE-AVX2', '-fma']
    else:
      flags = ['-mavx2', '-mfma']
  elif cpu == 'knc':
    flags = ['-mmic', '-fma']
  elif cpu == 'knl':
    if compiler == 'intel':
      flags = ['-xMIC-AVX512', '-fma']
    else:
      flags = ['-mavx512f', '-mavx512cd', '-mavx512pf', '-mavx512er', '-mfma']
  elif cpu == 'skx':
    if compiler == 'intel':
      flags = ['-xCORE-AVX512', '-fma']
    else:
      flags = ['-march=skylake-avx512']
  else:
    flags = []
  
  # enable interproc. opts for small cores
  if cpu in ['knc', 'knl', 'skx']:
    flags.extend(['-ip'])
              
  return flags

def getDefines(architecture):
  alignment = 'ALIGNMENT={}'.format(getAlignment(architecture))
  precision = 'REAL_SIZE={}'.format(getRealSize(architecture))
  return [alignment, precision]
