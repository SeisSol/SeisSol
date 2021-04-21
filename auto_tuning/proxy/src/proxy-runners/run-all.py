import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import argparse
import proxy_bindings as pb
from progress.bar import Bar
import json


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cells', default=50000, type=int, help="num cells in a time cluster")
parser.add_argument('-t', '--timesteps', default=30, type=int, help="num time steps/repeats")
parser.add_argument('-o', '--output_type', default='c', type=str, help="console(c), json (j)")
args = parser.parse_args()

if not args.output_type in ['j', 'c']:
  print(f'error: expected `j` of `c` as output type, given `{args.output_type}`')
  parser.print_help()
  sys.exit(-1)

kernels = [pb.Kernel.all,
           pb.Kernel.local,
           pb.Kernel.neigh,
           pb.Kernel.ader,
           pb.Kernel.localwoader,
           pb.Kernel.neigh_dr,
           pb.Kernel.godunov_dr]

output = {'time': [],
          'HW GFlops': [],
          'NZ GFlops': [],
          'kernel_str': []}

config = pb.ProxyConfig()
config.cells = args.cells
config.timesteps = args.timesteps
config.verbose = False

with Bar('proxy...    ', max=len(kernels)) as bar:
  for kernel in kernels:
    config.kernel = kernel
    result = pb.run_proxy(config)

    output['time'].append(result.time)
    output['HW GFlops'].append(result.hardware_gflops)
    output['NZ GFlops'].append(result.non_zero_gflops)
    output['kernel_str'].append(pb.Aux.kernel_to_string(kernel))
    bar.next()

if args.output_type == 'j':
  json_output = json.dumps(output)
  print(json_output)
else:
  for key, value in output.items():
    print(f'{key}:\t', value)
