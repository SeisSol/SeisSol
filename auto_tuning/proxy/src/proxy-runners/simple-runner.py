import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import argparse
import proxy_bindings as pb

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cells', default=100000, type=int, help="num cells in a time cluster")
parser.add_argument('-t', '--timesteps', default=20, type=int, help="num time steps/repeats")
parser.add_argument('-k', '--kernel', default='all', type=str, help="kernel types")
args = parser.parse_args()

config = pb.ProxyConfig()
config.cells = args.cells
config.timesteps = args.timesteps
config.kernel = pb.Aux.str_to_kernel(args.kernel)

output = pb.run_proxy(config)
pb.Aux.display_output(output, args.kernel)
