import os, sys, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import argparse
import proxy_bindings as pb
from progress.bar import Bar
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import time
import os


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cells', default=50000, type=int, help="num cells in a time cluster")
parser.add_argument('-t', '--timesteps', default=30, type=int, help="num time steps/repeats")
parser.add_argument('--output_type', default='c', type=str, help="console(p), json (j), matplotlib (m)")
parser.add_argument('-o', '--output_dir', default='.', type=str, help="relative output directory")
args = parser.parse_args()

if not args.output_type in ['j', 'c', 'p']:
    print(f'notice: output only to terminal')

kernels = [pb.Kernel.all,
           pb.Kernel.local,
           pb.Kernel.neigh,
           pb.Kernel.ader,
           pb.Kernel.localwoader,
           pb.Kernel.neigh_dr,
           pb.Kernel.godunov_dr]

config = pb.ProxyConfig()
config.cells = args.cells
config.timesteps = args.timesteps
config.verbose = False

df = pd.DataFrame(columns=['kernel type', 'time', 'HW GFLOPS', 'NZ GFLOPS'])
with Bar('proxy...    ', max=len(kernels)) as bar:
    for index, kernel in enumerate(kernels):
        config.kernel = kernel
        result = pb.run_proxy(config)
        df.loc[index] = [pb.Aux.kernel_to_string(kernel),
                         result.time,
                         result.hardware_gflops,
                         result.non_zero_gflops]
        bar.next()

# prepare a unique file suffix
time_obj = time.localtime()
now = time.strftime("%B-%H-%M-%S_%Y", time_obj)

# prepare a directory where to write data
output_dir = os.path.join(os.getcwd(), args.output_dir)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# write data to files
if args.output_type == 'j':
    with open(f'{output_dir}/proxy_{now}.json', 'w') as file:
        file.write(df.to_json())
elif args.output_type == 'm':
    mpl.style.use('ggplot')
    df.plot(kind='bar',x='kernel type',y='HW GFLOPS', rot=45)
    plt.savefig(f'{output_dir}/proxy_hw-gflops_{now}.png', dpi=300, bbox_inches='tight')

    df.plot(kind='bar',x='kernel type',y='NZ GFLOPS', rot=45)
    plt.savefig(f'{output_dir}/proxy_nz-gflops_{now}.png', dpi=300, bbox_inches='tight')

    df.plot(kind='bar',x='kernel type',y='time', rot=45)
    plt.savefig(f'{output_dir}/proxy_time_{now}.png', dpi=300, bbox_inches='tight')
elif args.output_type == 'p':
    with open(f'{output_dir}/proxy_{now}.pd', 'w') as file:
        file.write(df.to_csv())

# print out to console
print(df)
