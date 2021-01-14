import argparse
import subprocess
import re
import numpy as np

argparser = argparse.ArgumentParser()
argparser.add_argument('equation')
args = argparser.parse_args()

binary = 'build/SeisSol_Release_dhsw_6_{}'.format(args.equation)
parameter_file = 'src/tests/convergence/parameters_{}_4.par'.format(args.equation)

SeisSol_cmd = "{} {}".format(binary, parameter_file)

result = subprocess.check_output(SeisSol_cmd, shell=True)

s = open('src/tests/convergence/expected_errors.txt', 'r').read()
expected_errors = eval(s)

for line in result.splitlines():
    err = re.search(
        "(LInf|L2|L1)\s*,\s+var\[\s*(\d+)\s*\]\s*=\s*([0-9\.eE+-]+)",
        line.decode('utf-8'),
        )
    if err:
        norm = err.groups()[0]
        quantity = err.groups()[1]
        e = float(err.groups()[2])
        exp_e = expected_errors[args.equation]['{},{}'.format(norm,quantity)]
        assert(exp_e == e, "expected {}-error of quantity {} to be {}, but was {}.".format(
            norm, quantity, exp_e, e))
    
