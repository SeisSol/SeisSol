#!/usr/bin/env python

import opentuner
from opentuner import ConfigurationManipulator
from opentuner import SwitchParameter
from opentuner import MeasurementInterface
from opentuner import Result
from tune import MemoryLayout
import os
import re

argparser = opentuner.default_argparser()
argparser.add_argument( '--equations', required=True,
                        help='Elastic, viscoelastic, or viscoelastic2')
argparser.add_argument('--order', required=True, type=int,
                        help='Convergence order of ADER-DG scheme')
argparser.add_argument( '--numberOfMechanisms', type=int, default=0,
                        help='Required for viscoelastic')
argparser.add_argument( '--arch', required=True,
                        help='{s,d} x {snb,hsw,knl}')
argparser.add_argument( '--nelem', default=10000, type=int,
                        help='Number of elements used in measurement')
argparser.add_argument( '--ntimesteps', default=100, type=int,
                        help='Number of timesteps used in measurement')
argparser.add_argument('--ncompileJobs', default=4, type=int,
                        help='Number of threads used for compilation')

class MemoryLayoutTuner(MeasurementInterface):
  def __init__(self, args):
    super(MemoryLayoutTuner, self).__init__(args)
    if self.args.equations == 'elastic':
      self.memoryLayouts = MemoryLayout.getElasticMemoryLayouts(self.args.order, self.args.arch)
    elif self.args.equations == 'viscoelastic':
      self.memoryLayouts = MemoryLayout.getViscoelasticMemoryLayouts(self.args.order, self.args.numberOfMechanisms, self.args.arch)
    elif self.args.equations == 'viscoelastic2':
      self.memoryLayouts = MemoryLayout.getViscoelastic2MemoryLayouts(self.args.order, self.args.arch)

    options = {
      'equations':          args.equations,
      'order':              args.order,
      'numberOfMechanisms': args.numberOfMechanisms,
      'arch':               args.arch,
      'compileMode':        'release',
      '--jobs':             args.ncompileJobs
    }
    self.options = ' '.join(['{}={}'.format(key, value) for key, value in options.items()])
  
  def manipulator(self):
    """
    Define the search space by creating a
    ConfigurationManipulator
    """
    manipulator = ConfigurationManipulator()
    for name, options in self.memoryLayouts[1].iteritems():
      manipulator.add_parameter(
        SwitchParameter(name, 1+len(options))
      )
    return manipulator

  def run(self, desired_result, input, limit):
    """
    Compile and run a given configuration then
    return performance
    """
    layoutFileName = 'layout.xml'

    cfg = desired_result.configuration.data
    self.cfg_to_layout(desired_result.configuration.data, layoutFileName)
    
    cwd = os.getcwd()
    os.chdir(os.path.dirname(os.path.realpath(__file__)) + '/../proxy/')
    sconsCmd = 'scons {} memLayout={} buildDir={} tools=none'.format(
                  self.options,
                  os.path.join(cwd, layoutFileName),
                  cwd )
    compile_result = self.call_program(sconsCmd)
    os.chdir(cwd)

    run_cmd = './bin/seissol_proxy {} {} all'.format(self.args.nelem, self.args.ntimesteps)

    run_result = self.call_program(run_cmd)
    assert run_result['returncode'] == 0
    time = re.search(r'time for seissol proxy\W*:\W*([0-9\.\+e]+)', run_result['stdout'])
    assert time != None

    return Result(time=float(time.group(1)))

  def save_final_config(self, configuration):
    """called at the end of tuning"""
    print "Optimal block size written to tuned_layout.xml."
    self.cfg_to_layout(configuration.data, 'tuned_layout.xml')

  def cfg_to_layout(self, cfg, layoutFileName):
    cfgMinusDense = {key: value-1 for key, value in cfg.iteritems() if value > 0}
    MemoryLayout.generateLayoutFile(layoutFileName, cfgMinusDense, self.memoryLayouts)


if __name__ == '__main__':
  MemoryLayoutTuner.main(argparser.parse_args())
