import numpy as np
from yateto import Tensor
from yateto.input import parseJSONMatrixFile
from multSim import OptionalDimTensor
from common import generate_kernel_name_prefix
from yateto import Tensor, Scalar, simpleParameterSpace

def addKernels(generator, aderdg, matricesDir, targets=['cpu']):
    for target in targets:
        name_prefix = generate_kernel_name_prefix(targets)

        vtkbase = parseJSONMatrixFile(f'{matricesDir}/vtkbase.json')
        vtko = parseJSONMatrixFile(f'{matricesDir}/vtko{aderdg.order}.json')

        # the following is due to a shortcut in Yateto where 1-column matrices are interpreted as rank-1 vectors

        qb = Tensor('qb', (aderdg.numberOf3DBasisFunctions(),))
        xv = [Tensor(f'xv({i})', (((i+1)*(i+2)*(i+3))//6,)) for i in range(8)]
        xf = [Tensor(f'xf({i})', (((i+1)*(i+2))//2,)) for i in range(8)]

        generator.addFamily(f'{name_prefix}projectBasisToVtkVolume',
                  simpleParameterSpace(8),
                  lambda i: xv[i]['p'] <= vtko.byName(f'coll{aderdg.order}vd{i}')['pb'] * qb['b'],
                  target=target)
        
        generator.addFamily(f'{name_prefix}projectBasisToVtkFace',
                    simpleParameterSpace(8, 4),
                  lambda i,j: xf[i]['p'] <= vtko.byName(f'coll{aderdg.order}f{j}d{i}')['pb'] * qb['b'],
                  target=target)

def includeTensors(matricesDir, includeTensors):
    vtkbase = parseJSONMatrixFile(f'{matricesDir}/vtkbase.json')
    for x in vtkbase.__dict__:
        includeTensors.add(vtkbase.__dict__[x])
