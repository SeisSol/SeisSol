from tensorforge import Tensor
from tensorforge.input import parseXMLMatrixFile, parseJSONMatrixFile

def addStiffnessTensor(generator):
    stiffnessTensor = Tensor('stiffnessTensor', (3, 3, 3, 3))
    direction = Tensor('direction', (3,))
    christoffel = Tensor('christoffel', (3,3))

    computeChristoffel = christoffel['ik'] <= stiffnessTensor['ijkl'] * direction['j'] * direction['l']
    generator.add('computeChristoffel', computeChristoffel)

def includeMatrices(matricesDir):
    matrices = parseJSONMatrixFile('{}/sampling_directions.json'.format(matricesDir))
    return set([matrices.samplingDirections])
