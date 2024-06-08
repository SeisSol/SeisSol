import tensorforge.functions as forge
from tensorforge import Tensor

def test1(generator):
    A = Tensor('Atest', (32,32))
    B = Tensor('Btest', (32,32))
    C = Tensor('Ctest', (32,32))
    D = Tensor('Dtest', (32,32))
    generator.add(name='test1', ast=forge.scalarblock([
        forge.conditional(C['ij'], [forge.loop(forge.tensor(A['ij']) <= forge.tensor(C['ij']),
            [forge.assign(A['ij'], forge.abs(forge.exp(B['ij']) * forge.tensor(D['ij'])))]
        )]),
    ]), target='gpu')
