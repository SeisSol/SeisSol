from yateto.input import parseXMLMatrixFile, parseJSONMatrixFile

def main():
    alignStride = lambda name: True
    data = parseXMLMatrixFile('{}/matrices_{}.xml'.format('matrices', 56), alignStride=alignStride)

    limits = [56, 35, 20, 10, 4, 1]

    for o in range(5):
        print(f"""__device__ __forceinline__ void dgkernel{limits[o]}(real* __restrict__ temp, const real* __restrict__ dQ) {{
for (int j = 0; j < Quantities; ++j) {{
    real column[{limits[o]}];
    #pragma unroll
    for (int i = 0; i < {limits[o]}; ++i) {{
        column[i] = dQ[MATRIX(i, j, {limits[o]}, Quantities)];
    }}
""")
        for d in range(3):
            derivative = data.kDivMT[d]
            values = derivative.values()
            for k in range(limits[o+1]):
                entries = []
                for v in values:
                    if k == v[0] and v[1] < limits[o]:
                        entries += [f'{values[v]}f * column[{v[1]}]']
                if limits[o+1] >= 10:
                    print(f'    temp[MATRIX({k}, j + Quantities * {d}, {limits[o+1]}, Quantities * 3)] = {" + ".join(entries)};')
                else:
                    print(f'    temp[{k} * Quantities * 3 + (j + Quantities * {d})] = {" + ".join(entries)};')
        print('    }')
        print('}')


if __name__ == '__main__':
    main()
