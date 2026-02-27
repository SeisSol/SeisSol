# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from yateto import Tensor
from yateto.memory import DenseMemoryLayout


def generate_kernel_name_prefix(target):
    return f"{target}_" if target == "gpu" else ""


def tensor_to_numpy(tensor):
    np_tensor = np.zeros(tensor.shape())
    for indices in tensor.values():
        np_tensor[indices] = np.float64(tensor._values[indices])
    return np_tensor


def numpy_to_tensor(
    name,
    np_array,
    memoryLayoutClass=DenseMemoryLayout,
    alignStride=False,
    namespace=None,
):

    spp = {}
    nnz_count = np.count_nonzero(np_array)
    nnz_indices = np.array(np.nonzero(np_array))
    for i in range(nnz_count):
        index = tuple(nnz_indices[:, i])
        spp[index] = str(np_array[index])

    return Tensor(name, np_array.shape, spp, memoryLayoutClass, alignStride, namespace)
