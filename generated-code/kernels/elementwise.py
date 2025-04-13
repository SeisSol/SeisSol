# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

from kernels.common import generate_kernel_name_prefix
from yateto import Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile


def addKernels(generator, aderdg, matricesDir, targets=["cpu"]):
    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)

        elemwise = parseJSONMatrixFile(
            f"{matricesDir}/elemwise-collocate-p{aderdg.order}.json"
        )

        ew_extraweight_M = Tensor(
            f"ew_extraweight_M", elemwise.ew_quad_weights_v.shape()
        )
        ew_extraweight_k = Tensor(
            f"ew_extraweight_k", (*elemwise.ew_quad_weights_v.shape(), 3, 3)
        )
        ew_extraweight_r = Tensor(
            f"ew_extraweight_r", (*elemwise.ew_quad_weights_f.shape(), 4)
        )

        ew_M = Tensor(
            f"ew_M",
            (aderdg.numberOf3DBasisFunctions(), aderdg.numberOf3DBasisFunctions()),
        )
        ew_k = Tensor(
            f"ew_k",
            (aderdg.numberOf3DBasisFunctions(), aderdg.numberOf3DBasisFunctions(), 3),
        )
        ew_kT = Tensor(
            f"ew_kT",
            (aderdg.numberOf3DBasisFunctions(), aderdg.numberOf3DBasisFunctions(), 3),
        )
        ew_r = Tensor(
            f"ew_r",
            (aderdg.numberOf3DBasisFunctions(), aderdg.numberOf2DBasisFunctions(), 4),
        )

        generator.add(
            f"{name_prefix}bootstrap",
            [
                ew_M["ij"]
                <= elemwise.ew_collocate_f_vv["wi"]
                * elemwise.ew_collocate_f_vv["wj"]
                * ew_extraweight_M["w"]
                * elemwise.ew_quad_weights_v["w"],
                ew_k["ijd"]
                <= elemwise.ew_collocate_df_vv["wie"]
                * elemwise.ew_collocate_f_vv["wj"]
                * ew_extraweight_k["wed"]
                * elemwise.ew_quad_weights_v["w"],
                ew_kT["ijd"]
                <= (-1)
                * elemwise.ew_collocate_f_vv["wi"]
                * elemwise.ew_collocate_df_vv["wje"]
                * ew_extraweight_k["wed"]
                * elemwise.ew_quad_weights_v["w"],
                ew_r["ijf"]
                <= elemwise.ew_collocate_f_vf["wif"]
                * elemwise.ew_collocate_f_ff["wj"]
                * ew_extraweight_r["wf"]
                * elemwise.ew_quad_weights_f["w"],
            ],
        )


def includeTensors(includeTensors, aderdg, matricesDir):
    vtkbase = parseJSONMatrixFile(
        f"{matricesDir}/elemwise-collocate-p{aderdg.order}.json"
    )
    for x in vtkbase.__dict__:
        if isinstance(vtkbase.__dict__[x], dict):
            for y in vtkbase.__dict__[x]:
                includeTensors.add(vtkbase.__dict__[x][y])
        else:
            includeTensors.add(vtkbase.__dict__[x])
