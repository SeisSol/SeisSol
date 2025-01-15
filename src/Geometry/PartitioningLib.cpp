// SPDX-FileCopyrightText: 2024 SeisSol Group
// SPDX-FileCopyrightText: 2023 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "PartitioningLib.h"
#include "Common/Fnv1a.h"
#include <PUML/Partition.h>
#include <string_view>

using PUML::PartitionerType;
using namespace std::literals;
using namespace seissol::literals;

namespace seissol {

PartitionerType toPartitionerType(std::string_view partitioningLib) {
  switch (fnv1a(partitioningLib)) {
  case "Default"_fnv1a:
#if defined(USE_PARMETIS)
  case "Parmetis"_fnv1a:
    return PartitionerType::Parmetis;
  case "ParmetisGeometric"_fnv1a:
    return PartitionerType::ParmetisGeometric;
#endif
#if defined(USE_PTSCOTCH)
  case "PtScotch"_fnv1a:
    return PartitionerType::PtScotch;
  case "PtScotchQuality"_fnv1a:
    return PartitionerType::PtScotchQuality;
  case "PtScotchBalance"_fnv1a:
    return PartitionerType::PtScotchBalance;
  case "PtScotchBalanceQuality"_fnv1a:
    return PartitionerType::PtScotchBalanceQuality;
  case "PtScotchSpeed"_fnv1a:
    return PartitionerType::PtScotchSpeed;
  case "PtScotchBalanceSpeed"_fnv1a:
    return PartitionerType::PtScotchBalanceSpeed;
#endif
#if defined(USE_PARHIP)
  case "ParHIPUltrafastMesh"_fnv1a:
    return PartitionerType::ParHIPUltrafastMesh;
  case "ParHIPFastMesh"_fnv1a:
    return PartitionerType::ParHIPFastMesh;
  case "ParHIPEcoMesh"_fnv1a:
    return PartitionerType::ParHIPEcoMesh;
  case "ParHIPUltrafastSocial"_fnv1a:
    return PartitionerType::ParHIPUltrafastSocial;
  case "ParHIPFastSocial"_fnv1a:
    return PartitionerType::ParHIPFastSocial;
  case "ParHIPEcoSocial"_fnv1a:
    return PartitionerType::ParHIPEcoSocial;
#endif
  default:
    break;
  }
  return PartitionerType::None;
}

std::string_view toStringView(PartitionerType type) {
  switch (type) {
#if defined(USE_PARMETIS)
  case PartitionerType::Parmetis:
    return "Parmetis"sv;
  case PartitionerType::ParmetisGeometric:
    return "ParmetisGeometric"sv;
#endif
#if defined(USE_PTSCOTCH)
  case PartitionerType::PtScotch:
    return "PtScotch"sv;
  case PartitionerType::PtScotchQuality:
    return "PtScotchQuality"sv;
  case PartitionerType::PtScotchBalance:
    return "PtScotchBalance"sv;
  case PartitionerType::PtScotchBalanceQuality:
    return "PtScotchBalanceQuality"sv;
  case PartitionerType::PtScotchSpeed:
    return "PtScotchSpeed"sv;
  case PartitionerType::PtScotchBalanceSpeed:
    return "PtScotchBalanceSpeed"sv;
#endif
#if defined(USE_PARHIP)
  case PartitionerType::ParHIPUltrafastMesh:
    return "ParHIPUltrafastMesh"sv;
  case PartitionerType::ParHIPFastMesh:
    return "ParHIPFastMesh"sv;
  case PartitionerType::ParHIPEcoMesh:
    return "ParHIPEcoMesh"sv;
  case PartitionerType::ParHIPUltrafastSocial:
    return "ParHIPUltrafastSocial"sv;
  case PartitionerType::ParHIPFastSocial:
    return "ParHIPFastSocial"sv;
  case PartitionerType::ParHIPEcoSocial:
    return "ParHIPEcoSocial"sv;
#endif
  default:
    break;
  }
  return "None"sv;
}

} // namespace seissol
