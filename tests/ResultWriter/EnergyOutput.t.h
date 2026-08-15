// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Energy.h"
#include "ResultWriter/EnergyOutput.h"

#include <array>
#include <cstddef>
#include <string_view>

namespace seissol::unit_test {

using seissol::model::EnergyDescriptor;
using seissol::model::EnergyUnit;
using seissol::writer::EnergiesStorage;

TEST_CASE("EnergiesStorage handles") {
  EnergiesStorage storage;
  storage.setSimcount(1);

  const auto first = storage.addEnergy(EnergyDescriptor{"first"});
  const auto second = storage.addEnergy(EnergyDescriptor{"second"});
  const auto third = storage.addEnergy(EnergyDescriptor{"third"});

  SUBCASE("handles are consecutive and follow registration order") {
    REQUIRE(first == 0);
    REQUIRE(second == 1);
    REQUIRE(third == 2);
    REQUIRE(storage.descriptors().size() == 3);
    // registration order, not the alphabetical order of the name -> handle map
    REQUIRE(storage.descriptors()[0].name == "first");
    REQUIRE(storage.descriptors()[1].name == "second");
    REQUIRE(storage.descriptors()[2].name == "third");
  }

  SUBCASE("lookup by name and by handle agree") {
    storage.energy(second, 0) = 17.0;
    REQUIRE(storage.energy("second", 0) == 17.0);
    storage.energy("third", 0) = -3.5;
    REQUIRE(storage.energy(third, 0) == -3.5);
  }

  SUBCASE("has() reports registration, and does not create entries") {
    REQUIRE(storage.has("first"));
    REQUIRE_FALSE(storage.has("frist"));
    REQUIRE(storage.descriptors().size() == 3);
  }

  SUBCASE("reset clears the values but keeps the handles") {
    storage.energy(first, 0) = 1.0;
    storage.energy(second, 0) = 2.0;
    storage.reset();
    REQUIRE(storage.energy(first, 0) == 0.0);
    REQUIRE(storage.energy(second, 0) == 0.0);
    REQUIRE(storage.descriptors().size() == 3);
    REQUIRE(storage.has("first"));
  }
}

TEST_CASE("EnergiesStorage multi-simulation indexing") {
  constexpr std::size_t Simcount = 4;
  EnergiesStorage storage;
  storage.setSimcount(Simcount);

  const auto alpha = storage.addEnergy(EnergyDescriptor{"alpha"});
  const auto beta = storage.addEnergy(EnergyDescriptor{"beta"});

  for (std::size_t sim = 0; sim < Simcount; ++sim) {
    storage.energy(alpha, sim) = 10.0 + static_cast<double>(sim);
    storage.energy(beta, sim) = 100.0 + static_cast<double>(sim);
  }

  SUBCASE("simulations do not alias each other") {
    for (std::size_t sim = 0; sim < Simcount; ++sim) {
      REQUIRE(storage.energy(alpha, sim) == 10.0 + static_cast<double>(sim));
      REQUIRE(storage.energy("beta", sim) == 100.0 + static_cast<double>(sim));
    }
  }

  SUBCASE("the value buffer is sized handles x simulations") {
    REQUIRE(storage.values().size() == 2 * Simcount);
  }

  SUBCASE("reset clears every simulation") {
    storage.reset();
    for (std::size_t sim = 0; sim < Simcount; ++sim) {
      REQUIRE(storage.energy(alpha, sim) == 0.0);
      REQUIRE(storage.energy(beta, sim) == 0.0);
    }
  }
}

TEST_CASE("Energy descriptors of the configured material are well formed") {
  using Compute = seissol::model::EnergyCompute<seissol::model::MaterialT>;

  SUBCASE("EnergyCount matches the descriptor list") {
    REQUIRE(Compute::EnergyCount == Compute::Energies.size());
    REQUIRE(Compute::EnergyCount > 0);
  }

  SUBCASE("every descriptor is named and unique") {
    for (std::size_t i = 0; i < Compute::Energies.size(); ++i) {
      REQUIRE_FALSE(Compute::Energies[i].name.empty());
      for (std::size_t j = i + 1; j < Compute::Energies.size(); ++j) {
        REQUIRE(Compute::Energies[i].name != Compute::Energies[j].name);
      }
    }
  }

  SUBCASE("every group carries exactly one heading") {
    for (const auto& descriptor : Compute::Energies) {
      if (descriptor.group.empty()) {
        REQUIRE(descriptor.groupLabel.empty());
        REQUIRE(descriptor.shortLabel.empty());
        continue;
      }
      std::size_t headings = 0;
      for (const auto& other : Compute::Energies) {
        if (other.group == descriptor.group && !other.groupLabel.empty()) {
          ++headings;
        }
      }
      REQUIRE(headings == 1);
    }
  }

  SUBCASE("members of a group share a unit") {
    // a group is printed as one total plus per-member shares, so summing its
    // members has to make physical sense
    for (const auto& descriptor : Compute::Energies) {
      if (descriptor.group.empty()) {
        continue;
      }
      for (const auto& other : Compute::Energies) {
        if (other.group == descriptor.group) {
          REQUIRE(other.unit == descriptor.unit);
        }
      }
    }
  }

  SUBCASE("registering the material's energies yields distinct handles") {
    EnergiesStorage storage;
    storage.setSimcount(1);
    for (const auto& descriptor : Compute::Energies) {
      storage.addEnergy(descriptor);
    }
    REQUIRE(storage.descriptors().size() == Compute::EnergyCount);
    for (const auto& descriptor : Compute::Energies) {
      REQUIRE(storage.has(descriptor.name));
    }
  }
}

TEST_CASE("The stored potential energy is fully accounted for") {
  using Compute = seissol::model::EnergyCompute<seissol::model::MaterialT>;

  const auto groupOf = [](std::string_view name) -> std::string_view {
    for (const auto& descriptor : Compute::Energies) {
      if (descriptor.name == name) {
        return descriptor.group;
      }
    }
    return {};
  };
  const auto declares = [](std::string_view name) {
    for (const auto& descriptor : Compute::Energies) {
      if (descriptor.name == name) {
        return true;
      }
    }
    return false;
  };

  // Everything that is part of the stored potential must sit in the same group
  // as the kinetic energy it is reported against. For a viscoelastic material
  // that includes the Maxwell branch springs: leaving them out would make the
  // printed potential share understate the actual potential.
  if (declares("elastic_kinetic_energy") && declares("viscoelastic_energy")) {
    REQUIRE(groupOf("viscoelastic_energy") == groupOf("elastic_kinetic_energy"));
    REQUIRE(groupOf("elastic_energy") == groupOf("elastic_kinetic_energy"));
  }

  // A dissipation rate is a flux, not a stored energy, and must never be summed
  // into an energy group.
  for (const auto& descriptor : Compute::Energies) {
    if (descriptor.unit == seissol::model::EnergyUnit::Power) {
      for (const auto& other : Compute::Energies) {
        if (other.group == descriptor.group) {
          REQUIRE(other.unit == seissol::model::EnergyUnit::Power);
        }
      }
    }
  }
}

} // namespace seissol::unit_test
