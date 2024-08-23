#ifndef EQUATIONS_DATASTRUCTURES_H_
#define EQUATIONS_DATASTRUCTURES_H_

// Gather all datastructure Headers here
#include "Equations/anisotropic/Model/datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"
#include "Equations/poroelastic/Model/datastructures.hpp"
#include "Equations/viscoelastic2/Model/datastructures.hpp"

#include "Equations/anisotropic/Model/integrationData.hpp"
#include "Equations/elastic/Model/integrationData.hpp"
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/integrationData.hpp"
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/integrationData.hpp"
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/integrationData.hpp"
#endif

namespace seissol::model {
#if defined(USE_ANISOTROPIC)
using MaterialT = AnisotropicMaterial;
#elif defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
using MaterialT = ViscoElasticMaterial;
#elif defined(USE_ELASTIC)
using MaterialT = ElasticMaterial;
#elif defined(USE_POROELASTIC)
using MaterialT = PoroElasticMaterial;
#else
#error "Material class unknown."
#endif
} // namespace seissol::model

#endif
