#ifndef EQUATIONS_DATASTRUCTURES_H_
#define EQUATIONS_DATASTRUCTURES_H_

//Gather all datastructure Headers here
#include "Equations/anisotropic/Model/datastructures.hpp"
#include "Equations/damaged-elastic/Model/datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"
#include "Equations/poroelastic/Model/datastructures.hpp"
#include "Equations/viscoelastic2/Model/datastructures.hpp"

#include "Equations/anisotropic/Model/integrationData.hpp"
#include "Equations/damaged-elastic/Model/integrationData.hpp"
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
using Material_t = AnisotropicMaterial;
#elif defined(USE_VISCOELASTIC) || defined (USE_VISCOELASTIC2)
using Material_t = ViscoElasticMaterial;
#elif defined(USE_ELASTIC)
using Material_t = ElasticMaterial;
#elif defined(USE_POROELASTIC)
using Material_t = PoroElasticMaterial;
#elif defined(USE_DAMAGEDELASTIC)
using Material_t = DamagedElasticMaterial;
#else
#error "Material class unknown."
#endif
}

#endif
