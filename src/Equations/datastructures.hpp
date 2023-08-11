#ifndef EQUATIONS_DATASTRUCTURES_H_
#define EQUATIONS_DATASTRUCTURES_H_

#include "Common/constants.hpp"

//Gather all datastructure Headers here
#include "Equations/anisotropic/Model/datastructures.hpp"
#include "Equations/elastic/Model/datastructures.hpp"
#include "Equations/poroelastic/Model/datastructures.hpp"
#include "Equations/viscoelastic/Model/datastructures.hpp"

#include "Equations/anisotropic/Model/integrationData.hpp"
#include "Equations/elastic/Model/integrationData.hpp"
#include "Equations/poroelastic/Model/integrationData.hpp"
#include "Equations/viscoelastic/Model/integrationData.hpp"

namespace seissol::model {
#if defined(USE_ANISOTROPIC)
using Material_t = AnisotropicMaterial;
#elif defined(USE_VISCOELASTIC) || defined (USE_VISCOELASTIC2)
using Material_t = ViscoElasticMaterial;
#elif defined(USE_ELASTIC)
using Material_t = ElasticMaterial;
#elif defined(USE_POROELASTIC)
using Material_t = PoroElasticMaterial;
#else
#error "Material class unknown."
#endif
}

#endif
