#ifndef EQUATIONS_SETUP_H_
#define EQUATIONS_SETUP_H_

//Gather all Setup Headers here
#include "Equations/anisotropic/Model/AnisotropicSetup.h"
#include "Equations/elastic/Model/ElasticSetup.h"
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/PoroelasticSetup.h"
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/ViscoelasticSetup.h"
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/ViscoelasticSetup.h"
#endif

#endif
