#ifndef EQUATIONS_SETUP_H_
#define EQUATIONS_SETUP_H_

// Gather all Setup Headers here
#include "Equations/anisotropic/Model/AnisotropicSetup.h" // IWYU pragma: keep
#include "Equations/elastic/Model/ElasticSetup.h"         // IWYU pragma: keep
#include "Equations/acoustic/Model/AcousticSetup.h"       // IWYU pragma: keep
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/PoroelasticSetup.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/ViscoelasticSetup.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/ViscoelasticSetup.h" // IWYU pragma: keep
#endif

#endif
