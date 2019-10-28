#ifndef PHYSICS_INITIALFIELD_H
#define PHYSICS_INITIALFIELD_H

#include <vector>
#include <array>
#include <complex>
#include <Kernels/precision.hpp>
#include <generated_code/init.h>
#ifdef USE_POROELASTIC
#include <Equations/poroelastic/Model/datastructures.hpp>
#endif

namespace seissol {
  namespace physics {
    class InitialField {
    public:
      virtual ~InitialField() {}
      virtual void evaluate(  double time,
                              std::vector<std::array<double, 3>> const& points,
                              yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const = 0;
    };

    class ZeroField : public InitialField {
    public:
      void evaluate(  double,
                      std::vector<std::array<double, 3>> const&,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const
      {
        dofsQP.setZero();
      }
    };

    class Planarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
#ifdef USE_POROELASTIC
      Planarwave(seissol::model::Material material, real phase = 0.0);
#elif
      Planarwave(real phase = 0.0);
#endif

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      std::complex<real> m_eigenvectors[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
      int const                                             m_setVar;
      std::vector<int>                                      m_varField;
      std::vector<std::complex<real>>                       m_ampField;
      std::array<double, 3>                                 m_kVec;
      std::array<std::complex<real>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      real                                                  m_phase;
    };
  }
}

#endif // PHYSICS_INITIALFIELD_H
