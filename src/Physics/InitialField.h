#ifndef PHYSICS_INITIALFIELD_H
#define PHYSICS_INITIALFIELD_H

#include <vector>
#include <array>
#include <complex>
#include <Kernels/precision.hpp>
#include <generated_code/init.h>

namespace seissol {
  namespace physics {
    class InitialField {
    public:
      virtual void evaluate(  double time,
                              std::vector<std::array<double, 3>> const& points,
                              init::dofsQP::view::type& dofsQP ) const = 0;
    };

    class ZeroField : public InitialField {
    public:
      void evaluate(  double,
                      std::vector<std::array<double, 3>> const&,
                      init::dofsQP::view::type& dofsQP ) const
      {
        dofsQP.setZero();
      }
    };

    class Planarwave : public InitialField {
    public:
      Planarwave();

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      init::dofsQP::view::type& dofsQP ) const;
    private:
      std::complex<real> m_eigenvectors[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
      int const                                             m_setVar;
      std::vector<int>                                      m_varField;
      std::vector<std::complex<real>>                       m_ampField;
      std::array<double, 3>                                 m_kVec;
      std::array<std::complex<real>, NUMBER_OF_QUANTITIES>  m_lambdaA;
    };
  }
}

#endif // PHYSICS_INITIALFIELD_H
