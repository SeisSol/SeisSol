#ifndef PHYSICS_INITIALFIELD_H
#define PHYSICS_INITIALFIELD_H

#include <vector>
#include <array>
#include <Numerical_aux/MatrixView.h>

namespace seissol {
  namespace physics {
    class InitialField {
    public:
      virtual void evaluate(  double time,
                              std::vector<std::array<double, 3>> const& points,
                              MatrixView dofsQP ) const = 0;
    };

    class ZeroField : public InitialField {
    public:
      void evaluate(  double,
                      std::vector<std::array<double, 3>> const&,
                      MatrixView dofsQP ) const
      {
        dofsQP.setZero();
      }
    };

    class PlanarwaveElastic : public InitialField {
    public:
      PlanarwaveElastic();

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      MatrixView dofsQP ) const;
    private:
      real m_eigenvectors[9*9];
      int m_setVar;
      std::array<int, 2> m_varField;
      std::array<double, 2> m_ampField;
      std::array<double, 3> m_kVec;
      std::array<double, 9> m_lambdaA;
      double m_kVecNorm;
    };
  }
}

#endif // PHYSICS_INITIALFIELD_H
