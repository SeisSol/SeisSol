#ifndef PHYSICS_INITIALFIELD_H
#define PHYSICS_INITIALFIELD_H

#include <vector>
#include <array>
#include <complex>
#include "Initializer/typedefs.hpp"
#include <Kernels/precision.hpp>
#include <generated_code/init.h>

namespace seissol {
  namespace physics {
    class InitialField {
    public:
      virtual ~InitialField() {}
      virtual void evaluate(double time,
                            std::vector<std::array<double, 3>> const& points,
                            const CellMaterialData& materialData,
                            yateto::DenseTensorView<2,real,unsigned>& dofsQP) const = 0;
    };

    class ZeroField : public InitialField {
    public:
      void evaluate(double,
                    std::vector<std::array<double, 3>> const&,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const {
        dofsQP.setZero();
      }
    };

    class Planarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      Planarwave(model::Material material, double phase = 0.0);

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      const CellMaterialData& materialData,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      const std::vector<int>                                        m_varField;
      const std::vector<std::complex<double>>                       m_ampField;
      const double                                                  m_phase;
      const std::array<double, 3>                                   m_kVec;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      std::complex<double> m_eigenvectors[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
    };

    class ScholteWave : public InitialField {
    public:
      ScholteWave() {

      }
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const;
    };
    class SnellsLaw : public InitialField {
    public:
      SnellsLaw() {

      }
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const;
    };
      class Ocean : public InitialField {
      public:
          Ocean() {

          }
          void evaluate(double time,
                        std::vector<std::array<double, 3>> const& points,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2,real,unsigned>& dofsQP) const;
      };
  }
}

#endif // PHYSICS_INITIALFIELD_H
