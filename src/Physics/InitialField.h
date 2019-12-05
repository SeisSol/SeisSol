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
      Planarwave(double phase = 0.0, std::array<double, 3> kVec = {M_PI, M_PI, M_PI});

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      const CellMaterialData& materialData,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      std::complex<double> m_eigenvectors[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
      int                                                     m_setVar;
      std::vector<int>                                        m_varField;
      std::vector<std::complex<double>>                       m_ampField;
      std::array<double, 3>                                   m_kVec;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      double                                                  m_phase;
    };

    class AnisotropicPlanarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      AnisotropicPlanarwave(real phase = 0.0);

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      const CellMaterialData& materialData,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      std::array<std::array<double, 3>, 3>  m_kVec;
      std::array<Planarwave, 3>             m_pw;
      double                                m_phase;
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
