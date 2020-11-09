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

    //A planar wave travelling in direction kVec
    class Planarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
    Planarwave(const CellMaterialData& materialData, double phase = 0.0, std::array<double, 3> kVec = {M_PI, M_PI, M_PI});

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      const CellMaterialData& materialData,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      const std::vector<int>                            m_varField;
      const std::vector<std::complex<double>>           m_ampField;
      const double                                      m_phase;
      const std::array<double, 3>                       m_kVec;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES> m_eigenvectors;
    };

    //superimpose three planar waves travelling into different directions
    class SuperimposedPlanarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      SuperimposedPlanarwave(const CellMaterialData& materialData, real phase = 0.0);

      void evaluate(  double time,
                      std::vector<std::array<double, 3>> const& points,
                      const CellMaterialData& materialData,
                      yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    private:
      const std::array<std::array<double, 3>, 3>  m_kVec;
      const double                                m_phase;
      std::array<Planarwave, 3>             m_pw;
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
