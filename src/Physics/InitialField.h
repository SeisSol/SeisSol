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
      virtual ~InitialField() = default;
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
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override {
        dofsQP.setZero();
      }
    };

    //A planar wave travelling in direction kVec
    class Planarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
    Planarwave(const CellMaterialData& materialData, 
               double phase = 0.0, 
               std::array<double, 3> kVec = {M_PI, M_PI, M_PI},
#ifndef USE_POROELASTIC
               //Elastic materials have the following wave modes:
               //-P, -S2, -S1, N, N, N, S1, S2, P
               //Here we impose the -S2 and P mode
               std::vector<int> varField = {1,8}, 
               std::vector<std::complex<double>> ampField = {1.0, 1.0}
#else
               //Poroelastic materials have the following wave modes:
               //-P, -S2, -S1, -Ps, N, N, N, N, N, Ps, S1, S2, P
               //Here we impose -S1, -Ps and P
               std::vector<int> varField = {2,3,12}, 
               std::vector<std::complex<double>> ampField = {1.0, 1.0, 1.0}
#endif
               );

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const CellMaterialData& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const override;
    protected:
      const std::vector<int>                                        m_varField;
      const std::vector<std::complex<double>>                       m_ampField;
      const double                                                  m_phase;
      const std::array<double, 3>                                   m_kVec;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES> m_eigenvectors;
    };

    //superimpose three planar waves travelling into different directions
    class SuperimposedPlanarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      SuperimposedPlanarwave(const CellMaterialData& materialData, real phase = 0.0);

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const CellMaterialData& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const override;
    private:
      const std::array<std::array<double, 3>, 3>  m_kVec;
      const double                                m_phase;
      std::array<Planarwave, 3>                   m_pw;
    };

    //A part of a planar wave travelling in one direction
    class TravellingWave : public Planarwave{
    public:
      TravellingWave(const CellMaterialData& materialData, const TravellingWaveParameters& travellingWaveParameters);

      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
      private:
      std::array<double, 3> m_origin;
    };

    class ScholteWave : public InitialField {
    public:
      ScholteWave() = default;
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
    };
    class SnellsLaw : public InitialField {
    public:
      SnellsLaw() = default;
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const CellMaterialData& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
    };
    /*
     * From
     * Abrahams, L. S., Krenz, L., Dunham, E. M., & Gabriel, A. A. (2019, December).
     * Verification of a 3D fully-coupled earthquake and tsunami model.
     * In AGU Fall Meeting Abstracts (Vol. 2019, pp. NH43F-1000).
     * A 3D extension of the 2D scenario in
     * Lotto, G. C., & Dunham, E. M. (2015).
     * High-order finite difference modeling of tsunami generation in a compressible ocean from offshore earthquakes.
     * Computational Geosciences, 19(2), 327-340.
     */
      class Ocean : public InitialField {
      private:
        const int mode;
        const double gravitationalAcceleration;
      public:
          Ocean(int mode, double gravitationalAcceleration);
          void evaluate(double time,
                        std::vector<std::array<double, 3>> const& points,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
      };
  }
}

#endif // PHYSICS_INITIALFIELD_H
