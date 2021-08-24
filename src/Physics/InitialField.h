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
    Planarwave(const CellMaterialData& materialData, 
               double phase = 0.0, 
               std::array<double, 3> kVec = {M_PI, M_PI, M_PI},
               std::vector<int> varField = {1, 8}, 
               std::vector<std::complex<double>> ampField = {1.0, 1.0});

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const CellMaterialData& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
    protected:
      const double                                      m_phase;
      const std::array<double, 3>                       m_kVec;
      const std::vector<int>                            m_varField;
      const std::vector<std::complex<double>>           m_ampField;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      std::complex<double> m_eigenvectors[NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES];
    };

    //superimpose three planar waves travelling into different directions
    class SuperimposedPlanarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      SuperimposedPlanarwave(const CellMaterialData& materialData, real phase = 0.0);

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const CellMaterialData& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const;
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
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const;
      private:
      std::array<double, 3> m_origin;
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
                        yateto::DenseTensorView<2,real,unsigned>& dofsQP) const;
      };
  }
}

#endif // PHYSICS_INITIALFIELD_H
