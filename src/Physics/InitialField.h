#ifndef PHYSICS_INITIALFIELD_H
#define PHYSICS_INITIALFIELD_H

#include <vector>
#include <array>
#include <complex>
#include "Initializer/typedefs.hpp"
#include <Kernels/precision.hpp>
#include <generated_code/init.h>
#include "Equations/datastructures.hpp"

namespace seissol {
  namespace physics {
    class InitialField {
    public:
      virtual ~InitialField() = default;
      virtual void evaluate(double time,
                            std::vector<std::array<double, 3>> const& points,
                            const seissol::model::Material_t& materialData,
                            yateto::DenseTensorView<2,real,unsigned>& dofsQP) const = 0;
    };

    class ZeroField : public InitialField {
    public:
      void evaluate(double,
                    std::vector<std::array<double, 3>> const&,
                    const seissol::model::Material_t& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override {
        dofsQP.setZero();
      }
    };

    //A planar wave travelling in direction kVec
    class Planarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
    Planarwave(const seissol::model::Material_t& materialData, 
               double phase,
               std::array<double, 3> kVec,
               std::vector<int> varField,
               std::vector<std::complex<double>> ampField
               );
    explicit Planarwave(const seissol::model::Material_t& materialData,
                        double phase = 0.0,
                        std::array<double, 3> kVec = {M_PI, M_PI, M_PI});

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const seissol::model::Material_t& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const override;
    protected:
      std::vector<int> m_varField;
      std::vector<std::complex<double>> m_ampField;
      const double m_phase;
      const std::array<double, 3> m_kVec;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES>  m_lambdaA;
      std::array<std::complex<double>, NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES> m_eigenvectors;
  private:
      void init(const seissol::model::Material_t& materialData);
    };

    //superimpose three planar waves travelling into different directions
    class SuperimposedPlanarwave : public InitialField {
    public:
      //! Choose phase in [0, 2*pi]
      SuperimposedPlanarwave(const seissol::model::Material_t& materialData, real phase = 0.0);

      void evaluate( double time,
                     std::vector<std::array<double, 3>> const& points,
                     const seissol::model::Material_t& materialData,
                     yateto::DenseTensorView<2,real,unsigned>& dofsQP ) const override;
    private:
      const std::array<std::array<double, 3>, 3>  m_kVec;
      const double                                m_phase;
      std::array<Planarwave, 3>                   m_pw;
    };

    //A part of a planar wave travelling in one direction
    class TravellingWave : public Planarwave{
    public:
      TravellingWave(const seissol::model::Material_t& materialData, const TravellingWaveParameters& travellingWaveParameters);

      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const seissol::model::Material_t& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
      private:
      std::array<double, 3> m_origin;
    };

    class ScholteWave : public InitialField {
    public:
      ScholteWave() = default;
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const seissol::model::Material_t& materialData,
                    yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
    };
    class SnellsLaw : public InitialField {
    public:
      SnellsLaw() = default;
      void evaluate(double time,
                    std::vector<std::array<double, 3>> const& points,
                    const seissol::model::Material_t& materialData,
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
                        const seissol::model::Material_t& materialData,
                        yateto::DenseTensorView<2,real,unsigned>& dofsQP) const override;
      };
  }
}

#endif // PHYSICS_INITIALFIELD_H
