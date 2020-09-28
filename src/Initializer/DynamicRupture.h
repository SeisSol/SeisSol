/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/
 
#ifndef INITIALIZER_DR_H_
#define INITIALIZER_DR_H_

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>

namespace seissol {
  namespace initializers {
    struct DynamicRupture;
    struct DR_FL_2;
    struct DR_FL_3;
    struct DR_FL_6;
    struct DR_FL_33;
    struct DR_FL_103;
    struct DR_FL_103_Thermal;
  }
}

//TODO: remove space independent parameters:  rs_f0, rs_b, rs_sr0



struct seissol::initializers::DynamicRupture {
  //assert(init::QInterpolated::Start[0] == 0); ?
protected:
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
public:
  virtual ~DynamicRupture() {}
  bool IsFaultParameterizedByTraction;        //true if Traction T_n , T_s, T_d is stored in iniBulk_XX, iniShearXY, iniShearXZ
  Variable<real*>                                                   timeDerivativePlus;
  Variable<real*>                                                   timeDerivativeMinus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStatePlus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStateMinus;
  Variable<DRGodunovData>                                           godunovData;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverPlus;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverMinus;
  Variable<DRFaceInformation>                                       faceInformation;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsPlus;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsMinus;

  Variable<ImpedancesAndEta>                            impAndEta;
  //size padded for vectorization

  Variable<real[ numOfPointsPadded ]>                   iniBulkXX;
  Variable<real[ numOfPointsPadded ]>                   iniBulkYY;
  Variable<real[ numOfPointsPadded ]>                   iniBulkZZ;
  Variable<real[ numOfPointsPadded ]>                   iniShearXY;
  Variable<real[ numOfPointsPadded ]>                   iniShearYZ;
  Variable<real[ numOfPointsPadded ]>                   iniShearXZ;

  Variable<real[numOfPointsPadded][6]>                  initialStressInFaultCS;
  Variable<real[ numOfPointsPadded ]>                   cohesion;
  Variable<real[ numOfPointsPadded ]>                   mu;
  Variable<real[ numOfPointsPadded ]>                   slip;
  Variable<real[ numOfPointsPadded ]>                   slip1;
  Variable<real[ numOfPointsPadded ]>                   slip2;
  Variable<real[ numOfPointsPadded ]>                   slipRate1;
  Variable<real[ numOfPointsPadded ]>                   slipRate2;
  Variable<real[ numOfPointsPadded ]>                   rupture_time;
  Variable<bool[ numOfPointsPadded ]>                   RF;
  Variable<real[ numOfPointsPadded ]>                   peakSR;
  Variable<real[ numOfPointsPadded ]>                   tracXY;
  Variable<real[ numOfPointsPadded ]>                   tracXZ;

  virtual void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      timeDerivativePlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(     timeDerivativeMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(        imposedStatePlus,             mask,     PAGESIZE_HEAP,      seissol::memory::Standard );
    tree.addVar(       imposedStateMinus,             mask,     PAGESIZE_HEAP,      seissol::memory::Standard );
    tree.addVar(             godunovData,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          fluxSolverPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         fluxSolverMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         faceInformation,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          waveSpeedsPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         waveSpeedsMinus,             mask,                 1,      seissol::memory::Standard );

    tree.addVar(      impAndEta,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniBulkXX,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniBulkYY,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniBulkZZ,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniShearXY,                     mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniShearYZ,                     mask,                 1,      seissol::memory::Standard );
    tree.addVar(      iniShearXZ,                     mask,                 1,      seissol::memory::Standard );
    tree.addVar(      initialStressInFaultCS,         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      cohesion,                       mask,                 1,      seissol::memory::Standard );
    tree.addVar(      rupture_time,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RF,                             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      mu,                             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slip,                           mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slip1,                          mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slip2,                          mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slipRate1,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slipRate2,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      peakSR,                         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      tracXY,                         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      tracXZ,                         mask,                 1,      seissol::memory::Standard );
  }
};


struct seissol::initializers::DR_FL_2 : public seissol::initializers::DynamicRupture {
    Variable<real[ numOfPointsPadded ]>                   d_c;
    Variable<real[ numOfPointsPadded ]>                   mu_S;
    Variable<real[ numOfPointsPadded ]>                   mu_D;
    Variable<real[ numOfPointsPadded ]>                   forced_rupture_time;
    Variable<bool[ numOfPointsPadded ]>                   DS;
    Variable<real>                                        averaged_Slip;
    Variable<real[ numOfPointsPadded ]>                   dynStress_time;

    virtual void addTo(initializers::LTSTree& tree) {
        seissol::initializers::DynamicRupture::addTo(tree);
        LayerMask mask = LayerMask(Ghost);
        tree.addVar(      d_c,                              mask,                 1,      seissol::memory::Standard );
        tree.addVar(      mu_S,                             mask,                 1,      seissol::memory::Standard );
        tree.addVar(      mu_D,                             mask,                 1,      seissol::memory::Standard );
        tree.addVar(      forced_rupture_time,              mask,                 1,      seissol::memory::Standard );
        tree.addVar(      DS,                               mask,                 1,      seissol::memory::Standard );
        tree.addVar(      averaged_Slip,                    mask,                 1,      seissol::memory::Standard );
        tree.addVar(      dynStress_time,                   mask,                 1,      seissol::memory::Standard );
    }
};

struct seissol::initializers::DR_FL_3 : public seissol::initializers::DynamicRupture {
  Variable<real>                                                  RS_f0;                      //face independent
  Variable<real>                                                  RS_a;                       //face independent
  Variable<real>                                                  RS_b;                       //face independent
  Variable<real>                                                  RS_sl0;                     //face independent
  Variable<real>                                                  RS_sr0;                     //face independent
  Variable<real[ numOfPointsPadded ]>                             StateVar;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      RS_f0,            mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_a,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_b,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_sl0,           mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_sr0,           mask,                 1,      seissol::memory::Standard );
    tree.addVar(      StateVar,         mask,                 1,      seissol::memory::Standard );
  }
};


struct seissol::initializers::DR_FL_33 : public seissol::initializers::DynamicRupture {
  Variable<real[numOfPointsPadded][6]>                            nucleationStressInFaultCS;
  Variable<real>                                                  averaged_Slip;

    virtual void addTo(initializers::LTSTree& tree) {
      seissol::initializers::DynamicRupture::addTo(tree);
      LayerMask mask = LayerMask(Ghost);
      tree.addVar(      nucleationStressInFaultCS,  mask,                 1,      seissol::memory::Standard );
      tree.addVar(      averaged_Slip,              mask,                 1,      seissol::memory::Standard );
    }
};

struct seissol::initializers::DR_FL_103 : public seissol::initializers::DynamicRupture {
  Variable<real[numOfPointsPadded][6]>                            nucleationStressInFaultCS;
  Variable<real[ numOfPointsPadded ]>                             RS_sl0_array;
  Variable<real[ numOfPointsPadded ]>                             RS_a_array;
  Variable<real[ numOfPointsPadded ]>                             RS_srW_array;
  Variable<bool[ numOfPointsPadded ]>                             DS;
  Variable<real>                                                  averaged_Slip;
  Variable<real[ numOfPointsPadded ]>                             stateVar;
  Variable<real[ numOfPointsPadded ]>                             dynStress_time;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      nucleationStressInFaultCS,  mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_sl0_array,               mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_a_array,                 mask,                 1,      seissol::memory::Standard );
    tree.addVar(      RS_srW_array,               mask,                 1,      seissol::memory::Standard );
    tree.addVar(      DS,                         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      stateVar,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      averaged_Slip,              mask,                 1,      seissol::memory::Standard );
    tree.addVar(      dynStress_time,             mask,                 1,      seissol::memory::Standard );
  }
};

struct seissol::initializers::DR_FL_103_Thermal : public seissol::initializers::DR_FL_103 {

  static constexpr unsigned int TP_grid_nz = 60;  //todo: make this global?
  Variable<real[numOfPointsPadded][2]>                            TP;
  Variable<real[numOfPointsPadded][TP_grid_nz]>                   TP_Theta;
  Variable<real[numOfPointsPadded][TP_grid_nz]>                   TP_sigma;
  Variable<real[numOfPointsPadded]>                               TP_half_width_shear_zone;
  Variable<real[numOfPointsPadded]>                               alpha_hy;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DR_FL_103::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      TP,                         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_Theta,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_sigma,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_half_width_shear_zone,   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      alpha_hy,                   mask,                 1,      seissol::memory::Standard );
  }
};


struct seissol::initializers::DR_FL_6 : public seissol::initializers::DR_FL_2 {

  Variable<real[numOfPointsPadded]>                               strengthData;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DR_FL_2::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      strengthData,                   mask,                 1,      seissol::memory::Standard );
  }
};

#endif
