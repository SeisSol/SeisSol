//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_INITIALIZER_BASE_H
#define SEISSOL_DR_INITIALIZER_BASE_H

#include <c++/8.3.0/iostream>
#include <c++/8.3.0/unordered_map>
#include <Solver/Interoperability.h>

namespace seissol {
  namespace dr {
    namespace initializer {
      struct Base;
      struct FL_2;
      struct FL_3; //aging law
      struct FL_33;
    }
  }
}


class seissol::dr::initializer::Base {
protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
public:
  virtual ~Base() {}

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
        initializers::LTSTree* dynRupTree,
        std::unordered_map<std::string,
        double*> faultParameters,
        unsigned* ltsFaceToMeshFace,
        seissol::Interoperability &e_interoperability) {
    //TODO: add base initialization of DynamicRupture here
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real  (*mu)[ numOfPointsPadded ]            = it->var(dynRup->mu);  //fortran
      real  (*slip)[ numOfPointsPadded ]          = it->var(dynRup->slip);
      real  (*slip1)[numOfPointsPadded ]          = it->var(dynRup->slip1);
      real  (*slip2)[ numOfPointsPadded ]         = it->var(dynRup->slip2);
      real  (*slipRate1)[ numOfPointsPadded ]     = it->var(dynRup->slipRate1);    //fortran
      real  (*slipRate2)[numOfPointsPadded ]      = it->var(dynRup->slipRate2);    //fortran
      real  (*rupture_time)[ numOfPointsPadded ]  = it->var(dynRup->rupture_time);
      bool  (*RF)[ numOfPointsPadded ]            = it->var(dynRup->RF);                  //fortran
      real (*peakSR)[ numOfPointsPadded ]         = it->var(dynRup->peakSR);
      real  (*tracXY)[ numOfPointsPadded ]        = it->var(dynRup->tracXY);
      real  (*tracXZ)[ numOfPointsPadded ]        = it->var(dynRup->tracXZ);

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
        for (unsigned iBndGP = 0; iBndGP < init::QInterpolated::Stop[0]; ++iBndGP) {    //loop includes padded elements
          slip[ltsFace][iBndGP] = 0.0;
          slip1[ltsFace][iBndGP] = 0.0;
          slip2[ltsFace][iBndGP] = 0.0;
          rupture_time[ltsFace][iBndGP] = 0.0;
          peakSR[ltsFace][iBndGP] = 0.0;
          tracXY[ltsFace][iBndGP] = 0.0;
          tracXZ[ltsFace][iBndGP] = 0.0;
        }
        //get initial values from fortran
        e_interoperability.getDynRupParameters(ltsFace, meshFace, mu, slipRate1, slipRate2, RF);
      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
    std::cout << "init DR for Base\n";
  }
};
class seissol::dr::initializer::FL_2 : public seissol::dr::initializer::Base {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
        initializers::LTSTree* dynRupTree,
        std::unordered_map<std::string,
        double*> faultParameters,
        unsigned* ltsFaceToMeshFace,
        seissol::Interoperability &e_interoperability) override {
    Base::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_2 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_2 *>(dynRup);
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;


    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real (*initialStressInFaultCS)[numOfPointsPadded][6] = it->var(ConcreteLts->initialStressInFaultCS);
      real (*cohesion)[numOfPointsPadded]                  = it->var(ConcreteLts->cohesion);
      real (*d_c)[numOfPointsPadded]                       = it->var(ConcreteLts->d_c);
      real (*mu_S)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_S);
      real (*mu_D)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_D);
      real (*forced_rupture_time)[numOfPointsPadded]       = it->var(ConcreteLts->forced_rupture_time);
      bool *inst_healing                                   = it->var(ConcreteLts->inst_healing);       //fortran
      real *t_0                                            = it->var(ConcreteLts->t_0);                //fortran
      bool *magnitude_out                                  = it->var(ConcreteLts->magnitude_out);      //fortran
      bool (*DS)[numOfPointsPadded]                        = it->var(ConcreteLts->DS);                 //fortran


      real *averaged_Slip                                  = it->var(ConcreteLts->averaged_Slip);
      real (*dynStress_time)[numOfPointsPadded]            = it->var(ConcreteLts->dynStress_time);

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          cohesion[ltsFace][iBndGP]               = static_cast<real>( faultParameters["cohesion"][meshFace * numberOfPoints] );
          d_c[ltsFace][iBndGP]                    = static_cast<real>( faultParameters["d_c"][meshFace * numberOfPoints] );
          mu_S[ltsFace][iBndGP]                   = static_cast<real>( faultParameters["mu_s"][meshFace * numberOfPoints] );
          mu_D[ltsFace][iBndGP]                   = static_cast<real>( faultParameters["mu_d"][meshFace * numberOfPoints] );
          if(faultParameters["forced_rupture_time"] != NULL ){
              forced_rupture_time[ltsFace][iBndGP]    = static_cast<real>( faultParameters["forced_rupture_time"][meshFace * numberOfPoints] );
          }else{
              forced_rupture_time[ltsFace][iBndGP]    = 0.0;
          }
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          cohesion[ltsFace][iBndGP]               = 0.0;
          d_c[ltsFace][iBndGP]                    = 0.0;
          mu_S[ltsFace][iBndGP]                   = 0.0;
          mu_D[ltsFace][iBndGP]                   = 0.0;
          forced_rupture_time[ltsFace][iBndGP]    = 0.0;
        }
        averaged_Slip[ltsFace]= 0.0;

        //get initial values from fortran
        //TODO: get intial initialStressInFaultCS;
        e_interoperability.getDynRupFL_2(ltsFace, meshFace, initialStressInFaultCS, t_0, magnitude_out, DS, inst_healing);

        for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {    //loop includes padded elements
          dynStress_time[ltsFace][iBndGP] = 0.0;
        }
      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
  }
};


class seissol::dr::initializer::FL_3 : public seissol::dr::initializer::Base {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          initializers::LTSTree* dynRupTree,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override {
    Base::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_3 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_3 *>(dynRup);
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real *RS_f0                                               = it->var(ConcreteLts->RS_f0);
      real *RS_a                                                = it->var(ConcreteLts->RS_a);
      real *RS_b                                                = it->var(ConcreteLts->RS_b);
      real *RS_sl0                                              = it->var(ConcreteLts->RS_sl0);
      real *RS_sr0                                              = it->var(ConcreteLts->RS_sr0);
      real (*StateVar)[numOfPointsPadded]                       = it->var(ConcreteLts->StateVar);

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];

        //get initial values from fortran
        //TODO: write this function_
        //e_interoperability.getDynRupFL_3(ltsFace, meshFace, RS_f0, RS_a, RS_b, RS_sl0, RS_sr0,StateVar);

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop
  }
};



class seissol::dr::initializer::FL_33 : public seissol::dr::initializer::Base {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
          initializers::LTSTree* dynRupTree,
          std::unordered_map<std::string,
          double*> faultParameters,
          unsigned* ltsFaceToMeshFace,
          seissol::Interoperability &e_interoperability) override {
    Base::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_33 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_33 *>(dynRup);

  }
};

#endif //SEISSOL_DR_INITIALIZER_BASE_H
