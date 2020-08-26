//
// Created by adrian on 09.07.20.
//

#ifndef SEISSOL_DR_INITIALIZER_BASE_H
#define SEISSOL_DR_INITIALIZER_BASE_H

#include <c++/8.3.0/iostream>
#include <c++/8.3.0/unordered_map>
#include <Solver/Interoperability.h>
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace dr {
    namespace initializer {
      struct Base;
      struct FL_2;
      struct FL_3; //aging law
      struct FL_33;
      struct FL_103;
    }
  }
}


class seissol::dr::initializer::Base {
protected:
  static constexpr int numberOfPoints = tensor::QInterpolated::Shape[0];
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
  YAML::Node m_InputParam;
public:
  virtual ~Base() {}
  void setInputParam(const YAML::Node& Param) {m_InputParam = Param;}

  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
        initializers::LTSTree* dynRupTree,
        std::unordered_map<std::string,
        double*> faultParameters,
        unsigned* ltsFaceToMeshFace,
        seissol::Interoperability &e_interoperability) {
    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
      real  (*initialStressInFaultCS)[numOfPointsPadded][6] = it->var(dynRup->initialStressInFaultCS);  //get from fortran  EQN%InitialStressInFaultCS
      real  (*cohesion)[numOfPointsPadded]                  = it->var(dynRup->cohesion);                //get from faultParameters
      real  (*mu)[ numOfPointsPadded ]            = it->var(dynRup->mu);                                //get from fortran  EQN%IniMu(:,:)
      real  (*slip)[ numOfPointsPadded ]          = it->var(dynRup->slip);                              // = 0
      real  (*slip1)[numOfPointsPadded ]          = it->var(dynRup->slip1);                             // = 0
      real  (*slip2)[ numOfPointsPadded ]         = it->var(dynRup->slip2);                             // = 0
      real  (*slipRate1)[ numOfPointsPadded ]     = it->var(dynRup->slipRate1);                         //get from fortran  EQN%IniSlipRate1
      real  (*slipRate2)[numOfPointsPadded ]      = it->var(dynRup->slipRate2);                         //get from fortran  EQN%IniSlipRate2
      real  (*rupture_time)[ numOfPointsPadded ]  = it->var(dynRup->rupture_time);                      // = 0
      bool  (*RF)[ numOfPointsPadded ]            = it->var(dynRup->RF);                                //get from fortran
      real  (*peakSR)[ numOfPointsPadded ]        = it->var(dynRup->peakSR);                            // = 0
      real  (*tracXY)[ numOfPointsPadded ]        = it->var(dynRup->tracXY);                            // = 0
      real  (*tracXZ)[ numOfPointsPadded ]        = it->var(dynRup->tracXZ);                            // = 0

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
        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          if(faultParameters["cohesion"] != NULL ){
            cohesion[ltsFace][iBndGP] = static_cast<real>( faultParameters["cohesion"][meshFace * numberOfPoints] );
          }else{
            //TODO: maybe not log it = too much spam?
            //std::cout << "DR_initializer_base: cohesion set to 0, not found from faultParameters";
            cohesion[ltsFace][iBndGP] = 0;
          }
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          cohesion[ltsFace][iBndGP]               = 0.0;
        }
        e_interoperability.getDynRupParameters(ltsFace, meshFace, initialStressInFaultCS, mu, slipRate1, slipRate2, RF);

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
      real (*d_c)[numOfPointsPadded]                       = it->var(ConcreteLts->d_c);                 //from faultParameters
      real (*mu_S)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_S);                //from faultParameters
      real (*mu_D)[numOfPointsPadded]                      = it->var(ConcreteLts->mu_D);                //from faultParameters
      real (*forced_rupture_time)[numOfPointsPadded]       = it->var(ConcreteLts->forced_rupture_time); //from faultParameters
      bool *inst_healing                                   = it->var(ConcreteLts->inst_healing);        //from parameter file
      real *t_0                                            = it->var(ConcreteLts->t_0);                 //from parameter file
      bool *magnitude_out                                  = it->var(ConcreteLts->magnitude_out);       //from parameter file
      bool (*DS)[numOfPointsPadded]                        = it->var(ConcreteLts->DS);                  //from parameter file
      real *averaged_Slip                                  = it->var(ConcreteLts->averaged_Slip);       // = 0
      real (*dynStress_time)[numOfPointsPadded]            = it->var(ConcreteLts->dynStress_time);      // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];
        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
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
          d_c[ltsFace][iBndGP]                    = 0.0;
          mu_S[ltsFace][iBndGP]                   = 0.0;
          mu_D[ltsFace][iBndGP]                   = 0.0;
          forced_rupture_time[ltsFace][iBndGP]    = 0.0;
        }
        averaged_Slip[ltsFace]= 0.0;

        for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {    //loop includes padded elements
          dynStress_time[ltsFace][iBndGP] = 0.0;
          DS[ltsFace][iBndGP] = (m_InputParam["ds_output_on"]) ? true : false;
        }

        inst_healing[ltsFace] = (m_InputParam["inst_healing"]) ? true : false;
        t_0[ltsFace] = (m_InputParam["t_0"]) ? m_InputParam["t_0"].as<real>() : 0;
        magnitude_out[ltsFace] = (m_InputParam["magnitude_output_on"]) ? true : false;

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
        e_interoperability.getDynRupFL_3(ltsFace, meshFace, RS_f0, RS_a, RS_b, RS_sl0, RS_sr0, StateVar);

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

class seissol::dr::initializer::FL_103 : public seissol::dr::initializer::Base {
public:
  virtual void initializeFrictionMatrices(seissol::initializers::DynamicRupture *dynRup,
                                          initializers::LTSTree* dynRupTree,
                                          std::unordered_map<std::string,
                                              double*> faultParameters,
                                          unsigned* ltsFaceToMeshFace,
                                          seissol::Interoperability &e_interoperability) override {
    Base::initializeFrictionMatrices(dynRup, dynRupTree, faultParameters, ltsFaceToMeshFace, e_interoperability);
    seissol::initializers::DR_FL_103 *ConcreteLts = dynamic_cast<seissol::initializers::DR_FL_103 *>(dynRup);

    unsigned* layerLtsFaceToMeshFace = ltsFaceToMeshFace;

    for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {


      real  (*nucleationStressInFaultCS)[numOfPointsPadded][6]  = it->var(ConcreteLts->nucleationStressInFaultCS); //get from fortran
      bool *magnitude_out                                       = it->var(ConcreteLts->magnitude_out);      //par file
      real *t_0                                                 = it->var(ConcreteLts->t_0);                //par file
      real *RS_f0                                               = it->var(ConcreteLts->RS_f0);              //par file
      real *RS_b                                                = it->var(ConcreteLts->RS_b);               //par file
      real *RS_sr0                                              = it->var(ConcreteLts->RS_sr0);             //par file
      real *Mu_w                                                = it->var(ConcreteLts->Mu_w);               //par file
      real (*RS_sl0_array)[numOfPointsPadded]                   = it->var(ConcreteLts->RS_sl0_array);       //get from faultParameters
      real (*RS_a_array)[numOfPointsPadded]                     = it->var(ConcreteLts->RS_a_array);         //get from faultParameters
      real (*RS_srW_array)[numOfPointsPadded]                   = it->var(ConcreteLts->RS_srW_array);       //get from faultParameters
      bool (*DS)[numOfPointsPadded]                             = it->var(ConcreteLts->DS);                 //par file
      real *averaged_Slip                                       = it->var(ConcreteLts->averaged_Slip);      // = 0
      real (*stateVar)[numOfPointsPadded]                       = it->var(ConcreteLts->stateVar);           //get from Fortran = EQN%IniStateVar
      real (*dynStress_time)[numOfPointsPadded]                 = it->var(ConcreteLts->dynStress_time);     // = 0

      for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        unsigned meshFace = layerLtsFaceToMeshFace[ltsFace];


        e_interoperability.getDynRupFL_103(ltsFace, meshFace, nucleationStressInFaultCS, stateVar);

        t_0[ltsFace]      = m_InputParam["t_0"] ?     m_InputParam["t_0"].as<real>()    : 0;
        RS_f0[ltsFace]    = m_InputParam["rs_f0"] ?   m_InputParam["rs_f0"].as<real>()  : 0;
        RS_b[ltsFace]     = m_InputParam["rs_b"] ?    m_InputParam["rs_b"].as<real>()   : 0;
        RS_sr0[ltsFace]   = m_InputParam["rs_sr0"] ?  m_InputParam["rs_sr0"].as<real>() : 0;
        Mu_w[ltsFace]     = m_InputParam["mu_w"] ?    m_InputParam["mu_w"].as<real>()   : 0;

        for (unsigned iBndGP = 0; iBndGP < numOfPointsPadded; ++iBndGP) {    //loop includes padded elements
          dynStress_time[ltsFace][iBndGP] = 0.0;
          DS[ltsFace][iBndGP] = (m_InputParam["ds_output_on"]) ? true : false;
        }
        averaged_Slip[ltsFace]= 0.0;
        magnitude_out[ltsFace] = (m_InputParam["magnitude_output_on"]) ? true : false;

        for (unsigned iBndGP = 0; iBndGP < numberOfPoints; ++iBndGP) {
          RS_a_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["rs_a"][meshFace * numberOfPoints] );
          RS_srW_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["rs_srW"][meshFace * numberOfPoints] );
          RS_sl0_array[ltsFace][iBndGP] = static_cast<real>( faultParameters["RS_sl0"][meshFace * numberOfPoints] );
        }
        //initialize padded elements for vectorization
        for (unsigned iBndGP = numberOfPoints; iBndGP < numOfPointsPadded; ++iBndGP) {
          RS_a_array[ltsFace][iBndGP] = 0.0;
          RS_srW_array[ltsFace][iBndGP] = 0.0;
          RS_sl0_array[ltsFace][iBndGP] = 0.0;
        }

      }//lts-face loop
      layerLtsFaceToMeshFace += it->getNumberOfCells();
    }//leaf_iterator loop

  }
};

#endif //SEISSOL_DR_INITIALIZER_BASE_H
