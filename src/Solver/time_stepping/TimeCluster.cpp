/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2017, SeisSol Group
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
 * LTS cluster in SeisSol.
 **/

#include "Parallel/MPI.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SeisSol.h"
#include "TimeCluster.h"
#include <Solver/Interoperability.h>
#include <SourceTerm/PointSource.h>
#include <Kernels/TimeCommon.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Receiver.h>
#include <Monitoring/FlopCounter.hpp>


#include <cassert>
#include <cstring>


#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
extern volatile unsigned int* volatile g_handleRecvs;
extern volatile unsigned int* volatile g_handleSends;
#endif

//! fortran interoperability
extern seissol::Interoperability e_interoperability;


seissol::time_stepping::TimeCluster::TimeCluster( unsigned int                   i_clusterId,
                                                  unsigned int                   i_globalClusterId,
                                                  struct MeshStructure          *i_meshStructure,
                                                  struct GlobalData             *i_globalData,
                                                  seissol::initializers::TimeCluster* i_clusterData,
                                                  seissol::initializers::TimeCluster* i_dynRupClusterData,
                                                  seissol::initializers::LTS*         i_lts,
                                                  seissol::initializers::DynamicRupture* i_dynRup,
                                                  LoopStatistics*                        i_loopStatistics ):
 // cluster ids
 m_clusterId(               i_clusterId                ),
 m_globalClusterId(         i_globalClusterId          ),
 // mesh structure
 m_meshStructure(           i_meshStructure            ),
 // global data
 m_globalData(              i_globalData               ),
 m_clusterData(             i_clusterData              ),
 m_dynRupClusterData(       i_dynRupClusterData        ),
 m_lts(                     i_lts                      ),
 m_dynRup(                  i_dynRup                   ),
 // cells
 m_cellToPointSources(      NULL                       ),
 m_numberOfCellToPointSourcesMappings(0                ),
 m_pointSources(            NULL                       ),

 m_loopStatistics(          i_loopStatistics           ),
 m_receiverCluster(          nullptr                   ),
 //Code added by Adrian:
 m_friction_data(tensor::QInterpolated::Shape[0], e_interoperability.getnSide() )
{
    // assert all pointers are valid
    assert( m_meshStructure                            != NULL );
    assert( m_globalData                               != NULL );
    assert( m_clusterData                              != NULL );

  // default: no updates are allowed
  m_updatable.localCopy           = false;
  m_updatable.neighboringCopy     = false;
  m_updatable.localInterior       = false;
  m_updatable.neighboringInterior = false;
#ifdef USE_MPI
  m_sendLtsBuffers                = false;
#endif
  m_resetLtsBuffers               = false;
  // set timings to zero
  m_numberOfTimeSteps             = 0;
  m_receiverTime                  = 0;
  m_timeStepWidth                 = 0;
  m_subTimeStart                  = 0;
  m_numberOfFullUpdates           = 0;
  m_fullUpdateTime                = 0;
  m_predictionTime                = 0;


//Code added by adrian
//e_interoperability.getFrictionData(tensor::QInterpolated::Shape[0], m_friction_data);



  m_dynamicRuptureFaces = (i_dynRupClusterData->child<Ghost>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Copy>().getNumberOfCells() > 0)
	|| (i_dynRupClusterData->child<Interior>().getNumberOfCells() > 0);
  
  m_timeKernel.setGlobalData(m_globalData);
  m_localKernel.setGlobalData(m_globalData);
  m_localKernel.setInitConds(&e_interoperability.getInitialConditions());
  m_neighborKernel.setGlobalData(m_globalData);
  m_dynamicRuptureKernel.setGlobalData(m_globalData);

  computeFlops();

  m_regionComputeLocalIntegration = m_loopStatistics->getRegion("computeLocalIntegration");
  m_regionComputeNeighboringIntegration = m_loopStatistics->getRegion("computeNeighboringIntegration");
  m_regionComputeDynamicRupture = m_loopStatistics->getRegion("computeDynamicRupture");
}

seissol::time_stepping::TimeCluster::~TimeCluster() {
#ifndef NDEBUG
  logInfo() << "#(time steps):" << m_numberOfTimeSteps;
#endif
}

void seissol::time_stepping::TimeCluster::setPointSources( sourceterm::CellToPointSourcesMapping const* i_cellToPointSources,
                                                           unsigned i_numberOfCellToPointSourcesMappings,
                                                           sourceterm::PointSources const* i_pointSources )
{
  m_cellToPointSources = i_cellToPointSources;
  m_numberOfCellToPointSourcesMappings = i_numberOfCellToPointSourcesMappings;
  m_pointSources = i_pointSources;
}

void seissol::time_stepping::TimeCluster::writeReceivers() {
  SCOREP_USER_REGION( "writeReceivers", SCOREP_USER_REGION_TYPE_FUNCTION )

  if (m_receiverCluster != nullptr) {
    m_receiverTime = m_receiverCluster->calcReceivers(m_receiverTime, m_fullUpdateTime, m_timeStepWidth);
  }
}

void seissol::time_stepping::TimeCluster::computeSources() {
  SCOREP_USER_REGION( "computeSources", SCOREP_USER_REGION_TYPE_FUNCTION )

  // Return when point sources not initialised. This might happen if there
  // are no point sources on this rank.
  if (m_numberOfCellToPointSourcesMappings != 0) {
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
    for (unsigned mapping = 0; mapping < m_numberOfCellToPointSourcesMappings; ++mapping) {
      unsigned startSource = m_cellToPointSources[mapping].pointSourcesOffset;
      unsigned endSource = m_cellToPointSources[mapping].pointSourcesOffset + m_cellToPointSources[mapping].numberOfPointSources;
      if (m_pointSources->mode == sourceterm::PointSources::NRF) {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceNRF( m_pointSources->mInvJInvPhisAtSources[source],
                                                       m_pointSources->tensor[source],
                                                       m_pointSources->A[source],
                                                       m_pointSources->stiffnessTensor[source],
                                                       m_pointSources->slipRates[source],
                                                       m_fullUpdateTime,
                                                       m_fullUpdateTime + m_timeStepWidth,
                                                       *m_cellToPointSources[mapping].dofs );
        }
      } else {
        for (unsigned source = startSource; source < endSource; ++source) {
          sourceterm::addTimeIntegratedPointSourceFSRM( m_pointSources->mInvJInvPhisAtSources[source],
                                                        m_pointSources->tensor[source],
                                                        m_pointSources->slipRates[source][0],
                                                        m_fullUpdateTime,
                                                        m_fullUpdateTime + m_timeStepWidth,
                                                        *m_cellToPointSources[mapping].dofs );
        }
      }
    }
  }
}

void seissol::time_stepping::TimeCluster::computeDynamicRupture( seissol::initializers::Layer&  layerData ) {
  SCOREP_USER_REGION( "computeDynamicRupture", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeDynamicRupture);

  DRFaceInformation*                    faceInformation                                                   = layerData.var(m_dynRup->faceInformation);
  DRGodunovData*                        godunovData                                                       = layerData.var(m_dynRup->godunovData);
  real**                                timeDerivativePlus                                                = layerData.var(m_dynRup->timeDerivativePlus);
  real**                                timeDerivativeMinus                                               = layerData.var(m_dynRup->timeDerivativeMinus);
  real                                (*imposedStatePlus)[tensor::QInterpolated::size()]                  = layerData.var(m_dynRup->imposedStatePlus);
  real                                (*imposedStateMinus)[tensor::QInterpolated::size()]                 = layerData.var(m_dynRup->imposedStateMinus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                                    = layerData.var(m_dynRup->waveSpeedsPlus);
  seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                                   = layerData.var(m_dynRup->waveSpeedsMinus);

  alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
  alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][tensor::QInterpolated::size()];

  //Code added by ADRIAN

  //TODO: outsource this to initialisation:
  /*
  int nSide_Fortran = 0;
  e_interoperability.getnSide(nSide_Fortran);
  const size_t nSide = nSide_Fortran;
    struct seissol::physics::FrictionData friction_data(numberOfPoints, nSide);
    e_interoperability.getFrictionData(numberOfPoints, friction_data);
*/
  //std::cout << "before initializing fric data" << std::endl;

    const size_t numberOfPoints = tensor::QInterpolated::Shape[0];
    if(m_friction_data.initialized == false){
        e_interoperability.getFrictionData(m_friction_data);
        m_friction_data.initialized = true;
    }
    m_friction_data.function_call++;


    //double resampleMatrix[]  = {};
  //resampleMatrix =    init::resample::Values;
  //yateto::DenseTensorView<2, double> view = init::resample::view::create(const_cast<double *>(init::resample::Values));


/*
    const size_t nSide = e_interoperability.getnSide();

    std::cout << "l_slip before:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (iSide: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            std::cout << m_friction_data.getSlip(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "l_slip1 before:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (iSide: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            std::cout << m_friction_data.getSlip1(i,j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "l_slip2 before:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (iSide: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            std::cout << m_friction_data.getSlip2(i,j) << " ";
        }
        std::cout << std::endl;
    }


    std::cout << "l_rupture_time before:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (iSide: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            std::cout << m_friction_data.getRupture_time(i,j) << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "l_averaged_Slip before:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (iSide: " << j << ") ";
        std::cout << m_friction_data.averaged_Slip[j] << " ";
    }
    std::cout << std::endl;
/*

    for(int k = 0; k < nSide ; k++){
        std::cout << " (k: " << k << ") ";
        for(int j = 0; j < 6; j++){
            for (int i = 0; i < numberOfPoints; i++){
                std::cout << friction_data.initialStressInFaultCS[i + j* numberOfPoints+ k*6*numberOfPoints] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }



    for(int j = 0; j < nSide; j++){
        //std::cout << " (j: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            //std::cout << friction_data.mu_S[j][i] << " ";
        }
        //std::cout << std::endl;
    }
    std::cout << "cohesion(21,43): " << friction_data.getCohesion(21,43) << std::endl;
    std::cout << "initialStressInFaultCS(10,1,32): " << friction_data.getInitialStressInFaultCS(10,1,32) << std::endl;

    std::cout << "cohesion:" << std::endl;
    for(int j = 0; j < nSide; j++){
        std::cout << " (j: " << j << ") ";
        for (int i = 0; i < numberOfPoints; i++){
            std::cout << friction_data.cohesion[i + j* numberOfPoints] << " ";
        }
        std::cout << std::endl;
    }


   //*/

    struct seissol::physics::FrictionData fortran_data_before(m_friction_data.numberOfPoints, m_friction_data.nFace);
    e_interoperability.getFrictionData(fortran_data_before);
    m_friction_data.isEqualToFortran(fortran_data_before);

//TODO: change back to parallel
#ifdef _OPENMP
  #pragma omp parallel for schedule(static) private(QInterpolatedPlus,QInterpolatedMinus)
#endif
  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    unsigned prefetchFace = (face < layerData.getNumberOfCells()-1) ? face+1 : face;
    m_dynamicRuptureKernel.spaceTimeInterpolation(  faceInformation[face],
                                                    m_globalData,
                                                   &godunovData[face],
                                                    timeDerivativePlus[face],
                                                    timeDerivativeMinus[face],
                                                    QInterpolatedPlus,
                                                    QInterpolatedMinus,
                                                    timeDerivativePlus[prefetchFace],
                                                    timeDerivativeMinus[prefetchFace] );




    //Code added by ADRIAN
    //insert c++ evaluate_friction_law here:

/*
    // legacy code:
    e_interoperability.evaluateFrictionLaw( static_cast<int>(faceInformation[face].meshFace),
                                            QInterpolatedPlus,
                                            QInterpolatedMinus,
                                            imposedStatePlus[face],
                                            imposedStateMinus[face],
                                            m_fullUpdateTime,
                                            m_dynamicRuptureKernel.timePoints,
                                            m_dynamicRuptureKernel.timeWeights,
                                            waveSpeedsPlus[face],
                                            waveSpeedsMinus[face] );

      std::cout << "Fortran imposedStatePlusView: " << std::endl;
      auto imposedStatePlusView3 = init::QInterpolated::view::create(imposedStatePlus[face]);
      for(int j = 0; j < 9; j++){
          std::cout << " (j: " << j << ") " <<  "(face: " << face << ") ";
          for (int i = 0; i < numberOfPoints; i++) {
              std::cout << imposedStatePlusView3(i,j) << " ";
          }
          std::cout << std::endl;
      }
      //*/


    double TractionGP_XY[numberOfPoints][CONVERGENCE_ORDER] = {{}};
    double TractionGP_XZ[numberOfPoints][CONVERGENCE_ORDER] = {{}};
    double NorStressGP[numberOfPoints][CONVERGENCE_ORDER] = {{}};
    double XYStressGP[numberOfPoints][CONVERGENCE_ORDER]= {{}};
    double XZStressGP[numberOfPoints][CONVERGENCE_ORDER]= {{}};

    double *TractionGP_XY2[numberOfPoints];
    double *TractionGP_XZ2[numberOfPoints];
    double *NorStressGP2[numberOfPoints];
    double *XYStressGP2[numberOfPoints];
    double *XZStressGP2[numberOfPoints];
    for (size_t i = 0; i < numberOfPoints; ++i){
      TractionGP_XY2[i] = TractionGP_XY[i];
      TractionGP_XZ2[i] = TractionGP_XZ[i];
      NorStressGP2[i] = NorStressGP[i];
      XYStressGP2[i] = XYStressGP[i];
      XZStressGP2[i] = XZStressGP[i];
    }
    int iFace;
    double Zp_inv, Zp_neig_inv, Zs_inv, Zs_neig_inv, eta_p, eta_s;
    double time = m_fullUpdateTime;
    double *timePoints = &m_dynamicRuptureKernel.timePoints[0];
    double rho = waveSpeedsPlus->density;
    double rho_neig= waveSpeedsMinus->density;
    double w_speed[3] = {waveSpeedsPlus->pWaveVelocity, waveSpeedsPlus->sWaveVelocity, waveSpeedsPlus->sWaveVelocity};
    double w_speed_neig[3] = {waveSpeedsMinus->pWaveVelocity, waveSpeedsMinus->sWaveVelocity, waveSpeedsMinus->sWaveVelocity};


    iFace = static_cast<int>(faceInformation[face].meshFace);
    //int godunovLd = init::QInterpolated::Stop[0] - init::QInterpolated::Start[0];
    //int iSide = friction_data.side[iFace];     //= l_domain%MESH%Fault%Face(i_face,2,1)          ! iElem denotes "+" side
    //int iElem = friction_data.elem[iFace];     //= l_domain%MESH%Fault%Face(i_face,1,1)
    Zp_inv = 1.0 / (rho * w_speed[0]);       //(rho * waveSpeedsPlus->pWaveVelocity);
    Zp_neig_inv = 1.0 / (rho_neig *w_speed_neig[0]);     //(rho_neig * waveSpeedsMinus->pWaveVelocity);
    Zs_inv = 1.0 / (rho * w_speed[1]);               //(rho * waveSpeedsPlus->sWaveVelocity);
    Zs_neig_inv = 1.0 / (rho_neig *w_speed_neig[1]); //(rho_neig * waveSpeedsMinus->sWaveVelocity);
    eta_p = 1.0 / (Zp_inv + Zp_neig_inv);
    eta_s = 1.0 / (Zs_inv + Zs_neig_inv);

    //int numberOfPointest = tensor::QInterpolated::Shape[0];
    //int BasisfuncStart = init::QInterpolated::Start[0];
    //int BasisfuncStop = init::QInterpolated::Stop[0];


    for(int j = 0; j < CONVERGENCE_ORDER; j++){
        auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
        auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      for(int i = 0; i < numberOfPoints; i++){
          NorStressGP[i][j] = eta_p * (QInterpolatedMinusView(i,6) - QInterpolatedPlusView(i,6) + QInterpolatedPlusView(i,0) * Zp_inv + QInterpolatedMinusView(i,0) * Zp_neig_inv);
          XYStressGP[i][j]  = eta_s * (QInterpolatedMinusView(i,7) - QInterpolatedPlusView(i,7) + QInterpolatedPlusView(i,3) * Zs_inv + QInterpolatedMinusView(i,3) * Zs_neig_inv);
          XZStressGP[i][j] = eta_s * (QInterpolatedMinusView(i,8) - QInterpolatedPlusView(i,8) + QInterpolatedPlusView(i,5) * Zs_inv + QInterpolatedMinusView(i,5) * Zs_neig_inv);
      }
    }
    static_assert(tensor::QInterpolated::Shape[0] == tensor::resample::Shape[0], "Different number of quadrature points?");

    /*
    e_interoperability.evaluateFrictionLaw( static_cast<int>(faceInformation[face].meshFace),
                                          QInterpolatedPlus,
                                          QInterpolatedMinus,
                                          imposedStatePlus[face],
                                          imposedStateMinus[face],
                                          m_fullUpdateTime,
                                          m_dynamicRuptureKernel.timePoints,
                                          m_dynamicRuptureKernel.timeWeights,
                                          waveSpeedsPlus[face],
                                          waveSpeedsMinus[face] );
    //*/


    //std::cout << "(face: " << face << ") ";
    seissol::physics::Evaluate_friction_law evaluateFriction;
    evaluateFriction.Eval_friction_law(TractionGP_XY2, TractionGP_XZ2, NorStressGP2, XYStressGP2, XZStressGP2,
                                       iFace, m_friction_data.side[iFace],
                                       m_friction_data.elem[iFace], time, timePoints, rho, rho_neig, w_speed, w_speed_neig,
                                       const_cast<double *>(init::resample::Values),
                                       m_friction_data );


    auto imposedStatePlusView = init::QInterpolated::view::create(imposedStatePlus[face]);
    auto imposedStateMinusView = init::QInterpolated::view::create(imposedStateMinus[face]);
    //initialize to 0
    imposedStateMinusView.setZero();
    imposedStatePlusView.setZero();

    for(int j = 0; j < CONVERGENCE_ORDER; j++) {
        auto QInterpolatedPlusView = init::QInterpolated::view::create(QInterpolatedPlus[j]);
        auto QInterpolatedMinusView = init::QInterpolated::view::create(QInterpolatedMinus[j]);
      for (int i = 0; i < numberOfPoints; i++) {
          imposedStateMinusView(i,0) += m_dynamicRuptureKernel.timeWeights[j] * NorStressGP[i][j];
          imposedStateMinusView(i,3) += m_dynamicRuptureKernel.timeWeights[j]  * TractionGP_XY[i][j];
          imposedStateMinusView(i,5) += m_dynamicRuptureKernel.timeWeights[j]  * TractionGP_XZ[i][j];
          imposedStateMinusView(i,6) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedMinusView(i,6) - Zp_neig_inv * (NorStressGP[i][j]-QInterpolatedMinusView(i,0)));
          imposedStateMinusView(i,7) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedMinusView(i,7) - Zs_neig_inv * (TractionGP_XY[i][j]-QInterpolatedMinusView(i,3)));
          imposedStateMinusView(i,8) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedMinusView(i,8) - Zs_neig_inv * (TractionGP_XZ[i][j]-QInterpolatedMinusView(i,5)));

          imposedStatePlusView(i,0) += m_dynamicRuptureKernel.timeWeights[j]  * NorStressGP[i][j];
          imposedStatePlusView(i,3) += m_dynamicRuptureKernel.timeWeights[j]  * TractionGP_XY[i][j];
          imposedStatePlusView(i,5) += m_dynamicRuptureKernel.timeWeights[j]  * TractionGP_XZ[i][j];
          imposedStatePlusView(i,6) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedPlusView(i,6) + Zp_inv * (NorStressGP[i][j]-QInterpolatedPlusView(i,0)));
          imposedStatePlusView(i,7) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedPlusView(i,7) + Zs_inv * (TractionGP_XY[i][j]-QInterpolatedPlusView(i,3)));
          imposedStatePlusView(i,8) += m_dynamicRuptureKernel.timeWeights[j]  * (QInterpolatedPlusView(i,8) + Zs_inv * (TractionGP_XZ[i][j]-QInterpolatedPlusView(i,5)));
      } //End numberOfPoints-loop
    } //End CONVERGENCE_ORDER-loop

    e_interoperability.setFrictionOutput( m_friction_data, iFace);
      //*/




    //std::cout << " finished: " <<  "(face: " << face << ") "  << std::endl;
/*
      std::cout << "C++ imposedStatePlusView: " << std::endl;
      auto imposedStatePlusView2 = init::QInterpolated::view::create(imposedStatePlus[face]);
      for(int j = 0; j < 9; j++){
          std::cout << " (j: " << j << ") " <<  "(face: " << face << ") ";
          for (int i = 0; i < numberOfPoints; i++) {
              std::cout << imposedStatePlusView2(i,j) << " ";
          }
          std::cout << std::endl;
      }

      //*/
  } //End layerData.getNumberOfCells()-loop


    /*
    std::cout << "C++ imposedStatePlusView: " << std::endl;
    auto imposedStatePlusView2 = init::QInterpolated::view::create(imposedStatePlus[4]);
    for(int j = 0; j < 9; j++){
        std::cout << " (j: " << j << ") " <<  "(face: " << 4 << ") ";
        for (int i = 0; i < numberOfPoints; i++) {
            std::cout << imposedStatePlusView2(i,j) << " ";
        }
        std::cout << std::endl;
    }

    double getRupture_time[m_friction_data.numberOfPoints][m_friction_data.nFace];
    for(int j = 0; j < m_friction_data.nFace; j++){
        for (int i = 0; i < numberOfPoints; i++){
            getRupture_time[i][j] = m_friction_data.getRupture_time(i, j);
        }
    }
    m_friction_data.getRupture_time(0,0) = 100;
    e_interoperability.setFrictionOutput(m_friction_data);
    m_friction_data.getRupture_time(0,0) = -3;
    m_friction_data.getRupture_time(1,1) = -1;
    e_interoperability.getFrictionData(m_friction_data);
    double getRupture_time_old[m_friction_data.numberOfPoints][m_friction_data.nFace];
    for(int j = 0; j < m_friction_data.nFace; j++){
        for (int i = 0; i < numberOfPoints; i++){
            getRupture_time_old[i][j] = m_friction_data.getRupture_time(i, j);
        }
    }


     //*/
  struct seissol::physics::FrictionData fortran_data(m_friction_data.numberOfPoints, m_friction_data.nFace);
  e_interoperability.getFrictionData(fortran_data);
  m_friction_data.isEqualToFortran(fortran_data);

  m_loopStatistics->end(m_regionComputeDynamicRupture, layerData.getNumberOfCells());
}


void seissol::time_stepping::TimeCluster::computeDynamicRuptureFlops( seissol::initializers::Layer& layerData,
                                                                      long long&                    nonZeroFlops,
                                                                      long long&                    hardwareFlops )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;

  DRFaceInformation* faceInformation = layerData.var(m_dynRup->faceInformation);

  for (unsigned face = 0; face < layerData.getNumberOfCells(); ++face) {
    long long faceNonZeroFlops, faceHardwareFlops;
    m_dynamicRuptureKernel.flopsGodunovState( faceInformation[face], faceNonZeroFlops, faceHardwareFlops);

    nonZeroFlops += faceNonZeroFlops;
    hardwareFlops += faceHardwareFlops;
  }
}

#ifdef USE_MPI
/*
 * MPI-Communication during the simulation; exchange of DOFs.
 */
void seissol::time_stepping::TimeCluster::receiveGhostLayer(){
  SCOREP_USER_REGION( "receiveGhostLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  /*
   * Receive data of the ghost regions
   */
  for( unsigned int l_region = 0; l_region < m_meshStructure->numberOfRegions; l_region++ ) {
    // continue only if the cluster qualifies for communication
    if( m_resetLtsBuffers || m_meshStructure->neighboringClusters[l_region][1] <= static_cast<int>(m_globalClusterId) ) {
      // post receive request
      MPI_Irecv(   m_meshStructure->ghostRegions[l_region],                // initial address
                   m_meshStructure->ghostRegionSizes[l_region],            // number of elements in the receive buffer
                   MPI_C_REAL,                                               // datatype of each receive buffer element
                   m_meshStructure->neighboringClusters[l_region][0],      // rank of source
                   timeData+m_meshStructure->receiveIdentifiers[l_region], // message tag
                   seissol::MPI::mpi.comm(),                               // communicator
                   m_meshStructure->receiveRequests + l_region             // communication request
               );

      // add receive request to list of receives
      m_receiveQueue.push_back( m_meshStructure->receiveRequests + l_region );
    }
  }
}

void seissol::time_stepping::TimeCluster::sendCopyLayer(){
  SCOREP_USER_REGION( "sendCopyLayer", SCOREP_USER_REGION_TYPE_FUNCTION )

  /*
   * Send data of the copy regions
   */
  for( unsigned int l_region = 0; l_region < m_meshStructure->numberOfRegions; l_region++ ) {
    if( m_sendLtsBuffers || m_meshStructure->neighboringClusters[l_region][1] <= static_cast<int>(m_globalClusterId) ) {
      // post send request
      MPI_Isend(   m_meshStructure->copyRegions[l_region],              // initial address
                   m_meshStructure->copyRegionSizes[l_region],          // number of elements in the send buffer
                   MPI_C_REAL,                                            // datatype of each send buffer element
                   m_meshStructure->neighboringClusters[l_region][0],   // rank of destination
                   timeData+m_meshStructure->sendIdentifiers[l_region], // message tag
                   seissol::MPI::mpi.comm(),                            // communicator
                   m_meshStructure->sendRequests + l_region             // communication request
               );

      // add send request to list of sends
      m_sendQueue.push_back(m_meshStructure->sendRequests + l_region );
    }
  }
}

bool seissol::time_stepping::TimeCluster::testForGhostLayerReceives(){
  SCOREP_USER_REGION( "testForGhostLayerReceives", SCOREP_USER_REGION_TYPE_FUNCTION )

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  bool l_return;
  if (g_handleRecvs[m_clusterId] == 0) {
    l_return = true;
  } else {
    l_return = false;
  }
  return l_return;
#else
  // iterate over all pending receives
  for( std::list<MPI_Request*>::iterator l_receive = m_receiveQueue.begin(); l_receive != m_receiveQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the receive is complete
    MPI_Test( *l_receive, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending receives if completed
    if( l_mpiStatus == 1 )   l_receive = m_receiveQueue.erase( l_receive );
    // continue otherwise
    else                   ++l_receive;
  }

  // return true if the communication is finished
  return m_receiveQueue.empty();
#endif
}

bool seissol::time_stepping::TimeCluster::testForCopyLayerSends(){
  SCOREP_USER_REGION( "testForCopyLayerSends", SCOREP_USER_REGION_TYPE_FUNCTION )

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  bool l_return;
  if (g_handleSends[m_clusterId] == 0) {
    l_return = true;
  } else {
    l_return = false;
  }
  return l_return;
#else
  for( std::list<MPI_Request*>::iterator l_send = m_sendQueue.begin(); l_send != m_sendQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the send is complete
    MPI_Test( *l_send, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending sends if completed
    if( l_mpiStatus == 1 )   l_send = m_sendQueue.erase( l_send );
    // continue otherwise
    else                   ++l_send;
  }

  // return true if the communication is finished
  return m_sendQueue.empty();
#endif
}

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeCluster::initReceiveGhostLayer(){
  g_handleRecvs[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::initSendCopyLayer(){
  g_handleSends[m_clusterId] = 1;
}

void seissol::time_stepping::TimeCluster::waitForInits() {
  while( g_handleRecvs[m_clusterId] == 1 || g_handleSends[m_clusterId] == 1 );
}
#endif

#endif

void seissol::time_stepping::TimeCluster::computeLocalIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeLocalIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeLocalIntegration);

  // local integration buffer
  real l_integrationBuffer[tensor::I::size()] __attribute__((aligned(ALIGNMENT)));

  // pointer for the call of the ADER-function
  real* l_bufferPointer;

  real** buffers = i_layerData.var(m_lts->buffers);
  real** derivatives = i_layerData.var(m_lts->derivatives);
  real** displacements = i_layerData.var(m_lts->displacements);
  CellMaterialData* materialData = i_layerData.var(m_lts->material);

  kernels::LocalData::Loader loader;
  loader.load(*m_lts, i_layerData);
  kernels::LocalTmp tmp;

#ifdef _OPENMP
  #pragma omp parallel for private(l_bufferPointer, l_integrationBuffer, tmp) schedule(static)
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    auto data = loader.entry(l_cell);
    // overwrite cell buffer
    // TODO: Integrate this step into the kernel

    bool l_buffersProvided = (data.cellInformation.ltsSetup >> 8)%2 == 1; // buffers are provided
    bool l_resetBuffers = l_buffersProvided && ( (data.cellInformation.ltsSetup >> 10) %2 == 0 || m_resetLtsBuffers ); // they should be reset

    if (l_resetBuffers) {
      // assert presence of the buffer
      assert(buffers[l_cell] != nullptr);

      l_bufferPointer = buffers[l_cell];
    } else {
      // work on local buffer
      l_bufferPointer = l_integrationBuffer;
    }

    m_timeKernel.computeAder(m_timeStepWidth,
                             data,
                             tmp,
                             l_bufferPointer,
                             derivatives[l_cell]);

    // Compute local integrals (including some boundary conditions)
    CellBoundaryMapping (*boundaryMapping)[4] = i_layerData.var(m_lts->boundaryMapping);
    m_localKernel.computeIntegral(l_bufferPointer,
                                  data,
                                  tmp,
                                  &materialData[l_cell],
                                  &boundaryMapping[l_cell],
                                  m_fullUpdateTime,
                                  m_timeStepWidth
    );
    
    // Update displacement
    if (displacements[l_cell] != nullptr) {
      kernel::addVelocity krnl;
      krnl.I = l_bufferPointer;
      krnl.selectVelocity = init::selectVelocity::Values;
      krnl.displacement = displacements[l_cell];
      krnl.execute();
    }

    // update lts buffers if required
    // TODO: Integrate this step into the kernel
    if (!l_resetBuffers && l_buffersProvided) {
      assert (buffers[l_cell] != nullptr);

      for (unsigned int l_dof = 0; l_dof < tensor::I::size(); ++l_dof) {
        buffers[l_cell][l_dof] += l_integrationBuffer[l_dof];
      }
    }
  }

  m_loopStatistics->end(m_regionComputeLocalIntegration, i_layerData.getNumberOfCells());
}

void seissol::time_stepping::TimeCluster::computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData ) {
  SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

  m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

  real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
  CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
  CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
#ifdef USE_PLASTICITY
  PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
  real (*pstrain)[7] = i_layerData.var(m_lts->pstrain);
  unsigned numberOTetsWithPlasticYielding = 0;
#endif

  kernels::NeighborData::Loader loader;
  loader.load(*m_lts, i_layerData);

  real *l_timeIntegrated[4];
  real *l_faceNeighbors_prefetch[4];

#ifdef _OPENMP
#ifdef USE_PLASTICITY
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch) reduction(+:numberOTetsWithPlasticYielding)
#else
  #pragma omp parallel for schedule(static) private(l_timeIntegrated, l_faceNeighbors_prefetch)
#endif
#endif
  for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
    auto data = loader.entry(l_cell);
    seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                   data.cellInformation.ltsSetup,
                                                   data.cellInformation.faceTypes,
                                                   m_subTimeStart,
                                                   m_timeStepWidth,
                                                   faceNeighbors[l_cell],
#ifdef _OPENMP
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalData->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else
                                                   *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif
                                                   l_timeIntegrated);

#ifdef ENABLE_MATRIX_PREFETCH
#pragma message("the current prefetch structure (flux matrices and tDOFs is tuned for higher order and shouldn't be harmful for lower orders")
    l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][1] :
      drMapping[l_cell][1].godunov;
    l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][2] :
      drMapping[l_cell][2].godunov;
    l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != FaceType::dynamicRupture) ?
      faceNeighbors[l_cell][3] :
      drMapping[l_cell][3].godunov;

    // fourth face's prefetches
    if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
      l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != FaceType::dynamicRupture) ?
	faceNeighbors[l_cell+1][0] :
	drMapping[l_cell+1][0].godunov;
    } else {
      l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
    }
#endif

    m_neighborKernel.computeNeighborsIntegral( data,
                                               drMapping[l_cell],
#ifdef ENABLE_MATRIX_PREFETCH
                                               l_timeIntegrated, l_faceNeighbors_prefetch
#else
                                               l_timeIntegrated
#endif
                                               );

#ifdef USE_PLASTICITY
  numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_relaxTime,
                                                                                     m_timeStepWidth,
                                                                                     m_globalData,
                                                                                     &plasticity[l_cell],
                                                                                     data.dofs,
                                                                                     pstrain[l_cell] );
#endif
#ifdef INTEGRATE_QUANTITIES
  seissol::SeisSol::main.postProcessor().integrateQuantities( m_timeStepWidth,
                                                              i_layerData,
                                                              l_cell,
                                                              dofs[l_cell] );
#endif // INTEGRATE_QUANTITIES
  }

  #ifdef USE_PLASTICITY
  g_SeisSolNonZeroFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_nonZero[PlasticityCheck] + numberOTetsWithPlasticYielding * m_flops_nonZero[PlasticityYield];
  g_SeisSolHardwareFlopsPlasticity += i_layerData.getNumberOfCells() * m_flops_hardware[PlasticityCheck] + numberOTetsWithPlasticYielding * m_flops_hardware[PlasticityYield];
  #endif

  m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells());
}

#ifdef USE_MPI
bool seissol::time_stepping::TimeCluster::computeLocalCopy(){
  SCOREP_USER_REGION( "computeLocalCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localCopy ) {
    logError() << "Invalid call of computeLocalCopy, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // continue only if copy layer sends are complete
  if( !testForCopyLayerSends() ) return false;

  // post receive requests
#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  initReceiveGhostLayer();
#else
  receiveGhostLayer();
#endif

  // MPI checks for receiver writes receivers either in the copy layer or interior
  if( m_updatable.localInterior ) {
    writeReceivers();
  }

  // integrate copy layer locally
  computeLocalIntegration( m_clusterData->child<Copy>() );

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[LocalCopy];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[LocalCopy];

#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  initSendCopyLayer();
#else
  sendCopyLayer();
#endif

#ifndef USE_COMM_THREAD
  // continue with communication
  testForGhostLayerReceives();
#endif

  // compute sources, update simulation time
  if( !m_updatable.localInterior ) {
    computeSources();
    m_predictionTime += m_timeStepWidth;
  }

  // update finished
  m_updatable.localCopy  = false;

  // wait until communication thread finished initializing the receives
#if defined(_OPENMP) && defined(USE_COMM_THREAD)
  waitForInits();
#endif

  return true;
}
#endif

void seissol::time_stepping::TimeCluster::computeLocalInterior(){
  SCOREP_USER_REGION( "computeLocalInterior", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.localInterior ) {
    logError() << "Invalid call of computeLocalInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // MPI checks for receiver writes receivers either in the copy layer or interior
#ifdef USE_MPI
  if( m_updatable.localCopy ) {
    writeReceivers();
  }
#else
  // non-MPI checks for write in the interior
  writeReceivers();
#endif

  // integrate interior cells locally
  computeLocalIntegration( m_clusterData->child<Interior>() );

  g_SeisSolNonZeroFlopsLocal += m_flops_nonZero[LocalInterior];
  g_SeisSolHardwareFlopsLocal += m_flops_hardware[LocalInterior];

#ifdef USE_MPI
#ifndef USE_COMM_THREAD
  // continue with communication
  testForGhostLayerReceives();
  testForCopyLayerSends();
#endif
#endif

  // compute sources, update simulation time
  if( !m_updatable.localCopy ) {
    computeSources();
    m_predictionTime += m_timeStepWidth;
  }

  // update finished
  m_updatable.localInterior = false;
}

#ifdef USE_MPI
bool seissol::time_stepping::TimeCluster::computeNeighboringCopy() {
  SCOREP_USER_REGION( "computeNeighboringCopy", SCOREP_USER_REGION_TYPE_FUNCTION )

  // ensure a valid call
  if( !m_updatable.neighboringCopy ) {
    logError() << "Invalid call of computeNeighboringCopy aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  // continue only of ghost layer receives are complete
  if( !testForGhostLayerReceives() ) return false;

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  if (m_dynamicRuptureFaces == true) {
    if (m_updatable.neighboringInterior) {
      computeDynamicRupture(m_dynRupClusterData->child<Interior>());
      g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawInterior];
      g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawInterior];
    }

    computeDynamicRupture(m_dynRupClusterData->child<Copy>());
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawCopy];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawCopy];
  }

  computeNeighboringIntegration( m_clusterData->child<Copy>() );

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[NeighborCopy];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[NeighborCopy];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRNeighborCopy];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRNeighborCopy];

#ifndef USE_COMM_THREAD
  // continue with communication
  testForCopyLayerSends();
#endif

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringInterior ) {
    // First cluster calls fault receiver output
    // TODO: Change from iteration based to time based
    if (m_clusterId == 0) {
      e_interoperability.faultOutput( m_fullUpdateTime, m_timeStepWidth );
    }

    m_fullUpdateTime      += m_timeStepWidth;
    m_subTimeStart        += m_timeStepWidth;
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringCopy = false;

  return true;
}
#endif

void seissol::time_stepping::TimeCluster::computeNeighboringInterior() {
  SCOREP_USER_REGION( "computeNeighboringInterior", SCOREP_USER_REGION_TYPE_FUNCTION )
  // ensure a valid call
  if( !m_updatable.neighboringInterior ) {
    logError() << "Invalid call of computeNeighboringInterior, aborting:"
      << this             << m_clusterId      << m_globalClusterId << m_numberOfTimeSteps
      << m_fullUpdateTime << m_predictionTime << m_timeStepWidth   << m_subTimeStart      << m_resetLtsBuffers;
  }

  if (m_dynamicRuptureFaces == true && m_updatable.neighboringCopy == true) {
    computeDynamicRupture(m_dynRupClusterData->child<Interior>());
    g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRFrictionLawInterior];
    g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRFrictionLawInterior];
  }

  // Update all cells in the interior with the neighboring boundary contribution.
  computeNeighboringIntegration( m_clusterData->child<Interior>() );

  g_SeisSolNonZeroFlopsNeighbor += m_flops_nonZero[NeighborInterior];
  g_SeisSolHardwareFlopsNeighbor += m_flops_hardware[NeighborInterior];
  g_SeisSolNonZeroFlopsDynamicRupture += m_flops_nonZero[DRNeighborInterior];
  g_SeisSolHardwareFlopsDynamicRupture += m_flops_hardware[DRNeighborInterior];

  // compute dynamic rupture, update simulation time and statistics
  if( !m_updatable.neighboringCopy ) {
    // First cluster calls fault receiver output
    // TODO: Change from iteration based to time based
    if (m_clusterId == 0) {
      e_interoperability.faultOutput( m_fullUpdateTime, m_timeStepWidth );
    }

    m_fullUpdateTime      += m_timeStepWidth;
    m_subTimeStart        += m_timeStepWidth;
    m_numberOfFullUpdates += 1;
    m_numberOfTimeSteps   += 1;
  }

  // update finished
  m_updatable.neighboringInterior = false;
}

void seissol::time_stepping::TimeCluster::computeLocalIntegrationFlops( unsigned                    numberOfCells,
                                                                        CellLocalInformation const* cellInformation,
                                                                        long long&                  nonZeroFlops,
                                                                        long long&                  hardwareFlops  )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;

  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    unsigned cellNonZero, cellHardware;
    // TODO(Lukas) Maybe include avg. displacement computation here at some point.
    m_timeKernel.flopsAder(cellNonZero, cellHardware);
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
    m_localKernel.flopsIntegral(cellInformation[cell].faceTypes, cellNonZero, cellHardware);
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
  }
}

void seissol::time_stepping::TimeCluster::computeNeighborIntegrationFlops(  unsigned                    numberOfCells,
                                                                            CellLocalInformation const* cellInformation,
                                                                            CellDRMapping const       (*drMapping)[4],
                                                                            long long&                  nonZeroFlops,
                                                                            long long&                  hardwareFlops,
                                                                            long long&                  drNonZeroFlops,
                                                                            long long&                  drHardwareFlops )
{
  nonZeroFlops = 0;
  hardwareFlops = 0;
  drNonZeroFlops = 0;
  drHardwareFlops = 0;

  for (unsigned cell = 0; cell < numberOfCells; ++cell) {
    unsigned cellNonZero, cellHardware;
    long long cellDRNonZero, cellDRHardware;
    m_neighborKernel.flopsNeighborsIntegral(  cellInformation[cell].faceTypes,
                                              cellInformation[cell].faceRelations,
                                              drMapping[cell],
                                              cellNonZero,
                                              cellHardware,
                                              cellDRNonZero,
                                              cellDRHardware );
    nonZeroFlops += cellNonZero;
    hardwareFlops += cellHardware;
    drNonZeroFlops += cellDRNonZero;
    drHardwareFlops += cellDRHardware;

    /// \todo add lts time integration
    /// \todo add plasticity
  }
}

void seissol::time_stepping::TimeCluster::computeFlops()
{
#ifdef USE_MPI
  computeLocalIntegrationFlops( m_meshStructure->numberOfCopyCells,
                                m_clusterData->child<Copy>().var(m_lts->cellInformation),
                                m_flops_nonZero[LocalCopy],
                                m_flops_hardware[LocalCopy] );
#endif

  computeLocalIntegrationFlops( m_meshStructure->numberOfInteriorCells,
                                m_clusterData->child<Interior>().var(m_lts->cellInformation),
                                m_flops_nonZero[LocalInterior],
                                m_flops_hardware[LocalInterior] );

#ifdef USE_MPI
  computeNeighborIntegrationFlops(  m_meshStructure->numberOfCopyCells,
                                    m_clusterData->child<Copy>().var(m_lts->cellInformation),
                                    m_clusterData->child<Copy>().var(m_lts->drMapping),
                                    m_flops_nonZero[NeighborCopy],
                                    m_flops_hardware[NeighborCopy],
                                    m_flops_nonZero[DRNeighborCopy],
                                    m_flops_hardware[DRNeighborCopy] );
#endif

  computeNeighborIntegrationFlops(  m_meshStructure->numberOfInteriorCells,
                                    m_clusterData->child<Interior>().var(m_lts->cellInformation),
                                    m_clusterData->child<Interior>().var(m_lts->drMapping),
                                    m_flops_nonZero[NeighborInterior],
                                    m_flops_hardware[NeighborInterior],
                                    m_flops_nonZero[DRNeighborInterior],
                                    m_flops_hardware[DRNeighborInterior] );

  computeDynamicRuptureFlops( m_dynRupClusterData->child<Copy>(), m_flops_nonZero[DRFrictionLawCopy], m_flops_hardware[DRFrictionLawCopy] );
  computeDynamicRuptureFlops( m_dynRupClusterData->child<Interior>(), m_flops_nonZero[DRFrictionLawInterior], m_flops_hardware[DRFrictionLawInterior] );

  seissol::kernels::Plasticity::flopsPlasticity(  m_flops_nonZero[PlasticityCheck],
                                                  m_flops_hardware[PlasticityCheck],
                                                  m_flops_nonZero[PlasticityYield],
                                                  m_flops_hardware[PlasticityYield] );
}

#if defined(_OPENMP) && defined(USE_MPI) && defined(USE_COMM_THREAD)
void seissol::time_stepping::TimeCluster::pollForCopyLayerSends(){
  for( std::list<MPI_Request*>::iterator l_send = m_sendQueue.begin(); l_send != m_sendQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the send is complete
    MPI_Test( *l_send, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending sends if completed
    if( l_mpiStatus == 1 )   l_send = m_sendQueue.erase( l_send );
    // continue otherwise
    else                   ++l_send;
  }

  if (m_sendQueue.empty()) {
    g_handleSends[m_clusterId] = 0;
  }
}

void seissol::time_stepping::TimeCluster::pollForGhostLayerReceives(){
  // iterate over all pending receives
  for( std::list<MPI_Request*>::iterator l_receive = m_receiveQueue.begin(); l_receive != m_receiveQueue.end(); ) {
    int l_mpiStatus = 0;

    // check if the receive is complete
    MPI_Test( *l_receive, &l_mpiStatus, MPI_STATUS_IGNORE );

    // remove from list of pending receives if completed
    if( l_mpiStatus == 1 )   l_receive = m_receiveQueue.erase( l_receive );
    // continue otherwise
    else                   ++l_receive;
  }

  if (m_receiveQueue.empty()) {
    g_handleRecvs[m_clusterId] = 0;
  }
}

void seissol::time_stepping::TimeCluster::startReceiveGhostLayer() {
  receiveGhostLayer();
}

void seissol::time_stepping::TimeCluster::startSendCopyLayer() {
  sendCopyLayer();
}
#endif

