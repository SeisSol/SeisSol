//
// Created by adrian on 02.06.20.
//
#ifndef TMPFRICTIONDATASTRUCT_CPP_
#define TMPFRICTIONDATASTRUCT_CPP_


#include <Initializer/typedefs.hpp>

namespace seissol {
        namespace physics {
            struct FrictionData{
                const size_t numberOfPoints;
                const size_t nFace;

                int *elem;
                int *side;
                real ***initialStressInFaultCS;
                real **cohesion;
                real **D_C;
                real **mu_S;
                real **mu_D;
                int inst_healing;
                double t_0;
                int FL;
                real **forced_rupture_time;
                bool *magnitude_out;
                real **mu;
                real **slip;
                real **slip1;
                real **slip2;
                real **slipRate1;
                real **slipRate2;
                real **rupture_time;
                bool **RF;
                bool **DS;
                real **PeakSR;
                real **dynStress_time;
                real **TracXY;
                real **TracXZ;

                real *averaged_Slip;

                FrictionData(const size_t numberOfPoints_in, const size_t nFace_in): numberOfPoints(numberOfPoints_in), nFace(nFace_in){

                    int inst_healing= -1;
                    double t_0= -1;
                    int FL= -1;

                    initialStressInFaultCS = new real **[nFace];
                    for (int i = 0; i < nFace; i++) {
                        initialStressInFaultCS[i] = new real *[6]();
                        for (int j = 0; j < 6; j++)
                            initialStressInFaultCS[i][j] = new real[numberOfPoints]();
                    }
                    magnitude_out = new bool [nFace];
                    elem = new int [nFace];
                    side = new int [nFace];
                    averaged_Slip = new double [nFace];
                    for (int i = 0; i < nFace; i++) {
                        averaged_Slip[i] = 0;
                    }

                    cohesion = new double *[nFace];
                    D_C = new double *[nFace];
                    mu_S = new double *[nFace];
                    mu_D = new double *[nFace];
                    forced_rupture_time = new double *[nFace];
                    mu = new double *[nFace];
                    slip = new double *[nFace];
                    slip1 = new double *[nFace];
                    slip2 = new double *[nFace];
                    slipRate1 = new double *[nFace];
                    slipRate2 = new double *[nFace];
                    rupture_time = new double *[nFace];
                    RF = new bool *[nFace];
                    DS = new bool *[nFace];
                    PeakSR = new double *[nFace];
                    dynStress_time = new double *[nFace];
                    TracXY = new double *[nFace];
                    TracXZ = new double *[nFace];

                    for (int i = 0; i < nFace; i++) {
                        cohesion[i] = new double[numberOfPoints]();
                        D_C[i] = new double[numberOfPoints]();
                        mu_S[i] = new double[numberOfPoints]();
                        mu_D[i] = new double[numberOfPoints]();
                        forced_rupture_time[i] = new double[numberOfPoints]();
                        mu[i] = new double[numberOfPoints]();
                        slip[i] = new double[numberOfPoints]();
                        slip1[i] = new double[numberOfPoints]();
                        slip2[i] = new double[numberOfPoints]();
                        slipRate1[i] = new double[numberOfPoints]();
                        slipRate2[i] = new double[numberOfPoints]();
                        rupture_time[i] = new double[numberOfPoints]();
                        RF[i] = new bool[numberOfPoints]();
                        DS[i] = new bool[numberOfPoints]();
                        PeakSR[i] = new double[numberOfPoints]();
                        dynStress_time[i] = new double[numberOfPoints]();
                        TracXY[i] = new double[numberOfPoints]();
                        TracXZ[i] = new double[numberOfPoints]();
                    }

                }
                ~FrictionData() {
                    //TODO: write destructor
                    /*
                    for (int i = 0; i < nFace; i++){
                        for (int j = 0; j < 6; j++)
                            delete[] initialStressInFaultCS[i][j];
                        delete[] initialStressInFaultCS[i];
                    }
                    delete[] initialStressInFaultCS;

                     */
                }

            };

    }
}


//in Initializer/typedefs.hpp:PlasticityData
/*
// plasticity information per cell
struct PlasticityData {
    // initial loading (stress tensor)
    real initialLoading[6];
    real cohesionTimesCosAngularFriction;
    real sinAngularFriction;
    real mufactor;
};
 */

/*
/src/Initializer/
DISC%DynRup%SlipRate1     = EQN%IniSlipRate1
DISC%DynRup%SlipRate2     = EQN%IniSlipRate2
DISC%DynRup%Slip          = 0.0D0
DISC%DynRup%Slip1         = 0.0D0
DISC%DynRup%Slip2         = 0.0D0
DISC%DynRup%TracXY        = 0.0D0
DISC%DynRup%TracXZ        = 0.0D0
DISC%DynRup%Mu(:,:)       = EQN%IniMu(:,:)
DISC%DynRup%StateVar(:,:) = EQN%IniStateVar
        DISC%DynRup%PeakSR        = 0.0D0
DISC%DynRup%rupture_time  = 0.0D0
DISC%DynRup%dynStress_time = 0.0D0

 */

#endif