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
                const size_t nsize;
                bool initialized = false;

                int inst_healing;
                double t_0;
                int FL;

                //size nFace
                int *elem;
                int *side;
                bool *magnitude_out;
                real *averaged_Slip;

                //size [nBndGP][6][nFace]
                real *initialStressInFaultCS;

                //all size: [nBndGP][nFace]
                real *cohesion;
                real *D_C;
                real *mu_S;
                real *mu_D;
                real *forced_rupture_time;
                real *mu;
                real *slip;
                real *slip1;
                real *slip2;
                real *slipRate1;
                real *slipRate2;
                real *rupture_time;
                bool *RF;
                bool *DS;
                real *peakSR;
                real *dynStress_time;
                real *tracXY;
                real *tracXZ;


                FrictionData(const size_t numberOfPoints_in, const size_t nFace_in):
                    numberOfPoints(numberOfPoints_in),
                    nFace(nFace_in),
                    nsize(numberOfPoints_in * nFace_in){
                    if(!initialized){

                        initialized = true;
                        int inst_healing= -1;
                        double t_0= -1;
                        int FL= -1;


                        //size nFace
                        elem = new int [nFace];
                        side = new int [nFace];
                        magnitude_out = new bool [nFace];
                        averaged_Slip = new double [nFace];

                        cohesion = new double [nsize];
                        D_C = new double [nsize];
                        mu_S = new double [nsize];
                        mu_D = new double [nsize];
                        forced_rupture_time = new double [nsize];
                        mu = new double [nsize];
                        slip = new double [nsize];
                        slip1 = new double [nsize];
                        slip2 = new double [nsize];
                        slipRate1 = new double [nsize];
                        slipRate2 = new double [nsize];
                        rupture_time = new double [nsize];
                        RF = new bool [nsize];
                        DS = new bool [nsize];
                        peakSR = new double [nsize];
                        dynStress_time = new double [nsize];
                        tracXY = new double [nsize];
                        tracXZ = new double [nsize];

                        initialStressInFaultCS = new real [nsize*6];




                        /*
                         *
                        for (int i = 0; i < nFace; i++) {
                            averaged_Slip[i] = 0;
                            side
                            elem
                            magnitude_out

                        }

                         for (int i = 0; i < nFace; i++) {
                            initialStressInFaultCS[i] = new real *[6]();
                            for (int j = 0; j < 6; j++) {
                                initialStressInFaultCS[i][j] = new real[numberOfPoints]();
                            }
                        }
                        cohesion = new double *[nFace];
                        for (int i = 0; i < nFace; i++) {
                            //cohesion[i] = new double[numberOfPoints]();
                            D_C[i] = new double[numberOfPoints]();
                            mu_S[i] = new double[nFace]();      //testing
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

                         */

                    }


                }

                FrictionData(const FrictionData& that) : nFace(0), nsize(0), numberOfPoints(0) {
                    assert(nFace == 0);
                    assert(nFace != 0);
                }

                FrictionData& operator=(const FrictionData& that) {
                    assert(nFace == 0);
                    assert(nFace != 0);
                    return *this;
                }

                ~FrictionData() {
                    delete [] initialStressInFaultCS;
                    delete [] cohesion;
                    delete [] D_C;
                    delete [] mu_S;
                    delete [] mu_D;
                    delete [] forced_rupture_time;
                    delete [] mu;
                    delete [] slip;
                    delete [] slip1;
                    delete [] slip2;
                    delete [] slipRate1;
                    delete [] slipRate2;
                    delete [] rupture_time;
                    delete [] RF;
                    delete [] DS;
                    delete [] peakSR;
                    delete [] dynStress_time;
                    delete [] tracXY;
                    delete [] tracXZ;
                    delete [] magnitude_out;
                    delete [] elem;
                    delete [] side;
                    delete [] averaged_Slip;

                    //TODO: write destructor and copy constructor and copy assignment operator

                }
                //*/

                real getInitialStressInFaultCS(int iBndGP, int i,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    assert(i < 6);
                    return initialStressInFaultCS[iBndGP + i* numberOfPoints+ iFace*6*numberOfPoints];
                }
                real getCohesion(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return cohesion[iBndGP + iFace* numberOfPoints];
                }
                real getD_C(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return D_C[iBndGP + iFace* numberOfPoints];
                }
                real getMu_S(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return mu_S[iBndGP + iFace* numberOfPoints];
                }
                real getMu_D(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return mu_D[iBndGP + iFace* numberOfPoints];
                }
                real getforced_rupture_time(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return forced_rupture_time[iBndGP + iFace* numberOfPoints];
                }
                real& getMu(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return mu[iBndGP + iFace* numberOfPoints];
                }
                real setMu(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return mu[iBndGP + iFace* numberOfPoints];
                }
                real& getSlip(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip[iBndGP + iFace* numberOfPoints];
                }
                real setSlip(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip[iBndGP + iFace* numberOfPoints];
                }
                real& getSlip1(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip1[iBndGP + iFace* numberOfPoints];
                }
                real setSlip1(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip1[iBndGP + iFace* numberOfPoints];
                }
                real& getSlip2(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip2[iBndGP + iFace* numberOfPoints];
                }
                real setSlip2(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slip2[iBndGP + iFace* numberOfPoints];
                }
                real& getSlipRate1(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slipRate1[iBndGP + iFace* numberOfPoints];
                }
                real setSlipRate1(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slipRate1[iBndGP + iFace* numberOfPoints];
                }
                real& getSlipRate2(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slipRate2[iBndGP + iFace* numberOfPoints];
                }
                real setSlipRate2(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return slipRate2[iBndGP + iFace* numberOfPoints];
                }
                real& getRupture_time(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return rupture_time[iBndGP + iFace* numberOfPoints];
                }
                real setRupture_time(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return rupture_time[iBndGP + iFace* numberOfPoints];
                }
                bool getRF(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return RF[iBndGP + iFace* numberOfPoints];
                }
                void setRF(bool value, int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    RF[iBndGP + iFace* numberOfPoints] = value;
                }
                bool& getDS(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return DS[iBndGP + iFace* numberOfPoints];
                }
                void setDS(bool value, int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    RF[iBndGP + iFace* numberOfPoints];
                }
                real& getPeakSR(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return peakSR[iBndGP + iFace* numberOfPoints];
                }
                real setPeakSR(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return peakSR[iBndGP + iFace* numberOfPoints];
                }

                real& getDynStress_time(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return dynStress_time[iBndGP + iFace* numberOfPoints];
                }
                real setDynStress_time(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return dynStress_time[iBndGP + iFace* numberOfPoints];
                }
                real& getTracXY(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return tracXY[iBndGP + iFace* numberOfPoints];
                }
                real setTracXY(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return tracXY[iBndGP + iFace* numberOfPoints];
                }
                real& getTracXZ(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return tracXZ[iBndGP + iFace* numberOfPoints];
                }
                real setTracXZ(int iBndGP,int iFace){
                    assert(iBndGP < numberOfPoints);
                    assert(iFace < nFace);
                    return tracXZ[iBndGP + iFace* numberOfPoints];
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