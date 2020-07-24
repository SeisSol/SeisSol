//
// Created by adrian on 02.06.20.
//
#ifndef TMPFRICTIONDATASTRUCT_CPP_
#define TMPFRICTIONDATASTRUCT_CPP_

#include <Initializer/typedefs.hpp>
#include <c++/8.3.0/iostream>



namespace seissol {
        namespace physics {
            struct TmpFrictionData{
                //TODO make members private
                const size_t numberOfPoints;
                const size_t nFace;
                const size_t nsize;
                bool initialized = false;
                bool allocated = false;
                bool tmpFrictionOnly = false;
                int function_call = 0;


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

                TmpFrictionData(const size_t numberOfPoints_in, const size_t nFace_in):
                    numberOfPoints(numberOfPoints_in),
                    nFace(nFace_in),
                    nsize(numberOfPoints_in * nFace_in){
                    if(!allocated){
                        allocated = true;
                        inst_healing= -1;
                        t_0= -1;
                        FL= -1;


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
                    }
                }

                TmpFrictionData(const TmpFrictionData& that) : nFace(0), nsize(0), numberOfPoints(0) {
                    //dont use copy constructor
                    assert(nFace == 0);
                    assert(nFace != 0);
                }

                TmpFrictionData& operator=(const TmpFrictionData& that) {
                    //dont assign
                    assert(nFace == 0);
                    assert(nFace != 0);
                    return *this;
                }

                ~TmpFrictionData() {
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
                }

                //Debugger Function
                bool isEqualToFortran(struct seissol::physics::TmpFrictionData &fortran_data){


                    bool b_initialStressInFaultCS = true;
                    bool b_cohesion = true;
                    bool b_D_C = true;
                    bool b_mu_S = true;
                    bool b_mu_D = true;
                    bool b_inst_healing = true;
                    bool b_t_0 = true;
                    bool b_FL = true;
                    bool b_forced_rupture_time = true;
                    bool b_magnitude_out = true;
                    bool b_mu = true;
                    bool b_slip = true;
                    bool b_slip1 = true;
                    bool b_slip2 = true;
                    bool b_slipRate1 = true;
                    bool b_slipRate2 = true;
                    bool b_rupture_time = true;
                    bool b_RF = true;
                    bool b_DS = true;
                    bool b_peakSR = true;

                    bool b_averaged_Slip = true;
                    bool b_dynStress_time = true;
                    bool b_tracXY = true;
                    bool b_tracXZ = true;

                    bool inputs = true;
                    bool outputs = true;


                    b_FL = fortran_data.FL == this->FL;
                    b_inst_healing = fortran_data.inst_healing == this->inst_healing;
                    b_t_0 = fortran_data.t_0 == this->t_0;
                    for (int i = 0; i < nFace; i++) {
                        if( fortran_data.magnitude_out[i] != magnitude_out[i])
                            b_magnitude_out = false;
                        if(magnitude_out[i] == true){
                            if(  fortran_data.averaged_Slip[i] != averaged_Slip[i]){
                                b_averaged_Slip = false;
                            }
                        }
                    }

                    for (int i = 0; i < nsize; i++) {
                        //inputs
                        if (fortran_data.initialStressInFaultCS[i] != initialStressInFaultCS[i])
                            b_initialStressInFaultCS = false;
                        if (fortran_data.cohesion[i] != cohesion[i])
                            b_cohesion = false;
                        if (fortran_data.D_C[i] != D_C[i])
                            b_D_C = false;
                        if (fortran_data.mu_S[i] != mu_S[i])
                            b_mu_S = false;
                        if (fortran_data.mu_D[i] != mu_D[i])
                            b_mu_D = false;
                        if (fortran_data.forced_rupture_time[i] != forced_rupture_time[i])
                            b_forced_rupture_time = false;

                        //outputs
                        if (fortran_data.mu[i] != mu[i]) {
                            b_mu = false;
                            //int nSides = i / numberOfPoints;
                            //std::cout << " mu not equal: ( iBnGP: " << i%numberOfPoints << " , iFace: " << nSides << " ) " << std::endl;
                            //std::cout << " mu c++ equal: " << mu[i] << " ,  mu fortran: " << fortran_data.mu[i] << " ) " << std::endl;
                        }

                        if (fortran_data.slip[i] != slip[i]){
                            b_slip = false;
                            int nSides = i / numberOfPoints;
                            //std::cout << " slip not equal: ( iBnGP: " << i % numberOfPoints << " , iFace: " << nSides << " ) " << std::endl;
                            //std::cout << " slip c++ equal: " << slip[i] << " ,  slip fortran: " << fortran_data.slip[i] << " ) " << std::endl;
                        }
                        if(  fortran_data.slip1[i] != slip1[i] ){
                            b_slip1 = false;
                            //int nSides = i / numberOfPoints;
                            //std::cout << " slip1 not equal: ( iBnGP: " << i%numberOfPoints << " , iFace: " << nSides << " ) " << std::endl;
                            //std::cout << " slip1 c++ equal: " << slip1[i] << " ,  slip1 fortran: " << fortran_data.slip1[i] << " ) " << std::endl;
                        }
                        if(  fortran_data.slip2[i] != slip2[i] )
                            b_slip2 = false;
                        if(  fortran_data.slipRate1[i] != slipRate1[i] )
                            b_slipRate1 = false;
                        if(  fortran_data.slipRate2[i] != slipRate2[i] )
                            b_slipRate2 = false;
                        if(  fortran_data.rupture_time[i] != rupture_time[i] )
                            b_rupture_time = false;
                        if(  fortran_data.RF[i] != RF[i] )
                            b_RF = false;
                        if(  fortran_data.DS[i] != DS[i] )
                            b_DS = false;
                        if(  fortran_data.peakSR[i] != peakSR[i] )
                            b_peakSR = false;
                        if(  fortran_data.dynStress_time[i] != dynStress_time[i] )
                            b_dynStress_time = false;
                        if(  fortran_data.tracXY[i] != tracXY[i] )
                            b_tracXY = false;
                        if(  fortran_data.tracXZ[i] != tracXZ[i] )
                            b_tracXZ = false;
                    }

                    inputs = (b_initialStressInFaultCS && b_cohesion && b_D_C && b_mu_S && b_mu_D && b_inst_healing &&
                            b_t_0 && b_FL && b_forced_rupture_time && b_magnitude_out);
                    outputs = (b_mu && b_slip  && b_slip1  && b_slip2  && b_slipRate1  && b_slipRate2  && b_rupture_time
                            && b_RF  && b_DS  && b_peakSR  && b_averaged_Slip &&
                            b_dynStress_time  && b_tracXY && b_tracXZ);

                    //assert( (inputs && outputs) == true );
                    return (inputs && outputs);
                }

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
                real* getCohesionFace(int iFace){
                    assert(iFace < nFace);
                    return &cohesion[iFace* numberOfPoints];
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

#endif