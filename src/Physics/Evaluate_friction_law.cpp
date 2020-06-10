//
// Created by adrian on 14.05.20.
//

#include "Evaluate_friction_law.h"

//test
//protected
/*
m_data[0] = mu;
m_data[1] = slipRate1;
m_data[2] = slipRate2;
m_data[3] = slip;
m_data[4] = slip1;
m_data[5] = slip2;
m_data[6] = state;
m_data[7] = strength;
*/
//fault.data()
/*
 *     CellMaterialData* material = io_ltsTree->var(i_lts->material);
for (initializers::LTSTree::leaf_iterator it = dynRupTree->beginLeaf(initializers::LayerMask(Ghost)); it != dynRupTree->endLeaf(); ++it) {
    for (unsigned ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
        DRGodunovData*                        godunovData                                               = it->var(dynRup->godunovData);
        seissol::model::IsotropicWaveSpeeds*  waveSpeedsPlus                                            = it->var(dynRup->waveSpeedsPlus);
        seissol::model::IsotropicWaveSpeeds*  waveSpeedsMinus                                           = it->var(dynRup->waveSpeedsMinus);

    }
}
*/

/*
//Constructor
seissol::physics::Evaluate_friction_law ::Evaluate_friction_law()
{
}
//Destructor
seissol::physics::Evaluate_friction_law::~Evaluate_friction_law() {
}
*/


void seissol::physics::Evaluate_friction_law::Eval_friction_law(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int &iFace, int &iSide, int &iElem, double &time, double *timePoints,  // IN: element ID, time, inv Trafo
        double &rho, double &rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double resampleMatrix[],                                         //
        seissol::physics::FrictionData &friction_data
                                                   //global variables
){
    //required input
    int FL;     //EQN%FL
    int nBndGP;// DISC%Galerkin%nBndGP
    int nTimeGP; // DISC%Galerkin%nTimeGP
    void*EQN, *DISC, *MESH, *MPI,  *IO, *BND;

    //local variables
    double DeltaT[nTimeGP];

    DeltaT[0]=timePoints[0];
    for(int iTimeGP = 1; iTimeGP< nTimeGP; iTimeGP++ ){
        DeltaT[iTimeGP] = timePoints[iTimeGP]-timePoints[iTimeGP-1];
    }
    DeltaT[nTimeGP] = DeltaT[nTimeGP] + DeltaT[0];  // to fill last segment of Gaussian integration

    if(FL == 0){
        seissol::physics::Evaluate_friction_law::no_fault(XYStressGP,  XZStressGP,  TractionGP_XY, TractionGP_XZ);
    }else if(FL == 2 || FL == 16 ){
        seissol::physics::Evaluate_friction_law::Linear_slip_weakening_TPV1617(   TractionGP_XY, TractionGP_XZ, NorStressGP, XYStressGP, XZStressGP, iFace, iSide, iElem,
                nBndGP, nTimeGP, rho, rho_neig, w_speed, w_speed_neig, time, DeltaT, resampleMatrix, friction_data);
    }else if(FL == 6){
        seissol::physics::Evaluate_friction_law::Linear_slip_weakening_bimaterial( TractionGP_XY, TractionGP_XZ, NorStressGP, XYStressGP, XZStressGP, iFace, iSide, iElem,
                                                                                   nBndGP, nTimeGP, rho, rho_neig, w_speed, w_speed_neig, time, DeltaT, EQN, DISC, MESH, MPI, IO);
    }else if(FL == 33){

    }else if(FL == 3 || FL == 4){

    }else if(FL == 7){

    }else if(FL == 101){

    }else if(FL == 103){

    }else{
        //TODO:  logError(*) 'ERROR in friction.f90: friction law case',EQN%FL,' not implemented!'
    }
}

void seissol::physics::Evaluate_friction_law::no_fault(
        double **XYStressGP, double **XZStressGP,                                //INPUT
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
){
    TractionGP_XY = XYStressGP; //TODO: copy value not pointer value
    TractionGP_XZ = XZStressGP;
}

/*
 * Output: TractionGP_XY, TractionGP_XZ, DISC
 */

void seissol::physics::Evaluate_friction_law::Linear_slip_weakening_bimaterial(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,         //
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO
        //initializers::LTSTree* io_ltsTree, initializers::LTS*  i_lts, initializers::LTSTree* dynRupTree, initializers::DynamicRupture* dynRup,  //data structs
        //checkpoint::Fault fault
){
    //dummy:
    int nFace;

    //***********************************
    // GET THESE FROM required DATA STRUCT
    // input:
    double InitialStressInFaultCS[nBndGP][6][nFace]; // = EQN%InitialStressInFaultCS[i][1][iFace] //TODO: check right size and maybe all indecis need to be shifted by 1
    double Cohesion[nBndGP][nFace];   //DISC%DynRup%cohesion(nBndGP,iFace)
    double D_C[nBndGP][nFace];         // DISC%DynRup%D_C(nBndGP,iFace)
    double Mu_S[nBndGP][nFace]; //DISC%DynRup%Mu_S(nBndGP,iFace)
    double Mu_D[nBndGP][nFace]; //DISC%DynRup%Mu_D(nBndGP,iFace)
    bool inst_healing;    //DISC%DynRup%inst_healing //TODO is actually a bool?
    double v_star;          //DISC%DynRup%v_star         // !< reference velocity of prakash-cliff regularization
    double prakash_length;   //DISC%DynRup%L  reference length of prakash-cliff regularization

    //in and output:
    double Mu[nBndGP][nFace];           //DISC%DynRup%Mu(nBndGP,iFace)
    double Slip[nBndGP][nFace];   //DISC%DynRup%Slip(nBndGP,iFace)
    double Slip1[nBndGP][nFace];   //DISC%DynRup%Slip1(nBndGP,iFace)
    double Slip2[nBndGP][nFace];   //DISC%DynRup%Slip2(nBndGP,iFace)
    double SlipRate1[nBndGP][nFace];  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double SlipRate2[nBndGP][nFace];  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    bool DS[nBndGP][nFace];             //DISC%DynRup%DS(:,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double StrengthData[nBndGP][nFace];        //DISC%DynRup%Strength(iBndGP,iFace)

    //only output
    double dynStress_time[nBndGP][nFace]; //DISC%DynRup%dynStress_time(:,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    //***********************************


    //initialize local variables
    double LocMu, LocD_C, LocSlip, LocSlip1, LocSlip2, LocP;
    double LocMu_S, LocMu_D;
    double LocSR1, LocSR2;
    double cohesion;
    double Strength_exp;
    double P_0;
    double time_inc;
    double LocSR;
    double sigma;
    double ShTest;
    double Strength;
    double LocTracXY, LocTracXZ;


    for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++){ //loop over all points

        LocMu     = Mu[iBndGP][iFace];     // Current friction coefficient at given fault node
        LocMu_S   = Mu_S[iBndGP][iFace];  //DISC%DynRup%Mu_S(iBndGP,iFace)               //!< Static friction coefficient at given fault node - seisolxx.f90->ini_seisol.f90->readpar.f90->parameters.par found in fault.yaml
        LocMu_D   = Mu_D[iBndGP][iFace]; //DISC%DynRup%Mu_D(iBndGP,iFace)               //!< Dynamic friction coefficient at given fault node - found in fault.yaml
        LocD_C    = D_C[iBndGP][iFace]; //DISC%DynRup%D_C(iBndGP,iFace)               //!< Critical slip at given fault node - found in fault.yaml
        LocSlip   = Slip[iBndGP][iFace]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
        LocSlip1   = Slip1[iBndGP][iFace]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
        LocSlip2   = Slip2[iBndGP][iFace]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
        LocSR1    = SlipRate1[iBndGP][iFace]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSR2    = SlipRate2[iBndGP][iFace]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
        cohesion  = Cohesion[iBndGP][iFace]; //DISC%DynRup%cohesion(iBndGP,iFace)          // !< cohesion at given fault node  (should be negative since negative normal stress is compression)
        P_0       = InitialStressInFaultCS[iBndGP][0][iFace]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];
        Strength_exp = StrengthData[iBndGP][iFace]; //DISC%DynRup%Strength[iBndGP][iFace];         //!< save strength since it is used for bimaterial

        for(int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++){ //loop over time steps

            LocP   = NorStressGP[iBndGP][iTimeGP];
            time_inc = DeltaT[iTimeGP];
            //  modify strength according to prakash clifton
            LocSR = std::sqrt(LocSR1*LocSR1 + LocSR2*LocSR1);
            sigma = LocP+P_0;
            prak_clif_mod(Strength_exp, sigma, LocSR, v_star, prakash_length, LocMu, time_inc);
            ShTest = sqrt(std::pow(InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP], 2) + std::pow(InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP],2) );

            if(ShTest > Strength){
                // 1 evaluate friction
                LocTracXY = ((InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP])/ShTest)*Strength;
                LocTracXZ = ((InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP])/ShTest)*Strength;

                // 2 update stress change
                LocTracXY = LocTracXY - InitialStressInFaultCS[iBndGP][3][iFace];
                LocTracXZ = LocTracXZ - InitialStressInFaultCS[iBndGP][5][iFace];
            }else{
                LocTracXY = XYStressGP[iBndGP][iTimeGP];
                LocTracXZ = XZStressGP[iBndGP][iTimeGP];
            }
            //!Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
            LocSR1     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXY-XYStressGP[iBndGP][iTimeGP]);
            LocSR2     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXZ-XZStressGP[iBndGP][iTimeGP]);
            LocSR      = sqrt(LocSR1*LocSR1 + LocSR2*LocSR2);
            // Update slip
            LocSlip1 = LocSlip1 + LocSR1*time_inc;
            LocSlip2 = LocSlip2 + LocSR2*time_inc;
            LocSlip = LocSlip + LocSR*time_inc;
            if(abs(LocSlip) < LocD_C){
                LocMu = LocMu_S - (LocMu_S-LocMu_D)/LocD_C*abs(LocSlip);
            }else{
                LocMu = LocMu_D;
            }

            // instantaneous healing
            if(inst_healing == true){
                if(LocSR < u_0){
                    LocMu = LocMu_S;
                    // reset slip history for LSW
                    LocSlip = 0.0; //0.0D0
                }
            }

            //Save traction for flux computation
            TractionGP_XY[iBndGP][iTimeGP] = LocTracXY;
            TractionGP_XZ[iBndGP][iTimeGP] = LocTracXZ;

        } //End of time loop

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if(RF[iBndGP][iFace] && LocSR > 0.001){
            rupture_time[iBndGP][iFace]=time;
            RF[iBndGP][iFace] = false;
        }

        //output time when shear stress is equal to the dynamic stress after rupture arrived
        //currently only for linear slip weakening
        if ( (rupture_time[iBndGP][iFace] > 0.0) && (rupture_time[iBndGP][iFace] <= time)){
            if(DS[iBndGP][iFace] && abs(LocSlip) >= LocD_C){
                dynStress_time[iBndGP][iFace] = time;
                DS[iBndGP][iFace] = false;
            }
        }
        if(LocSR > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR;
        }

        Mu[iBndGP][iFace]       = LocMu;
        SlipRate1[iBndGP][iFace] = LocSR1;
        SlipRate2[iBndGP][iFace] = LocSR2;
        Slip[iBndGP][iFace]     = LocSlip;
        Slip1[iBndGP][iFace]     = LocSlip1;
        Slip2[iBndGP][iFace]     = LocSlip2;
        TracXY[iBndGP][iFace]    = LocTracXY;
        TracXZ[iBndGP][iFace]    = LocTracXZ;
        StrengthData[iBndGP][iFace]  = Strength_exp;

    }//End loop over all points
}


/*
 *
 */
void seissol::physics::Evaluate_friction_law::Linear_slip_weakening_TPV1617(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,
        double resampleMatrix[],
        seissol::physics::FrictionData &friction_data
){
    //dummy:
    int nFace;
    double matmul[nBndGP];
    double sum_tmpSlip;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    auto resampleMatrixView = init::resample::view::create(resampleMatrix);
    double t_0 = friction_data.t_0; //= DISC%DynRup%t_0
    double ***InitialStressInFaultCS = friction_data.initialStressInFaultCS; // = EQN%InitialStressInFaultCS[i][1][iFace] //TODO: check right size and maybe all indecis need to be shifted by 1
    double **cohesion = friction_data.cohesion;   //DISC%DynRup%cohesion(nBndGP,iFace)
    double **D_C = friction_data.D_C;         // DISC%DynRup%D_C(nBndGP,iFace)
    int FL = friction_data.FL;                  //EQN%FL
    double **forced_rupture_time = friction_data.forced_rupture_time; //DISC%DynRup%forced_rupture_time(nBndGP,iFace)
    double **Mu_S = friction_data.mu_S; //DISC%DynRup%Mu_S(nBndGP,iFace)
    double **Mu_D = friction_data.mu_D; //DISC%DynRup%Mu_D(nBndGP,iFace)
    bool inst_healing = friction_data.inst_healing;    //DISC%DynRup%inst_healing //TODO is actually a bool?
    bool *magnitude_out = friction_data.magnitude_out;        //DISC%DynRup%magnitude_out(iFace) //TODO is this a bool?

    //in and output:
    double **Mu = friction_data.mu;           //DISC%DynRup%Mu(nBndGP,iFace)
    double **Slip1 = friction_data.slip1; //DISC%DynRup%Slip1(nBndGP,iFace)
    double **Slip2 = friction_data.slip2;   //DISC%DynRup%Slip2(nBndGP,iFace)
    double **Slip = friction_data.slip;   //DISC%DynRup%Slip(nBndGP,iFace)
    bool **RF = friction_data.RF;            //DISC%DynRup%RF(nBndGP,iFace)
    double **rupture_time = friction_data.rupture_time;    //DISC%DynRup%rupture_time(nBndGP,iFace)
    bool **DS = friction_data.DS;             //DISC%DynRup%DS(:,iFace)
    double **PeakSR = friction_data.PeakSR;       //DISC%DynRup%PeakSR(:,iFace)
    double *averaged_Slip = friction_data.averaged_Slip;    //DISC%DynRup%averaged_Slip(iFace)

    //only output
    double **dynStress_time = friction_data.dynStress_time; //DISC%DynRup%dynStress_time(:,iFace)
    double **TracXY = friction_data.TracXY;           //DISC%DynRup%TracXY(:,iFace)
    double **TracXZ = friction_data.TracXY;       //DISC%DynRup%TracXZ(:,iFace)
    double **SlipRate1 = friction_data.slipRate1;  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double **SlipRate2 = friction_data.slipRate2;  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2
    //***********************************


    //initialize local variables
    double tmpSlip[nBndGP];
    for(int i = 0; i < nBndGP; i++){
        tmpSlip[i] = 0.0; //D0
    }
    double Z = rho * w_speed[2];
    double Z_neig = rho_neig * w_speed_neig[2];
    double eta = Z*Z_neig / (Z+Z_neig);
    double tn = time;
    double time_inc = 0;
    double P[nBndGP];
    double Strength[nBndGP][iFace];
    double ShTest[nBndGP];
    double LocSR[nBndGP];
    double LocSR1[nBndGP];
    double LocSR2[nBndGP];
    double LocTracXY[nBndGP];
    double LocTracXZ[nBndGP];
    double f1[nBndGP];
    double f2[nBndGP];

    for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {  //loop over time steps

            time_inc = DeltaT[iTimeGP];
            tn = tn + time_inc;


            P[iBndGP] = InitialStressInFaultCS[iBndGP][0][iFace] + NorStressGP[iBndGP][iTimeGP];
            Strength[iBndGP][iFace] = -cohesion[iBndGP][iFace] - Mu[iBndGP][iFace] * std::min(P[iBndGP], 0.0);
            ShTest[iBndGP] = std::sqrt(
                    std::pow(InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP], 2) +
                    std::pow(InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP], 2));
            LocSR[iBndGP] = std::max(0.0, (ShTest[iBndGP] - Strength[iBndGP][iFace]) / eta);
            LocSR1[iBndGP] = LocSR[iBndGP] * (InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP]) /
                             (Strength[iBndGP][iFace] + eta * LocSR[iBndGP]);
            LocSR2[iBndGP] = LocSR[iBndGP] * (InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP]) /
                             (Strength[iBndGP][iFace] + eta * LocSR[iBndGP]);
            LocTracXY[iBndGP] = XYStressGP[iBndGP][iTimeGP] - eta * LocSR1[iBndGP];
            LocTracXZ[iBndGP] = XZStressGP[iBndGP][iTimeGP] - eta * LocSR2[iBndGP];

            //Update slip
            Slip1[iBndGP][iFace] = Slip1[iBndGP][iFace] + LocSR1[iBndGP] * time_inc;
            Slip2[iBndGP][iFace] = Slip2[iBndGP][iFace] + LocSR2[iBndGP] * time_inc;

            //Resample slip-rate, such that the state (Slip) lies in the same polynomial space as the degrees of freedom
            //resampleMatrix first projects LocSR on the two-dimensional basis on the reference triangle with
            //degree less or equal than CONVERGENCE_ORDER-1, and then evaluates the polynomial at the quadrature points
            matmul[iBndGP] = 0;

            for (int j = 0; j < nBndGP; j++) {
                //TODO: deck if resampleMatrix is symmetric  = dimension (nBndGP,nBndGP)
                matmul[iBndGP] += resampleMatrixView(iBndGP,j) * LocSR[j];
            }
            Slip[iBndGP][iFace] = Slip[iBndGP][iFace] + matmul[iBndGP] * time_inc;
            tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSR[iBndGP] * time_inc;

            //Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
            f1[iBndGP] = std::min(std::abs(Slip[iBndGP][iFace]) / D_C[iBndGP][iFace], 1.0);

            if (FL == 16) {
                if (t_0 == 0) {
                    if (tn >= forced_rupture_time[iBndGP][iFace]) {
                        f2[iBndGP] = 1.0;
                    } else {
                        f2[iBndGP] = 0.0;
                    }
                } else {
                    //TODO: why unreachable?
                    f2[iBndGP] = std::max(0.0, std::min((time - forced_rupture_time[iBndGP][iFace]) / t_0, 1.0));
                }
            } else {
                f2[iBndGP] = 0.0;
            }

            Mu[iBndGP][iFace] = Mu_S[iBndGP][iFace] -
                                (Mu_S[iBndGP][iFace] - Mu_D[iBndGP][iFace]) * std::max(f1[iBndGP], f2[iBndGP]);

            //instantaneous healing
            if (inst_healing == true) {
                if (LocSR[iBndGP] < u_0) {
                    Mu[iBndGP][iFace] = Mu_S[iBndGP][iFace];
                    Slip[iBndGP][iFace] = 0.0;
                }
            }

            TractionGP_XY[iBndGP][iTimeGP] = LocTracXY[iBndGP];
            TractionGP_XZ[iBndGP][iTimeGP] = LocTracXZ[iBndGP];
        }     //end iTimeGP loop

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR[iBndGP] > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }
        //output time when shear stress is equal to the dynamic stress after rupture arrived
        //currently only for linear slip weakening
        if (rupture_time[iBndGP][iFace] > 0.0 &&
                rupture_time[iBndGP][iFace] <= time &&
                DS[iBndGP][iFace] &&
                std::abs(Slip[iBndGP][iFace] >= D_C[iBndGP][iFace]))
        {
            dynStress_time[iBndGP][iFace] = time;
            DS[iBndGP][iFace] = false;
        }

        if(LocSR[iBndGP] > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR[iBndGP];
        }

        TracXY[iBndGP][iFace] = LocTracXY[iBndGP];
        TracXZ[iBndGP][iFace] = LocTracXZ[iBndGP];
        SlipRate1[iBndGP][iFace] = LocSR1[iBndGP];
        SlipRate2[iBndGP][iFace] = LocSR2[iBndGP];

    } //end i < nBndGP loop

    //---compute and store slip to determine the magnitude of an earthquake ---
    //    to this end, here the slip is computed and averaged per element
    //    in calc_seissol.f90 this value will be multiplied by the element surface
    //    and an output happened once at the end of the simulation
    if(magnitude_out[iFace] ){
        for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
            sum_tmpSlip += tmpSlip[iBndGP];
        }
        averaged_Slip[iFace] = averaged_Slip[iFace] + sum_tmpSlip / nBndGP;
    }
}



/*
 * !< T. Ulrich 27.07.17
 * !< This friction law allows imposing a slip rate on the DR boundary
 */
void seissol::physics::Evaluate_friction_law::ImposedSlipRateOnDRBoundary(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,         //
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO
){
    //dummy:
    int nFace;
    double sum_tmpSlip;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    double Tnuc; //DISC%DynRup%t_0
    double NucleationStressInFaultCS[nBndGP][6][nFace]; //EQN%NucleationStressInFaultCS(iBndGP,6?,iFace)
    double Slip1[nBndGP][nFace]; //DISC%DynRup%Slip1(iBndGP,iFace)
    double Slip2[nBndGP][nFace];  //DISC%DynRup%Slip2(iBndGP,iFace)
    double Slip[nBndGP][nFace];  //DISC%DynRup%Slip(iBndGP,iFace)
    bool magnitude_out[nFace];        //DISC%DynRup%magnitude_out(iFace) //TODO is this a bool?

    //in and output:
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double averaged_Slip[nFace];    //DISC%DynRup%averaged_Slip(iFace)

    //only output
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    //***********************************


    //initialize local variables
    double dt = 0.0;
    double tmpSlip[nBndGP];
    for(int i = 0; i < nBndGP; i++){
        tmpSlip[i] = 0.0; //D0
    }
    double Gnuc;
    double time_inc;
    double LocTracXY[nBndGP];
    double LocTracXZ[nBndGP];
    double SlipRate1[nBndGP][nFace];
    double SlipRate2[nBndGP][nFace];
    double LocSR[nBndGP];

    double eta = (w_speed[2]*rho*w_speed_neig[2]*rho_neig) / (w_speed[2]*rho + w_speed_neig[2]*rho_neig);
    double tn = time;
    for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
        dt += DeltaT[iTimeGP];
    }

    for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
            time_inc = DeltaT[iTimeGP];
            tn=tn + time_inc;
            Gnuc = Calc_SmoothStepIncrement(tn, Tnuc, time_inc)/time_inc;

            // EQN%NucleationStressInFaultCS (1 and 2) contains the slip in FaultCS
            LocTracXY[iBndGP]  = XYStressGP[iBndGP][iTimeGP] - eta * NucleationStressInFaultCS[iBndGP][1][iFace]*Gnuc;
            LocTracXZ[iBndGP] =  XZStressGP[iBndGP][iTimeGP]  - eta * NucleationStressInFaultCS[iBndGP][2][iFace]*Gnuc;
            SlipRate1[iBndGP][iFace]     = NucleationStressInFaultCS[iBndGP][1][iFace]*Gnuc;
            SlipRate2[iBndGP][iFace]    = NucleationStressInFaultCS[iBndGP][2][iFace]*Gnuc;
            LocSR[nBndGP]    = std::sqrt(std::pow(SlipRate1[iBndGP][iFace],2) + std::pow(SlipRate2[iBndGP][iFace],2));
            // Update slip
            Slip1[iBndGP][iFace] += SlipRate1[iBndGP][iFace]*time_inc;
            Slip2[iBndGP][iFace] += SlipRate2[iBndGP][iFace]*time_inc;
            Slip[iBndGP][iFace] += LocSR[iBndGP]*time_inc;
            tmpSlip[iBndGP] = tmpSlip[iBndGP] + LocSR[iBndGP]*time_inc;

            TractionGP_XY[iBndGP][iTimeGP] = LocTracXY[iBndGP];
            TractionGP_XZ[iBndGP][iTimeGP] = LocTracXZ[iBndGP];
        }   //end time loop
        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR[iBndGP] > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }

        if(LocSR[iBndGP] > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR[iBndGP];
        }

        TracXY[iBndGP][iFace] = LocTracXY[iBndGP];
        TracXZ[iBndGP][iFace] = LocTracXZ[iBndGP];

        //---compute and store slip to determine the magnitude of an earthquake ---
        //   to this end, here the slip is computed and averaged per element
        //   in calc_seissol.f90 this value will be multiplied by the element surface
        //   and an output happened once at the end of the simulation
        if(magnitude_out[iFace] ){
            for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
                sum_tmpSlip += tmpSlip[iBndGP];
            }
            averaged_Slip[iFace] = averaged_Slip[iFace] + sum_tmpSlip / nBndGP;
        }

    }// end of iBndGP-loop
}

/*
 *     !> friction case 3,4: rate and state friction
 *     !> aging (3) and slip law (4)
 */
void seissol::physics::Evaluate_friction_law::rate_and_state(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO
){
    //dummy:
    int nFace;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    double InitialStressInFaultCS[nBndGP][6][nFace]; //EQN%InitialStressInFaultCS(iBndGP,6?,iFace)
    double Cohesion[nBndGP][nFace];   //DISC%DynRup%cohesion(nBndGP,iFace)
    double RS_f0;  //DISC%DynRup%RS_f0  !< Reference friction coefficient
    double RS_a;    //DISC%DynRup%RS_a  !< RS constitutive parameter "a"
    double RS_b;       //DISC%DynRup%RS_b  !< RS constitutive parameter "b"
    double RS_sl0;  //DISC%DynRup%RS_sl0     !< Reference slip
    double RS_sr0;  //DISC%DynRup%RS_sr0     !< Reference slip rate
    int FL;     //EQN%FL

    //in and output:
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double StateVar[nBndGP][nFace];   //DISC%DynRup%StateVar(iBndGP,iFace)
    double Slip1[nBndGP][nFace]; //DISC%DynRup%Slip1(iBndGP,iFace)
    double Slip2[nBndGP][nFace];  //DISC%DynRup%Slip2(iBndGP,iFace)
    double Slip[nBndGP][nFace];  //DISC%DynRup%Slip(iBndGP,iFace)
    double SlipRate1[nBndGP][nFace];  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double SlipRate2[nBndGP][nFace];  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2

    //only output
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    double Mu[nBndGP][nFace];         //DISC%DynRup%Mu(iBndGP,iFace)
    //***********************************

    //initialize local variables
    double LocSlip, LocSlip1, LocSlip2, LocSR, LocSR1, LocSR2, LocSV;
    double cohesion;
    double P;
    double P_0;
    double LocP;
    double time_inc;
    double ShTest;
    double SV0;
    int nSRupdates, nSVupdates;
    double tmp;
    double tmp2;
    double SRtest;
    double NR;
    double dNR;
    double LocMu;
    double LocTracXY;
    double LocTracXZ;


    for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {

        LocSlip   = Slip[iBndGP][iFace]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
        LocSlip1   = Slip1[iBndGP][iFace]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
        LocSlip2   = Slip2[iBndGP][iFace]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
        LocSR1    = SlipRate1[iBndGP][iFace]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSR2    = SlipRate2[iBndGP][iFace]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSV   = StateVar[iBndGP][iFace];     //DISC%DynRup%StateVar(iBndGP,iFace)
        cohesion  = Cohesion[iBndGP][iFace]; //DISC%DynRup%cohesion(iBndGP,iFace)          // !< cohesion at given fault node  (should be negative since negative normal stress is compression)
        P_0       = InitialStressInFaultCS[iBndGP][0][iFace]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];
        for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
            LocP   = NorStressGP[iBndGP][iTimeGP];
            time_inc = DeltaT[iTimeGP];

            //SignSR1   = SIGN(1.0,LocSR1)                    ! Gets the sign of the slip rate
            //SignSR2   = SIGN(1.0,LocSR2)                    ! Gets the sign of the slip rate

            // load traction and normal stress
            P = LocP+P_0;
            ShTest = std::sqrt(
                    std::pow(InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP], 2) +
                    std::pow(InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP], 2));

            // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996) //TODO: look up
            // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

            SV0=LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

            // The following process is adapted from that described by Kaneko et al. (2008) TODO: look up
            nSRupdates = 5; //TODO: can be put outside of loop
            nSVupdates = 2;

            LocSR      = std::sqrt(std::pow(LocSR1,2) + std::pow(LocSR2,2));
            tmp        = abs(LocSR);
            for(int j = 0; j < nSVupdates; j++){ //!This loop corrects SV values
                LocSR=abs(LocSR);
                if(FL == 3){    //aging law
                    LocSV=SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
                }else if(FL == 4){  //slip law
                    LocSV=RS_sl0/tmp*pow(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));    //TODO: loop up if pow is right
                }
                // Newton-Raphson algorithm to determine the value of the slip rate.
                // We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
                //  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
                // In our case we equalize the values of the traction for two equations:
                //             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
                //             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
                //               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))

                SRtest=LocSR;   // We use as first guess the SR value of the previous time step
                for(int i = 0; i < nSRupdates; i++){   //!This loop corrects SR values
                    tmp          = 0.5/RS_sr0* exp( (RS_f0+RS_b*log(RS_sr0*LocSV/RS_sl0) ) /RS_a);
                    tmp2         = tmp*SRtest;
                    NR           = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) *
                            (abs(P)*RS_a*log(tmp2+sqrt(pow(tmp2,2)+1.0))-ShTest)-SRtest;              //!TODO:not sure if ShTest should be + or -...
                    dNR          = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) *
                            (abs(P)*RS_a/sqrt(1+pow(tmp2,2))*tmp)-1.0;
                    SRtest = abs(SRtest-NR/dNR);             // no ABS needed around NR/dNR at least for aging law
                }   // End
                tmp=0.5*(LocSR+abs(SRtest));  //! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
                LocSR=abs(SRtest);
            }   // End SV-Loop

            if(FL == 3){    //aging law
                LocSV=SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
            }else if(FL == 4){  //slip law
                LocSV=RS_sl0/tmp*pow(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));    //TODO: loop up if pow is right
            }
            tmp  = 0.5 * (LocSR)/RS_sr0 * exp((RS_f0 + RS_b*log(RS_sr0*LocSV/RS_sl0)) / RS_a);
            LocMu    = RS_a * log(tmp + sqrt(pow(tmp,2) + 1.0));
            // 2D:
            // LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
            // LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
            // update stress change
            LocTracXY = -((InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP])/ShTest)*(LocMu*P+abs(cohesion));
            LocTracXZ = -((InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP])/ShTest)*(LocMu*P+abs(cohesion));
            LocTracXY = LocTracXY - InitialStressInFaultCS[iBndGP][3][iFace];
            LocTracXZ = LocTracXZ - InitialStressInFaultCS[iBndGP][5][iFace];

            // Compute slip
            LocSlip   = LocSlip  + (LocSR)*time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

            //Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
            LocSR1     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXY-XYStressGP[iBndGP][iTimeGP]);
            LocSR2     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXZ-XZStressGP[iBndGP][iTimeGP]);

            LocSlip1   = LocSlip1  + (LocSR1)*time_inc;
            LocSlip2   = LocSlip2  + (LocSR2)*time_inc;

            //LocSR1     = SignSR1*ABS(LocSR1)
            //LocSR2     = SignSR2*ABS(LocSR2)

            //Save traction for flux computation
            TractionGP_XY[iBndGP,iTimeGP] = &LocTracXY;
            TractionGP_XZ[iBndGP,iTimeGP] = &LocTracXZ;
        } //End of iTimeGP loop

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }
        if(LocSR > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR;
        }

        Mu[iBndGP][iFace]       = LocMu;
        SlipRate1[iBndGP][iFace]  = LocSR1;
        SlipRate2[iBndGP][iFace]  = LocSR2;
        Slip[iBndGP][iFace]       = LocSlip;
        Slip1[iBndGP][iFace]     = LocSlip1;
        Slip2[iBndGP][iFace]      = LocSlip2;
        StateVar[iBndGP][iFace]   = LocSV;
        TracXY[iBndGP][iFace] = LocTracXY;
        TracXZ[iBndGP][iFace] = LocTracXZ;
    }//End of iBndGP-loop
}

/*
 *     !> friction case 7: severe velocity weakening rate and state friction
 *     !< after Ampuero and Ben-Zion 2008
 */
void seissol::physics::Evaluate_friction_law::rate_and_state_vw(
        double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO
) {
    //dummy:
    int nFace;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    double InitialStressInFaultCS[nBndGP][6][nFace]; //EQN%InitialStressInFaultCS(iBndGP,6?,iFace)
    //double Cohesion[nBndGP][nFace];   //DISC%DynRup%cohesion(nBndGP,iFace)
    double RS_f0;  //DISC%DynRup%RS_f0  !< Reference friction coefficient, equivalent to static friction coefficient
    double RS_a;    //DISC%DynRup%RS_a  !< RS constitutive parameter "a", direct effect
    double RS_b;       //DISC%DynRup%RS_b  !< RS constitutive parameter "b", evolution effect
    double RS_sl0;  //DISC%DynRup%RS_sl0     !< Reference slip  , Dc, char. lengt scale
    double RS_sr0;  //DISC%DynRup%RS_sr0     !< Reference slip rate, Vc, char. velocity scale
    //int FL;     //EQN%FL

    //in and output:
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double StateVar[nBndGP][nFace];   //DISC%DynRup%StateVar(iBndGP,iFace)
    double Slip1[nBndGP][nFace]; //DISC%DynRup%Slip1(iBndGP,iFace)
    double Slip2[nBndGP][nFace];  //DISC%DynRup%Slip2(iBndGP,iFace)
    double Slip[nBndGP][nFace];  //DISC%DynRup%Slip(iBndGP,iFace)
    double SlipRate1[nBndGP][nFace];  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double SlipRate2[nBndGP][nFace];  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2

    //only output
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    double Mu[nBndGP][nFace];         //DISC%DynRup%Mu(iBndGP,iFace)
    //***********************************

    //initialize local variables
    double LocSlip, LocSlip1, LocSlip2, LocSR1, LocSR2, LocSR;
    double LocSV;
    double P_0;
    double LocP;
    double time_inc;
    double P;
    double ShTest;
    double SV0;
    double Tc;
    double coeft;
    int nSRupdates, nSVupdates;
    double SRtest;
    double tmp;
    double NR, dNR;
    double LocMu;
    double LocTracXY, LocTracXZ;

    for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {

        LocSlip = Slip[iBndGP][iFace]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
        LocSlip1 = Slip1[iBndGP][iFace]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
        LocSlip2 = Slip2[iBndGP][iFace]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
        LocSR1 = SlipRate1[iBndGP][iFace]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSR2 = SlipRate2[iBndGP][iFace]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSV = StateVar[iBndGP][iFace];     //DISC%DynRup%StateVar(iBndGP,iFace)
        P_0 = InitialStressInFaultCS[iBndGP][0][iFace]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];
        for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
            LocP = NorStressGP[iBndGP][iTimeGP];
            time_inc = DeltaT[iTimeGP];


            // load traction and normal stress
            P = LocP + P_0;
            ShTest = std::sqrt(
                    std::pow(InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP], 2) +
                    std::pow(InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP], 2));

            // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996) //TODO: look up
            // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

            SV0 = LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

            // The following process is adapted from that described by Kaneko et al. (2008) TODO: look up
            nSRupdates = 5; //TODO: can be put outside of loop
            nSVupdates = 2;

            LocSR = std::sqrt(std::pow(LocSR1, 2) + std::pow(LocSR2, 2)); //can be put outside of the loop

            //charact. time scale Tc
            Tc = RS_sl0 / RS_sr0;
            // exponent
            coeft= exp(-time_inc / Tc);

            for (int j = 0; j < nSVupdates; j++) { //!This loop corrects SV values
                LocSR = abs(LocSR);
                //exact integration assuming constant V in this loop
                LocSV=Tc*LocSR*(1.0-coeft) + coeft*SV0;
                //! Newton-Raphson algorithm to determine the value of the slip rate.
                //! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
                //!  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
                //! In our case we equalize the values of the traction for two equations:
                //!             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
                //!             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
                //!             where mu = mu_s + a V/(V+Vc) - b SV/(SV + Vc)
                //!
                SRtest = LocSR;   // We use as first guess the SR value of the previous time step
                for (int i = 0; i < nSRupdates; i++) {   //!This loop corrects SR values
                    tmp = RS_f0+RS_a*SRtest/(SRtest+RS_sr0)-RS_b*LocSV/(LocSV+RS_sl0); //=mu
                    NR = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) * (abs(P)*tmp-ShTest)-SRtest;
                    dNR          = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) *
                            (abs(P)*(RS_a/(SRtest+RS_sr0)-RS_a*SRtest/pow(SRtest+RS_sr0,2))) -1.0;
                    SRtest = SRtest-NR/dNR;
                }   // End nSRupdates-Loop
                tmp=0.5*(LocSR+abs(SRtest));  // For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
                LocSR=abs(SRtest);
            }   // End nSVupdates-Loop -  This loop corrects SV values

            LocSV    = Tc*tmp*(1-coeft) + coeft*SV0;
            tmp = 0.5 * (LocSR)/RS_sr0 * exp((RS_f0 + RS_b*log(RS_sr0*LocSV/RS_sl0)) / RS_a);

            //! Ampuero and Ben-Zion 2008 (eq. 1):
            // LocMu = friction coefficient (mu_f)
            // RS_f0 = static coefficient (mu_s)
            // RS_a = positive coefficient, quantifying  the direct effect (alpha)
            // LocSR = slip rate (V)
            // RS_sr0 = characteristic velocity scale (V_c)
            // RS_b = positive coefficient, quantifying  the evolution effect (beta)
            // RS_sl0 = characteristic velocity scale (V_c)
            LocMu = RS_f0+RS_a*LocSR/(LocSR+RS_sr0)-RS_b*LocSV/(LocSV+RS_sl0);


            // update stress change
            LocTracXY = -((InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP]) / ShTest) *LocMu * P;
            LocTracXZ = -((InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP]) / ShTest) *LocMu * P;
            LocTracXY = LocTracXY - InitialStressInFaultCS[iBndGP][3][iFace];
            LocTracXZ = LocTracXZ - InitialStressInFaultCS[iBndGP][5][iFace];

            // Compute slip
            LocSlip = LocSlip + LocSR * time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

            //Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
            LocSR1 = -(1.0 / (w_speed[2] * rho) + 1.0 / (w_speed_neig[2] * rho_neig)) * (LocTracXY - XYStressGP[iBndGP][iTimeGP]);
            LocSR2 = -(1.0 / (w_speed[2] * rho) + 1.0 / (w_speed_neig[2] * rho_neig)) * (LocTracXZ - XZStressGP[iBndGP][iTimeGP]);

            LocSlip1 = LocSlip1 + (LocSR1) * time_inc;
            LocSlip2 = LocSlip2 + (LocSR2) * time_inc;

            //Save traction for flux computation
            TractionGP_XY[iBndGP, iTimeGP] = &LocTracXY;
            TractionGP_XZ[iBndGP, iTimeGP] = &LocTracXZ;
        } //End of iTimeGP loop

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }
        if (LocSR > PeakSR[iBndGP][iFace]) {
            PeakSR[iBndGP][iFace] = LocSR;
        }

        Mu[iBndGP][iFace] = LocMu;
        SlipRate1[iBndGP][iFace] = LocSR1;
        SlipRate2[iBndGP][iFace] = LocSR2;
        Slip[iBndGP][iFace] = LocSlip;
        Slip1[iBndGP][iFace] = LocSlip1;
        Slip2[iBndGP][iFace] = LocSlip2;
        StateVar[iBndGP][iFace] = LocSV;
        TracXY[iBndGP][iFace] = LocTracXY;
        TracXZ[iBndGP][iFace] = LocTracXZ;
    }//End of iBndGP-loop
}


/*
 *     !> special friction case for SCEC TPV101: rate and state friction
 *     !> aging law
 *     !< with time and space dependent nucleation
 */
void seissol::physics::Evaluate_friction_law::rate_and_state_nuc101(
        double **TractionGP_XY, double **TractionGP_XZ,                     // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT, double iT[],                           //IN: time, inv Trafo
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO, void *BND
){
    //dummy:
    int nFace;
    int unkown;
    int nObject;
    int nSide = 4;
    int nElem;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    auto resampleMatrixView = init::resample::view::create(iT);
    int Face[nFace][unkown][unkown];              //MESH%Fault%Face(iFace,1,1)  !<Assigns each element's side in domain a fault-plane-index
    int BoundaryToObject[unkown][unkown];    //MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor) //!<Mapping from element and side index to boundary object nr
    int MPINumber[unkown][unkown];      //MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)   !<Index into MPI communication structure
    double NeighborCoords[unkown][unkown][unkown][nObject];  //BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex);     !< Vertex coordinates of elements adjacent to MPI boundary (outside)
    double xyNode[unkown][unkown];      //MESH%VRTX%xyNode(1,Vertex(1:4,iElem))     !<Coordinates of the nodes
    int Vertex[4][nElem];           //MESH%ELEM%Vertex(1:4,iElem)   !<Connection index to element's vertices
    double BndGP_Tri[unkown][nBndGP]; //MESH%ELEM%BndGP_Tri(1,iBndGP)    !<GaussWeights in 2D boundary
    double ShearXY_0;               //EQN%ShearXY_0      !< Initial shear stress
    int GlobalElemType;             //MESH%GlobalElemType    !<Triangle=3, quads=4, tets = 4, hex = 6, mixed = 7
    double IniBulk_xx[nBndGP][nFace];    //EQN%IniBulk_xx(iBndGP,iFace)  !< Initial bulk stress at fault
    double IniBulk_yy[nBndGP][nFace]; //EQN%IniBulk_yy(:,iFace)
    double IniBulk_zz[nBndGP][nFace]; //EQN%IniBulk_zz(:,iFace)
    double IniShearYZ[nBndGP][nFace]; //EQN%IniShearYZ(:,iFace)
    double IniShearXZ[nBndGP][nFace]; //EQN%IniShearXZ(:,iFace)
    double RS_f0;  //DISC%DynRup%RS_f0  !< Reference friction coefficient, equivalent to static friction coefficient
    //double RS_a;    //DISC%DynRup%RS_a  !< RS constitutive parameter "a", direct effect
    double RS_b;       //DISC%DynRup%RS_b  !< RS constitutive parameter "b", evolution effect
    double RS_sl0;  //DISC%DynRup%RS_sl0     !< Reference slip  , Dc, char. lengt scale
    double RS_sr0;  //DISC%DynRup%RS_sr0     !< Reference slip rate, Vc, char. velocity scale
    double RS_a_array[nFace][nBndGP]; //DISC%DynRup%RS_a_array(iFace,iBndGP)
    int FL;     //EQN%FL

    //in and output:
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double StateVar[nBndGP][nFace];   //DISC%DynRup%StateVar(iBndGP,iFace)
    double Slip1[nBndGP][nFace]; //DISC%DynRup%Slip1(iBndGP,iFace)
    double Slip2[nBndGP][nFace];  //DISC%DynRup%Slip2(iBndGP,iFace)
    double Slip[nBndGP][nFace];  //DISC%DynRup%Slip(iBndGP,iFace)
    double SlipRate1[nBndGP][nFace];  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double SlipRate2[nBndGP][nFace];  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2
    double IniShearXY[nFace][nBndGP];   // EQN%IniShearXY(iFace,iBndGP)      !< Initial shear stress at fault

    //only output
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    double Mu[nBndGP][nFace];         //DISC%DynRup%Mu(iBndGP,iFace)

    //unkown origin???
    int VertexSide[4][3];
    //double iT[6][6];            function argument // inverse Transformation matrix   = resampleMatrix
    //***********************************

    //initialize local variables
    bool nodewise;
    double Rnuc;
    double Tnuc;
    int iNeighbor, iLocalNeighborSide;
    int iObject;
    int MPIIndex;
    double xV[4];
    double yV[4];
    double zV[4];
    double chi, tau, xi,eta,zeta;
    double xGP, yGP,zGP;
    double radius;
    double Fnuc;
    double Gnuc;
    double xp[GlobalElemType], yp[GlobalElemType], zp[GlobalElemType];
    double maxvalX = 0.0;
    double maxvalZ = 0.0;
    double Stress[6][nBndGP];
    double matmul;
    double LocSlip, LocSlip1, LocSlip2, LocSR1, LocSR2, LocSV, P_0;
    double LocP, time_inc;
    double RS_a;
    double P, ShTest;
    double SV0;
    int nSRupdates, nSVupdates;
    double LocSR, tmp, tmp2, NR, dNR;
    double SRtest;
    double LocMu, LocTracXY, LocTracXZ;




    // switch for Gauss node wise stress assignment
    nodewise; // = true;    //TODO maybe external config?

    //Apply time dependent nucleation at global time step not sub time steps for simplicity
    //initialize time and space dependent nucleation
    Rnuc=3000.0;
    Tnuc=1.0;

    if(time <= Tnuc){
        if(nodewise){
            // Gauss node coordinate definition and stress assignment
            // get vertices of complete tet
            if(Face[iFace][0][0] == 0){
                // iElem is in the neighbor domain
                // The neighbor element belongs to a different MPI domain
                iNeighbor = Face[iFace][0][1];          // iNeighbor denotes "-" side
                iLocalNeighborSide  = Face[iFace][1][1];
                iObject  = BoundaryToObject[iLocalNeighborSide][iNeighbor];
                MPIIndex = MPINumber[iLocalNeighborSide][iNeighbor];
                for(int i = 0; i < 4; i++){
                    xV[i] = NeighborCoords[0][i][MPIIndex][iObject];    //indexing shifted
                    yV[i] = NeighborCoords[1][i][MPIIndex][iObject];
                    zV[i] = NeighborCoords[2][i][MPIIndex][iObject];
                }
            }else{
                // get vertices
                for(int i = 0; i < 4; i++){
                    int vtx = Vertex[i][iElem];
                    xV[i] = xyNode[0][vtx]; //indexing shifted
                    yV[i] = xyNode[1][vtx];
                    zV[i] = xyNode[2][vtx];
                }
            }
            for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
                // Transformation of boundary GP's into XYZ coordinate system
                chi  = BndGP_Tri[0][iBndGP];
                tau  = BndGP_Tri[1][iBndGP];
                TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0);
                TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV);

                //radial distance to hypocenter
                radius = sqrt(pow(xGP,2)+pow(zGP+7500.0,2));
                // Inside nucleation patch add shear stress perturbation of 25 MPa along strike

                if(radius < Rnuc){
                    Fnuc=exp(  pow(radius,2) / ( pow(radius,2)-pow(Rnuc,2) )  );
                    if(time < 0.0){
                        Gnuc=exp( pow(time-Tnuc,2)/(time*(time-2.0*Tnuc)) );
                    }else{
                        Gnuc=0.0;
                    }
                    IniShearXY[iFace][iBndGP] = ShearXY_0 + 25.0e6*Fnuc*Gnuc;
                }
            } //End iBndGP-loop
        }else{  //if not nodewise
            // get coordinates needed for nucleation zone
            if(iElem != 0){
                for(int j = 0; j < 3; j++){
                    xp[j] = xyNode[0][Vertex[ VertexSide[iSide][j] ] [iElem] ]; //indexing shifted
                    yp[j] = xyNode[1][Vertex[ VertexSide[iSide][j] ] [iElem] ];
                    zp[j] = xyNode[2][Vertex[ VertexSide[iSide][j] ] [iElem] ];
                }
            }else if(iElem == 0){   //in case "+" element is not present in the local domain
                iLocalNeighborSide = Face[iFace][1][1]; //shifting index
                for(int j = 0; j < 3; j++){
                    xp[j] = xyNode[0][   Vertex[ VertexSide[iLocalNeighborSide][j] ][ Face[iFace][0][1] ]  ];
                    yp[j] = xyNode[1][   Vertex[ VertexSide[iLocalNeighborSide][j] ][ Face[iFace][0][1] ]  ];
                    zp[j] = xyNode[2][   Vertex[ VertexSide[iLocalNeighborSide][j] ][ Face[iFace][0][1] ]  ];
                }
            }

            //max radial distance to hypocenter
            //TODO: max value calc: (maybe as external function)
            for(int j = 0; j < 3; j++){
                maxvalX = std::max(maxvalX, xp[j]);
                maxvalZ = std::max(maxvalZ, zp[j]);
            }
            radius = sqrt( pow( maxvalX,2) + pow((maxvalZ+7500.0),2) );
            //   Inside nucleation patch add shear stress perturbation of 25 MPa along strike
            if(radius < Rnuc){
                Fnuc = exp(pow(radius,2)/(pow(radius,2)-pow(Rnuc,2)) );
                if(time > 0.0){
                    Gnuc=exp( pow(time-Tnuc,2) / (time*(time-2.0*Tnuc)) );
                }else{
                    Gnuc=0.0;
                }
                for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
                    IniShearXY[iFace][iBndGP] = ShearXY_0 + 25.0e6*Fnuc*Gnuc;
                }
            }   //end if Rnuc
        } //end if nodewise
    } //end if time <= Tnuc

    //Background stress rotation to face's reference system
    for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        Stress[0][iBndGP]= IniBulk_xx[iBndGP][iFace];   //shifted index
        Stress[1][iBndGP]= IniBulk_yy[iBndGP][iFace];
        Stress[2][iBndGP]= IniBulk_zz[iBndGP][iFace];
        Stress[3][iBndGP]= IniShearXY[iBndGP][iFace];
        Stress[4][iBndGP]= IniShearYZ[iBndGP][iFace];
        Stress[5][iBndGP]= IniShearXZ[iBndGP][iFace];
        matmul = 0;
        for(int j = 0; j <6; j++){
            for(int k = 0; k <6; k++){
                matmul += resampleMatrixView(j,k) * Stress[k][iBndGP];
            }
            Stress[j][iBndGP] = matmul;
        }
    }

    for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        LocSlip = Slip[iBndGP][iFace];      //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
        LocSlip1 = Slip1[iBndGP][iFace];    //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
        LocSlip2 = Slip2[iBndGP][iFace];    //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
        LocSR1 = SlipRate1[iBndGP][iFace];  //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSR2 = SlipRate2[iBndGP][iFace];      //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
        LocSV = StateVar[iBndGP][iFace];     //DISC%DynRup%StateVar(iBndGP,iFace)
        P_0 = Stress[0][iBndGP]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];

        for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
            LocP   = NorStressGP[iBndGP][iTimeGP];
            time_inc = DeltaT[iTimeGP];
            // spatially variabel a
            RS_a   = RS_a_array[iFace][iBndGP];

            P = LocP + P_0;
            ShTest = std::sqrt(
                    std::pow(Stress[iBndGP][3] + XYStressGP[iBndGP][iTimeGP], 2) +
                    std::pow(Stress[iBndGP][5]+ XZStressGP[iBndGP][iTimeGP], 2));

            // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
            // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

            SV0=LocSV;    // Careful, the SV must always be corrected using SV0 and not LocSV!

            // The following process is adapted from that described by Kaneko et al. (2008) TODO: look up
            nSRupdates = 5; //TODO: can be put outside of loop
            nSVupdates = 2;

            LocSR      = std::sqrt(std::pow(LocSR1,2) + std::pow(LocSR2,2));
            tmp        = abs(LocSR);

            for(int j = 0; j < nSVupdates; j++){ //!This loop corrects SV values
                LocSR=abs(LocSR);
                if(FL == 102){    //aging law
                    LocSV=SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
                }else if(FL == 104){  //slip law
                    LocSV=RS_sl0/tmp*pow(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));    //TODO: loop up if pow is right
                }
                // Newton-Raphson algorithm to determine the value of the slip rate.
                // We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
                //  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
                // In our case we equalize the values of the traction for two equations:
                //             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
                //             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
                //               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))

                SRtest=LocSR;  // We use as first guess the SR value of the previous time step

                for(int i = 0; i < nSRupdates; i++){   //!This loop corrects SR values
                    tmp          = 0.5/RS_sr0* exp( (RS_f0+RS_b*log(RS_sr0*LocSV/RS_sl0) ) /RS_a);
                    tmp2         = tmp*SRtest;
                    NR           = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) *                     //TODO: shift index
                                   (abs(P)*RS_a*log(tmp2+sqrt(pow(tmp2,2)+1.0))-ShTest)-SRtest;              //!TODO:not sure if ShTest should be + or -...
                    dNR          = -(1.0/w_speed[2]/rho+1.0/w_speed_neig[2]/rho_neig) *
                                   (abs(P)*RS_a/sqrt(1+pow(tmp2,2))*tmp)-1.0;
                    SRtest = SRtest-NR/dNR;             // no ABS needed around NR/dNR at least for aging law
                }   // End nSRupdates-loop
                tmp=0.5*(LocSR+abs(SRtest));  //! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
                LocSR=abs(SRtest);
            }   // End SV-Loop  j=1,nSVupdates   !This loop corrects SV values



            if(FL == 102){    //aging law
                LocSV=SV0*exp(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-exp(-tmp*time_inc/RS_sl0));
            }else if(FL == 104){  //slip law
                LocSV=RS_sl0/tmp*pow(tmp*SV0/RS_sl0, exp(-tmp*time_inc/RS_sl0));    //TODO: loop up if pow is right
            }
            tmp  = 0.5 * (LocSR)/RS_sr0 * exp((RS_f0 + RS_b*log(RS_sr0*LocSV/RS_sl0)) / RS_a);

            LocMu    = RS_a * log(tmp + sqrt(pow(tmp,2) + 1.0));
            // 2D:
            // LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
            // LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results

            // update stress change
            LocTracXY = -((Stress[iBndGP][3] + XYStressGP[iBndGP][iTimeGP])/ShTest)*(LocMu*P);
            LocTracXZ = -((Stress[iBndGP][5] + XZStressGP[iBndGP][iTimeGP])/ShTest)*(LocMu*P);
            LocTracXY = LocTracXY - Stress[iBndGP][3];
            LocTracXZ = LocTracXZ - Stress[iBndGP][5];

            // Compute slip
            LocSlip   = LocSlip  + (LocSR)*time_inc; // ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

            //Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
            LocSR1     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXY-XYStressGP[iBndGP][iTimeGP]);        //TODO: shift index
            LocSR2     = -(1.0/(w_speed[2]*rho)+1.0/(w_speed_neig[2]*rho_neig))*(LocTracXZ-XZStressGP[iBndGP][iTimeGP]);

            LocSlip1   = LocSlip1  + (LocSR1)*time_inc;
            LocSlip2   = LocSlip2  + (LocSR2)*time_inc;
            //LocSR1     = SignSR1*ABS(LocSR1)
            //LocSR2     = SignSR2*ABS(LocSR2)

            //Save traction for flux computation
            TractionGP_XY[iBndGP][iTimeGP] = LocTracXY;
            TractionGP_XZ[iBndGP][iTimeGP] = LocTracXZ;

        } // End of iTimeGP loop

        // output rupture front
        // outside of iTimeGP loop in order to safe an 'if' in a loop
        // this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }
        if(LocSR > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR;
        }

        Mu[iBndGP][iFace]       = LocMu;
        SlipRate1[iBndGP][iFace]  = LocSR1;
        SlipRate2[iBndGP][iFace]  = LocSR2;
        Slip[iBndGP][iFace]       = LocSlip;
        Slip1[iBndGP][iFace]     = LocSlip1;
        Slip2[iBndGP][iFace]      = LocSlip2;
        StateVar[iBndGP][iFace]   = LocSV;
        TracXY[iBndGP][iFace] = LocTracXY;
        TracXZ[iBndGP][iFace] = LocTracXZ;

    }// End of iBndGP loop
}


void seissol::physics::Evaluate_friction_law::rate_and_state_nuc103(
        double **TractionGP_XY, double **TractionGP_XZ,                     // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
        double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
        int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
        double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
        double time, double *DeltaT,                                        //IN: time
        double resampleMatrix[],
        void *EQN, void *DISC, void *MESH, void *MPI, void *IO, void *BND
){
    //dummy:
    int nFace;
    int unkown;
//    int nObject;
//    int nSide = 4;
//    int nElem;

    //***********************************
    // GET THESE FROM DATA STRUCT
    //required input:
    auto resampleMatrixView = init::resample::view::create(resampleMatrix);
    int TP_grid_nz; //DISC%dynRup%TP_grid_nz;  !< number of grid points to solve advection for TP in z-direction
    double t_0 = 0; //= DISC%DynRup%t_0
    double InitialStressInFaultCS[nBndGP][6][nFace]; //EQN%InitialStressInFaultCS(iBndGP,6?,iFace)
    double NucleationStressInFaultCS[nBndGP][6][nFace]; //EQN%NucleationStressInFaultCS(iBndGP,6?,iFace)
    double RS_f0;  //DISC%DynRup%RS_f0  !< Reference friction coefficient, equivalent to static friction coefficient
    double RS_b;       //DISC%DynRup%RS_b  !< RS constitutive parameter "b", evolution effect
    double RS_sl0;  //DISC%DynRup%RS_sl0     !< Reference slip  , Dc, char. lengt scale
    double RS_sr0;  //DISC%DynRup%RS_sr0     !< Reference slip rate, Vc, char. velocity scale
    double RS_a_array[nFace][nBndGP]; //DISC%DynRup%RS_a_array(iFace,iBndGP)
    double Mu_w; //= DISC%DynRup%Mu_w      ! mu_w, weakening friction coefficient
    double RS_srW_array[nFace][nBndGP]; //DISC%DynRup%RS_srW_array(:,iFace)    ! Vw, weakening sliding velocity, space dependent
    int ThermalPress; //DISC%DynRup%ThermalPress   !< thermal pressurization switch
    double TP[nBndGP][nFace][unkown];     //DISC%DynRup%TP(:,iFace,2)    !< Temperature and Pressure for TP along each fault point
    double TP_Theta[nBndGP][nFace][TP_grid_nz]; //DISC%DynRup%TP_Theta(iBndGP, iFace,:) !< Fourier transformed pressure
    double TP_sigma[nBndGP][nFace][TP_grid_nz]; //DISC%DynRup%TP_sigma(iBndGP, iFace,:) !< Fourier transformed temperature
    double TP_half_width_shear_zone[nBndGP][nFace];  //DISC%DynRup%TP_half_width_shear_zone(iBndGP,iFace)    !< spatial dependent half width of the shearing layer for TP
    double alpha_th;        //DISC%DynRup%alpha_th   !< thermal diffusion parameter for TP
    double alpha_hy[nBndGP][nFace];    //DISC%DynRup%alpha_hy(iBndGP,iFace)    !< spatial dependent hydraulic diffusion parameter for TP
    double rho_c;        //DISC%DynRup%rho_c !< heat capacity for TP
    double TP_Lambda;          // DISC%DynRup%TP_Lambda    !< pore pressure increase per unit increase
    double TP_grid[TP_grid_nz];    //DISC%DynRup%TP_grid   !< grid for TP
    double TP_DFinv[TP_grid_nz]; //DISC%DynRup%TP_DFinv  !< inverse Fourier coefficients
    double temp_0;          //EQN%Temp_0     !< Initial temperature for TP
    double pressure_0;            //EQN%Pressure_0           !< Initial pressure for TP
    bool magnitude_out[nFace];        //DISC%DynRup%magnitude_out(iFace) //TODO is this a bool?

    //in and output:
    bool DS[nBndGP][nFace];             //DISC%DynRup%DS(:,iFace)
    bool RF[nBndGP][nFace];            //DISC%DynRup%RF(nBndGP,iFace)
    double PeakSR[nBndGP][nFace];       //DISC%DynRup%PeakSR(:,iFace)
    double StateVar[nBndGP][nFace];   //DISC%DynRup%StateVar(iBndGP,iFace)
    double Slip1[nBndGP][nFace]; //DISC%DynRup%Slip1(iBndGP,iFace)
    double Slip2[nBndGP][nFace];  //DISC%DynRup%Slip2(iBndGP,iFace)
    double Slip[nBndGP][nFace];  //DISC%DynRup%Slip(iBndGP,iFace)
    double SlipRate1[nBndGP][nFace];  //DISC%DynRup%SlipRate1(:,iFace) = LocSR1
    double SlipRate2[nBndGP][nFace];  //DISC%DynRup%SlipRate2(:,iFace) = LocSR2

    //only output
    double rupture_time[nBndGP][nFace];     //DISC%DynRup%rupture_time(nBndGP,iFace)
    double TracXY[nBndGP][nFace];           //DISC%DynRup%TracXY(:,iFace)
    double TracXZ[nBndGP][nFace];       //DISC%DynRup%TracXZ(:,iFace)
    double Mu[nBndGP][nFace];         //DISC%DynRup%Mu(iBndGP,iFace)
    double dynStress_time[nBndGP][nFace]; //DISC%DynRup%dynStress_time(:,iFace)
    double averaged_Slip[nFace];    //DISC%DynRup%averaged_Slip(iFace)

    //***********************************

    //initialize local variables
    double Tnuc;
    int nSRupdates, nSVupdates;
    double dt;
    double Gnuc;
    double LocSlip[nBndGP], LocSlip1[nBndGP], LocSlip2[nBndGP], LocSR[nBndGP], LocSR1[nBndGP], LocSR2[nBndGP], LocSV[nBndGP], LocMu[nBndGP], P_0[nBndGP];
    double LocP[nBndGP], time_inc;
    double RS_fw, RS_srW[nBndGP], RS_a[nBndGP] ;
    double P[nBndGP], ShTest[nBndGP] , SV0[nBndGP];
    double SR_tmp[nBndGP], invZ;
    double P_f[nBndGP];
    double tmp[nBndGP], tmp2[nBndGP];
    double S[nBndGP];
    double Theta_tmp[TP_grid_nz], Sigma_tmp[TP_grid_nz];
    double n_stress[nBndGP] ;
    double SRtest[nBndGP];
    bool has_converged;
    double tmpSlip[nBndGP];
    double LocTracXY[nBndGP], LocTracXZ[nBndGP];
    double matmul;
    double sum_tmpSlip;


    // switch for Gauss node wise stress assignment
    bool nodewise; //= true;    //TODO: configureable?

    for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        tmpSlip[iBndGP] = 0.0;
    }


    //Apply time dependent nucleation at global time step not sub time steps for simplicity
    //initialize time and space dependent nucleation
    Tnuc = t_0;

    //!TU 7.07.16: if the SR is too close to zero, we will have problems (NaN)
    //!as a consequence, the SR is affected the AlmostZero value when too small
    double AlmostZero = 1e-45; //d-45;

    //!PARAMETERS of THE optimisation loops
    //!absolute tolerance on the function to be optimzed
    //! This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be the most adapted
    //! Number of iteration in the loops
    nSRupdates = 60;
    nSVupdates = 2;

    dt = 0;
    for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
        dt += DeltaT[iTimeGP];
    }
    if (time <= Tnuc){
        Gnuc = Calc_SmoothStepIncrement(time, Tnuc, dt);

        //DISC%DynRup%NucBulk_** is already in fault coordinate system

        for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
            for (int i = 0; i < 6; i++) {
                InitialStressInFaultCS[iBndGP][i][iFace]=InitialStressInFaultCS[iBndGP][i][iFace]+NucleationStressInFaultCS[iBndGP][i][iFace]*Gnuc;
            }
        }

    } //end If-Tnuc

    for (int iTimeGP = 0; iTimeGP < nTimeGP; iTimeGP++) {
        for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {

            LocMu[iBndGP] = Mu[iBndGP][iFace];     // Current friction coefficient at given fault node
            LocSlip[iBndGP] = Slip[iBndGP][iFace]; //DISC%DynRup%Slip(iBndGP,iFace)              //!< Slip path at given fault node
            LocSlip1[iBndGP] = Slip1[iBndGP][iFace]; //DISC%DynRup%Slip1(iBndGP,iFace)            //!< Slip at given fault node along loc dir 1
            LocSlip2[iBndGP] = Slip2[iBndGP][iFace]; //DISC%DynRup%Slip2(iBndGP,iFace)            // !< Slip at given fault node along loc dir 2
            LocSR1[iBndGP] = SlipRate1[iBndGP][iFace]; //DISC%DynRup%SlipRate1(iBndGP,iFace)         // !< Slip Rate at given fault node
            LocSR2[iBndGP] = SlipRate2[iBndGP][iFace]; //DISC%DynRup%SlipRate2(iBndGP,iFace)         // !< Slip Rate at given fault node
            P_0[iBndGP] = InitialStressInFaultCS[iBndGP][0][iFace]; //EQN%InitialStressInFaultCS[iBndGP][1][iFace];
            LocSV[iBndGP] = StateVar[iBndGP][iFace];     //DISC%DynRup%StateVar(iBndGP,iFace)



            // friction develops as                    mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
            // state variable SV develops as        dSV / dt = -(V - L) * (SV - SV_ss)
            //                                      SV_ss = a * ln[ 2*V0/V * sinh(mu_ss/a) ]
            //                                      mu_ss = mu_w + [mu_lv - mu_w] / [ 1 + (V/Vw)^8 ] ^ (1/8) ]
            //                                      mu_lv = mu_0 - (b-a) ln (V/V0)

            LocP[iBndGP] = NorStressGP[iBndGP][iTimeGP];
            time_inc = DeltaT[iTimeGP];

            //RS_f0  //= DISC%DynRup%RS_f0     ! mu_0, reference friction coefficient
            //RS_sr0 //= DISC%DynRup%RS_sr0    ! V0, reference velocity scale
            RS_fw = Mu_w;              //= DISC%DynRup%Mu_w      ! mu_w, weakening friction coefficient
            RS_srW[iBndGP] = RS_srW_array[iBndGP][iFace];   //! Vw, weakening sliding velocity, space dependent
            RS_a[iBndGP] = RS_a_array[iBndGP][iFace];  //! a, direct effect, space dependent
            //RS_b   //= DISC%DynRup%RS_b       ! b, evolution effect
            //RS_sl0 //= DISC%DynRup%RS_sl0_array(:,iFace)     ! L, char. length scale


            // load traction and normal stress
            P[iBndGP] = LocP[iBndGP] + P_0[iBndGP];
            ShTest[iBndGP] = std::sqrt(
                    std::pow(InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP], 2) +
                    std::pow(InitialStressInFaultCS[iBndGP][5][iFace] + XZStressGP[iBndGP][iTimeGP], 2));

            // We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996) //TODO: look up
            // ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )

            SV0[iBndGP] = LocSV[iBndGP];    // Careful, the SV must always be corrected using SV0 and not LocSV!

            // The following process is adapted from that described by Kaneko et al. (2008)
            LocSR[iBndGP] = std::sqrt(std::pow(LocSR1[iBndGP], 2) + std::pow(LocSR2[iBndGP], 2));
            LocSR[iBndGP] = std::max(AlmostZero, LocSR[iBndGP]);

            SR_tmp[iBndGP] = LocSR[iBndGP];
            invZ = (1.0 / w_speed[2] / rho + 1.0 / w_speed_neig[2] / rho_neig);

            if (ThermalPress == 1) {
                P_f[iBndGP] = TP[iBndGP][iFace][1];
            } else {
                P_f[iBndGP] = 0.0;
            }

        }// End of iBndGP-loop
        for(int j = 1; j < nSVupdates; j++){

            for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
            //fault strength using LocMu and P_f from previous timestep/iteration
            //1.update SV using Vold from the previous time step
            updateStateVariable(nBndGP, RS_f0, RS_b, RS_a[iBndGP]  , RS_sr0, RS_fw, RS_srW[iBndGP], RS_sl0, SV0[iBndGP], time_inc, SR_tmp[iBndGP], LocSV[iBndGP] );
            if(ThermalPress == 1){
                S[iBndGP]  = -LocMu[iBndGP]*(P[iBndGP] - P_f[iBndGP]);

                for (int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++){
                    //!recover original values as it gets overwritten in the ThermalPressure routine
                    Theta_tmp[iTP_grid_nz] = TP_Theta[iBndGP][iFace][iTP_grid_nz];
                    Sigma_tmp[iTP_grid_nz] = TP_sigma[iBndGP][iFace][iTP_grid_nz];
                }
                Calc_ThermalPressure(temp_0, pressure_0, time_inc, TP_grid_nz, TP_half_width_shear_zone[iBndGP][iFace],
                                     alpha_th, alpha_hy[iBndGP][iFace],rho_c, TP_Lambda, Theta_tmp, Sigma_tmp, S[iBndGP], LocSR[iBndGP]  , TP_grid,
                                     TP_DFinv,TP[iBndGP][iFace][0], TP[iBndGP][iFace][1] );
                P_f[iBndGP] = TP[iBndGP][iFace][1];
            }
            //2. solve for Vnew , applying the Newton-Raphson algorithm
            //effective normal stress including initial stresses and pore fluid pressure
            n_stress[iBndGP]   = P[iBndGP] - P_f[iBndGP];

            }// End of iBndGP-loop

            IterativelyInvertSR (nBndGP, nSRupdates, LocSR, RS_sr0, LocSV, RS_a, n_stress, ShTest, invZ, SRtest, has_converged);

            for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {

                // 3. update theta, now using V=(Vnew+Vold)/2
                SR_tmp[iBndGP]=0.5*(LocSR[iBndGP]+abs(SRtest[iBndGP])) ; // For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)

                // 4. solve again for Vnew
                LocSR[iBndGP]=abs(SRtest[iBndGP]);
                //!update LocMu
                tmp[iBndGP]  = 0.5 /RS_sr0 * exp(LocSV[iBndGP]/RS_a[iBndGP]);
                tmp2[iBndGP]  = LocSR[iBndGP]*tmp[iBndGP] ;
                // mu from LocSR
                LocMu[iBndGP]  = RS_a[iBndGP]*log(tmp2[iBndGP] +sqrt(pow(tmp2[iBndGP] ,2)+1.0));
            }// End of iBndGP-loop

        } //End nSVupdates-loop   j=1,nSVupdates   !This loop corrects SV values


        if (!has_converged){
            //!logError(*) 'nonConvergence RS Newton', time
            if (tmp[0] != tmp[0]) {
                //TODO: error logging : logError(*) 'NaN detected', time
                return;
            }
        }

        for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
            updateStateVariable (nBndGP, RS_f0, RS_b, RS_a[iBndGP], RS_sr0, RS_fw, RS_srW[iBndGP], RS_sl0, SV0[iBndGP], time_inc, SR_tmp[iBndGP], LocSV[iBndGP]);

            //! 5. get final theta, mu, traction and slip
            //! SV from mean slip rate in tmp

            if(ThermalPress == 1){
                S[iBndGP] = -LocMu[iBndGP]*(P[iBndGP] - P_f[iBndGP]);

                for (int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++){
                    Theta_tmp[iTP_grid_nz] = TP_Theta[iBndGP][iFace][iTP_grid_nz];
                    Sigma_tmp[iTP_grid_nz] = TP_sigma[iBndGP][iFace][iTP_grid_nz];
                    //!use Theta/Sigma from last call in this update, dt/2 and new SR from NS

                    Calc_ThermalPressure(temp_0, pressure_0, time_inc, TP_grid_nz, TP_half_width_shear_zone[iBndGP][iFace],
                                         alpha_th, alpha_hy[iBndGP][iFace],rho_c, TP_Lambda, Theta_tmp, Sigma_tmp, S[iBndGP], LocSR[iBndGP]  , TP_grid,
                                         TP_DFinv,TP[iBndGP][iFace][0], TP[iBndGP][iFace][1] );

                    P_f[iBndGP] =TP[iBndGP][iFace][1];
                    TP_Theta[iBndGP][iFace][iTP_grid_nz] = Theta_tmp[iTP_grid_nz];
                    TP_sigma[iBndGP][iFace][iTP_grid_nz] = Sigma_tmp[iTP_grid_nz];
                }
            }

            //!update LocMu for next strength determination, only needed for last update
            //! X in Asinh(x) for mu calculation
            tmp[iBndGP] = 0.5/RS_sr0 * exp(LocSV[iBndGP]/RS_a[iBndGP]);
            tmp2[iBndGP] = LocSR[iBndGP]*tmp[iBndGP];
            //! mu from LocSR
            LocMu[iBndGP]  = RS_a[iBndGP]*log(tmp2[iBndGP]+sqrt(pow(tmp2[iBndGP],2)+1.0));

            //! update stress change
            LocTracXY[iBndGP] = -((InitialStressInFaultCS[iBndGP][3][iFace] + XYStressGP[iBndGP][iTimeGP])/ShTest[iBndGP])*LocMu[iBndGP]*(P[iBndGP]-P_f[iBndGP]);
            LocTracXZ[iBndGP] = -((InitialStressInFaultCS[iBndGP][5][iFace]  + XZStressGP[iBndGP][iTimeGP]) /ShTest[iBndGP])*LocMu[iBndGP]*(P[iBndGP]-P_f[iBndGP]);
            LocTracXY[iBndGP] = LocTracXY[iBndGP] - InitialStressInFaultCS[iBndGP][3][iFace];
            LocTracXZ[iBndGP] = LocTracXZ[iBndGP] - InitialStressInFaultCS[iBndGP][5][iFace];

            //Compute slip
            LocSlip[iBndGP]   = LocSlip[iBndGP]  + (LocSR[iBndGP])*time_inc; //! ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening

            //!Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
            LocSR1[iBndGP]     = -invZ*(LocTracXY[iBndGP]-XYStressGP[iBndGP][iTimeGP]);
            LocSR2[iBndGP]     = -invZ*(LocTracXZ[iBndGP]-XZStressGP[iBndGP][iTimeGP]);

            //!TU 07.07.16: correct LocSR1_2 to avoid numerical errors
            tmp[iBndGP]  = sqrt( pow(LocSR1[iBndGP],2)+ pow(LocSR2[iBndGP],2) );
            if( tmp != 0){
                LocSR1[iBndGP] = LocSR[iBndGP]*LocSR1[iBndGP]/tmp[iBndGP];
                LocSR2[iBndGP] = LocSR[iBndGP]*LocSR2[iBndGP]/tmp[iBndGP];
            }

            tmpSlip[iBndGP]  = tmpSlip[iBndGP]  + tmp[iBndGP] *time_inc;

            LocSlip1[iBndGP]   = LocSlip1[iBndGP]  + (LocSR1[iBndGP])*time_inc;
            LocSlip2[iBndGP]   = LocSlip2[iBndGP]  + (LocSR2[iBndGP])*time_inc;
            //LocSR1     = SignSR1*ABS(LocSR1)
            //LocSR2     = SignSR2*ABS(LocSR2)

            //!Save traction for flux computation
            TractionGP_XY[iBndGP][iTimeGP] = LocTracXY[iBndGP];
            TractionGP_XZ[iBndGP][iTimeGP] = LocTracXZ[iBndGP];
        }

    } // End of iTimeGP-loop

    for (int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
        //!output rupture front
        //! outside of iTimeGP loop in order to safe an 'if' in a loop
        //! this way, no subtimestep resolution possible
        if (RF[iBndGP][iFace] && LocSR[iBndGP] > 0.001) {
            rupture_time[iBndGP][iFace] = time;
            RF[iBndGP][iFace] = false;
        }
        if(LocSR[iBndGP] > PeakSR[iBndGP][iFace]){
            PeakSR[iBndGP][iFace] = LocSR[iBndGP];
        }

        //!output time when shear stress is equal to the dynamic stress after rupture arrived
        //!currently only for linear slip weakening
        if( (rupture_time[iBndGP][iFace] > 0.0) &&
            (rupture_time[iBndGP][iFace] <= time) &&
                DS[iBndGP][iFace] &&
                (Mu[iBndGP][iFace] <= (RS_fw+0.05*(RS_f0-RS_fw) ))){
            dynStress_time[iBndGP][iFace] = time;
            DS[iBndGP][iFace] = false;
        }


        Mu[iBndGP][iFace]       = LocMu[iBndGP];
        SlipRate1[iBndGP][iFace]  = LocSR1[iBndGP];
        SlipRate2[iBndGP][iFace]  = LocSR2[iBndGP];
        Slip[iBndGP][iFace]       = LocSlip[iBndGP];
        Slip1[iBndGP][iFace]     = LocSlip1[iBndGP];
        Slip2[iBndGP][iFace]      = LocSlip2[iBndGP];
        TracXY[iBndGP][iFace] = LocTracXY[iBndGP];
        TracXZ[iBndGP][iFace] = LocTracXZ[iBndGP];


        matmul = 0.0;
        for (int j = 0; j < nBndGP; j++) {
            matmul += resampleMatrixView(iBndGP,j) * (LocSV[j] - StateVar[j][iFace]);
        }
        StateVar[iBndGP][iFace]   += matmul;


        if(magnitude_out[iFace] ){
            sum_tmpSlip = 0.0;
            for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++) {
                sum_tmpSlip += tmpSlip[iBndGP];
            }
            averaged_Slip[iFace] = averaged_Slip[iFace] + sum_tmpSlip / nBndGP;
        }

    }
}




/*
 * src/Solver/prak_clif_mod.f90
 *
module prak_clif_mod

contains

        elemental subroutine prakash_cliff_fric(strength,sigma,V,Vstar,L,mu,dt)
implicit none
real, intent(inout):: strength !< shear strength (or friction)
real, intent(in)   :: V        !< slip velocity
real, intent(in)   :: sigma    !< normal traction
real, intent(in)   :: Vstar    !< referenz Velocity
real, intent(in)   :: L        !< referenz length
real, intent(in)   :: mu       !< friction coefficient
real, intent(in)   :: dt       !< time increment
real :: expterm
        expterm=exp(-(abs(V) + Vstar)*dt/L)
strength =   strength*expterm - max(0.,-mu*sigma)*(expterm-1.)
end subroutine prakash_cliff_fric

        end module prak_clif_mod
*/

void seissol::physics::Evaluate_friction_law::prak_clif_mod(double &strength, double &sigma, double &V, double &Vstar, double &L, double &mu, double &dt){
    double expterm;
    expterm = std::exp(-(std::abs(V) + Vstar)*dt/L);
    strength =  strength* expterm - std::max(0.,-mu*sigma)*(expterm-1.);
}

/*
 * Function in NucleationFunctions_mod
 */
double seissol::physics::Evaluate_friction_law::Calc_SmoothStepIncrement(double time, double Tnuc, double dt){
    double Gnuc;
    double prevtime;
    if(time > 0.0 && time <= Tnuc){
        Gnuc = Calc_SmoothStep(time, Tnuc);
        prevtime = time - dt;
        if(prevtime > 0.0){
            Gnuc = Gnuc - Calc_SmoothStep(prevtime, Tnuc);
        }
    }else{
        Gnuc = 0.0;
    }
    return Gnuc;
}

/*
 * Function in NucleationFunctions_mod
 */
double seissol::physics::Evaluate_friction_law::Calc_SmoothStep(double time, double Tnuc){
    double Gnuc;
    if (time <= 0){
        Gnuc=0.0;
    }else{
        if (time < Tnuc){
            Gnuc = std::exp( std::pow(time-Tnuc,2)/(time*(time-2.0*Tnuc)));
        }else{
            Gnuc=1.0;
        }
    }
    return Gnuc;
}

void seissol::physics::Evaluate_friction_law::TrafoChiTau2XiEtaZeta(double &xi, double &eta, double &zeta, double chi, double tau, int iSide, int iNeighborVertex){
    double chi1, tau1;
    switch(iNeighborVertex){
        // Inside the element itself
        case 0:
            switch(iSide){
                case 1:
                    xi   = tau;
                    eta  = chi;
                    zeta = 0.0;
                    break;
                case 2:
                    xi   = chi;
                    eta  = 0.0;
                    zeta = tau;
                    break;
                case 3:
                    xi   = 0.0;
                    eta  = tau;
                    zeta = chi;
                    break;
                case 4:
                    xi   = 1.0-chi-tau;
                    eta  = chi;
                    zeta = tau;
                    break;
            }
            return;
        // Inside the neighbor
        case 1:
            chi1 = tau;
            tau1 = chi;
            break;
        case 2:
            chi1 = 1.-chi-tau;
            tau1 = tau;
            break;
        case 3:
            chi1 = chi;
            tau1 = 1.-chi-tau;
            break;
    }
    switch(iSide){
        case 1:
            xi   = tau1;
            eta  = chi1;
            zeta = 0.0;
            break;
        case 2:
            xi   = chi1;
            eta  = 0.0;
            zeta = tau1;
            break;
        case 3:
            xi   = 0.0;
            eta  = tau1;
            zeta = chi1;
            break;
        case 4:
            xi   = 1.-chi1-tau1;
            eta  = chi1;
            zeta = tau1;
            break;
    }
}

/*
 * Shifted indexing by 1 towards zero
 */
void seissol::physics::Evaluate_friction_law::TetraTrafoXiEtaZeta2XYZ(double &xP, double &yP, double &zP, double xi, double eta, double zeta,double x[4], double y[4], double z[4]){
    xP = x[0] + (x[1]-x[0])*xi + (x[2]-x[0])*eta + (x[3]-x[0])*zeta;
    yP = y[0] + (y[1]-y[0])*xi + (y[2]-y[0])*eta + (y[3]-y[0])*zeta;
    zP = z[0] + (z[1]-z[0])*xi + (z[2]-z[0])*eta + (z[3]-z[0])*zeta;
}

void seissol::physics::Evaluate_friction_law::updateStateVariable(int nBndGP, double RS_f0, double RS_b, double RS_a, double RS_sr0,
        double RS_fw, double RS_srW, double RS_sl0, double SV0, double time_inc, double SR_tmp, double &LocSV){
    double flv, fss, SVss;

    // low-velocity steady state friction coefficient
    flv = RS_f0 - (RS_b-RS_a)* log(SR_tmp/RS_sr0);
    // steady state friction coefficient
    fss = RS_fw + (flv - RS_fw)/pow(1.0+pow(SR_tmp/RS_srW,8) ,1.0/8.0);
    // steady-state state variable
    // For compiling reasons we write SINH(X)=(EXP(X)-EXP(-X))/2
    SVss = RS_a * log(2.0*RS_sr0/SR_tmp * (exp(fss/RS_a)-exp(-fss/RS_a))/2.0);

    // exact integration of dSV/dt DGL, assuming constant V over integration step
    LocSV = SVss*(1.0-exp(-SR_tmp*time_inc/RS_sl0))+exp(-SR_tmp*time_inc/RS_sl0)*SV0;

    /*  //TODO log error NaN detected
    if (ANY(IsNaN(LocSV)) == true){
        logError(*) 'NaN detected'
    }
     */
}

void seissol::physics::Evaluate_friction_law::Calc_ThermalPressure(double temp_0, double pressure_0, double dt,
        int TP_grid_nz, double TP_half_width_shear_zone, double alpha_th, double alpha_hy, double rho_c,
        double Lambda, double theta[], double sigma[], double Sh, double SR,
        double Dwn[], double DFinv[], double temp, double pressure){

    double tauV,Lambda_prime, tmp, theta_current, sigma_current, T, p;

    double omega;

    T = 0.0;
    p = 0.0;

    for (int iTP_grid_nz = 0; iTP_grid_nz < TP_grid_nz; iTP_grid_nz++){
        tauV = Sh*SR; //!fault strenght*slip rate
        Lambda_prime = Lambda*alpha_th/(alpha_hy-alpha_th);
        tmp = pow(Dwn[iTP_grid_nz]/TP_half_width_shear_zone,2);
        //!1. Calculate diffusion of the field at previous timestep

        //!temperature
        theta_current = theta[iTP_grid_nz]*exp(-alpha_th*dt*tmp);
        //!pore pressure + lambda'*temp
        sigma_current = sigma[iTP_grid_nz]*exp(-alpha_hy*dt*tmp);

        //!2. Add current contribution and get new temperature
        heat_source(TP_half_width_shear_zone,alpha_th,dt,Dwn[iTP_grid_nz],TP_grid_nz,omega);
        theta[iTP_grid_nz] = theta_current + (tauV/rho_c)*omega;
        heat_source(TP_half_width_shear_zone,alpha_hy,dt,Dwn[iTP_grid_nz],TP_grid_nz,omega);
        sigma[iTP_grid_nz] = sigma_current + ((Lambda+Lambda_prime)*tauV)/(rho_c)*omega;

        //!3. Recover temperature and pressure using inverse Fourier
        //! transformation with the calculated fourier coefficients


        //!new contribution
        T += (DFinv[iTP_grid_nz]/TP_half_width_shear_zone)*theta[iTP_grid_nz] ;   //TODO redo loop
        p += (DFinv[iTP_grid_nz]/TP_half_width_shear_zone)*sigma[iTP_grid_nz];
    }

    //Update pore pressure change (sigma = pore pressure + lambda'*temp)
    //In the BIEM code (Lapusta) they use T without initial value
    p = p - Lambda_prime*T;

    //Temp and pore pressure change at single GP on the fault + initial values
    temp = T + temp_0;
    pressure = -p + pressure_0;
}

void seissol::physics::Evaluate_friction_law::heat_source(double TP_half_width_shear_zone, double alpha, double dt, double Dwn, int TP_grid_nz, double &omega){
    double pi=3.141592653589793;
    double tmp;
    //!Gaussian shear zone in spectral domain, normalized by w
    tmp = pow(Dwn/TP_half_width_shear_zone,2);
    //!original function in spatial domain
    //!omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_half_width_shear_zone).^2);
    //!function in the wavenumber domain *including additional factors in front of the heat source function*
    //!omega = 1/(*alpha*Dwn**2**(sqrt(2.0*pi))*exp(-0.5*(Dwn*TP_half_width_shear_zone)**2)*(1-exp(-alpha**dt**tmp))
    //!inserting Dwn/TP_half_width_shear_zone (scaled) for Dwn cancels out TP_half_width_shear_zone
    omega = 1.0/(alpha*tmp*(sqrt(2.0*pi))) * exp(-0.5*pow(Dwn,2)) * (1.0 - exp(-alpha*dt*tmp)) ;
}

void seissol::physics::Evaluate_friction_law::IterativelyInvertSR (int nBndGP, int nSRupdates, double LocSR[], double RS_sr0,
        double LocSV[], double RS_a[], double n_stress[], double sh_stress[], double invZ, double SRtest[], bool &has_converged){

    double tmp, tmp2, tmp3, mu_f, dmu_f, NR[nBndGP], dNR;
    double aTolF = 1e-8;
    double AlmostZero = 1e-45;


    //!solve for Vnew = SR , applying the Newton-Raphson algorithm
    //!SR fulfills g(SR)=f(SR)
    //!-> find root of NR=f-g using a Newton-Raphson algorithm with dNR = d(NR)/d(SR)
    //!SR_{i+1}=SR_i-( NR_i / dNR_i )
    //!
    //        !equalize:
    //!         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
    //!         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
    //!  where mu = friction coefficient, dependening on the RSF law used

    for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++){
        //! first guess = SR value of the previous step
        SRtest[iBndGP] = LocSR[iBndGP];
        tmp   =  0.5 / RS_sr0 *exp(LocSV[iBndGP]/RS_a[iBndGP]);
    }



    has_converged = false;

    for(int i = 0; i < nSRupdates; i++){
        for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++){


            //!f = ( tmp2 * ABS(LocP+P_0)- ABS(S_0))*(S_0)/ABS(S_0)
            //!g = SRtest * 1.0/(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) + ABS(ShTest)
            //!for compiling reasons ASINH(X)=LOG(X+SQRT(X^2+1))

            //!calculate friction coefficient
            tmp2  = tmp*SRtest[iBndGP];
            mu_f  = RS_a[iBndGP] * log(tmp2+sqrt(pow(tmp2,2)+1.0));
            dmu_f = RS_a[iBndGP] / sqrt(1.0+pow(tmp2,2))*tmp;
            NR[iBndGP]    = -invZ * (abs(n_stress[iBndGP])*mu_f-sh_stress[iBndGP])-SRtest[iBndGP];
        }
        has_converged = true;
        for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++){
            if (abs(NR[iBndGP]) >= aTolF ){
                has_converged = false;
                break;
            }
        }
        if(has_converged){
            return;
        }
        for(int iBndGP = 0; iBndGP < nBndGP; iBndGP++){
            //!derivative of NR
            dNR   = -invZ * (abs(n_stress[iBndGP])*dmu_f) - 1.0;
            //!ratio
            tmp3 = NR[iBndGP]/dNR;

            //!update SRtest
            SRtest[iBndGP] = std::max(AlmostZero,SRtest[iBndGP]-tmp3);
        }
    }
}