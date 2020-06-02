//
// Created by adrian on 14.05.20.
//

#ifndef SEISSOL_EVALUATE_FRICTION_LAW_H
#define SEISSOL_EVALUATE_FRICTION_LAW_H

#endif //SEISSOL_EVALUATE_FRICTION_LAW_H

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/LTS.h>
#include <Initializer/DynamicRupture.h>
#include <Initializer/tree/Layer.hpp>
#include <Checkpoint/Fault.h>


namespace seissol {
    namespace physics {
        class Evaluate_friction_law;
    }
}
class seissol::physics::Evaluate_friction_law {
private:
    //TODO: Interface Eval_friction_law??
    float const u_0 = 10e-14; //slip rate is considered as being zero for instaneous healing
    double const ZERO = 0.0;; // CHECK: 0.0D0 -> double precision exponent

    //Methods:
    /*
    PRIVATE :: updateStateVariable
    PRIVATE :: IterativelyInvertSR

    PRIVATE :: ImposedSlipRateOnDRBoundary
    PRIVATE :: rate_and_state
    PRIVATE :: rate_and_state_vw
    PRIVATE :: rate_and_state_nuc101
    PRIVATE :: rate_and_state_nuc103
    */

    /*
     * actually from another file: src/Solver/prak_clif_mod.f90
     */
    void prak_clif_mod(double &strength, double &sigma, double &V, double &Vstar, double &L, double &mu, double &dt);

    /*
     * Function in NucleationFunctions_mod
     */
    double Calc_SmoothStepIncrement(double time, double Tnuc, double dt);

    /*
     * Function in NucleationFunctions_mod
     */
    double Calc_SmoothStep(double time, double Tnuc);

    //PRIVATE :: no_fault
    //PURE SUBROUTINE no_fault(TractionGP_XY,TractionGP_XZ,XYStressGP,XZStressGP)
    void no_fault(
            double **XYStressGP, double **XZStressGP,                            //INPUT
            double **TractionGP_XY,                                             //OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
    );
    /*
     * from src/Numerical_aux/dgbasis.f90 - TrafoChiTau2XiEtaZeta
     * aps side-local chi-tau coordinates to reference xi-eta-zeta coordinates
     */
    void TrafoChiTau2XiEtaZeta(double &xi, double &eta, double &zeta, double chi, double tau, int iSide, int iNeighborVertex);

    /*
     * from src/Numerical_aux/dgbasis.f90
     * maps xi-eta-zeta coords. to xyz coords. in a tetrahedron
     */
    void TetraTrafoXiEtaZeta2XYZ(double &xP, double &yP, double &zP, double xi, double eta, double zeta,double x[4], double y[4], double z[4]);

    /*
     * Method from Evaluate_friction_law.f90
     */
    void updateStateVariable(int nBndGP, double RS_f0, double RS_b, double RS_a, double RS_sr0, double RS_fw, double RS_srW,
            double RS_sl0, double SV0, double time_inc, double SR_tmp, double &LocSV);

    /*
     * form src/Physics/thermalpressure.f90
     *
     * theta[TP_grid_nz], sigma[TP_grid_nz], Dwn[TP_grid_nz]
     */
    void Calc_ThermalPressure(double temp_0, double pressure_0, double dt,
                              int TP_grid_nz, double TP_half_width_shear_zone, double alpha_th, double alpha_hy, double rho_c,
                              double Lambda, double theta[], double sigma[], double Sh, double SR,
                              double Dwn[], double DFinv[], double temp, double pressure);

    /*
     * form src/Physics/thermalpressure.f90
     */
    void heat_source(double TP_half_width_shear_zone, double alpha, double dt, double Dwn, int TP_grid_nz, double &omega);

    /*
     * Method from Evaluate_friction_law.f90
     * LocSR[nBndGP], LocSV[nBndGP], RS_a[nBndGP], n_stress[nBndGP], sh_stress[nBndGP], SRtest[nBndGP]
     */
    void IterativelyInvertSR (int nBndGP, int nSRupdates, double LocSR[], double RS_sr0,
            double LocSV[], double RS_a[], double n_stress[], double sh_stress[],
            double invZ, double SRtest[], bool &has_converged);


    /*
     * Special friction case 6: linear slip weakening with Prakash-Clifton regularization
     */
    void Linear_slip_weakening_bimaterial(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,
            real const **resampleMatrix,                                         //
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO
            //initializers::LTSTree* io_ltsTree, initializers::LTS*  i_lts, initializers::LTSTree* dynRupTree,  initializers::DynamicRupture* dynRup,             //data structs
            //checkpoint::Fault fault
    );

    /*
     *     !> friction case 16,17
     *     !> Specific conditions for SCEC TPV16/17
     *     !> basically, introduction of a time dependent forced rupture
     */
    void Linear_slip_weakening_TPV1617(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,
            real const **resampleMatrix,                                         //
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO
    );

    /*
     * !< T. Ulrich 27.07.17
     * !< This friction law allows imposing a slip rate on the DR boundary
     */
    void ImposedSlipRateOnDRBoundary(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO
    );

    /*
     *     !> friction case 3,4: rate and state friction
     *     !> aging (3) and slip law (4)
     */
    void rate_and_state(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO
    );


    /*
     *     !> friction case 7: severe velocity weakening rate and state friction
     *     !< after Ampuero and Ben-Zion 2008
     */
    void rate_and_state_vw(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO
    );

    /*
     *     !> special friction case for SCEC TPV101: rate and state friction
     *     !> aging law
     *     !< with time and space dependent nucleation
     */
    void rate_and_state_nuc101(
            double **TractionGP_XY, double **TractionGP_XZ,                     // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT, double iT[6][6],                           //IN: time, inv Trafo
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO, void *BND
    );



    /*
     *     !> special friction case for SCEC TPV103: rate and state friction
     *     !> slip law with strong weakening
     *     !< with time and space dependent nucleation
     */
    void rate_and_state_nuc103(
            double **TractionGP_XY, double **TractionGP_XZ,                     // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int iFace, int iSide, int iElem, int nBndGP, int nTimeGP,           // IN: element ID, nBndGP = Nr of boundary Gausspoints, nTimeGP = Nr of time Gausspoints
            double rho, double rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            double time, double *DeltaT,                                        //IN: time
            real const **resampleMatrix,
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO, void *BND
    );



public:
    /*
     * output:
     *
     * updated Traction:
     *  TractionGP_XY, TractionGP_XZ: 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
     *
     * input:
     *
     * Godunov status:
     *  NorStressGP, XYStressGP, XZStressGP: 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
     *
     * element ID, time, inv Trafo:
     *  int iFace = int i_face                              (Interoperability)
     *  int iSide = l_domain%MESH%Fault%Face(i_face,2,1)    (f_ctof_bind_interoperability)
     *      with l_domain = void* i_domain                  (Interoperability)
     *  int iElem = l_domain%MESH%Fault%Face(i_face,1,1)    (f_ctof_bind_interoperability) Remark: iElem denotes "+" side
     *      with l_domain = void* i_domain                  (Interoperability)
     *  double* time = l_time = i_time                      (f_ctof_bind_interoperability)
     *      with i_time = double* i_fullUpdateTime          (Interoperability)
     *  double *timePoints = double* timePoints             (Interoperability)
     *      with  size: (CONVERGENCE_ORDER)                 (f_ctof_bind_interoperability)
     *  background values:
     *  double rho = double densityPlus                     (Interoperability)
     *  double rho_neig = double densityMinus               (Interoperability)
     *  double *w_speed = [pWaveVelocityPlus, sWaveVelocityPlus, sWaveVelocityPlus]             (f_ctof_bind_interoperability)
     *      with array size = [3]
     *      with double pWaveVelocityPlus                                                       (Interoperability)
     *      with double sWaveVelocityPlus                                                       (Interoperability)
     *  double *w_speed_neig = [pWaveVelocityMinus, sWaveVelocityMinus, sWaveVelocityMinus]     (f_ctof_bind_interoperability)
     *      with array size = [3]
     *      with double pWaveVelocityMinus                                                      (Interoperability)
     *      with double sWaveVelocityMinus                                                      (Interoperability)
     *
     *  global variables:
     *  EQN = l_domain%eqn
     *      with type tEquations                        (/Numerical_aux/typesdef.f90)
     *      EQN%InitialStressInFaultCS
     *      EQN%FL                      actually not really needed bc of restructuring the the code
     *      EQN%IniBulk_xx(:,iFace)
     *      EQN%IniBulk_yy(:,iFace)
     *      EQN%IniBulk_zz(:,iFace)
     *      EQN%IniShearXY(:,iFace)
     *      EQN%IniShearYZ(:,iFace)
     *      EQN%IniShearXZ(:,iFace)
     *      EQN%NucleationStressInFaultCS
     *  DISC = l_domain%disc
     *      with type tDiscretization                   (/Numerical_aux/typesdef.f90)
     *  MESH = l_domain%mesh
     *      with type tUnstructMesh                     (/Numerical_aux/typesdef.f90)
     *  MPI =l_domain%mpi
     *      with type tMPI                              (/Numerical_aux/typesdef.f90)
     *  IO = l_domain%io
     *      with type tInputOutput                      (/Numerical_aux/typesdef.f90)
     *  BND = l_domain%bnd
     *      with l_domain = tUnstructDomainDescript     (/Numerical_aux/typesdef.f90)
     *
     * depending on the EQN.FL setting this function calls different friction models
     *
     */
     void Eval_friction_law(
            double **TractionGP_XY,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **TractionGP_XZ,                                              // OUT: updated Traction 2D array with size [1:i_numberOfPoints, CONVERGENCE_ORDER]
            double **NorStressGP, double **XYStressGP, double **XZStressGP,        // IN: Godunov status
            int &iFace, int &iSide, int &iElem, double &time, double *timePoints,  // IN: element ID, time, inv Trafo
            double &rho, double &rho_neig, double *w_speed, double *w_speed_neig, // IN: background values
            real const **resampleMatrix,                                         //
            void *EQN, void *DISC, void *MESH, void *MPI, void *IO, void *BND                                           //data structs
    );


};
