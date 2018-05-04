/*!============================================================
  |                                                           |
  |         micromorphic_linear_elasticity_voigt.cpp          |
  |                                                           |
  -------------------------------------------------------------
  | The source file for the definition of a                   |
  | micromorphic linear elasticity using voigt notation.      |
  -------------------------------------------------------------
  | Notes: Micromorphic constitutive models should be         |
  |        developed in the namespace micro_material          |
  |        and have the function get_stress. This             |
  |        function should read in the right Cauchy           |
  |        green deformation tensor, Psi, Gamma, and          |
  |        write the PK2 stress, the symmetric stress         |
  |        in the reference configuration, and the            |
  |        higher order couple stress in the reference        |
  |        configuration. (ADDITIONAL VALUES WILL BE          |
  |        ADDED OVER TIME).                                  |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#include <iostream>  
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <deformation_measures.h>
#include <micromorphic_linear_elasticity_voigt.h>

namespace micro_material{

    void LinearElasticity::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                          const double (&grad_u)[3][3],           const double (&phi)[9],
                                          const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                          const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                          Vector_9 &cauchy, Vector_9 &s, Vector_27 &m, std::vector<Eigen::VectorXd> &ADD_TERMS){
        /*!
        =======================
        |    evaluate_model   |
        =======================
        
        Evaluate the constitutive model 
        from the general incoming values.
        
        Only returns stresses and additional 
        terms.
        
        */
        
        //Extract the time
        double t  = time[0];
        double dt = time[1];

        //Extract the parameters        
        double params[18];
        if(fparams.size() == 18){
            for(int i=0; i<18; i++){
                params[i] = fparams[i];
            }
        }
        else{std::cout << "Error: Material parameters incorrectly specified\n";}
        
        //Compute the required deformation measures
        Matrix_3x3 F;
        Matrix_3x3 chi;
        Matrix_3x9 grad_chi;
        get_deformation_measures(grad_u, phi, grad_phi, F, chi, grad_chi);
        
        //Compute the stresses
        Vector_9  PK2;
        Vector_9  SIGMA;
        Vector_27 M;
        
        get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M);
        deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);
        
        return;
    }
                                
    void LinearElasticity::evaluate_model(const std::vector<double> &time,        const std::vector<double> (&fparams),
                                          const double (&grad_u)[3][3],           const double (&phi)[9],
                                          const double (&grad_phi)[9][3],         std::vector<double> &SDVS,
                                          const std::vector<double> &ADD_DOF,     const std::vector<std::vector<double>> &ADD_grad_DOF,
                                          Vector_9    &cauchy,    Vector_9    &s,           Vector_27    &m,
                                          Matrix_9x9  &DcauchyDgrad_u, Matrix_9x9  &DcauchyDphi, Matrix_9x27  &DcauchyDgrad_phi,
                                          Matrix_9x9  &DsDgrad_u,      Matrix_9x9  &DsDphi,      Matrix_9x27  &DsDgrad_phi,
                                          Matrix_27x9 &DmDgrad_u,      Matrix_27x9 &DmDphi,      Matrix_27x27 &DmDgrad_phi,
                                          std::vector<Eigen::VectorXd> &ADD_TERMS,               std::vector<Eigen::MatrixXd> &ADD_JACOBIANS){
        /*!
        ========================
        |    evaluate_model    |
        ========================
        
        Evaluate the constitutive model 
        from the general incoming values.
        
        Returns stresses, additional 
        terms, and their jacobians.
        
        */

        //Extract the time
        double t  = time[0];
        double dt = time[1];

        //Extract the parameters        
        double params[18];
        if(fparams.size() == 18){
            for(int i=0; i<18; i++){
                params[i] = fparams[i];
            }
        }
        else{std::cout << "Error: Material parameters incorrectly specified\n";assert(-21==-20);}

        //Compute the required deformation measures
        Matrix_3x3 F;
        Matrix_3x3 chi;
        Matrix_3x9 grad_chi;
        get_deformation_measures(grad_u, phi, grad_phi, F, chi, grad_chi);
        
        //Compute the stresses and jacobians
        Vector_9  PK2;
        Vector_9  SIGMA;
        Vector_27 M;
        
        Matrix_9x9   dPK2dF;
        Matrix_9x9   dPK2dchi;
        Matrix_9x27  dPK2dgrad_chi;
        Matrix_9x9   dSIGMAdF;
        Matrix_9x9   dSIGMAdchi;
        Matrix_9x27  dSIGMAdgrad_chi;
        Matrix_27x9  dMdF;
        Matrix_27x9  dMdchi;
        Matrix_27x27 dMdgrad_chi;
        
        //Note: D(x)Dphi = d(x)dchi
        Matrix_9x9   dcauchydF;
        Matrix_9x27  dcauchydgrad_chi;
        Matrix_9x9   dsdF;
        Matrix_9x27  dsdgrad_chi;
        Matrix_27x9  dmdF;
        Matrix_27x27 dmdgrad_chi;

        //assert(13==14);

        get_stress(t, dt, params, F, chi, grad_chi, SDVS, PK2, SIGMA, M,
                   dPK2dF,   dPK2dchi,   dPK2dgrad_chi,
                   dSIGMAdF, dSIGMAdchi, dSIGMAdgrad_chi,
                   dMdF,     dMdchi,     dMdgrad_chi);
                   
        //assert(14==15);
        deformation_measures::map_stresses_to_current_configuration(F, chi, PK2, SIGMA, M, cauchy, s, m);

        //assert(15==16);
        
        deformation_measures::map_jacobians_to_current_configuration(F,      chi,      PK2,           SIGMA,     M,           cauchy, s,   m,
                                                                     dPK2dF, dPK2dchi, dPK2dgrad_chi, dSIGMAdF,  dSIGMAdchi,  dSIGMAdgrad_chi,
                                                                     dMdF,   dMdchi,   dMdgrad_chi,   dcauchydF, DcauchyDphi, dcauchydgrad_chi,
                                                                     dsdF,   DsDphi,   dsdgrad_chi,   dmdF,      DmDphi,      dmdgrad_chi);
        //assert(16==17);                                                                    
        Matrix_3x9 _grad_phi;
        Vector_27  _grad_phi_v;
        Matrix_3x3 eye = Matrix_3x3::Identity();
        deformation_measures::assemble_grad_chi(grad_phi, eye, _grad_phi); //Put grad_phi into an Eigen Matrix
        deformation_measures::voigt_3x9_tensor(_grad_phi,_grad_phi_v);     //Put grad_phi into voigt notation
        deformation_measures::compute_total_derivatives(F, _grad_phi_v,
                                                        dcauchydF,      dcauchydgrad_chi, dsdF,      dsdgrad_chi, dmdF,      dmdgrad_chi,
                                                        DcauchyDgrad_u, DcauchyDgrad_phi, DsDgrad_u, DsDgrad_phi, DmDgrad_u, DmDgrad_phi);
        //assert(17==18);
        return;
    }
    
    void LinearElasticity::get_deformation_measures(const double (&grad_u)[3][3], const double (&phi)[9], const double (&grad_phi)[9][3],
                                                    Matrix_3x3 &F,                Matrix_3x3 &chi,        Matrix_3x9 &grad_chi){
        /*!
        ==================================
        |    get_deformation_measures    |
        ==================================
        
        Compute the deformation measures from the degrees of freedom and their 
        gradients.
        
        */
        
        deformation_measures::get_deformation_gradient(grad_u, F);

        deformation_measures::assemble_chi(phi, chi);

        deformation_measures::assemble_grad_chi(grad_phi, F, grad_chi);
        
        return;
    }

    void get_stress(const double &t,     const double &dt,      const double (&params)[18],
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS,    Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M){
    
        //Extract the parameters
        double lambda  = params[ 0];
        double mu      = params[ 1];
        double eta     = params[ 2];
        double tau     = params[ 3];
        double kappa   = params[ 4];
        double nu      = params[ 5];
        double sigma   = params[ 6];
        double tau1    = params[ 7];
        double tau2    = params[ 8];
        double tau3    = params[ 9];
        double tau4    = params[10];
        double tau5    = params[11];
        double tau6    = params[12];
        double tau7    = params[13];
        double tau8    = params[14];
        double tau9    = params[15];
        double tau10   = params[16];
        double tau11   = params[17];

        //Initialize the stiffness matrices
        //SpMat A( 9, 9);
        //SpMat B( 9, 9);
        //SpMat C(27,27);
        //SpMat D( 9, 9);

        Matrix_9x9   A;
        Matrix_9x9   B;
        Matrix_27x27 C;
        Matrix_9x9   D;

        //Populate the stiffness matrices
        compute_A_voigt(lambda,mu,A);
        compute_B_voigt(eta,kappa,nu,sigma,tau,B);
        compute_C_voigt(tau1,tau2,tau3, tau4, tau5,tau6,
                        tau7,tau8,tau9,tau10,tau11,C);
        compute_D_voigt(sigma,tau,D);

        //Compute the deformation measures
        Matrix_3x3 RCG; //The right cauchy green deformation tensor
        deformation_measures::get_right_cauchy_green(F,RCG);
        Matrix_3x3 RCGinv = RCG.inverse(); //The inverse of the right cauchy green deformation tensor
        Matrix_3x3 Psi; //The second order micromorphic deformation measure
        deformation_measures::get_psi(F,chi,Psi);
        Matrix_3x9 Gamma; //The higher order micromorphic deformation measure
        deformation_measures::get_gamma(F,grad_chi,Gamma);

        //Compute the strain measures
        Matrix_3x3 E;
        Matrix_3x3 E_micro;
        deformation_measures::get_lagrange_strain(F,E);
        deformation_measures::get_micro_strain(Psi,E_micro);
        
        //std::cout << "F:\n" << F << "\n";
        //std::cout << "chi:\n" << chi << "\n";
        //std::cout << "grad chi:\n" << grad_chi << "\n";
        //std::cout << "Psi:\n" << Psi << "\n";
        //std::cout << "Gamma:\n" << Gamma << "\n";
        //std::cout << "E:\n" << E << "\n";
        //std::cout << "E_micro:\n" << E_micro << "\n";

        //Put the strain measures in voigt notation
        Vector_9  E_voigt;
        Vector_9  E_micro_voigt;
        Vector_27 Gamma_voigt;

        deformation_measures::voigt_3x3_tensor(E,       E_voigt);
        deformation_measures::voigt_3x3_tensor(E_micro, E_micro_voigt);
        deformation_measures::voigt_3x9_tensor(Gamma,   Gamma_voigt);

        //Compute the stress measures
        compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      PK2);
                           
        compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      SIGMA);
                           
        compute_higher_order_stress(Gamma_voigt, C, M);

        return;
    }
    
    void get_stress(const double &t,     const double &dt,      const double (&params)[18],
                    const Matrix_3x3 &F, const Matrix_3x3 &chi, const Matrix_3x9 &grad_chi,
                    std::vector<double> &SDVS,    Vector_9 &PK2, Vector_9 &SIGMA, Vector_27 &M,
                    Matrix_9x9  &dPK2dF,   Matrix_9x9  &dPK2dchi,   Matrix_9x27  &dPK2dgrad_chi,
                    Matrix_9x9  &dSIGMAdF, Matrix_9x9  &dSIGMAdchi, Matrix_9x27  &dSIGMAdgrad_chi,
                    Matrix_27x9 &dMdF,     Matrix_27x9 &dMdchi,     Matrix_27x27 &dMdgrad_chi){
        /*!=================
        |    get_stress    |
        ====================
        
        Computes the stress measures and their jacobians.
        
        */
        //Extract the parameters
        double lambda  = params[ 0];
        double mu      = params[ 1];
        double eta     = params[ 2];
        double tau     = params[ 3];
        double kappa   = params[ 4];
        double nu      = params[ 5];
        double sigma   = params[ 6];
        double tau1    = params[ 7];
        double tau2    = params[ 8];
        double tau3    = params[ 9];
        double tau4    = params[10];
        double tau5    = params[11];
        double tau6    = params[12];
        double tau7    = params[13];
        double tau8    = params[14];
        double tau9    = params[15];
        double tau10   = params[16];
        double tau11   = params[17];

        //Initialize the stiffness matrices
        //SpMat A( 9, 9);
        //SpMat B( 9, 9);
        //SpMat C(27,27);
        //SpMat D( 9, 9);

        Matrix_9x9   A;
        Matrix_9x9   B;
        Matrix_27x27 C;
        Matrix_9x9   D;

        //assert(101==102);
        //Populate the stiffness matrices
        compute_A_voigt(lambda,mu,A);
        compute_B_voigt(eta,kappa,nu,sigma,tau,B);
        compute_C_voigt(tau1,tau2,tau3, tau4, tau5,tau6,
                        tau7,tau8,tau9,tau10,tau11,C);
        compute_D_voigt(sigma,tau,D);

        //Compute the deformation measures
        Matrix_3x3 RCG; //The right cauchy green deformation tensor
        deformation_measures::get_right_cauchy_green(F,RCG);
        Matrix_3x3 RCGinv = RCG.inverse(); //The inverse of the right cauchy green deformation tensor
        Matrix_3x3 Psi; //The second order micromorphic deformation measure
        deformation_measures::get_psi(F,chi,Psi);
        Matrix_3x9 Gamma; //The higher order micromorphic deformation measure
        deformation_measures::get_gamma(F,grad_chi,Gamma);

        //Compute the strain measures
        Matrix_3x3 E;
        Matrix_3x3 E_micro;
        deformation_measures::get_lagrange_strain(F,E);
        deformation_measures::get_micro_strain(Psi,E_micro);

        //Put the strain measures in voigt notation
        Vector_9  E_voigt;
        Vector_9  E_micro_voigt;
        Vector_27 Gamma_voigt;

        deformation_measures::voigt_3x3_tensor(E,       E_voigt);
        deformation_measures::voigt_3x3_tensor(E_micro, E_micro_voigt);
        deformation_measures::voigt_3x9_tensor(Gamma,   Gamma_voigt);

        //Compute the stress measures
        compute_PK2_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      PK2);
                           
        compute_symmetric_stress(E_voigt, E_micro_voigt, Gamma_voigt, RCGinv, Psi, Gamma,
                           A,       B,             C,           D,      SIGMA);
                           
        compute_higher_order_stress(Gamma_voigt, C, M);
        
        //Compute the jacobians w.r.t. the derived deformation measures
        Matrix_9x9  dPK2dRCG;
        Matrix_9x9  dPK2dPsi;
        Matrix_9x27 dPK2dGamma;
        
        Matrix_9x9  dSIGMAdRCG;
        Matrix_9x9  dSIGMAdPsi;
        Matrix_9x27 dSIGMAdGamma;
        
        Matrix_27x27 dMdGamma; //Note: other gradients are zero for M in this form.
        
        //Terms to speed up computation.
        Matrix_9x9  dPK2dRCGterms[4];
        Matrix_9x9  dPK2dPsiterms[3];
        Matrix_9x27 dPK2dGammaterms[2];
        
        //Gradients of the PK2 stress
        compute_dPK2dRCG(RCG, RCGinv, Gamma, Gamma_voigt, E,             E_micro, E_voigt, E_micro_voigt,
                         A,   B,      C,     D,           dPK2dRCGterms, dPK2dRCG);
                        
        compute_dPK2dPsi(RCGinv, E_micro, E_voigt,       E_micro_voigt,
                          B,      D,       dPK2dPsiterms, dPK2dPsi);
                          
        compute_dPK2dGamma(RCGinv, Gamma,           Gamma_voigt,
                           C,      dPK2dGammaterms, dPK2dGamma);
        
        //Gradients of the symmetric stress (reference configuration)
        compute_dSIGMAdRCG(dPK2dRCGterms, dSIGMAdRCG);
    
        compute_dSIGMAdPsi(dPK2dPsiterms, dSIGMAdPsi);
    
        compute_dSIGMAdGamma(dPK2dGammaterms, dSIGMAdGamma);
        
        compute_dMdGamma(Matrix_27x27(C), dMdGamma);
        
        //Gradients of the derived measures
        //Note: Replaced sparse matricies with dense matrices
        //      This is not the most efficient but it seems required
        //      for use in MOOSE.
        Matrix_9x9   dRCGdF;
        
        Matrix_9x9   dPsidF;
        Matrix_9x9   dPsidchi;
        
        Matrix_27x9  dGammadF;
        Matrix_27x27 dGammadgrad_chi;

        deformation_measures::compute_dRCGdF(F,dRCGdF);

        deformation_measures::compute_dPsidF(chi,dPsidF);
        deformation_measures::compute_dPsidchi(F,dPsidchi);
        
        Vector_27 grad_chi_voigt;
        deformation_measures::voigt_3x9_tensor(grad_chi,grad_chi_voigt);
        
        deformation_measures::compute_dGammadF(grad_chi_voigt,dGammadF);
        deformation_measures::compute_dGammadgrad_chi(F, dGammadgrad_chi);

        //Compute the jacobians of the stresses w.r.t. the fundamental deformation measures.
        dPK2dF   = dPK2dRCG*dRCGdF   + dPK2dPsi*dPsidF   + dPK2dGamma*dGammadF;
        dSIGMAdF = dSIGMAdRCG*dRCGdF + dSIGMAdPsi*dPsidF + dSIGMAdGamma*dGammadF;
        dMdF     = dMdGamma*dGammadF; //Note: all other derivatives are zero.
        
        dPK2dchi   = dPK2dPsi*dPsidchi;
        dSIGMAdchi = dSIGMAdPsi*dPsidchi;
        dMdchi     = Matrix_27x9::Zero(); //Note: M is independent of the magnitude of chi (not so for grad_chi)
        
        dPK2dgrad_chi   = dPK2dGamma*dGammadgrad_chi;
        dSIGMAdgrad_chi = dSIGMAdGamma*dGammadgrad_chi;
        dMdgrad_chi     = dMdGamma*dGammadgrad_chi;
        
        return;    
    }

    void compute_A_voigt(const double &lambda,const double &mu, SpMat &A) {
        /*!=========================
           |    compute_A_voigt    |
           =========================
           
           Compute the A stiffness matrix in voigt notation.
        */

        std::vector<T> tripletList;
        tripletList.reserve(21);

        tripletList.push_back(T(0,0,lambda + 2*mu));
        tripletList.push_back(T(0,1,lambda));
        tripletList.push_back(T(0,2,lambda));
        tripletList.push_back(T(1,0,lambda));
        tripletList.push_back(T(1,1,lambda + 2*mu));
        tripletList.push_back(T(1,2,lambda));
        tripletList.push_back(T(2,0,lambda));
        tripletList.push_back(T(2,1,lambda));
        tripletList.push_back(T(2,2,lambda + 2*mu));
        tripletList.push_back(T(3,3,mu));
        tripletList.push_back(T(3,6,mu));
        tripletList.push_back(T(4,4,mu));
        tripletList.push_back(T(4,7,mu));
        tripletList.push_back(T(5,5,mu));
        tripletList.push_back(T(5,8,mu));
        tripletList.push_back(T(6,3,mu));
        tripletList.push_back(T(6,6,mu));
        tripletList.push_back(T(7,4,mu));
        tripletList.push_back(T(7,7,mu));
        tripletList.push_back(T(8,5,mu));
        tripletList.push_back(T(8,8,mu));

        A.setFromTriplets(tripletList.begin(), tripletList.end());

        return;
    }

    void compute_A_voigt(const double &lambda,const double &mu, Matrix_9x9 &A) {
        /*!=========================
           |    compute_A_voigt    |
           =========================
           
           Compute the A stiffness matrix in voigt notation.
        */

        A = Matrix_9x9::Zero();

        A(0,0) = lambda + 2*mu;
        A(0,1) = lambda;
        A(0,2) = lambda;
        A(1,0) = lambda;
        A(1,1) = lambda + 2*mu;
        A(1,2) = lambda;
        A(2,0) = lambda;
        A(2,1) = lambda;
        A(2,2) = lambda + 2*mu;
        A(3,3) = mu;
        A(3,6) = mu;
        A(4,4) = mu;
        A(4,7) = mu;
        A(5,5) = mu;
        A(5,8) = mu;
        A(6,3) = mu;
        A(6,6) = mu;
        A(7,4) = mu;
        A(7,7) = mu;
        A(8,5) = mu;
        A(8,8) = mu;

        return;
    }
    
    void compute_B_voigt(const double &eta,   const double &kappa, const double &nu,
                         const double &sigma, const double &tau,   SpMat &B) {
        /*!=========================
           |    compute_B_voigt    |
           =========================
           
           Compute the B stiffness matrix in voigt notation.
        */
        
        std::vector<T> tripletList;
        tripletList.reserve(21);
        
        tripletList.push_back(T(0,0,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(0,1,eta - tau));
        tripletList.push_back(T(0,2,eta - tau));
        tripletList.push_back(T(1,0,eta - tau));
        tripletList.push_back(T(1,1,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(1,2,eta - tau));
        tripletList.push_back(T(2,0,eta - tau));
        tripletList.push_back(T(2,1,eta - tau));
        tripletList.push_back(T(2,2,eta + kappa + nu - 2*sigma - tau));
        tripletList.push_back(T(3,3,kappa - sigma));
        tripletList.push_back(T(3,6,nu - sigma));
        tripletList.push_back(T(4,4,kappa - sigma));
        tripletList.push_back(T(4,7,nu - sigma));
        tripletList.push_back(T(5,5,kappa - sigma));
        tripletList.push_back(T(5,8,nu - sigma));
        tripletList.push_back(T(6,3,nu - sigma));
        tripletList.push_back(T(6,6,kappa - sigma));
        tripletList.push_back(T(7,4,nu - sigma));
        tripletList.push_back(T(7,7,kappa - sigma));
        tripletList.push_back(T(8,5,nu - sigma));
        tripletList.push_back(T(8,8,kappa - sigma));
        
        B.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_B_voigt(const double &eta,   const double &kappa, const double &nu,
                         const double &sigma, const double &tau,   Matrix_9x9 &B) {
        /*!=========================
           |    compute_B_voigt    |
           =========================
           
           Compute the B stiffness matrix in voigt notation.
        */

        B = Matrix_9x9::Zero();
        
        B(0,0) = eta + kappa + nu - 2*sigma - tau;
        B(0,1) = eta - tau;
        B(0,2) = eta - tau;
        B(1,0) = eta - tau;
        B(1,1) = eta + kappa + nu - 2*sigma - tau;
        B(1,2) = eta - tau;
        B(2,0) = eta - tau;
        B(2,1) = eta - tau;
        B(2,2) = eta + kappa + nu - 2*sigma - tau;
        B(3,3) = kappa - sigma;
        B(3,6) = nu - sigma;
        B(4,4) = kappa - sigma;
        B(4,7) = nu - sigma;
        B(5,5) = kappa - sigma;
        B(5,8) = nu - sigma;
        B(6,3) = nu - sigma;
        B(6,6) = kappa - sigma;
        B(7,4) = nu - sigma;
        B(7,7) = kappa - sigma;
        B(8,5) = nu - sigma;
        B(8,8) = kappa - sigma;
        
        return;
    }

    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, SpMat &C) {
        /*!=========================
           |    compute_C_voigt    |
           =========================
           
        Compute the C stiffness tensor in voigt 
        format.
        
        */
        std::vector<T> tripletList;
        tripletList.reserve(183);

        tripletList.push_back(T(0,0,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(0,1,tau1 + tau4 + tau5));
        tripletList.push_back(T(0,2,tau1 + tau4 + tau5));
        tripletList.push_back(T(0,14,tau2 + tau5 + tau6));
        tripletList.push_back(T(0,17,tau1 + tau2 + tau3));
        tripletList.push_back(T(0,22,tau2 + tau5 + tau6));
        tripletList.push_back(T(0,25,tau1 + tau2 + tau3));
        tripletList.push_back(T(1,0,tau1 + tau4 + tau5));
        tripletList.push_back(T(1,1,tau4 + tau7 + tau9));
        tripletList.push_back(T(1,2,tau4));
        tripletList.push_back(T(1,14,tau10 + tau5 + tau8));
        tripletList.push_back(T(1,17,tau1 + tau11 + tau8));
        tripletList.push_back(T(1,22,tau5));
        tripletList.push_back(T(1,25,tau1));
        tripletList.push_back(T(2,0,tau1 + tau4 + tau5));
        tripletList.push_back(T(2,1,tau4));
        tripletList.push_back(T(2,2,tau4 + tau7 + tau9));
        tripletList.push_back(T(2,14,tau5));
        tripletList.push_back(T(2,17,tau1));
        tripletList.push_back(T(2,22,tau10 + tau5 + tau8));
        tripletList.push_back(T(2,25,tau1 + tau11 + tau8));
        tripletList.push_back(T(3,3,tau7));
        tripletList.push_back(T(3,6,tau9));
        tripletList.push_back(T(3,13,tau10));
        tripletList.push_back(T(3,16,tau8));
        tripletList.push_back(T(3,23,tau8));
        tripletList.push_back(T(3,26,tau11));
        tripletList.push_back(T(4,4,tau10 + tau3 + tau7));
        tripletList.push_back(T(4,7,tau2 + tau8 + tau9));
        tripletList.push_back(T(4,12,tau3));
        tripletList.push_back(T(4,15,tau2));
        tripletList.push_back(T(4,18,tau1 + tau11 + tau8));
        tripletList.push_back(T(4,19,tau1));
        tripletList.push_back(T(4,20,tau1 + tau2 + tau3));
        tripletList.push_back(T(5,5,tau10 + tau3 + tau7));
        tripletList.push_back(T(5,8,tau2 + tau8 + tau9));
        tripletList.push_back(T(5,9,tau1 + tau11 + tau8));
        tripletList.push_back(T(5,10,tau1 + tau2 + tau3));
        tripletList.push_back(T(5,11,tau1));
        tripletList.push_back(T(5,21,tau2));
        tripletList.push_back(T(5,24,tau3));
        tripletList.push_back(T(6,3,tau9));
        tripletList.push_back(T(6,6,tau7));
        tripletList.push_back(T(6,13,tau8));
        tripletList.push_back(T(6,16,tau11));
        tripletList.push_back(T(6,23,tau10));
        tripletList.push_back(T(6,26,tau8));
        tripletList.push_back(T(7,4,tau2 + tau8 + tau9));
        tripletList.push_back(T(7,7,tau11 + tau6 + tau7));
        tripletList.push_back(T(7,12,tau2));
        tripletList.push_back(T(7,15,tau6));
        tripletList.push_back(T(7,18,tau10 + tau5 + tau8));
        tripletList.push_back(T(7,19,tau5));
        tripletList.push_back(T(7,20,tau2 + tau5 + tau6));
        tripletList.push_back(T(8,5,tau2 + tau8 + tau9));
        tripletList.push_back(T(8,8,tau11 + tau6 + tau7));
        tripletList.push_back(T(8,9,tau10 + tau5 + tau8));
        tripletList.push_back(T(8,10,tau2 + tau5 + tau6));
        tripletList.push_back(T(8,11,tau5));
        tripletList.push_back(T(8,21,tau6));
        tripletList.push_back(T(8,24,tau2));
        tripletList.push_back(T(9,5,tau1 + tau11 + tau8));
        tripletList.push_back(T(9,8,tau10 + tau5 + tau8));
        tripletList.push_back(T(9,9,tau4 + tau7 + tau9));
        tripletList.push_back(T(9,10,tau1 + tau4 + tau5));
        tripletList.push_back(T(9,11,tau4));
        tripletList.push_back(T(9,21,tau5));
        tripletList.push_back(T(9,24,tau1));
        tripletList.push_back(T(10,5,tau1 + tau2 + tau3));
        tripletList.push_back(T(10,8,tau2 + tau5 + tau6));
        tripletList.push_back(T(10,9,tau1 + tau4 + tau5));
        tripletList.push_back(T(10,10,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(10,11,tau1 + tau4 + tau5));
        tripletList.push_back(T(10,21,tau2 + tau5 + tau6));
        tripletList.push_back(T(10,24,tau1 + tau2 + tau3));
        tripletList.push_back(T(11,5,tau1));
        tripletList.push_back(T(11,8,tau5));
        tripletList.push_back(T(11,9,tau4));
        tripletList.push_back(T(11,10,tau1 + tau4 + tau5));
        tripletList.push_back(T(11,11,tau4 + tau7 + tau9));
        tripletList.push_back(T(11,21,tau10 + tau5 + tau8));
        tripletList.push_back(T(11,24,tau1 + tau11 + tau8));
        tripletList.push_back(T(12,4,tau3));
        tripletList.push_back(T(12,7,tau2));
        tripletList.push_back(T(12,12,tau10 + tau3 + tau7));
        tripletList.push_back(T(12,15,tau2 + tau8 + tau9));
        tripletList.push_back(T(12,18,tau1));
        tripletList.push_back(T(12,19,tau1 + tau11 + tau8));
        tripletList.push_back(T(12,20,tau1 + tau2 + tau3));
        tripletList.push_back(T(13,3,tau10));
        tripletList.push_back(T(13,6,tau8));
        tripletList.push_back(T(13,13,tau7));
        tripletList.push_back(T(13,16,tau9));
        tripletList.push_back(T(13,23,tau11));
        tripletList.push_back(T(13,26,tau8));
        tripletList.push_back(T(14,0,tau2 + tau5 + tau6));
        tripletList.push_back(T(14,1,tau10 + tau5 + tau8));
        tripletList.push_back(T(14,2,tau5));
        tripletList.push_back(T(14,14,tau11 + tau6 + tau7));
        tripletList.push_back(T(14,17,tau2 + tau8 + tau9));
        tripletList.push_back(T(14,22,tau6));
        tripletList.push_back(T(14,25,tau2));
        tripletList.push_back(T(15,4,tau2));
        tripletList.push_back(T(15,7,tau6));
        tripletList.push_back(T(15,12,tau2 + tau8 + tau9));
        tripletList.push_back(T(15,15,tau11 + tau6 + tau7));
        tripletList.push_back(T(15,18,tau5));
        tripletList.push_back(T(15,19,tau10 + tau5 + tau8));
        tripletList.push_back(T(15,20,tau2 + tau5 + tau6));
        tripletList.push_back(T(16,3,tau8));
        tripletList.push_back(T(16,6,tau11));
        tripletList.push_back(T(16,13,tau9));
        tripletList.push_back(T(16,16,tau7));
        tripletList.push_back(T(16,23,tau8));
        tripletList.push_back(T(16,26,tau10));
        tripletList.push_back(T(17,0,tau1 + tau2 + tau3));
        tripletList.push_back(T(17,1,tau1 + tau11 + tau8));
        tripletList.push_back(T(17,2,tau1));
        tripletList.push_back(T(17,14,tau2 + tau8 + tau9));
        tripletList.push_back(T(17,17,tau10 + tau3 + tau7));
        tripletList.push_back(T(17,22,tau2));
        tripletList.push_back(T(17,25,tau3));
        tripletList.push_back(T(18,4,tau1 + tau11 + tau8));
        tripletList.push_back(T(18,7,tau10 + tau5 + tau8));
        tripletList.push_back(T(18,12,tau1));
        tripletList.push_back(T(18,15,tau5));
        tripletList.push_back(T(18,18,tau4 + tau7 + tau9));
        tripletList.push_back(T(18,19,tau4));
        tripletList.push_back(T(18,20,tau1 + tau4 + tau5));
        tripletList.push_back(T(19,4,tau1));
        tripletList.push_back(T(19,7,tau5));
        tripletList.push_back(T(19,12,tau1 + tau11 + tau8));
        tripletList.push_back(T(19,15,tau10 + tau5 + tau8));
        tripletList.push_back(T(19,18,tau4));
        tripletList.push_back(T(19,19,tau4 + tau7 + tau9));
        tripletList.push_back(T(19,20,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,4,tau1 + tau2 + tau3));
        tripletList.push_back(T(20,7,tau2 + tau5 + tau6));
        tripletList.push_back(T(20,12,tau1 + tau2 + tau3));
        tripletList.push_back(T(20,15,tau2 + tau5 + tau6));
        tripletList.push_back(T(20,18,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,19,tau1 + tau4 + tau5));
        tripletList.push_back(T(20,20,2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9));
        tripletList.push_back(T(21,5,tau2));
        tripletList.push_back(T(21,8,tau6));
        tripletList.push_back(T(21,9,tau5));
        tripletList.push_back(T(21,10,tau2 + tau5 + tau6));
        tripletList.push_back(T(21,11,tau10 + tau5 + tau8));
        tripletList.push_back(T(21,21,tau11 + tau6 + tau7));
        tripletList.push_back(T(21,24,tau2 + tau8 + tau9));
        tripletList.push_back(T(22,0,tau2 + tau5 + tau6));
        tripletList.push_back(T(22,1,tau5));
        tripletList.push_back(T(22,2,tau10 + tau5 + tau8));
        tripletList.push_back(T(22,14,tau6));
        tripletList.push_back(T(22,17,tau2));
        tripletList.push_back(T(22,22,tau11 + tau6 + tau7));
        tripletList.push_back(T(22,25,tau2 + tau8 + tau9));
        tripletList.push_back(T(23,3,tau8));
        tripletList.push_back(T(23,6,tau10));
        tripletList.push_back(T(23,13,tau11));
        tripletList.push_back(T(23,16,tau8));
        tripletList.push_back(T(23,23,tau7));
        tripletList.push_back(T(23,26,tau9));
        tripletList.push_back(T(24,5,tau3));
        tripletList.push_back(T(24,8,tau2));
        tripletList.push_back(T(24,9,tau1));
        tripletList.push_back(T(24,10,tau1 + tau2 + tau3));
        tripletList.push_back(T(24,11,tau1 + tau11 + tau8));
        tripletList.push_back(T(24,21,tau2 + tau8 + tau9));
        tripletList.push_back(T(24,24,tau10 + tau3 + tau7));
        tripletList.push_back(T(25,0,tau1 + tau2 + tau3));
        tripletList.push_back(T(25,1,tau1));
        tripletList.push_back(T(25,2,tau1 + tau11 + tau8));
        tripletList.push_back(T(25,14,tau2));
        tripletList.push_back(T(25,17,tau3));
        tripletList.push_back(T(25,22,tau2 + tau8 + tau9));
        tripletList.push_back(T(25,25,tau10 + tau3 + tau7));
        tripletList.push_back(T(26,3,tau11));
        tripletList.push_back(T(26,6,tau8));
        tripletList.push_back(T(26,13,tau8));
        tripletList.push_back(T(26,16,tau10));
        tripletList.push_back(T(26,23,tau9));
        tripletList.push_back(T(26,26,tau7));
        
        C.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }
    
    void compute_C_voigt(const double &tau1,  const double &tau2,  const double &tau3,
                         const double &tau4,  const double &tau5,  const double &tau6,
                         const double &tau7,  const double &tau8,  const double &tau9,
                         const double &tau10, const double &tau11, Matrix_27x27 &C) {
        /*!=========================
           |    compute_C_voigt    |
           =========================
           
        Compute the C stiffness tensor in voigt 
        format.
        
        */
        C = Matrix_27x27::Zero();

        C( 0, 0) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C( 0, 1) = tau1 + tau4 + tau5;
        C( 0, 2) = tau1 + tau4 + tau5;
        C( 0,14) = tau2 + tau5 + tau6;
        C( 0,17) = tau1 + tau2 + tau3;
        C( 0,22) = tau2 + tau5 + tau6;
        C( 0,25) = tau1 + tau2 + tau3;
        C( 1, 0) = tau1 + tau4 + tau5;
        C( 1, 1) = tau4 + tau7 + tau9;
        C( 1, 2) = tau4;
        C( 1,14) = tau10 + tau5 + tau8;
        C( 1,17) = tau1 + tau11 + tau8;
        C( 1,22) = tau5;
        C( 1,25) = tau1;
        C( 2, 0) = tau1 + tau4 + tau5;
        C( 2, 1) = tau4;
        C( 2, 2) = tau4 + tau7 + tau9;
        C( 2,14) = tau5;
        C( 2,17) = tau1;
        C( 2,22) = tau10 + tau5 + tau8;
        C( 2,25) = tau1 + tau11 + tau8;
        C( 3, 3) = tau7;
        C( 3, 6) = tau9;
        C( 3,13) = tau10;
        C( 3,16) = tau8;
        C( 3,23) = tau8;
        C( 3,26) = tau11;
        C( 4, 4) = tau10 + tau3 + tau7;
        C( 4, 7) = tau2 + tau8 + tau9;
        C( 4,12) = tau3;
        C( 4,15) = tau2;
        C( 4,18) = tau1 + tau11 + tau8;
        C( 4,19) = tau1;
        C( 4,20) = tau1 + tau2 + tau3;
        C( 5, 5) = tau10 + tau3 + tau7;
        C( 5, 8) = tau2 + tau8 + tau9;
        C( 5, 9) = tau1 + tau11 + tau8;
        C( 5,10) = tau1 + tau2 + tau3;
        C( 5,11) = tau1;
        C( 5,21) = tau2;
        C( 5,24) = tau3;
        C( 6, 3) = tau9;
        C( 6, 6) = tau7;
        C( 6,13) = tau8;
        C( 6,16) = tau11;
        C( 6,23) = tau10;
        C( 6,26) = tau8;
        C( 7, 4) = tau2 + tau8 + tau9;
        C( 7, 7) = tau11 + tau6 + tau7;
        C( 7,12) = tau2;
        C( 7,15) = tau6;
        C( 7,18) = tau10 + tau5 + tau8;
        C( 7,19) = tau5;
        C( 7,20) = tau2 + tau5 + tau6;
        C( 8, 5) = tau2 + tau8 + tau9;
        C( 8, 8) = tau11 + tau6 + tau7;
        C( 8, 9) = tau10 + tau5 + tau8;
        C( 8,10) = tau2 + tau5 + tau6;
        C( 8,11) = tau5;
        C( 8,21) = tau6;
        C( 8,24) = tau2;
        C( 9, 5) = tau1 + tau11 + tau8;
        C( 9, 8) = tau10 + tau5 + tau8;
        C( 9, 9) = tau4 + tau7 + tau9;
        C( 9,10) = tau1 + tau4 + tau5;
        C( 9,11) = tau4;
        C( 9,21) = tau5;
        C( 9,24) = tau1;
        C(10, 5) = tau1 + tau2 + tau3;
        C(10, 8) = tau2 + tau5 + tau6;
        C(10, 9) = tau1 + tau4 + tau5;
        C(10,10) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C(10,11) = tau1 + tau4 + tau5;
        C(10,21) = tau2 + tau5 + tau6;
        C(10,24) = tau1 + tau2 + tau3;
        C(11, 5) = tau1;
        C(11, 8) = tau5;
        C(11, 9) = tau4;
        C(11,10) = tau1 + tau4 + tau5;
        C(11,11) = tau4 + tau7 + tau9;
        C(11,21) = tau10 + tau5 + tau8;
        C(11,24) = tau1 + tau11 + tau8;
        C(12, 4) = tau3;
        C(12, 7) = tau2;
        C(12,12) = tau10 + tau3 + tau7;
        C(12,15) = tau2 + tau8 + tau9;
        C(12,18) = tau1;
        C(12,19) = tau1 + tau11 + tau8;
        C(12,20) = tau1 + tau2 + tau3;
        C(13, 3) = tau10;
        C(13, 6) = tau8;
        C(13,13) = tau7;
        C(13,16) = tau9;
        C(13,23) = tau11;
        C(13,26) = tau8;
        C(14, 0) = tau2 + tau5 + tau6;
        C(14, 1) = tau10 + tau5 + tau8;
        C(14, 2) = tau5;
        C(14,14) = tau11 + tau6 + tau7;
        C(14,17) = tau2 + tau8 + tau9;
        C(14,22) = tau6;
        C(14,25) = tau2;
        C(15, 4) = tau2;
        C(15, 7) = tau6;
        C(15,12) = tau2 + tau8 + tau9;
        C(15,15) = tau11 + tau6 + tau7;
        C(15,18) = tau5;
        C(15,19) = tau10 + tau5 + tau8;
        C(15,20) = tau2 + tau5 + tau6;
        C(16, 3) = tau8;
        C(16, 6) = tau11;
        C(16,13) = tau9;
        C(16,16) = tau7;
        C(16,23) = tau8;
        C(16,26) = tau10;
        C(17, 0) = tau1 + tau2 + tau3;
        C(17, 1) = tau1 + tau11 + tau8;
        C(17, 2) = tau1;
        C(17,14) = tau2 + tau8 + tau9;
        C(17,17) = tau10 + tau3 + tau7;
        C(17,22) = tau2;
        C(17,25) = tau3;
        C(18, 4) = tau1 + tau11 + tau8;
        C(18, 7) = tau10 + tau5 + tau8;
        C(18,12) = tau1;
        C(18,15) = tau5;
        C(18,18) = tau4 + tau7 + tau9;
        C(18,19) = tau4;
        C(18,20) = tau1 + tau4 + tau5;
        C(19, 4) = tau1;
        C(19, 7) = tau5;
        C(19,12) = tau1 + tau11 + tau8;
        C(19,15) = tau10 + tau5 + tau8;
        C(19,18) = tau4;
        C(19,19) = tau4 + tau7 + tau9;
        C(19,20) = tau1 + tau4 + tau5;
        C(20, 4) = tau1 + tau2 + tau3;
        C(20, 7) = tau2 + tau5 + tau6;
        C(20,12) = tau1 + tau2 + tau3;
        C(20,15) = tau2 + tau5 + tau6;
        C(20,18) = tau1 + tau4 + tau5;
        C(20,19) = tau1 + tau4 + tau5;
        C(20,20) = 2*tau1 + tau10 + tau11 + 2*tau2 + tau3 + tau4 + 2*tau5 + tau6 + tau7 + 2*tau8 + tau9;
        C(21, 5) = tau2;
        C(21, 8) = tau6;
        C(21, 9) = tau5;
        C(21,10) = tau2 + tau5 + tau6;
        C(21,11) = tau10 + tau5 + tau8;
        C(21,21) = tau11 + tau6 + tau7;
        C(21,24) = tau2 + tau8 + tau9;
        C(22, 0) = tau2 + tau5 + tau6;
        C(22, 1) = tau5;
        C(22, 2) = tau10 + tau5 + tau8;
        C(22,14) = tau6;
        C(22,17) = tau2;
        C(22,22) = tau11 + tau6 + tau7;
        C(22,25) = tau2 + tau8 + tau9;
        C(23, 3) = tau8;
        C(23, 6) = tau10;
        C(23,13) = tau11;
        C(23,16) = tau8;
        C(23,23) = tau7;
        C(23,26) = tau9;
        C(24, 5) = tau3;
        C(24, 8) = tau2;
        C(24, 9) = tau1;
        C(24,10) = tau1 + tau2 + tau3;
        C(24,11) = tau1 + tau11 + tau8;
        C(24,21) = tau2 + tau8 + tau9;
        C(24,24) = tau10 + tau3 + tau7;
        C(25, 0) = tau1 + tau2 + tau3;
        C(25, 1) = tau1;
        C(25, 2) = tau1 + tau11 + tau8;
        C(25,14) = tau2;
        C(25,17) = tau3;
        C(25,22) = tau2 + tau8 + tau9;
        C(25,25) = tau10 + tau3 + tau7;
        C(26, 3) = tau11;
        C(26, 6) = tau8;
        C(26,13) = tau8;
        C(26,16) = tau10;
        C(26,23) = tau9;
        C(26,26) = tau7;

        return;
    }
    void compute_D_voigt(const double &sigma, const double &tau, SpMat &D){
        /*!=========================
           |    compute_D_voigt    |
           =========================
           
        Compute the D stiffness tensor in voigt 
        format.
        
        */
        std::vector<T> tripletList;
        tripletList.reserve(21);

        tripletList.push_back(T(0,0,2*sigma + tau));
        tripletList.push_back(T(0,1,tau));
        tripletList.push_back(T(0,2,tau));
        tripletList.push_back(T(1,0,tau));
        tripletList.push_back(T(1,1,2*sigma + tau));
        tripletList.push_back(T(1,2,tau));
        tripletList.push_back(T(2,0,tau));
        tripletList.push_back(T(2,1,tau));
        tripletList.push_back(T(2,2,2*sigma + tau));
        tripletList.push_back(T(3,3,sigma));
        tripletList.push_back(T(3,6,sigma));
        tripletList.push_back(T(4,4,sigma));
        tripletList.push_back(T(4,7,sigma));
        tripletList.push_back(T(5,5,sigma));
        tripletList.push_back(T(5,8,sigma));
        tripletList.push_back(T(6,3,sigma));
        tripletList.push_back(T(6,6,sigma));
        tripletList.push_back(T(7,4,sigma));
        tripletList.push_back(T(7,7,sigma));
        tripletList.push_back(T(8,5,sigma));
        tripletList.push_back(T(8,8,sigma));       
        
        D.setFromTriplets(tripletList.begin(), tripletList.end());
        return;
    }


    void compute_D_voigt(const double &sigma, const double &tau, Matrix_9x9 &D){
        /*!=========================
           |    compute_D_voigt    |
           =========================
           
        Compute the D stiffness tensor in voigt 
        format.
        
        */
        D = Matrix_9x9::Zero();

        D( 0, 0) = 2*sigma + tau;
        D( 0, 1) = tau;
        D( 0, 2) = tau;
        D( 1, 0) = tau;
        D( 1, 1) = 2*sigma + tau;
        D( 1, 2) = tau;
        D( 2, 0) = tau;
        D( 2, 1) = tau;
        D( 2, 2) = 2*sigma + tau;
        D( 3, 3) = sigma;
        D( 3, 6) = sigma;
        D( 4, 4) = sigma;
        D( 4, 7) = sigma;
        D( 5, 5) = sigma;
        D( 5, 8) = sigma;
        D( 6, 3) = sigma;
        D( 6, 6) = sigma;
        D( 7, 4) = sigma;
        D( 7, 7) = sigma;
        D( 8, 5) = sigma;
        D( 8, 8) = sigma;

        return;
    }
    
    void compute_PK2_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const SpMat &A, const SpMat &B,    const SpMat &C,
                            const SpMat &D, Vector_9 &PK2){
        /*!============================
        |    compute_PK2_stress    |
        ============================
           
        Compute the second piola kirchoff stress.
           
        */
        
        PK2 = A*E_voigt;        //Compute the first terms
        PK2 += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Vector_9 term3_4_voigt;
        deformation_measures::voigt_3x3_tensor(Temp1*(RCGinv*Psi).transpose(),term3_4_voigt);
        
        PK2 += term3_4_voigt;
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        Vector_9 term5_voigt;
        deformation_measures::voigt_3x3_tensor(Temp2*(RCGinv*Gamma).transpose(),term5_voigt);
        PK2 += term5_voigt;
        return;
    }
    
    void compute_PK2_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const Matrix_9x9 &A,        const Matrix_9x9 &B,           const Matrix_27x27 &C,
                            const Matrix_9x9 &D,        Vector_9 &PK2){
        /*!============================
        |    compute_PK2_stress    |
        ============================
           
        Compute the second piola kirchoff stress.
           
        */
        
        PK2 = A*E_voigt;        //Compute the first terms
        PK2 += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Vector_9 term3_4_voigt;
        deformation_measures::voigt_3x3_tensor(Temp1*(RCGinv*Psi).transpose(),term3_4_voigt);
        
        PK2 += term3_4_voigt;
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        Vector_9 term5_voigt;
        deformation_measures::voigt_3x3_tensor(Temp2*(RCGinv*Gamma).transpose(),term5_voigt);
        PK2 += term5_voigt;
        return;
    }

    void compute_symmetric_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                            const Matrix_3x3 &RCGinv,     const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                            const SpMat &A, const SpMat &B,    const SpMat &C,
                            const SpMat &D, Vector_9 &SIGMA){
        /*!=====================================
           |    compute_symmetric_stress    |
           ==================================
           
           Compute the symmetric stress in the reference configuration.
           
        */
        
        SIGMA = A*E_voigt;        //Compute the first terms
        SIGMA += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Matrix_3x3 symmetric_part = Temp1*(RCGinv*Psi).transpose();
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        symmetric_part += Temp2*(RCGinv*Gamma).transpose();
        Vector_9 vector_symm_part;
        deformation_measures::voigt_3x3_tensor((symmetric_part + symmetric_part.transpose()),vector_symm_part);
        
        SIGMA += vector_symm_part;
        
        return;
    }
    
    void compute_symmetric_stress(const Vector_9 &E_voigt,    const Vector_9 &E_micro_voigt, const Vector_27 &Gamma_voigt,
                                  const Matrix_3x3 &RCGinv,   const Matrix_3x3 &Psi,         const Matrix_3x9 &Gamma,
                                  const Matrix_9x9 &A,        const Matrix_9x9 &B,           const Matrix_27x27 &C,
                                  const Matrix_9x9 &D,        Vector_9 &SIGMA){
        /*!=====================================
           |    compute_symmetric_stress    |
           ==================================
           
           Compute the symmetric stress in the reference configuration.
           
        */
        
        SIGMA = A*E_voigt;        //Compute the first terms
        SIGMA += D*E_micro_voigt;
        
        //Compute the middle terms
        Matrix_3x3 Temp1;
        deformation_measures::undo_voigt_3x3_tensor(B*E_micro_voigt+D*E_voigt,Temp1);
        Matrix_3x3 symmetric_part = Temp1*(RCGinv*Psi).transpose();
        
        //Compute the end terms
        Matrix_3x9 Temp2;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,Temp2);
        symmetric_part += Temp2*(RCGinv*Gamma).transpose();
        Vector_9 vector_symm_part;
        deformation_measures::voigt_3x3_tensor((symmetric_part + symmetric_part.transpose()),vector_symm_part);
        
        SIGMA += vector_symm_part;
        
        return;
    }

    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const SpMat &C, Vector_27 &M){
        /*!=====================================
        |    compute_higher_order_stress    |
        =====================================
          
        Compute the higher order stress in the reference configuration.
          
        */
        
        M = C*Gamma_voigt; //Compute the stress (requires positive permutation)
        deformation_measures::perform_right_positive_cyclic_permutation(M); //Perform the permutation
    }
    
    void compute_higher_order_stress(const Vector_27 &Gamma_voigt, const Matrix_27x27 &C, Vector_27 &M){
        /*!=====================================
        |    compute_higher_order_stress    |
        =====================================
          
        Compute the higher order stress in the reference configuration.
          
        */
        
        M = C*Gamma_voigt; //Compute the stress (requires positive permutation)
        deformation_measures::perform_right_positive_cyclic_permutation(M); //Perform the permutation
    }

    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        Matrix_9x9 term1;
        term1 = 0.5*A;
        
        //Compute term2
        Matrix_9x9 term2;
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,term2);
        
        //Compute term3
        Matrix_9x9 term3;
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,term3);
        
        //Compute term4
        Matrix_9x9 term4;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,term4);
        
        //Assemble the derivative
        dPK2dRCG = (term1+term2+term3+term4);
        
        return;
    }
    
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,  const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,   const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,   const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        Matrix_9x9 term1;
        term1 = 0.5*A;
        
        //Compute term2
        Matrix_9x9 term2;
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,term2);
        
        //Compute term3
        Matrix_9x9 term3;
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,term3);
        
        //Compute term4
        Matrix_9x9 term4;
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,term4);
        
        //Assemble the derivative
        dPK2dRCG = (term1+term2+term3+term4);
        
        return;
    }

    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &A,      const SpMat &B,        const SpMat &C,  const SpMat &D, Matrix_9x9 (&terms)[4], Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        terms[0] = 0.5*A;
        
        //Compute term2
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,terms[1]);
        
        //Compute term3
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,terms[2]);
        
        //Compute term4
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,terms[3]);
        
        //Assemble the derivative
        dPK2dRCG = (terms[0] + terms[1] + terms[2] + terms[3]);
        
        return;
    }
    
    void compute_dPK2dRCG(const Matrix_3x3 &RCG, const Matrix_3x3 &RCGinv,  const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                          const Matrix_3x3 &E,   const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &A,   const Matrix_9x9 &B,       const Matrix_27x27 &C,   const Matrix_9x9 &D, Matrix_9x9 (&terms)[4], Matrix_9x9 &dPK2dRCG){
        /*!==========================
        |    compute_dPK2dRCG    |
        ==========================
        
        Compute the derivative of the PK2 stress w.r.t. 
        the deformation gradient.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        Matrix_9x9 dRCGinvdRCG;
        deformation_measures::compute_dAinvdA(RCGinv,dRCGinvdRCG);
        
        //Compute term1
        terms[0] = 0.5*A;
        
        //Compute term2
        T1 = RCGinv*(E_micro+Matrix_3x3::Identity());
        deformation_measures::dot_2ot_4ot(1,1,T1,0.5*D,terms[1]);
        
        //Compute term3
        V1 = (B*E_micro_voigt+D*E_voigt);
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::dot_2ot_4ot(1,0,T1*(E_micro+Matrix_3x3::Identity()).transpose(),dRCGinvdRCG,terms[2]);
        
        //Compute term4
        deformation_measures::undo_voigt_3x9_tensor(C*Gamma_voigt,T3);
        T1 = T3*Gamma.transpose();
        deformation_measures::dot_2ot_4ot(1,0,T1,dRCGinvdRCG,terms[3]);
        
        //Assemble the derivative
        dPK2dRCG = (terms[0] + terms[1] + terms[2] + terms[3]);
        
        return;
    }

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        dPK2dPsi = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        dPK2dPsi += T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        dPK2dPsi += T2;
        
        return;
    }
    
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D,       Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        dPK2dPsi = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        dPK2dPsi += T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        dPK2dPsi += T2;
        
        return;
    }

    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const SpMat &B, const SpMat &D, Matrix_9x9 (&terms)[3], Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        terms[0] = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        terms[1] = T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        terms[2] = T2;
        
        dPK2dPsi = terms[0] + terms[1] + terms[2];
        
        return;
    }
    
    void compute_dPK2dPsi(const Matrix_3x3 &RCGinv, const Matrix_3x3 &E_micro, const Vector_9 &E_voigt, const Vector_9 &E_micro_voigt,
                          const Matrix_9x9 &B,      const Matrix_9x9 &D,       Matrix_9x9 (&terms)[3],  Matrix_9x9 &dPK2dPsi){
        /*!==========================
        |    compute_dPK2dPsi    |
        ==========================
        
        Compute the derivative of the second piola kirchoff 
        stress w.r.t. the deformation measure Psi.
        
        */
        
        //Initialize temporary terms
        Vector_9   V1;
        Vector_27  V2;
        Matrix_3x3 T1;
        Matrix_9x9 T2;
        Matrix_3x9 T3;
        
        //Add term1
        terms[0] = D;
        
        //Add term2
        deformation_measures::dot_2ot_4ot(1, 1, RCGinv*(E_micro + Matrix_3x3::Identity()), B, T2);
        terms[1] = T2;
        
        //Add term3
        V1 = B*E_micro_voigt + D*E_voigt;
        deformation_measures::undo_voigt_3x3_tensor(V1,T1);
        deformation_measures::two_sot_to_fot(2,T1,RCGinv,T2);
        terms[2] = T2;
        
        dPK2dPsi = terms[0] + terms[1] + terms[2];
        
        return;
    }

    void compute_dPK2dGamma_term2(const Vector_27 &term1, const Matrix_27x27 &C, Matrix_9x27 &term2){
        /*!==================================
        |    compute_dPK2dGamma_term2    |
        ==================================
        
        Compute term2 of dPK2dGamma. Note that this is identical as 
        term2 for dSIGMAdGamma if symmetrized.
        
        */
        
        //Extract term1
        double term1111 = term1(0);
        double term1122 = term1(1);
        double term1133 = term1(2);
        double term1123 = term1(3);
        double term1113 = term1(4);
        double term1112 = term1(5);
        double term1132 = term1(6);
        double term1131 = term1(7);
        double term1121 = term1(8);
        double term1211 = term1(9);
        double term1222 = term1(10);
        double term1233 = term1(11);
        double term1223 = term1(12);
        double term1213 = term1(13);
        double term1212 = term1(14);
        double term1232 = term1(15);
        double term1231 = term1(16);
        double term1221 = term1(17);
        double term1311 = term1(18);
        double term1322 = term1(19);
        double term1333 = term1(20);
        double term1323 = term1(21);
        double term1313 = term1(22);
        double term1312 = term1(23);
        double term1332 = term1(24);
        double term1331 = term1(25);
        double term1321 = term1(26);

        //Extract C
        double C111111 = C(0,0);
        double C111122 = C(0,1);
        double C111133 = C(0,2);
        double C111123 = C(0,3);
        double C111113 = C(0,4);
        double C111112 = C(0,5);
        double C111132 = C(0,6);
        double C111131 = C(0,7);
        double C111121 = C(0,8);
        double C111211 = C(0,9);
        double C111222 = C(0,10);
        double C111233 = C(0,11);
        double C111223 = C(0,12);
        double C111213 = C(0,13);
        double C111212 = C(0,14);
        double C111232 = C(0,15);
        double C111231 = C(0,16);
        double C111221 = C(0,17);
        double C111311 = C(0,18);
        double C111322 = C(0,19);
        double C111333 = C(0,20);
        double C111323 = C(0,21);
        double C111313 = C(0,22);
        double C111312 = C(0,23);
        double C111332 = C(0,24);
        double C111331 = C(0,25);
        double C111321 = C(0,26);
        double C122111 = C(1,0);
        double C122122 = C(1,1);
        double C122133 = C(1,2);
        double C122123 = C(1,3);
        double C122113 = C(1,4);
        double C122112 = C(1,5);
        double C122132 = C(1,6);
        double C122131 = C(1,7);
        double C122121 = C(1,8);
        double C122211 = C(1,9);
        double C122222 = C(1,10);
        double C122233 = C(1,11);
        double C122223 = C(1,12);
        double C122213 = C(1,13);
        double C122212 = C(1,14);
        double C122232 = C(1,15);
        double C122231 = C(1,16);
        double C122221 = C(1,17);
        double C122311 = C(1,18);
        double C122322 = C(1,19);
        double C122333 = C(1,20);
        double C122323 = C(1,21);
        double C122313 = C(1,22);
        double C122312 = C(1,23);
        double C122332 = C(1,24);
        double C122331 = C(1,25);
        double C122321 = C(1,26);
        double C133111 = C(2,0);
        double C133122 = C(2,1);
        double C133133 = C(2,2);
        double C133123 = C(2,3);
        double C133113 = C(2,4);
        double C133112 = C(2,5);
        double C133132 = C(2,6);
        double C133131 = C(2,7);
        double C133121 = C(2,8);
        double C133211 = C(2,9);
        double C133222 = C(2,10);
        double C133233 = C(2,11);
        double C133223 = C(2,12);
        double C133213 = C(2,13);
        double C133212 = C(2,14);
        double C133232 = C(2,15);
        double C133231 = C(2,16);
        double C133221 = C(2,17);
        double C133311 = C(2,18);
        double C133322 = C(2,19);
        double C133333 = C(2,20);
        double C133323 = C(2,21);
        double C133313 = C(2,22);
        double C133312 = C(2,23);
        double C133332 = C(2,24);
        double C133331 = C(2,25);
        double C133321 = C(2,26);
        double C123111 = C(3,0);
        double C123122 = C(3,1);
        double C123133 = C(3,2);
        double C123123 = C(3,3);
        double C123113 = C(3,4);
        double C123112 = C(3,5);
        double C123132 = C(3,6);
        double C123131 = C(3,7);
        double C123121 = C(3,8);
        double C123211 = C(3,9);
        double C123222 = C(3,10);
        double C123233 = C(3,11);
        double C123223 = C(3,12);
        double C123213 = C(3,13);
        double C123212 = C(3,14);
        double C123232 = C(3,15);
        double C123231 = C(3,16);
        double C123221 = C(3,17);
        double C123311 = C(3,18);
        double C123322 = C(3,19);
        double C123333 = C(3,20);
        double C123323 = C(3,21);
        double C123313 = C(3,22);
        double C123312 = C(3,23);
        double C123332 = C(3,24);
        double C123331 = C(3,25);
        double C123321 = C(3,26);
        double C113111 = C(4,0);
        double C113122 = C(4,1);
        double C113133 = C(4,2);
        double C113123 = C(4,3);
        double C113113 = C(4,4);
        double C113112 = C(4,5);
        double C113132 = C(4,6);
        double C113131 = C(4,7);
        double C113121 = C(4,8);
        double C113211 = C(4,9);
        double C113222 = C(4,10);
        double C113233 = C(4,11);
        double C113223 = C(4,12);
        double C113213 = C(4,13);
        double C113212 = C(4,14);
        double C113232 = C(4,15);
        double C113231 = C(4,16);
        double C113221 = C(4,17);
        double C113311 = C(4,18);
        double C113322 = C(4,19);
        double C113333 = C(4,20);
        double C113323 = C(4,21);
        double C113313 = C(4,22);
        double C113312 = C(4,23);
        double C113332 = C(4,24);
        double C113331 = C(4,25);
        double C113321 = C(4,26);
        double C112111 = C(5,0);
        double C112122 = C(5,1);
        double C112133 = C(5,2);
        double C112123 = C(5,3);
        double C112113 = C(5,4);
        double C112112 = C(5,5);
        double C112132 = C(5,6);
        double C112131 = C(5,7);
        double C112121 = C(5,8);
        double C112211 = C(5,9);
        double C112222 = C(5,10);
        double C112233 = C(5,11);
        double C112223 = C(5,12);
        double C112213 = C(5,13);
        double C112212 = C(5,14);
        double C112232 = C(5,15);
        double C112231 = C(5,16);
        double C112221 = C(5,17);
        double C112311 = C(5,18);
        double C112322 = C(5,19);
        double C112333 = C(5,20);
        double C112323 = C(5,21);
        double C112313 = C(5,22);
        double C112312 = C(5,23);
        double C112332 = C(5,24);
        double C112331 = C(5,25);
        double C112321 = C(5,26);
        double C132111 = C(6,0);
        double C132122 = C(6,1);
        double C132133 = C(6,2);
        double C132123 = C(6,3);
        double C132113 = C(6,4);
        double C132112 = C(6,5);
        double C132132 = C(6,6);
        double C132131 = C(6,7);
        double C132121 = C(6,8);
        double C132211 = C(6,9);
        double C132222 = C(6,10);
        double C132233 = C(6,11);
        double C132223 = C(6,12);
        double C132213 = C(6,13);
        double C132212 = C(6,14);
        double C132232 = C(6,15);
        double C132231 = C(6,16);
        double C132221 = C(6,17);
        double C132311 = C(6,18);
        double C132322 = C(6,19);
        double C132333 = C(6,20);
        double C132323 = C(6,21);
        double C132313 = C(6,22);
        double C132312 = C(6,23);
        double C132332 = C(6,24);
        double C132331 = C(6,25);
        double C132321 = C(6,26);
        double C131111 = C(7,0);
        double C131122 = C(7,1);
        double C131133 = C(7,2);
        double C131123 = C(7,3);
        double C131113 = C(7,4);
        double C131112 = C(7,5);
        double C131132 = C(7,6);
        double C131131 = C(7,7);
        double C131121 = C(7,8);
        double C131211 = C(7,9);
        double C131222 = C(7,10);
        double C131233 = C(7,11);
        double C131223 = C(7,12);
        double C131213 = C(7,13);
        double C131212 = C(7,14);
        double C131232 = C(7,15);
        double C131231 = C(7,16);
        double C131221 = C(7,17);
        double C131311 = C(7,18);
        double C131322 = C(7,19);
        double C131333 = C(7,20);
        double C131323 = C(7,21);
        double C131313 = C(7,22);
        double C131312 = C(7,23);
        double C131332 = C(7,24);
        double C131331 = C(7,25);
        double C131321 = C(7,26);
        double C121111 = C(8,0);
        double C121122 = C(8,1);
        double C121133 = C(8,2);
        double C121123 = C(8,3);
        double C121113 = C(8,4);
        double C121112 = C(8,5);
        double C121132 = C(8,6);
        double C121131 = C(8,7);
        double C121121 = C(8,8);
        double C121211 = C(8,9);
        double C121222 = C(8,10);
        double C121233 = C(8,11);
        double C121223 = C(8,12);
        double C121213 = C(8,13);
        double C121212 = C(8,14);
        double C121232 = C(8,15);
        double C121231 = C(8,16);
        double C121221 = C(8,17);
        double C121311 = C(8,18);
        double C121322 = C(8,19);
        double C121333 = C(8,20);
        double C121323 = C(8,21);
        double C121313 = C(8,22);
        double C121312 = C(8,23);
        double C121332 = C(8,24);
        double C121331 = C(8,25);
        double C121321 = C(8,26);
        double C211111 = C(9,0);
        double C211122 = C(9,1);
        double C211133 = C(9,2);
        double C211123 = C(9,3);
        double C211113 = C(9,4);
        double C211112 = C(9,5);
        double C211132 = C(9,6);
        double C211131 = C(9,7);
        double C211121 = C(9,8);
        double C211211 = C(9,9);
        double C211222 = C(9,10);
        double C211233 = C(9,11);
        double C211223 = C(9,12);
        double C211213 = C(9,13);
        double C211212 = C(9,14);
        double C211232 = C(9,15);
        double C211231 = C(9,16);
        double C211221 = C(9,17);
        double C211311 = C(9,18);
        double C211322 = C(9,19);
        double C211333 = C(9,20);
        double C211323 = C(9,21);
        double C211313 = C(9,22);
        double C211312 = C(9,23);
        double C211332 = C(9,24);
        double C211331 = C(9,25);
        double C211321 = C(9,26);
        double C222111 = C(10,0);
        double C222122 = C(10,1);
        double C222133 = C(10,2);
        double C222123 = C(10,3);
        double C222113 = C(10,4);
        double C222112 = C(10,5);
        double C222132 = C(10,6);
        double C222131 = C(10,7);
        double C222121 = C(10,8);
        double C222211 = C(10,9);
        double C222222 = C(10,10);
        double C222233 = C(10,11);
        double C222223 = C(10,12);
        double C222213 = C(10,13);
        double C222212 = C(10,14);
        double C222232 = C(10,15);
        double C222231 = C(10,16);
        double C222221 = C(10,17);
        double C222311 = C(10,18);
        double C222322 = C(10,19);
        double C222333 = C(10,20);
        double C222323 = C(10,21);
        double C222313 = C(10,22);
        double C222312 = C(10,23);
        double C222332 = C(10,24);
        double C222331 = C(10,25);
        double C222321 = C(10,26);
        double C233111 = C(11,0);
        double C233122 = C(11,1);
        double C233133 = C(11,2);
        double C233123 = C(11,3);
        double C233113 = C(11,4);
        double C233112 = C(11,5);
        double C233132 = C(11,6);
        double C233131 = C(11,7);
        double C233121 = C(11,8);
        double C233211 = C(11,9);
        double C233222 = C(11,10);
        double C233233 = C(11,11);
        double C233223 = C(11,12);
        double C233213 = C(11,13);
        double C233212 = C(11,14);
        double C233232 = C(11,15);
        double C233231 = C(11,16);
        double C233221 = C(11,17);
        double C233311 = C(11,18);
        double C233322 = C(11,19);
        double C233333 = C(11,20);
        double C233323 = C(11,21);
        double C233313 = C(11,22);
        double C233312 = C(11,23);
        double C233332 = C(11,24);
        double C233331 = C(11,25);
        double C233321 = C(11,26);
        double C223111 = C(12,0);
        double C223122 = C(12,1);
        double C223133 = C(12,2);
        double C223123 = C(12,3);
        double C223113 = C(12,4);
        double C223112 = C(12,5);
        double C223132 = C(12,6);
        double C223131 = C(12,7);
        double C223121 = C(12,8);
        double C223211 = C(12,9);
        double C223222 = C(12,10);
        double C223233 = C(12,11);
        double C223223 = C(12,12);
        double C223213 = C(12,13);
        double C223212 = C(12,14);
        double C223232 = C(12,15);
        double C223231 = C(12,16);
        double C223221 = C(12,17);
        double C223311 = C(12,18);
        double C223322 = C(12,19);
        double C223333 = C(12,20);
        double C223323 = C(12,21);
        double C223313 = C(12,22);
        double C223312 = C(12,23);
        double C223332 = C(12,24);
        double C223331 = C(12,25);
        double C223321 = C(12,26);
        double C213111 = C(13,0);
        double C213122 = C(13,1);
        double C213133 = C(13,2);
        double C213123 = C(13,3);
        double C213113 = C(13,4);
        double C213112 = C(13,5);
        double C213132 = C(13,6);
        double C213131 = C(13,7);
        double C213121 = C(13,8);
        double C213211 = C(13,9);
        double C213222 = C(13,10);
        double C213233 = C(13,11);
        double C213223 = C(13,12);
        double C213213 = C(13,13);
        double C213212 = C(13,14);
        double C213232 = C(13,15);
        double C213231 = C(13,16);
        double C213221 = C(13,17);
        double C213311 = C(13,18);
        double C213322 = C(13,19);
        double C213333 = C(13,20);
        double C213323 = C(13,21);
        double C213313 = C(13,22);
        double C213312 = C(13,23);
        double C213332 = C(13,24);
        double C213331 = C(13,25);
        double C213321 = C(13,26);
        double C212111 = C(14,0);
        double C212122 = C(14,1);
        double C212133 = C(14,2);
        double C212123 = C(14,3);
        double C212113 = C(14,4);
        double C212112 = C(14,5);
        double C212132 = C(14,6);
        double C212131 = C(14,7);
        double C212121 = C(14,8);
        double C212211 = C(14,9);
        double C212222 = C(14,10);
        double C212233 = C(14,11);
        double C212223 = C(14,12);
        double C212213 = C(14,13);
        double C212212 = C(14,14);
        double C212232 = C(14,15);
        double C212231 = C(14,16);
        double C212221 = C(14,17);
        double C212311 = C(14,18);
        double C212322 = C(14,19);
        double C212333 = C(14,20);
        double C212323 = C(14,21);
        double C212313 = C(14,22);
        double C212312 = C(14,23);
        double C212332 = C(14,24);
        double C212331 = C(14,25);
        double C212321 = C(14,26);
        double C232111 = C(15,0);
        double C232122 = C(15,1);
        double C232133 = C(15,2);
        double C232123 = C(15,3);
        double C232113 = C(15,4);
        double C232112 = C(15,5);
        double C232132 = C(15,6);
        double C232131 = C(15,7);
        double C232121 = C(15,8);
        double C232211 = C(15,9);
        double C232222 = C(15,10);
        double C232233 = C(15,11);
        double C232223 = C(15,12);
        double C232213 = C(15,13);
        double C232212 = C(15,14);
        double C232232 = C(15,15);
        double C232231 = C(15,16);
        double C232221 = C(15,17);
        double C232311 = C(15,18);
        double C232322 = C(15,19);
        double C232333 = C(15,20);
        double C232323 = C(15,21);
        double C232313 = C(15,22);
        double C232312 = C(15,23);

        double C232332 = C(15,24);
        double C232331 = C(15,25);
        double C232321 = C(15,26);
        double C231111 = C(16,0);
        double C231122 = C(16,1);
        double C231133 = C(16,2);
        double C231123 = C(16,3);
        double C231113 = C(16,4);
        double C231112 = C(16,5);
        double C231132 = C(16,6);
        double C231131 = C(16,7);
        double C231121 = C(16,8);
        double C231211 = C(16,9);
        double C231222 = C(16,10);
        double C231233 = C(16,11);
        double C231223 = C(16,12);
        double C231213 = C(16,13);
        double C231212 = C(16,14);
        double C231232 = C(16,15);
        double C231231 = C(16,16);
        double C231221 = C(16,17);
        double C231311 = C(16,18);
        double C231322 = C(16,19);
        double C231333 = C(16,20);
        double C231323 = C(16,21);
        double C231313 = C(16,22);
        double C231312 = C(16,23);
        double C231332 = C(16,24);
        double C231331 = C(16,25);
        double C231321 = C(16,26);
        double C221111 = C(17,0);
        double C221122 = C(17,1);
        double C221133 = C(17,2);
        double C221123 = C(17,3);
        double C221113 = C(17,4);
        double C221112 = C(17,5);
        double C221132 = C(17,6);
        double C221131 = C(17,7);
        double C221121 = C(17,8);
        double C221211 = C(17,9);
        double C221222 = C(17,10);
        double C221233 = C(17,11);
        double C221223 = C(17,12);
        double C221213 = C(17,13);
        double C221212 = C(17,14);
        double C221232 = C(17,15);
        double C221231 = C(17,16);
        double C221221 = C(17,17);
        double C221311 = C(17,18);
        double C221322 = C(17,19);
        double C221333 = C(17,20);
        double C221323 = C(17,21);
        double C221313 = C(17,22);
        double C221312 = C(17,23);
        double C221332 = C(17,24);
        double C221331 = C(17,25);
        double C221321 = C(17,26);
        double C311111 = C(18,0);
        double C311122 = C(18,1);
        double C311133 = C(18,2);
        double C311123 = C(18,3);
        double C311113 = C(18,4);
        double C311112 = C(18,5);
        double C311132 = C(18,6);
        double C311131 = C(18,7);
        double C311121 = C(18,8);
        double C311211 = C(18,9);
        double C311222 = C(18,10);
        double C311233 = C(18,11);
        double C311223 = C(18,12);
        double C311213 = C(18,13);
        double C311212 = C(18,14);
        double C311232 = C(18,15);
        double C311231 = C(18,16);
        double C311221 = C(18,17);
        double C311311 = C(18,18);
        double C311322 = C(18,19);
        double C311333 = C(18,20);
        double C311323 = C(18,21);
        double C311313 = C(18,22);
        double C311312 = C(18,23);
        double C311332 = C(18,24);
        double C311331 = C(18,25);
        double C311321 = C(18,26);
        double C322111 = C(19,0);
        double C322122 = C(19,1);
        double C322133 = C(19,2);
        double C322123 = C(19,3);
        double C322113 = C(19,4);
        double C322112 = C(19,5);
        double C322132 = C(19,6);
        double C322131 = C(19,7);
        double C322121 = C(19,8);
        double C322211 = C(19,9);
        double C322222 = C(19,10);
        double C322233 = C(19,11);
        double C322223 = C(19,12);
        double C322213 = C(19,13);
        double C322212 = C(19,14);
        double C322232 = C(19,15);
        double C322231 = C(19,16);
        double C322221 = C(19,17);
        double C322311 = C(19,18);
        double C322322 = C(19,19);
        double C322333 = C(19,20);
        double C322323 = C(19,21);
        double C322313 = C(19,22);
        double C322312 = C(19,23);
        double C322332 = C(19,24);
        double C322331 = C(19,25);
        double C322321 = C(19,26);
        double C333111 = C(20,0);
        double C333122 = C(20,1);
        double C333133 = C(20,2);
        double C333123 = C(20,3);
        double C333113 = C(20,4);
        double C333112 = C(20,5);
        double C333132 = C(20,6);
        double C333131 = C(20,7);
        double C333121 = C(20,8);
        double C333211 = C(20,9);
        double C333222 = C(20,10);
        double C333233 = C(20,11);
        double C333223 = C(20,12);
        double C333213 = C(20,13);
        double C333212 = C(20,14);
        double C333232 = C(20,15);
        double C333231 = C(20,16);
        double C333221 = C(20,17);
        double C333311 = C(20,18);
        double C333322 = C(20,19);
        double C333333 = C(20,20);
        double C333323 = C(20,21);
        double C333313 = C(20,22);
        double C333312 = C(20,23);
        double C333332 = C(20,24);
        double C333331 = C(20,25);
        double C333321 = C(20,26);
        double C323111 = C(21,0);
        double C323122 = C(21,1);
        double C323133 = C(21,2);
        double C323123 = C(21,3);
        double C323113 = C(21,4);
        double C323112 = C(21,5);
        double C323132 = C(21,6);
        double C323131 = C(21,7);
        double C323121 = C(21,8);
        double C323211 = C(21,9);
        double C323222 = C(21,10);
        double C323233 = C(21,11);
        double C323223 = C(21,12);
        double C323213 = C(21,13);
        double C323212 = C(21,14);
        double C323232 = C(21,15);
        double C323231 = C(21,16);
        double C323221 = C(21,17);
        double C323311 = C(21,18);
        double C323322 = C(21,19);
        double C323333 = C(21,20);
        double C323323 = C(21,21);
        double C323313 = C(21,22);
        double C323312 = C(21,23);
        double C323332 = C(21,24);
        double C323331 = C(21,25);
        double C323321 = C(21,26);
        double C313111 = C(22,0);
        double C313122 = C(22,1);
        double C313133 = C(22,2);
        double C313123 = C(22,3);
        double C313113 = C(22,4);
        double C313112 = C(22,5);
        double C313132 = C(22,6);
        double C313131 = C(22,7);
        double C313121 = C(22,8);
        double C313211 = C(22,9);
        double C313222 = C(22,10);
        double C313233 = C(22,11);
        double C313223 = C(22,12);
        double C313213 = C(22,13);
        double C313212 = C(22,14);
        double C313232 = C(22,15);
        double C313231 = C(22,16);
        double C313221 = C(22,17);
        double C313311 = C(22,18);
        double C313322 = C(22,19);
        double C313333 = C(22,20);
        double C313323 = C(22,21);
        double C313313 = C(22,22);
        double C313312 = C(22,23);
        double C313332 = C(22,24);
        double C313331 = C(22,25);
        double C313321 = C(22,26);
        double C312111 = C(23,0);
        double C312122 = C(23,1);
        double C312133 = C(23,2);
        double C312123 = C(23,3);
        double C312113 = C(23,4);
        double C312112 = C(23,5);
        double C312132 = C(23,6);
        double C312131 = C(23,7);
        double C312121 = C(23,8);
        double C312211 = C(23,9);
        double C312222 = C(23,10);
        double C312233 = C(23,11);
        double C312223 = C(23,12);
        double C312213 = C(23,13);
        double C312212 = C(23,14);
        double C312232 = C(23,15);
        double C312231 = C(23,16);
        double C312221 = C(23,17);
        double C312311 = C(23,18);
        double C312322 = C(23,19);
        double C312333 = C(23,20);
        double C312323 = C(23,21);
        double C312313 = C(23,22);
        double C312312 = C(23,23);
        double C312332 = C(23,24);
        double C312331 = C(23,25);
        double C312321 = C(23,26);
        double C332111 = C(24,0);
        double C332122 = C(24,1);
        double C332133 = C(24,2);
        double C332123 = C(24,3);
        double C332113 = C(24,4);
        double C332112 = C(24,5);
        double C332132 = C(24,6);
        double C332131 = C(24,7);
        double C332121 = C(24,8);
        double C332211 = C(24,9);
        double C332222 = C(24,10);
        double C332233 = C(24,11);
        double C332223 = C(24,12);
        double C332213 = C(24,13);
        double C332212 = C(24,14);
        double C332232 = C(24,15);
        double C332231 = C(24,16);
        double C332221 = C(24,17);
        double C332311 = C(24,18);
        double C332322 = C(24,19);
        double C332333 = C(24,20);
        double C332323 = C(24,21);
        double C332313 = C(24,22);
        double C332312 = C(24,23);
        double C332332 = C(24,24);
        double C332331 = C(24,25);
        double C332321 = C(24,26);
        double C331111 = C(25,0);
        double C331122 = C(25,1);
        double C331133 = C(25,2);
        double C331123 = C(25,3);
        double C331113 = C(25,4);
        double C331112 = C(25,5);
        double C331132 = C(25,6);
        double C331131 = C(25,7);
        double C331121 = C(25,8);
        double C331211 = C(25,9);
        double C331222 = C(25,10);
        double C331233 = C(25,11);
        double C331223 = C(25,12);
        double C331213 = C(25,13);
        double C331212 = C(25,14);
        double C331232 = C(25,15);
        double C331231 = C(25,16);
        double C331221 = C(25,17);
        double C331311 = C(25,18);
        double C331322 = C(25,19);
        double C331333 = C(25,20);
        double C331323 = C(25,21);
        double C331313 = C(25,22);
        double C331312 = C(25,23);
        double C331332 = C(25,24);
        double C331331 = C(25,25);
        double C331321 = C(25,26);
        double C321111 = C(26,0);
        double C321122 = C(26,1);
        double C321133 = C(26,2);
        double C321123 = C(26,3);
        double C321113 = C(26,4);
        double C321112 = C(26,5);
        double C321132 = C(26,6);
        double C321131 = C(26,7);
        double C321121 = C(26,8);
        double C321211 = C(26,9);
        double C321222 = C(26,10);
        double C321233 = C(26,11);
        double C321223 = C(26,12);
        double C321213 = C(26,13);
        double C321212 = C(26,14);
        double C321232 = C(26,15);
        double C321231 = C(26,16);
        double C321221 = C(26,17);
        double C321311 = C(26,18);
        double C321322 = C(26,19);
        double C321333 = C(26,20);
        double C321323 = C(26,21);
        double C321313 = C(26,22);
        double C321312 = C(26,23);
        double C321332 = C(26,24);
        double C321331 = C(26,25);
        double C321321 = C(26,26);

         //Assemble term2
        term2(0,0) = C111111*term1111 + C112111*term1112 + C113111*term1113 + C121111*term1121 + C122111*term1122 + C123111*term1123 + C131111*term1131 + C132111*term1132 + C133111*term1133;
        term2(0,1) = C111122*term1111 + C112122*term1112 + C113122*term1113 + C121122*term1121 + C122122*term1122 + C123122*term1123 + C131122*term1131 + C132122*term1132 + C133122*term1133;
        term2(0,2) = C111133*term1111 + C112133*term1112 + C113133*term1113 + C121133*term1121 + C122133*term1122 + C123133*term1123 + C131133*term1131 + C132133*term1132 + C133133*term1133;
        term2(0,3) = C111123*term1111 + C112123*term1112 + C113123*term1113 + C121123*term1121 + C122123*term1122 + C123123*term1123 + C131123*term1131 + C132123*term1132 + C133123*term1133;
        term2(0,4) = C111113*term1111 + C112113*term1112 + C113113*term1113 + C121113*term1121 + C122113*term1122 + C123113*term1123 + C131113*term1131 + C132113*term1132 + C133113*term1133;
        term2(0,5) = C111112*term1111 + C112112*term1112 + C113112*term1113 + C121112*term1121 + C122112*term1122 + C123112*term1123 + C131112*term1131 + C132112*term1132 + C133112*term1133;
        term2(0,6) = C111132*term1111 + C112132*term1112 + C113132*term1113 + C121132*term1121 + C122132*term1122 + C123132*term1123 + C131132*term1131 + C132132*term1132 + C133132*term1133;
        term2(0,7) = C111131*term1111 + C112131*term1112 + C113131*term1113 + C121131*term1121 + C122131*term1122 + C123131*term1123 + C131131*term1131 + C132131*term1132 + C133131*term1133;
        term2(0,8) = C111121*term1111 + C112121*term1112 + C113121*term1113 + C121121*term1121 + C122121*term1122 + C123121*term1123 + C131121*term1131 + C132121*term1132 + C133121*term1133;
        term2(0,9) = C111211*term1111 + C112211*term1112 + C113211*term1113 + C121211*term1121 + C122211*term1122 + C123211*term1123 + C131211*term1131 + C132211*term1132 + C133211*term1133;
        term2(0,10) = C111222*term1111 + C112222*term1112 + C113222*term1113 + C121222*term1121 + C122222*term1122 + C123222*term1123 + C131222*term1131 + C132222*term1132 + C133222*term1133;
        term2(0,11) = C111233*term1111 + C112233*term1112 + C113233*term1113 + C121233*term1121 + C122233*term1122 + C123233*term1123 + C131233*term1131 + C132233*term1132 + C133233*term1133;
        term2(0,12) = C111223*term1111 + C112223*term1112 + C113223*term1113 + C121223*term1121 + C122223*term1122 + C123223*term1123 + C131223*term1131 + C132223*term1132 + C133223*term1133;
        term2(0,13) = C111213*term1111 + C112213*term1112 + C113213*term1113 + C121213*term1121 + C122213*term1122 + C123213*term1123 + C131213*term1131 + C132213*term1132 + C133213*term1133;
        term2(0,14) = C111212*term1111 + C112212*term1112 + C113212*term1113 + C121212*term1121 + C122212*term1122 + C123212*term1123 + C131212*term1131 + C132212*term1132 + C133212*term1133;
        term2(0,15) = C111232*term1111 + C112232*term1112 + C113232*term1113 + C121232*term1121 + C122232*term1122 + C123232*term1123 + C131232*term1131 + C132232*term1132 + C133232*term1133;
        term2(0,16) = C111231*term1111 + C112231*term1112 + C113231*term1113 + C121231*term1121 + C122231*term1122 + C123231*term1123 + C131231*term1131 + C132231*term1132 + C133231*term1133;
        term2(0,17) = C111221*term1111 + C112221*term1112 + C113221*term1113 + C121221*term1121 + C122221*term1122 + C123221*term1123 + C131221*term1131 + C132221*term1132 + C133221*term1133;
        term2(0,18) = C111311*term1111 + C112311*term1112 + C113311*term1113 + C121311*term1121 + C122311*term1122 + C123311*term1123 + C131311*term1131 + C132311*term1132 + C133311*term1133;
        term2(0,19) = C111322*term1111 + C112322*term1112 + C113322*term1113 + C121322*term1121 + C122322*term1122 + C123322*term1123 + C131322*term1131 + C132322*term1132 + C133322*term1133;
        term2(0,20) = C111333*term1111 + C112333*term1112 + C113333*term1113 + C121333*term1121 + C122333*term1122 + C123333*term1123 + C131333*term1131 + C132333*term1132 + C133333*term1133;
        term2(0,21) = C111323*term1111 + C112323*term1112 + C113323*term1113 + C121323*term1121 + C122323*term1122 + C123323*term1123 + C131323*term1131 + C132323*term1132 + C133323*term1133;
        term2(0,22) = C111313*term1111 + C112313*term1112 + C113313*term1113 + C121313*term1121 + C122313*term1122 + C123313*term1123 + C131313*term1131 + C132313*term1132 + C133313*term1133;
        term2(0,23) = C111312*term1111 + C112312*term1112 + C113312*term1113 + C121312*term1121 + C122312*term1122 + C123312*term1123 + C131312*term1131 + C132312*term1132 + C133312*term1133;
        term2(0,24) = C111332*term1111 + C112332*term1112 + C113332*term1113 + C121332*term1121 + C122332*term1122 + C123332*term1123 + C131332*term1131 + C132332*term1132 + C133332*term1133;
        term2(0,25) = C111331*term1111 + C112331*term1112 + C113331*term1113 + C121331*term1121 + C122331*term1122 + C123331*term1123 + C131331*term1131 + C132331*term1132 + C133331*term1133;
        term2(0,26) = C111321*term1111 + C112321*term1112 + C113321*term1113 + C121321*term1121 + C122321*term1122 + C123321*term1123 + C131321*term1131 + C132321*term1132 + C133321*term1133;
        term2(1,0) = C211111*term1211 + C212111*term1212 + C213111*term1213 + C221111*term1221 + C222111*term1222 + C223111*term1223 + C231111*term1231 + C232111*term1232 + C233111*term1233;
        term2(1,1) = C211122*term1211 + C212122*term1212 + C213122*term1213 + C221122*term1221 + C222122*term1222 + C223122*term1223 + C231122*term1231 + C232122*term1232 + C233122*term1233;
        term2(1,2) = C211133*term1211 + C212133*term1212 + C213133*term1213 + C221133*term1221 + C222133*term1222 + C223133*term1223 + C231133*term1231 + C232133*term1232 + C233133*term1233;
        term2(1,3) = C211123*term1211 + C212123*term1212 + C213123*term1213 + C221123*term1221 + C222123*term1222 + C223123*term1223 + C231123*term1231 + C232123*term1232 + C233123*term1233;
        term2(1,4) = C211113*term1211 + C212113*term1212 + C213113*term1213 + C221113*term1221 + C222113*term1222 + C223113*term1223 + C231113*term1231 + C232113*term1232 + C233113*term1233;
        term2(1,5) = C211112*term1211 + C212112*term1212 + C213112*term1213 + C221112*term1221 + C222112*term1222 + C223112*term1223 + C231112*term1231 + C232112*term1232 + C233112*term1233;
        term2(1,6) = C211132*term1211 + C212132*term1212 + C213132*term1213 + C221132*term1221 + C222132*term1222 + C223132*term1223 + C231132*term1231 + C232132*term1232 + C233132*term1233;
        term2(1,7) = C211131*term1211 + C212131*term1212 + C213131*term1213 + C221131*term1221 + C222131*term1222 + C223131*term1223 + C231131*term1231 + C232131*term1232 + C233131*term1233;
        term2(1,8) = C211121*term1211 + C212121*term1212 + C213121*term1213 + C221121*term1221 + C222121*term1222 + C223121*term1223 + C231121*term1231 + C232121*term1232 + C233121*term1233;
        term2(1,9) = C211211*term1211 + C212211*term1212 + C213211*term1213 + C221211*term1221 + C222211*term1222 + C223211*term1223 + C231211*term1231 + C232211*term1232 + C233211*term1233;
        term2(1,10) = C211222*term1211 + C212222*term1212 + C213222*term1213 + C221222*term1221 + C222222*term1222 + C223222*term1223 + C231222*term1231 + C232222*term1232 + C233222*term1233;
        term2(1,11) = C211233*term1211 + C212233*term1212 + C213233*term1213 + C221233*term1221 + C222233*term1222 + C223233*term1223 + C231233*term1231 + C232233*term1232 + C233233*term1233;
        term2(1,12) = C211223*term1211 + C212223*term1212 + C213223*term1213 + C221223*term1221 + C222223*term1222 + C223223*term1223 + C231223*term1231 + C232223*term1232 + C233223*term1233;
        term2(1,13) = C211213*term1211 + C212213*term1212 + C213213*term1213 + C221213*term1221 + C222213*term1222 + C223213*term1223 + C231213*term1231 + C232213*term1232 + C233213*term1233;
        term2(1,14) = C211212*term1211 + C212212*term1212 + C213212*term1213 + C221212*term1221 + C222212*term1222 + C223212*term1223 + C231212*term1231 + C232212*term1232 + C233212*term1233;
        term2(1,15) = C211232*term1211 + C212232*term1212 + C213232*term1213 + C221232*term1221 + C222232*term1222 + C223232*term1223 + C231232*term1231 + C232232*term1232 + C233232*term1233;
        term2(1,16) = C211231*term1211 + C212231*term1212 + C213231*term1213 + C221231*term1221 + C222231*term1222 + C223231*term1223 + C231231*term1231 + C232231*term1232 + C233231*term1233;
        term2(1,17) = C211221*term1211 + C212221*term1212 + C213221*term1213 + C221221*term1221 + C222221*term1222 + C223221*term1223 + C231221*term1231 + C232221*term1232 + C233221*term1233;
        term2(1,18) = C211311*term1211 + C212311*term1212 + C213311*term1213 + C221311*term1221 + C222311*term1222 + C223311*term1223 + C231311*term1231 + C232311*term1232 + C233311*term1233;
        term2(1,19) = C211322*term1211 + C212322*term1212 + C213322*term1213 + C221322*term1221 + C222322*term1222 + C223322*term1223 + C231322*term1231 + C232322*term1232 + C233322*term1233;
        term2(1,20) = C211333*term1211 + C212333*term1212 + C213333*term1213 + C221333*term1221 + C222333*term1222 + C223333*term1223 + C231333*term1231 + C232333*term1232 + C233333*term1233;
        term2(1,21) = C211323*term1211 + C212323*term1212 + C213323*term1213 + C221323*term1221 + C222323*term1222 + C223323*term1223 + C231323*term1231 + C232323*term1232 + C233323*term1233;
        term2(1,22) = C211313*term1211 + C212313*term1212 + C213313*term1213 + C221313*term1221 + C222313*term1222 + C223313*term1223 + C231313*term1231 + C232313*term1232 + C233313*term1233;
        term2(1,23) = C211312*term1211 + C212312*term1212 + C213312*term1213 + C221312*term1221 + C222312*term1222 + C223312*term1223 + C231312*term1231 + C232312*term1232 + C233312*term1233;
        term2(1,24) = C211332*term1211 + C212332*term1212 + C213332*term1213 + C221332*term1221 + C222332*term1222 + C223332*term1223 + C231332*term1231 + C232332*term1232 + C233332*term1233;
        term2(1,25) = C211331*term1211 + C212331*term1212 + C213331*term1213 + C221331*term1221 + C222331*term1222 + C223331*term1223 + C231331*term1231 + C232331*term1232 + C233331*term1233;
        term2(1,26) = C211321*term1211 + C212321*term1212 + C213321*term1213 + C221321*term1221 + C222321*term1222 + C223321*term1223 + C231321*term1231 + C232321*term1232 + C233321*term1233;
        term2(2,0) = C311111*term1311 + C312111*term1312 + C313111*term1313 + C321111*term1321 + C322111*term1322 + C323111*term1323 + C331111*term1331 + C332111*term1332 + C333111*term1333;
        term2(2,1) = C311122*term1311 + C312122*term1312 + C313122*term1313 + C321122*term1321 + C322122*term1322 + C323122*term1323 + C331122*term1331 + C332122*term1332 + C333122*term1333;
        term2(2,2) = C311133*term1311 + C312133*term1312 + C313133*term1313 + C321133*term1321 + C322133*term1322 + C323133*term1323 + C331133*term1331 + C332133*term1332 + C333133*term1333;
        term2(2,3) = C311123*term1311 + C312123*term1312 + C313123*term1313 + C321123*term1321 + C322123*term1322 + C323123*term1323 + C331123*term1331 + C332123*term1332 + C333123*term1333;
        term2(2,4) = C311113*term1311 + C312113*term1312 + C313113*term1313 + C321113*term1321 + C322113*term1322 + C323113*term1323 + C331113*term1331 + C332113*term1332 + C333113*term1333;
        term2(2,5) = C311112*term1311 + C312112*term1312 + C313112*term1313 + C321112*term1321 + C322112*term1322 + C323112*term1323 + C331112*term1331 + C332112*term1332 + C333112*term1333;
        term2(2,6) = C311132*term1311 + C312132*term1312 + C313132*term1313 + C321132*term1321 + C322132*term1322 + C323132*term1323 + C331132*term1331 + C332132*term1332 + C333132*term1333;
        term2(2,7) = C311131*term1311 + C312131*term1312 + C313131*term1313 + C321131*term1321 + C322131*term1322 + C323131*term1323 + C331131*term1331 + C332131*term1332 + C333131*term1333;
        term2(2,8) = C311121*term1311 + C312121*term1312 + C313121*term1313 + C321121*term1321 + C322121*term1322 + C323121*term1323 + C331121*term1331 + C332121*term1332 + C333121*term1333;
        term2(2,9) = C311211*term1311 + C312211*term1312 + C313211*term1313 + C321211*term1321 + C322211*term1322 + C323211*term1323 + C331211*term1331 + C332211*term1332 + C333211*term1333;
        term2(2,10) = C311222*term1311 + C312222*term1312 + C313222*term1313 + C321222*term1321 + C322222*term1322 + C323222*term1323 + C331222*term1331 + C332222*term1332 + C333222*term1333;
        term2(2,11) = C311233*term1311 + C312233*term1312 + C313233*term1313 + C321233*term1321 + C322233*term1322 + C323233*term1323 + C331233*term1331 + C332233*term1332 + C333233*term1333;
        term2(2,12) = C311223*term1311 + C312223*term1312 + C313223*term1313 + C321223*term1321 + C322223*term1322 + C323223*term1323 + C331223*term1331 + C332223*term1332 + C333223*term1333;
        term2(2,13) = C311213*term1311 + C312213*term1312 + C313213*term1313 + C321213*term1321 + C322213*term1322 + C323213*term1323 + C331213*term1331 + C332213*term1332 + C333213*term1333;
        term2(2,14) = C311212*term1311 + C312212*term1312 + C313212*term1313 + C321212*term1321 + C322212*term1322 + C323212*term1323 + C331212*term1331 + C332212*term1332 + C333212*term1333;
        term2(2,15) = C311232*term1311 + C312232*term1312 + C313232*term1313 + C321232*term1321 + C322232*term1322 + C323232*term1323 + C331232*term1331 + C332232*term1332 + C333232*term1333;
        term2(2,16) = C311231*term1311 + C312231*term1312 + C313231*term1313 + C321231*term1321 + C322231*term1322 + C323231*term1323 + C331231*term1331 + C332231*term1332 + C333231*term1333;
        term2(2,17) = C311221*term1311 + C312221*term1312 + C313221*term1313 + C321221*term1321 + C322221*term1322 + C323221*term1323 + C331221*term1331 + C332221*term1332 + C333221*term1333;
        term2(2,18) = C311311*term1311 + C312311*term1312 + C313311*term1313 + C321311*term1321 + C322311*term1322 + C323311*term1323 + C331311*term1331 + C332311*term1332 + C333311*term1333;
        term2(2,19) = C311322*term1311 + C312322*term1312 + C313322*term1313 + C321322*term1321 + C322322*term1322 + C323322*term1323 + C331322*term1331 + C332322*term1332 + C333322*term1333;
        term2(2,20) = C311333*term1311 + C312333*term1312 + C313333*term1313 + C321333*term1321 + C322333*term1322 + C323333*term1323 + C331333*term1331 + C332333*term1332 + C333333*term1333;
        term2(2,21) = C311323*term1311 + C312323*term1312 + C313323*term1313 + C321323*term1321 + C322323*term1322 + C323323*term1323 + C331323*term1331 + C332323*term1332 + C333323*term1333;
        term2(2,22) = C311313*term1311 + C312313*term1312 + C313313*term1313 + C321313*term1321 + C322313*term1322 + C323313*term1323 + C331313*term1331 + C332313*term1332 + C333313*term1333;
        term2(2,23) = C311312*term1311 + C312312*term1312 + C313312*term1313 + C321312*term1321 + C322312*term1322 + C323312*term1323 + C331312*term1331 + C332312*term1332 + C333312*term1333;
        term2(2,24) = C311332*term1311 + C312332*term1312 + C313332*term1313 + C321332*term1321 + C322332*term1322 + C323332*term1323 + C331332*term1331 + C332332*term1332 + C333332*term1333;
        term2(2,25) = C311331*term1311 + C312331*term1312 + C313331*term1313 + C321331*term1321 + C322331*term1322 + C323331*term1323 + C331331*term1331 + C332331*term1332 + C333331*term1333;
        term2(2,26) = C311321*term1311 + C312321*term1312 + C313321*term1313 + C321321*term1321 + C322321*term1322 + C323321*term1323 + C331321*term1331 + C332321*term1332 + C333321*term1333;
        term2(3,0) = C211111*term1311 + C212111*term1312 + C213111*term1313 + C221111*term1321 + C222111*term1322 + C223111*term1323 + C231111*term1331 + C232111*term1332 + C233111*term1333;
        term2(3,1) = C211122*term1311 + C212122*term1312 + C213122*term1313 + C221122*term1321 + C222122*term1322 + C223122*term1323 + C231122*term1331 + C232122*term1332 + C233122*term1333;
        term2(3,2) = C211133*term1311 + C212133*term1312 + C213133*term1313 + C221133*term1321 + C222133*term1322 + C223133*term1323 + C231133*term1331 + C232133*term1332 + C233133*term1333;
        term2(3,3) = C211123*term1311 + C212123*term1312 + C213123*term1313 + C221123*term1321 + C222123*term1322 + C223123*term1323 + C231123*term1331 + C232123*term1332 + C233123*term1333;
        term2(3,4) = C211113*term1311 + C212113*term1312 + C213113*term1313 + C221113*term1321 + C222113*term1322 + C223113*term1323 + C231113*term1331 + C232113*term1332 + C233113*term1333;
        term2(3,5) = C211112*term1311 + C212112*term1312 + C213112*term1313 + C221112*term1321 + C222112*term1322 + C223112*term1323 + C231112*term1331 + C232112*term1332 + C233112*term1333;
        term2(3,6) = C211132*term1311 + C212132*term1312 + C213132*term1313 + C221132*term1321 + C222132*term1322 + C223132*term1323 + C231132*term1331 + C232132*term1332 + C233132*term1333;
        term2(3,7) = C211131*term1311 + C212131*term1312 + C213131*term1313 + C221131*term1321 + C222131*term1322 + C223131*term1323 + C231131*term1331 + C232131*term1332 + C233131*term1333;
        term2(3,8) = C211121*term1311 + C212121*term1312 + C213121*term1313 + C221121*term1321 + C222121*term1322 + C223121*term1323 + C231121*term1331 + C232121*term1332 + C233121*term1333;
        term2(3,9) = C211211*term1311 + C212211*term1312 + C213211*term1313 + C221211*term1321 + C222211*term1322 + C223211*term1323 + C231211*term1331 + C232211*term1332 + C233211*term1333;
        term2(3,10) = C211222*term1311 + C212222*term1312 + C213222*term1313 + C221222*term1321 + C222222*term1322 + C223222*term1323 + C231222*term1331 + C232222*term1332 + C233222*term1333;
        term2(3,11) = C211233*term1311 + C212233*term1312 + C213233*term1313 + C221233*term1321 + C222233*term1322 + C223233*term1323 + C231233*term1331 + C232233*term1332 + C233233*term1333;
        term2(3,12) = C211223*term1311 + C212223*term1312 + C213223*term1313 + C221223*term1321 + C222223*term1322 + C223223*term1323 + C231223*term1331 + C232223*term1332 + C233223*term1333;
        term2(3,13) = C211213*term1311 + C212213*term1312 + C213213*term1313 + C221213*term1321 + C222213*term1322 + C223213*term1323 + C231213*term1331 + C232213*term1332 + C233213*term1333;
        term2(3,14) = C211212*term1311 + C212212*term1312 + C213212*term1313 + C221212*term1321 + C222212*term1322 + C223212*term1323 + C231212*term1331 + C232212*term1332 + C233212*term1333;
        term2(3,15) = C211232*term1311 + C212232*term1312 + C213232*term1313 + C221232*term1321 + C222232*term1322 + C223232*term1323 + C231232*term1331 + C232232*term1332 + C233232*term1333;
        term2(3,16) = C211231*term1311 + C212231*term1312 + C213231*term1313 + C221231*term1321 + C222231*term1322 + C223231*term1323 + C231231*term1331 + C232231*term1332 + C233231*term1333;
        term2(3,17) = C211221*term1311 + C212221*term1312 + C213221*term1313 + C221221*term1321 + C222221*term1322 + C223221*term1323 + C231221*term1331 + C232221*term1332 + C233221*term1333;
        term2(3,18) = C211311*term1311 + C212311*term1312 + C213311*term1313 + C221311*term1321 + C222311*term1322 + C223311*term1323 + C231311*term1331 + C232311*term1332 + C233311*term1333;
        term2(3,19) = C211322*term1311 + C212322*term1312 + C213322*term1313 + C221322*term1321 + C222322*term1322 + C223322*term1323 + C231322*term1331 + C232322*term1332 + C233322*term1333;
        term2(3,20) = C211333*term1311 + C212333*term1312 + C213333*term1313 + C221333*term1321 + C222333*term1322 + C223333*term1323 + C231333*term1331 + C232333*term1332 + C233333*term1333;
        term2(3,21) = C211323*term1311 + C212323*term1312 + C213323*term1313 + C221323*term1321 + C222323*term1322 + C223323*term1323 + C231323*term1331 + C232323*term1332 + C233323*term1333;
        term2(3,22) = C211313*term1311 + C212313*term1312 + C213313*term1313 + C221313*term1321 + C222313*term1322 + C223313*term1323 + C231313*term1331 + C232313*term1332 + C233313*term1333;
        term2(3,23) = C211312*term1311 + C212312*term1312 + C213312*term1313 + C221312*term1321 + C222312*term1322 + C223312*term1323 + C231312*term1331 + C232312*term1332 + C233312*term1333;
        term2(3,24) = C211332*term1311 + C212332*term1312 + C213332*term1313 + C221332*term1321 + C222332*term1322 + C223332*term1323 + C231332*term1331 + C232332*term1332 + C233332*term1333;
        term2(3,25) = C211331*term1311 + C212331*term1312 + C213331*term1313 + C221331*term1321 + C222331*term1322 + C223331*term1323 + C231331*term1331 + C232331*term1332 + C233331*term1333;
        term2(3,26) = C211321*term1311 + C212321*term1312 + C213321*term1313 + C221321*term1321 + C222321*term1322 + C223321*term1323 + C231321*term1331 + C232321*term1332 + C233321*term1333;
        term2(4,0) = C111111*term1311 + C112111*term1312 + C113111*term1313 + C121111*term1321 + C122111*term1322 + C123111*term1323 + C131111*term1331 + C132111*term1332 + C133111*term1333;
        term2(4,1) = C111122*term1311 + C112122*term1312 + C113122*term1313 + C121122*term1321 + C122122*term1322 + C123122*term1323 + C131122*term1331 + C132122*term1332 + C133122*term1333;
        term2(4,2) = C111133*term1311 + C112133*term1312 + C113133*term1313 + C121133*term1321 + C122133*term1322 + C123133*term1323 + C131133*term1331 + C132133*term1332 + C133133*term1333;
        term2(4,3) = C111123*term1311 + C112123*term1312 + C113123*term1313 + C121123*term1321 + C122123*term1322 + C123123*term1323 + C131123*term1331 + C132123*term1332 + C133123*term1333;
        term2(4,4) = C111113*term1311 + C112113*term1312 + C113113*term1313 + C121113*term1321 + C122113*term1322 + C123113*term1323 + C131113*term1331 + C132113*term1332 + C133113*term1333;
        term2(4,5) = C111112*term1311 + C112112*term1312 + C113112*term1313 + C121112*term1321 + C122112*term1322 + C123112*term1323 + C131112*term1331 + C132112*term1332 + C133112*term1333;
        term2(4,6) = C111132*term1311 + C112132*term1312 + C113132*term1313 + C121132*term1321 + C122132*term1322 + C123132*term1323 + C131132*term1331 + C132132*term1332 + C133132*term1333;
        term2(4,7) = C111131*term1311 + C112131*term1312 + C113131*term1313 + C121131*term1321 + C122131*term1322 + C123131*term1323 + C131131*term1331 + C132131*term1332 + C133131*term1333;
        term2(4,8) = C111121*term1311 + C112121*term1312 + C113121*term1313 + C121121*term1321 + C122121*term1322 + C123121*term1323 + C131121*term1331 + C132121*term1332 + C133121*term1333;
        term2(4,9) = C111211*term1311 + C112211*term1312 + C113211*term1313 + C121211*term1321 + C122211*term1322 + C123211*term1323 + C131211*term1331 + C132211*term1332 + C133211*term1333;
        term2(4,10) = C111222*term1311 + C112222*term1312 + C113222*term1313 + C121222*term1321 + C122222*term1322 + C123222*term1323 + C131222*term1331 + C132222*term1332 + C133222*term1333;
        term2(4,11) = C111233*term1311 + C112233*term1312 + C113233*term1313 + C121233*term1321 + C122233*term1322 + C123233*term1323 + C131233*term1331 + C132233*term1332 + C133233*term1333;
        term2(4,12) = C111223*term1311 + C112223*term1312 + C113223*term1313 + C121223*term1321 + C122223*term1322 + C123223*term1323 + C131223*term1331 + C132223*term1332 + C133223*term1333;
        term2(4,13) = C111213*term1311 + C112213*term1312 + C113213*term1313 + C121213*term1321 + C122213*term1322 + C123213*term1323 + C131213*term1331 + C132213*term1332 + C133213*term1333;
        term2(4,14) = C111212*term1311 + C112212*term1312 + C113212*term1313 + C121212*term1321 + C122212*term1322 + C123212*term1323 + C131212*term1331 + C132212*term1332 + C133212*term1333;
        term2(4,15) = C111232*term1311 + C112232*term1312 + C113232*term1313 + C121232*term1321 + C122232*term1322 + C123232*term1323 + C131232*term1331 + C132232*term1332 + C133232*term1333;
        term2(4,16) = C111231*term1311 + C112231*term1312 + C113231*term1313 + C121231*term1321 + C122231*term1322 + C123231*term1323 + C131231*term1331 + C132231*term1332 + C133231*term1333;
        term2(4,17) = C111221*term1311 + C112221*term1312 + C113221*term1313 + C121221*term1321 + C122221*term1322 + C123221*term1323 + C131221*term1331 + C132221*term1332 + C133221*term1333;
        term2(4,18) = C111311*term1311 + C112311*term1312 + C113311*term1313 + C121311*term1321 + C122311*term1322 + C123311*term1323 + C131311*term1331 + C132311*term1332 + C133311*term1333;
        term2(4,19) = C111322*term1311 + C112322*term1312 + C113322*term1313 + C121322*term1321 + C122322*term1322 + C123322*term1323 + C131322*term1331 + C132322*term1332 + C133322*term1333;
        term2(4,20) = C111333*term1311 + C112333*term1312 + C113333*term1313 + C121333*term1321 + C122333*term1322 + C123333*term1323 + C131333*term1331 + C132333*term1332 + C133333*term1333;
        term2(4,21) = C111323*term1311 + C112323*term1312 + C113323*term1313 + C121323*term1321 + C122323*term1322 + C123323*term1323 + C131323*term1331 + C132323*term1332 + C133323*term1333;
        term2(4,22) = C111313*term1311 + C112313*term1312 + C113313*term1313 + C121313*term1321 + C122313*term1322 + C123313*term1323 + C131313*term1331 + C132313*term1332 + C133313*term1333;
        term2(4,23) = C111312*term1311 + C112312*term1312 + C113312*term1313 + C121312*term1321 + C122312*term1322 + C123312*term1323 + C131312*term1331 + C132312*term1332 + C133312*term1333;
        term2(4,24) = C111332*term1311 + C112332*term1312 + C113332*term1313 + C121332*term1321 + C122332*term1322 + C123332*term1323 + C131332*term1331 + C132332*term1332 + C133332*term1333;
        term2(4,25) = C111331*term1311 + C112331*term1312 + C113331*term1313 + C121331*term1321 + C122331*term1322 + C123331*term1323 + C131331*term1331 + C132331*term1332 + C133331*term1333;
        term2(4,26) = C111321*term1311 + C112321*term1312 + C113321*term1313 + C121321*term1321 + C122321*term1322 + C123321*term1323 + C131321*term1331 + C132321*term1332 + C133321*term1333;
        term2(5,0) = C111111*term1211 + C112111*term1212 + C113111*term1213 + C121111*term1221 + C122111*term1222 + C123111*term1223 + C131111*term1231 + C132111*term1232 + C133111*term1233;
        term2(5,1) = C111122*term1211 + C112122*term1212 + C113122*term1213 + C121122*term1221 + C122122*term1222 + C123122*term1223 + C131122*term1231 + C132122*term1232 + C133122*term1233;
        term2(5,2) = C111133*term1211 + C112133*term1212 + C113133*term1213 + C121133*term1221 + C122133*term1222 + C123133*term1223 + C131133*term1231 + C132133*term1232 + C133133*term1233;
        term2(5,3) = C111123*term1211 + C112123*term1212 + C113123*term1213 + C121123*term1221 + C122123*term1222 + C123123*term1223 + C131123*term1231 + C132123*term1232 + C133123*term1233;
        term2(5,4) = C111113*term1211 + C112113*term1212 + C113113*term1213 + C121113*term1221 + C122113*term1222 + C123113*term1223 + C131113*term1231 + C132113*term1232 + C133113*term1233;
        term2(5,5) = C111112*term1211 + C112112*term1212 + C113112*term1213 + C121112*term1221 + C122112*term1222 + C123112*term1223 + C131112*term1231 + C132112*term1232 + C133112*term1233;
        term2(5,6) = C111132*term1211 + C112132*term1212 + C113132*term1213 + C121132*term1221 + C122132*term1222 + C123132*term1223 + C131132*term1231 + C132132*term1232 + C133132*term1233;
        term2(5,7) = C111131*term1211 + C112131*term1212 + C113131*term1213 + C121131*term1221 + C122131*term1222 + C123131*term1223 + C131131*term1231 + C132131*term1232 + C133131*term1233;
        term2(5,8) = C111121*term1211 + C112121*term1212 + C113121*term1213 + C121121*term1221 + C122121*term1222 + C123121*term1223 + C131121*term1231 + C132121*term1232 + C133121*term1233;
        term2(5,9) = C111211*term1211 + C112211*term1212 + C113211*term1213 + C121211*term1221 + C122211*term1222 + C123211*term1223 + C131211*term1231 + C132211*term1232 + C133211*term1233;
        term2(5,10) = C111222*term1211 + C112222*term1212 + C113222*term1213 + C121222*term1221 + C122222*term1222 + C123222*term1223 + C131222*term1231 + C132222*term1232 + C133222*term1233;
        term2(5,11) = C111233*term1211 + C112233*term1212 + C113233*term1213 + C121233*term1221 + C122233*term1222 + C123233*term1223 + C131233*term1231 + C132233*term1232 + C133233*term1233;
        term2(5,12) = C111223*term1211 + C112223*term1212 + C113223*term1213 + C121223*term1221 + C122223*term1222 + C123223*term1223 + C131223*term1231 + C132223*term1232 + C133223*term1233;
        term2(5,13) = C111213*term1211 + C112213*term1212 + C113213*term1213 + C121213*term1221 + C122213*term1222 + C123213*term1223 + C131213*term1231 + C132213*term1232 + C133213*term1233;
        term2(5,14) = C111212*term1211 + C112212*term1212 + C113212*term1213 + C121212*term1221 + C122212*term1222 + C123212*term1223 + C131212*term1231 + C132212*term1232 + C133212*term1233;
        term2(5,15) = C111232*term1211 + C112232*term1212 + C113232*term1213 + C121232*term1221 + C122232*term1222 + C123232*term1223 + C131232*term1231 + C132232*term1232 + C133232*term1233;
        term2(5,16) = C111231*term1211 + C112231*term1212 + C113231*term1213 + C121231*term1221 + C122231*term1222 + C123231*term1223 + C131231*term1231 + C132231*term1232 + C133231*term1233;
        term2(5,17) = C111221*term1211 + C112221*term1212 + C113221*term1213 + C121221*term1221 + C122221*term1222 + C123221*term1223 + C131221*term1231 + C132221*term1232 + C133221*term1233;
        term2(5,18) = C111311*term1211 + C112311*term1212 + C113311*term1213 + C121311*term1221 + C122311*term1222 + C123311*term1223 + C131311*term1231 + C132311*term1232 + C133311*term1233;
        term2(5,19) = C111322*term1211 + C112322*term1212 + C113322*term1213 + C121322*term1221 + C122322*term1222 + C123322*term1223 + C131322*term1231 + C132322*term1232 + C133322*term1233;
        term2(5,20) = C111333*term1211 + C112333*term1212 + C113333*term1213 + C121333*term1221 + C122333*term1222 + C123333*term1223 + C131333*term1231 + C132333*term1232 + C133333*term1233;
        term2(5,21) = C111323*term1211 + C112323*term1212 + C113323*term1213 + C121323*term1221 + C122323*term1222 + C123323*term1223 + C131323*term1231 + C132323*term1232 + C133323*term1233;
        term2(5,22) = C111313*term1211 + C112313*term1212 + C113313*term1213 + C121313*term1221 + C122313*term1222 + C123313*term1223 + C131313*term1231 + C132313*term1232 + C133313*term1233;
        term2(5,23) = C111312*term1211 + C112312*term1212 + C113312*term1213 + C121312*term1221 + C122312*term1222 + C123312*term1223 + C131312*term1231 + C132312*term1232 + C133312*term1233;
        term2(5,24) = C111332*term1211 + C112332*term1212 + C113332*term1213 + C121332*term1221 + C122332*term1222 + C123332*term1223 + C131332*term1231 + C132332*term1232 + C133332*term1233;
        term2(5,25) = C111331*term1211 + C112331*term1212 + C113331*term1213 + C121331*term1221 + C122331*term1222 + C123331*term1223 + C131331*term1231 + C132331*term1232 + C133331*term1233;
        term2(5,26) = C111321*term1211 + C112321*term1212 + C113321*term1213 + C121321*term1221 + C122321*term1222 + C123321*term1223 + C131321*term1231 + C132321*term1232 + C133321*term1233;
        term2(6,0) = C311111*term1211 + C312111*term1212 + C313111*term1213 + C321111*term1221 + C322111*term1222 + C323111*term1223 + C331111*term1231 + C332111*term1232 + C333111*term1233;
        term2(6,1) = C311122*term1211 + C312122*term1212 + C313122*term1213 + C321122*term1221 + C322122*term1222 + C323122*term1223 + C331122*term1231 + C332122*term1232 + C333122*term1233;
        term2(6,2) = C311133*term1211 + C312133*term1212 + C313133*term1213 + C321133*term1221 + C322133*term1222 + C323133*term1223 + C331133*term1231 + C332133*term1232 + C333133*term1233;
        term2(6,3) = C311123*term1211 + C312123*term1212 + C313123*term1213 + C321123*term1221 + C322123*term1222 + C323123*term1223 + C331123*term1231 + C332123*term1232 + C333123*term1233;
        term2(6,4) = C311113*term1211 + C312113*term1212 + C313113*term1213 + C321113*term1221 + C322113*term1222 + C323113*term1223 + C331113*term1231 + C332113*term1232 + C333113*term1233;
        term2(6,5) = C311112*term1211 + C312112*term1212 + C313112*term1213 + C321112*term1221 + C322112*term1222 + C323112*term1223 + C331112*term1231 + C332112*term1232 + C333112*term1233;
        term2(6,6) = C311132*term1211 + C312132*term1212 + C313132*term1213 + C321132*term1221 + C322132*term1222 + C323132*term1223 + C331132*term1231 + C332132*term1232 + C333132*term1233;
        term2(6,7) = C311131*term1211 + C312131*term1212 + C313131*term1213 + C321131*term1221 + C322131*term1222 + C323131*term1223 + C331131*term1231 + C332131*term1232 + C333131*term1233;
        term2(6,8) = C311121*term1211 + C312121*term1212 + C313121*term1213 + C321121*term1221 + C322121*term1222 + C323121*term1223 + C331121*term1231 + C332121*term1232 + C333121*term1233;
        term2(6,9) = C311211*term1211 + C312211*term1212 + C313211*term1213 + C321211*term1221 + C322211*term1222 + C323211*term1223 + C331211*term1231 + C332211*term1232 + C333211*term1233;
        term2(6,10) = C311222*term1211 + C312222*term1212 + C313222*term1213 + C321222*term1221 + C322222*term1222 + C323222*term1223 + C331222*term1231 + C332222*term1232 + C333222*term1233;
        term2(6,11) = C311233*term1211 + C312233*term1212 + C313233*term1213 + C321233*term1221 + C322233*term1222 + C323233*term1223 + C331233*term1231 + C332233*term1232 + C333233*term1233;
        term2(6,12) = C311223*term1211 + C312223*term1212 + C313223*term1213 + C321223*term1221 + C322223*term1222 + C323223*term1223 + C331223*term1231 + C332223*term1232 + C333223*term1233;
        term2(6,13) = C311213*term1211 + C312213*term1212 + C313213*term1213 + C321213*term1221 + C322213*term1222 + C323213*term1223 + C331213*term1231 + C332213*term1232 + C333213*term1233;
        term2(6,14) = C311212*term1211 + C312212*term1212 + C313212*term1213 + C321212*term1221 + C322212*term1222 + C323212*term1223 + C331212*term1231 + C332212*term1232 + C333212*term1233;
        term2(6,15) = C311232*term1211 + C312232*term1212 + C313232*term1213 + C321232*term1221 + C322232*term1222 + C323232*term1223 + C331232*term1231 + C332232*term1232 + C333232*term1233;
        term2(6,16) = C311231*term1211 + C312231*term1212 + C313231*term1213 + C321231*term1221 + C322231*term1222 + C323231*term1223 + C331231*term1231 + C332231*term1232 + C333231*term1233;
        term2(6,17) = C311221*term1211 + C312221*term1212 + C313221*term1213 + C321221*term1221 + C322221*term1222 + C323221*term1223 + C331221*term1231 + C332221*term1232 + C333221*term1233;
        term2(6,18) = C311311*term1211 + C312311*term1212 + C313311*term1213 + C321311*term1221 + C322311*term1222 + C323311*term1223 + C331311*term1231 + C332311*term1232 + C333311*term1233;
        term2(6,19) = C311322*term1211 + C312322*term1212 + C313322*term1213 + C321322*term1221 + C322322*term1222 + C323322*term1223 + C331322*term1231 + C332322*term1232 + C333322*term1233;
        term2(6,20) = C311333*term1211 + C312333*term1212 + C313333*term1213 + C321333*term1221 + C322333*term1222 + C323333*term1223 + C331333*term1231 + C332333*term1232 + C333333*term1233;
        term2(6,21) = C311323*term1211 + C312323*term1212 + C313323*term1213 + C321323*term1221 + C322323*term1222 + C323323*term1223 + C331323*term1231 + C332323*term1232 + C333323*term1233;
        term2(6,22) = C311313*term1211 + C312313*term1212 + C313313*term1213 + C321313*term1221 + C322313*term1222 + C323313*term1223 + C331313*term1231 + C332313*term1232 + C333313*term1233;
        term2(6,23) = C311312*term1211 + C312312*term1212 + C313312*term1213 + C321312*term1221 + C322312*term1222 + C323312*term1223 + C331312*term1231 + C332312*term1232 + C333312*term1233;
        term2(6,24) = C311332*term1211 + C312332*term1212 + C313332*term1213 + C321332*term1221 + C322332*term1222 + C323332*term1223 + C331332*term1231 + C332332*term1232 + C333332*term1233;
        term2(6,25) = C311331*term1211 + C312331*term1212 + C313331*term1213 + C321331*term1221 + C322331*term1222 + C323331*term1223 + C331331*term1231 + C332331*term1232 + C333331*term1233;
        term2(6,26) = C311321*term1211 + C312321*term1212 + C313321*term1213 + C321321*term1221 + C322321*term1222 + C323321*term1223 + C331321*term1231 + C332321*term1232 + C333321*term1233;
        term2(7,0) = C311111*term1111 + C312111*term1112 + C313111*term1113 + C321111*term1121 + C322111*term1122 + C323111*term1123 + C331111*term1131 + C332111*term1132 + C333111*term1133;
        term2(7,1) = C311122*term1111 + C312122*term1112 + C313122*term1113 + C321122*term1121 + C322122*term1122 + C323122*term1123 + C331122*term1131 + C332122*term1132 + C333122*term1133;
        term2(7,2) = C311133*term1111 + C312133*term1112 + C313133*term1113 + C321133*term1121 + C322133*term1122 + C323133*term1123 + C331133*term1131 + C332133*term1132 + C333133*term1133;
        term2(7,3) = C311123*term1111 + C312123*term1112 + C313123*term1113 + C321123*term1121 + C322123*term1122 + C323123*term1123 + C331123*term1131 + C332123*term1132 + C333123*term1133;
        term2(7,4) = C311113*term1111 + C312113*term1112 + C313113*term1113 + C321113*term1121 + C322113*term1122 + C323113*term1123 + C331113*term1131 + C332113*term1132 + C333113*term1133;
        term2(7,5) = C311112*term1111 + C312112*term1112 + C313112*term1113 + C321112*term1121 + C322112*term1122 + C323112*term1123 + C331112*term1131 + C332112*term1132 + C333112*term1133;
        term2(7,6) = C311132*term1111 + C312132*term1112 + C313132*term1113 + C321132*term1121 + C322132*term1122 + C323132*term1123 + C331132*term1131 + C332132*term1132 + C333132*term1133;
        term2(7,7) = C311131*term1111 + C312131*term1112 + C313131*term1113 + C321131*term1121 + C322131*term1122 + C323131*term1123 + C331131*term1131 + C332131*term1132 + C333131*term1133;
        term2(7,8) = C311121*term1111 + C312121*term1112 + C313121*term1113 + C321121*term1121 + C322121*term1122 + C323121*term1123 + C331121*term1131 + C332121*term1132 + C333121*term1133;
        term2(7,9) = C311211*term1111 + C312211*term1112 + C313211*term1113 + C321211*term1121 + C322211*term1122 + C323211*term1123 + C331211*term1131 + C332211*term1132 + C333211*term1133;
        term2(7,10) = C311222*term1111 + C312222*term1112 + C313222*term1113 + C321222*term1121 + C322222*term1122 + C323222*term1123 + C331222*term1131 + C332222*term1132 + C333222*term1133;
        term2(7,11) = C311233*term1111 + C312233*term1112 + C313233*term1113 + C321233*term1121 + C322233*term1122 + C323233*term1123 + C331233*term1131 + C332233*term1132 + C333233*term1133;
        term2(7,12) = C311223*term1111 + C312223*term1112 + C313223*term1113 + C321223*term1121 + C322223*term1122 + C323223*term1123 + C331223*term1131 + C332223*term1132 + C333223*term1133;
        term2(7,13) = C311213*term1111 + C312213*term1112 + C313213*term1113 + C321213*term1121 + C322213*term1122 + C323213*term1123 + C331213*term1131 + C332213*term1132 + C333213*term1133;
        term2(7,14) = C311212*term1111 + C312212*term1112 + C313212*term1113 + C321212*term1121 + C322212*term1122 + C323212*term1123 + C331212*term1131 + C332212*term1132 + C333212*term1133;
        term2(7,15) = C311232*term1111 + C312232*term1112 + C313232*term1113 + C321232*term1121 + C322232*term1122 + C323232*term1123 + C331232*term1131 + C332232*term1132 + C333232*term1133;
        term2(7,16) = C311231*term1111 + C312231*term1112 + C313231*term1113 + C321231*term1121 + C322231*term1122 + C323231*term1123 + C331231*term1131 + C332231*term1132 + C333231*term1133;
        term2(7,17) = C311221*term1111 + C312221*term1112 + C313221*term1113 + C321221*term1121 + C322221*term1122 + C323221*term1123 + C331221*term1131 + C332221*term1132 + C333221*term1133;
        term2(7,18) = C311311*term1111 + C312311*term1112 + C313311*term1113 + C321311*term1121 + C322311*term1122 + C323311*term1123 + C331311*term1131 + C332311*term1132 + C333311*term1133;
        term2(7,19) = C311322*term1111 + C312322*term1112 + C313322*term1113 + C321322*term1121 + C322322*term1122 + C323322*term1123 + C331322*term1131 + C332322*term1132 + C333322*term1133;
        term2(7,20) = C311333*term1111 + C312333*term1112 + C313333*term1113 + C321333*term1121 + C322333*term1122 + C323333*term1123 + C331333*term1131 + C332333*term1132 + C333333*term1133;
        term2(7,21) = C311323*term1111 + C312323*term1112 + C313323*term1113 + C321323*term1121 + C322323*term1122 + C323323*term1123 + C331323*term1131 + C332323*term1132 + C333323*term1133;
        term2(7,22) = C311313*term1111 + C312313*term1112 + C313313*term1113 + C321313*term1121 + C322313*term1122 + C323313*term1123 + C331313*term1131 + C332313*term1132 + C333313*term1133;
        term2(7,23) = C311312*term1111 + C312312*term1112 + C313312*term1113 + C321312*term1121 + C322312*term1122 + C323312*term1123 + C331312*term1131 + C332312*term1132 + C333312*term1133;
        term2(7,24) = C311332*term1111 + C312332*term1112 + C313332*term1113 + C321332*term1121 + C322332*term1122 + C323332*term1123 + C331332*term1131 + C332332*term1132 + C333332*term1133;
        term2(7,25) = C311331*term1111 + C312331*term1112 + C313331*term1113 + C321331*term1121 + C322331*term1122 + C323331*term1123 + C331331*term1131 + C332331*term1132 + C333331*term1133;
        term2(7,26) = C311321*term1111 + C312321*term1112 + C313321*term1113 + C321321*term1121 + C322321*term1122 + C323321*term1123 + C331321*term1131 + C332321*term1132 + C333321*term1133;
        term2(8,0) = C211111*term1111 + C212111*term1112 + C213111*term1113 + C221111*term1121 + C222111*term1122 + C223111*term1123 + C231111*term1131 + C232111*term1132 + C233111*term1133;
        term2(8,1) = C211122*term1111 + C212122*term1112 + C213122*term1113 + C221122*term1121 + C222122*term1122 + C223122*term1123 + C231122*term1131 + C232122*term1132 + C233122*term1133;
        term2(8,2) = C211133*term1111 + C212133*term1112 + C213133*term1113 + C221133*term1121 + C222133*term1122 + C223133*term1123 + C231133*term1131 + C232133*term1132 + C233133*term1133;
        term2(8,3) = C211123*term1111 + C212123*term1112 + C213123*term1113 + C221123*term1121 + C222123*term1122 + C223123*term1123 + C231123*term1131 + C232123*term1132 + C233123*term1133;
        term2(8,4) = C211113*term1111 + C212113*term1112 + C213113*term1113 + C221113*term1121 + C222113*term1122 + C223113*term1123 + C231113*term1131 + C232113*term1132 + C233113*term1133;
        term2(8,5) = C211112*term1111 + C212112*term1112 + C213112*term1113 + C221112*term1121 + C222112*term1122 + C223112*term1123 + C231112*term1131 + C232112*term1132 + C233112*term1133;
        term2(8,6) = C211132*term1111 + C212132*term1112 + C213132*term1113 + C221132*term1121 + C222132*term1122 + C223132*term1123 + C231132*term1131 + C232132*term1132 + C233132*term1133;
        term2(8,7) = C211131*term1111 + C212131*term1112 + C213131*term1113 + C221131*term1121 + C222131*term1122 + C223131*term1123 + C231131*term1131 + C232131*term1132 + C233131*term1133;
        term2(8,8) = C211121*term1111 + C212121*term1112 + C213121*term1113 + C221121*term1121 + C222121*term1122 + C223121*term1123 + C231121*term1131 + C232121*term1132 + C233121*term1133;
        term2(8,9) = C211211*term1111 + C212211*term1112 + C213211*term1113 + C221211*term1121 + C222211*term1122 + C223211*term1123 + C231211*term1131 + C232211*term1132 + C233211*term1133;
        term2(8,10) = C211222*term1111 + C212222*term1112 + C213222*term1113 + C221222*term1121 + C222222*term1122 + C223222*term1123 + C231222*term1131 + C232222*term1132 + C233222*term1133;
        term2(8,11) = C211233*term1111 + C212233*term1112 + C213233*term1113 + C221233*term1121 + C222233*term1122 + C223233*term1123 + C231233*term1131 + C232233*term1132 + C233233*term1133;
        term2(8,12) = C211223*term1111 + C212223*term1112 + C213223*term1113 + C221223*term1121 + C222223*term1122 + C223223*term1123 + C231223*term1131 + C232223*term1132 + C233223*term1133;
        term2(8,13) = C211213*term1111 + C212213*term1112 + C213213*term1113 + C221213*term1121 + C222213*term1122 + C223213*term1123 + C231213*term1131 + C232213*term1132 + C233213*term1133;
        term2(8,14) = C211212*term1111 + C212212*term1112 + C213212*term1113 + C221212*term1121 + C222212*term1122 + C223212*term1123 + C231212*term1131 + C232212*term1132 + C233212*term1133;
        term2(8,15) = C211232*term1111 + C212232*term1112 + C213232*term1113 + C221232*term1121 + C222232*term1122 + C223232*term1123 + C231232*term1131 + C232232*term1132 + C233232*term1133;
        term2(8,16) = C211231*term1111 + C212231*term1112 + C213231*term1113 + C221231*term1121 + C222231*term1122 + C223231*term1123 + C231231*term1131 + C232231*term1132 + C233231*term1133;
        term2(8,17) = C211221*term1111 + C212221*term1112 + C213221*term1113 + C221221*term1121 + C222221*term1122 + C223221*term1123 + C231221*term1131 + C232221*term1132 + C233221*term1133;
        term2(8,18) = C211311*term1111 + C212311*term1112 + C213311*term1113 + C221311*term1121 + C222311*term1122 + C223311*term1123 + C231311*term1131 + C232311*term1132 + C233311*term1133;
        term2(8,19) = C211322*term1111 + C212322*term1112 + C213322*term1113 + C221322*term1121 + C222322*term1122 + C223322*term1123 + C231322*term1131 + C232322*term1132 + C233322*term1133;
        term2(8,20) = C211333*term1111 + C212333*term1112 + C213333*term1113 + C221333*term1121 + C222333*term1122 + C223333*term1123 + C231333*term1131 + C232333*term1132 + C233333*term1133;
        term2(8,21) = C211323*term1111 + C212323*term1112 + C213323*term1113 + C221323*term1121 + C222323*term1122 + C223323*term1123 + C231323*term1131 + C232323*term1132 + C233323*term1133;
        term2(8,22) = C211313*term1111 + C212313*term1112 + C213313*term1113 + C221313*term1121 + C222313*term1122 + C223313*term1123 + C231313*term1131 + C232313*term1132 + C233313*term1133;
        term2(8,23) = C211312*term1111 + C212312*term1112 + C213312*term1113 + C221312*term1121 + C222312*term1122 + C223312*term1123 + C231312*term1131 + C232312*term1132 + C233312*term1133;
        term2(8,24) = C211332*term1111 + C212332*term1112 + C213332*term1113 + C221332*term1121 + C222332*term1122 + C223332*term1123 + C231332*term1131 + C232332*term1132 + C233332*term1133;
        term2(8,25) = C211331*term1111 + C212331*term1112 + C213331*term1113 + C221331*term1121 + C222331*term1122 + C223331*term1123 + C231331*term1131 + C232331*term1132 + C233331*term1133;
        term2(8,26) = C211321*term1111 + C212321*term1112 + C213321*term1113 + C221321*term1121 + C222321*term1122 + C223321*term1123 + C231321*term1131 + C232321*term1132 + C233321*term1133;
        return;

    }
    
    void compute_dPK2dGamma_term3(const Vector_27 &term1, const Matrix_3x3 &RCGinv, Matrix_9x27 &term3){
        /*!==================================
        |    compute_dPK2dGamma_term3    |
        ==================================
        
        Compute term3 of dPK2dGamma. Note, also the same as dSIGMAdGamma if symmetrized.
        
        */
        
        //Extract term1
        double term1111 = term1(0);
        double term1122 = term1(1);
        double term1133 = term1(2);
        double term1123 = term1(3);
        double term1113 = term1(4);
        double term1112 = term1(5);
        double term1132 = term1(6);
        double term1131 = term1(7);
        double term1121 = term1(8);
        double term1211 = term1(9);
        double term1222 = term1(10);
        double term1233 = term1(11);
        double term1223 = term1(12);
        double term1213 = term1(13);
        double term1212 = term1(14);
        double term1232 = term1(15);
        double term1231 = term1(16);
        double term1221 = term1(17);
        double term1311 = term1(18);
        double term1322 = term1(19);
        double term1333 = term1(20);
        double term1323 = term1(21);
        double term1313 = term1(22);
        double term1312 = term1(23);
        double term1332 = term1(24);
        double term1331 = term1(25);
        double term1321 = term1(26);

        //Extract RCGinv
        double RCGinv11 = RCGinv(0,0);
        double RCGinv12 = RCGinv(0,1);
        double RCGinv13 = RCGinv(0,2);
        double RCGinv21 = RCGinv(1,0);
        double RCGinv22 = RCGinv(1,1);
        double RCGinv23 = RCGinv(1,2);
        double RCGinv31 = RCGinv(2,0);
        double RCGinv32 = RCGinv(2,1);
        double RCGinv33 = RCGinv(2,2);

         //Assemble term3
        term3(0,0) = RCGinv11*term1111;
        term3(0,1) = RCGinv11*term1122;
        term3(0,2) = RCGinv11*term1133;
        term3(0,3) = RCGinv11*term1123;
        term3(0,4) = RCGinv11*term1113;
        term3(0,5) = RCGinv11*term1112;
        term3(0,6) = RCGinv11*term1132;
        term3(0,7) = RCGinv11*term1131;
        term3(0,8) = RCGinv11*term1121;
        term3(0,9) = RCGinv12*term1111;
        term3(0,10) = RCGinv12*term1122;
        term3(0,11) = RCGinv12*term1133;
        term3(0,12) = RCGinv12*term1123;
        term3(0,13) = RCGinv12*term1113;
        term3(0,14) = RCGinv12*term1112;
        term3(0,15) = RCGinv12*term1132;
        term3(0,16) = RCGinv12*term1131;
        term3(0,17) = RCGinv12*term1121;
        term3(0,18) = RCGinv13*term1111;
        term3(0,19) = RCGinv13*term1122;
        term3(0,20) = RCGinv13*term1133;
        term3(0,21) = RCGinv13*term1123;
        term3(0,22) = RCGinv13*term1113;
        term3(0,23) = RCGinv13*term1112;
        term3(0,24) = RCGinv13*term1132;
        term3(0,25) = RCGinv13*term1131;
        term3(0,26) = RCGinv13*term1121;
        term3(1,0) = RCGinv21*term1211;
        term3(1,1) = RCGinv21*term1222;
        term3(1,2) = RCGinv21*term1233;
        term3(1,3) = RCGinv21*term1223;
        term3(1,4) = RCGinv21*term1213;
        term3(1,5) = RCGinv21*term1212;
        term3(1,6) = RCGinv21*term1232;
        term3(1,7) = RCGinv21*term1231;
        term3(1,8) = RCGinv21*term1221;
        term3(1,9) = RCGinv22*term1211;
        term3(1,10) = RCGinv22*term1222;
        term3(1,11) = RCGinv22*term1233;
        term3(1,12) = RCGinv22*term1223;
        term3(1,13) = RCGinv22*term1213;
        term3(1,14) = RCGinv22*term1212;
        term3(1,15) = RCGinv22*term1232;
        term3(1,16) = RCGinv22*term1231;
        term3(1,17) = RCGinv22*term1221;
        term3(1,18) = RCGinv23*term1211;
        term3(1,19) = RCGinv23*term1222;
        term3(1,20) = RCGinv23*term1233;
        term3(1,21) = RCGinv23*term1223;
        term3(1,22) = RCGinv23*term1213;
        term3(1,23) = RCGinv23*term1212;
        term3(1,24) = RCGinv23*term1232;
        term3(1,25) = RCGinv23*term1231;
        term3(1,26) = RCGinv23*term1221;
        term3(2,0) = RCGinv31*term1311;
        term3(2,1) = RCGinv31*term1322;
        term3(2,2) = RCGinv31*term1333;
        term3(2,3) = RCGinv31*term1323;
        term3(2,4) = RCGinv31*term1313;
        term3(2,5) = RCGinv31*term1312;
        term3(2,6) = RCGinv31*term1332;
        term3(2,7) = RCGinv31*term1331;
        term3(2,8) = RCGinv31*term1321;
        term3(2,9) = RCGinv32*term1311;
        term3(2,10) = RCGinv32*term1322;
        term3(2,11) = RCGinv32*term1333;
        term3(2,12) = RCGinv32*term1323;
        term3(2,13) = RCGinv32*term1313;
        term3(2,14) = RCGinv32*term1312;
        term3(2,15) = RCGinv32*term1332;
        term3(2,16) = RCGinv32*term1331;
        term3(2,17) = RCGinv32*term1321;
        term3(2,18) = RCGinv33*term1311;
        term3(2,19) = RCGinv33*term1322;
        term3(2,20) = RCGinv33*term1333;
        term3(2,21) = RCGinv33*term1323;
        term3(2,22) = RCGinv33*term1313;
        term3(2,23) = RCGinv33*term1312;
        term3(2,24) = RCGinv33*term1332;
        term3(2,25) = RCGinv33*term1331;
        term3(2,26) = RCGinv33*term1321;
        term3(3,0) = RCGinv31*term1211;
        term3(3,1) = RCGinv31*term1222;
        term3(3,2) = RCGinv31*term1233;
        term3(3,3) = RCGinv31*term1223;
        term3(3,4) = RCGinv31*term1213;
        term3(3,5) = RCGinv31*term1212;
        term3(3,6) = RCGinv31*term1232;
        term3(3,7) = RCGinv31*term1231;
        term3(3,8) = RCGinv31*term1221;
        term3(3,9) = RCGinv32*term1211;
        term3(3,10) = RCGinv32*term1222;
        term3(3,11) = RCGinv32*term1233;
        term3(3,12) = RCGinv32*term1223;
        term3(3,13) = RCGinv32*term1213;
        term3(3,14) = RCGinv32*term1212;
        term3(3,15) = RCGinv32*term1232;
        term3(3,16) = RCGinv32*term1231;
        term3(3,17) = RCGinv32*term1221;
        term3(3,18) = RCGinv33*term1211;
        term3(3,19) = RCGinv33*term1222;
        term3(3,20) = RCGinv33*term1233;
        term3(3,21) = RCGinv33*term1223;
        term3(3,22) = RCGinv33*term1213;
        term3(3,23) = RCGinv33*term1212;
        term3(3,24) = RCGinv33*term1232;
        term3(3,25) = RCGinv33*term1231;
        term3(3,26) = RCGinv33*term1221;
        term3(4,0) = RCGinv31*term1111;
        term3(4,1) = RCGinv31*term1122;
        term3(4,2) = RCGinv31*term1133;
        term3(4,3) = RCGinv31*term1123;
        term3(4,4) = RCGinv31*term1113;
        term3(4,5) = RCGinv31*term1112;
        term3(4,6) = RCGinv31*term1132;
        term3(4,7) = RCGinv31*term1131;
        term3(4,8) = RCGinv31*term1121;
        term3(4,9) = RCGinv32*term1111;
        term3(4,10) = RCGinv32*term1122;
        term3(4,11) = RCGinv32*term1133;
        term3(4,12) = RCGinv32*term1123;
        term3(4,13) = RCGinv32*term1113;
        term3(4,14) = RCGinv32*term1112;
        term3(4,15) = RCGinv32*term1132;
        term3(4,16) = RCGinv32*term1131;
        term3(4,17) = RCGinv32*term1121;
        term3(4,18) = RCGinv33*term1111;
        term3(4,19) = RCGinv33*term1122;
        term3(4,20) = RCGinv33*term1133;
        term3(4,21) = RCGinv33*term1123;
        term3(4,22) = RCGinv33*term1113;
        term3(4,23) = RCGinv33*term1112;
        term3(4,24) = RCGinv33*term1132;
        term3(4,25) = RCGinv33*term1131;
        term3(4,26) = RCGinv33*term1121;
        term3(5,0) = RCGinv21*term1111;
        term3(5,1) = RCGinv21*term1122;
        term3(5,2) = RCGinv21*term1133;
        term3(5,3) = RCGinv21*term1123;
        term3(5,4) = RCGinv21*term1113;
        term3(5,5) = RCGinv21*term1112;
        term3(5,6) = RCGinv21*term1132;
        term3(5,7) = RCGinv21*term1131;
        term3(5,8) = RCGinv21*term1121;
        term3(5,9) = RCGinv22*term1111;
        term3(5,10) = RCGinv22*term1122;
        term3(5,11) = RCGinv22*term1133;
        term3(5,12) = RCGinv22*term1123;
        term3(5,13) = RCGinv22*term1113;
        term3(5,14) = RCGinv22*term1112;
        term3(5,15) = RCGinv22*term1132;
        term3(5,16) = RCGinv22*term1131;
        term3(5,17) = RCGinv22*term1121;
        term3(5,18) = RCGinv23*term1111;
        term3(5,19) = RCGinv23*term1122;
        term3(5,20) = RCGinv23*term1133;
        term3(5,21) = RCGinv23*term1123;
        term3(5,22) = RCGinv23*term1113;
        term3(5,23) = RCGinv23*term1112;
        term3(5,24) = RCGinv23*term1132;
        term3(5,25) = RCGinv23*term1131;
        term3(5,26) = RCGinv23*term1121;
        term3(6,0) = RCGinv21*term1311;
        term3(6,1) = RCGinv21*term1322;
        term3(6,2) = RCGinv21*term1333;
        term3(6,3) = RCGinv21*term1323;
        term3(6,4) = RCGinv21*term1313;
        term3(6,5) = RCGinv21*term1312;
        term3(6,6) = RCGinv21*term1332;
        term3(6,7) = RCGinv21*term1331;
        term3(6,8) = RCGinv21*term1321;
        term3(6,9) = RCGinv22*term1311;
        term3(6,10) = RCGinv22*term1322;
        term3(6,11) = RCGinv22*term1333;
        term3(6,12) = RCGinv22*term1323;
        term3(6,13) = RCGinv22*term1313;
        term3(6,14) = RCGinv22*term1312;
        term3(6,15) = RCGinv22*term1332;
        term3(6,16) = RCGinv22*term1331;
        term3(6,17) = RCGinv22*term1321;
        term3(6,18) = RCGinv23*term1311;
        term3(6,19) = RCGinv23*term1322;
        term3(6,20) = RCGinv23*term1333;
        term3(6,21) = RCGinv23*term1323;
        term3(6,22) = RCGinv23*term1313;
        term3(6,23) = RCGinv23*term1312;
        term3(6,24) = RCGinv23*term1332;
        term3(6,25) = RCGinv23*term1331;
        term3(6,26) = RCGinv23*term1321;
        term3(7,0) = RCGinv11*term1311;
        term3(7,1) = RCGinv11*term1322;
        term3(7,2) = RCGinv11*term1333;
        term3(7,3) = RCGinv11*term1323;
        term3(7,4) = RCGinv11*term1313;
        term3(7,5) = RCGinv11*term1312;
        term3(7,6) = RCGinv11*term1332;
        term3(7,7) = RCGinv11*term1331;
        term3(7,8) = RCGinv11*term1321;
        term3(7,9) = RCGinv12*term1311;
        term3(7,10) = RCGinv12*term1322;
        term3(7,11) = RCGinv12*term1333;
        term3(7,12) = RCGinv12*term1323;
        term3(7,13) = RCGinv12*term1313;
        term3(7,14) = RCGinv12*term1312;
        term3(7,15) = RCGinv12*term1332;
        term3(7,16) = RCGinv12*term1331;
        term3(7,17) = RCGinv12*term1321;
        term3(7,18) = RCGinv13*term1311;
        term3(7,19) = RCGinv13*term1322;
        term3(7,20) = RCGinv13*term1333;
        term3(7,21) = RCGinv13*term1323;
        term3(7,22) = RCGinv13*term1313;
        term3(7,23) = RCGinv13*term1312;
        term3(7,24) = RCGinv13*term1332;
        term3(7,25) = RCGinv13*term1331;
        term3(7,26) = RCGinv13*term1321;
        term3(8,0) = RCGinv11*term1211;
        term3(8,1) = RCGinv11*term1222;
        term3(8,2) = RCGinv11*term1233;
        term3(8,3) = RCGinv11*term1223;
        term3(8,4) = RCGinv11*term1213;
        term3(8,5) = RCGinv11*term1212;
        term3(8,6) = RCGinv11*term1232;
        term3(8,7) = RCGinv11*term1231;
        term3(8,8) = RCGinv11*term1221;
        term3(8,9) = RCGinv12*term1211;
        term3(8,10) = RCGinv12*term1222;
        term3(8,11) = RCGinv12*term1233;
        term3(8,12) = RCGinv12*term1223;
        term3(8,13) = RCGinv12*term1213;
        term3(8,14) = RCGinv12*term1212;
        term3(8,15) = RCGinv12*term1232;
        term3(8,16) = RCGinv12*term1231;
        term3(8,17) = RCGinv12*term1221;
        term3(8,18) = RCGinv13*term1211;
        term3(8,19) = RCGinv13*term1222;
        term3(8,20) = RCGinv13*term1233;
        term3(8,21) = RCGinv13*term1223;
        term3(8,22) = RCGinv13*term1213;
        term3(8,23) = RCGinv13*term1212;
        term3(8,24) = RCGinv13*term1232;
        term3(8,25) = RCGinv13*term1231;
        term3(8,26) = RCGinv13*term1221;
        return;
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_9x27 term2;
        Matrix_27x27 _C = C; //Copy the sparse matrix to a dense matrix.
        compute_dPK2dGamma_term2(term1,_C,term2);
        
        //Compute term3
        Matrix_9x27 term3;
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,term3);
        
        dPK2dGamma = term2 + term3;
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C, Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_9x27 term2;
        compute_dPK2dGamma_term2(term1, C,term2);
        
        //Compute term3
        Matrix_9x27 term3;
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,term3);
        
        dPK2dGamma = term2 + term3;
    }

    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            SpMat &C, Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        Matrix_27x27 _C = C; //Copy the sparse matrix to a dense matrix.
        compute_dPK2dGamma_term2(term1,_C,terms[0]);
        
        //Compute term3
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,terms[1]);
        
        dPK2dGamma = terms[0] + terms[1];
    }
    
    void compute_dPK2dGamma(const Matrix_3x3 &RCGinv, const Matrix_3x9 &Gamma, const Vector_27 &Gamma_voigt,
                            Matrix_27x27 &C, Matrix_9x27 (&terms)[2], Matrix_9x27 &dPK2dGamma){
        /*!============================
        |    compute_dPK2dGamma    |
        ============================
        
        Compute the derivative of the second piola kirchoff stress 
        w.r.t. the deformation measure Gamma.
        
        */
        
        //Compute term1
        Vector_27 term1;
        deformation_measures::voigt_3x9_tensor(RCGinv*Gamma,term1);
        
        //Compute term2
        compute_dPK2dGamma_term2(term1,C,terms[0]);
        
        //Compute term3
        compute_dPK2dGamma_term3(C*Gamma_voigt,RCGinv,terms[1]);
        
        dPK2dGamma = terms[0] + terms[1];
    }

    void compute_dSIGMAdRCG(Matrix_9x9 (&terms)[4], Matrix_9x9 &dSIGMAdRCG){
        /*!=========================
        |    compute_dSIGMAdRCG    |
        ============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the 
        computation of dPK2dRCG.
        
        This is the preferred method since the computation has 
        already been done for the PK2 derivative.
        
        */
        
        dSIGMAdRCG = terms[0];
        
        Matrix_9x9 temp;

        temp = terms[1] + terms[2] + terms[3]; //TODO: Just use two terms. No need for four
        dSIGMAdRCG        += temp;
        dSIGMAdRCG.row(0) += temp.row(0);
        dSIGMAdRCG.row(1) += temp.row(1);
        dSIGMAdRCG.row(2) += temp.row(2);
        
        dSIGMAdRCG.row(3) += temp.row(6);
        dSIGMAdRCG.row(4) += temp.row(7);
        dSIGMAdRCG.row(5) += temp.row(8);
        dSIGMAdRCG.row(6) += temp.row(3);
        dSIGMAdRCG.row(7) += temp.row(4);
        dSIGMAdRCG.row(8) += temp.row(5);
 
//        for (int i=1; i<4; i++){
//            temp = terms[i];
//            
//            temp.row(0) *= 2;
//            temp.row(1) *= 2;
//            temp.row(2) *= 2;
//            temp.row(3) += terms[i].row(6);
//            temp.row(4) += terms[i].row(7);
//            temp.row(5) += terms[i].row(8);
//            temp.row(6) += terms[i].row(3);
//            temp.row(7) += terms[i].row(4);
//            temp.row(8) += terms[i].row(5);
//            
//            dSIGMAdRCG += temp;
//        }
        return;
    }
    
    void compute_dSIGMAdPsi(Matrix_9x9 (&terms)[3], Matrix_9x9 &dSIGMAdPsi){
        /*!============================
        |    compute_dSIGMAdPsi    |
        ============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the 
        computation of dPK2dPsi.
        
        This is the preferred method since the computation has already 
        been done for the PK2 derivative.
        
        */
        
        dSIGMAdPsi = terms[0];
        
        Matrix_9x9 temp;
        
        for (int i=1; i<3; i++){
            temp = terms[i];
            
            temp.row(0) *= 2;
            temp.row(1) *= 2;
            temp.row(2) *= 2;
            
            temp.row(3) += terms[i].row(6);
            temp.row(4) += terms[i].row(7);
            temp.row(5) += terms[i].row(8);
            temp.row(6) += terms[i].row(3);
            temp.row(7) += terms[i].row(4);
            temp.row(8) += terms[i].row(5);
            
            dSIGMAdPsi += temp;
        }
        return;
    }
    
    void compute_dSIGMAdGamma(Matrix_9x27 (&terms)[2], Matrix_9x27 &dSIGMAdGamma){
        /*!===========================
        |    compute_dSIGMAdGamma    |
        ==============================
        
        Compute the derivative of the symmetric stress in the 
        reference configuration using the terms from the computation 
        of dPK2dGamma.
        
        This is the preferred method since the computation has already
        been done for the PK2 derivative.
        
        */
        
        dSIGMAdGamma = Matrix_9x27::Zero();
        
        Matrix_9x27 temp;
        
        for (int i=0; i<2; i++){
            temp = terms[i];
            
            temp.row(0) *= 2.;
            temp.row(1) *= 2.;
            temp.row(2) *= 2.;
            temp.row(3) += terms[i].row(6);
            temp.row(4) += terms[i].row(7);
            temp.row(5) += terms[i].row(8);
            temp.row(6) += terms[i].row(3);
            temp.row(7) += terms[i].row(4);
            temp.row(8) += terms[i].row(5);
            
            dSIGMAdGamma += temp;
        }
        return;
    }
    
    void compute_dMdGamma(const Matrix_27x27 &C, Matrix_27x27 &dMdGamma){
        /*!==========================
        |    compute_dMdGamma    |
        ==========================
        
        Compute the derivative of the higher order stress tensor w.r.t. Gamma.
        
        */
        
        //Extract C
        double C111111 = C(0,0);
        double C111122 = C(0,1);
        double C111133 = C(0,2);
        double C111123 = C(0,3);
        double C111113 = C(0,4);
        double C111112 = C(0,5);
        double C111132 = C(0,6);
        double C111131 = C(0,7);
        double C111121 = C(0,8);
        double C111211 = C(0,9);
        double C111222 = C(0,10);
        double C111233 = C(0,11);
        double C111223 = C(0,12);
        double C111213 = C(0,13);
        double C111212 = C(0,14);
        double C111232 = C(0,15);
        double C111231 = C(0,16);
        double C111221 = C(0,17);
        double C111311 = C(0,18);
        double C111322 = C(0,19);
        double C111333 = C(0,20);
        double C111323 = C(0,21);
        double C111313 = C(0,22);
        double C111312 = C(0,23);
        double C111332 = C(0,24);
        double C111331 = C(0,25);
        double C111321 = C(0,26);
        double C122111 = C(1,0);
        double C122122 = C(1,1);
        double C122133 = C(1,2);
        double C122123 = C(1,3);
        double C122113 = C(1,4);
        double C122112 = C(1,5);
        double C122132 = C(1,6);
        double C122131 = C(1,7);
        double C122121 = C(1,8);
        double C122211 = C(1,9);
        double C122222 = C(1,10);
        double C122233 = C(1,11);
        double C122223 = C(1,12);
        double C122213 = C(1,13);
        double C122212 = C(1,14);
        double C122232 = C(1,15);
        double C122231 = C(1,16);
        double C122221 = C(1,17);
        double C122311 = C(1,18);
        double C122322 = C(1,19);
        double C122333 = C(1,20);
        double C122323 = C(1,21);
        double C122313 = C(1,22);
        double C122312 = C(1,23);
        double C122332 = C(1,24);
        double C122331 = C(1,25);
        double C122321 = C(1,26);
        double C133111 = C(2,0);
        double C133122 = C(2,1);
        double C133133 = C(2,2);
        double C133123 = C(2,3);
        double C133113 = C(2,4);
        double C133112 = C(2,5);
        double C133132 = C(2,6);
        double C133131 = C(2,7);
        double C133121 = C(2,8);
        double C133211 = C(2,9);
        double C133222 = C(2,10);
        double C133233 = C(2,11);
        double C133223 = C(2,12);
        double C133213 = C(2,13);
        double C133212 = C(2,14);
        double C133232 = C(2,15);
        double C133231 = C(2,16);
        double C133221 = C(2,17);
        double C133311 = C(2,18);
        double C133322 = C(2,19);
        double C133333 = C(2,20);
        double C133323 = C(2,21);
        double C133313 = C(2,22);
        double C133312 = C(2,23);
        double C133332 = C(2,24);
        double C133331 = C(2,25);
        double C133321 = C(2,26);
        double C123111 = C(3,0);
        double C123122 = C(3,1);
        double C123133 = C(3,2);
        double C123123 = C(3,3);
        double C123113 = C(3,4);
        double C123112 = C(3,5);
        double C123132 = C(3,6);
        double C123131 = C(3,7);
        double C123121 = C(3,8);
        double C123211 = C(3,9);
        double C123222 = C(3,10);
        double C123233 = C(3,11);
        double C123223 = C(3,12);
        double C123213 = C(3,13);
        double C123212 = C(3,14);
        double C123232 = C(3,15);
        double C123231 = C(3,16);
        double C123221 = C(3,17);
        double C123311 = C(3,18);
        double C123322 = C(3,19);
        double C123333 = C(3,20);
        double C123323 = C(3,21);
        double C123313 = C(3,22);
        double C123312 = C(3,23);
        double C123332 = C(3,24);
        double C123331 = C(3,25);
        double C123321 = C(3,26);
        double C113111 = C(4,0);
        double C113122 = C(4,1);
        double C113133 = C(4,2);
        double C113123 = C(4,3);
        double C113113 = C(4,4);
        double C113112 = C(4,5);
        double C113132 = C(4,6);
        double C113131 = C(4,7);
        double C113121 = C(4,8);
        double C113211 = C(4,9);
        double C113222 = C(4,10);
        double C113233 = C(4,11);
        double C113223 = C(4,12);
        double C113213 = C(4,13);
        double C113212 = C(4,14);
        double C113232 = C(4,15);
        double C113231 = C(4,16);
        double C113221 = C(4,17);
        double C113311 = C(4,18);
        double C113322 = C(4,19);
        double C113333 = C(4,20);
        double C113323 = C(4,21);
        double C113313 = C(4,22);
        double C113312 = C(4,23);
        double C113332 = C(4,24);
        double C113331 = C(4,25);
        double C113321 = C(4,26);
        double C112111 = C(5,0);
        double C112122 = C(5,1);
        double C112133 = C(5,2);
        double C112123 = C(5,3);
        double C112113 = C(5,4);
        double C112112 = C(5,5);
        double C112132 = C(5,6);
        double C112131 = C(5,7);
        double C112121 = C(5,8);
        double C112211 = C(5,9);
        double C112222 = C(5,10);
        double C112233 = C(5,11);
        double C112223 = C(5,12);
        double C112213 = C(5,13);
        double C112212 = C(5,14);
        double C112232 = C(5,15);
        double C112231 = C(5,16);
        double C112221 = C(5,17);
        double C112311 = C(5,18);
        double C112322 = C(5,19);
        double C112333 = C(5,20);
        double C112323 = C(5,21);
        double C112313 = C(5,22);
        double C112312 = C(5,23);
        double C112332 = C(5,24);
        double C112331 = C(5,25);
        double C112321 = C(5,26);
        double C132111 = C(6,0);
        double C132122 = C(6,1);
        double C132133 = C(6,2);
        double C132123 = C(6,3);
        double C132113 = C(6,4);
        double C132112 = C(6,5);
        double C132132 = C(6,6);
        double C132131 = C(6,7);
        double C132121 = C(6,8);
        double C132211 = C(6,9);
        double C132222 = C(6,10);
        double C132233 = C(6,11);
        double C132223 = C(6,12);
        double C132213 = C(6,13);
        double C132212 = C(6,14);
        double C132232 = C(6,15);
        double C132231 = C(6,16);
        double C132221 = C(6,17);
        double C132311 = C(6,18);
        double C132322 = C(6,19);
        double C132333 = C(6,20);
        double C132323 = C(6,21);
        double C132313 = C(6,22);
        double C132312 = C(6,23);
        double C132332 = C(6,24);
        double C132331 = C(6,25);
        double C132321 = C(6,26);
        double C131111 = C(7,0);
        double C131122 = C(7,1);
        double C131133 = C(7,2);
        double C131123 = C(7,3);
        double C131113 = C(7,4);
        double C131112 = C(7,5);
        double C131132 = C(7,6);
        double C131131 = C(7,7);
        double C131121 = C(7,8);
        double C131211 = C(7,9);
        double C131222 = C(7,10);
        double C131233 = C(7,11);
        double C131223 = C(7,12);
        double C131213 = C(7,13);
        double C131212 = C(7,14);
        double C131232 = C(7,15);
        double C131231 = C(7,16);
        double C131221 = C(7,17);
        double C131311 = C(7,18);
        double C131322 = C(7,19);
        double C131333 = C(7,20);
        double C131323 = C(7,21);
        double C131313 = C(7,22);
        double C131312 = C(7,23);
        double C131332 = C(7,24);
        double C131331 = C(7,25);
        double C131321 = C(7,26);
        double C121111 = C(8,0);
        double C121122 = C(8,1);
        double C121133 = C(8,2);
        double C121123 = C(8,3);
        double C121113 = C(8,4);
        double C121112 = C(8,5);
        double C121132 = C(8,6);
        double C121131 = C(8,7);
        double C121121 = C(8,8);
        double C121211 = C(8,9);
        double C121222 = C(8,10);
        double C121233 = C(8,11);
        double C121223 = C(8,12);
        double C121213 = C(8,13);
        double C121212 = C(8,14);
        double C121232 = C(8,15);
        double C121231 = C(8,16);
        double C121221 = C(8,17);
        double C121311 = C(8,18);
        double C121322 = C(8,19);
        double C121333 = C(8,20);
        double C121323 = C(8,21);
        double C121313 = C(8,22);
        double C121312 = C(8,23);
        double C121332 = C(8,24);
        double C121331 = C(8,25);
        double C121321 = C(8,26);
        double C211111 = C(9,0);
        double C211122 = C(9,1);
        double C211133 = C(9,2);
        double C211123 = C(9,3);
        double C211113 = C(9,4);
        double C211112 = C(9,5);
        double C211132 = C(9,6);
        double C211131 = C(9,7);
        double C211121 = C(9,8);
        double C211211 = C(9,9);
        double C211222 = C(9,10);
        double C211233 = C(9,11);
        double C211223 = C(9,12);
        double C211213 = C(9,13);
        double C211212 = C(9,14);
        double C211232 = C(9,15);
        double C211231 = C(9,16);
        double C211221 = C(9,17);
        double C211311 = C(9,18);
        double C211322 = C(9,19);
        double C211333 = C(9,20);
        double C211323 = C(9,21);
        double C211313 = C(9,22);
        double C211312 = C(9,23);
        double C211332 = C(9,24);
        double C211331 = C(9,25);
        double C211321 = C(9,26);
        double C222111 = C(10,0);
        double C222122 = C(10,1);
        double C222133 = C(10,2);
        double C222123 = C(10,3);
        double C222113 = C(10,4);
        double C222112 = C(10,5);
        double C222132 = C(10,6);
        double C222131 = C(10,7);
        double C222121 = C(10,8);
        double C222211 = C(10,9);
        double C222222 = C(10,10);
        double C222233 = C(10,11);
        double C222223 = C(10,12);
        double C222213 = C(10,13);
        double C222212 = C(10,14);
        double C222232 = C(10,15);
        double C222231 = C(10,16);
        double C222221 = C(10,17);
        double C222311 = C(10,18);
        double C222322 = C(10,19);
        double C222333 = C(10,20);
        double C222323 = C(10,21);
        double C222313 = C(10,22);
        double C222312 = C(10,23);
        double C222332 = C(10,24);
        double C222331 = C(10,25);
        double C222321 = C(10,26);
        double C233111 = C(11,0);
        double C233122 = C(11,1);
        double C233133 = C(11,2);
        double C233123 = C(11,3);
        double C233113 = C(11,4);
        double C233112 = C(11,5);
        double C233132 = C(11,6);
        double C233131 = C(11,7);
        double C233121 = C(11,8);
        double C233211 = C(11,9);
        double C233222 = C(11,10);
        double C233233 = C(11,11);
        double C233223 = C(11,12);
        double C233213 = C(11,13);
        double C233212 = C(11,14);
        double C233232 = C(11,15);
        double C233231 = C(11,16);
        double C233221 = C(11,17);
        double C233311 = C(11,18);
        double C233322 = C(11,19);
        double C233333 = C(11,20);
        double C233323 = C(11,21);
        double C233313 = C(11,22);
        double C233312 = C(11,23);
        double C233332 = C(11,24);
        double C233331 = C(11,25);
        double C233321 = C(11,26);
        double C223111 = C(12,0);
        double C223122 = C(12,1);
        double C223133 = C(12,2);
        double C223123 = C(12,3);
        double C223113 = C(12,4);
        double C223112 = C(12,5);
        double C223132 = C(12,6);
        double C223131 = C(12,7);
        double C223121 = C(12,8);
        double C223211 = C(12,9);
        double C223222 = C(12,10);
        double C223233 = C(12,11);
        double C223223 = C(12,12);
        double C223213 = C(12,13);
        double C223212 = C(12,14);
        double C223232 = C(12,15);
        double C223231 = C(12,16);
        double C223221 = C(12,17);
        double C223311 = C(12,18);
        double C223322 = C(12,19);
        double C223333 = C(12,20);
        double C223323 = C(12,21);
        double C223313 = C(12,22);
        double C223312 = C(12,23);
        double C223332 = C(12,24);
        double C223331 = C(12,25);
        double C223321 = C(12,26);
        double C213111 = C(13,0);
        double C213122 = C(13,1);
        double C213133 = C(13,2);
        double C213123 = C(13,3);
        double C213113 = C(13,4);
        double C213112 = C(13,5);
        double C213132 = C(13,6);
        double C213131 = C(13,7);
        double C213121 = C(13,8);
        double C213211 = C(13,9);
        double C213222 = C(13,10);
        double C213233 = C(13,11);
        double C213223 = C(13,12);
        double C213213 = C(13,13);
        double C213212 = C(13,14);
        double C213232 = C(13,15);
        double C213231 = C(13,16);
        double C213221 = C(13,17);
        double C213311 = C(13,18);
        double C213322 = C(13,19);
        double C213333 = C(13,20);
        double C213323 = C(13,21);
        double C213313 = C(13,22);
        double C213312 = C(13,23);
        double C213332 = C(13,24);
        double C213331 = C(13,25);
        double C213321 = C(13,26);
        double C212111 = C(14,0);
        double C212122 = C(14,1);
        double C212133 = C(14,2);
        double C212123 = C(14,3);
        double C212113 = C(14,4);
        double C212112 = C(14,5);
        double C212132 = C(14,6);
        double C212131 = C(14,7);
        double C212121 = C(14,8);
        double C212211 = C(14,9);
        double C212222 = C(14,10);
        double C212233 = C(14,11);
        double C212223 = C(14,12);
        double C212213 = C(14,13);
        double C212212 = C(14,14);
        double C212232 = C(14,15);
        double C212231 = C(14,16);
        double C212221 = C(14,17);
        double C212311 = C(14,18);
        double C212322 = C(14,19);
        double C212333 = C(14,20);
        double C212323 = C(14,21);
        double C212313 = C(14,22);
        double C212312 = C(14,23);
        double C212332 = C(14,24);
        double C212331 = C(14,25);
        double C212321 = C(14,26);
        double C232111 = C(15,0);
        double C232122 = C(15,1);
        double C232133 = C(15,2);
        double C232123 = C(15,3);
        double C232113 = C(15,4);
        double C232112 = C(15,5);
        double C232132 = C(15,6);
        double C232131 = C(15,7);
        double C232121 = C(15,8);
        double C232211 = C(15,9);
        double C232222 = C(15,10);
        double C232233 = C(15,11);
        double C232223 = C(15,12);
        double C232213 = C(15,13);
        double C232212 = C(15,14);
        double C232232 = C(15,15);
        double C232231 = C(15,16);
        double C232221 = C(15,17);
        double C232311 = C(15,18);
        double C232322 = C(15,19);
        double C232333 = C(15,20);
        double C232323 = C(15,21);
        double C232313 = C(15,22);
        double C232312 = C(15,23);
        double C232332 = C(15,24);
        double C232331 = C(15,25);
        double C232321 = C(15,26);
        double C231111 = C(16,0);
        double C231122 = C(16,1);
        double C231133 = C(16,2);
        double C231123 = C(16,3);
        double C231113 = C(16,4);
        double C231112 = C(16,5);
        double C231132 = C(16,6);
        double C231131 = C(16,7);
        double C231121 = C(16,8);
        double C231211 = C(16,9);
        double C231222 = C(16,10);
        double C231233 = C(16,11);
        double C231223 = C(16,12);
        double C231213 = C(16,13);
        double C231212 = C(16,14);
        double C231232 = C(16,15);
        double C231231 = C(16,16);
        double C231221 = C(16,17);
        double C231311 = C(16,18);
        double C231322 = C(16,19);
        double C231333 = C(16,20);
        double C231323 = C(16,21);
        double C231313 = C(16,22);
        double C231312 = C(16,23);
        double C231332 = C(16,24);
        double C231331 = C(16,25);
        double C231321 = C(16,26);
        double C221111 = C(17,0);
        double C221122 = C(17,1);
        double C221133 = C(17,2);
        double C221123 = C(17,3);
        double C221113 = C(17,4);
        double C221112 = C(17,5);
        double C221132 = C(17,6);
        double C221131 = C(17,7);
        double C221121 = C(17,8);
        double C221211 = C(17,9);
        double C221222 = C(17,10);
        double C221233 = C(17,11);
        double C221223 = C(17,12);
        double C221213 = C(17,13);
        double C221212 = C(17,14);
        double C221232 = C(17,15);
        double C221231 = C(17,16);
        double C221221 = C(17,17);
        double C221311 = C(17,18);
        double C221322 = C(17,19);
        double C221333 = C(17,20);
        double C221323 = C(17,21);
        double C221313 = C(17,22);
        double C221312 = C(17,23);
        double C221332 = C(17,24);
        double C221331 = C(17,25);
        double C221321 = C(17,26);
        double C311111 = C(18,0);
        double C311122 = C(18,1);
        double C311133 = C(18,2);
        double C311123 = C(18,3);
        double C311113 = C(18,4);
        double C311112 = C(18,5);
        double C311132 = C(18,6);
        double C311131 = C(18,7);
        double C311121 = C(18,8);
        double C311211 = C(18,9);
        double C311222 = C(18,10);
        double C311233 = C(18,11);
        double C311223 = C(18,12);
        double C311213 = C(18,13);
        double C311212 = C(18,14);
        double C311232 = C(18,15);
        double C311231 = C(18,16);
        double C311221 = C(18,17);
        double C311311 = C(18,18);
        double C311322 = C(18,19);
        double C311333 = C(18,20);
        double C311323 = C(18,21);
        double C311313 = C(18,22);
        double C311312 = C(18,23);
        double C311332 = C(18,24);
        double C311331 = C(18,25);
        double C311321 = C(18,26);
        double C322111 = C(19,0);
        double C322122 = C(19,1);
        double C322133 = C(19,2);
        double C322123 = C(19,3);
        double C322113 = C(19,4);
        double C322112 = C(19,5);
        double C322132 = C(19,6);
        double C322131 = C(19,7);
        double C322121 = C(19,8);
        double C322211 = C(19,9);
        double C322222 = C(19,10);
        double C322233 = C(19,11);
        double C322223 = C(19,12);
        double C322213 = C(19,13);
        double C322212 = C(19,14);
        double C322232 = C(19,15);
        double C322231 = C(19,16);
        double C322221 = C(19,17);
        double C322311 = C(19,18);
        double C322322 = C(19,19);
        double C322333 = C(19,20);
        double C322323 = C(19,21);
        double C322313 = C(19,22);
        double C322312 = C(19,23);
        double C322332 = C(19,24);
        double C322331 = C(19,25);
        double C322321 = C(19,26);
        double C333111 = C(20,0);
        double C333122 = C(20,1);
        double C333133 = C(20,2);
        double C333123 = C(20,3);
        double C333113 = C(20,4);
        double C333112 = C(20,5);
        double C333132 = C(20,6);
        double C333131 = C(20,7);
        double C333121 = C(20,8);
        double C333211 = C(20,9);
        double C333222 = C(20,10);
        double C333233 = C(20,11);
        double C333223 = C(20,12);
        double C333213 = C(20,13);
        double C333212 = C(20,14);
        double C333232 = C(20,15);
        double C333231 = C(20,16);
        double C333221 = C(20,17);
        double C333311 = C(20,18);
        double C333322 = C(20,19);
        double C333333 = C(20,20);
        double C333323 = C(20,21);
        double C333313 = C(20,22);
        double C333312 = C(20,23);
        double C333332 = C(20,24);
        double C333331 = C(20,25);
        double C333321 = C(20,26);
        double C323111 = C(21,0);
        double C323122 = C(21,1);
        double C323133 = C(21,2);
        double C323123 = C(21,3);
        double C323113 = C(21,4);
        double C323112 = C(21,5);
        double C323132 = C(21,6);
        double C323131 = C(21,7);
        double C323121 = C(21,8);
        double C323211 = C(21,9);
        double C323222 = C(21,10);
        double C323233 = C(21,11);
        double C323223 = C(21,12);
        double C323213 = C(21,13);
        double C323212 = C(21,14);
        double C323232 = C(21,15);
        double C323231 = C(21,16);
        double C323221 = C(21,17);
        double C323311 = C(21,18);
        double C323322 = C(21,19);
        double C323333 = C(21,20);
        double C323323 = C(21,21);
        double C323313 = C(21,22);
        double C323312 = C(21,23);
        double C323332 = C(21,24);
        double C323331 = C(21,25);
        double C323321 = C(21,26);
        double C313111 = C(22,0);
        double C313122 = C(22,1);
        double C313133 = C(22,2);
        double C313123 = C(22,3);
        double C313113 = C(22,4);
        double C313112 = C(22,5);
        double C313132 = C(22,6);
        double C313131 = C(22,7);
        double C313121 = C(22,8);
        double C313211 = C(22,9);
        double C313222 = C(22,10);
        double C313233 = C(22,11);
        double C313223 = C(22,12);
        double C313213 = C(22,13);
        double C313212 = C(22,14);
        double C313232 = C(22,15);
        double C313231 = C(22,16);
        double C313221 = C(22,17);
        double C313311 = C(22,18);
        double C313322 = C(22,19);
        double C313333 = C(22,20);
        double C313323 = C(22,21);
        double C313313 = C(22,22);
        double C313312 = C(22,23);
        double C313332 = C(22,24);
        double C313331 = C(22,25);
        double C313321 = C(22,26);
        double C312111 = C(23,0);
        double C312122 = C(23,1);
        double C312133 = C(23,2);
        double C312123 = C(23,3);
        double C312113 = C(23,4);
        double C312112 = C(23,5);
        double C312132 = C(23,6);
        double C312131 = C(23,7);
        double C312121 = C(23,8);
        double C312211 = C(23,9);
        double C312222 = C(23,10);
        double C312233 = C(23,11);
        double C312223 = C(23,12);
        double C312213 = C(23,13);
        double C312212 = C(23,14);
        double C312232 = C(23,15);
        double C312231 = C(23,16);
        double C312221 = C(23,17);
        double C312311 = C(23,18);
        double C312322 = C(23,19);
        double C312333 = C(23,20);
        double C312323 = C(23,21);
        double C312313 = C(23,22);
        double C312312 = C(23,23);
        double C312332 = C(23,24);
        double C312331 = C(23,25);
        double C312321 = C(23,26);
        double C332111 = C(24,0);
        double C332122 = C(24,1);
        double C332133 = C(24,2);
        double C332123 = C(24,3);
        double C332113 = C(24,4);
        double C332112 = C(24,5);
        double C332132 = C(24,6);
        double C332131 = C(24,7);
        double C332121 = C(24,8);
        double C332211 = C(24,9);
        double C332222 = C(24,10);
        double C332233 = C(24,11);
        double C332223 = C(24,12);
        double C332213 = C(24,13);
        double C332212 = C(24,14);
        double C332232 = C(24,15);
        double C332231 = C(24,16);
        double C332221 = C(24,17);
        double C332311 = C(24,18);
        double C332322 = C(24,19);
        double C332333 = C(24,20);
        double C332323 = C(24,21);
        double C332313 = C(24,22);
        double C332312 = C(24,23);
        double C332332 = C(24,24);
        double C332331 = C(24,25);
        double C332321 = C(24,26);
        double C331111 = C(25,0);
        double C331122 = C(25,1);
        double C331133 = C(25,2);
        double C331123 = C(25,3);
        double C331113 = C(25,4);
        double C331112 = C(25,5);
        double C331132 = C(25,6);
        double C331131 = C(25,7);
        double C331121 = C(25,8);
        double C331211 = C(25,9);
        double C331222 = C(25,10);
        double C331233 = C(25,11);
        double C331223 = C(25,12);
        double C331213 = C(25,13);
        double C331212 = C(25,14);
        double C331232 = C(25,15);
        double C331231 = C(25,16);
        double C331221 = C(25,17);
        double C331311 = C(25,18);
        double C331322 = C(25,19);
        double C331333 = C(25,20);
        double C331323 = C(25,21);
        double C331313 = C(25,22);
        double C331312 = C(25,23);
        double C331332 = C(25,24);
        double C331331 = C(25,25);
        double C331321 = C(25,26);
        double C321111 = C(26,0);
        double C321122 = C(26,1);
        double C321133 = C(26,2);
        double C321123 = C(26,3);
        double C321113 = C(26,4);
        double C321112 = C(26,5);
        double C321132 = C(26,6);
        double C321131 = C(26,7);
        double C321121 = C(26,8);
        double C321211 = C(26,9);
        double C321222 = C(26,10);
        double C321233 = C(26,11);
        double C321223 = C(26,12);
        double C321213 = C(26,13);
        double C321212 = C(26,14);
        double C321232 = C(26,15);
        double C321231 = C(26,16);
        double C321221 = C(26,17);
        double C321311 = C(26,18);
        double C321322 = C(26,19);
        double C321333 = C(26,20);
        double C321323 = C(26,21);
        double C321313 = C(26,22);
        double C321312 = C(26,23);
        double C321332 = C(26,24);
        double C321331 = C(26,25);
        double C321321 = C(26,26);

        //Assemble dMdGamma
        dMdGamma(0,0) = C111111;
        dMdGamma(0,1) = C111122;
        dMdGamma(0,2) = C111133;
        dMdGamma(0,3) = C111123;
        dMdGamma(0,4) = C111113;
        dMdGamma(0,5) = C111112;
        dMdGamma(0,6) = C111132;
        dMdGamma(0,7) = C111131;
        dMdGamma(0,8) = C111121;
        dMdGamma(0,9) = C111211;
        dMdGamma(0,10) = C111222;
        dMdGamma(0,11) = C111233;
        dMdGamma(0,12) = C111223;
        dMdGamma(0,13) = C111213;
        dMdGamma(0,14) = C111212;
        dMdGamma(0,15) = C111232;
        dMdGamma(0,16) = C111231;
        dMdGamma(0,17) = C111221;
        dMdGamma(0,18) = C111311;
        dMdGamma(0,19) = C111322;
        dMdGamma(0,20) = C111333;
        dMdGamma(0,21) = C111323;
        dMdGamma(0,22) = C111313;
        dMdGamma(0,23) = C111312;
        dMdGamma(0,24) = C111332;
        dMdGamma(0,25) = C111331;
        dMdGamma(0,26) = C111321;
        dMdGamma(1,0) = C221111;
        dMdGamma(1,1) = C221122;
        dMdGamma(1,2) = C221133;
        dMdGamma(1,3) = C221123;
        dMdGamma(1,4) = C221113;
        dMdGamma(1,5) = C221112;
        dMdGamma(1,6) = C221132;
        dMdGamma(1,7) = C221131;
        dMdGamma(1,8) = C221121;
        dMdGamma(1,9) = C221211;
        dMdGamma(1,10) = C221222;
        dMdGamma(1,11) = C221233;
        dMdGamma(1,12) = C221223;
        dMdGamma(1,13) = C221213;
        dMdGamma(1,14) = C221212;
        dMdGamma(1,15) = C221232;
        dMdGamma(1,16) = C221231;
        dMdGamma(1,17) = C221221;
        dMdGamma(1,18) = C221311;
        dMdGamma(1,19) = C221322;
        dMdGamma(1,20) = C221333;
        dMdGamma(1,21) = C221323;
        dMdGamma(1,22) = C221313;
        dMdGamma(1,23) = C221312;
        dMdGamma(1,24) = C221332;
        dMdGamma(1,25) = C221331;
        dMdGamma(1,26) = C221321;
        dMdGamma(2,0) = C331111;
        dMdGamma(2,1) = C331122;
        dMdGamma(2,2) = C331133;
        dMdGamma(2,3) = C331123;
        dMdGamma(2,4) = C331113;
        dMdGamma(2,5) = C331112;
        dMdGamma(2,6) = C331132;
        dMdGamma(2,7) = C331131;
        dMdGamma(2,8) = C331121;
        dMdGamma(2,9) = C331211;
        dMdGamma(2,10) = C331222;
        dMdGamma(2,11) = C331233;
        dMdGamma(2,12) = C331223;
        dMdGamma(2,13) = C331213;
        dMdGamma(2,14) = C331212;
        dMdGamma(2,15) = C331232;
        dMdGamma(2,16) = C331231;
        dMdGamma(2,17) = C331221;
        dMdGamma(2,18) = C331311;
        dMdGamma(2,19) = C331322;
        dMdGamma(2,20) = C331333;
        dMdGamma(2,21) = C331323;
        dMdGamma(2,22) = C331313;
        dMdGamma(2,23) = C331312;
        dMdGamma(2,24) = C331332;
        dMdGamma(2,25) = C331331;
        dMdGamma(2,26) = C331321;
        dMdGamma(3,0) = C231111;
        dMdGamma(3,1) = C231122;
        dMdGamma(3,2) = C231133;
        dMdGamma(3,3) = C231123;
        dMdGamma(3,4) = C231113;
        dMdGamma(3,5) = C231112;
        dMdGamma(3,6) = C231132;
        dMdGamma(3,7) = C231131;
        dMdGamma(3,8) = C231121;
        dMdGamma(3,9) = C231211;
        dMdGamma(3,10) = C231222;
        dMdGamma(3,11) = C231233;
        dMdGamma(3,12) = C231223;
        dMdGamma(3,13) = C231213;
        dMdGamma(3,14) = C231212;
        dMdGamma(3,15) = C231232;
        dMdGamma(3,16) = C231231;
        dMdGamma(3,17) = C231221;
        dMdGamma(3,18) = C231311;
        dMdGamma(3,19) = C231322;
        dMdGamma(3,20) = C231333;
        dMdGamma(3,21) = C231323;
        dMdGamma(3,22) = C231313;
        dMdGamma(3,23) = C231312;
        dMdGamma(3,24) = C231332;
        dMdGamma(3,25) = C231331;
        dMdGamma(3,26) = C231321;
        dMdGamma(4,0) = C131111;
        dMdGamma(4,1) = C131122;
        dMdGamma(4,2) = C131133;
        dMdGamma(4,3) = C131123;
        dMdGamma(4,4) = C131113;
        dMdGamma(4,5) = C131112;
        dMdGamma(4,6) = C131132;
        dMdGamma(4,7) = C131131;
        dMdGamma(4,8) = C131121;
        dMdGamma(4,9) = C131211;
        dMdGamma(4,10) = C131222;
        dMdGamma(4,11) = C131233;
        dMdGamma(4,12) = C131223;
        dMdGamma(4,13) = C131213;
        dMdGamma(4,14) = C131212;
        dMdGamma(4,15) = C131232;
        dMdGamma(4,16) = C131231;
        dMdGamma(4,17) = C131221;
        dMdGamma(4,18) = C131311;
        dMdGamma(4,19) = C131322;
        dMdGamma(4,20) = C131333;
        dMdGamma(4,21) = C131323;
        dMdGamma(4,22) = C131313;
        dMdGamma(4,23) = C131312;
        dMdGamma(4,24) = C131332;
        dMdGamma(4,25) = C131331;
        dMdGamma(4,26) = C131321;
        dMdGamma(5,0) = C121111;
        dMdGamma(5,1) = C121122;
        dMdGamma(5,2) = C121133;
        dMdGamma(5,3) = C121123;
        dMdGamma(5,4) = C121113;
        dMdGamma(5,5) = C121112;
        dMdGamma(5,6) = C121132;
        dMdGamma(5,7) = C121131;
        dMdGamma(5,8) = C121121;
        dMdGamma(5,9) = C121211;
        dMdGamma(5,10) = C121222;
        dMdGamma(5,11) = C121233;
        dMdGamma(5,12) = C121223;
        dMdGamma(5,13) = C121213;
        dMdGamma(5,14) = C121212;
        dMdGamma(5,15) = C121232;
        dMdGamma(5,16) = C121231;
        dMdGamma(5,17) = C121221;
        dMdGamma(5,18) = C121311;
        dMdGamma(5,19) = C121322;
        dMdGamma(5,20) = C121333;
        dMdGamma(5,21) = C121323;
        dMdGamma(5,22) = C121313;
        dMdGamma(5,23) = C121312;
        dMdGamma(5,24) = C121332;
        dMdGamma(5,25) = C121331;
        dMdGamma(5,26) = C121321;
        dMdGamma(6,0) = C321111;
        dMdGamma(6,1) = C321122;
        dMdGamma(6,2) = C321133;
        dMdGamma(6,3) = C321123;
        dMdGamma(6,4) = C321113;
        dMdGamma(6,5) = C321112;
        dMdGamma(6,6) = C321132;
        dMdGamma(6,7) = C321131;
        dMdGamma(6,8) = C321121;
        dMdGamma(6,9) = C321211;
        dMdGamma(6,10) = C321222;
        dMdGamma(6,11) = C321233;
        dMdGamma(6,12) = C321223;
        dMdGamma(6,13) = C321213;
        dMdGamma(6,14) = C321212;
        dMdGamma(6,15) = C321232;
        dMdGamma(6,16) = C321231;
        dMdGamma(6,17) = C321221;
        dMdGamma(6,18) = C321311;
        dMdGamma(6,19) = C321322;
        dMdGamma(6,20) = C321333;
        dMdGamma(6,21) = C321323;
        dMdGamma(6,22) = C321313;
        dMdGamma(6,23) = C321312;
        dMdGamma(6,24) = C321332;
        dMdGamma(6,25) = C321331;
        dMdGamma(6,26) = C321321;
        dMdGamma(7,0) = C311111;
        dMdGamma(7,1) = C311122;
        dMdGamma(7,2) = C311133;
        dMdGamma(7,3) = C311123;
        dMdGamma(7,4) = C311113;
        dMdGamma(7,5) = C311112;
        dMdGamma(7,6) = C311132;
        dMdGamma(7,7) = C311131;
        dMdGamma(7,8) = C311121;
        dMdGamma(7,9) = C311211;
        dMdGamma(7,10) = C311222;
        dMdGamma(7,11) = C311233;
        dMdGamma(7,12) = C311223;
        dMdGamma(7,13) = C311213;
        dMdGamma(7,14) = C311212;
        dMdGamma(7,15) = C311232;
        dMdGamma(7,16) = C311231;
        dMdGamma(7,17) = C311221;
        dMdGamma(7,18) = C311311;
        dMdGamma(7,19) = C311322;
        dMdGamma(7,20) = C311333;
        dMdGamma(7,21) = C311323;
        dMdGamma(7,22) = C311313;
        dMdGamma(7,23) = C311312;
        dMdGamma(7,24) = C311332;
        dMdGamma(7,25) = C311331;
        dMdGamma(7,26) = C311321;
        dMdGamma(8,0) = C211111;
        dMdGamma(8,1) = C211122;
        dMdGamma(8,2) = C211133;
        dMdGamma(8,3) = C211123;
        dMdGamma(8,4) = C211113;
        dMdGamma(8,5) = C211112;
        dMdGamma(8,6) = C211132;
        dMdGamma(8,7) = C211131;
        dMdGamma(8,8) = C211121;
        dMdGamma(8,9) = C211211;
        dMdGamma(8,10) = C211222;
        dMdGamma(8,11) = C211233;
        dMdGamma(8,12) = C211223;
        dMdGamma(8,13) = C211213;
        dMdGamma(8,14) = C211212;
        dMdGamma(8,15) = C211232;
        dMdGamma(8,16) = C211231;
        dMdGamma(8,17) = C211221;
        dMdGamma(8,18) = C211311;
        dMdGamma(8,19) = C211322;
        dMdGamma(8,20) = C211333;
        dMdGamma(8,21) = C211323;
        dMdGamma(8,22) = C211313;
        dMdGamma(8,23) = C211312;
        dMdGamma(8,24) = C211332;
        dMdGamma(8,25) = C211331;
        dMdGamma(8,26) = C211321;
        dMdGamma(9,0) = C112111;
        dMdGamma(9,1) = C112122;
        dMdGamma(9,2) = C112133;
        dMdGamma(9,3) = C112123;
        dMdGamma(9,4) = C112113;
        dMdGamma(9,5) = C112112;
        dMdGamma(9,6) = C112132;
        dMdGamma(9,7) = C112131;
        dMdGamma(9,8) = C112121;
        dMdGamma(9,9) = C112211;
        dMdGamma(9,10) = C112222;
        dMdGamma(9,11) = C112233;
        dMdGamma(9,12) = C112223;
        dMdGamma(9,13) = C112213;
        dMdGamma(9,14) = C112212;
        dMdGamma(9,15) = C112232;
        dMdGamma(9,16) = C112231;
        dMdGamma(9,17) = C112221;
        dMdGamma(9,18) = C112311;
        dMdGamma(9,19) = C112322;
        dMdGamma(9,20) = C112333;
        dMdGamma(9,21) = C112323;
        dMdGamma(9,22) = C112313;
        dMdGamma(9,23) = C112312;
        dMdGamma(9,24) = C112332;
        dMdGamma(9,25) = C112331;
        dMdGamma(9,26) = C112321;
        dMdGamma(10,0) = C222111;
        dMdGamma(10,1) = C222122;
        dMdGamma(10,2) = C222133;
        dMdGamma(10,3) = C222123;
        dMdGamma(10,4) = C222113;
        dMdGamma(10,5) = C222112;
        dMdGamma(10,6) = C222132;
        dMdGamma(10,7) = C222131;
        dMdGamma(10,8) = C222121;
        dMdGamma(10,9) = C222211;
        dMdGamma(10,10) = C222222;
        dMdGamma(10,11) = C222233;
        dMdGamma(10,12) = C222223;
        dMdGamma(10,13) = C222213;
        dMdGamma(10,14) = C222212;
        dMdGamma(10,15) = C222232;
        dMdGamma(10,16) = C222231;
        dMdGamma(10,17) = C222221;
        dMdGamma(10,18) = C222311;
        dMdGamma(10,19) = C222322;
        dMdGamma(10,20) = C222333;
        dMdGamma(10,21) = C222323;
        dMdGamma(10,22) = C222313;
        dMdGamma(10,23) = C222312;
        dMdGamma(10,24) = C222332;
        dMdGamma(10,25) = C222331;
        dMdGamma(10,26) = C222321;
        dMdGamma(11,0) = C332111;
        dMdGamma(11,1) = C332122;
        dMdGamma(11,2) = C332133;
        dMdGamma(11,3) = C332123;
        dMdGamma(11,4) = C332113;
        dMdGamma(11,5) = C332112;
        dMdGamma(11,6) = C332132;
        dMdGamma(11,7) = C332131;
        dMdGamma(11,8) = C332121;
        dMdGamma(11,9) = C332211;
        dMdGamma(11,10) = C332222;
        dMdGamma(11,11) = C332233;
        dMdGamma(11,12) = C332223;
        dMdGamma(11,13) = C332213;
        dMdGamma(11,14) = C332212;
        dMdGamma(11,15) = C332232;
        dMdGamma(11,16) = C332231;
        dMdGamma(11,17) = C332221;
        dMdGamma(11,18) = C332311;
        dMdGamma(11,19) = C332322;
        dMdGamma(11,20) = C332333;
        dMdGamma(11,21) = C332323;
        dMdGamma(11,22) = C332313;
        dMdGamma(11,23) = C332312;
        dMdGamma(11,24) = C332332;
        dMdGamma(11,25) = C332331;
        dMdGamma(11,26) = C332321;
        dMdGamma(12,0) = C232111;
        dMdGamma(12,1) = C232122;
        dMdGamma(12,2) = C232133;
        dMdGamma(12,3) = C232123;
        dMdGamma(12,4) = C232113;
        dMdGamma(12,5) = C232112;
        dMdGamma(12,6) = C232132;
        dMdGamma(12,7) = C232131;
        dMdGamma(12,8) = C232121;
        dMdGamma(12,9) = C232211;
        dMdGamma(12,10) = C232222;
        dMdGamma(12,11) = C232233;
        dMdGamma(12,12) = C232223;
        dMdGamma(12,13) = C232213;
        dMdGamma(12,14) = C232212;
        dMdGamma(12,15) = C232232;
        dMdGamma(12,16) = C232231;
        dMdGamma(12,17) = C232221;
        dMdGamma(12,18) = C232311;
        dMdGamma(12,19) = C232322;
        dMdGamma(12,20) = C232333;
        dMdGamma(12,21) = C232323;
        dMdGamma(12,22) = C232313;
        dMdGamma(12,23) = C232312;
        dMdGamma(12,24) = C232332;
        dMdGamma(12,25) = C232331;
        dMdGamma(12,26) = C232321;
        dMdGamma(13,0) = C132111;
        dMdGamma(13,1) = C132122;
        dMdGamma(13,2) = C132133;
        dMdGamma(13,3) = C132123;
        dMdGamma(13,4) = C132113;
        dMdGamma(13,5) = C132112;
        dMdGamma(13,6) = C132132;
        dMdGamma(13,7) = C132131;
        dMdGamma(13,8) = C132121;
        dMdGamma(13,9) = C132211;
        dMdGamma(13,10) = C132222;
        dMdGamma(13,11) = C132233;
        dMdGamma(13,12) = C132223;
        dMdGamma(13,13) = C132213;
        dMdGamma(13,14) = C132212;
        dMdGamma(13,15) = C132232;
        dMdGamma(13,16) = C132231;
        dMdGamma(13,17) = C132221;
        dMdGamma(13,18) = C132311;
        dMdGamma(13,19) = C132322;
        dMdGamma(13,20) = C132333;
        dMdGamma(13,21) = C132323;
        dMdGamma(13,22) = C132313;
        dMdGamma(13,23) = C132312;
        dMdGamma(13,24) = C132332;
        dMdGamma(13,25) = C132331;
        dMdGamma(13,26) = C132321;
        dMdGamma(14,0) = C122111;
        dMdGamma(14,1) = C122122;
        dMdGamma(14,2) = C122133;
        dMdGamma(14,3) = C122123;
        dMdGamma(14,4) = C122113;
        dMdGamma(14,5) = C122112;
        dMdGamma(14,6) = C122132;
        dMdGamma(14,7) = C122131;
        dMdGamma(14,8) = C122121;
        dMdGamma(14,9) = C122211;
        dMdGamma(14,10) = C122222;
        dMdGamma(14,11) = C122233;
        dMdGamma(14,12) = C122223;
        dMdGamma(14,13) = C122213;
        dMdGamma(14,14) = C122212;
        dMdGamma(14,15) = C122232;
        dMdGamma(14,16) = C122231;
        dMdGamma(14,17) = C122221;
        dMdGamma(14,18) = C122311;
        dMdGamma(14,19) = C122322;
        dMdGamma(14,20) = C122333;
        dMdGamma(14,21) = C122323;
        dMdGamma(14,22) = C122313;
        dMdGamma(14,23) = C122312;
        dMdGamma(14,24) = C122332;
        dMdGamma(14,25) = C122331;
        dMdGamma(14,26) = C122321;
        dMdGamma(15,0) = C322111;
        dMdGamma(15,1) = C322122;
        dMdGamma(15,2) = C322133;
        dMdGamma(15,3) = C322123;
        dMdGamma(15,4) = C322113;
        dMdGamma(15,5) = C322112;
        dMdGamma(15,6) = C322132;
        dMdGamma(15,7) = C322131;
        dMdGamma(15,8) = C322121;
        dMdGamma(15,9) = C322211;
        dMdGamma(15,10) = C322222;
        dMdGamma(15,11) = C322233;
        dMdGamma(15,12) = C322223;
        dMdGamma(15,13) = C322213;
        dMdGamma(15,14) = C322212;
        dMdGamma(15,15) = C322232;
        dMdGamma(15,16) = C322231;
        dMdGamma(15,17) = C322221;
        dMdGamma(15,18) = C322311;
        dMdGamma(15,19) = C322322;
        dMdGamma(15,20) = C322333;
        dMdGamma(15,21) = C322323;
        dMdGamma(15,22) = C322313;
        dMdGamma(15,23) = C322312;
        dMdGamma(15,24) = C322332;
        dMdGamma(15,25) = C322331;
        dMdGamma(15,26) = C322321;
        dMdGamma(16,0) = C312111;
        dMdGamma(16,1) = C312122;
        dMdGamma(16,2) = C312133;
        dMdGamma(16,3) = C312123;
        dMdGamma(16,4) = C312113;
        dMdGamma(16,5) = C312112;
        dMdGamma(16,6) = C312132;
        dMdGamma(16,7) = C312131;
        dMdGamma(16,8) = C312121;
        dMdGamma(16,9) = C312211;
        dMdGamma(16,10) = C312222;
        dMdGamma(16,11) = C312233;
        dMdGamma(16,12) = C312223;
        dMdGamma(16,13) = C312213;
        dMdGamma(16,14) = C312212;
        dMdGamma(16,15) = C312232;
        dMdGamma(16,16) = C312231;
        dMdGamma(16,17) = C312221;
        dMdGamma(16,18) = C312311;
        dMdGamma(16,19) = C312322;
        dMdGamma(16,20) = C312333;
        dMdGamma(16,21) = C312323;
        dMdGamma(16,22) = C312313;
        dMdGamma(16,23) = C312312;
        dMdGamma(16,24) = C312332;
        dMdGamma(16,25) = C312331;
        dMdGamma(16,26) = C312321;
        dMdGamma(17,0) = C212111;
        dMdGamma(17,1) = C212122;
        dMdGamma(17,2) = C212133;
        dMdGamma(17,3) = C212123;
        dMdGamma(17,4) = C212113;
        dMdGamma(17,5) = C212112;
        dMdGamma(17,6) = C212132;
        dMdGamma(17,7) = C212131;
        dMdGamma(17,8) = C212121;
        dMdGamma(17,9) = C212211;
        dMdGamma(17,10) = C212222;
        dMdGamma(17,11) = C212233;
        dMdGamma(17,12) = C212223;
        dMdGamma(17,13) = C212213;
        dMdGamma(17,14) = C212212;
        dMdGamma(17,15) = C212232;
        dMdGamma(17,16) = C212231;
        dMdGamma(17,17) = C212221;
        dMdGamma(17,18) = C212311;
        dMdGamma(17,19) = C212322;
        dMdGamma(17,20) = C212333;
        dMdGamma(17,21) = C212323;
        dMdGamma(17,22) = C212313;
        dMdGamma(17,23) = C212312;
        dMdGamma(17,24) = C212332;
        dMdGamma(17,25) = C212331;
        dMdGamma(17,26) = C212321;
        dMdGamma(18,0) = C113111;
        dMdGamma(18,1) = C113122;
        dMdGamma(18,2) = C113133;
        dMdGamma(18,3) = C113123;
        dMdGamma(18,4) = C113113;
        dMdGamma(18,5) = C113112;
        dMdGamma(18,6) = C113132;
        dMdGamma(18,7) = C113131;
        dMdGamma(18,8) = C113121;
        dMdGamma(18,9) = C113211;
        dMdGamma(18,10) = C113222;
        dMdGamma(18,11) = C113233;
        dMdGamma(18,12) = C113223;
        dMdGamma(18,13) = C113213;
        dMdGamma(18,14) = C113212;
        dMdGamma(18,15) = C113232;
        dMdGamma(18,16) = C113231;
        dMdGamma(18,17) = C113221;
        dMdGamma(18,18) = C113311;
        dMdGamma(18,19) = C113322;
        dMdGamma(18,20) = C113333;
        dMdGamma(18,21) = C113323;
        dMdGamma(18,22) = C113313;
        dMdGamma(18,23) = C113312;
        dMdGamma(18,24) = C113332;
        dMdGamma(18,25) = C113331;
        dMdGamma(18,26) = C113321;
        dMdGamma(19,0) = C223111;
        dMdGamma(19,1) = C223122;
        dMdGamma(19,2) = C223133;
        dMdGamma(19,3) = C223123;
        dMdGamma(19,4) = C223113;
        dMdGamma(19,5) = C223112;
        dMdGamma(19,6) = C223132;
        dMdGamma(19,7) = C223131;
        dMdGamma(19,8) = C223121;
        dMdGamma(19,9) = C223211;
        dMdGamma(19,10) = C223222;
        dMdGamma(19,11) = C223233;
        dMdGamma(19,12) = C223223;
        dMdGamma(19,13) = C223213;
        dMdGamma(19,14) = C223212;
        dMdGamma(19,15) = C223232;
        dMdGamma(19,16) = C223231;
        dMdGamma(19,17) = C223221;
        dMdGamma(19,18) = C223311;
        dMdGamma(19,19) = C223322;
        dMdGamma(19,20) = C223333;
        dMdGamma(19,21) = C223323;
        dMdGamma(19,22) = C223313;
        dMdGamma(19,23) = C223312;
        dMdGamma(19,24) = C223332;
        dMdGamma(19,25) = C223331;
        dMdGamma(19,26) = C223321;
        dMdGamma(20,0) = C333111;
        dMdGamma(20,1) = C333122;
        dMdGamma(20,2) = C333133;
        dMdGamma(20,3) = C333123;
        dMdGamma(20,4) = C333113;
        dMdGamma(20,5) = C333112;
        dMdGamma(20,6) = C333132;
        dMdGamma(20,7) = C333131;
        dMdGamma(20,8) = C333121;
        dMdGamma(20,9) = C333211;
        dMdGamma(20,10) = C333222;
        dMdGamma(20,11) = C333233;
        dMdGamma(20,12) = C333223;
        dMdGamma(20,13) = C333213;
        dMdGamma(20,14) = C333212;
        dMdGamma(20,15) = C333232;
        dMdGamma(20,16) = C333231;
        dMdGamma(20,17) = C333221;
        dMdGamma(20,18) = C333311;
        dMdGamma(20,19) = C333322;
        dMdGamma(20,20) = C333333;
        dMdGamma(20,21) = C333323;
        dMdGamma(20,22) = C333313;
        dMdGamma(20,23) = C333312;
        dMdGamma(20,24) = C333332;
        dMdGamma(20,25) = C333331;
        dMdGamma(20,26) = C333321;
        dMdGamma(21,0) = C233111;
        dMdGamma(21,1) = C233122;
        dMdGamma(21,2) = C233133;
        dMdGamma(21,3) = C233123;
        dMdGamma(21,4) = C233113;
        dMdGamma(21,5) = C233112;
        dMdGamma(21,6) = C233132;
        dMdGamma(21,7) = C233131;
        dMdGamma(21,8) = C233121;
        dMdGamma(21,9) = C233211;
        dMdGamma(21,10) = C233222;
        dMdGamma(21,11) = C233233;
        dMdGamma(21,12) = C233223;
        dMdGamma(21,13) = C233213;
        dMdGamma(21,14) = C233212;
        dMdGamma(21,15) = C233232;
        dMdGamma(21,16) = C233231;
        dMdGamma(21,17) = C233221;
        dMdGamma(21,18) = C233311;
        dMdGamma(21,19) = C233322;
        dMdGamma(21,20) = C233333;
        dMdGamma(21,21) = C233323;
        dMdGamma(21,22) = C233313;
        dMdGamma(21,23) = C233312;
        dMdGamma(21,24) = C233332;
        dMdGamma(21,25) = C233331;
        dMdGamma(21,26) = C233321;
        dMdGamma(22,0) = C133111;
        dMdGamma(22,1) = C133122;
        dMdGamma(22,2) = C133133;
        dMdGamma(22,3) = C133123;
        dMdGamma(22,4) = C133113;
        dMdGamma(22,5) = C133112;
        dMdGamma(22,6) = C133132;
        dMdGamma(22,7) = C133131;
        dMdGamma(22,8) = C133121;
        dMdGamma(22,9) = C133211;
        dMdGamma(22,10) = C133222;
        dMdGamma(22,11) = C133233;
        dMdGamma(22,12) = C133223;
        dMdGamma(22,13) = C133213;
        dMdGamma(22,14) = C133212;
        dMdGamma(22,15) = C133232;
        dMdGamma(22,16) = C133231;
        dMdGamma(22,17) = C133221;
        dMdGamma(22,18) = C133311;
        dMdGamma(22,19) = C133322;
        dMdGamma(22,20) = C133333;
        dMdGamma(22,21) = C133323;
        dMdGamma(22,22) = C133313;
        dMdGamma(22,23) = C133312;
        dMdGamma(22,24) = C133332;
        dMdGamma(22,25) = C133331;
        dMdGamma(22,26) = C133321;
        dMdGamma(23,0) = C123111;
        dMdGamma(23,1) = C123122;
        dMdGamma(23,2) = C123133;
        dMdGamma(23,3) = C123123;
        dMdGamma(23,4) = C123113;
        dMdGamma(23,5) = C123112;
        dMdGamma(23,6) = C123132;
        dMdGamma(23,7) = C123131;
        dMdGamma(23,8) = C123121;
        dMdGamma(23,9) = C123211;
        dMdGamma(23,10) = C123222;
        dMdGamma(23,11) = C123233;
        dMdGamma(23,12) = C123223;
        dMdGamma(23,13) = C123213;
        dMdGamma(23,14) = C123212;
        dMdGamma(23,15) = C123232;
        dMdGamma(23,16) = C123231;
        dMdGamma(23,17) = C123221;
        dMdGamma(23,18) = C123311;
        dMdGamma(23,19) = C123322;
        dMdGamma(23,20) = C123333;
        dMdGamma(23,21) = C123323;
        dMdGamma(23,22) = C123313;
        dMdGamma(23,23) = C123312;
        dMdGamma(23,24) = C123332;
        dMdGamma(23,25) = C123331;
        dMdGamma(23,26) = C123321;
        dMdGamma(24,0) = C323111;
        dMdGamma(24,1) = C323122;
        dMdGamma(24,2) = C323133;
        dMdGamma(24,3) = C323123;
        dMdGamma(24,4) = C323113;
        dMdGamma(24,5) = C323112;
        dMdGamma(24,6) = C323132;
        dMdGamma(24,7) = C323131;
        dMdGamma(24,8) = C323121;
        dMdGamma(24,9) = C323211;
        dMdGamma(24,10) = C323222;
        dMdGamma(24,11) = C323233;
        dMdGamma(24,12) = C323223;
        dMdGamma(24,13) = C323213;
        dMdGamma(24,14) = C323212;
        dMdGamma(24,15) = C323232;
        dMdGamma(24,16) = C323231;
        dMdGamma(24,17) = C323221;
        dMdGamma(24,18) = C323311;
        dMdGamma(24,19) = C323322;
        dMdGamma(24,20) = C323333;
        dMdGamma(24,21) = C323323;
        dMdGamma(24,22) = C323313;
        dMdGamma(24,23) = C323312;
        dMdGamma(24,24) = C323332;
        dMdGamma(24,25) = C323331;
        dMdGamma(24,26) = C323321;
        dMdGamma(25,0) = C313111;
        dMdGamma(25,1) = C313122;
        dMdGamma(25,2) = C313133;
        dMdGamma(25,3) = C313123;
        dMdGamma(25,4) = C313113;
        dMdGamma(25,5) = C313112;
        dMdGamma(25,6) = C313132;
        dMdGamma(25,7) = C313131;
        dMdGamma(25,8) = C313121;
        dMdGamma(25,9) = C313211;
        dMdGamma(25,10) = C313222;
        dMdGamma(25,11) = C313233;
        dMdGamma(25,12) = C313223;
        dMdGamma(25,13) = C313213;
        dMdGamma(25,14) = C313212;
        dMdGamma(25,15) = C313232;
        dMdGamma(25,16) = C313231;
        dMdGamma(25,17) = C313221;
        dMdGamma(25,18) = C313311;
        dMdGamma(25,19) = C313322;
        dMdGamma(25,20) = C313333;
        dMdGamma(25,21) = C313323;
        dMdGamma(25,22) = C313313;
        dMdGamma(25,23) = C313312;
        dMdGamma(25,24) = C313332;
        dMdGamma(25,25) = C313331;
        dMdGamma(25,26) = C313321;
        dMdGamma(26,0) = C213111;
        dMdGamma(26,1) = C213122;
        dMdGamma(26,2) = C213133;
        dMdGamma(26,3) = C213123;
        dMdGamma(26,4) = C213113;
        dMdGamma(26,5) = C213112;
        dMdGamma(26,6) = C213132;
        dMdGamma(26,7) = C213131;
        dMdGamma(26,8) = C213121;
        dMdGamma(26,9) = C213211;
        dMdGamma(26,10) = C213222;
        dMdGamma(26,11) = C213233;
        dMdGamma(26,12) = C213223;
        dMdGamma(26,13) = C213213;
        dMdGamma(26,14) = C213212;
        dMdGamma(26,15) = C213232;
        dMdGamma(26,16) = C213231;
        dMdGamma(26,17) = C213221;
        dMdGamma(26,18) = C213311;
        dMdGamma(26,19) = C213322;
        dMdGamma(26,20) = C213333;
        dMdGamma(26,21) = C213323;
        dMdGamma(26,22) = C213313;
        dMdGamma(26,23) = C213312;
        dMdGamma(26,24) = C213332;
        dMdGamma(26,25) = C213331;
        dMdGamma(26,26) = C213321;
        return;
    }
}
