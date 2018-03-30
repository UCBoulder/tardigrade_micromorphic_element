/*!===========================================================================
   |                                                                         |
   |                        micromorphic_measures.cpp                        |
   |                                                                         |
   ===========================================================================
   | The source file for a wrapper that converts variables and their         |
   | gradients generated by MOOSE into deformation measures which can be     |
   | used to compute the micromorphic stress measures and their tangents.    |
   | This is done to avoid any possible assumptions of symmetry which could  |
   | be present in the Tensor Mechanics physics module.                      |
   ===========================================================================
   | Dependencies:                                                           |
   |     Eigen: A matrix library available at eigen.tuxfamily.org            |
   ===========================================================================
   */

#include<micromorphic_measures.h>

namespace micromorphic_measures

{

    void get_deformation_gradient(const double (&_grad_u)[3][3], Matrix_3x3 &F){

        /*!==================================
           |    get_deformation_gradient    |
           ==================================

           Compute the deformation gradient from the 
           gradients of the displacements. The deformation 
           gradient can be computed via

           Finv(i,j) = eye(i,j) - _grad_u[i][j]

           where Finv is the inverse of the deformation gradient.

        */

        for (int i=0; i<3; i++){

            for (int j=0; j<3; j++){

                F(i,j) = -_grad_u[i][j];

            }

        }

        for (int i=0; i<3; i++){F(i,i) += 1;}

        F = F.inverse().eval();

        return;
    }

    void assemble_chi(const double (&_phi)[9], Matrix_3x3 &chi){

        /*!======================
           |    assemble_chi    |
           ======================

           Assemble the micro-deformation 
           tensor chi from the phi degrees of 
           freedom.

           _phi is assumed to be organized in 
           Voigt notation form i.e.

           phi_11, phi_22, phi_33, phi_23, phi_13, phi_12, phi_32, phi_31, phi_21

        */

        chi(0,0) = 1. + _phi[0];
        chi(1,1) = 1. + _phi[1];
        chi(2,2) = 1. + _phi[2];
        chi(1,2) = _phi[3];
        chi(0,2) = _phi[4];
        chi(0,1) = _phi[5];
        chi(2,1) = _phi[6];
        chi(2,0) = _phi[7];
        chi(1,0) = _phi[8];

    }

    void assemble_grad_chi(const double (&_grad_phi)[9][3], Matrix_3x9 &grad_chi){

        /*!==========================
          |    assemble_grad_chi    |
          ===========================

          Assemble the gradient of chi w.r.t. the 
          current coordinates.

          _grad_phi is assumed to be organized

          _grad_phi[I][k] = phi_ij,k

          where I is a ``super'' index which has the form

               0  1  2  3  4  5  6  7  8
          I = 11,22,33,23,13,12,32,31,21

          and k ranges from

              0  1  2
          k = 1, 2, 3

        */

        grad_chi(0,0) = _grad_phi[0][0]; //111
        grad_chi(0,1) = _grad_phi[5][1]; //122
        grad_chi(0,2) = _grad_phi[4][2]; //133
        grad_chi(0,3) = _grad_phi[5][2]; //123
        grad_chi(0,4) = _grad_phi[0][2]; //113
        grad_chi(0,5) = _grad_phi[0][1]; //112
        grad_chi(0,6) = _grad_phi[4][1]; //132
        grad_chi(0,7) = _grad_phi[4][0]; //131
        grad_chi(0,8) = _grad_phi[5][0]; //121
        
        grad_chi(1,0) = _grad_phi[8][0]; //211
        grad_chi(1,1) = _grad_phi[1][1]; //222
        grad_chi(1,2) = _grad_phi[3][2]; //233
        grad_chi(1,3) = _grad_phi[1][2]; //223
        grad_chi(1,4) = _grad_phi[8][2]; //213
        grad_chi(1,5) = _grad_phi[8][1]; //212
        grad_chi(1,6) = _grad_phi[3][1]; //232
        grad_chi(1,7) = _grad_phi[3][0]; //231
        grad_chi(1,8) = _grad_phi[1][0]; //221

        grad_chi(2,0) = _grad_phi[7][0]; //311
        grad_chi(2,1) = _grad_phi[6][1]; //322
        grad_chi(2,2) = _grad_phi[2][2]; //333
        grad_chi(2,3) = _grad_phi[6][2]; //323
        grad_chi(2,4) = _grad_phi[7][2]; //313
        grad_chi(2,5) = _grad_phi[7][1]; //312
        grad_chi(2,6) = _grad_phi[2][1]; //332
        grad_chi(2,7) = _grad_phi[2][0]; //331
        grad_chi(2,8) = _grad_phi[6][0]; //321

        return;

    }

    void get_right_cauchy_green(const Matrix_3x3 &F, Matrix_3x3 &RCG){
        /*!================================
           |    get_right_cauchy_green    |
           ================================

           Get the right Cauchy-Green deformation 
           tensor.

           C_IJ = F_iI F_iJ

        */

        RCG = F.transpose()*F;
        return;
    }

    void get_left_cauchy_green(const Matrix_3x3 &F, Matrix_3x3 &LCG){

        /*!===============================
           |    get_left_cauchy_green    |
           ===============================

           Get the left Cauchy-Green deformation
           tensor.

           b_ij = F_iI F_jI

        */

        LCG = F*F.transpose();
        return;
    }

    void get_lagrange_strain(const Matrix_3x3 &F, Matrix_3x3& E){

        /*!=============================
           |    get_lagrange_strain    |
           =============================

           Get the lagrange strain tensor.

           E_IJ = 0.5*(F_iI F_iJ - I_IJ)

        */

        E = F.transpose()*F;
        for (int i=0; i<3; i++){E(i,i) -= 1.;}
        return;
    }

    void get_almansi_strain(const Matrix_3x3 &F, Matrix_3x3 &e){
        /*!============================
           |    get_almansi_strain    |
           ============================

           Compute the almansi strain.

           e_ij = 0.5*(I_ij - (F_iI F_jI)**-1)

        */

        e = -(F*F.transpose()).inverse();
        for (int i=0; i<3; i++){e(i,i) += 1.;}
        return;
    }

    void get_small_strain(const double (&_grad_u)[3][3], Matrix_3x3 &epsilon){
         /*!=========================
           |    get_small_strain    |
           ==========================

           Compute the small strain tensor.

           epsilon(i,j) = 0.5*(_grad_u[i][j] + _grad_u[j][i])

         */

         for (int i=0; i<3; i++){

             for (int j=0; j<3; j++){

                 epsilon(i,j) = 0.5*(_grad_u[i][j] + _grad_u[j][i]);

             }
         }
             
         return;
    }

}
