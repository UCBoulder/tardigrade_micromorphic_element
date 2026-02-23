/*!============================================================
  |                                                           |
  |                    balance_equations.h                    |
  |                                                           |
  -------------------------------------------------------------
  | The header file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  -------------------------------------------------------------
  | Note: N is the test function and eta is the interpolation |
  |       function. We use this notation since these symbols  |
  |       are not used by the micromorphic formulation.       |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

// Required to fix linking problems
//__asm__(".symver memcpy,memcpy@GLIBC_2.14");

#ifndef BALANCE_EQUATIONS_H
#define BALANCE_EQUATIONS_H

#define USE_EIGEN
#include <vector>

namespace balance_equations {

    typedef double                      variableType;
    typedef std::vector<variableType>   variableVector;
    typedef std::vector<variableVector> variableMatrix;

    /*============================================
    | Forces from the balance of linear momentum |
    ============================================*/

    int compute_internal_force(const double (&dNdX)[3], const variableVector &F, const variableVector &PK2,
                               double (&fint)[3]);

    int compute_internal_force(const unsigned int &i, const double (&dNdX)[3], const variableVector &F,
                               const variableVector &PK2, double &fint_i);

    void compute_body_force(const double &N, const double &density, const double (&b)[3], double (&fb)[3]);

    void compute_body_force(const unsigned int &i, const double &N, const double &density, const double (&b)[3],
                            double &fb);

    int compute_inertia_force(const double &N, const double &density, const double (&a)[3], double (&finertia)[3]);

    int compute_inertia_force(const unsigned int &i, const double &N, const double &density, const double (&a)[3],
                              double &finertia_i);

    /*===========================================================
    | Stresses from the balance of the first moment of momentum |
    ===========================================================*/

    // Stresses from the balance of first moment of momentum
    int compute_internal_couple(const double &N, const double (&dNdX)[3], const variableVector &F,
                                const variableVector &chi, const variableVector &PK2, const variableVector &SIGMA,
                                const variableVector &M, double (&cint)[9]);

    int compute_internal_couple(const unsigned int &i, const unsigned int &j, const double &N, const double (&dNdX)[3],
                                const variableVector &F, const variableVector &chi, const variableVector &PK2,
                                const variableVector &SIGMA, const variableVector &M, double &cint_ij);

    void compute_body_couple(const double &N, const double &density, const double (&l)[9], double (&cb)[9]);

    void compute_body_couple(const unsigned int &i, const unsigned int &j, const double &N, const double &density,
                             const double (&l)[9], double &cb_ij);

    void compute_inertia_couple(const double &N, const double &density, const double (&omega)[9],
                                double (&cinertia)[9]);

    void compute_inertia_couple(const unsigned int &i, const unsigned int &j, const double &N, const double &density,
                                const double (&omega)[9], double &cinertia_ij);

    int compute_inertia_couple(const double &N, const double &density, const double (&chi)[9],
                               const double (&D2ChiDt2)[9], const double (&referenceInertia)[9], double (&cinertia)[9]);

    int compute_inertia_couple(const unsigned int &i, const unsigned int &j, const double &N, const double &density,
                               const double (&chi)[9], const double (&D2ChiDt2)[9], const double (&referenceInertia)[9],
                               double &cinertia_ij);

    /*=======================================================
    | The Jacobians of the balance of linear momentum terms |
    =======================================================*/

    int compute_internal_force_jacobian(const double &N, const double (&dNdX)[3], const double &eta,
                                        const double (&detadX)[3], const variableVector &F, const variableVector &PK2,
                                        const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                        const variableMatrix &DPK2Dgrad_phi, variableMatrix &DfintDU);

    int compute_internal_force_jacobian(const unsigned int &i, const unsigned int &j, const double &N,
                                        const double (&dNdX)[3], const double &eta, const double (&detadX)[3],
                                        const variableVector &F, const variableVector &PK2,
                                        const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                        const variableMatrix &DPK2Dgrad_phi, variableType &DfintDU_ij);

    int compute_inertia_force_jacobian(const double &N, const double &eta, const double &density, const double (&a)[3],
                                       const variableMatrix &DaDu, variableMatrix &DfinertiaDU);

    int compute_inertia_force_jacobian(const unsigned int &i, const unsigned int &j, const double &N, const double &eta,
                                       const double &density, const double (&a)[3], const variableVector &DaDu_i,
                                       variableType &DfinertiaDU_ij);

    /*====================================================================
    | The Jacobians of the balance of the first moment of momentum terms |
    ====================================================================*/

    int compute_internal_couple_jacobian(const double &N, const double (&dNdX)[3], const double &eta,
                                         const double (&detadX)[3], const variableVector &F, const variableVector &chi,
                                         const variableVector &PK2, const variableVector &SIGMA,
                                         const variableVector &M, const variableMatrix &DPK2Dgrad_u,
                                         const variableMatrix &DPK2Dphi, const variableMatrix &DPK2Dgrad_phi,
                                         const variableMatrix &DSIGMADgrad_u, const variableMatrix &DSIGMADphi,
                                         const variableMatrix &DSIGMADgrad_phi, const variableMatrix &DMDgrad_u,
                                         const variableMatrix &DMDphi, const variableMatrix &DMDgrad_phi,
                                         variableMatrix &DcintDU);

    int compute_internal_couple_jacobian(const unsigned int &i, const unsigned int &j, const double &N,
                                         const double (&dNdX)[3], const double &eta, const double (&detadX)[3],
                                         const variableVector &F, const variableVector &chi, const variableVector &PK2,
                                         const variableVector &SIGMA, const variableVector &M,
                                         const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                         const variableMatrix &DPK2Dgrad_phi, const variableMatrix &DSIGMADgrad_u,
                                         const variableMatrix &DSIGMADphi, const variableMatrix &DSIGMADgrad_phi,
                                         const variableMatrix &DMDgrad_u, const variableMatrix &DMDphi,
                                         const variableMatrix &DMDgrad_phi, variableType &DcintDU_ij);

    int compute_inertia_couple_jacobian(const double &N, const double &eta, const double &density,
                                        const double (&chi)[9], const double (&D2ChiDt2)[9],
                                        const variableMatrix &D3ChiDt2dChi, const double (&referenceInertia)[9],
                                        variableMatrix       &DcinertiaDU);

    int compute_inertia_couple_jacobian(const unsigned int &i, const unsigned int &j, const double &N,
                                        const double &eta, const double &density, const double (&chi)[9],
                                        const double (&D2ChiDt2)[9], const variableVector &D3ChiDt2dChi_j,
                                        const double (&referenceInertia)[9], double       &DcinertiaDU_ij);
}  // namespace balance_equations

#endif
