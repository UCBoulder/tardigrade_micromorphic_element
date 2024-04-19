/*!============================================================
  |                                                           |
  |                   balance_equations.cpp                   |
  |                                                           |
  -------------------------------------------------------------
  | The source file for the definition of the balance         |
  | equations used for standard as well as micromorphic       |
  | continuum mechanics.                                      |
  |                                                           |
  | Note: Compile with SKIP_ERROR_HANDLING to ignore the      |
  |       error handling which adds some overhead. This is    |
  |       only recommended if the code is known to be         |
  |       functioning correctly and performance is the        |
  |       primary concern.                                    |
  =============================================================
  | Dependencies:                                             |
  | Eigen: Open source matrix library available at            |
  |        eigen.tuxfamily.org.                               |
  =============================================================*/

#include <balance_equations.h>

namespace balance_equations{

    int compute_internal_force( const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &PK2, double ( &fint )[ 3 ]){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress.
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * Note: This has already accounted for the volume change term.
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, and 2 if PK2 is the wrong size.
         *
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double ( &fint )[ 3 ]: The internal force.
         */

        int errorCode;
        for ( unsigned int i = 0; i < 3; i++ ){
            errorCode = compute_internal_force( i, dNdX, F, PK2, fint[ i ] );

#ifndef SKIP_ERROR_HANDLING
            if ( errorCode != 0 ){
                return errorCode;
            }
#endif
        }

        return 0;
    }
    
    int compute_internal_force( const unsigned int &i,
                                const double (&dNdX)[3], const variableVector &F, const variableVector &PK2, double &fint_i ){
        /*!
         * Compute the internal force given the gradient of the shape function in the reference configuration,
         * the deformation gradient, and the PK2 stress on a single component.
         *
         * fint_i = - N_{ , I } * PK2_{ I J } * F_{ i J }
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, and 2 if PK2 is the wrong size.
         *
         * :param const unsigned int &i: The component to compute the stress on.
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape-functions in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The PK2 stress.
         * :param double fint_i: The i'th component of the internal force.
         */

        //Assume 3D
        const unsigned int dim = 3;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( PK2.size() != dim * dim ){
            return 2;
        }
#endif

        fint_i = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                fint_i -= dNdX[ I ] * PK2[ dim * I + J ] * F[ dim * i + J ];
            }
        }

        return 0;
    }

    void compute_body_force( const double &N, const double &density, const double ( &b )[ 3 ], double ( &fb )[ 3 ] ){
        /*!
         * Compute the body force
         *
         * fb_i = N * \rho * b_i
         *
         * :param const double &N: The value of the shape function
         * :param const double &density: The density in the current configuration
         * :param const double ( &b )[ 3 ]: The body force in the current configuration
         * :param double ( &fb )[ 3 ]: The body force in the current configuration
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        for ( unsigned int i = 0; i < dim; i++ ){
            compute_body_force( i, N, density, b, fb[ i ] );
        }

        return;
    }
    
    void compute_body_force( const unsigned int &i, const double &N, const double &density, const double ( &b )[ 3 ], double &fb_i ){
        /*!
         * Compute a single component of the body force
         *
         * fb_i = N * \rho * b_i
         *
         * :param const unsigned int &i: The component to compute the body force on
         * :param const double &N: The value of the shape function
         * :param const double &density: The density in the current configuration
         * :param const double ( &b )[ 3 ]: The body force in the current configuration
         * :param double fb_i: The body force in the current configuration in direction i
         */
        
        fb_i = N * density * b[ i ];

        return;
    }

    int compute_inertia_force( const double &N, const double &density, const double ( &a )[ 3 ], double ( &finertia )[ 3 ] ){
        /*!
         * Compute the inertia force
         *
         * finertia_i = -N * \rho * a_i
         *
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param double ( &finertia )[ 3 ]: The inertia force
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        for ( unsigned int i = 0; i < dim; i++ ){
            compute_inertia_force( i, N, density, a, finertia[ i ] );
        }

        return 0;
    }

    int compute_inertia_force( const unsigned int &i,
                               const double &N, const double &density, const double ( &a )[ 3 ], double &finertia_i ){
        /*!
         * Compute the inertia force for the ith component.
         *
         * finertia_i = -N * \rho * a_i
         *
         * :param const unsigned int &i: The component to compute the inertia force on
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param double finertia_i: The inertia force in the current configuration in direction i
         */
       
        finertia_i = -N * density * a[ i ]; 

        return 0;
    }
    
    int compute_internal_couple( const double &N, const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &chi,
                                 const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                 double ( &cint )[ 9 ] ){
        /*!
         * Compute the internal couple defined as
         * cint_{ ij } = N F_{ jJ } ( PK2_{ JI } - SIGMA_{ JI } ) F_{ iI } - N_{ ,K } F_{ iI } \chi_{ jJ } M_{ KIJ }
         *
         * :param const double &N: The shape-function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference configuration
         * :param const variableVector &F: The deformation gradient
         * :param const variableVector &chi: The micro-deformation
         * :param const variableVector &PK2: The second Piola-Kirchoff stress
         * :param const variableVector &SIGMA: The symmetric micro stress in the reference configuration.
         * :param const variableVector &M: The higher order stress in the reference configuration.
         * :param double ( &cint )[ 9 ] ): The internal couple stored as
         *     [ cint_{ 11 }, cint_{ 12 }, cint_{ 13 }, cint_{ 21 }, cint_{ 22 }, cint_{ 23 }, cint_{ 31 }, cint_{ 32 }, cint_{ 33 } ]
         */

        //Assume 3d
        const unsigned int dim = 3;

        int errorCode;
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                errorCode = compute_internal_couple( i, j, N, dNdX, F, chi, PK2, SIGMA, M, cint[ dim * i + j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_couple( const unsigned int &i, const unsigned int &j,
                                 const double &N, const double ( &dNdX )[ 3 ], const variableVector &F, const variableVector &chi,
                                 const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                 double &cint_ij ){
        /*!
         * Compute the internal couple at index ij defined as
         * cint_{ ij } = N F_{ jJ } ( PK2_{ JI } - SIGMA_{ JI } ) F_{ iI } - N_{ ,K } F_{ iI } \chi_{ jJ } M_{ KIJ }
         *
         * Function returns 0 for no errors, 1 if F is the wrong size, 2 if chi is the wrong size,
         * 3 if PK2 is the wrong size, 4 if SIGMA is the wrong size, and 5 if M is the wrong size
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference configuration
         * :param const variableVector &F: The deformation gradient
         * :param const variableVector &chi: The micro-deformation
         * :param const variableVector &PK2: The second Piola-Kirchoff stress
         * :param const variableVector &SIGMA: The symmetric micro stress in the reference configuration.
         * :param const variableVector &M: The higher order stress in the reference configuration.
         * :param double cint_ij: The ij'th internal couple
         */

        //Assume 3D
        const unsigned int dim = 3;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( chi.size() != dim * dim ){
            return 2;
        }

        if ( PK2.size() != dim * dim ){
            return 3;
        }

        if ( SIGMA.size() != dim * dim ){
            return 4;
        }

        if ( M.size() != dim * dim * dim ){
            return 5;
        }
#endif

        cint_ij = 0;
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                cint_ij += N * F[ dim * j + J ] * ( PK2[ dim * J + I ] - SIGMA[ dim * J + I ] ) * F[ dim * i + I ];

                for ( unsigned int K = 0; K < dim; K++ ){
                    cint_ij -= dNdX[ K ] * F[ dim * i + I ] * chi[ dim * j + J ] * M[ dim * dim * K + dim * I + J ];
                }
            }
        }

        return 0;
    }
    
    void compute_body_couple( const double &N, const double &density, const double ( &l )[ 9 ], double ( &cb )[ 9 ] ){
        /*!
         * Compute the body couple term
         *
         * couple_body_{ ij } = N * density * l_{ ij }
         *
         * :param const double &N: The shape-function value.
         * :param const double &density: The density in the current configuration.
         * :param const double ( %l )[ 9 ]: The body-force couple in the current configuration.
         * :param const double ( &cb )[ 9 ]: The body force couple term.
         */
        
        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_body_couple( i, j, N, density, l, cb[ dim * i + j ] );
            }
        }
        
        return;
    }

    void compute_body_couple( const unsigned int &i, const unsigned int &j,
                              const double &N, const double &density, const double ( &l )[ 9 ], double &cb_ij ){
        /*!
         * Compute the body couple term for the indices i and j
         *
         * couple_body_{ ij } = N * density * l_{ ij }
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value.
         * :param const double &density: The density in the current configuration.
         * :param const double ( %l )[ 9 ]: The body-force couple in the current configuration.
         * :param const double &cb_ij: The ij'th body force couple term.
         */

        //Assume 3D
        const unsigned int dim = 3;
        
        cb_ij = N * density * l[ dim * i + j ];
        
        return;
    }
    
    void compute_inertia_couple( const double &N, const double &density, const double ( &omega )[ 9 ], double ( &cinertia )[ 9 ] ){
        /*!
         * Compute the inertia couple in the current configuration
         *
         * cinertia_{ ij } = -N * density * omega_{ ij }
         *
         * :param const double &N: The shape-function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &omega )[ 9 ]: The micro-inertia tensor in the current configuration.
         * :param const double ( &cinertia )[ 9 ]: The inertia couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                compute_inertia_couple( i, j, N, density, omega, cinertia[ dim * i + j ] );
            }
        }
        
        return;
    }
    
    void compute_inertia_couple( const unsigned int &i, const unsigned int &j,
                                 const double &N, const double &density, const double ( &omega )[ 9 ], double &cinertia_ij ){
        /*!
         * Compute the inertia couple in the current configuration for the indices i and j
         *
         * cinertia_{ ij } = -N * density * omega_{ ij }
         *
         * :param const unsigned int &i: The first index.
         * :param const unsigned int &j: The second index.
         * :param const double &N: The shape-function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &omega )[ 9 ]: The micro-inertia tensor in the current configuration.
         * :param const double &cinternal_ij: The ij'th inertia couple.
         */

        //Assume 3D
        const unsigned int dim = 3;

        cinertia_ij = - N * density * omega[ dim * i + j ];
        
        return;
    }

    int compute_inertia_couple( const double &N, const double &density,
                                const double ( &chi )[ 9 ], const double ( &D2ChiDt2 )[ 9 ],
                                const double ( &referenceInertia )[ 9 ], double ( &cinertia )[ 9 ] ){
        /*!
         * Compute the inertia couple in the reference configuration
         *
         * integrand_ij = -N \rho_0 I_{ IJ } \ddot{ \chi }_{iI} \chi_{jJ}
         *
         * :param const double &N: The shape function value
         * :param const double &density: The density in the reference configuration.
         * :param const double ( &chi )[ 9 ]: The micro deformation tensor
         * :param const double ( &D2ChiDt2 )[ 9 ]: The second temporal derivative of the micro deformation tensor
         * :param const double ( &referenceInertia )[ 9 ]: The moment of inertia in the reference configuration
         * :param double ( &cinertia )[ 9 ]: The inertia couple
         */

        //Assume 3D
        const unsigned int dim = 3;

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                compute_inertia_couple( i, j, N, density, chi, D2ChiDt2, referenceInertia, cinertia[ dim * i + j ] );

            }

        }

        return 0;

    }

    int compute_inertia_couple( const unsigned int &i, const unsigned int &j, const double &N, const double &density,
                                const double ( &chi )[ 9 ], const double ( &D2ChiDt2 )[ 9 ], const double ( &referenceInertia )[ 9 ],
                                double &cinertia_ij ){
        /*!
         * Compute the ij index of the inertia couple
         *
         * integrand_ij = -N \rho_0 I_{ IJ } \ddot{ \chi }_{iI} \chi_{jJ}
         *
         * :param const unsigned int &i: The row index of the inertia couple
         * :param const unsigned int &j: The column index of the inertia couple
         * :param const double &N: The shape function value
         * :param const double &density: The density in the reference configuration.
         * :param const double ( &chi )[ 9 ]: The micro deformation tensor
         * :param const double ( &D2ChiDt2 )[ 9 ]: The second temporal derivative of the micro deformation tensor
         * :param const double ( &referenceInertia )[ 9 ]: The moment of inertia in the reference configuration
         * :param double ( &cinertia_ij )[ 9 ]: The inertia couple
         */

        //Assume 3D
        const unsigned int dim = 3;

        cinertia_ij = 0;

        for ( unsigned int I = 0; I < dim; I++ ){

            for ( unsigned int J = 0; J < dim; J++ ){

                cinertia_ij -= N * D2ChiDt2[ dim * i + I ] * chi[ dim * j + J ] * density * referenceInertia[ dim * I + J ];

            }

        }

        return 0;

    }

    int compute_internal_force_jacobian( const double &N, const double ( &dNdX )[ 3 ], const double &eta, const double ( &detadX )[ 3 ],
                                         const variableVector &F, const variableVector &PK2,
                                         const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                         const variableMatrix &DPK2Dgrad_phi,
                                         variableMatrix &DfintDU ){
        /*!
         * Compute the jacobian of the internal force
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * returns 0 if there are no errors, 1 if F is not the right size, 2 if PK2 is not the right size, 3 if DPK2Dgrad_u is
         * not the right size, 4 if DPK2Dphi is not the right size, and 5 if DPK2Dgrad_phi is not the right size.
         *
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. X in the reference configuration.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. X in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The second Piola Kirchoff stress.
         * :param const variableMatrix &DPK2Dgrad_u: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient of the 
         *     macro displacement w.r.t. X in the reference configuration.
         * :param const variableVector &DPK2Dphi: The Jacobian of the second Piola Kirchoff stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient 
         *     of the micro displacement.
         * :param variableMatrix DfintDU: The Jacobian of the internal force w.r.t. the degree of freedom vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 degrees of freedom
        const unsigned int NDOF = 12;

        int errorCode;
        DfintDU = variableMatrix( dim, variableVector( NDOF, 0 ) );
        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < NDOF; j++ ){
                errorCode = compute_internal_force_jacobian( i, j, N, dNdX, eta, detadX, F, PK2,
                                                             DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                                             DfintDU[ i ][ j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_force_jacobian( const unsigned int &i, const unsigned int &j,
                                         const double &N, const double ( &dNdX )[ 3 ], const double &eta, const double ( &detadX )[ 3 ],
                                         const variableVector &F, const variableVector &PK2,
                                         const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                         const variableMatrix &DPK2Dgrad_phi,
                                         variableType &DfintDU_ij ){
        /*!
         * Compute the jacobian of the internal force
         *
         * fint_i = - N_{ , I } PK2_{ I J } F_{ i J }
         *
         * returns 0 if there are no errors, 1 if F is not the right size, 2 if PK2 is not the right size, 3 if DPK2Dgrad_u is
         * not the right size, 4 if DPK2Dphi is not the right size, 5 if DPK2Dgrad_phi is not the right size, 6 if 
         * i is out of range, 7 if j is out of range.
         *
         * :param const unsigned int &i: The row index of the Jacobian ( the internal force vector index )
         * :param const unsigned int &j: The column index of the Jacobian ( the degree of freedom vector index )
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. X in the reference configuration.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. X in the reference configuration.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &PK2: The second Piola Kirchoff stress.
         * :param const variableMatrix &DPK2Dgrad_u: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient of the 
         *     macro displacement w.r.t. X in the reference configuration.
         * :param const variableVector &DPK2Dphi: The Jacobian of the second Piola Kirchoff stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The Jacobian of the second Piola Kirchoff stress w.r.t. the gradient 
         *     of the micro displacement.
         * :param variableType DfintDU_ij: The ij'th component of the internal force jacobian w.r.t. the degree of freedom vector
         *     which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         *     where it's position in this vector is j
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 DOF
        const unsigned int NDOF = 12;

#ifndef SKIP_ERROR_HANDLING
        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( PK2.size() != dim * dim ){
            return 2;
        }

        if ( DPK2Dgrad_u.size() != dim * dim ){
            return 3;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_u.size(); i++ ){
            if ( DPK2Dgrad_u[ i ].size() != dim * dim ){
                return 3;
            }
        }

        if ( DPK2Dphi.size() != dim * dim ){
            return 4;
        }

        for ( unsigned int i = 0; i < DPK2Dphi.size(); i++ ){
            if ( DPK2Dphi[ i ].size() != dim * dim ){
                return 4;
            }
        }

        if ( DPK2Dgrad_phi.size() != dim * dim ){
            return 5;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_phi.size(); i++ ){
            if ( DPK2Dgrad_phi[ i ].size() != dim * dim * dim ){
                return 5;
            }
        }

        if ( i >= dim ){
            return 6;
        }

        if ( j >= NDOF ){
            return 7;
        }
#endif

        DfintDU_ij = 0;

        //Assemble the Jacobians w.r.t. the macro displacement
        if ( j < 3 ){
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        DfintDU_ij -= dNdX[ I ] * DPK2Dgrad_u[ dim * I + J ][ dim * j + K ] * F[ dim * i + J ] * detadX[ K ];
                    }

                    if ( i == j ){
                        DfintDU_ij -= dNdX[ I ] * PK2[ dim * I + J ] * detadX[ J ];
                    }
                }
            }
        }
        else if ( ( j < 12 ) && ( j >= 3 ) ){

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    DfintDU_ij -= dNdX[ I ] * DPK2Dphi[ dim * I + J ][ j - 3 ] * eta * F[ dim * i + J ];

                    for ( unsigned int K = 0; K < dim; K++ ){
                        DfintDU_ij -= dNdX[ I ] * DPK2Dgrad_phi[ dim * I + J ][ dim * ( j - 3 ) + K ] * detadX[ K ] * F[ dim * i + J ];
                    }
                }
            }
        }

        return 0;
    }

    int compute_internal_couple_jacobian( const double &N, const double ( &dNdX )[ 3 ],
                                          const double &eta, const double ( &detadX )[ 3 ],
                                          const variableVector &F, const variableVector &chi,
                                          const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                          const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                          const variableMatrix &DPK2Dgrad_phi,
                                          const variableMatrix &DSIGMADgrad_u, const variableMatrix &DSIGMADphi,
                                          const variableMatrix &DSIGMADgrad_phi,
                                          const variableMatrix &DMDgrad_u, const variableMatrix &DMDphi,
                                          const variableMatrix &DMDgrad_phi,
                                          variableMatrix &DcintDU ){
        /*!
         * Compute the jacobian of the internal couple
         * cint_{ ij } = N F_{ iI } ( PK2_{ JI } - SIGMA_{ JI } ) F_{ jJ } - N_{ ,K } F_{ iI } \chi_{ jJ } M_{ KIJ }
         *
         * Returns 0 if there are no errors, 1 if F has an incorrect size, 2 if chi has an incorrect size,
         * 3 if PK2 has an incorrect size, 4 if SIGMA has an incorrect size, 5 if M has an incorrect size,
         * 6 if DPK2Dgrad_u has an incorrect size, 7 if DPK2Dphi has an incorrect size, 8 if DPK2Dgrad_phi
         * has an incorrect size, 9 if DSIGMADgrad_u has an incorrect size, 10 if DSIGMADphi has an incorrect
         * size, 11 if DSIGMADgrad_phi has an incorrect size, 12 if DMDgrad_u has an incorrect size, 13 if
         * DMDphi has an incorrect size, 14 if DMDgrad_phi has an incorrect size, 15 if i is out of range
         * and 16 if j is out of range.
         *
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference coordinates.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. the reference coordinates.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &chi: The micro deformation tensor.
         * :param const variableVector &PK2: The second Piola Kirchoff stress tensor.
         * :param const variableVector &SIGMA: The symmetric micro stress tensor in the reference configuration.
         * :param const variableVector &M: The higher order stress tensor in the reference configuration.
         * :param const variableMatrix &DPK2Dgrad_u: The derivative of the PK2 stress w.r.t. the gradient of the
         *     macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DPK2Dphi: The derivative of the PK2 stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The derivative of the PK2 stress w.r.t. the gradient of the 
         *     micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADgrad_u: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADphi: The derivative of the symmetric micro stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DSIGMADgrad_phi: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDgrad_u: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDphi: The derivative of the higher order stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DMDgrad_phi: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableMatrix &DcintDU: The Jacobian of the internal couple stress w.r.t. the degree of freedom
         *     vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        const unsigned int dim = 3;

        //Assume 12 degrees of freedom
        const unsigned int NDOF = 12;

        int errorCode;

        //Compute the Jacobian terms
        DcintDU = variableMatrix( dim * dim, variableVector( NDOF, 0 ) );
        for ( unsigned int i = 0; i < dim * dim; i++ ){
            for ( unsigned int j = 0; j < NDOF; j++ ){
                errorCode = compute_internal_couple_jacobian( i, j, N, dNdX, eta, detadX, F, chi, PK2, SIGMA, M,
                                                              DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                                              DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                                              DMDgrad_u, DMDphi, DMDgrad_phi,
                                                              DcintDU[ i ][ j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){
                    return errorCode;
                }
#endif
            }
        }

        return 0;
    }

    int compute_internal_couple_jacobian( const unsigned int &i, const unsigned int &j,
                                          const double &N, const double ( &dNdX )[ 3 ],
                                          const double &eta, const double ( &detadX )[ 3 ],
                                          const variableVector &F, const variableVector &chi,
                                          const variableVector &PK2, const variableVector &SIGMA, const variableVector &M,
                                          const variableMatrix &DPK2Dgrad_u, const variableMatrix &DPK2Dphi,
                                          const variableMatrix &DPK2Dgrad_phi,
                                          const variableMatrix &DSIGMADgrad_u, const variableMatrix &DSIGMADphi,
                                          const variableMatrix &DSIGMADgrad_phi,
                                          const variableMatrix &DMDgrad_u, const variableMatrix &DMDphi,
                                          const variableMatrix &DMDgrad_phi,
                                          variableType &DcintDU_ij ){
        /*!
         * Compute the jacobian of the internal couple
         * cint_{ ij } = N F_{ iI } ( PK2_{ JI } - SIGMA_{ JI } ) F_{ jJ } - N_{ ,K } F_{ iI } \chi_{ jJ } M_{ KIJ }
         *
         * Returns 0 if there are no errors, 1 if F has an incorrect size, 2 if chi has an incorrect size,
         * 3 if PK2 has an incorrect size, 4 if SIGMA has an incorrect size, 5 if M has an incorrect size,
         * 6 if DPK2Dgrad_u has an incorrect size, 7 if DPK2Dphi has an incorrect size, 8 if DPK2Dgrad_phi
         * has an incorrect size, 9 if DSIGMADgrad_u has an incorrect size, 10 if DSIGMADphi has an incorrect
         * size, 11 if DSIGMADgrad_phi has an incorrect size, 12 if DMDgrad_u has an incorrect size, 13 if
         * DMDphi has an incorrect size, 14 if DMDgrad_phi has an incorrect size, 15 if i is out of range
         * and 16 if j is out of range.
         *
         * :param const unsigned int &i: The row index of the Jacobian ( the internal couple vector index )
         * :param const unsigned int &j: The column index of the Jacobian ( the degree of freedom vector index )
         * :param const double &N: The shape function value
         * :param const double ( &dNdX )[ 3 ]: The gradient of the shape function w.r.t. the reference coordinates.
         * :param const double &eta: The interpolation function value
         * :param const double ( &detadX )[ 3 ]: The gradient of the interpolation function w.r.t. the reference coordinates.
         * :param const variableVector &F: The deformation gradient.
         * :param const variableVector &chi: The micro deformation tensor.
         * :param const variableVector &PK2: The second Piola Kirchoff stress tensor.
         * :param const variableVector &SIGMA: The symmetric micro stress tensor in the reference configuration.
         * :param const variableVector &M: The higher order stress tensor in the reference configuration.
         * :param const variableMatrix &DPK2Dgrad_u: The derivative of the PK2 stress w.r.t. the gradient of the
         *     macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DPK2Dphi: The derivative of the PK2 stress w.r.t. the micro displacement.
         * :param const variableMatrix &DPK2Dgrad_phi: The derivative of the PK2 stress w.r.t. the gradient of the 
         *     micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADgrad_u: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DSIGMADphi: The derivative of the symmetric micro stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DSIGMADgrad_phi: The derivative of the symmetric micro stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDgrad_u: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the macro displacement w.r.t. the reference configuration.
         * :param const variableMatrix &DMDphi: The derivative of the higher order stress in the reference 
         *     configuration w.r.t. the micro displacement.
         * :param const variableMatrix &DMDgrad_phi: The derivative of the higher order stress in the reference
         *     configuration w.r.t. the gradient of the micro displacement w.r.t. the reference configuration.
         * :param variableMatrix &DcintDU: The Jacobian of the internal couple stress w.r.t. the degree of freedom
         *     vector which is organized:
         *     [ u1, u2, u3, phi_11, phi_12, phi_13, phi_21, phi_22, phi_23, phi_31, phi_32, phi_33 ]
         */

        //Assume 3D
        unsigned int dim = 3;

        //Assume 12 DOF
        unsigned int NDOF = 12;

#ifndef SKIP_ERROR_HANDLING

        if ( F.size() != dim * dim ){
            return 1;
        }

        if ( chi.size() != dim * dim ){
            return 2;
        }

        if ( PK2.size() != dim * dim ){
            return 3;
        }

        if ( SIGMA.size() != dim * dim ){
            return 4;
        }

        if ( M.size() != dim * dim * dim ){
            return 5;
        }

        if ( DPK2Dgrad_u.size() != dim * dim ){
            return 6;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_u.size(); i++ ){
            if ( DPK2Dgrad_u[ i ].size() != dim * dim ){
                return 6;
            }
        }

        if ( DPK2Dphi.size() != dim * dim ){
            return 7;
        }

        for ( unsigned int i = 0; i < DPK2Dphi.size(); i++ ){
            if ( DPK2Dphi[ i ].size() != dim * dim ){
                return 7;
            }
        }

        if ( DPK2Dgrad_phi.size() != dim * dim ){
            return 8;
        }

        for ( unsigned int i = 0; i < DPK2Dgrad_phi.size(); i++ ){
            if ( DPK2Dgrad_phi[ i ].size() != dim * dim * dim ){
                return 8;
            }
        }

        if ( DSIGMADgrad_u.size() != dim * dim ){
            return 9;
        }

        for ( unsigned int i = 0; i < DSIGMADgrad_u.size(); i++ ){
            if ( DSIGMADgrad_u[ i ].size() != dim * dim ){
                return 9;
            }
        }

        if ( DSIGMADphi.size() != dim * dim ){
            return 10;
        }

        for ( unsigned int i = 0; i < DSIGMADphi.size(); i++ ){
            if ( DSIGMADphi[ i ].size() != dim * dim ){
                return 10;
            }
        }

        if ( DSIGMADgrad_phi.size() != dim * dim ){
            return 11;
        }

        for ( unsigned int i = 0; i < DSIGMADgrad_phi.size(); i++ ){
            if ( DSIGMADgrad_phi[ i ].size() != dim * dim * dim ){
                return 11;
            }
        }

        if ( DMDgrad_u.size() != dim * dim * dim ){
            return 12;
        }

        for ( unsigned int i = 0; i < DMDgrad_u.size(); i++ ){
            if ( DMDgrad_u[ i ].size() != dim * dim ){
                return 12;
            }
        }

        if ( DMDphi.size() != dim * dim * dim ){
            return 13;
        }

        for ( unsigned int i = 0; i < DMDphi.size(); i++ ){
            if ( DMDphi[ i ].size() != dim * dim ){
                return 13;
            }
        }

        if ( DMDgrad_phi.size() != dim * dim * dim ){
            return 14;
        }

        for ( unsigned int i = 0; i < DMDgrad_phi.size(); i++ ){
            if ( DMDgrad_phi[ i ].size() != dim * dim * dim ){
                return 14;
            }
        }

        if ( i > dim * dim ){
            return 15;
        }

        if ( j > NDOF ){
            return 16;
        }

#endif

        DcintDU_ij = 0;
        if ( j < 3 ){

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( DPK2Dgrad_u[ dim * J + I ][ dim * j + K ]
                                                                     - DSIGMADgrad_u[ dim * J + I ][ dim * j + K ]
                                                                     ) * detadX[ K ] * F[ dim * ( i % 3 ) + J ];

                        for ( unsigned int L = 0; L < dim; L++ ){
                            DcintDU_ij -= dNdX[ K ] * F[ dim * ( i / 3 ) + I ] * chi[ dim * ( i % 3 ) + J ]
                                        * DMDgrad_u[ dim * dim * K + dim * I + J ][ dim * j + L ] * detadX[ L ];
                        }
                        
                        if ( ( i / 3) == j ){
                            DcintDU_ij -= dNdX[ K ] * detadX[ I ] * chi[ dim * ( i % 3 ) + J ] * M[ dim * dim * K + dim * I + J ];
                        }
                    }

                    if ( ( i / 3 ) == j ){
                        DcintDU_ij += N * detadX[ I ] * ( PK2[ dim * J + I ] - SIGMA[ dim * J + I ] ) * F[ dim * ( i % 3 ) + J ];
                    }

                    if ( ( i % 3 ) == j ){
                        DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( PK2[ dim * J + I ] - SIGMA[ dim * J + I ] ) * detadX[ J ];
                    }
                }
            }
        }
        else if ( ( j >= 3 ) && ( j < 12 ) ){

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( DPK2Dphi[ dim * J + I ][ j - 3 ]
                                                                 - DSIGMADphi[ dim * J + I ][ j - 3 ]
                                                                 ) * eta * F[ dim * ( i % 3 ) + J ];
                    for ( unsigned int K = 0; K < dim; K++ ){
                        DcintDU_ij += N * F[ dim * ( i / 3 ) + I ] * ( DPK2Dgrad_phi[ dim * J + I ][ dim * ( j - 3 ) + K ]
                                                                     - DSIGMADgrad_phi[ dim * J + I ][ dim * ( j - 3 ) + K ]
                                                                     ) * detadX[ K ] * F[ dim * ( i % 3 ) + J ];

                        DcintDU_ij -= dNdX[ K ] * F[ dim * ( i / 3 ) + I ] * chi[ dim * ( i % 3 ) + J ]
                                    * DMDphi[ dim * dim * K + dim * I + J ][ j - 3 ] * eta;

                        for ( unsigned int L = 0; L < dim; L++ ){
                            DcintDU_ij -= dNdX[ K ] * F[ dim * ( i / 3 ) + I ] * chi[ dim * ( i % 3 ) + J ] * DMDgrad_phi[ dim * dim * K + dim * I + J ][ dim * ( j - 3 ) + L ] * detadX[ L ];
                        }

                        if ( ( dim * ( i % 3 ) + J ) == ( j - 3 ) ){
                            DcintDU_ij -= dNdX[ K ] * F[ dim * ( i / 3 ) + I ] * eta * M[ dim * dim * K + dim * I + J ];
                        }
                    }
                }
            }

        }

        return 0;
    }

    int compute_inertia_couple_jacobian( const double &N, const double &eta, const double &density, const double ( &chi )[ 9 ],
                                         const double ( &D2ChiDt2 )[ 9 ], const variableMatrix &D3ChiDt2DChi,
                                         const double ( &referenceInertia )[ 9 ], variableMatrix &DcinertiaDU ){
        /*!
         * Compute the jacobian of the inertia couple term
         *
         * integrand_kl = -N \rho_0 I_{KL } \ddot{ \chi }_{kK} \chi_{lL}
         *
         * returns:
         * 0 if no errors
         * 1 if D3ChiDt2DChi doesn't have the proper number of rows
         * 2 if D3ChiDt2DChi doesn't have the proper number of columns
         *
         * :param const double &N: The shapefunction value
         * :param const double &eta: The interpolation function value
         * :param const double &density: The density in the reference configuration
         * :param const double ( &chi )[ 9 ]: The micro-deformation
         * :param const double ( &D2ChiDt2 )[ 9 ]: The second temporal derivative of the micro-deformation
         * :param const variableMatrix &D3ChiDt2DChi: The derivative of the second temporal derivative of 
         *     the micro-deformation w.r.t. the micro-deformation
         * :param const double ( &referenceInertia )[ 9 ]: The mass moment of inertia in the reference configuration
         * :param variableMatrix &DcinertiaDU: The Jacobian of the inertia couple term w.r.t. the degree of freedom vector
         */

        //Assume 3d
        const unsigned int dim = 3;

        //Error handling
#ifndef SKIP_ERROR_HANDLING
        if ( D3ChiDt2DChi.size( ) != dim * dim ){

            //Return the error code
            return 1;

        }

        for ( unsigned int i = 0; i < D3ChiDt2DChi.size( ); i++ ){

            if ( D3ChiDt2DChi[ i ].size( ) != dim * dim ){

                return 2;

            }

        }
#endif

        DcinertiaDU = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

        variableVector D3ChiDt2DChi_j( dim * dim, 0 );

        int errorCode;

        for ( unsigned int i = 0; i < dim * dim; i++ ){

            for ( unsigned int j = 0; j < dim * dim; j++ ){

                for ( unsigned int _i = 0; _i < dim * dim; _i++ ){

                    D3ChiDt2DChi_j[ _i ] = D3ChiDt2DChi[ _i ][ j ];

                }

                errorCode = compute_inertia_couple_jacobian( i, j, N, eta, density, chi, D2ChiDt2,
                                                             D3ChiDt2DChi_j, referenceInertia, DcinertiaDU[ i ][ j ] );

#ifndef SKIP_ERROR_HANDLING
                if ( errorCode != 0 ){

                    return errorCode;

                }
#endif

            }

        }

        return 0;
    }

    int compute_inertia_couple_jacobian( const unsigned int &i, const unsigned int &j,
                                         const double &N, const double &eta, const double &density, const double ( &chi )[ 9 ],
                                         const double ( &D2ChiDt2 )[ 9 ], const variableVector &D3ChiDt2DChi_j,
                                         const double ( &referenceInertia )[ 9 ], double &DcinertiaDU_ij ){
        /*
         * Compute the jacobian of the inertia couple term
         *
         * integrand_kl = -N \rho_0 I_{KL } \ddot{ \chi }_{kK} \chi_{lL}
         *
         * returns:
         * 0 if no errors
         * 2 if D3ChiDt2DChi doesn't have the proper number of columns
         *
         * :param const unsigned int &i: The row index of the Jacobian required the Jacobian rows are stored as
         *     [ 11, 12, 13, 21, 22, 23, 31, 32, 33 ] where i refers to the index of the preceeding array
         * :param const unsigned int &j: The column index of the Jacobian required which is stored as
         *     [ 11, 12, 13, 21, 22, 23, 31, 32, 33 ] where j refers to the index of the preceeding array
         * :param const double &N: The shapefunction value
         * :param const double &eta: The interpolation function value
         * :param const double &density: The density in the reference configuration
         * :param const double ( &chi )[ 9 ]: The micro-deformation
         * :param const double ( &D2ChiDt2 )[ 9 ]: The second temporal derivative of the micro-deformation
         * :param const variableVector &D3ChiDt2DChi_j: The derivative of the second temporal derivative of 
         *     the micro-deformations w.r.t. the micro-deformation ( i.e. a column of the full Jacobian )
         * :param const double ( &referenceInertia )[ 9 ]: The mass moment of inertia in the reference configuration
         * :param double &DcinertiaDU_ij: The Jacobian of the inertia couple term w.r.t. the degree of freedom vector
         *     for the indices i and j
         */

        //Assume 3D
        const unsigned int dim = 3;

        const unsigned int k = i / dim;
        const unsigned int l = i % dim;

        const unsigned int m = j / dim;
        const unsigned int M = j % dim;

        //Error handling
#ifndef SKIP_ERROR_HANDLING
        if ( D3ChiDt2DChi_j.size( ) != dim * dim ){

            return 2;

        }
#endif

        DcinertiaDU_ij = 0;
        for ( unsigned int K = 0; K < dim; K++ ){

            for ( unsigned int L = 0; L < dim; L++ ){

                DcinertiaDU_ij -= N * eta * density * referenceInertia[ dim * K + L ]
                                * D3ChiDt2DChi_j[ dim * k + K ] * chi[ dim * l + L ];

                if ( ( l == m ) && ( L == M ) ){

                    DcinertiaDU_ij -= N * eta * density * referenceInertia[ dim * K + L ] * D2ChiDt2[ dim * k + K ];

                }

            }

        }

        return 0;
    }

    int compute_inertia_force_jacobian( const double &N, const double &eta, const double &density, const double ( &a )[ 3 ],
                                        const variableMatrix &DaDu, variableMatrix &DfinertiaDU ){
        /*!
         * Compute the inertia force jacobian
         *
         * finertia_i = -N * \rho * a_i
         *
         * return:
         * 0 if no errors
         * 1 if DaDu is not consistent with the dimension
         *
         * :param const double &N: The shape function value
         * :param const double &eta: The interpolation function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const variableMatrix &DaDu: The jacobian of the acceleration w.r.t. the degree of freedom vector
         *     [ u_1, u_2, u_3 ] is the ordering of the degree of freedom vector
         * :param variableMatrix &DfinertiaDU: The jacobian of the inertia force
         */

        //Assume 3D
        const unsigned int dim = 3;

        if ( dim != DaDu.size( ) ){

            return 1;

        }

        int errorCode = 0;

        DfinertiaDU = variableMatrix( dim, variableVector( dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int j = 0; j < dim; j++ ){

                errorCode = compute_inertia_force_jacobian( i, j, N, eta, density, a, DaDu[ i ], DfinertiaDU[ i ][ j ] );

                if ( errorCode != 0 ){

                    return errorCode;

                }

            }

        }

        return 0;
    }

    int compute_inertia_force_jacobian( const unsigned int &i, const unsigned int &j,
                                        const double &N, const double &eta, const double &density, const double ( &a )[ 3 ],
                                        const variableVector &DaDu_i, variableType &DfinertiaDU_ij ){
        /*!
         * Compute the inertia force jacobian for the ith component in the reference configuration
         *
         * finertia_i = -N * \rho * a_i
         *
         * \frac{ \partial f^{inertia}_i }{ \partial U_j } = -N * \rho * \frac{\partial a_i }{ \partial U_j }
         *
         * return:
         * 0 if no errors
         * 2 if j is larger than the ith row of DaDu
         *
         * :param const unsigned int &i: The component to compute the inertia force on
         * :param const double &N: The shape function value
         * :param const double &density: The density in the current configuration
         * :param const double ( &a )[ 3 ]: The acceleration in the current configuration
         * :param const variableVector &DaDu_i: The ith row of the Jacobian of the acceleration w.r.t. the deformation.
         * :param const double finertia_i: The inertia force in the current configuration in direction i
         */

        if ( j >= DaDu_i.size( ) ){

            return 2;

        }

        DfinertiaDU_ij = -N * density * eta * DaDu_i[ j ];

        return 0;

    }

}
