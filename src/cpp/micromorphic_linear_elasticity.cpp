/*
 * micromorphic_linear_elasticity.cpp
 *
 * An implimentation of linear elasticity in the micromorphic context.
 *
 * Based around a quadratic form of the Helmholtz free energy:
 * \rho \psi = \frac{1}{2} E_{IJ} A_{IJKL} E_{KL} + \frac{1}{2} \mathcal{E}_{IJ} B_{IJKL} \mathcal{E}_{KL} 
 *           + \frac{1}{2} \Gamma_{IJK} C_{IJKLMN} \Gamma_{LMN} + E_{IJ} D_{IJKL} \mathcal{E}_{KL}
 */

#include<micromorphic_linear_elasticity.h>

namespace micromorphicLinearElasticity{

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
         * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         */

        //Compute the required deformation measures
        variableVector RCG, Psi, Gamma;
        errorOut error = computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                     RCG, Psi, Gamma );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReference",
                                             "Error in the computation of the deformation measures" );
            result->addNext( error );
            return result;
        }

        error = linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                          PK2Stress, referenceMicroStress,
                                                          referenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReference",
                                             "Error in the computation of the reference stresses" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation,
                                                       const variableVector &Psi, const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
         * :param const variableVector &Psi: The micro-deformation measure Psi
         * :param const variableVector &Gamma: The higher order deformation measure Gamma
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
         * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the strain measures
        variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
        variableVector microStrain   = Psi - eye;

        //Compute the higher order stress
        errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                             "Error in computation of higher-order stress" );
            result->addNext( error );
            return result;
        }

        //Compute the first common term for the PK2 and symmetric micro-stress
        variableVector term1;
        error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                             "Error in computation of term 1" );
            result->addNext( error );
            return result;
        }

        //Compute the second common term for the PK2 and symmetric micro-stress
        variableVector invRCGPsi;
        error = computeInvRCGPsi( invRCG, Psi, invRCGPsi );
        
        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                             "Error in computation of invRCG Psi product" );
            result->addNext( error );
            return result;
        }

        variableVector term2;
        error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2 );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMesures",
                                             "Error in computation of term 2" );
            result->addNext( error );
            return result;
        }

        //Compute the third common term for the PK2 and symmetric micro-stress
        variableVector invRCGGamma;
        error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                             "Error in computation of invRCG Gamma product" );
            result->addNext(error);
            return result;
        }

        variableVector term3;
        error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3 );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                             "Error in computation of term 3" );
            result->addNext( error );
            return result;
        }

        //Construct the PK2 and reference symmetric stresses
        PK2Stress            = term1 + term2 + term3;

        variableVector symmTerm2Term3;
        error = constitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3 );
        referenceMicroStress = term1 + 2 * symmTerm2Term3;

        return NULL;
    }

    errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                        const variableVector &gradientMicroDeformation,
                                        const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                        const parameterVector &D,
                                        variableVector &PK2Stress, variableVector &referenceMicroStress,
                                        variableVector &referenceHigherOrderStress,
                                        variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdXi, variableMatrix &dPK2StressdGradXi,
                                        variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdXi,
                                        variableMatrix &dReferenceMicroStressdGradXi, variableMatrix &dMdF, variableMatrix &dMdGradXi ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * Also computes the Jacobians
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
         * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         * :param variableMatrix &dPK2StressdF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
         * :param variableMatrix &dPK2StressdXi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
         * :param variableMatrix &dPK2StressdGradXi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
         * :param variableMatrix &dReferenceMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient.
         * :param variableMatrix &dReferenceMicroStressdXi: The Jacobian of the Micro stress w.r.t. the micro deformation.
         * :param variableMatrix &dReferenceStressdGradXi: The Jacobian of the Micro stress w.r.t. the gradient of the micro deformation.
         * :param variableMatrix &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
         * :param variableMatrix &dMdGradXi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
         */

        //Assume 3d
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        //Compute the required deformation measures
        variableVector RCG, Psi, Gamma;
        variableMatrix dRCGdF, dPsidF, dPsidXi, dGammadF, dGammadGradXi;
        errorOut error = computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                     RCG, Psi, Gamma, dRCGdF, dPsidF, dPsidXi, dGammadF, dGammadGradXi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReference (jacobian)",
                                             "Error in the computation of the deformation measures" );
            result->addNext( error );
            return result;
        }

        variableVector invRCG = vectorTools::inverse( RCG, dim, dim );

        //Compute the strain measures
        variableVector greenLagrangeStrain = 0.5 * ( RCG - eye );
        variableVector microStrain   = Psi - eye;

        variableMatrix dGreenLagrangeStraindF = 0.5 * dRCGdF;
        variableMatrix dMicroStraindF  = dPsidF;
        variableMatrix dMicroStraindXi = dPsidXi;

        //Compute the higher order stress
        variableMatrix dMdGamma;
        error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress, dMdGamma );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity (jacobian)",
                                             "Error in computation of higher-order stress" );
            result->addNext( error );
            return result;
        }

        dMdF = vectorTools::dot( dMdGamma, dGammadF );
        dMdGradXi = vectorTools::dot( dMdGamma, dGammadGradXi );

        //Compute the first common term for the PK2 and symmetric micro-stress
        variableVector term1;

        variableMatrix dTerm1dGreenLagrangeStrain, dTerm1dMicroStrain;
        error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1,
                                           dTerm1dGreenLagrangeStrain, dTerm1dMicroStrain );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity (jacobian)",
                                             "Error in computation of term 1" );
            result->addNext( error );
            return result;
        }

        //Assemble term1 jacobians w.r.t. F and Xi
        variableMatrix dTerm1dF = vectorTools::dot( dTerm1dGreenLagrangeStrain, dGreenLagrangeStraindF )
                                + vectorTools::dot( dTerm1dMicroStrain, dMicroStraindF );

        variableMatrix dTerm1dXi = vectorTools::dot( dTerm1dMicroStrain, dMicroStraindXi );

        //Compute the second common term for the PK2 and symmetric micro-stress
        variableVector invRCGPsi;
        variableMatrix dInvRCGPsidRCG, dInvRCGPsidPsi;

        error = computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity (jacobian)",
                                             "Error in computation of invRCG Psi product" );
            result->addNext( error );
            return result;
        }

        //Assemble InvCPsi jacobians w.r.t. F and Xi
        variableMatrix dInvRCGPsidF = vectorTools::dot( dInvRCGPsidRCG, dRCGdF )
                                    + vectorTools::dot( dInvRCGPsidPsi, dPsidF );
        variableMatrix dInvRCGPsidXi = vectorTools::dot( dInvRCGPsidPsi, dPsidXi );

        variableVector term2;
        variableMatrix dTerm2dGreenLagrangeStrain, dTerm2dMicroStrain, dTerm2dInvRCGPsi;
        error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                           dTerm2dGreenLagrangeStrain, dTerm2dMicroStrain, dTerm2dInvRCGPsi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in computation of term 2" );
            result->addNext( error );
            return result;
        }
        
        //Assemble term2 jacobians w.r.t. F and Xi
        variableMatrix dTerm2dF = vectorTools::dot( dTerm2dGreenLagrangeStrain, dGreenLagrangeStraindF )
                                + vectorTools::dot( dTerm2dMicroStrain, dMicroStraindF )
                                + vectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidF );

        variableMatrix dTerm2dXi = vectorTools::dot( dTerm2dMicroStrain, dMicroStraindXi )
                                 + vectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidXi );

        //Compute the third common term for the PK2 and symmetric micro-stress
        variableVector invRCGGamma;
        variableMatrix dInvRCGGammadRCG, dInvRCGGammadGamma;

        error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma );

        if ( error ){
            errorOut result = new errorNode( "linear elasticity (jacobian)", "Error in computation of invRCG Gamma product" );
            result->addNext( error );
            return result;
        }

        variableMatrix dInvCGammadF = vectorTools::dot( dInvRCGGammadRCG, dRCGdF )
                                    + vectorTools::dot( dInvRCGGammadGamma, dGammadF );
        variableMatrix dInvCGammadGradXi = vectorTools::dot( dInvRCGGammadGamma, dGammadGradXi );

        variableVector term3;
        variableMatrix dTerm3dInvCGamma, dTerm3dM;
        error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvCGamma, dTerm3dM );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in computation of term 3" );
            result->addNext( error );
            return result;
        }

        variableMatrix dTerm3dF = vectorTools::dot( dTerm3dInvCGamma, dInvCGammadF )
                                + vectorTools::dot( dTerm3dM, dMdF );

        variableMatrix dTerm3dGradXi = vectorTools::dot( dTerm3dInvCGamma, dInvCGammadGradXi )
                                     + vectorTools::dot( dTerm3dM, dMdGradXi );

        //Construct the PK2 and reference symmetric stresses
        PK2Stress            = term1 + term2 + term3;

        dPK2StressdF = dTerm1dF + dTerm2dF + dTerm3dF;
        dPK2StressdXi = dTerm1dXi + dTerm2dXi;
        dPK2StressdGradXi = dTerm3dGradXi;

        variableVector symmTerm2Term3;
        variableMatrix dSymmTerm2Term3dTerm2Term3;
        error = constitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 );
        referenceMicroStress = term1 + 2 * symmTerm2Term3;

        dReferenceMicroStressdF = dTerm1dF + 2 * ( vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dF )
                                                 + vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dF ) );

        dReferenceMicroStressdXi = dTerm1dXi + 2 * vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dXi );
        dReferenceMicroStressdGradXi = 2 * vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dGradXi );

        return NULL;
    }

    errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                       const variableVector &Gamma,
                                                       const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                       const parameterVector &D,
                                                       variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                       variableVector &referenceHigherOrderStress,
                                                       variableMatrix &dPK2StressdRCG, variableMatrix &dPK2StressdPsi,
                                                       variableMatrix &dPK2StressdGamma,
                                                       variableMatrix &dReferenceMicroStressdRCG,
                                                       variableMatrix &dReferenceMicroStressdPsi,
                                                       variableMatrix &dReferenceMicroStressdGamma,
                                                       variableMatrix &dMdGamma ){
        /*!
         * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
         * of a quadratic decomposition of the energy.
         *
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
         * :param const variableVector &Psi: The micro-deformation measure Psi
         * :param const variableVector &Gamma: The higher order deformation measure Gamma
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
         * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
         *     reference configuration.
         * :param variableMatrix &dPK2StressdRCG: The Jacobian of the PK2 stress w.r.t. the right Cauchy-Green
         *     deformation metric.
         * :param variableMatrix &dPK2StressdPsi: The Jacobian of the PK2 stress w.r.t. the micro deformation 
         *     metric.
         * :param variableMatrix &dPK2StressdGamma: The Jacobian of the PK2 stress w.r.t. the higher order 
         *     deformation measure.
         * :param variableMatrix &dReferenceMicroStressdRCG: The Jacobian of the reference micro stress w.r.t. the 
         *     right Cacuhy-Green deformation metric.
         * :param variableMatrix &dReferenceMicroStressdPsi: The Jacobian of the reference micro stress w.r.t. the 
         *     micro deformation measure.
         * :param variableMatrix &dReferenceMicroStrssdGamma: The Jacobian of the reference micro stress w.r.t. the 
         *     higher order deformation measure.
         * :param variableMatrix &dMdGamma: The Jacobian of the reference higher order stress w.r.t. 
         *     the higher order deformation measure.
         */

        //Assume 3d
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the strain measures
        variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
        variableVector microStrain   = Psi - eye;

        //Compute the higher order stress
        errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress, dMdGamma );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityRefereceDerivedMetrics (jacobian)",
                                             "Error in computation of higher-order stress" );
            result->addNext( error );
            return result;
        }

        //Compute the first common term for the PK2 and symmetric micro-stress
        variableVector term1;

        variableMatrix dTerm1dRCG, dTerm1dPsi;
        error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1,
                                           dTerm1dRCG, dTerm1dPsi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityDerivedMetrics (jacobian)",
                                             "Error in computation of term 1" );
            result->addNext( error );
            return result;
        }

        //Assemble term1 jacobians w.r.t. F and Xi
        dTerm1dRCG *= 0.5;

        //Compute the second common term for the PK2 and symmetric micro-stress
        variableVector invRCGPsi;
        variableMatrix dInvRCGPsidRCG, dInvRCGPsidPsi;

        error = computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                             "Error in computation of invRCG Psi product" );
            result->addNext( error );
            return result;
        }

        variableVector term2;
        variableMatrix dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi;
        error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                           dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in computation of term 2" );
            result->addNext( error );
            return result;
        }

        dTerm2dRCG *= 0.5;
        dTerm2dRCG += vectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidRCG );

        dTerm2dPsi += vectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidPsi );

        //Compute the third common term for the PK2 and symmetric micro-stress
        variableVector invRCGGamma;
        variableMatrix dInvRCGGammadRCG, dInvRCGGammadGamma;

        error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma );

        if ( error ){
            errorOut result = new errorNode( "linear elasticity (jacobian)", "Error in computation of invRCG Gamma product" );
            result->addNext( error );
            return result;
        }

        variableVector term3;
        variableMatrix dTerm3dInvRCGGamma, dTerm3dM;
        error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvRCGGamma, dTerm3dM );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in computation of term 3" );
            result->addNext( error );
            return result;
        }

        variableMatrix dTerm3dRCG = vectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadRCG );
        variableMatrix dTerm3dGamma = vectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadGamma )
                                    + vectorTools::dot( dTerm3dM, dMdGamma );

        //Construct the PK2 and reference symmetric stresses
        PK2Stress            = term1 + term2 + term3;

        dPK2StressdRCG    = dTerm1dRCG + dTerm2dRCG + dTerm3dRCG;
        dPK2StressdPsi    = dTerm1dPsi + dTerm2dPsi;
        dPK2StressdGamma  = dTerm3dGamma;

        variableVector symmTerm2Term3;
        variableMatrix dSymmTerm2Term3dTerm2Term3;
        error = constitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 );
        referenceMicroStress = term1 + 2 * symmTerm2Term3;

        dReferenceMicroStressdRCG = dTerm1dRCG + 2 * ( vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dRCG )
                                                     + vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dRCG ) );

        dReferenceMicroStressdPsi = dTerm1dPsi + 2 * vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dPsi );
        dReferenceMicroStressdGamma = 2 * vectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dGamma );

        return NULL;
    }

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma ){
        /*!
         * Compute the deformation measures
         * C_{IJ} = F_{iI} F_{iJ}
         * Psi_{IJ} = F_{iI} \Xi_{iJ}
         * \Gamma_{IJK} = F_{iI} \Xi_{iJ, K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor
         * :param variableVector &Psi: The micro-deformation measure
         * :param variableVector &Gamma: The gradient micro-deformation measure
         */

        errorOut error = constitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen );

        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures",
                                             "Error in the computation of the right Cauchy-Green Deformation measure" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( deformationGradient, microDeformation, Psi );
        
        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures",
                                             "Error in the computation of Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma );

        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures",
                                             "Error in the computation of Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;

    }

    errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &gradientMicroDeformation,
                                         variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                         variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidXi,
                                         variableMatrix &dGammadF, variableMatrix &dGammadGradXi ){
        /*!
         * Compute the deformation measures
         * C_{IJ} = F_{iI} F_{iJ}
         * Psi_{IJ} = F_{iI} \Xi_{iJ}
         * \Gamma_{IJK} = F_{iI} \Xi_{iJ, K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor
         * :param variableVector &Psi: The micro-deformation measure
         * :param variableVector &Gamma: The gradient micro-deformation measure
         * :param variableMatrix &dCdF: The gradient of the right Cauchy green deformation tensor w.r.t. 
         *     the deformation gradient.
         * :param variableMatrix &dPsidF: The gradient of Psi w.r.t. the deformation gradient.
         * :param variableMatrix &dPsidXi: The gradient of Psi w.r.t. the microDeformation.
         * :param variableMatrix &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * :param variableMatrix &dGammadGradXi: The gradient of Gamma w.r.t. the spatial gradient of Xi
         */

        errorOut error = constitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen, dCdF );

        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                             "Error in the computation of the right Cauchy-Green Deformation measure" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computePsi( deformationGradient, microDeformation, Psi, dPsidF, dPsidXi );
        
        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                             "Error in the computation of Psi" );
            result->addNext( error );
            return result;
        }

        error = micromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma, dGammadF, dGammadGradXi );

        if ( error ){
            errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                             "Error in the computation of Gamma" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain, 
                                        const parameterVector &A, const parameterVector &D, variableVector &term1 ){
        /*!
         * Compute the first term for the linear elastic model
         * term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL}
         *
         * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain.
         * :param const variableVector &microStrain: The micro-strain
         * :param const parameterVector &A: The A stiffness matrix
         * :param const parameterVEctor &D: The D stiffness matrix
         * :param variableVector &term1: The first term.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( greenLagrangeStrain.size() != dim * dim ){
            return new errorNode( "computeLinearElasticTerm1",
                                  "The green lagrange strain must have a length of 9" );
        }

        if ( microStrain.size() != dim * dim ){
            return new errorNode( "computeLinearElasticTerm1",
                                  "The micro-strain must have a length of 9" );
        }

        if ( A.size() != dim * dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm1",
                                  "A must have a size of 3**4" );
        }

        if ( D.size() != dim * dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm1",
                                  "D must have a size of 3**4" );
        }

        //Compute the first common term for the PK2 and symmetric micro-stress
        term1 = variableVector( dim * dim, 0 );
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        term1[ dim * I + J ] += A[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * greenLagrangeStrain[ dim * K + L ]
                                              + D[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * microStrain[ dim * K + L ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain, 
                                        const parameterVector &A, const parameterVector &D, variableVector &term1,
                                        variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain ){
        /*!
         * Compute the first term for the linear elastic model
         * term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL}
         *
         * Also return the Jacobian
         * \frac{\partial term^1_{IJ} }{ E_{MN} } = A_{IJMN}
         * \frac{\partial term^1_{IJ} }{ \mathcal{E}_{MN} } = D_{IJMN}
         *
         * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain.
         * :param const variableVector &microStrain: The micro-strain
         * :param const parameterVector &A: The A stiffness matrix
         * :param const parameterVEctor &D: The D stiffness matrix
         * :param variableVector &term1: The first term.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 );

        if ( error ){
            errorOut result = new errorNode( "computeLinearElasticTerm1 (jacobian)",
                                             "Error in computation of linear elastic term1" );
            result->addNext( error );
            return result;
        }

        //Compute the first common term for the PK2 and symmetric micro-stress
        dTerm1dGreenLagrangeStrain = variableMatrix( term1.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
        dTerm1dMicroStrain = variableMatrix( term1.size(), variableVector( microStrain.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        dTerm1dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] = A[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                        dTerm1dMicroStrain[ dim * I + J ][ dim * M + N ] = D[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2 ){
        /*!
         * Compute the second term from the linear elastic constitutive model
         *
         * term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}
         *
         * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)
         * :param const variableVector &microStrain: The micro-strain \mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}
         * :param const variableVector &invCPsi: The product C_{JR}^{-1} \Psi_{RQ}
         * :param const variableVector &B: The B stiffness matrix
         * :param const variableVector &D: The D stiffness matrix
         * :param variableVector &term2: The second term.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( greenLagrangeStrain.size() != dim * dim ){
            return new errorNode( "computeLinearElasticTerm2",
                                  "The green lagrange strain must have a length of 9" );
        }

        if ( microStrain.size() != dim * dim ){
            return new errorNode( "computeLinearElasticTerm2",
                                  "The micro-strain must have a length of 9" );
        }

        if ( invCPsi.size() != dim * dim ){
            return new errorNode( "computeLinearElasticTerm2",
                                  "invCPsi must have a size of 9" );
        }

        if ( B.size() != dim * dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm2",
                                  "B must have a size of 3**4" );
        }

        if ( D.size() != dim * dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm2",
                                  "D must have a size of 3**4" );
        }

        term2 = variableVector( greenLagrangeStrain.size(), 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int Q = 0; Q < dim; Q++ ){
                            term2[ dim * I + J] += ( B[ dim * dim * dim * I + dim * dim * Q + dim * K + L ] * microStrain[ dim * K + L ]
                                                 + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + Q ] )
                                                 * invCPsi[ dim * J + Q ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                        const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                        variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                        variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi ){
        /*!
         * Compute the second term from the linear elastic constitutive model
         *
         * term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}
         *
         * Also return the Jacobians
         * \frac{ \partial term^2_{IJ} }{ \partial E_{MN} } = D_{MNIK} C_{JR}^{-1} \Psi_{RK}
         * \frac{ \partial term^2_{IJ} }{ \partial \mathcal{E}_{MN} } = B_{IKMN} C_{JR}^{-1} \Psi_{RK}
         * \frac{ \partial term^2_{IJ} }{ \partial C_{MO}^{-1} \Psi_{ON} } = \left( B_{INKL} \mathcal{E}_{KL} + E_{KL} D_{KLIN} \right) \delta_{JM} 
         *
         * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)
         * :param const variableVector &microStrain: The micro-strain \mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}
         * :param const variableVector &invCPsi: The product C_{JR}^{-1} \Psi_{RQ}
         * :param const variableVector &B: The B stiffness matrix
         * :param const variableVector &D: The D stiffness matrix
         * :param variableVector &term2: The second term.
         * :param variableMatrix &dTerm2dGreenLagrangeStrain: The jacobian of term 2 w.r.t. the Green-Lagrange strain.
         * :param variableMatrix &dTerm2dMicroStrain: The jacobian of term 2 w.r.t. the microStrain.
         * :param variableMatrix &dTerm2dInvCPsi: The jacobian of term 2 w.r.t. C_{JR}^{-1} \Psi_{RQ}
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invCPsi, B, D, term2 );

        if ( error ){
            errorOut result = new errorNode("computeLinearElasticTerm2 (jacobian)",
                                            "Error in computation of term 2" );
            result->addNext( error );
            return result;
        }

        //Compute the Jacobians
        constantVector eye( dim * dim );
        vectorTools::eye( eye );
        dTerm2dGreenLagrangeStrain = variableMatrix( term2.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
        dTerm2dMicroStrain         = variableMatrix( term2.size(), variableVector( microStrain.size(), 0 ) );
        dTerm2dInvCPsi             = variableMatrix( term2.size(), variableVector( invCPsi.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int K = 0; K < dim; K++ ){
                            dTerm2dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] += D[ dim * dim * dim * M + dim * dim * N + dim * I + K] * invCPsi[ dim * J + K ];
                            dTerm2dMicroStrain[ dim * I + J ][ dim * M + N ] += B[ dim * dim * dim * I + dim * dim * K + dim * M + N] * invCPsi[ dim * J + K ];
                            for ( unsigned int L = 0; L < dim; L++ ){
                                dTerm2dInvCPsi[ dim * I + J ][ dim * M + N ] += ( B[ dim * dim * dim * I + dim * dim * N + dim * K + L ] * microStrain[ dim * K + L ] + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + N ] ) * eye[ dim * J + M ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress ){
        /*!
         * Compute the higher order stress in the reference configuration.
         * M_{IJK} = C_{JKILMN} Gamma_{LMN}
         *
         * :param const variableVector &Gamma: The micro-gradient deformation measure.
         * :param const parameterVector &C: The C stiffness tensor.
         * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference 
         *     configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( Gamma.size() != dim * dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStress",
                                  "Gamma must have a length of 27" );
        }

        if ( C.size() != dim * dim * dim * dim * dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStress",
                                  "The C stiffness tensor have a length of 3**6" );
        }

        referenceHigherOrderStress = variableVector( dim * dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                referenceHigherOrderStress[ dim * dim * I + dim * J + K ] += C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * L + dim * M + N ] * Gamma[ dim * dim * L + dim * M + N ];
                            }
                        }
                    }
                }
            }
        }
        return NULL;
    }

    errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C, 
                                                variableVector &referenceHigherOrderStress,
                                                variableMatrix &dReferenceHigherOrderStressdGamma ){
        /*!
         * Compute the higher order stress in the reference configuration.
         * M_{IJK} = C_{JKILMN} Gamma_{LMN}
         *
         * Also compute the Jacobian
         * \frac{ \partial M_{IJK} }{\partial \Gamma_{OPQ} } = C_{JKIOPQ}
         *
         * :param const variableVector &Gamma: The micro-gradient deformation measure.
         * :param const parameterVector &C: The C stiffness tensor.
         * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference 
         *     configuration.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceHigherOrderStress (jacobian)",
                                             "Error in computation of higher order stress" );
            result->addNext( error );
            return result;
        }

        //Assemble the Jacobian
        dReferenceHigherOrderStressdGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int O = 0; O < dim; O++ ){
                        for ( unsigned int P = 0; P < dim; P++ ){
                            for ( unsigned int Q = 0; Q < dim; Q++ ){
                                dReferenceHigherOrderStressdGamma[ dim * dim * I + dim * J + K ][ dim * dim * O + dim * P + Q ] +=
                                    C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * O + dim * P + Q ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm3( const variableVector &invCGamma, 
                                        const variableVector &referenceHigherOrderStress, variableVector &term3 ){
        /*!
         * Compute the value of the third term in the micromorphic linear elasticity formulation.
         * term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR}
         *
         * :param const variableVector &invCGamma: C_{JS}^{-1} \Gamma_{SQR}
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param variableVector &term3: The third term in the linear elastic equation.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( invCGamma.size() != dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm3",
                                  "invCGamma must have a size of 27" );
        }

        if ( referenceHigherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeLinearElasticTerm3",
                                  "The referenceHigherOrder stress must have a size of 27" );
        }

        term3 = variableVector( dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        term3[ dim * I + J ] += referenceHigherOrderStress[ dim * dim * I + dim * Q + R ] * invCGamma[ dim * dim * J + dim * Q + R ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                        const variableVector &referenceHigherOrderStress, variableVector &term3,
                                        variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress ){
        /*!
         * Compute the value of the third term in the micromorphic linear elasticity formulation.
         * term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR}
         *
         * Also returns the Jacobians
         * \frac{ \partial term3_{IJ} }{ \partial M_{TUV} } = \delta_{IT} C_{JS}^{-1} \Gamma_{SUV}
         * \frac{ \partial term3_{IJ} }{ \partial C_{TW}^{-1} \Gamma_{WUV} = M_{IUV} \delta_{JT}
         *
         * :param const variableVector &invCGamma: C_{JS}^{-1} \Gamma_{SQR}
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param variableVector &term3: The third term in the linear elastic equation.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress, term3 );

        if ( error ){
            errorOut result = new errorNode( "computeLinearElasticTerm3 (jacobian)",
                                             "Error in computation of term 3" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dTerm3dInvCGamma = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );
        dTerm3dReferenceHigherOrderStress = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int T = 0; T < dim; T++ ){
                    for ( unsigned int U = 0; U < dim; U++ ){
                        for ( unsigned int V = 0; V < dim; V++ ){
                            dTerm3dInvCGamma[ dim * I + J ][ dim * dim * T + dim * U + V ] = referenceHigherOrderStress[ dim * dim * I + dim * U + V ] * eye[ dim * J + T ];
                            dTerm3dReferenceHigherOrderStress[ dim * I + J ][ dim * dim * T + dim * U + V ] = eye[ dim * I + T ] * invCGamma[ dim * dim * J + dim * U + V ];
                        }
                    }
                }
            }
        }
        return NULL;
    }

    errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi ){
        /*!
         * Compute the product C_{IK}^{-1} \Psi_{KJ}
         *
         * :param const variableVector &invRCG: The inverse of the right cauchy green deformation tensor.
         * :param const variableVector &Psi: The micro-deformation measure.
         * :param variableVector &invRCGPsi: the product.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( invRCG.size() != dim * dim ){
            return new errorNode( "computeInvRCGGamma", "invRCG has an improper dimension" );
        }

        if ( Psi.size() != dim * dim ){
            return new errorNode( "computeInvRCGGamma", "Psi has an improper dimension" );
        }

        invRCGPsi = vectorTools::matrixMultiply( invRCG, Psi, dim, dim, dim, dim );

        return NULL;
    }

    errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                               variableMatrix &dInvRCGPsidRCG, variableMatrix &dInvRCGPsidPsi ){
        /*!
         * Compute the product C_{IK}^{-1} \Psi_{KJ}
         *
         * Also compute the Jacobians
         * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial C_{KL} } = -C_{IK}^{-1} C_{LO}^{-1} \Psi_{OJ}
         * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial \Psi_{KL} } = C_{IK}^{-1} \delta_{JL}
         *
         * :param const variableVector &invRCG: The inverse of the right cauchy green deformation tensor.
         * :param const variableVector &Psi: The micro-deformation measure.
         * :param variableVector &invRCGPsi: the product.
         * :param variableMatrix &dInvRCGPsidRCG: The Jacobian of the product w.r.t. the right cauchy green
         *     deformation tensor.
         * :param variableMatrix &dInvRCGPsidPsi: The Jacobian of the product w.r.t. the micro-deformation measure.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeInvRCGPsi( invRCG, Psi, invRCGPsi );

        if ( error ){
            errorOut result = new errorNode( "computeInvRCGPsi (jacobian)", "Error in computation of invRCG Psi product" );
            result->addNext( error );
            return result;
        }

        //Construct the jacobians
        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        dInvRCGPsidRCG = variableMatrix( invRCGPsi.size(), variableVector( invRCG.size(), 0 ) );
        dInvRCGPsidPsi = variableMatrix( invRCGPsi.size(), variableVector( Psi.size(), 0 ) );

        for ( unsigned int I = 0; I < 3; I++ ){
            for ( unsigned int J = 0; J < 3; J++ ){
                for ( unsigned int K = 0; K < 3; K++ ){
                    for ( unsigned int L = 0; L < 3; L++ ){
                        dInvRCGPsidRCG[ dim * I + J ][ dim * K + L ] = -invRCG[ dim * I + K ] * invRCGPsi[ dim * L + J ];
                        dInvRCGPsidPsi[ dim * I + J ][ dim * K + L ] = invRCG[ dim * I + K ] * eye[ dim * J + L ];
                    }
                }
            }
        }
        
        return NULL;
    }

    errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma ){
        /*!
         * Compute the product C_{IS}^{-1} \Gamma_{SQR}
         *
         * :param const variableVector &invRCG: The inverse of the right Cauchy Green deformation tensor.
         * :param const variableVector &Gamma: The gradient of the micro-deformation deformation tensor.
         * :param variableVector &invRCGGamma: The product.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( invRCG.size() != dim * dim ){
            return new errorNode( "computeInvRCGGamma", "invRCG has an improper dimension" );
        }

        if ( Gamma.size() != dim * dim * dim ){
            return new errorNode( "computeInvRCGGamma", "Gamma has an improper dimension" );
        }

        invRCGGamma = variableVector( dim * dim * dim, 0 );
        for ( unsigned int J = 0; J < dim; J++ ){
            for ( unsigned int Q = 0; Q < dim; Q++ ){
                for ( unsigned int R = 0; R < dim; R++ ){
                    for ( unsigned int S = 0; S < dim; S++ ){
                        invRCGGamma[ dim * dim * J + dim * Q + R ] += invRCG[ dim * J + S ] * Gamma[ dim * dim * S + dim * Q + R ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                 variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma ){
        /*!
         * Compute the product C_{IS}^{-1} \Gamma_{SQR}
         *
         * Also compute the Jacobians
         * \frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial C_{TU} } = -C_{JT}^{-1} C_{US}^{-1} \Gamma_{SQR}
          \frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial \Gamma_{TUV} } = C_{JT}^{-1} \delta_{QU} \delta_{RV}
         *
         * :param const variableVector &invRCG: The inverse of the right Cauchy Green deformation tensor.
         * :param const variableVector &Gamma: The gradient of the micro-deformation deformation tensor.
         * :param variableVector &invRCGGamma: The product.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma );

        if ( error ){
            errorOut result = new errorNode( "computeInvRCGGamma (jacobian)", "Error in computation of invRCG Gamma product" );
            result->addNext( error );
            return result;
        }

        //Assemble jacobians of invCGamma w.r.t. C and Gamma
        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        dInvRCGGammadRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dInvRCGGammadGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int J = 0; J < dim; J++ ){
            for ( unsigned int Q = 0; Q < dim; Q++ ){
                for ( unsigned int R = 0; R < dim; R++ ){
                    for ( unsigned int T = 0; T < dim; T++ ){
                        for ( unsigned int U = 0; U < dim; U++ ){
                            dInvRCGGammadRCG[ dim * dim * J + dim * Q + R ][ dim * T + U ] 
                                = -invRCG[ dim * J + T] * invRCGGamma[ dim * dim * U + dim * Q + R ];
                            for ( unsigned int V = 0; V < dim; V++ ){
                                dInvRCGGammadGamma[ dim * dim * J + dim * Q + R ][ dim * dim * T + dim * U + V]
                                    = invRCG[ dim * J + T ] * eye[ dim * Q + U ] * eye[ dim * R + V ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress ){
        /*!
         * Map the stress measures in the reference configuration to the current configuration.
         *
         * :param const variableVector &deformationGradient: The deformation gradient between the 
         *     reference configuration and the current configuration.
         * :param const variableVector &microDeformation: The micro-deformation map between the 
         *     reference configuration and the current configuration.
         * :param const variableVector &PK2Stress: The Second Piola-Kirchoff stress.
         * :param const variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in 
         *     the reference configuration.
         * :param variableVector &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
         * :param variableVector &microStress: The symmetric micro-stress in the current configuration.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Map the PK2 stress to the Cauchy stress
        errorOut error = micromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the PK2 stress to the Cauchy stress" );
            result->addNext( error );
            return result;
        }

        //Map the symmetric micro stress to the current configuration
        error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the micro-stress to the current configuation" );
            result->addNext( error );
            return result;
        }

        //Map the higher order stress to the current configuration
        error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                 microDeformation, higherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the higher-order stress to the current configuation" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                         const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                         const variableVector &referenceHigherOrderStress,
                                         variableVector &cauchyStress, variableVector &microStress,
                                         variableVector &higherOrderStress,
                                         variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                         variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                         variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdXi,
                                         variableMatrix &dHigherOrderStressdReferenceHigherOrderStress ){
        /*!
         * Map the stress measures in the reference configuration to the current configuration.
         *
         * Also computes the Jacobians
         *
         * :param const variableVector &deformationGradient: The deformation gradient between the 
         *     reference configuration and the current configuration.
         * :param const variableVector &microDeformation: The micro-deformation map between the 
         *     reference configuration and the current configuration.
         * :param const variableVector &PK2Stress: The Second Piola-Kirchoff stress.
         * :param const variableVector &referenceMicroStress: The symmetric micro-stress in the 
         *     reference configuration.
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in 
         *     the reference configuration.
         * :param variableVector &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
         * :param variableVector &microStress: The symmetric micro-stress in the current configuration.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param variableMatrix &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the 
         *     deformation gradient.
         * :param variableMatrix &dCauchyStressdPK2Stress: The Jacobian of the Cauchy stress w.r.t. the 
         *     PK2 stress.
         * :param variableMatrix &dMicroStressdF: The Jacobian of the micro stress w.r.t. the 
         *     deformation gradient.
         * :param variableMatrix &dMicroStressdReferenceMicroStress: The Jacobian of the micro-stress 
         *     in the current configuration w.r.t. the micro-stress in the reference configuration.
         * :param variableMatrix &dHigherOrderStressdF: The Jacobian of the higher-order stress w.r.t.
         *     the deformation gradient.
         * :param variableMatrix &dHigherOrderStressdXi: The Jacobian of the higher-order stress 
         *     w.r.t. the micro-deformation.
         * :param variableMatrix &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
         *     higher-order stress w.r.t. the higher order stress in the reference configuration.
         */

        //Map the PK2 stress to the Cauchy stress
        errorOut error = micromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress,
                                                                  dCauchyStressdPK2Stress, dCauchyStressdF );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the PK2 stress to the Cauchy stress" );
            result->addNext( error );
            return result;
        }

        //Map the symmetric micro stress to the current configuration
        error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress,
                                                                    dMicroStressdReferenceMicroStress, dMicroStressdF );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the micro-stress to the current configuation" );
            result->addNext( error );
            return result;
        }

        //Map the higher order stress to the current configuration
        error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                 microDeformation, higherOrderStress,
                                                                 dHigherOrderStressdReferenceHigherOrderStress,
                                                                 dHigherOrderStressdF,
                                                                 dHigherOrderStressdXi );

        if ( error ){
            errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                             "Error in the map of the higher-order stress to the current configuation" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }

    errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                               const variableVector &gradientMicroDeformation,
                               const parameterVector &A, const parameterVector &B, const parameterVector &C,
                               const parameterVector &D,
                               variableVector &cauchyStress, variableVector &microStress,
                               variableVector &higherOrderStress ){
        /*!
         * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &CauchyStress: The Cauchy stress.
         * :param variableVector &microStress: The symmetric micro-stress.
         * :param variableVector &higherOrderStress: The higher-order stress.
         */

        variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
        errorOut error = linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                    A, B, C, D,
                                                    PK2Stress, referenceMicroStress, referenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in the computation of the stresses in the reference configuration" );
            result->addNext( error );
            return result;
        }

        error = mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                            referenceMicroStress, referenceHigherOrderStress,
                                            cauchyStress, microStress, higherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in mapping the reference stresses to the current configuration" );
            result->addNext( error );
            return result;
        }

        return NULL;
    }
    
    errorOut linearElasticity(  const variableVector &deformationGradient, const variableVector &microDeformation,
                                const variableVector &gradientMicroDeformation,
                                const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                const parameterVector &D,
                                variableVector &cauchyStress, variableVector &microStress,
                                variableVector &higherOrderStress,
                                variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdXi, variableMatrix &dCauchyStressdGradXi,
                                variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdXi, variableMatrix &dMicroStressdGradXi,
                                variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdXi,
                                variableMatrix &dHigherOrderStressdGradXi ){
        /*!
         * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
         * of a quadratic decomposition of the energy.
         *
         * Also compute the Jacobians
         *
         * :param const variableVector &deformationGradient: The deformation gradient
         * :param const variableVector &microDeformation: The micro-deformation
         * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
         * :param const parameterVector &A: The A stiffness matrix.
         * :param const parameterVector &B: The B stiffness matrix.
         * :param const parameterVector &C: The C stiffness matrix.
         * :param const parameterVector &D: The D stiffness matrix.
         * :param variableVector &CauchyStress: The Cauchy stress.
         * :param variableVector &microStress: The symmetric micro-stress.
         * :param variableVector &higherOrderStress: The higher-order stress.
         * :param variableMatrix &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the deformation gradient
         * :param variableMatrix &dCauchyStressdXi: The Jacobian of the Cauchy stress w.r.t. the micro deformation.
         * :param variableMatrix &dCauchyStressdGradXi: The Jacobian of the Cauchy stress w.r.t. the gradient of the 
         *     micro-deformation.
         * :param variableMatrix &dMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient
         * :param variableMatrix &dMicroStressdXi: The Jacobian of the Micro stress w.r.t. the micro deformation.
         * :param variableMatrix &dMicroStressdGradXi: The Jacobian of the Micro stress w.r.t. the gradient of the 
         *     micro-deformation.
         * :param variableMatrix &dHigherOrderStressdF: The Jacobian of the Higher Order stress w.r.t. the deformation gradient
         * :param variableMatrix &dHigherOrderStressdXi: The Jacobian of the Higher Order stress w.r.t. the micro deformation.
         * :param variableMatrix &dHigherOrderStressdGradXi: The Jacobian of the Higher Order stress w.r.t. the gradient of the 
         *     micro-deformation.
         */

        variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;

        variableMatrix dPK2StressdF, dPK2StressdXi, dPK2StressdGradXi;
        variableMatrix dReferenceMicroStressdF, dReferenceMicroStressdXi, dReferenceMicroStressdGradXi;
        variableMatrix dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradXi;

        errorOut error = linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                    A, B, C, D,
                                                    PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                    dPK2StressdF, dPK2StressdXi, dPK2StressdGradXi,
                                                    dReferenceMicroStressdF, dReferenceMicroStressdXi, dReferenceMicroStressdGradXi,
                                                    dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradXi );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in the computation of the stresses in the reference configuration" );
            result->addNext( error );
            return result;
        }

        variableMatrix dCauchyStressdPK2Stress, dMicroStressdReferenceMicroStress, dHigherOrderStressdReferenceHigherOrderStress;

        error = mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                            referenceMicroStress, referenceHigherOrderStress,
                                            cauchyStress, microStress, higherOrderStress,
                                            dCauchyStressdF, dCauchyStressdPK2Stress,
                                            dMicroStressdF, dMicroStressdReferenceMicroStress,
                                            dHigherOrderStressdF, dHigherOrderStressdXi,
                                            dHigherOrderStressdReferenceHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "linearElasticity",
                                             "Error in mapping the reference stresses to the current configuration" );
            result->addNext( error );
            return result;
        }

        //Assemble the jacobians of the Cauchy stress
        dCauchyStressdF += vectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdF );
        dCauchyStressdXi = vectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdXi );
        dCauchyStressdGradXi = vectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdGradXi );

        //Assemble the jacobians of the symmetric micro-stress
        dMicroStressdF += vectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdF );
        dMicroStressdXi = vectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdXi );
        dMicroStressdGradXi = vectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdGradXi );

        //Assemble the jacobians of the higher-order stress
        dHigherOrderStressdF += vectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                  dReferenceHigherOrderStressdF );
        dHigherOrderStressdGradXi = vectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                      dReferenceHigherOrderStressdGradXi );

        return NULL;
    }

    errorOut formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A ){
        /*!
         * Form the isotropic A stiffness tensor.
         * A_{KLMN} = \lambda \delta_{KL} \delta_{MN} + \mu \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )
         *
         * :param const parameterType &lambda: The micromorphic lambda parameter.
         * :param const parameterType &mu: The micromorphic mu parameter.
         * :param parameterVector &A: The isotropic A stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        A = parameterVector( dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        A[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = lambda * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                               + mu * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                      + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return NULL;
    }

    errorOut formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                             const parameterType &nu,  const parameterType &sigma, parameterVector &B ){
        /*!
         * Form the isotropic B stiffness tensor.
         * B_{KLMN} = ( eta - tau ) \delta_{KL} \delta_{MN} + \kappa \delta_{KM} \delta_{LN} + \nu \delta_{KN} \delta_{LM}
         *          - \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right)
         *
         * :param const parameterType &eta: The micromorphic eta parameter.
         * :param const parameterType &tau: The micromorphic tau parameter.
         * :param const parameterType &kappa: The micromorphic kappa parameter.
         * :param const parameterType &nu: The micromorphic nu parameter.
         * :param const parameterType &sigma: The micromorphic sigma parameter
         * :param parameterVector &B: The isotropic B stiffnes tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        B = parameterVector( dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        B[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = ( eta - tau ) * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                               + kappa * eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                               + nu * eye[ dim * K + N ] * eye[ dim * L + M ]
                                                                               - sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                         + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return NULL;
    }

    errorOut formIsotropicC( const parameterVector &taus, parameterVector &C ){
        /*!
         * Form the isotropic C stiffness tensor.
         * C_{KLMNPQ} = \tau_1 \left( \delta_{KL} \delta_{MN} \delta_{PQ} + \delta_{KQ} \delta_{LM} \delta_{NP} \right) 
         *            + \tau_2 \left( \delta_{KL} \delta_{MP} \delta_{NQ} + \delta_{KM} \delta_{LQ} \delta_{NP} \right)
         *            + \tau_3 \delta_{KL} \delta_{MQ} \delta_{NP}
         *            + \tau_4 \delta_{KN} \delta_{LM} \delta_{PQ}
         *            + \tau_5 \left( \delta_{KM} \delta_{LN} \delta_{PQ} + \delta_{KP} \delta_{LM} \delta_{NQ} )
         *            + \tau_6 \delta_{KM} \delta_{LP} \delta_{NQ}
         *            + \tau_7 \delta_{KN} \delta_{LP} \delta_{MQ}
         *            + \tau_8 \left( \delta_{KP} \delta_{LQ} \delta_{MN} + \delta_{KQ} \delta_{LN} \delta_{MP} )
         *            + \tau_9 \delta_{KN} \delta_{LQ} \delta_{MP}
         *            + \tau_{10} \delta_{KP} \delta_{LN} \delta_{MQ}
         *            + \tau_{11} \delta_{KQ} \delta_{LP} \delta_{MN}
         *
         * :param const parameterVector &taus: The moduli (11 independent terms)
         * :param parameterVector &C: The isotropic C stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( taus.size() != 11 ){
            return new errorNode( "formIsotropicC", "11 moduli required to form C" );
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        C = parameterVector( dim * dim * dim * dim * dim * dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        for ( unsigned int P = 0; P < dim; P++ ){
                            for ( unsigned int Q = 0; Q < dim; Q++ ){
                                C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M 
                                 + dim * dim * N + dim * P + Q ] = taus[0] * ( eye[ dim * K + L ] * eye[ dim * M + N ] * eye[ dim * P + Q ]
                                                                             + eye[ dim * K + Q ] * eye[ dim * L + M ] * eye[ dim * N + P ] )
                                                                 + taus[1] * ( eye[ dim * K + L ] * eye[ dim * M + P ] * eye[ dim * N + Q ]
                                                                             + eye[ dim * K + M ] * eye[ dim * L + Q ] * eye[ dim * N + P ] )
                                                                 + taus[2] * eye[ dim * K + L ] * eye[ dim * M + Q ] * eye[ dim * N + P]
                                                                 + taus[3] * eye[ dim * K + N ] * eye[ dim * L + M ] * eye[ dim * P + Q]
                                                                 + taus[4] * ( eye[ dim * K + M ] * eye[ dim * L + N ] * eye[ dim * P + Q ]
                                                                             + eye[ dim * K + P ] * eye[ dim * L + M ] * eye[ dim * N + Q ] )
                                                                 + taus[5] * eye[ dim * K + M ] * eye[ dim * L + P ] * eye[ dim * N + Q ]
                                                                 + taus[6] * eye[ dim * K + N ] * eye[ dim * L + P ] * eye[ dim * M + Q ]
                                                                 + taus[7] * ( eye[ dim * K + P ] * eye[ dim * L + Q ] * eye[ dim * M + N ]
                                                                             + eye[ dim * K + Q ] * eye[ dim * L + N ] * eye[ dim * M + P ] )
                                                                 + taus[8] * eye[ dim * K + N ] * eye[ dim * L + Q ] * eye[ dim * M + P ]
                                                                 + taus[9] * eye[ dim * K + P ] * eye[ dim * L + N ] * eye[ dim * M + Q ]
                                                                 + taus[10] * eye[ dim * K + Q ] * eye[ dim * L + P ] * eye[ dim * M + N ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D ) {
        /*!
         * Form the isotropic tensor D.
         * D_{KLMN} = \tau \delta_{KL} \delta_{MN} + \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} )
         *
         * :param const parameterType &tau: The micromorphic tau parameter.
         * :param const parameterType &sigma: The micromorphic sigma parameter.
         * :param parameterVector &D: The D stiffness tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        D = parameterVector( dim * dim * dim * dim, 0 );
        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int L = 0; L < dim; L++ ){
                for ( unsigned int M = 0; M < dim; M++ ){
                    for ( unsigned int N = 0; N < dim; N++ ){
                        D[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = tau * eye[ dim * K + L ] * eye[ dim * M + N ]
                            + sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ] + eye[ dim * K + N ] * eye[ dim * L + M ] );
                    }
                }
            }
        }

        return NULL;
    }
}
