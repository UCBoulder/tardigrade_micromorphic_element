#ifndef TARDIGRADE_MICROMORPHIC_LINEAR_ELASTICITY_INTERFACE_H
#define TARDIGRADE_MICROMORPHIC_LINEAR_ELASTICITY_INTERFACE_H

#include<tardigrade_micromorphic_linear_elasticity.h>
#include<micromorphic_material_library.h>
namespace tardigradeMicromorphicLinearElasticity{
    class LinearElasticity: public micromorphic_material_library::IMaterial{
        /*!
         * The class which is called when evaluating a
         * linear elastic micromorphic constitutive model.
         *
         * This class registers linear elasticity into a library of models which can then
         * be called by name.
         */

        public:
            int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                std::vector< double > &SDVS,
                                const std::vector< double > &current_ADD_DOF,
                                const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                const std::vector< double > &previous_ADD_DOF,
                                const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                                std::vector< std::vector< double > > &ADD_TERMS,
                                std::string &output_message
#ifdef DEBUG_MODE
                                , std::map< std::string, std::map< std::string, std::map< std::string, std::vector< double > > > > &DEBUG
#endif
                              )
            {
                return tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
                                         time, fparams,
                                         current_grad_u, current_phi, current_grad_phi,
                                         previous_grad_u, previous_phi, previous_grad_phi,
                                         SDVS,
                                         current_ADD_DOF, current_ADD_grad_DOF,
                                         previous_ADD_DOF, previous_ADD_grad_DOF,
                                         PK2, SIGMA, M,
                                         ADD_TERMS,
                                         output_message
                                         );
            }

            int evaluate_model( const std::vector< double > &time,            const std::vector< double > ( &fparams ),
                                const double ( &current_grad_u )[ 3 ][ 3 ],   const double ( &current_phi )[ 9 ],
                                const double ( &current_grad_phi )[ 9 ][ 3 ], const double ( &previous_grad_u )[ 3 ][ 3 ],
                                const double ( &previous_phi )[ 9 ],          const double ( &previous_grad_phi )[ 9 ][ 3 ],
                                std::vector< double > &SDVS,
                                const std::vector< double > &current_ADD_DOF,
                                const std::vector< std::vector< double > > &current_ADD_grad_DOF,
                                const std::vector< double > &previous_ADD_DOF,
                                const std::vector< std::vector< double > > &previous_ADD_grad_DOF,
                                std::vector< double > &PK2, std::vector< double > &SIGMA, std::vector< double > &M,
                                std::vector< std::vector< double > > &DPK2Dgrad_u,
                                std::vector< std::vector< double > > &DPK2Dphi,
                                std::vector< std::vector< double > > &DPK2Dgrad_phi,
                                std::vector< std::vector< double > > &DSIGMADgrad_u,
                                std::vector< std::vector< double > > &DSIGMADphi,
                                std::vector< std::vector< double > > &DSIGMADgrad_phi,
                                std::vector< std::vector< double > > &DMDgrad_u,
                                std::vector< std::vector< double > > &DMDphi,
                                std::vector< std::vector< double > > &DMDgrad_phi,
                                std::vector< std::vector< double > > &ADD_TERMS,
                                std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS,
                                std::string &output_message
#ifdef DEBUG_MODE
                                , std::map< std::string, std::map< std::string, std::map< std::string, std::vector< double > > > > &DEBUG
#endif
                              )
            {
                return tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
                                         time, fparams,
                                         current_grad_u, current_phi, current_grad_phi,
                                         previous_grad_u, previous_phi, previous_grad_phi,
                                         SDVS,
                                         current_ADD_DOF, current_ADD_grad_DOF,
                                         previous_ADD_DOF, previous_ADD_grad_DOF,
                                         PK2, SIGMA, M,
                                         DPK2Dgrad_u, DPK2Dphi, DPK2Dgrad_phi,
                                         DSIGMADgrad_u, DSIGMADphi, DSIGMADgrad_phi,
                                         DMDgrad_u, DMDphi, DMDgrad_phi,
                                         ADD_TERMS, ADD_JACOBIANS,
                                         output_message
                                         );
            }
    };

    REGISTER_MATERIAL(LinearElasticity)
}

#endif
