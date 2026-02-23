/*!
=====================================================================
|                 micromorphic_material_library.cpp                 |
=====================================================================
| A header file which defines a class which registers all of the    |
| available micromorphic material models and their interface. The   |
| library is implemented in such a way that someone can register a  |
| new micromorphic constitutive model and have it available for use |
| merely by calling it by name. This allows us to re-use code as    |
| much as possible without needing to do time consuming rewrites of |
| already existing code.                                            |
---------------------------------------------------------------------
| Note: Registration approach taken from stackoverflow question     |
|       compile time plugin system 2                                |
=====================================================================
*/

#include "micromorphic_material_library.h"

#include <iostream>

namespace micromorphic_material_library {

    // IMaterial::~Imaterial(){}

    // IMaterialRegistrar::~IMaterialRegistrar(){}

    int IMaterial::evaluate_model_numeric_gradients(
        const std::vector<double> &time, const std::vector<double>(&fparams), const double (&current_grad_u)[3][3],
        const double (&current_phi)[9], const double (&current_grad_phi)[9][3], const double (&previous_grad_u)[3][3],
        const double (&previous_phi)[9], const double (&previous_grad_phi)[9][3], std::vector<double> &SDVS,
        const std::vector<double> &current_ADD_DOF, const std::vector<std::vector<double> > &current_ADD_grad_DOF,
        const std::vector<double> &previous_ADD_DOF, const std::vector<std::vector<double> > &previous_ADD_grad_DOF,
        std::vector<double> &PK2, std::vector<double> &SIGMA, std::vector<double> &M,
        std::vector<std::vector<double> > &DPK2Dgrad_u, std::vector<std::vector<double> > &DPK2Dphi,
        std::vector<std::vector<double> > &DPK2Dgrad_phi, std::vector<std::vector<double> > &DSIGMADgrad_u,
        std::vector<std::vector<double> > &DSIGMADphi, std::vector<std::vector<double> > &DSIGMADgrad_phi,
        std::vector<std::vector<double> > &DMDgrad_u, std::vector<std::vector<double> > &DMDphi,
        std::vector<std::vector<double> > &DMDgrad_phi, std::vector<std::vector<double> > &ADD_TERMS,
        std::vector<std::vector<std::vector<double> > > &ADD_JACOBIANS, std::string &output_message,
#ifdef DEBUG_MODE
        std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > &DEBUG,

#endif
        double delta) {

        /*!
         * Evaluate the jacobian of the material model using a numeric gradient
         * TODO: Add Jacobians of the additional terms
         *
         * :param const std::vector< double > &time: The current time and the timestep
         *     [ current_t, dt ]
         * :param const std::vector< double > ( &fparams ): The parameters for the constitutive model
         * :param const double ( &current_grad_u )[ 3 ][ 3 ]: The current displacement gradient
         *     Assumed to be of the form [ [ u_{1,1}, u_{1,2}, u_{1,3} ],
         *                                 [ u_{2,1}, u_{2,2}, u_{2,3} ],
         *                                 [ u_{3,1}, u_{3,2}, u_{3,3} ] ]
         * :param const double ( &current_phi )[ 9 ]: The current micro displacment values.
         *     Assumed to be of the form [ \phi_{11}, \phi_{12}, \phi_{13}, \phi_{21}, \phi_{22}, \phi_{23}, \phi_{31},
         * \phi_{32}, \phi_{33} ] :param const double ( &current_grad_phi )[ 9 ][ 3 ]: The current micro displacement
         * gradient Assumed to be of the form [ [ \phi_{11,1}, \phi_{11,2}, \phi_{11,3} ], [ \phi_{12,1}, \phi_{12,2},
         * \phi_{12,3} ], [ \phi_{13,1}, \phi_{13,2}, \phi_{13,3} ], [ \phi_{21,1}, \phi_{21,2}, \phi_{21,3} ], [
         * \phi_{22,1}, \phi_{22,2}, \phi_{22,3} ], [ \phi_{23,1}, \phi_{23,2}, \phi_{23,3} ], [ \phi_{31,1},
         * \phi_{31,2}, \phi_{31,3} ], [ \phi_{32,1}, \phi_{32,2}, \phi_{32,3} ], [ \phi_{33,1}, \phi_{33,2},
         * \phi_{33,3} ] ] :param const double ( &previous_grad_u )[ 3 ][ 3 ]: The previous displacement gradient.
         * :param const double ( &previous_phi )[ 9 ]: The previous micro displacement.
         * :param const double ( &previous_grad_phi )[ 9 ][ 3 ]: The previous micro displacement gradient.
         * :param std::vector< double > &SDVS: The previously converged values of the state variables
         * :param std::vector< double > &current_ADD_DOF: The current values of the additional degrees of freedom
         * :param std::vector< std::vector< double > > &current_ADD_grad_DOF: The current values of the gradients of the
         *     additional degrees of freedom
         * :param std::vector< double > &current_PK2: The current value of the second Piola Kirchhoff stress tensor. The
         * format is [ S_{11}, S_{12}, S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ] :param std::vector<
         * double > &current_SIGMA: The current value of the reference micro stress. The format is [ S_{11}, S_{12},
         * S_{13}, S_{21}, S_{22}, S_{23}, S_{31}, S_{32}, S_{33} ] :param std::vector< double > &current_M: The current
         * value of the reference higher order stress. The format is [ M_{111}, M_{112}, M_{113}, M_{121}, M_{122},
         * M_{123}, M_{131}, M_{132}, M_{133}, M_{211}, M_{212}, M_{213}, M_{221}, M_{222}, M_{223}, M_{231}, M_{232},
         * M_{233}, M_{311}, M_{312}, M_{313}, M_{321}, M_{322}, M_{323}, M_{331}, M_{332}, M_{333} ] :param
         * std::vector< std::vector< double > > &DPK2Dgrad_u: The Jacobian of the PK2 stress w.r.t. the gradient of
         * macro displacement. :param std::vector< std::vector< double > > &DPK2Dphi: The Jacobian of the PK2 stress
         * w.r.t. the micro displacement. :param std::vector< std::vector< double > > &DPK2Dgrad_phi: The Jacobian of
         * the PK2 stress w.r.t. the gradient of the micro displacement. :param std::vector< std::vector< double > >
         * &DSIGMAdgrad_u: The Jacobian of the reference symmetric micro stress w.r.t. the gradient of the macro
         * displacement. :param std::vector< std::vector< double > > &DSIGMAdphi: The Jacobian of the reference
         * symmetric micro stress w.r.t. the micro displacement. :param std::vector< std::vector< double > >
         * &DSIGMAdgrad_phi: The Jacobian of the reference symmetric micro stress w.r.t. the gradient of the micro
         * displacement. :param std::vector< std::vector< double > > &DMDgrad_u: The Jacobian of the reference higher
         * order stress w.r.t. the gradient of the macro displacement. :param std::vector< std::vector< double > >
         * &DMDphi: The Jacobian of the reference higher order stress w.r.t. the micro displacement. :param std::vector<
         * std::vector< double > > &DMDgrad_phi: The Jacobian of the reference higher order stress w.r.t. the gradient
         * of the micro displacement. :param std::vector< std::vector< double > > &ADD_TERMS: Additional terms :param
         * std::vector< std::vector< std::vector< double > > > &ADD_JACOBIANS: The jacobians of the additional terms
         * w.r.t. the deformation :param std::string &output_message: The output message string. :param double delta =
         * 1e-6: The perturbation to be applied to the incoming degrees of freedom.
         *
         * Returns:
         *     0: No errors. Solution converged.
         *     1: Convergence Error. Request timestep cutback.
         *     2: Fatal Errors encountered. Terminate the simulation.
         */

        // Evaluate the model at the set point
        std::vector<double> SDVS_previous = SDVS;
        int errorCode = evaluate_model(time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u,
                                       previous_phi, previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF,
                                       previous_ADD_DOF, previous_ADD_grad_DOF, PK2, SIGMA, M, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                                       ,
                                       DEBUG
#endif
        );

#ifdef DEBUG_MODE
        DEBUG.clear();
#endif

        if (errorCode > 0) {
            return errorCode;
        }

        // Assemble the Jacobians w.r.t. the macro displacement gradient
        DPK2Dgrad_u   = std::vector<std::vector<double> >(PK2.size(), std::vector<double>(9, 0));
        DSIGMADgrad_u = std::vector<std::vector<double> >(SIGMA.size(), std::vector<double>(9, 0));
        DMDgrad_u     = std::vector<std::vector<double> >(M.size(), std::vector<double>(9, 0));
        for (unsigned int i = 0; i < 9; i++) {
            std::vector<std::vector<double> > deltaMatrix(3, std::vector<double>(3, 0));
            deltaMatrix[i / 3][i % 3] = delta * fabs(current_grad_u[i / 3][i % 3]) + delta;

            double current_grad_u_P[3][3] = {
                {current_grad_u[0][0] + deltaMatrix[0][0], current_grad_u[0][1] + deltaMatrix[0][1],
                 current_grad_u[0][2] + deltaMatrix[0][2]},
                {current_grad_u[1][0] + deltaMatrix[1][0], current_grad_u[1][1] + deltaMatrix[1][1],
                 current_grad_u[1][2] + deltaMatrix[1][2]},
                {current_grad_u[2][0] + deltaMatrix[2][0], current_grad_u[2][1] + deltaMatrix[2][1],
                 current_grad_u[2][2] + deltaMatrix[2][2]}
            };

            double current_grad_u_M[3][3] = {
                {current_grad_u[0][0] - deltaMatrix[0][0], current_grad_u[0][1] - deltaMatrix[0][1],
                 current_grad_u[0][2] - deltaMatrix[0][2]},
                {current_grad_u[1][0] - deltaMatrix[1][0], current_grad_u[1][1] - deltaMatrix[1][1],
                 current_grad_u[1][2] - deltaMatrix[1][2]},
                {current_grad_u[2][0] - deltaMatrix[2][0], current_grad_u[2][1] - deltaMatrix[2][1],
                 current_grad_u[2][2] - deltaMatrix[2][2]}
            };

            std::vector<double> PK2_P, SIGMA_P, M_P, PK2_M, SIGMA_M, M_M;

            std::vector<double> SDVS_P, SDVS_M;
            SDVS_P = SDVS_previous;
            SDVS_M = SDVS_previous;

#ifdef DEBUG_MODE
            std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > DEBUG_P,
                DEBUG_M;
#endif

            errorCode =
                evaluate_model(time, fparams, current_grad_u_P, current_phi, current_grad_phi, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_P, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_P, SIGMA_P, M_P, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_P
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            errorCode =
                evaluate_model(time, fparams, current_grad_u_M, current_phi, current_grad_phi, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_M, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_M, SIGMA_M, M_M, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_M
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            for (unsigned int j = 0; j < PK2.size(); j++) {
                DPK2Dgrad_u[j][i]   = (PK2_P[j] - PK2_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
                DSIGMADgrad_u[j][i] = (SIGMA_P[j] - SIGMA_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
            }
            for (unsigned int j = 0; j < M.size(); j++) {
                DMDgrad_u[j][i] = (M_P[j] - M_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
            }
        }

        // Assemble the Jacobians w.r.t. the micro displacement
        DPK2Dphi   = std::vector<std::vector<double> >(9, std::vector<double>(9, 0));
        DSIGMADphi = std::vector<std::vector<double> >(9, std::vector<double>(9, 0));
        DMDphi     = std::vector<std::vector<double> >(27, std::vector<double>(9, 0));
        for (unsigned int i = 0; i < 9; i++) {
            std::vector<double> deltaVector(9, 0);
            deltaVector[i] = delta * fabs(current_phi[i]) + delta;

            double current_phi_P[9] = {
                current_phi[0] + deltaVector[0], current_phi[1] + deltaVector[1], current_phi[2] + deltaVector[2],
                current_phi[3] + deltaVector[3], current_phi[4] + deltaVector[4], current_phi[5] + deltaVector[5],
                current_phi[6] + deltaVector[6], current_phi[7] + deltaVector[7], current_phi[8] + deltaVector[8]};

            double current_phi_M[9] = {
                current_phi[0] - deltaVector[0], current_phi[1] - deltaVector[1], current_phi[2] - deltaVector[2],
                current_phi[3] - deltaVector[3], current_phi[4] - deltaVector[4], current_phi[5] - deltaVector[5],
                current_phi[6] - deltaVector[6], current_phi[7] - deltaVector[7], current_phi[8] - deltaVector[8]};

            std::vector<double> PK2_P, SIGMA_P, M_P, PK2_M, SIGMA_M, M_M;

            std::vector<double> SDVS_P, SDVS_M;
            SDVS_P = SDVS_previous;
            SDVS_M = SDVS_previous;

#ifdef DEBUG_MODE
            std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > DEBUG_P,
                DEBUG_M;
#endif

            errorCode =
                evaluate_model(time, fparams, current_grad_u, current_phi_P, current_grad_phi, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_P, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_P, SIGMA_P, M_P, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_P
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            errorCode =
                evaluate_model(time, fparams, current_grad_u, current_phi_M, current_grad_phi, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_M, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_M, SIGMA_M, M_M, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_M
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            for (unsigned int j = 0; j < PK2.size(); j++) {
                DPK2Dphi[j][i]   = (PK2_P[j] - PK2_M[j]) / (2 * deltaVector[i]);
                DSIGMADphi[j][i] = (SIGMA_P[j] - SIGMA_M[j]) / (2 * deltaVector[i]);
            }
            for (unsigned int j = 0; j < M.size(); j++) {
                DMDphi[j][i] = (M_P[j] - M_M[j]) / (2 * deltaVector[i]);
            }
        }

        // Assemble the Jacobians w.r.t. the gradient of the micro displacement
        DPK2Dgrad_phi   = std::vector<std::vector<double> >(PK2.size(), std::vector<double>(27, 0));
        DSIGMADgrad_phi = std::vector<std::vector<double> >(SIGMA.size(), std::vector<double>(27, 0));
        DMDgrad_phi     = std::vector<std::vector<double> >(M.size(), std::vector<double>(27, 0));

        for (unsigned int i = 0; i < 27; i++) {
            std::vector<std::vector<double> > deltaMatrix(9, std::vector<double>(3, 0));
            deltaMatrix[i / 3][i % 3] = delta * fabs(current_grad_phi[i / 3][i % 3]) + delta;

            double current_grad_phi_P[9][3] = {
                {current_grad_phi[0][0] + deltaMatrix[0][0], current_grad_phi[0][1] + deltaMatrix[0][1],
                 current_grad_phi[0][2] + deltaMatrix[0][2]},
                {current_grad_phi[1][0] + deltaMatrix[1][0], current_grad_phi[1][1] + deltaMatrix[1][1],
                 current_grad_phi[1][2] + deltaMatrix[1][2]},
                {current_grad_phi[2][0] + deltaMatrix[2][0], current_grad_phi[2][1] + deltaMatrix[2][1],
                 current_grad_phi[2][2] + deltaMatrix[2][2]},
                {current_grad_phi[3][0] + deltaMatrix[3][0], current_grad_phi[3][1] + deltaMatrix[3][1],
                 current_grad_phi[3][2] + deltaMatrix[3][2]},
                {current_grad_phi[4][0] + deltaMatrix[4][0], current_grad_phi[4][1] + deltaMatrix[4][1],
                 current_grad_phi[4][2] + deltaMatrix[4][2]},
                {current_grad_phi[5][0] + deltaMatrix[5][0], current_grad_phi[5][1] + deltaMatrix[5][1],
                 current_grad_phi[5][2] + deltaMatrix[5][2]},
                {current_grad_phi[6][0] + deltaMatrix[6][0], current_grad_phi[6][1] + deltaMatrix[6][1],
                 current_grad_phi[6][2] + deltaMatrix[6][2]},
                {current_grad_phi[7][0] + deltaMatrix[7][0], current_grad_phi[7][1] + deltaMatrix[7][1],
                 current_grad_phi[7][2] + deltaMatrix[7][2]},
                {current_grad_phi[8][0] + deltaMatrix[8][0], current_grad_phi[8][1] + deltaMatrix[8][1],
                 current_grad_phi[8][2] + deltaMatrix[8][2]}
            };

            double current_grad_phi_M[9][3] = {
                {current_grad_phi[0][0] - deltaMatrix[0][0], current_grad_phi[0][1] - deltaMatrix[0][1],
                 current_grad_phi[0][2] - deltaMatrix[0][2]},
                {current_grad_phi[1][0] - deltaMatrix[1][0], current_grad_phi[1][1] - deltaMatrix[1][1],
                 current_grad_phi[1][2] - deltaMatrix[1][2]},
                {current_grad_phi[2][0] - deltaMatrix[2][0], current_grad_phi[2][1] - deltaMatrix[2][1],
                 current_grad_phi[2][2] - deltaMatrix[2][2]},
                {current_grad_phi[3][0] - deltaMatrix[3][0], current_grad_phi[3][1] - deltaMatrix[3][1],
                 current_grad_phi[3][2] - deltaMatrix[3][2]},
                {current_grad_phi[4][0] - deltaMatrix[4][0], current_grad_phi[4][1] - deltaMatrix[4][1],
                 current_grad_phi[4][2] - deltaMatrix[4][2]},
                {current_grad_phi[5][0] - deltaMatrix[5][0], current_grad_phi[5][1] - deltaMatrix[5][1],
                 current_grad_phi[5][2] - deltaMatrix[5][2]},
                {current_grad_phi[6][0] - deltaMatrix[6][0], current_grad_phi[6][1] - deltaMatrix[6][1],
                 current_grad_phi[6][2] - deltaMatrix[6][2]},
                {current_grad_phi[7][0] - deltaMatrix[7][0], current_grad_phi[7][1] - deltaMatrix[7][1],
                 current_grad_phi[7][2] - deltaMatrix[7][2]},
                {current_grad_phi[8][0] - deltaMatrix[8][0], current_grad_phi[8][1] - deltaMatrix[8][1],
                 current_grad_phi[8][2] - deltaMatrix[8][2]}
            };

            std::vector<double> PK2_P, SIGMA_P, M_P, PK2_M, SIGMA_M, M_M;

            std::vector<double> SDVS_P, SDVS_M;
            SDVS_P = SDVS_previous;
            SDVS_M = SDVS_previous;

#ifdef DEBUG_MODE
            std::map<std::string, std::map<std::string, std::map<std::string, std::vector<double> > > > DEBUG_P,
                DEBUG_M;
#endif

            errorCode =
                evaluate_model(time, fparams, current_grad_u, current_phi, current_grad_phi_P, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_P, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_P, SIGMA_P, M_P, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_P
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            errorCode =
                evaluate_model(time, fparams, current_grad_u, current_phi, current_grad_phi_M, previous_grad_u,
                               previous_phi, previous_grad_phi, SDVS_M, current_ADD_DOF, current_ADD_grad_DOF,
                               previous_ADD_DOF, previous_ADD_grad_DOF, PK2_M, SIGMA_M, M_M, ADD_TERMS, output_message
#ifdef DEBUG_MODE
                               ,
                               DEBUG_P
#endif
                );

            if (errorCode > 0) {
                return errorCode;
            }

            for (unsigned int j = 0; j < PK2.size(); j++) {
                DPK2Dgrad_phi[j][i]   = (PK2_P[j] - PK2_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
                DSIGMADgrad_phi[j][i] = (SIGMA_P[j] - SIGMA_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
            }
            for (unsigned int j = 0; j < M.size(); j++) {
                DMDgrad_phi[j][i] = (M_P[j] - M_M[j]) / (2 * deltaMatrix[i / 3][i % 3]);
            }
        }

        return 0;
    }

    MaterialFactory &MaterialFactory::Instance() {
        static MaterialFactory instance;
        return instance;
    }

    void MaterialFactory::Register(IMaterialRegistrar *registrar, std::string name) { registry_[name] = registrar; }

    std::unique_ptr<IMaterial> MaterialFactory::GetMaterial(std::string name) {
        /* throws out_of_range if material unknown */
        try {
            IMaterialRegistrar *registrar = registry_.at(name);
            return registrar->GetMaterial();
        } catch (...) {
            std::string message = "Exception when attempting to access: " + name + "\n";
            message += "material models available are:\n";
            for (auto it = registry_.begin(); it != registry_.end(); it++) {
                message += it->first + "\n";
            }
            throw std::runtime_error(message);
        }
        return NULL;
    }

    void MaterialFactory::PrintMaterials() {
        /*! Prints all of the materials registered in the library*/
        std::string message = "Materials available in the library:\n";
        for (auto it = registry_.begin(); it != registry_.end(); it++) {
            message += "    " + it->first + "\n";
        }
        std::cerr << message;
    }

}  // namespace micromorphic_material_library
