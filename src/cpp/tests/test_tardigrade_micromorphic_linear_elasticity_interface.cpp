// Tests for tardigrade_constitutive_tools

#include <tardigrade_micromorphic_linear_elasticity.h>
#include <tardigrade_micromorphic_linear_elasticity_interface.h>

#include <fstream>
#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_linear_elasticity_interface
#include <boost/test/included/unit_test.hpp>

typedef tardigradeMicromorphicTools::constantType   constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType   parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType   variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

BOOST_AUTO_TEST_CASE(testMaterialLibraryInterface) {
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     * :param std::ofstream &results: The output file.
     */

    // Initialize the model
    std::string _model_name = "LinearElasticity";
    auto       &factory     = micromorphic_material_library::MaterialFactory::Instance();
    auto        material    = factory.GetMaterial(_model_name);

    // Set up the inputs
    const std::vector<double> time = {10, 2.7};

    const std::vector<double> fparams = {2,  1.7, 1.8, 5,  2.8, .76, .15, 9.8, 5.4, 11, 1.,  2.,
                                         3., 4.,  5.,  6., 7.,  8.,  9.,  10., 11., 2,  .76, 5.4};

    const double current_grad_u[3][3] = {
        {-1.07901185, -1.09656192, -0.04629144},
        {-0.77749189, -1.27877771, -0.82648234},
        {0.66484637,  -0.05552567, -1.65125738}
    };

    const double current_phi[9] = {-1.40391532, -0.42715691, 0.75393369, 0.2849511,  -2.06484257,
                                   -0.52190902, 1.07238446,  0.19155907, -0.39704566};

    const double current_grad_phi[9][3] = {
        {0.14940184, 0.12460812, 0.31971128},
        {0.67550862, 0.61095383, 0.87972732},
        {0.30872424, 0.32158187, 0.25480281},
        {0.45570006, 0.69090695, 0.72388584},
        {0.14880964, 0.67520596, 0.15106516},
        {0.77810545, 0.07641724, 0.09367471},
        {0.15905979, 0.0651695,  0.52150417},
        {0.91873444, 0.5622355,  0.50199447},
        {0.26729942, 0.89858519, 0.09043229}
    };

    const double previous_grad_u[3][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    const double previous_phi[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    const double previous_grad_phi[9][3] = {
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0},
        {0, 0, 0}
    };

    std::vector<double>                     SDVS;
    const std::vector<double>               current_ADD_DOF, previous_ADD_DOF;
    const std::vector<std::vector<double> > current_ADD_grad_DOF, previous_ADD_grad_DOF;

    std::vector<double> PK2_answer = {-26.78976487, 91.99831835,  135.04096376, -63.68792655, 149.68226149,
                                      186.67587146, -42.54105342, 125.2317492,  150.55767059};

    std::vector<double> SIGMA_answer = {-47.59920949, 20.84881327, 93.02392773, 20.84881327, 302.43209139,
                                        311.0104045,  93.02392773, 311.0104045, 312.60512922};

    std::vector<double> M_answer = {-50.37283054, -23.25778149, -37.92963077, -19.16962188, -32.97279228, -14.89104497,
                                    -33.4026237,  -15.47947779, -40.31460994, -16.29637436, -36.63942799, -18.22777296,
                                    -39.33546661, -86.69472439, -59.29150146, -15.76480164, -55.42039768, -35.09720118,
                                    -28.94394503, -17.96726082, -45.09734176, -16.46568416, -50.79898863, -39.19129183,
                                    -47.46372724, -42.98201472, -45.57864883};

    std::vector<double>               PK2_result, SIGMA_result, M_result;
    std::vector<std::vector<double> > ADD_TERMS;
    std::string                       output_message;

    // Evaluate the model
    int errorCode = material->evaluate_model(time, fparams, current_grad_u, current_phi, current_grad_phi,
                                             previous_grad_u, previous_phi, previous_grad_phi, SDVS, current_ADD_DOF,
                                             current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result,
                                             SIGMA_result, M_result, ADD_TERMS, output_message);

    BOOST_CHECK(errorCode <= 0);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(PK2_result, PK2_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(SIGMA_result, SIGMA_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(M_result, M_answer));

    // Check the Jacobian using the previously tested jacobian
    std::vector<std::vector<double> > DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer, DSIGMADgrad_u_answer,
        DSIGMADphi_answer, DSIGMADgrad_phi_answer, DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer;

    std::vector<std::vector<std::vector<double> > > ADD_JACOBIANS;

    errorCode = tardigradeMicromorphicLinearElasticity::evaluate_hydra_model(
        time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi,
        SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result, SIGMA_result,
        M_result, DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer, DSIGMADgrad_u_answer, DSIGMADphi_answer,
        DSIGMADgrad_phi_answer, DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer, ADD_TERMS, ADD_JACOBIANS,
        output_message);

    BOOST_CHECK(errorCode <= 0);

    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();

    std::vector<std::vector<double> > DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result, DSIGMADgrad_u_result,
        DSIGMADphi_result, DSIGMADgrad_phi_result, DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result;

    errorCode =
        material->evaluate_model(time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u,
                                 previous_phi, previous_grad_phi, SDVS, current_ADD_DOF, current_ADD_grad_DOF,
                                 previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result, SIGMA_result, M_result,
                                 DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result, DSIGMADgrad_u_result,
                                 DSIGMADphi_result, DSIGMADgrad_phi_result, DMDgrad_u_result, DMDphi_result,
                                 DMDgrad_phi_result, ADD_TERMS, ADD_JACOBIANS, output_message);

    BOOST_CHECK(errorCode == 0);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(PK2_result, PK2_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(SIGMA_result, SIGMA_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(M_result, M_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dgrad_u_result, DPK2Dgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dphi_result, DPK2Dphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADgrad_u_result, DSIGMADgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADphi_result, DSIGMADphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDgrad_u_result, DMDgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDphi_result, DMDphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDgrad_phi_result, DMDgrad_phi_answer));

    // Test the computed numeric Jacobian values
    errorCode = material->evaluate_model_numeric_gradients(
        time, fparams, current_grad_u, current_phi, current_grad_phi, previous_grad_u, previous_phi, previous_grad_phi,
        SDVS, current_ADD_DOF, current_ADD_grad_DOF, previous_ADD_DOF, previous_ADD_grad_DOF, PK2_result, SIGMA_result,
        M_result, DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result, DSIGMADgrad_u_result, DSIGMADphi_result,
        DSIGMADgrad_phi_result, DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result, ADD_TERMS, ADD_JACOBIANS,
        output_message, 1e-6);

    BOOST_CHECK(errorCode <= 0);

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(PK2_result, PK2_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(SIGMA_result, SIGMA_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(M_result, M_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dgrad_u_result, DPK2Dgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dphi_result, DPK2Dphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADgrad_u_result, DSIGMADgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADphi_result, DSIGMADphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDgrad_u_result, DMDgrad_u_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDphi_result, DMDphi_answer));

    BOOST_CHECK(tardigradeVectorTools::fuzzyEquals(DMDgrad_phi_result, DMDgrad_phi_answer));
}
