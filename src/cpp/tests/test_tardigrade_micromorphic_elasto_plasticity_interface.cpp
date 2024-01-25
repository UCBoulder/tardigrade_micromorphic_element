//Tests for tardigrade_constitutive_tools

#include<tardigrade_micromorphic_elasto_plasticity.h>
#include<micromorphic_material_library.h>
#include<sstream>
#include<fstream>
#include<iostream>
#include<iomanip>

#define BOOST_TEST_MODULE test_tardigrade_micromorphic_elasto_plasticity_interface
#include <boost/test/included/unit_test.hpp>

typedef tardigradeMicromorphicTools::constantType constantType;
typedef tardigradeMicromorphicTools::constantVector constantVector;
typedef tardigradeMicromorphicTools::constantMatrix constantMatrix;

typedef tardigradeMicromorphicTools::parameterType parameterType;
typedef tardigradeMicromorphicTools::parameterVector parameterVector;
typedef tardigradeMicromorphicTools::parameterMatrix parameterMatrix;

typedef tardigradeMicromorphicTools::variableType variableType;
typedef tardigradeMicromorphicTools::variableVector variableVector;
typedef tardigradeMicromorphicTools::variableMatrix variableMatrix;

typedef tardigradeMicromorphicTools::errorNode errorNode;
typedef tardigradeMicromorphicTools::errorOut errorOut;

BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    BOOST_CHECK( material );

    //Set up the inputs
    //Initialize the time
    std::vector< double > time = { 10., 2.5 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 2.4e2, 1.5e1,             //Macro hardening parameters
                                      2, 1.4e2, 2.0e1,             //Micro hardening parameters
                                      2, 2.0e0, 2.7e1,             //Micro gradient hardening parameters
                                      2, 0.56, 0.2,                //Macro flow parameters
                                      2, 0.15,-0.2,                //Micro flow parameters
                                      2, 0.82, 0.1,                //Micro gradient flow parameters
                                      2, 0.70, 0.3,                //Macro yield parameters
                                      2, 0.40,-0.3,                //Micro yield parameters
                                      2, 0.52, 0.4,                //Micro gradient yield parameters
                                      2, 696.47, 65.84,            //A stiffness tensor parameters
                                      5, -7.69, -51.92, 38.61, -27.31, 5.13,  //B stiffness tensor parameters
                                      11, 1.85, -0.19, -1.08, -1.57, 2.29, -0.61, 5.97, -2.02, 2.38, -0.32, -3.25, //C stiffness tensor parameters
                                      2, -51.92, 5.13,             //D stiffness tensor parameters
                                      0.4, 0.3, 0.35, 1e-8, 1e-8   //Integration parameters
                                    };

    //Initialize the gradient of the macro displacement
//    double current_grad_u[ 3 ][ 3 ] = { { -1.83182277, -0.66558173,  0.23458272 },
//                                        { -0.56632666, -0.21399259,  0.16367238 },
//                                        { -0.29129789, -0.22367825, -2.0632945  } };
//
//    double previous_grad_u[ 3 ][ 3 ] = { { -1.89906429,  0.20890208, -0.39814132 },
//                                         {  0.31303067, -1.23910631, -0.93837662 },
//                                         { -0.32571524, -0.95306342, -0.93025257 } };

    double current_grad_u[ 3 ][ 3 ] = { {0.200, 0.100, 0.000 },
                                        {0.100, 0.001, 0.000 },
                                        {0.000, 0.000, 0.000 } };

    double previous_grad_u[ 3 ][ 3 ] = { {0, 0, 0},
                                         {0, 0, 0},
                                         {0, 0, 0} };
    //Initialize the micro displacement
//    double current_phi[ 9 ] = { 0.84729289,  0.40617104,  0.59534561,  
//                                0.44195587,  0.34121966, -0.79098944, 
//                               -0.43965428,  0.88466225,  0.1684519 };
//
//    double previous_phi[ 9 ] = { -0.99935855, -0.21425717,  0.0668254 ,
//                                 -0.11111872, -0.07416114, -1.01048108,
//                                  0.1804018 , -1.01116291,  0.03248007 };

    double current_phi[ 9 ] = { 0.100, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000 };

    double previous_phi[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {  0.13890017, -0.3598602 , -0.08048856 },
                                          { -0.18572739,  0.06847269,  0.22931628 },
                                          { -0.01829735, -0.48731265, -0.25277529 },
                                          {  0.26626212,  0.4844646 , -0.31965177 },
                                          {  0.49197846,  0.19051656, -0.0365349  },
                                          { -0.06607774, -0.33526875, -0.15803078 },
                                          {  0.09738707, -0.49482218, -0.39584868 },
                                          { -0.45599864,  0.08585038, -0.09432794 },
                                          {  0.23055539,  0.07564162,  0.24051469 } };

//    double previous_grad_phi[ 9 ][ 3 ] = { { -0.47850242,  0.36472234,  0.37071411 },
//                                           {  0.00294417,  0.34480654, -0.34450988 },
//                                           {  0.21056511, -0.28113967, -0.45726839 },
//                                           { -0.26431286, -0.09985721,  0.47322301 },
//                                           { -0.18156887, -0.32226199, -0.37295847 },
//                                           {  0.15062371,  0.09439471,  0.09167948 },
//                                           { -0.46869859,  0.018301  ,  0.45013866 },
//                                           { -0.15455446,  0.40552715, -0.4216042  },
//                                           { -0.38930237,  0.10974753, -0.31188239 } };

//    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0},
//                                          {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

// Old approach
//    tardigradeSolverTools::floatVector PK2_answer = { 172.484,   15.3785,   -0.917177,
//                                             13.4848, 142.823,    -0.0214307,
//                                             -1.7635,   1.77719, 141.069 };
//
//    tardigradeSolverTools::floatVector SIGMA_answer = { 176.916,   15.8646,   -2.83731,
//                                               15.8646, 144.538,     1.85836,
//                                               -2.83731,  1.85836, 142.013 };
//
//    tardigradeSolverTools::floatVector M_answer = { 0.598283, -0.512218,  0.620664,    3.22636,   1.16682,
//                                          1.20593,   0.562825, -2.52317,     1.62616,  -2.61391,
//                                         -0.60994,  -1.02147,   0.668187,    0.49348,  -0.23916,
//                                         -2.77419,   0.760483,  1.71784,    -0.499389,  2.62828,
//                                         -0.761044,  1.23369,  -0.00778206, -2.25643,  -0.729551,
//                                          0.743204,  0.910521 };
//
//    tardigradeSolverTools::floatVector SDVS_answer = { -1.79592e-24,  0.0243222,    0.0822384,    0.0430345,   0.0435752,
//                                             -8.96006e-25,  0.00852191,   0.0465339,    0.0243507,   0.0246566,
//                                              0.00742998,   0.00500421,  -0.000296486,  0.00498757, -0.00260492,
//                                              0.000284355, -0.000367318,  0.000222511, -0.0015603,   0.00863313,
//                                              0.00537105,  -0.000347686,  0.00643802,  -0.00298667,  0.000297105,
//                                             -0.000422398,  0.000293946, -0.00181986,   0.0385522,  -0.0244021,
//                                             -0.00912035,   0.0105171,   -0.000328615,  0.0222843,   0.0185626,
//                                             -0.0118234,   -0.00785555,   0.0451085,    0.031607,    0.0212748,
//                                              0.0116981,    0.0161821,    0.00126031,   0.0147688,   0.00691805,
//                                             -0.0241431,    0.00942608,  -0.0350366,    0.0221571,  -0.0249697,
//                                              0.00935849,   0.0214931,    0.0169609,    0.0177352,   0.0203838 };

// Hydra approach
    tardigradeSolverTools::floatVector PK2_answer = { 1.72374955e+02,  1.53586446e+01, -9.18594711e-01,  1.34613772e+01,
        1.42758448e+02, -2.14353350e-02, -1.76155689e+00,  1.77580288e+00,
        1.41003910e+02 };

    tardigradeSolverTools::floatVector SIGMA_answer = { 176.85221932,  15.84248153,  -2.83823499,  15.84248153,
       144.5183467 ,   1.85778345,  -2.83823499,   1.85778345,
       141.99307014 };

    tardigradeSolverTools::floatVector M_answer = { 0.60016978, -0.51045338,  0.61980289,  3.23507063,  1.16925208,
        1.20665651,  0.56034359, -2.52042105,  1.62648849, -2.61882314,
       -0.61143851, -1.02227749,  0.67046998,  0.49701267, -0.23999964,
       -2.77670511,  0.75636495,  1.71897722, -0.49808019,  2.62569695,
       -0.76002078,  1.23462488, -0.00650376, -2.25591243, -0.73016414,
        0.74380723,  0.90861263 };

    tardigradeSolverTools::floatVector SDVS_answer = {  7.52790780e-03,  5.06936872e-03, -3.00562092e-04,  5.05251313e-03,
       -2.63763651e-03,  2.87928507e-04, -3.72441969e-04,  2.25230220e-04,
       -1.58168041e-03,  7.52714199e-03,  4.60154810e-03, -3.01710599e-04,
        5.56787320e-03, -2.63687070e-03,  2.58160226e-04, -3.65069825e-04,
        2.58160226e-04, -1.58168041e-03,  3.83196421e-02, -2.45387739e-02,
       -8.94156844e-03,  1.04125027e-02, -5.47145710e-04,  2.21268303e-02,
        1.85810005e-02, -1.14241984e-02, -7.62221572e-03,  4.54399478e-02,
        3.21567189e-02,  2.10892064e-02,  1.18343147e-02,  1.63154826e-02,
        1.01776270e-03,  1.47855827e-02,  7.10955065e-03, -2.40560608e-02,
        9.30082927e-03, -3.52878321e-02,  2.18504370e-02, -2.49577455e-02,
        9.18791192e-03,  2.13410715e-02,  1.69866062e-02,  1.77613753e-02,
        2.03937074e-02, -4.80897459e-21,  1.03588470e-02,  6.05863048e-02,
        3.16523418e-02,  3.20619023e-02, -7.99482390e-21,  2.11178609e-02,
        8.23638353e-02,  4.30296627e-02,  4.35864383e-02 };

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
                                            );

    BOOST_CHECK( errorCode == 0 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS, SDVS_answer ) );

    //Check the Jacobian using the previously tested jacobian
    std::vector< std::vector< double > > DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                         DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                         DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer;

    std::vector< std::vector< std::vector< double > > > ADD_JACOBIANS;

    SDVS = SDVSDefault;

//    errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_model(
    errorCode = tardigradeMicromorphicElastoPlasticity::evaluate_hydra_model(
                                time, fparams,
                                current_grad_u, current_phi, current_grad_phi,
                                previous_grad_u, previous_phi, previous_grad_phi,
                                SDVS,
                                current_ADD_DOF, current_ADD_grad_DOF,
                                previous_ADD_DOF, previous_ADD_grad_DOF,
                                PK2_result, SIGMA_result, M_result,
                                DPK2Dgrad_u_answer, DPK2Dphi_answer, DPK2Dgrad_phi_answer,
                                DSIGMADgrad_u_answer, DSIGMADphi_answer, DSIGMADgrad_phi_answer,
                                DMDgrad_u_answer, DMDphi_answer, DMDgrad_phi_answer,
                                ADD_TERMS, ADD_JACOBIANS,
                                output_message
                              );

    BOOST_CHECK( errorCode <= 0 );

    PK2_result.clear();
    SIGMA_result.clear();
    M_result.clear();

    SDVS = SDVSDefault;

    std::vector< std::vector< double > > DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                         DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                         DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result;

    errorCode = material->evaluate_model( time, fparams,
                                          current_grad_u, current_phi, current_grad_phi,
                                          previous_grad_u, previous_phi, previous_grad_phi,
                                          SDVS,
                                          current_ADD_DOF, current_ADD_grad_DOF,
                                          previous_ADD_DOF, previous_ADD_grad_DOF,
                                          PK2_result, SIGMA_result, M_result,
                                          DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                          DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                          DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                          ADD_TERMS, ADD_JACOBIANS,
                                          output_message
                                        );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS_answer, SDVS ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer ) );

    //Test the computed numeric Jacobian values
    SDVS = SDVSDefault;
    errorCode = material->evaluate_model_numeric_gradients( time, fparams,
                                                            current_grad_u, current_phi, current_grad_phi,
                                                            previous_grad_u, previous_phi, previous_grad_phi,
                                                            SDVS,
                                                            current_ADD_DOF, current_ADD_grad_DOF,
                                                            previous_ADD_DOF, previous_ADD_grad_DOF,
                                                            PK2_result, SIGMA_result, M_result,
                                                            DPK2Dgrad_u_result, DPK2Dphi_result, DPK2Dgrad_phi_result,
                                                            DSIGMADgrad_u_result, DSIGMADphi_result, DSIGMADgrad_phi_result,
                                                            DMDgrad_u_result, DMDphi_result, DMDgrad_phi_result,
                                                            ADD_TERMS, ADD_JACOBIANS,
                                                            output_message,
                                                            1e-6 );

    BOOST_CHECK( errorCode <= 0 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2_result, PK2_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMA_result, SIGMA_answer, 1e-5, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( M_result, M_answer, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVS, SDVS_answer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dgrad_u_result, DPK2Dgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dphi_result, DPK2Dphi_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DPK2Dgrad_phi_result, DPK2Dgrad_phi_answer, 1e-4, 1e-5 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADgrad_u_result, DSIGMADgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADphi_result, DSIGMADphi_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DSIGMADgrad_phi_result, DSIGMADgrad_phi_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDgrad_u_result, DMDgrad_u_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDphi_result, DMDphi_answer, 1e-4 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( DMDgrad_phi_result, DMDgrad_phi_answer, 1e-4 ) );

}

BOOST_AUTO_TEST_CASE( testMaterialLibraryInterface2 ){
    /*!
     * Test the interface to the linear elastic model
     * via the material library.
     *
     * NOTE: This function mostly exists to perform debugging
     *       on the implementation of the function into a 
     *       larger solver code.
     *
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    //Set up the inputs
    //Initialize the time
    std::vector< double > time = { 0.045, 0.01 };

    //Initialize the material parameters
    std::vector< double > fparams = { 2, 170, 15, 2, 140, 20, 2, 2, 27, 2, 0.56, 0.2, 2, 0.15, 0.3, 2, 0.82, 0.1, 2, 0.42, 0.3, 2, 0.05, 0.2, 2, 0.52, 0.4, 2, 29480, 25480, 5, 1000, 400, -1500, -1400, -3000, 11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 2, 400, -3000, 0.5, 0.5, 0.5, 1e-09, 1e-09 };

    //Initialize the gradient of the macro displacement
    double current_grad_u[ 3 ][ 3 ] =
    {
        { -0.00124343, -6.55319e-14, 3.99657e-13},
        { 0, 0.0045, 0},
        { -1.75135e-13, -1.35481e-13, -0.00124343 },
    };

    double previous_grad_u[ 3 ][ 3 ] =
    {
        { -0.00123858, -1.22379e-17, 5.04154e-18},
        { 0, 0.004, 0},
        { -1.47723e-18, 4.44523e-18, -0.00123858 },
    };

    //Initialize the micro displacement
    double current_phi[ 9 ] = { -0.00153489, -3.04626e-13, 5.16537e-13, 1.58771e-13, 0.00303407, 4.29828e-14, -4.38368e-13, -1.80694e-13, -0.00153489 };

    double previous_phi[ 9 ] = { -0.00164749, -2.63663e-17, 1.35603e-17, 8.65138e-19, 0.00325613, -2.13082e-20, -1.17433e-17, 2.24626e-18, -0.00164749 };

    //Initialize the gradient of the micro displacement
    double current_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };

    double previous_grad_phi[ 9 ][ 3 ] = { {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0},
                                           {0, 0, 0} };
                                           

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    //Evaluate the model
    int errorCode = material->evaluate_model( time, fparams,
                                              current_grad_u, current_phi, current_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
                                            );

    BOOST_CHECK( errorCode <= 0 );

}

BOOST_AUTO_TEST_CASE( testEvaluate_model_history ){
    /*!
     * Test the material model undergoing a time history.
     *
     */

    //Initialize the model
    std::string _model_name = "LinearElasticityDruckerPragerPlasticity";
    auto &factory = micromorphic_material_library::MaterialFactory::Instance();
    auto material = factory.GetMaterial(_model_name);

    std::vector< std::vector< double > > grad_u_0 = { { 0, 0, 0 },
                                                      { 0, 0, 0 },
                                                      { 0, 0, 0 } };

    std::vector< std::vector< double > > grad_u_f = { { 0.5, 0, 0 },
                                                      { 0.0, 0, 0 },
                                                      { 0.0, 0, 0 } };

    std::vector< double > phi_0 = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    std::vector< double > phi_f = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    std::vector< std::vector< double > > grad_phi_0 = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    std::vector< std::vector< double > > grad_phi_f = { { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0},
                                                        { 0, 0, 0} };

    double dt = 0.05;
    double t0 = 0.;
    double tf = 0.25;

    double t = t0;

    //Set up the model parameters
    std::vector< double > fparams = { 2, 1e3, 1e2,
                                      2, 7e2, 1e4,
                                      2, 1e3, 1e4,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 0., 0.0,
                                      2, 29480, 25480,
                                      5, 1000, 400, -1500, -1400, -3000,
                                      11, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0,
                                      2, 400, -3000,
                                      0.5, 0.5, 0.5, 1e-09, 1e-09 };
//                                      0.0, 0.0, 0.0, 1e-09, 1e-09 };

    //Initialize the state variable vector
    std::vector< double > SDVSDefault( 55, 0 );

    //Initialize the additional degree of freedom vectors
    std::vector< double > current_ADD_DOF;
    std::vector< std::vector< double > > current_ADD_grad_DOF;

    std::vector< double > previous_ADD_DOF;
    std::vector< std::vector< double > > previous_ADD_grad_DOF;

    //Initialize the stress measures
    std::vector< double > current_PK2( 9, 0 );

    std::vector< double > current_SIGMA( 9, 0 );

    std::vector< double > current_M( 27, 0 );

    //Initialize the additional terms vector
    std::vector< std::vector< double > > ADD_TERMS;

    //Initialize the output message string
    std::string output_message;

    std::vector< double > SDVS = SDVSDefault;

    std::vector< double > PK2_result, SIGMA_result, M_result;

    std::vector< double > PK2Answer   = { 5.16027494e+03, -7.22187000e-17,  9.00340000e-20, -7.17613000e-17,
                                          4.05236467e+03,  2.90161000e-19,  8.60690000e-20,  2.07349000e-19,
                                          4.05236467e+03 };

    std::vector< double > SIGMAAnswer = { 5.05182030e+03, -6.76134000e-17,  4.40726000e-20, -6.76136000e-17,
                                          4.11616405e+03,  1.97956000e-19,  4.40303000e-20,  1.98198000e-19,
                                          4.11616405e+03 };

    std::vector< double > MAnswer( 27, 0 );

    std::vector< double > SDVSAnswer = {  4.0658826e-02,  4.8225600e-21, -8.5640000e-24, -2.3326500e-21,
                                         -1.9704307e-02, -2.4829900e-22,  7.1380000e-24,  2.4050100e-22,
                                         -1.9704307e-02,  1.5088110e-02,  4.5280000e-23,  3.4898000e-23,
                                          6.0550000e-23, -7.4596400e-03,  6.4600000e-24, -9.9600000e-24,
                                          1.6290000e-23, -7.4596400e-03, -5.4500000e-25,  1.8890000e-24,
                                         -2.1205000e-23, -4.6510000e-24,  3.3450000e-24, -3.2660000e-26,
                                         -3.4010000e-26,  2.7576000e-25, -8.3000000e-27, -1.9900000e-26,
                                         -8.5180000e-27, -8.4260000e-26, -3.4000000e-28,  4.5730000e-26,
                                          2.1967000e-26,  2.2268000e-28,  6.3990000e-26, -5.1994900e-26,
                                         -4.9700000e-28,  8.2463000e-26, -4.4230000e-26, -4.1900000e-28,
                                         -1.0140000e-26,  7.3940000e-26, -6.2700000e-28,  1.7270000e-26,
                                         -4.8920000e-28,  4.4162700e-01,  0.0000000e+00,  0.0000000e+00,
                                          0.0000000e+00, -8.5310258e-29,  5.3643600e-02,  3.1683300e-02,
                                         -6.2561200e-25,  3.8434100e-25, -9.0882210e-25 };

    std::vector< std::vector< double > > grad_u_prev   = grad_u_0;
    std::vector< double > phi_prev                     = phi_0;
    std::vector< std::vector< double > > grad_phi_prev = grad_phi_0;

    std::vector< std::vector< double > > grad_u_curr;
    std::vector< double > phi_curr;
    std::vector< std::vector< double > > grad_phi_curr;

    std::vector< double > time;

    double current_grad_u[ 3 ][ 3 ], current_phi[ 9 ], current_grad_phi[ 9 ][ 3 ];
    double previous_grad_u[ 3 ][ 3 ], previous_phi[ 9 ], previous_grad_phi[ 9 ][ 3 ];

    //Initial state
    //Update the arrays
    for ( unsigned int i = 0; i < 3; i++ ){
        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
        }
    }

    for ( unsigned int i = 0; i < 9; i++ ){
        previous_phi[ i ] = phi_prev[ i ];

        for ( unsigned int j = 0; j < 3; j++ ){
            previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
        }
    }
    //Evaluate the model
    time = { 0., 0. };
    int errorCode = material->evaluate_model( time, fparams,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              previous_grad_u, previous_phi, previous_grad_phi,
                                              SDVS,
                                              current_ADD_DOF, current_ADD_grad_DOF,
                                              previous_ADD_DOF, previous_ADD_grad_DOF,
                                              PK2_result, SIGMA_result, M_result,
                                              ADD_TERMS,
                                              output_message
                                            );
   
    BOOST_CHECK( errorCode <= 0);

    //Begin iteration
    while ( t + dt < tf ){

        time = { t + dt, dt };

        //Increment the displacements
        grad_u_curr   = grad_u_prev   + dt * ( grad_u_f - grad_u_0 );
        phi_curr      = phi_prev      + dt * ( phi_f - phi_0 );
        grad_phi_curr = grad_phi_prev + dt * ( grad_phi_f - grad_phi_0 );

        //Update the arrays
        for ( unsigned int i = 0; i < 3; i++ ){
            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_u[ i ][ j ]  = grad_u_curr[ i ][ j ];
                previous_grad_u[ i ][ j ] = grad_u_prev[ i ][ j ];
            }
        }

        for ( unsigned int i = 0; i < 9; i++ ){
            current_phi[ i ] = phi_curr[ i ];
            previous_phi[ i ] = phi_prev[ i ];

            for ( unsigned int j = 0; j < 3; j++ ){
                current_grad_phi[ i ][ j ] = grad_phi_curr[ i ][ j ];
                previous_grad_phi[ i ][ j ] = grad_phi_prev[ i ][ j ];
            }
        }

        //Evaluate the model
        int errorCode = material->evaluate_model( time, fparams,
                                                  current_grad_u, current_phi, current_grad_phi,
                                                  previous_grad_u, previous_phi, previous_grad_phi,
                                                  SDVS,
                                                  current_ADD_DOF, current_ADD_grad_DOF,
                                                  previous_ADD_DOF, previous_ADD_grad_DOF,
                                                  PK2_result, SIGMA_result, M_result,
                                                  ADD_TERMS,
                                                  output_message
                                                );

        BOOST_CHECK( errorCode <= 0 );

        t += dt;

        grad_u_prev   = grad_u_curr;
        phi_prev      = phi_curr;
        grad_phi_prev = grad_phi_curr;

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SDVSAnswer, SDVS ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2Answer, PK2_result ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( SIGMAAnswer, SIGMA_result ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( MAnswer, M_result ) );
    
}
