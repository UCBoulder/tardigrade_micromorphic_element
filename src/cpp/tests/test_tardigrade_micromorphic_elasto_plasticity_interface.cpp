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
    tardigradeSolverTools::floatVector PK2_answer = { 176.85219663,  11.80698447,  -0.42826398,
                                                        9.69466718, 154.70890135,  -0.34866071,
                                                       -1.32423993,   1.55229825, 152.89259867 };

    tardigradeSolverTools::floatVector SIGMA_answer = { 179.90985149,  12.12517822,  -2.34267185,
                                                         12.12517822, 155.18595135,   1.57329507,
                                                         -2.34267185,   1.57329507, 152.67097908 };

    tardigradeSolverTools::floatVector M_answer = { 0.55375154, -0.53513182,  0.64105973,  3.02610427,  1.18082958,
                                                    1.15933657,  0.60697669, -2.62337981,  1.55437107, -2.49773592,
                                                   -0.69264839, -0.86588199,  0.70376012,  0.44225752, -0.21927807,
                                                   -2.77735624,  0.95102462,  1.70382866, -0.48471999,  2.69244538,
                                                   -0.62391716,  1.1710374 , -0.10182978, -2.27621329, -0.7875813 ,
                                                    0.752404  ,  0.99367537 };

    tardigradeSolverTools::floatVector SDVS_answer = { 0.00988018 ,  0.011861  , -0.000863237,  0.0110431  , -0.0133071 ,
                                                       0.000535357, -0.00150177,  0.00111378 , -0.0115157  ,  0.0419332 ,
                                                       0.0248495  , -0.00137643,  0.0300679  , -0.012865   ,  0.00100127,
                                                      -0.00166548 ,  0.00100127, -0.00863817 ,  0.0524155  , -0.0192113 ,
                                                      -0.0190793  ,  0.0137043 ,  0.0109677  ,  0.0343499  ,  0.0171499 ,
                                                      -0.0366952  , -0.0208241 ,  0.0284659  ,  1.45232e-05,  0.034603  ,
                                                       0.00526692 ,  0.0108779 ,  0.0137381  ,  0.0140525  , -0.00342129,
                                                      -0.0290165  ,  0.0155877 , -0.0185748  ,  0.0431296  , -0.0287399 ,
                                                       0.0195558  ,  0.0308138 ,  0.0160433  ,  0.0167955  ,  0.0198817 ,
                                                      -0.031206   ,  0.0553656 ,  0.0604235  ,  0.0343974  ,  0.0310536 ,
                                                      -0.0521235  ,  0.11287   ,  0.0821425  ,  0.0467614  ,  0.0422157 };

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

    std::vector< double > PK2Answer   = { 5.16350e+03, -3.00694e-17, -4.15205e-17,
                                         -2.94815e-17,  4.05388e+03,  3.19499e-17,
                                         -4.21351e-17,  3.07368e-17,  4.05388e+03 };

    std::vector< double > SIGMAAnswer = { 5.06571e+03, -4.04083e-17, -3.78865e-17,
                                         -3.90853e-17,  4.11088e+03,  7.82938e-18,
                                         -3.63730e-17,  1.43787e-17,  4.11088e+03 };

    std::vector< double > MAnswer( 27, 0 );

    std::vector< double > SDVSAnswer = { 4.03683e-02,  3.55329e-21,  7.12965e-22, -2.10963e-21,
                                        -1.95818e-02,  1.22340e-22,  7.98957e-22, -6.98456e-22,
                                        -1.95818e-02,  1.07605e-02,  1.46563e-21, -2.59163e-22,
                                         6.84998e-21, -5.30058e-03,  4.72773e-21, -4.80356e-22,
                                         5.58748e-21, -5.30058e-03,  4.21892e-22,  5.95521e-22,
                                         5.94044e-22,  7.77887e-22, -5.20839e-23,  3.50324e-24,
                                         7.23258e-24,  3.52995e-24, -4.23086e-24,  1.34296e-23,
                                         5.46166e-25,  6.86559e-24, -3.72983e-25, -6.95441e-25,
                                        -2.45526e-25,  5.40325e-27, -1.01136e-24,  4.06927e-26,
                                        -4.03804e-25, -8.94997e-25, -3.66759e-24, -3.66195e-25,
                                        -1.59872e-24, -2.52535e-24,  1.64963e-25,  6.20733e-25,
                                        -3.11079e-26, -2.39800e-01,  5.28270e-01, -1.47340e-25,
                                         1.14879e-23, -5.57978e-34,  6.22567e-02,  2.25813e-02,
                                        -1.36077e-24,  8.36356e-25,  3.56089e-26 };

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

        std::cout << "SDVS:\n"; tardigradeVectorTools::print( SDVS );

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
