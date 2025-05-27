#Set the tests for the balance equations
set(TEST_NAME "test_${BALANCE_EQUATION_LIBRARY}")
add_executable(${TEST_NAME} "test_${BALANCE_EQUATION_LIBRARY_FILENAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
target_compile_options(${TEST_NAME} PRIVATE "-lrt")
target_link_libraries(${TEST_NAME} PUBLIC ${BALANCE_EQUATION_LIBRARY})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${TEST_NAME} PUBLIC
                               ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                               ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                               ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH})
endif()

#Set the tests for the user subroutine interfaces
foreach(user_subroutine ${USER_SUBROUTINES}) 
    set(TEST_NAME "test_${user_subroutine}_interface")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    target_compile_options(${TEST_NAME} PRIVATE "-lrt")
#    target_link_libraries(${TEST_NAME} PUBLIC ${MATERIAL_MODEL_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC ${MATERIAL_MODEL_LIBRARY})
    target_link_libraries(${TEST_NAME} PUBLIC "$<LINK_LIBRARY:WHOLE_ARCHIVE,${user_subroutine}>")

    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(${TEST_NAME} PUBLIC
                                   ${Boost_INCLUDE_DIRS}  # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                                   ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${tardigrade_micromorphic_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${tardigrade_solver_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                                   ${user_subroutine_SOURCE_DIR}/${CPP_SRC_PATH})
    endif()
endforeach(user_subroutine)
