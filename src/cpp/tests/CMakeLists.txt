#Set the tests for the balance equations
set(TEST_NAME "test_${BALANCE_EQUATION_LIBRARY}")
add_executable(${TEST_NAME} "test_${BALANCE_EQUATION_LIBRARY_FILENAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
target_compile_options(${TEST_NAME} PRIVATE)
target_link_libraries(${TEST_NAME} PUBLIC ${BALANCE_EQUATION_LIBRARY})

get_target_property(target_sources "test_${BALANCE_EQUATION_LIBRARY}" SOURCES)
foreach(source ${target_sources})
    set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
endforeach(source)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${TEST_NAME}
        PUBLIC
            ${Boost_INCLUDE_DIRS} # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
            ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
            ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
    )
endif()

#Set the tests for the user subroutine interfaces
foreach(user_subroutine ${USER_SUBROUTINES})
    set(TEST_NAME "test_${user_subroutine}_interface")
    add_executable(${TEST_NAME} "${TEST_NAME}.cpp")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    target_compile_options(${TEST_NAME} PRIVATE)
    target_link_libraries(${TEST_NAME} PUBLIC ${MATERIAL_MODEL_LIBRARY} ${user_subroutine})

    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(
            ${TEST_NAME}
            PUBLIC
                ${Boost_INCLUDE_DIRS} # Required for MacOSX CMake builds. Not sure if it's a CMake or clang issue.
                ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_constitutive_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_micromorphic_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_solver_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${user_subroutine_SOURCE_DIR}/${CPP_SRC_PATH}
        )
    endif()
    get_target_property(target_sources ${TEST_NAME} SOURCES)
    foreach(source ${target_sources})
        set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    endforeach(source)
endforeach(user_subroutine)
