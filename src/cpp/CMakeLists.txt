# Project

# Build the micromorphic balance equations
add_library(
    ${BALANCE_EQUATION_LIBRARY}
    SHARED
    "${BALANCE_EQUATION_LIBRARY_FILENAME}.cpp"
    "${BALANCE_EQUATION_LIBRARY_FILENAME}.h"
)
set_target_properties(${BALANCE_EQUATION_LIBRARY} PROPERTIES PUBLIC_HEADER ${BALANCE_EQUATION_LIBRARY_FILENAME}.h)
target_link_libraries(${BALANCE_EQUATION_LIBRARY} tardigrade_error_tools)
target_compile_options(${BALANCE_EQUATION_LIBRARY} PUBLIC)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${BALANCE_EQUATION_LIBRARY}
        PUBLIC
            "${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
            "${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
    )
endif()

# Build the micromorphic material model library
set(LOCAL_BUILD_INCLUDES "")
foreach(user_subroutine ${USER_SUBROUTINES})
    add_library(
        "${user_subroutine}_interface"
        SHARED
        "${user_subroutine}_interface.h"
        "${user_subroutine}_interface.cpp"
    )
    set_target_properties("${user_subroutine}_interface" PROPERTIES PUBLIC_HEADER "${user_subroutine}_interface.h")
    target_compile_options("${user_subroutine}_interface" PUBLIC)
    target_link_libraries("${user_subroutine}_interface" tardigrade_error_tools ${user_subroutine})
    if(NOT cmake_build_type_lower STREQUAL "release")
        set(LOCAL_BUILD_INCLUDES ${LOCAL_BUILD_INCLUDES} "${${user_subroutine}_SOURCE_DIR}/${CPP_SRC_PATH}")
        target_include_directories(
            "${user_subroutine}_interface"
            PUBLIC "${${user_subroutine}_SOURCE_DIR}/${CPP_SRC_PATH}"
        )
    endif()
    get_target_property(target_sources ${user_subroutine} SOURCES)
    foreach(source ${target_sources})
        set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    endforeach(source)
endforeach(user_subroutine)

# Project
add_library(
    ${MATERIAL_MODEL_LIBRARY}
    SHARED
    "${MATERIAL_MODEL_LIBRARY_FILENAME}.cpp"
    "${MATERIAL_MODEL_LIBRARY_FILENAME}.h"
)
set_target_properties(
    ${MATERIAL_MODEL_LIBRARY}
    PROPERTIES PUBLIC_HEADER "${MATERIAL_MODEL_LIBRARY_FILENAME}.h" SUFFIX ".so"
)
target_compile_options(${MATERIAL_MODEL_LIBRARY} PUBLIC)
target_link_libraries(${MATERIAL_MODEL_LIBRARY} PUBLIC tardigrade_error_tools ${USER_INTERFACES})
foreach(user_interface ${USER_INTERFACES})
    target_link_options(
        ${MATERIAL_MODEL_LIBRARY}
        PUBLIC "SHELL:-Wl,--whole-archive $<TARGET_FILE:${user_interface}> -Wl,--no-whole-archive"
    )
endforeach(user_interface)

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(
        ${MATERIAL_MODEL_LIBRARY}
        PUBLIC
            "${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
            "${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
            ${LOCAL_BUILD_INCLUDES}
    )
endif()

install(
    TARGETS ${EXPORTED_LIBRARIES}
    EXPORT ${MATERIAL_MODEL_LIBRARY}_Targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

foreach(package ${BALANCE_EQUATION_LIBRARY} ${MATERIAL_MODEL_LIBRARY})
    install(
        TARGETS ${package}
        EXPORT ${package}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    get_target_property(target_sources ${package} SOURCES)
    foreach(source ${target_sources})
        set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
    endforeach(source)
endforeach(package)
