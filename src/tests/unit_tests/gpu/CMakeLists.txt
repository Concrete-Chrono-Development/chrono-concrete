# ------------------------------------------------------------------------------
# Additional include paths and libraries
# ------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${CH_GPU_INCLUDES})

SET(LIBRARIES
    Chrono_core
    Chrono_gpu
)

# ------------------------------------------------------------------------------
# List of all executables
# ------------------------------------------------------------------------------

SET(TESTS
    #utest_GPU_mini
    utest_GPU_frictionrolling
    utest_GPU_meshrolling
    utest_GPU_ballistic
    utest_GPU_stack
    utest_GPU_pyramid
)

# A hack to set the working directory in which to execute the CTest
# runs.  This is needed for tests that need to access the Chrono data
# directory (since we use a relative path to it)

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(MY_WORKING_DIR "${EXECUTABLE_OUTPUT_PATH}/Release")
else()
  set(MY_WORKING_DIR ${EXECUTABLE_OUTPUT_PATH})
endif()


# ------------------------------------------------------------------------------
# Add all executables
# ------------------------------------------------------------------------------

MESSAGE(STATUS "Add unit test programs for GPU module")

FOREACH(PROGRAM ${TESTS})
    # MESSAGE(STATUS "...add ${PROGRAM}")

    set(FILES ${PROGRAM}.cpp ut_gpu_utils.h)
    ADD_EXECUTABLE(${PROGRAM} ${FILES})
    SOURCE_GROUP("" FILES ${FILES})

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES FOLDER tests)
    SET_PROPERTY(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES} gtest_main)
    ADD_TEST(${PROGRAM} ${PROJECT_BINARY_DIR}/bin/${PROGRAM})

    INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})

    if(MSVC)
      set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
    endif()
    
    set_tests_properties(${PROGRAM} PROPERTIES WORKING_DIRECTORY ${MY_WORKING_DIR})
    
ENDFOREACH(PROGRAM)
