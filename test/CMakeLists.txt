# Add directories to the list of includes
include_directories (${phaistos_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/../src ${Boost_INCLUDE_DIRS}) 

# Output binaries to global bin dir
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/test)

# # add profiling flag
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")

# Create executable
add_executable(test_charmm test_charmm.cpp)

# Link executable with library
target_link_libraries(test_charmm libphaistos ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
