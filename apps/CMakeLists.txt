add_executable(RunARC RunARC.cpp)

target_link_libraries(RunARC PUBLIC ARC_Simulation_Reconstruction)
target_link_libraries(RunARC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

add_executable(OptimiseARC OptimiseARC.cpp)

target_link_libraries(OptimiseARC PUBLIC ARC_Simulation_Reconstruction)
target_link_libraries(OptimiseARC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree ROOT::Minuit2 OpenMP::OpenMP_CXX)

install(TARGETS OptimiseARC RunARC DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
