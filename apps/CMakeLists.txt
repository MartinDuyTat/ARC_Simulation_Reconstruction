add_executable(RunARC RunARC.cpp)

target_link_libraries(RunARC PUBLIC ARC_Simulation_Reconstruction)
target_link_libraries(RunARC PUBLIC ROOT::Physics ROOT::RIO ROOT::Tree)

install(TARGETS RunARC DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../bin)
