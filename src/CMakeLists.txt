add_library(ARC_Simulation_Reconstruction
	    ParticleTrack.cpp
	    Photon.cpp
	    RadiatorCell.cpp
	    SiPM.cpp
	    TrackingVolume.cpp)

target_include_directories(ARC_Simulation_Reconstruction PUBLIC ../include)

target_link_libraries(ARC_Simulation_Reconstruction PUBLIC ROOT::Physics ROOT::Tree ROOT::RooFit ROOT::Gpad)