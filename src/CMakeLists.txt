add_library(ARC_Simulation_Reconstruction
	    EventDisplay.cpp
	    ParticleTrack.cpp
	    ParticleMass.cpp
	    Photon.cpp
	    PhotonMapper.cpp
	    PhotonReconstructor.cpp
	    RadiatorArray.cpp
	    RadiatorCell.cpp
	    ReconstructedPhoton.cpp
	    Settings.cpp
	    SiPM.cpp
	    TrackingVolume.cpp)

target_include_directories(ARC_Simulation_Reconstruction PUBLIC ../include)

target_link_libraries(ARC_Simulation_Reconstruction PUBLIC ROOT::Physics ROOT::Tree ROOT::Gpad ROOT::MathMore)
