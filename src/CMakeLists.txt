add_library(ARC_Simulation_Reconstruction
	    ParticleTrack.cpp
	    ParticleMass.cpp
	    Photon.cpp
	    PhotonMapper.cpp
	    PhotonReconstructor.cpp
	    RadiatorCell.cpp
	    ReconstructedPhoton.cpp
	    SiPM.cpp
	    TrackingVolume.cpp)

target_include_directories(ARC_Simulation_Reconstruction PUBLIC ../include)

target_link_libraries(ARC_Simulation_Reconstruction PUBLIC ROOT::Physics ROOT::Tree ROOT::Gpad ROOT::MathMore)
