add_library(ARC_Simulation_Reconstruction
	    ARCVector.cpp
	    BarrelRadiatorArray.cpp
	    BarrelRadiatorCell.cpp
	    EndCapRadiatorArray.cpp
	    EndCapRadiatorCell.cpp
	    EventDisplay.cpp
	    HalfRadiatorCell.cpp
	    HelixFunctor.cpp
	    HelixPath.cpp
	    Quartic.cpp
	    Particle.cpp
	    ParticleTrack.cpp
	    ParticleMass.cpp
	    Photon.cpp
	    PhotonMapper.cpp
	    PhotonReconstructor.cpp
	    RadiatorArray.cpp
	    RadiatorCell.cpp
	    ReconstructedPhoton.cpp
	    ResolutionOptimizable.cpp
	    ResolutionUtilities.cpp
	    Settings.cpp
	    SiPM.cpp
	    TrackingVolume.cpp
	    Utilities.cpp)

target_include_directories(ARC_Simulation_Reconstruction PUBLIC ../include)

target_link_libraries(ARC_Simulation_Reconstruction PUBLIC ROOT::Physics ROOT::Tree ROOT::Gpad ROOT::MathMore ROOT::GenVector OpenMP::OpenMP_CXX)
