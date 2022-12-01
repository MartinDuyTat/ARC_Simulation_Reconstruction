// Martin Duy Tat 29th April 2022
/**
 * ParticleTrack represents a charged particle with position, direction and momentum inside the detector volume
 */

#ifndef PARTICLETRACK
#define PARTICLETRACK

#include<vector>
#include<memory>
#include"TLine.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Photon.h"
#include"TrackingVolume.h"
#include"RadiatorArray.h"
#include"Particle.h"
#include"HelixPath.h"

using Vector = ROOT::Math::XYZVector;

class ParticleTrack: public Particle {
 public:
  /**
   * Construct a charged particle with momentum and ID at the interaction point
   * @param ParticleID PDG particle ID convention
   * @param Momentum Particle momentum, in GeV
   * @param Position Particle initial position, by default it's the origin
   * @param BField Magnetic field that particle passes through
   */
  ParticleTrack(int ParticleID,
		const Vector &Momentum,
		std::size_t TrackNumber,
		double BField);
  /**
   * Need virtual constructor since it's a virtual class
   */
  ~ParticleTrack() = default;
  /**
   * Enum class with the location of the particle
   * TrackerVolume: Inside the tracker volume, where the IP is
   * EntranceWindow: Particle is at the entrance window to the radiator cell
   * MissedEntranceWindow: Particle didn't hit this radiator entrance window
   * MissedRadiator: Particle did not enter any radiator cells
   * Radiator: Inside the radiatorCell
   * Mirror: The particle has reached the mirror
   * MissedMirror: The particle missed the mirror
   */
  enum class Location{
    TrackerVolume,
    EntranceWindow,
    MissedEntranceWindow,
    Radiator,
    Mirror,
    MissedMirror};
  /**
   * Track particle through inner tracker with magnetic field
   * @param InnerTracker The tracking volume
   * @return Returns true if tracking through the magnetic field was successful
   */
  bool TrackThroughTracker(const TrackingVolume &InnerTracker);
  /**
   * Set radiator
   */
  void SetRadiator(const RadiatorCell *radiatorCell);
  /**
   * Convert to local radiator coordinates
   */
  virtual void ConvertToRadiatorCoordinates() override;
  /**
   * Convert back to global coordinates
   */
  virtual void ConvertBackToGlobalCoordinates() override;
  /**
   * Track particle through radiator cell
   * @return Returns true if tracking through the magnetic field was successful
   */
  bool TrackThroughAerogel();
  /**
   * Helper function to track particle through gas until it hits the correct mirror
   */
  void TrackThroughGasToMirror();
  /**
   * If particle misses the mirror in this cell, call this function to migrate to next cell
   */
  bool TrackToNextCell(const RadiatorArray &radiatorArray);
  /**
   * Generate Cherenkov photon from aerogel
   */
  Photon GeneratePhotonFromAerogel() const;
  /**
   * Generate Cherenkov photon from gas
   */
  Photon GeneratePhotonFromGas() const;
  /**
   * Generate Cherenkov photons from aerogel according to Frank-Tamm relation
   */
  std::vector<Photon> GeneratePhotonsFromAerogel() const;
  /**
   * Generate Cherenkov photons from gas according to Frank-Tamm relation
   */
  std::vector<Photon> GeneratePhotonsFromGas() const;
  /**
   * Generate Cherenkov photon
   * @param Entry_s Path length at entry of radiator
   * @param Exit_s Path length at exit of radiator
   * @param Radiator The medium the photon was emitted in, to determine the index of refraction
   */
  Photon GeneratePhoton(double Entry,
			double Exit,
			Photon::Radiator Radiator) const;
  /**
   * Get the particle speed, in units of c
   */
  double Beta() const;
  /**
   * Draw particle track and the photons
   */
  std::unique_ptr<TLine> DrawParticleTrack() const;
  /**
   * Get particle location
   */
  Location GetParticleLocation() const;
  /**
   * Get the position of the entrance window
   */
  Vector GetEntranceWindowPosition() const;
  /**
   * Rotate around the phi direction to map particle to a valid radiator cell
   * @param DeltaPhi Angle that is rotated in phi
   */
  virtual void MapPhi(double DeltaPhi) override;
  /**
   * Reflect everything in the z-direction since radiator cells are symmetric
   */
  virtual void ReflectZ() override;
  /**
   * Reflect everything in the y-direction
   */
  virtual void ReflectY() override;
  /**
   * Check if particle is at the radiator so that we can search for the radiator cell
   */
  virtual bool IsAtRadiator() const override;
 private:
  /**
   * Particle momentum, in GeV
   */
  ARCVector m_Momentum;
  /**
   * Initial particle position in the global coordinate system, in m
   */
  ARCVector m_InitialPosition;
  /**
   * Entry path length of aerogel
   */
  double m_AerogelEntry_s;
  /**
   * Exit path length of aerogel
   */
  double m_AerogelExit_s;
  /**
   * Entry path length of gas
   */
  double m_GasEntry_s;
  /**
   * Exit path length of gas
   */
  double m_GasExit_s;
  /**
   * Entrance window position
   */
  ARCVector m_EntranceWindowPosition;
  /**
   * Particle ID
   */
  int m_ParticleID;
  /**
   * Flag that keeps track of where the particle is
   */
  Location m_Location;
  /**
   * Flag that is true if emission point of photon is random
   */
  bool m_RandomEmissionPoint;
  /**
   * Flag that is true if chromatic dispersion is on
   */
  bool m_ChromaticDispersion;
  /**
   * Mass of particle
   */
  const double m_Mass;
  /**
   * Frank-Tamm relation for photon yield
   */
  double GetPhotonYield(double x, double Beta, double n) const;
  /**
   * The multiplication factor in the photon yield calculation
   */
  const double m_PhotonMultiplier;
  /**
   * Track number
   */
  std::size_t m_TrackNumber;
  /**
   * List of track numbers that we want to draw
   */
  std::vector<std::size_t> m_TracksToDraw;
  /**
   * Helix trajectory
   */
  HelixPath m_Helix;
  /**
   * Current path length
   */
  double m_PathLength;
  /**
   * Helper function that checks if this track should be drawn
   */
  bool IsTrackDrawn() const;
  /**
   * Helper function that determines the photon multiplier
   */
  static double GetPhotonMultiplier(const Vector &Momentum);
};

#endif
