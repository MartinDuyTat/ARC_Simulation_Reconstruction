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

using Vector = ROOT::Math::XYZVector;

class ParticleTrack {
 public:
  /**
   * Construct a charged particle with momentum and ID at the interaction point
   * @param ParticleID PDG particle ID convention
   * @param Momentum Particle momentum, in GeV
   * @param Position Particle initial position, by default it's the origin
   */
  ParticleTrack(int ParticleID,
		const Vector &Momentum,
		const Vector &Position = Vector(0.0, 0.0, 0.0));
  /**
   * Enum with the two coordinate systems used
   * GlobalDetector: z along symmetry axis, x is up, y is out of the plane (only used to generate charged tracks)
   * LocalRadiator: z axis in the radial direction, x in the theta direction and y in the phi direction
   */
  enum class CoordinateSystem{GlobalDetector, LocalRadiator};
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
    MissedRadiator,
    Radiator,
    Mirror,
    MissedMirror};
  /**
   * Track particle through inner tracker with magnetic field
   */
  void TrackThroughTracker(const TrackingVolume &InnerTracker);
  /**
   * Find the correct radiator that the particle goes through
   */
  bool FindRadiator(RadiatorArray &radiatorArray);
  /**
   * Set radiator
   */
  void SetRadiator(const RadiatorCell *radiatorCell);
  /**
   * Convert to local radiator coordinates
   */
  void ConvertToRadiatorCoordinates();
  /**
   * Track particle through radiator cell
   */
  void TrackThroughRadiatorCell();
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
   * Get index of fraction
   */
  double GetIndexRefraction(Photon::Radiator Radiator, double Energy) const;
  /**
   * Generate Cherenkov photon
   * @param Entry point of radiator
   * @param Exit point of ratiator
   * @param Radiator The medium the photon was emitted in, to determine the index of refraction
   */
  Photon GeneratePhoton(const Vector &Entry, const Vector &Exit, Photon::Radiator Radiator) const;
  /**
   * Get the particle speed, in units of c
   */
  double Beta() const;
  /**
   * Get momentum vector
   */
  const Vector& GetMomentum() const;
  /**
   * Get entry point of particle in radiator
   */
  const Vector& GetEntryPoint(Photon::Radiator Radiator) const;
  /**
   * Get exit point of particle in radiator
   */
  const Vector& GetExitPoint(Photon::Radiator Radiator) const;
  /**
   * Draw particle track
   */
  std::unique_ptr<TLine> DrawParticleTrack() const;
  /**
   * Get particle position
   */
  const Vector& GetPosition() const;
  /**
   * Get particle location
   */
  Location GetParticleLocation() const;
  /**
   * Get the position of the entrance window
   */
  const Vector GetEntranceWindowPosition() const;
  /**
   * Rotate around the phi direction to map particle to a valid radiator cell
   * @param DeltaPhi Angle that is rotated in phi
   */
  void MapPhi(double DeltaPhi);
  /**
   * Reflect everything in the z-direction since radiator cells are symmetric
   */
  void ReflectZ();
  /**
   * Reflect everything in the y-direction
   */
  void ReflectY();
  /**
   * Get the column number of the radiator
   */
  double GetRadiatorColumnNumber() const;
  /**
   * Get the row number of the radiator
   */
  double GetRadiatorRowNumber() const;
 private:
  /**
   * Particle momentum, in GeV
   */
  Vector m_Momentum;
  /**
   * Particle position, in m
   */
  Vector m_Position;
  /**
   * Initial particle position in the global coordinate system, in m
   */
  Vector m_InitialPosition;
  /**
   * Particle ID
   */
  int m_ParticleID;
  /**
   * Pointer to the radiator cell that track enters
   */
  const RadiatorCell *m_RadiatorCell;
  /**
   * Flag that keeps track of where the particle is
   */
  Location m_Location;
  /**
   * Entry point of aerogel
   */
  Vector m_AerogelEntry;
  /**
   * Exit point of aerogel
   */
  Vector m_AerogelExit;
  /**
   * Entry point of gas
   */
  Vector m_GasEntry;
  /**
   * Exit point of gas
   */
  Vector m_GasExit;
  /**
   * Entrance window position
   */
  Vector m_EntranceWindowPosition;
  /**
   * Which coordinate system the momentum and position is given in
   */
  CoordinateSystem m_CoordinateSystem;
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
   * Helper function to swap x and z directions (so that z points up in the dell
   */
  void SwapXZ(Vector &Vec) const;
  /**
   * Helper function to track particle through gas until it hits the correct mirror
   * If it ends up outside the cell in the theta direction, move to the next radiator cell
   */
  void TrackThroughGasToMirror();
  /**
   * Helper function to swap radiator cell if particle is outside in the theta direction
   * @return Returns true if swapping radiator cell was successful (no edges hit)
   */
  //bool SwapRadiatorCell();
  /**
   * Helper function to change coordinate origin
   */
  void ChangeCoordinateOrigin(const Vector &Shift);
};

#endif
