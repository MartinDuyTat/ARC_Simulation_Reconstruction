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
   * @param Momentum Particle momentum, in GeV
   * @param ParticleID PDG particle ID convention
   */
  ParticleTrack(const Vector &Momentum, int ParticleID);
  /**
   * Enum with the two coordinate systems used
   * GlobalDetector: z along symmetry axis, x is up, y is out of the plane (only used to generate charged tracks)
   * LocalRadiator: z axis in the radial direction, x in the theta direction and y in the phi direction
   */
  enum class CoordinateSystem{GlobalDetector, LocalRadiator};
  /**
   * Track particle through inner tracker with magnetic field
   */
  void TrackThroughTracker(const TrackingVolume &InnerTracker);
  /**
   * Convert to local radiator coordinates
   */
  void ConvertToRadiatorCoordinates(RadiatorArray &Cell);
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
   * Get photon hits in SiPM
   */
  const std::vector<PhotonHit>& GetPhotonHits() const;
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
   * Initial particle position, in m
   */
  const Vector m_InitialPosition;
  /**
   * Particle ID
   */
  int m_ParticleID;
  /**
   * Pointer to the radiator cell that track enters
   */
  RadiatorCell *m_RadiatorCell;
  /**
   * Flag that is true when track has been traced through inner tracker detector
   */
  bool m_TrackedThroughTracker;
  /**
   * Flag that is true when track has been traced through radiator cell
   */
  bool m_TrackedThroughRadiator;
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
   * Which coordinate system the momentum and position is given in
   */
  CoordinateSystem m_CoordinateSystem;
  /**
   * Frank-Tamm relation for photon yield
   */
  double GetPhotonYield(double x, double Beta, double n) const;
  /**
   * Helper function that maps the phi direction to the cell near phi = 0
   */
  void MapPhiBack(Vector &Vec) const;
  /**
   * Helper function to rotate around y axis in coordinate transformation
   */
  void RotateY(Vector &Vec) const;
};

#endif
