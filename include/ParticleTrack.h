// Martin Duy Tat 29th April 2022
/**
 * ParticleTrack represents a charged particle with position, direction and momentum inside the detector volume
 */

#ifndef PARTICLETRACK
#define PARTICLETRACK

#include<vector>
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Photon.h"
#include"TrackingVolume.h"
#include"RadiatorCell.h"

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
   */
  enum class CoordinateSystem{GlobalDetector, LocalRadiator};
  /**
   * Track particle through inner tracker with magnetic field
   */
  void TrackThroughTracker(const TrackingVolume &InnerTracker);
  /**
   * Convert to local radiator coordinates
   */
  void ConvertToRadiatorCoordinates(const RadiatorCell &Cell);
  /**
   * Track particle through radiator cell
   */
  void TrackThroughRadiatorCell(const RadiatorCell &Cell);
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
   * @param Entry point of radiator
   * @param Exit point of ratiator
   * @param n_phase Index of refraction for phase velocity
   */
  Photon GeneratePhoton(const Vector &Entry, const Vector &Exit, double n_phase) const;
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
   * Particle ID
   */
  int m_ParticleID;
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
