// Martin Duy Tat 29th April 2022
/**
 * ParticleTrack represents a charged particle with position, direction and momentum inside the detector volume
 */

#ifndef PARTICLETRACK
#define PARTICLETRACK

#include"Vector3Dfwd.h"
#include"Photon.h"
#include"TrackingVolume.h"

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
   * Generate Cherenkov photon
   * @param Entry point of radiator
   * @param Exit point of ratiator
   */
  Photon GeneratePhoton(const Vector &Entry, const Vector &Exit) const;
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
};

#endif
