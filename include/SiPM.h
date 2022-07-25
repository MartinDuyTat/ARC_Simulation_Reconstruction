// Martin Duy Tat 29th April 2022
/**
 * SiPM represents a single Silicon Photo-Multiplier, which is located on the detector plane of each radiator cell
 * All positions are in the local radiator coordinate system
 */

#ifndef SIPM
#define SIPM

#include<vector>
#include<string>
#include<memory>
#include"TObject.h"
#include"Math/Interpolator.h"
#include"Math/Vector3Dfwd.h"
#include"Photon.h"

using Vector = ROOT::Math::XYZVector;
using ROOT::Math::Interpolator;
namespace InterpolationType = ROOT::Math::Interpolation;

struct PhotonHit {
  /**
   * Constructor saving the photon hits
   */
  PhotonHit(double x, double y, const Photon *photon): x(x), y(y), m_Photon(photon) {}
  /**
   * Detector hit coordinates
   */
  const double x;
  const double y;
  /**
   * Pointer to the photon that caused this hit
   */
  const Photon *m_Photon;
};

class SiPM {
 public:
  /**
   * Constructor setting up the detector coordinates and detector size
   * @param xPosition Position in the x direction
   * @param yPosition Position in the y direction
   */
  SiPM(double xPosition, double yPosition);
  /**
   * Add a photon hit
   */
  void AddPhotonHit(Photon &photon);
  /**
   * Plot photon hits
   */
  void PlotHits(const std::string &Filename) const;
  /**
   * Get photon hits in SiPM
   */
  const std::vector<PhotonHit>& GetPhotonHits() const;
  /**
   * Draw SiPM
   */
  std::unique_ptr<TObject> DrawSiPM(const Vector &RadiatorPosition) const;
  /**
   * Check if photon hit the detector
   */
  bool IsDetectorHit(const Photon &photon) const;
  /**
   * Remove all saved photon hits to reset the detector
   */
  void ResetDetector();
 private:
  /**
   * Size of detector in x direction
   */
  const double m_DetectorSizeX;
  /**
   * Size of detector in y direction
   */
  const double m_DetectorSizeY;
  /**
   * Centre of detector in x direction
   */
  const double m_DetectorPositionX;
  /**
   * Centre of detector in y direction
   */
  const double m_DetectorPositionY;
  /**
   * Max PDE
   */
  const double m_MaxPDE;
  /**
   * Interpolator for PDE
   */
  std::unique_ptr<Interpolator> m_Interpolator;
  /**
   * Vector of photon hits
   */
  std::vector<PhotonHit> m_PhotonHits;
  /**
   * The wavelengths used to measure PDE in SiPM, in nm
   */
  static constexpr std::array<double, 16> m_Lambda{
    287.0, 299.0, 322.0, 340.0,
    367.0, 392.0, 402.0, 413.0,
    422.0, 437.0, 452.0, 467.0,
    503.0, 590.0, 700.0, 800.0};
  /**
   * The measured PDE
   */
  static constexpr std::array<double, 16> m_PDE{
    0.20, 0.41, 0.46, 0.47,
    0.52, 0.59, 0.60, 0.60,
    0.59, 0.57, 0.56, 0.53,
    0.51, 0.40, 0.26, 0.13};
};

#endif
