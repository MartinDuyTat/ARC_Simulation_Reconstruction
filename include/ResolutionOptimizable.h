// Martin Duy Tat 28th July 2022
/**
 * ResolutionOptimizable is a derived class from the interface IOptimizable
 * It takes care of all calculations of the resolution for a given radiator cell
 */

#ifndef RESOLUTIONOPTIMIZABLE
#define RESOLUTIONOPTIMIZABLE

#include<vector>
#include<string_view>
#include"RadiatorCell.h"
#include"Photon.h"
#include"DifferentialEvolution.h"
#include"ParticleTrack.h"
#include"RadiatorArray.h"

using Tracks = std::vector<ParticleTrack>;

class ResolutionOptimizable: public de::IOptimizable {
 public:
  /**
   * Constructor that saves the radiator configuration, tracks and photons
   * @param radiatorCell The radiator configuration
   * @param radiatorArray The full radiator array
   * @param ParticlesPhotons The generated particles and photons
   */
  ResolutionOptimizable(RadiatorCell &radiatorCell,
			const RadiatorArray &radiatorArray,
			const Tracks &Particles);
  /**
   * Function that evaluates the resolution
   * The free parameters are:
   * Mirror curvature
   * Shift in x position
   * Shift in z position
   * @param x Free parameters
   */
  double EvaluateCost(std::vector<double> x) const override;
  /**
   * Find the number of free parameters
   */
  std::size_t NumberOfParameters() const override;
  /**
   * Set up the parameter space
   */
  std::vector<Constraints> GetConstraints() const override;
  /**
   * Get the fixed parameters
   */
  std::unordered_map<std::size_t, double> GetFixedParameters() const;
 private:
  /**
   * The radiator cell we want to optimise
   */
  RadiatorCell *m_RadiatorCell;
  /**
   * The full radiator array
   */
  const RadiatorArray *m_RadiatorArray;
  /**
   * The generated tracks and photons used to calculate resolution
   */
  const Tracks *m_Particles;
  /**
   * Map of parameters and that are fixed and their values
   */
  std::unordered_map<std::size_t, double> m_FixedParameters;
  /**
   * The parameter space
   */
  std::vector<Constraints> m_ParameterSpace;
  /**
   * The seed used for each iteration
   */
  const std::size_t m_Seed;
  /**
   * Free parameter names
   */
  static constexpr std::array<std::string_view, 5> m_Parameters{
    "MirrorCurvature",
    "MirrorXPosition",
    "MirrorZPosition",
    "DetectorPosition",
    "DetectorTilt"};
};

#endif
