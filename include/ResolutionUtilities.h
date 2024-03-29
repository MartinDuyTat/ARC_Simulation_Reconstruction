// Martin Duy Tat 28th July 2022
/**
 * The ResolutionUtilites namespace contains the function to calculate the resolution from some generated tracks and their photons
 * OpenMP is used to parallelise the calculation
 */

#ifndef RESOLUTIONUTILITIES
#define RESOLUTIONUTILITIES

#include<vector>
#include<utility>
#include"Photon.h"
#include"RadiatorCell.h"
#include"ParticleTrack.h"
#include"RadiatorArray.h"

using Tracks = std::vector<ParticleTrack>;

namespace ResolutionUtilities {
  /**
   * Helper function where the (parallel) calculation takes place
   */
  double CalculateResolution(const Tracks &Particles,
			     const RadiatorCell *radiatorCell,
			     const RadiatorArray &radiatorArray,
			     bool IncludeCentrePenalty);
  /**
   * The "cost" function for minimisation
   */
  double fcn(double MirrorCurvature,
	     double MirrorXPosition,
	     double MirrorZPosition,
	     double DetectorPosition,
	     double DetectorTilt,
	     RadiatorCell &radiatorCell,
	     const RadiatorArray &radiatorArray,
	     const Tracks &Particles,
	     std::size_t Seed,
	     bool IncludeCentrePenalty);
  /**
   * Load fit results and plot the fit projections
   */
  void PlotProjections(RadiatorCell &radiatorCell,
		       const RadiatorArray &radiatorArray,
		       const Tracks &Particles);
  /**
   * Do the actual fit and save the results
   */
  void DoFit(RadiatorCell &radiatorCell,
	     const RadiatorArray &radiatorArray,
	     const Tracks &Particles,
	     std::size_t Column,
	     std::size_t Row);

}

#endif
