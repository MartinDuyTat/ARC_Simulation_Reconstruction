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

using TracksPhotons = std::vector<std::pair<ParticleTrack, std::vector<Photon>>>;

namespace ResolutionUtilities {
  /**
   * Helper function where the (parallel) calculation takes place
   */
  double CalculateResolution(const RadiatorCell &radiatorCell,
			     const TracksPhotons &ParticlesPhotons);
  /**
   * The "cost" function for minimisation
   */
  double fcn(double MirrorCurvature,
	     double MirrorXPosition,
	     double MirrorZPosition,
	     RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons);
  /**
   * Load fit results and plot the fit projections
   */
  void PlotProjections(RadiatorCell &radiatorCell,
		       const TracksPhotons &ParticlesPhotons);
  /**
   * Do the actual fit and save the results
   */
  void DoFit(RadiatorCell &radiatorCell,
	     const TracksPhotons &ParticlesPhotons);

}

#endif
