// Martin Duy Tat 3rd March 2022

#include<algorithm>
#include"PhotonReconstructor.h"
#include"Photon.h"
#include"ParticleTrack.h"
#include"RadiatorCell.h"
#include"ReconstructedPhoton.h"
#include"TMath.h"
#include"Math/Vector3Dfwd.h"
#include"Math/DisplacementVector3D.h"
#include"Math/Polynomial.h"

using Vector = ROOT::Math::XYZVector;
using Polynomial = ROOT::Math::Polynomial;

namespace PhotonReconstructor {

  double ReconstructCherenkovAngle(const ParticleTrack &Particle,
				   const PhotonHit &photonHit,
				   bool TrueEmissionPoint,
				   Photon::Radiator Radiator) {
    const RadiatorCell *radiatorCell = photonHit.m_Photon->m_RadiatorCell;
    auto MirrorCentre = radiatorCell->GetMirrorCentre();
    auto GetEmissionPoint = [=]() {
      if(TrueEmissionPoint) {
	return photonHit.m_Photon->m_EmissionPoint;
      } else {
	return (Particle.GetExitPoint(Radiator)
	      + Particle.GetEntryPoint(Radiator))*0.5;
      }
    };
    const Vector EmissionPoint = GetEmissionPoint() - MirrorCentre;
    const Vector DetectionPoint = Vector(photonHit.x, photonHit.y, 0.0)
                                - MirrorCentre;
    const double Curvature = radiatorCell->GetMirrorCurvature();
    // Distance between emission point and mirror centre,
    // in units of mirror curvature
    const double EmissionMirrorDist = TMath::Sqrt(EmissionPoint.Mag2())/Curvature;
    // Parallel component of distance between detection point and mirror centre,
    // in units of mirror curvature
    const double DetectionMirrorParaDist = EmissionPoint.Dot(DetectionPoint)
                                           /(EmissionMirrorDist*Curvature*Curvature);
    const double DetectionMirrorPerpDist = TMath::Sqrt(
					   DetectionPoint.Mag2()/(Curvature*Curvature)
					 - TMath::Power(DetectionMirrorParaDist, 2));
    auto quarticSolution = SolveQuartic(EmissionMirrorDist,
					DetectionMirrorParaDist,
					DetectionMirrorPerpDist);
    if(quarticSolution.m_DegenerateSolution) {
      return -1.0;
    }
    Vector ReflectionPoint = GetReflectionPoint(EmissionPoint,
						DetectionPoint,
						DetectionMirrorParaDist,
						DetectionMirrorPerpDist,
						EmissionMirrorDist,
						quarticSolution, 0);
    Vector OtherReflectionPoint = GetReflectionPoint(EmissionPoint,
						     DetectionPoint,
						     DetectionMirrorParaDist,
						     DetectionMirrorPerpDist,
						     EmissionMirrorDist,
						     quarticSolution, 1);
    if((OtherReflectionPoint - Particle.GetPosition()).Mag2()
       < (ReflectionPoint - Particle.GetPosition()).Mag2()) {
      std::swap(ReflectionPoint, OtherReflectionPoint);
    }
    const Vector ReflectionEmissionPoint = ReflectionPoint - EmissionPoint;
    const Vector Direction = Particle.GetMomentum().Unit();
    const double CosTheta = ReflectionEmissionPoint.Unit().Dot(Direction);
    return TMath::ACos(CosTheta);
  }
  
  ReconstructedPhoton ReconstructPhoton(const ParticleTrack &Particle,
					const PhotonHit &photonHit,
					Photon::Radiator Radiator) {
    ReconstructedPhoton reconstructedPhoton(*photonHit.m_Photon);
    reconstructedPhoton.m_CherenkovAngle_TrueEmissionPoint =
      ReconstructCherenkovAngle(Particle, photonHit, true, Radiator);
    reconstructedPhoton.m_CherenkovAngle =
      ReconstructCherenkovAngle(Particle, photonHit, false, Radiator);
    return reconstructedPhoton;
  }

  QuarticSolution SolveQuartic(double EmMirrorDist,
			       double DetMirrorParaDist,
			       double DetMirrorPerpDist) {
    // Em = Emission
    // Det = Detection
    const double DetMirrorDist2 = DetMirrorParaDist*DetMirrorParaDist
                                + DetMirrorPerpDist*DetMirrorPerpDist;
    const double Denominator = 4.0*EmMirrorDist*EmMirrorDist*DetMirrorDist2;
    const double a = -4.0*EmMirrorDist*EmMirrorDist*DetMirrorPerpDist/Denominator;
    const double b = (DetMirrorPerpDist*DetMirrorPerpDist
		    + TMath::Power(EmMirrorDist + DetMirrorParaDist, 2))/Denominator
                    - 1.0;
    const double c = 2.0*EmMirrorDist*DetMirrorPerpDist
                     *(EmMirrorDist - DetMirrorParaDist)/Denominator;
    const double d = DetMirrorPerpDist*DetMirrorPerpDist
                     *(EmMirrorDist*EmMirrorDist - 1.0)/Denominator;
    Polynomial polynomial(1.0, a, b, c, d);
    auto polySolutions = polynomial.FindRealRoots();
    QuarticSolution quarticSolution;
    if(polySolutions.size() != 2) {
      quarticSolution.m_DegenerateSolution = true;
    }
    for(std::size_t i = 0; i < 2; i++) {
      quarticSolution.m_SinBeta[i] = polySolutions[i];
      quarticSolution.m_CosBeta[i] = (EmMirrorDist + DetMirrorParaDist)*polySolutions[i];
      quarticSolution.m_CosBeta[i] += EmMirrorDist*DetMirrorPerpDist
	                              *(1.0 - 2.0*polySolutions[i]*polySolutions[i]);
      quarticSolution.m_CosBeta[i] /= DetMirrorPerpDist
	                            + 2.0*EmMirrorDist*DetMirrorParaDist*polySolutions[i];
    }
    return quarticSolution;
  }

  Vector GetReflectionPoint(const Vector &EmPoint,
			    const Vector &DetPoint,
			    double DetMirrorParaDist,
			    double DetMirrorPerpDist,
			    double EmMirrorDist,
			    QuarticSolution quarticSolution,
			    int SolutionNumber) {
    Vector Reflection = EmPoint*(quarticSolution.m_CosBeta[SolutionNumber]/EmMirrorDist);
    Reflection += (quarticSolution.m_SinBeta[SolutionNumber]/DetMirrorPerpDist)
                 *(DetPoint - EmPoint*(DetMirrorParaDist/EmMirrorDist));
    return Reflection;
  }

}
