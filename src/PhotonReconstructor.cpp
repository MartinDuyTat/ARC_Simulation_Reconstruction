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
	return (Particle.GetExitPoint(Radiator) + Particle.GetEntryPoint(Radiator))*0.5;
      }
    };
    const Vector EmissionPoint = GetEmissionPoint() - MirrorCentre;
    const Vector DetectionPoint = Vector(photonHit.x, photonHit.y, 0.0) - MirrorCentre;
    const double MirrorCurvature = radiatorCell->GetMirrorCurvature();
    // Distance between emission point and mirror centre, in units of mirror curvature
    const double EmissionMirrorDist = TMath::Sqrt(EmissionPoint.Mag2())/MirrorCurvature;
    // Parallel component of distance between detection point and mirror centre, in units of mirror curvature
    const double DetectionMirrorParaDist = EmissionPoint.Dot(DetectionPoint)/(EmissionMirrorDist*MirrorCurvature*MirrorCurvature);
    const double DetectionMirrorPerpDist = TMath::Sqrt(DetectionPoint.Mag2()/TMath::Power(MirrorCurvature, 2) - DetectionMirrorParaDist*DetectionMirrorParaDist);
    auto quarticSolution = SolveQuartic(EmissionMirrorDist, DetectionMirrorParaDist, DetectionMirrorPerpDist);
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
    if((OtherReflectionPoint - Particle.GetPosition()).Mag2() < (ReflectionPoint - Particle.GetPosition()).Mag2()) {
      std::swap(ReflectionPoint, OtherReflectionPoint);
    }
    const Vector ReflectionEmissionPoint = ReflectionPoint - EmissionPoint;
    const double CosTheta = ReflectionEmissionPoint.Unit().Dot(Particle.GetMomentum().Unit());
    return TMath::ACos(CosTheta);
  }
  
  ReconstructedPhoton ReconstructPhoton(const ParticleTrack &Particle,
					const PhotonHit &photonHit,
					Photon::Radiator Radiator) {
    ReconstructedPhoton reconstructedPhoton(*photonHit.m_Photon);
    reconstructedPhoton.m_CherenkovAngle_TrueEmissionPoint = ReconstructCherenkovAngle(Particle, photonHit, true, Radiator);
    reconstructedPhoton.m_CherenkovAngle = ReconstructCherenkovAngle(Particle, photonHit, false, Radiator);
    return reconstructedPhoton;
  }

  QuarticSolution SolveQuartic(double EmissionMirrorDist, double DetectionMirrorParaDist, double DetectionMirrorPerpDist) {
    const double DetectionMirrorDist2 = DetectionMirrorParaDist*DetectionMirrorParaDist + DetectionMirrorPerpDist*DetectionMirrorPerpDist;
    const double Denominator = 4.0*EmissionMirrorDist*EmissionMirrorDist*DetectionMirrorDist2;
    const double a = -4.0*EmissionMirrorDist*EmissionMirrorDist*DetectionMirrorPerpDist/Denominator;
    const double b = (DetectionMirrorPerpDist*DetectionMirrorPerpDist + TMath::Power(EmissionMirrorDist + DetectionMirrorParaDist, 2))/Denominator - 1.0;
    const double c = 2.0*EmissionMirrorDist*DetectionMirrorPerpDist*(EmissionMirrorDist - DetectionMirrorParaDist)/Denominator;
    const double d = DetectionMirrorPerpDist*DetectionMirrorPerpDist*(EmissionMirrorDist*EmissionMirrorDist - 1.0)/Denominator;
    Polynomial polynomial(1.0, a, b, c, d);
    auto polynomialSolutions = polynomial.FindRealRoots();
    QuarticSolution quarticSolution;
    if(polynomialSolutions.size() != 2) {
      quarticSolution.m_DegenerateSolution = true;
    }
    for(std::size_t i = 0; i < 2; i++) {
      quarticSolution.m_SinBeta[i] = polynomialSolutions[i];
      quarticSolution.m_CosBeta[i] = (EmissionMirrorDist + DetectionMirrorParaDist)*polynomialSolutions[i];
      quarticSolution.m_CosBeta[i] += EmissionMirrorDist*DetectionMirrorPerpDist*(1.0 - 2.0*polynomialSolutions[i]*polynomialSolutions[i]);
      quarticSolution.m_CosBeta[i] /= DetectionMirrorPerpDist + 2.0*EmissionMirrorDist*DetectionMirrorParaDist*polynomialSolutions[i];
    }
    return quarticSolution;
  }

  Vector GetReflectionPoint(const Vector &EmissionPoint,
			    const Vector &DetectionPoint,
			    double DetectionMirrorParaDist,
			    double DetectionMirrorPerpDist,
			    double EmissionMirrorDist,
			    QuarticSolution quarticSolution,
			    int SolutionNumber) {
    Vector ReflectionPoint = EmissionPoint*quarticSolution.m_CosBeta[SolutionNumber]/EmissionMirrorDist;
    ReflectionPoint += quarticSolution.m_SinBeta[SolutionNumber]*(DetectionPoint - EmissionPoint*DetectionMirrorParaDist/EmissionMirrorDist)/DetectionMirrorPerpDist;
    return ReflectionPoint;
  }

}
