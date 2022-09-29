// Martin Duy Tat 1st May 2022

#include<algorithm>
#include<memory>
#include"TMath.h"
#include"TRandom.h"
#include"PhotonMapper.h"
#include"Photon.h"
#include"Settings.h"
#include"SiPM.h"
#include"RadiatorArray.h"

namespace PhotonMapper {

  double PhotonMirrorDistance(const Photon &photon) {
    auto MirrorCentre = photon.GetRadiatorCell()->GetMirrorCentre();
    const double R = photon.GetRadiatorCell()->GetMirrorCurvature();
    const Vector MirrorCentreMinusPosition = MirrorCentre - photon.GetPosition();
    const double b = -photon.GetDirection().Dot(MirrorCentreMinusPosition);
    const double c = (MirrorCentreMinusPosition).Mag2() - R*R;
    if(b*b - c < 0.0) {
      return -100.0;
    }
    const double Discriminant = TMath::Sqrt(b*b - c);
    const double s1 = -b + Discriminant;
    const double s2 = -b - Discriminant;
    return std::max(s1, s2);
  }

  void TracePhotonToMirror(Photon &photon) {
    const double s = PhotonMirrorDistance(photon);
    photon.PropagatePhoton(s*photon.GetDirection());
    const double MaxHeight = photon.GetRadiatorCell()->GetRadiatorThickness()
                           - 2*photon.GetRadiatorCell()->GetVesselThickness()
                           - photon.GetRadiatorCell()->GetCoolingThickness();
    auto MirrorCentre = photon.GetRadiatorCell()->GetMirrorCentre();
    const auto Curvature = photon.GetRadiatorCell()->GetMirrorCurvature();
    if(photon.GetPosition().Z() > MaxHeight) {
      photon.UpdatePhotonStatus(Photon::Status::WallMiss);
    } else if(s < 0.0) {
      photon.UpdatePhotonStatus(Photon::Status::Backwards);
    } else if((photon.GetPosition() - MirrorCentre).Mag2() - Curvature*Curvature > 1e-6) {
      photon.UpdatePhotonStatus(Photon::Status::OutsideMirrorRadius);
    } else if(photon.GetRadiatorCell()->IsInsideCell(photon)) {
      photon.UpdatePhotonStatus(Photon::Status::MirrorHit);
      photon.RegisterMirrorHitPosition(photon.GetPosition());
    } else {
      photon.UpdatePhotonStatus(Photon::Status::MirrorMiss);
    }
  }

  PhotonHit TracePhotonToDetector(Photon &photon) {
    const Vector Mirror = photon.GetPosition();
    auto MirrorCentre = photon.GetRadiatorCell()->GetMirrorCentre();
    const auto Curvature = photon.GetRadiatorCell()->GetMirrorCurvature();
    const Vector Normal = (Mirror - MirrorCentre)/Curvature;
    photon.KickPhoton(-2*photon.GetDirection().Dot(Normal)*Normal);
    const auto Position = photon.GetPosition();
    const auto Direction = photon.GetDirection();
    auto Detector = photon.GetRadiatorCell()->GetDetector();
    const double Theta = Detector.GetDetectorTilt();
    const double CosTheta = TMath::Cos(Theta);
    const double SinTheta = (Theta < 0.0 ? -1.0 : +1.0)
                            *TMath::Sqrt(1.0 - CosTheta*CosTheta);
    // Find intersection between straight photon trajectory and tilted detector plane
    const double Det = SinTheta*Direction.X() - CosTheta*Direction.Z();
    const double s = ((Position.Z() - Detector.GetDetectorZPosition())*CosTheta
                    - (Position.X() - Detector.GetDetectorXPosition())*SinTheta)/Det;
    photon.PropagatePhoton(s*Direction);
    // Check if hit position is inside aerogel
    const double AerogelThickness = photon.GetRadiatorCell()->GetAerogelThickness();    
    if(photon.GetPosition().Z() < AerogelThickness) {
      const double Slope = TMath::Abs(Direction.Z());
      photon.AddAerogelTravelDistance(AerogelThickness/Slope);
    }
    // TODO: Account for photons generated in aerogel
    if(IsScatteredInAerogel(photon)) {
      photon.UpdatePhotonStatus(Photon::Status::AerogelScattered);
    }
    const PhotonHit photonHit = photon.GetRadiatorCell()->GetDetector().AddPhotonHit(photon);
    return photonHit;
  }

  std::unique_ptr<PhotonHit> TracePhoton(Photon &photon,
					 const RadiatorArray &radiatorArray) {
    TracePhotonToMirror(photon);
    std::size_t Counter = 0;
    while(photon.GetStatus() == Photon::Status::MirrorMiss && Counter < 10) {
      photon.UpdatePhotonStatus(Photon::Status::Emitted);
      photon.ConvertBackToGlobalCoordinates();
      if(photon.FindRadiator(radiatorArray)) {
	photon.ConvertToRadiatorCoordinates();
	photon.PutPhotonToEmissionPoint();
	TracePhotonToMirror(photon);
      } else {
	photon.UpdatePhotonStatus(Photon::Status::MirrorMiss);
	break;
      }
      Counter++;
    }
    if(photon.GetMirrorHitPosition()) {
      PhotonHit photonHit = TracePhotonToDetector(photon);
      return std::make_unique<PhotonHit>(photonHit);
    } else {
      return nullptr;
    }
  }

  bool IsScatteredInAerogel(const Photon &photon) {
    constexpr double T0 = 0.9679;
    constexpr double Clarity = 5.087e10;
    const double Lambda = 1239.8/photon.GetEnergy();
    const double Lambda2 = Lambda*Lambda;
    const double Lambda4 = Lambda2*Lambda2;
    const double Exponent = -Clarity*(photon.GetAerogelTravelDistance())/Lambda4;
    const double Transmission = T0*TMath::Exp(Exponent);
    return gRandom->Uniform(0.0, 1.0) > Transmission;
  }

}
