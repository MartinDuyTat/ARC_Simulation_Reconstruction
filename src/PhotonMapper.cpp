// Martin Duy Tat 1st May 2022

#include<algorithm>
#include<memory>
#include"TMath.h"
#include"TRandom.h"
#include"PhotonMapper.h"
#include"Photon.h"
#include"Settings.h"
#include"SiPM.h"

namespace PhotonMapper {

  double PhotonMirrorDistance(const Photon &photon) {
    auto MirrorCentre = photon.m_RadiatorCell->GetMirrorCentre();
    const double R = photon.m_RadiatorCell->GetMirrorCurvature();
    const Vector MirrorCentreMinusPosition = MirrorCentre - photon.m_Position;
    const double b = -photon.m_Direction.Dot(MirrorCentreMinusPosition);
    const double c = (MirrorCentreMinusPosition).Mag2() - R*R;
    const double Discriminant = TMath::Sqrt(b*b - c);
    const double s1 = -b + Discriminant;
    const double s2 = -b - Discriminant;
    return std::max(s1, s2);
  }

  void TracePhotonToMirror(Photon &photon) {
    const double s = PhotonMirrorDistance(photon);
    const Vector Mirror = photon.m_Position + s*photon.m_Direction;
    auto MirrorCentre = photon.m_RadiatorCell->GetMirrorCentre();
    const auto Curvature = photon.m_RadiatorCell->GetMirrorCurvature();
    const Vector Normal = (Mirror - MirrorCentre)/Curvature;
    const Vector Vout = photon.m_Direction - 2*photon.m_Direction.Dot(Normal)*Normal;
    photon.m_Position = Mirror;
    photon.m_Direction = Vout;
    const double MaxHeight = photon.m_RadiatorCell->GetRadiatorThickness()
                           - 2*photon.m_RadiatorCell->GetVesselThickness()
                           - photon.m_RadiatorCell->GetCoolingThickness();
    if(photon.m_Position.Z() > MaxHeight) {
      photon.m_Status = Photon::Status::WallMiss;
    } else if(s < 0.0) {
      photon.m_Status = Photon::Status::Backwards;
    } else if((photon.m_Position - MirrorCentre).Mag2() - Curvature*Curvature > 1e-6) {
      photon.m_Status = Photon::Status::MirrorMiss;
    } else if(photon.m_RadiatorCell->IsInsideCell(photon)) {
      photon.m_Status = Photon::Status::MirrorHit;
      photon.m_MirrorHitPosition = std::make_unique<Vector>(photon.m_Position);
    } else {
      photon.m_Status = Photon::Status::MirrorMiss;
    }
  }

  PhotonHit TracePhotonToDetector(Photon &photon) {
    const auto Position = photon.m_Position;
    const auto Direction = photon.m_Direction;
    auto Detector = photon.m_RadiatorCell->GetDetector();
    const double SinTheta = TMath::Sin(Detector.GetDetectorTilt());
    const double CosTheta = TMath::Cos(Detector.GetDetectorTilt());
    // Find intersection between straight photon trajectory and tilted detector plane
    const double Det = SinTheta*Direction.X() - CosTheta*Direction.Z();
    const double s = ((Position.Z() - Detector.GetDetectorZPosition())*CosTheta
                    - (Position.X() - Detector.GetDetectorXPosition())*SinTheta)/Det;
    photon.m_Position += s*Direction;
    // Check if hit position is inside aerogel
    const double AerogelThickness = photon.m_RadiatorCell->GetAerogelThickness();    
    if(photon.m_Position.Z() < AerogelThickness) {
      const double Slope = TMath::Abs(Direction.Z());
      photon.m_AerogelTravelDistance += AerogelThickness/Slope;
    }
    // TODO: Account for photons generated in aerogel
    if(IsScatteredInAerogel(photon)) {
      photon.m_Status = Photon::Status::AerogelScattered;
    }
    const PhotonHit photonHit = photon.m_RadiatorCell->GetDetector().AddPhotonHit(photon);
    return photonHit;
  }

  std::unique_ptr<PhotonHit> TracePhoton(Photon &photon) {
    TracePhotonToMirror(photon);
    if(photon.m_MirrorHitPosition) {
      PhotonHit photonHit = TracePhotonToDetector(photon);
      return std::make_unique<PhotonHit>(photonHit);
    } else {
      return nullptr;
    }
  }

  bool IsScatteredInAerogel(const Photon &photon) {
    constexpr double T0 = 0.9679;
    constexpr double Clarity = 5.087e10;
    const double Lambda = 1239.8/photon.m_Energy;
    const double Lambda2 = Lambda*Lambda;
    const double Lambda4 = Lambda2*Lambda2;
    const double Exponent = -Clarity*(photon.m_AerogelTravelDistance)/Lambda4;
    const double Transmission = T0*TMath::Exp(Exponent);
    return gRandom->Uniform(0.0, 1.0) > Transmission;
  }

}
