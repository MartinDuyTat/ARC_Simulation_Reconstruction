// Martin Duy Tat 1st May 2022

#include<algorithm>
#include"TMath.h"
#include"PhotonMapper.h"
#include"Photon.h"
#include"Settings.h"

namespace PhotonMapper {

  double PhotonMirrorDistance(const Photon &photon) {
    auto MirrorCentre = photon.m_RadiatorCell->GetMirrorCentre();
    const double R = photon.m_RadiatorCell->GetMirrorCurvature();
    const double a = 1.0;
    const double b = -2*photon.m_Direction.Dot(MirrorCentre - photon.m_Position);
    const double c = (MirrorCentre - photon.m_Position).Mag2() - R*R;
    const double s1 = (-b + TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    const double s2 = (-b - TMath::Sqrt(b*b - 4.0*a*c))/(2.0*a);
    return std::max(s1, s2);
  }

  void TracePhotonToMirror(Photon &photon) {
    const double s = PhotonMirrorDistance(photon);
    const Vector Mirror = photon.m_Position + s*photon.m_Direction;
    const Vector Normal = (Mirror - photon.m_RadiatorCell->GetMirrorCentre())/photon.m_RadiatorCell->GetMirrorCurvature();
    const Vector Vout = photon.m_Direction - 2*photon.m_Direction.Dot(Normal)*Normal;
    photon.m_Position = Mirror;
    photon.m_Direction = Vout;
    if(!photon.m_RadiatorCell->IsInsideThetaBoundary(photon)) {
      photon.m_Status = Photon::Status::MissedTheta;
    } else if(!photon.m_RadiatorCell->IsInsidePhiBoundary(photon)) {
      photon.m_Status = Photon::Status::MissedPhi;
    } else {
      photon.m_Status = Photon::Status::MirrorHit;
      photon.m_MirrorHitPosition = std::make_unique<Vector>(photon.m_Position);
    }
  }

  void TracePhotonToDetector(Photon &photon) {
    photon.m_Position += (photon.m_Position.Dot(Vector(0.0, 0.0, 1.0))/TMath::Abs(photon.m_Direction.Z()))*photon.m_Direction;
    photon.m_Status = photon.m_RadiatorCell->m_Detector.AddPhotonHit(photon) ? Photon::Status::DetectorHit : Photon::Status::DetectorMiss;
    if(!photon.m_RadiatorCell->IsInsideThetaBoundary(photon)) {
      photon.m_Status = Photon::Status::MissedDetectorPlane;
    }
  }

  void TracePhoton(Photon &photon) {
    TracePhotonToMirror(photon);
    if(photon.m_MirrorHitPosition) {
      TracePhotonToDetector(photon);
    }
  }

}
