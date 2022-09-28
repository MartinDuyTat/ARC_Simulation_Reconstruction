// Martin Duy Tat 29th April 2022

#include<string>
#include<iostream>
#include<memory>
#include<array>
#include<algorithm>
#include"TGraph.h"
#include"TCanvas.h"
#include"TLegend.h"
#include"TRandom.h"
#include"TPolyLine.h"
#include"TAxis.h"
#include"TPad.h"
#include"TMath.h"
#include"Math/Interpolator.h"
#include"SiPM.h"
#include"Photon.h"
#include"Settings.h"
#include"RadiatorCell.h"

SiPM::SiPM():
  m_DetectorSizeX(Settings::GetDouble("RadiatorCell/DetectorSize")),
  m_DetectorSizeY(m_DetectorSizeX),
  m_DetectorPositionX(0.0),
  m_DetectorPositionY(0.0),
  m_DetectorPositionZ(0.0),
  m_DetectorTilt(0.0),
  m_CoolingThickness(Settings::GetDouble("RadiatorCell/CoolingThickness")),
  m_MaxPDE(0.432),
  m_SwapXZ(Settings::GetString("General/BarrelOrEndcap") == "EndCap") {
}

PhotonHit SiPM::AddPhotonHit(Photon &photon) const {
  const Vector DetectorCentre(m_DetectorPositionX,
			      0.0,
			      m_DetectorPositionZ);
  const Vector PhotonXZPosition(photon.GetPosition().X(),
				0.0,
				photon.GetPosition().Z());
  PhotonHit photonHit{photon.GetPosition(),
                      &photon,
                      PhotonXZPosition - DetectorCentre};
  if(photon.GetStatus() != Photon::Status::MirrorHit) {
    return photonHit;
  }
  const double PhotonLambda = 1239.8/photon.GetEnergy();
  if(PhotonLambda < m_Lambda[0] || PhotonLambda > m_Lambda[15]) {
    std::cout << "Warning! Generated photon outside energy range of SiPM!\n";
  }
  if(!IsDetectorHit(photon)) {
    photon.UpdatePhotonStatus(Photon::Status::DetectorMiss);
    return photonHit;
  }
  const double RandomNumber = gRandom->Uniform(0.0, m_MaxPDE);
  const double Efficiency = m_Interpolator.Eval(PhotonLambda)*0.90*0.80;
  if(RandomNumber <= Efficiency) {
    photon.UpdatePhotonStatus(Photon::Status::DetectorHit);
  } else {
    photon.UpdatePhotonStatus(Photon::Status::EfficiencyMiss);
  }
  return photonHit;
}

void SiPM::PlotHits(const std::string &Filename,
		    const std::vector<PhotonHit> &photonHits) const {
  std::vector<double> xHits_Aerogel, yHits_Aerogel, xHits_Gas, yHits_Gas;
  for(const auto &PhotonHit : photonHits) {
    if(PhotonHit.m_Photon->GetRadiator() == Photon::Radiator::Aerogel) {
      xHits_Aerogel.push_back(PhotonHit.m_HitPosition.X()*100.0);
      yHits_Aerogel.push_back(PhotonHit.m_HitPosition.Y()*100.0);
    } else if(PhotonHit.m_Photon->GetRadiator() == Photon::Radiator::Gas) {
      xHits_Gas.push_back(PhotonHit.m_HitPosition.X()*100.0);
      yHits_Gas.push_back(PhotonHit.m_HitPosition.Y()*100.0);
    } else {
      std::cout << "Warning! Detector hit with unknown origin\n";
    }
  }
  TGraph Graph_Aerogel(static_cast<int>(xHits_Aerogel.size()),
		       xHits_Aerogel.data(),
		       yHits_Aerogel.data());
  TGraph Graph_Gas(static_cast<int>(xHits_Gas.size()),
		   xHits_Gas.data(),
		   yHits_Gas.data());
  TCanvas c("c", "", 1200, 1200);
  Graph_Aerogel.SetMarkerStyle(kFullDotLarge);
  Graph_Gas.SetMarkerStyle(kFullDotLarge);
  Graph_Aerogel.SetMarkerColor(kBlue);
  Graph_Gas.SetMarkerColor(kRed);
  Graph_Aerogel.SetTitle("Photon hits;x (cm); y (cm)");
  Graph_Aerogel.Draw("AP");
  gPad->SetLeftMargin(0.15f);
  gPad->SetBottomMargin(0.15f);
  Graph_Aerogel.GetXaxis()->SetLimits((m_DetectorPositionX
				       - m_DetectorSizeX/2.0)*100.0,
				      (m_DetectorPositionX
				       + m_DetectorSizeX/2.0)*100.0);
  Graph_Aerogel.SetMinimum((m_DetectorPositionY - m_DetectorSizeY/2.0)*100.0);
  Graph_Aerogel.SetMaximum((m_DetectorPositionY + m_DetectorSizeY/2.0)*100.0);
  Graph_Aerogel.Draw("AP");
  Graph_Gas.Draw("P SAME");
  TLegend Legend(0.7, 0.8, 0.9, 0.9);
  const int Momentum = Settings::GetInt("Particle/Momentum");
  const std::string AerogelLabel = "Aerogel " + std::to_string(Momentum) + " GeV";
  Legend.AddEntry(&Graph_Aerogel, AerogelLabel.c_str(), "P");
  const std::string GasLabel = "Gas " + std::to_string(Momentum) + " GeV";
  Legend.AddEntry(&Graph_Gas, GasLabel.c_str(), "P");
  Legend.Draw();
  c.Draw();
  c.SaveAs(Filename.c_str());
}

std::unique_ptr<TObject> SiPM::DrawSiPM(const Vector &RadiatorPosition) const {
  const double RadiatorZ = RadiatorPosition.Z();
  const double DetectorPositionX = m_DetectorPositionX + RadiatorPosition.X();
  std::array<double, 5> xPoints, zPoints;
  const double CosTheta = TMath::Cos(m_DetectorTilt);
  const double SinTheta = TMath::Sin(m_DetectorTilt);
  // Lower left corner
  xPoints[0] = DetectorPositionX - CosTheta*m_DetectorSizeX/2.0;
  zPoints[0] = RadiatorZ - 0.002 + m_DetectorPositionZ - SinTheta*m_DetectorSizeX/2.0;
  // Lower right corner
  xPoints[1] = DetectorPositionX + CosTheta*m_DetectorSizeX/2.0;
  zPoints[1] = RadiatorZ -0.002 + m_DetectorPositionZ + SinTheta*m_DetectorSizeX/2.0;
  // Upper right corner
  xPoints[2] = DetectorPositionX + CosTheta*m_DetectorSizeX/2.0;
  zPoints[2] = RadiatorZ - 0.0 + m_DetectorPositionZ + SinTheta*m_DetectorSizeX/2.0;
  // Upper left corner
  xPoints[3] = DetectorPositionX - CosTheta*m_DetectorSizeX/2.0;
  zPoints[3] = RadiatorZ - 0.0 + m_DetectorPositionZ - SinTheta*m_DetectorSizeX/2.0;
  // Back to left corner
  xPoints[4] = DetectorPositionX - CosTheta*m_DetectorSizeX/2.0;
  zPoints[4] = RadiatorZ - 0.002 + m_DetectorPositionZ - SinTheta*m_DetectorSizeX/2.0;
  if(m_SwapXZ) {
    std::swap(xPoints, zPoints);
  }
  TPolyLine DetectorLine(5, xPoints.data(), zPoints.data());
  DetectorLine.SetFillColor(kBlack);
  return std::make_unique<TPolyLine>(DetectorLine);
}

bool SiPM::IsDetectorHit(const Photon &photon) const {
  if(!photon.GetRadiatorCell()->IsInsideCell(photon)) {
    return false;
  }
  const double photon_x = photon.GetPosition().X();
  const double photon_z = photon.GetPosition().Z();
  const double XDist = TMath::Sqrt((photon_x - m_DetectorPositionX)*
				   (photon_x - m_DetectorPositionX)
				 + (photon_z - m_DetectorPositionZ)*
				   (photon_z - m_DetectorPositionZ));
  if(XDist > m_DetectorSizeX/2.0) {
    return false;
  } else if(TMath::Abs(photon.GetPosition().Y() - m_DetectorPositionY)
	    > m_DetectorSizeY/2.0) {
    return false;
  } else {
    return true;
  }
}

void SiPM::SetDetectorPosition(double x) {
  if(m_SwapXZ) {
    x = -x;
  }
  m_DetectorPositionX = x;
}

void SiPM::SetDetectorTilt(double Angle) {
  if(m_SwapXZ) {
    Angle = -Angle;
  }
  m_DetectorTilt = Angle;
  const double AbsAngle = TMath::Abs(Angle);
  const double DetectorSizeXOver2 = m_DetectorSizeX*0.5;
  if(AbsAngle > TMath::ASin(m_CoolingThickness/DetectorSizeXOver2)) {
    m_DetectorPositionZ = TMath::Sin(AbsAngle)*DetectorSizeXOver2;
    m_DetectorPositionZ -= m_CoolingThickness;
  } else {
    m_DetectorPositionZ = 0.0;
  }
}

const Interpolator SiPM::m_Interpolator{
  std::vector<double>(m_Lambda.begin(), m_Lambda.end()),
  std::vector<double>(m_PDE.begin(), m_PDE.end()),
  InterpolationType::kAKIMA
};

double SiPM::GetDetectorXPosition() const {
  return m_DetectorPositionX;
}

double SiPM::GetDetectorZPosition() const {
  return m_DetectorPositionZ;
}

double SiPM::GetDetectorTilt() const {
  return m_DetectorTilt;
}

double SiPM::GetDetectorSizeX() const {
  return m_DetectorSizeX;
}
