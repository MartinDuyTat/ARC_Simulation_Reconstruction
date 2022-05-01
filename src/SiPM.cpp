// Martin Duy Tat 29th April 2022

#include<string>
#include"TGraph.h"
#include"TCanvas.h"
#include"SiPM.h"
#include"Photon.h"

SiPM::SiPM(): m_DetectorSizeX(0.05),
	      m_DetectorSizeY(0.05),
	      m_DetectorPositionX(0.0),
	      m_DetectorPositionY(0.0) {
}

void SiPM::AddPhotonHit(const Photon &photon) {
  m_PhotonHits.emplace_back(photon.m_Position.X(), photon.m_Position.Y());
}

void SiPM::PlotHits(const std::string &Filename) const {
  std::vector<double> xHits(m_PhotonHits.size()), yHits(m_PhotonHits.size());
  for(std::size_t i = 0; i < m_PhotonHits.size(); i++) {
    xHits[i] = m_PhotonHits[i].x*100.0;
    yHits[i] = m_PhotonHits[i].y*100.0;
  }
  TGraph Graph(xHits.size(), xHits.data(), yHits.data());
  TCanvas c("c", "", 800, 800);
  Graph.SetMarkerStyle(kFullDotLarge);
  Graph.SetMarkerStyle(kFullDotLarge);
  Graph.SetMarkerColor(kBlue);
  Graph.SetTitle("Photon hits;x (cm); y (cm)");
  Graph.Draw("AP");
  //TLegend Legend(0.7, 0.8, 0.9, 0.9);
  //Legend.AddEntry(&Graph1, "Aerogel 5GeV", "P");
  //Legend.AddEntry(&Graph2, "Gas 5GeV", "P");
  //Legend.Draw();
  c.Draw();
  c.SaveAs(Filename.c_str());
}
