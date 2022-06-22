// Martin Duy Tat 27th May 2022

#include<string>
#include<memory>
#include"TCanvas.h"
#include"TPad.h"
#include"TLine.h"
#include"TStyle.h"
#include"EventDisplay.h"
#include"Settings.h"

void EventDisplay::DrawEventDisplay(const std::string &Filename) {
  gStyle->SetLineScalePS(0.01);
  TCanvas c("c", "", Settings::GetInt("EventDisplay/CanvasWidth"), Settings::GetInt("EventDisplay/CanvasHeight"));
  gPad->DrawFrame(-Settings::GetDouble("ARCGeometry/Length")/2.0,
		  0.0,
		  Settings::GetDouble("ARCGeometry/Length")/2.0,
		  Settings::GetDouble("ARCGeometry/Radius") + Settings::GetDouble("RadiatorCell/RadiatorThickness"));
  TLine InnerARC(-Settings::GetDouble("ARCGeometry/Length")/2.0,
		 Settings::GetDouble("ARCGeometry/Radius"),
		 Settings::GetDouble("ARCGeometry/Length")/2.0,
		 Settings::GetDouble("ARCGeometry/Radius"));
  TLine OuterARC(-Settings::GetDouble("ARCGeometry/Length")/2.0,
		 Settings::GetDouble("ARCGeometry/Radius") + Settings::GetDouble("RadiatorCell/RadiatorThickness"),
		 Settings::GetDouble("ARCGeometry/Length")/2.0,
		 Settings::GetDouble("ARCGeometry/Radius") + Settings::GetDouble("RadiatorCell/RadiatorThickness"));
  InnerARC.SetLineColor(kBlack);
  OuterARC.SetLineColor(kBlack);
  m_EventObjects.push_back(std::make_pair(std::make_unique<TLine>(InnerARC), ""));
  m_EventObjects.push_back(std::make_pair(std::make_unique<TLine>(OuterARC), ""));
  for(auto &Object : m_EventObjects) {
    Object.first->Draw(Object.second.c_str());
  }
  c.Draw();
  c.SaveAs(Filename.c_str());
}

void EventDisplay::AddObject(std::unique_ptr<TObject> Object, const std::string &Option) {
  m_EventObjects.push_back(std::make_pair(std::move(Object), Option));
}

void EventDisplay::AddObject(std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects) {
  for(auto &Object : Objects) {
    AddObject(std::move(Object.first), Object.second);
  }
}
