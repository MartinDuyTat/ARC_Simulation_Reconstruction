// Martin Duy Tat 23rd June 2022

#include<vector>
#include<stdexcept>
#include<string>
#include<utility>
#include<algorithm>
#include<array>
#include"Math/Vector3Dfwd.h"
#include"EndCapRadiatorArray.h"
#include"RadiatorCell.h"
#include"EndCapRadiatorCell.h"
#include"Settings.h"
#include"Particle.h"

EndCapRadiatorArray::EndCapRadiatorArray():
  RadiatorArray::RadiatorArray(),
  m_InnerRadius(Settings::GetDouble("ARCGeometry/EndCapInnerRadius")),
  m_OuterRadius(Settings::GetDouble("ARCGeometry/EndCapOuterRadius")) {
  constexpr std::size_t NumberCells = EndCapRadiatorCell::m_ValidCells.size();
  m_Cells.reserve(NumberCells);
  for(std::size_t i = 0; i < NumberCells; i++) {
    const std::size_t CellColumnNumber = EndCapRadiatorCell::m_ValidCells[i].first;
    const std::size_t CellRowNumber = EndCapRadiatorCell::m_ValidCells[i].second;
    m_Cells.emplace_back(std::make_unique<EndCapRadiatorCell>(CellColumnNumber,
							      CellRowNumber,
							      m_xHexDist));
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
EndCapRadiatorArray::DrawRadiatorArray() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorArrayObjects;
  for(const auto &Cell : m_Cells) {
    if(Cell->GetCellNumber().second == Settings::GetSizeT("EventDisplay/RowToDraw")) {
      auto RadiatorObjects = Cell->DrawRadiatorGeometry();
      RadiatorArrayObjects.insert(RadiatorArrayObjects.end(),
				  std::make_move_iterator(RadiatorObjects.begin()),
				  std::make_move_iterator(RadiatorObjects.end()));
    }
  }
  return RadiatorArrayObjects;
}

const RadiatorCell* EndCapRadiatorArray::FindRadiator(Particle &particle) const {
  if(!particle.IsAtRadiator()) {
    throw std::runtime_error("Cannot find radiator since track is not at entrance window");
  }
  const auto GlobalPosition = particle.GetPosition().GlobalVector();
  const double Radius = TMath::Sqrt(GlobalPosition.X()*GlobalPosition.X() +
				    GlobalPosition.Y()*GlobalPosition.Y());
  if(Radius < m_InnerRadius || Radius > m_OuterRadius) {
    return nullptr;
  }
  const double CosTheta_min = Settings::GetDouble("ARCGeometry/CosTheta_boundary");
  const double Theta = particle.GetPosition().GlobalVector().Theta();
  if(TMath::Cos(Theta) < CosTheta_min) {
    return nullptr;
  }
  // If track hits the other end cap, reflect
  if(particle.GetPosition().GlobalVector().Z() < 0.0) {
    particle.ReflectZ();
  }
  // If track hits the lower half, rotate by 180 degrees
  if(particle.GetPosition().GlobalVector().Phi() < 0.0) {
    particle.MapPhi(TMath::Pi());
  }
  // Rotate by 60 degrees until the track is in the interval [0, 60] degrees
  while(particle.GetPosition().GlobalVector().Phi() > TMath::Pi()/3.0) {
    particle.MapPhi(-TMath::Pi()/3.0);
  }
  bool IsReflected = false;
  // Reflect if particle is on the other side of the symmetry line
  if(particle.GetPosition().GlobalVector().Phi() > TMath::Pi()/6.0) {
    particle.MapPhi(-TMath::Pi()/6.0);
    particle.ReflectY();
    particle.MapPhi(TMath::Pi()/6.0);
    IsReflected = true;
  }
  // Loop from the bottom row until we find the correct row
  const double x = particle.GetPosition().GlobalVector().X();
  const double y = particle.GetPosition().GlobalVector().Y();
  std::size_t CellRowNumber = 1;
  while(true) {
    if(IsBelowThisRow(CellRowNumber + 1, x, y)) {
      break;
    } else {
      CellRowNumber++;
    }
  }
  std::size_t CellColumnNumber;
  if(CellRowNumber%2 == 1) {
    CellColumnNumber = static_cast<std::size_t>((x + m_xHexDist/2.0)/m_xHexDist);
  } else {
    CellColumnNumber = static_cast<std::size_t>(x/m_xHexDist) + 1;
  }
  if(IsReflected) {
    constexpr std::array<std::pair<std::size_t, std::size_t>, 4> ReflectCells{{
      {2, 2}, {3, 3}, {5, 4}, {6, 5}}};
    if(CellRowNumber == 1) {
      particle.ReflectY();
    } else if(std::find(ReflectCells.begin(),
			ReflectCells.end(),
			std::make_pair(CellColumnNumber, CellRowNumber))
	      != ReflectCells.end()) {
      particle.MapPhi(-TMath::Pi()/6.0);
      particle.ReflectY();
      particle.MapPhi(TMath::Pi()/6.0);
    }
  }
  return (*this)(CellColumnNumber, CellRowNumber);
}

bool EndCapRadiatorArray::IsBelowThisRow(std::size_t Row, double x, double y) const {
  // Shift y coordinate so that we're comparing with the main row
  const double y_shift = y - static_cast<double>(Row - 1)*m_yHexDist;
  if(Row%2 == 1) {
    return IsBelowOddRow(x, y_shift);
  } else {
    return IsBelowEvenRow(x, y_shift);
  }
}

bool EndCapRadiatorArray::IsBelowOddRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const int xUnits = static_cast<int>(x/(m_xHexDist/2.0));
  // Shift x coordinate
  const double x_shift = x - static_cast<double>(xUnits)*m_xHexDist/2.0;
  if(xUnits%2 == 0) {
    // If hexagon slope is upwards
    return x_shift > TMath::Sqrt(3)*y + m_xHexDist;
  } else {
    // If hexagon slope is downwards
    return x_shift < -m_xHexDist/2.0 - TMath::Sqrt(3)*y;
  }
}

bool EndCapRadiatorArray::IsBelowEvenRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const int xUnits = static_cast<int>(x/(m_xHexDist/2.0));
  // Shift x coordinate
  const double x_shift = x - static_cast<double>(xUnits)*m_xHexDist/2.0;
  if(xUnits%2 == 0) {
    // If hexagon slope is downwards
    return x_shift < -m_xHexDist/2.0 - TMath::Sqrt(3)*y;
  } else {
    // If hexagon slope is upwards
    return x_shift > TMath::Sqrt(3)*y + m_xHexDist;
  }
}

int EndCapRadiatorArray::FindRadiatorIndex(std::size_t i, std::size_t j) const {
  const auto iter = std::find(EndCapRadiatorCell::m_ValidCells.begin(),
			      EndCapRadiatorCell::m_ValidCells.end(),
			      std::make_pair(i, j));
  if(iter == EndCapRadiatorCell::m_ValidCells.end()) {
    const auto iter2 = std::find(EndCapRadiatorCell::m_NotValidCells.begin(),
				 EndCapRadiatorCell::m_NotValidCells.end(),
				 std::make_pair(i, j));
    if(iter2 != EndCapRadiatorCell::m_NotValidCells.end()) {
      return -1;
    } else {
      throw std::invalid_argument("Radiator ("
				  + std::to_string(i)
				  + ", "
				  + std::to_string(j)
				  + ") does not exist");
    }
  } else {
    const auto Index = iter - EndCapRadiatorCell::m_ValidCells.begin();
    return static_cast<int>(Index);
  }
}
