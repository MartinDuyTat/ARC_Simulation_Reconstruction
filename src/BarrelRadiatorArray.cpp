// Martin Duy Tat 23rd June 2022

#include<vector>
#include<stdexcept>
#include<string>
#include"Math/Vector3Dfwd.h"
#include"BarrelRadiatorArray.h"
#include"RadiatorCell.h"
#include"HalfRadiatorCell.h"
#include"Settings.h"
#include"Particle.h"
#include"BarrelRadiatorCell.h"

BarrelRadiatorArray::BarrelRadiatorArray():
  RadiatorArray::RadiatorArray(),
  m_FullArray(Settings::GetBool("General/FullArray")),
  m_BarrelRadius(Settings::GetDouble("ARCGeometry/Radius") +
		 Settings::GetDouble("RadiatorCell/CoolingThickness")),
  m_DeltaPhi(2.0*m_yHexDist/m_BarrelRadius) {
  if(m_FullArray) {
    m_Cells.reserve(2*m_NumberMainRowCells);
    for(std::size_t i = 0; i < m_NumberMainRowCells; i++) {
      // Main row
      m_Cells.emplace_back(std::make_unique<BarrelRadiatorCell>(i, 1, m_xHexDist));
    }
    for(std::size_t i = 0; i < m_NumberMainRowCells - 1; i++) {
      // Upper row
      m_Cells.emplace_back(std::make_unique<BarrelRadiatorCell>(i + 1, 2, m_xHexDist));
    }
    m_Cells.emplace_back(std::make_unique<HalfRadiatorCell>(m_NumberMainRowCells,
							    2, m_xHexDist));
  } else {
    m_Cells.emplace_back(std::make_unique<BarrelRadiatorCell>(0, 0, m_xHexDist));
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
BarrelRadiatorArray::DrawRadiatorArray() const {
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

const RadiatorCell* BarrelRadiatorArray::FindRadiator(Particle &particle) const {
  if(!particle.IsAtRadiator()) {
    throw std::runtime_error("Cannot find radiator since track is not at radiator");
  }
  if(m_FullArray) {
    auto Position = particle.GetPosition();
    if(TMath::Abs(Position.Phi()) > TMath::Pi()/2.0) {
      particle.MapPhi(TMath::Pi());
    }
    Position = particle.GetPosition();
    double x = Position.Z();
    const double Radius = TMath::Sqrt(Position.X()*Position.X() +
				      Position.Y()*Position.Y());
    // Check if particle hits outside half cell
    if(x > 8.5*m_xHexDist) {
      return nullptr;
    }
    // Account for curvature when calculating y coordinate (azimuthal)
    double ProjectedY = Position.Y()*m_BarrelRadius/Radius;
    double y = m_BarrelRadius*TMath::ASin(ProjectedY/m_BarrelRadius);
    // First check if track is on the correct side and not hitting endcap
    bool IsReflected = false;
    if(x < 0.0) {
      particle.ReflectZ();
      IsReflected = true;
      x = particle.GetPosition().Z();
    }
    // Make sure track is not below the main row
    while(IsBelowMainRow(x, y)) {
      particle.MapPhi(+m_DeltaPhi);
      Position = particle.GetPosition();
      ProjectedY = Position.Y()*m_BarrelRadius/Radius;
      y = m_BarrelRadius*TMath::ASin(ProjectedY/m_BarrelRadius);
    }
    // Make sure track is not above the upper row
    while(IsAboveUpperRow(x, y, m_BarrelRadius)) {
      particle.MapPhi(-m_DeltaPhi);
      Position = particle.GetPosition();
      ProjectedY = Position.Y()*m_BarrelRadius/Radius;
      y = m_BarrelRadius*TMath::ASin(ProjectedY/m_BarrelRadius);
    }
    if(IsAboveMainRow(x, y)) {
      // Particle hits the upper row
      const std::size_t xIndex = static_cast<std::size_t>(x/m_xHexDist) + 1;
      const std::size_t yIndex = 2;
      particle.MapPhi(-m_DeltaPhi/2.0);
      particle.SetPhiRotated(-m_DeltaPhi/2.0);
      return (*this)(xIndex, yIndex);
    } else {
      // Particle hits the main row
      const std::size_t xIndex = static_cast<std::size_t>((x + m_xHexDist/2.0)/m_xHexDist);
      const std::size_t yIndex = 1;
      if(xIndex == 0 && yIndex == 1 && IsReflected) {
	particle.ReflectZ();
      }
      return (*this)(xIndex, yIndex);
    }
  } else {
    return m_Cells.begin()->get();
  }
}

bool BarrelRadiatorArray::IsBelowMainRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const std::size_t xUnits = static_cast<std::size_t>(x/(m_xHexDist/2.0));
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

bool BarrelRadiatorArray::IsAboveUpperRow(double x, double y, double Radius) const {
  // The boundary on the top is just the boundary below, but shifted
  const double y_shift = y - m_DeltaPhi*Radius;
  return !IsBelowMainRow(x, y_shift);
}

bool BarrelRadiatorArray::IsAboveMainRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const std::size_t xUnits = static_cast<std::size_t>(x/(m_xHexDist/2.0));
  // Shift x coordinate
  const double x_shift = x - static_cast<double>(xUnits)*m_xHexDist/2.0;
  if(xUnits%2 == 0) {
    // If hexagon slope is downwards
    return x_shift > m_xHexDist - TMath::Sqrt(3)*y;
  } else {
    // If hexagon slope is downwards
    return x_shift < TMath::Sqrt(3)*y - m_xHexDist/2.0;
  }
}

int BarrelRadiatorArray::FindRadiatorIndex(std::size_t i, std::size_t j) const {
  if(m_FullArray) {
    if((j != 1 && j != 2) || (i == 0 && j == 2)) {
      throw std::invalid_argument("Radiator ("
				  + std::to_string(i)
				  + ", "
				  + std::to_string(j)
				  + ") does not exist");
    } else {
      if(j == 1) {
	if(i > 8) {
	  return -1;
	} else {
	  return static_cast<int>(i);
	}
      } else {
	if(i > 9) {
	  return -1;
	} else {
	  return static_cast<int>(m_NumberMainRowCells + i - 1);
	}
      }
    }
  } else {
    if(i != 0 && j != 0) {
      throw std::invalid_argument("Cannot have radiator cell index ("
				  + std::to_string(i)
				  + ", "
				  + std::to_string(j)
				  + " )");
    } else {
      return 0;
    }
  }
}
