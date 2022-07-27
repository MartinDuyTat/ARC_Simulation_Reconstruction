// Martin Duy Tat 23rd June 2022

#include<vector>
#include<stdexcept>
#include<string>
#include"Math/Vector3Dfwd.h"
#include"RadiatorArray.h"
#include"RadiatorCell.h"
#include"Settings.h"
#include"ParticleTrack.h"

RadiatorArray::RadiatorArray():
  m_FullArray(Settings::GetBool("General/FullArray")),
  m_NumberMainRowCells(Settings::GetInt("ARCGeometry/CellsPerRow")),
  m_xHexDist(Settings::GetDouble("ARCGeometry/Length")/(2*m_NumberMainRowCells - 1)),
  m_yHexDist(m_xHexDist*TMath::Sqrt(3)/2.0) {
  if(m_FullArray) {
    m_Cells.reserve(2*m_NumberMainRowCells);
    for(int i = 0; i < m_NumberMainRowCells; i++) {
      // Main row
      m_Cells.emplace_back(i, 1, m_xHexDist);
    }
    for(int i = 0; i < m_NumberMainRowCells; i++) {
      // Upper row
      m_Cells.emplace_back(i + 1, 2, m_xHexDist);
    }
  } else {
    m_Cells.emplace_back(0, 0, m_xHexDist);
  }
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
RadiatorArray::DrawRadiatorArray() const {
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> RadiatorArrayObjects;
  for(const auto &Cell : m_Cells) {
    if(Cell.GetCellNumber().second == Settings::GetInt("EventDisplay/RowToDraw")) {
      auto RadiatorObjects = Cell.DrawRadiatorGeometry();
      RadiatorArrayObjects.insert(RadiatorArrayObjects.end(),
				  std::make_move_iterator(RadiatorObjects.begin()),
				  std::make_move_iterator(RadiatorObjects.end()));
    }
  }
  return RadiatorArrayObjects;
}

RadiatorIter RadiatorArray::operator()(int i, int j) {
  if(m_FullArray) {
    if((j != 1 && j != 2) || (i == 0 && j == 2)) {
      throw std::invalid_argument("Radiator ("
				  + std::to_string(i)
				  + ", "
				  + std::to_string(j)
				  + ") does not exist");
    } else {
      if(j == 1) {
	return m_Cells.begin() + i;
      } else {
	return m_Cells.begin() + m_NumberMainRowCells + i - 1;
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
      return m_Cells.begin();
    }
  }
}

RadiatorIter RadiatorArray::FindRadiator(ParticleTrack &particleTrack) {
  if(particleTrack.GetParticleLocation() != ParticleTrack::Location::EntranceWindow) {
    throw std::runtime_error("Cannot find radiator since track is not at entrance window");
  }
  if(m_FullArray) {
    auto Position = particleTrack.GetPosition();
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    double x = Position.Z();
    double y = TMath::ATan2(Position.Y(), Position.X())*Radius;
    // First check if track is on the correct side and not hitting endcap
    bool IsReflected = false;
    if(x < 0.0) {
      particleTrack.ReflectZ();
      IsReflected = true;
      x = particleTrack.GetPosition().Z();
    } else if(x > Settings::GetDouble("ARCGeometry/Length")/2.0) {
      throw std::runtime_error("Particle hitting endcap");
    }
    // Make sure track is not below the main row
    const double DeltaPhi = 2*m_yHexDist/Radius;
    while(IsBelowMainRow(x, y)) {
      particleTrack.MapPhi(+DeltaPhi);
      Position = particleTrack.GetPosition();
      y = TMath::ATan2(Position.Y(), Position.X())*Radius;
    }
    // Make sure track is not above the upper row
    while(IsAboveUpperRow(x, y)) {
      particleTrack.MapPhi(-DeltaPhi);
      Position = particleTrack.GetPosition();
      y = TMath::ATan2(Position.Y(), Position.X())*Radius;
    }
    if(IsAboveMainRow(x, y)) {
      // Particle hits the upper row
      const int xIndex = (x/m_xHexDist) + 1;
      const int yIndex = 2;
      if(xIndex == 0 && yIndex == 1 && IsReflected) {
	particleTrack.ReflectZ();
      }
      particleTrack.MapPhi(-DeltaPhi/2.0);
      return (*this)(xIndex, yIndex);
    } else {
      // Particle hits the main row
      const int xIndex = (x + m_xHexDist/2.0)/m_xHexDist;
      const int yIndex = 1;
      if(xIndex == 0 && yIndex == 1 && IsReflected) {
	particleTrack.ReflectZ();
      }
      return (*this)(xIndex, yIndex);
    }
  } else {
    return m_Cells.begin();
  }
}

bool RadiatorArray::IsBelowMainRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const int xUnits = x/(m_xHexDist/2.0);
  // Shift x coordinate
  const double x_shift = x - xUnits*m_xHexDist/2.0;
  if(xUnits%2 == 0) {
    // If hexagon slope is upwards
    return x_shift > TMath::Sqrt(3)*y + m_xHexDist;
  } else {
    // If hexagon slope is downwards
    return x_shift < -m_xHexDist/2.0 - TMath::Sqrt(3)*y;
  }
}

bool RadiatorArray::IsAboveUpperRow(double x, double y) const {
  // The boundary on the top is just the boundary below, but shifted
  const double y_shift = y - 2*m_yHexDist;
  return !IsBelowMainRow(x, y_shift);
}

bool RadiatorArray::IsAboveMainRow(double x, double y) const {
  // Count half-hexagon lengths along x
  const int xUnits = x/(m_xHexDist/2.0);
  // Shift x coordinate
  const double x_shift = x - xUnits*m_xHexDist/2.0;
  if(xUnits%2 == 0) {
    // If hexagon slope is downwards
    return x_shift > m_xHexDist - TMath::Sqrt(3)*y;
  } else {
    // If hexagon slope is downwards
    return x_shift < TMath::Sqrt(3)*y - m_xHexDist/2.0;
  }
}
