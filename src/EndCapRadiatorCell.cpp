// Martin Duy Tat 1st September 2022

#include<array>
#include<utility>
#include"TMath.h"
#include"EndCapRadiatorCell.h"
#include"RadiatorCell.h"
#include"Settings.h"

EndCapRadiatorCell::EndCapRadiatorCell(int CellColumnNumber,
				       int CellRowNumber,
				       double HexagonSize):
  RadiatorCell::RadiatorCell(CellColumnNumber, CellRowNumber, HexagonSize),
  m_HexagonDistY(m_HexagonSize*TMath::Sqrt(3)/2.0),
  m_Position(GetCellPosition(CellColumnNumber, CellRowNumber)) {
}

bool EndCapRadiatorCell::IsInsideCell(const Vector &Position) const {
  // Get x and y coordinates after mapping everything to first quadrant
  const double x = TMath::Abs(Position.X());
  const double y = TMath::Abs(Position.Y());
  // First part is checking the sloped part, the other is the vertical part
  return x < std::min(m_HexagonSize - y*TMath::Sqrt(3.0), m_HexagonSize*0.5);
}

std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
EndCapRadiatorCell::DrawRadiatorGeometry() const {
  return {};
}

const Vector& EndCapRadiatorCell::GetRadiatorPosition() const {
  return m_Position;
}

Vector EndCapRadiatorCell::GetCellPosition(int CellColumnNumber,
					   int CellRowNumber) const {
  if(std::find(m_ValidCells.begin(),
	       m_ValidCells.end(),
	       std::make_pair(CellColumnNumber, CellRowNumber)) ==
     m_ValidCells.end()) {
    throw std::invalid_argument("Invalid cell column number: "
				+ std::to_string(CellColumnNumber));
  }
  const double ZPosition = Settings::GetDouble("ARCGeometry/BarrelZ")
                         + m_CoolingThickness;
  const double YPosition = m_HexagonDistY*(CellRowNumber - 1);
  if(CellRowNumber%2 == 1) {
    const double XPosition = m_HexagonSize*CellColumnNumber;
    return Vector(XPosition, YPosition, ZPosition);
  } else {
    const double XPosition = m_HexagonSize*(CellColumnNumber - 0.5);
    return Vector(XPosition, YPosition, ZPosition);
  }
}
