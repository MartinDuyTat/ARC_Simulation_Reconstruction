// Martin Duy Tat 1st September 2022
/**
 * EndCapRadiatorCell inherits from RadiatorCell
 * It is not curved and its position is in the end caps
 */

#ifndef ENDCAPRADIATORCELL
#define ENDCAPRADIATORCELL

#include<vector>
#include<utility>
#include<memory>
#include<array>
#include"TObject.h"
#include"RadiatorCell.h"
#include"ARCVector.h"

class EndCapRadiatorCell: public RadiatorCell {
 public:
  /**
   * Default constructor that calls parent constructor
   */
  EndCapRadiatorCell(std::size_t CellColumnNumber,
		     std::size_t CellRowNumber,
		     double HexagonSize);
  /**
   * Delete copy constructor
   */
  EndCapRadiatorCell(const EndCapRadiatorCell &radiatorCell) = delete;
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const ARCVector &Position) const override;
  /**
   * Function that checks if the position is inside the radiator
   */
  bool IsInsideCell(const Vector &Position) const;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const override;
  /**
   * Function that checks if the detector is contained inside the radiator cell
   */
  virtual bool IsDetectorInsideCell() const override;
  /**
   * All the cell numbers that are valid for the end cap
   */
  static constexpr std::array<std::pair<std::size_t, std::size_t>, 23>
  m_ValidCells{{
    /*{1, 1},*/ {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1},
    {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, /*{8, 2},*/
    {3, 3}, {4, 3}, {5, 3}, {6, 3}, {7, 3},
    {5, 4}, {6, 4}, {7, 4},
    {6, 5}
  }};
  /**
   * All the cell numbers that are invalid, but next to valid cells for the end cap
   */
  static constexpr std::array<std::pair<std::size_t, std::size_t>, 9>
  m_NotValidCells{{
    {8, 1}, {9, 2}, {8, 3}, {8, 4}, {7, 5},
    {0, 1}, {1, 1}, {1, 2},
    {8, 2} //remove
  }};
 private:
  /**
   * Inner radius of end cap
   */
  const double m_InnerRadius;
  /**
   * Outer radius of end cap
   */
  const double m_OuterRadius;
  /**
   * Helper function to get cell position based on cell number
   */
  static Vector GetCellPosition(std::size_t CellColumnNumber,
				std::size_t CellRowNumber,
				double HexagonSize);
  /**
   * Helper function to determine the rotation from global to local coordinates
   */
  static Rotation3D GetCellOrientation(std::size_t CellColumnNumber,
				       std::size_t CellRowNumber,
				       double HexagonSize);
  /**
   * Helper function that determines the rotation angle of the coordinate system
   * The coordinate system is oriented so that the x-axis is radial
   */
  static double GetCoordinateRotation(std::size_t CellColumnNumber,
				      std::size_t CellRowNumber,
				      double HexagonSize);
  /**
   * Helper function that checks if the particle is within the radial acceptance
   * @param Position Global position vector
   */
  bool IsInsideRadialAcceptance(const Vector &Position) const;
};

#endif
