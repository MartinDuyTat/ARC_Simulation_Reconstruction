// Martin Duy Tat 1st September 2022
/**
 * EndCapRadiatorCell inherits from RadiatorCell
 * It is not curved and its position is in the end caps
 */

#ifndef ENDCAPRADIATORCELL
#define ENDCAPRADIATORCELL

#include"RadiatorCell.h"

class EndCapRadiatorCell: public RadiatorCell {
 public:
  /**
   * Default constructor that calls parent constructor
   */
  EndCapRadiatorCell(int CellColumnNumber, int CellRowNumber, double HexagonSize);
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const Vector &Position) const override;
  /**
   * Get the radiator cell position in the global coordinates
   */
  virtual const Vector& GetRadiatorPosition() const override;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const override;
  /**
   * All the cell numbers that are valid for the end cap
   */
  static constexpr std::array<std::pair<int, int>, 21> m_ValidCells{{
    {2, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 1},
    {2, 2}, {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2},
    {3, 3}, {4, 3}, {5, 3}, {6, 3}, {7, 3},
    {5, 4}, {6, 4}, {7, 4},
    {6, 5}
  }};
  /**
   * All the cell numbers that are invalid, but next to valid cells for the end cap
   */
  static constexpr std::array<std::pair<int, int>, 21> m_NotValidCells{{
    {8, 1}, {8, 2}, {8, 3}, {8, 4}, {7, 5},
    {0, 1}, {1, 1}, {1, 2}
  }};
 private:
  /**
   * The distance between hexagons in the y-direction
   */
  const double m_HexagonDistY;
  /**
   * Position of radiator cell in global coordinates, but rotated so that x is along the theta direction and y is along the phi direction
   */
  const Vector m_Position;
  /**
   * Helper function to get cell position based on cell number
   */
  Vector GetCellPosition(int CellColumnNumber, int CellRowNumber) const;
};

#endif
