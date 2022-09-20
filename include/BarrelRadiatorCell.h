// Martin Duy Tat 6th September 2022
/**
 * BarrelRadiatorCell inherits from RadiatorCell and contains the position of the radiator
 */

#ifndef BARRELRADIATORCELL
#define BARRELRADIATORCELL

#include<vector>
#include<memory>
#include<string>
#include"TObject.h"
#include"RadiatorCell.h"

class BarrelRadiatorCell: public RadiatorCell {
 public:
  /**
   * Constructor that sets up the geometry and a default position of the barrel radiator cell
   * @param CellColumnNumber The column number of this cell
   * @param CellRowNumber The row number of this cell
   * @param HexagonSize The length between two opposide edges (not points) of the hexagons
   */
  BarrelRadiatorCell(std::size_t CellColumnNumber, std::size_t CellRowNumber, double HexagonSize);
  /**
   * Delete copy constructor
   */
  BarrelRadiatorCell(const BarrelRadiatorCell &radiatorCell) = delete;
  /**
   * Get the radiator cell position in the global coordinates
   */
  virtual const Vector& GetRadiatorPosition() const override;
  /**
   * Function that checks if the position is inside the radiator
   */
  virtual bool IsInsideCell(const Vector &Position) const override;
  /**
   * Draw radiator geometry
   */
  virtual std::vector<std::pair<std::unique_ptr<TObject>, std::string>>
  DrawRadiatorGeometry() const override;
 protected:
  /**
   * Position of radiator cell in global coordinates, but rotated so that x is along the theta direction and y is along the phi direction
   */
  const Vector m_Position;
 private:
  /**
   * Helper function to get cell position based on cell number
   */
  Vector GetCellPosition(std::size_t CellColumnNumber, std::size_t CellRowNumber) const;
};

#endif
