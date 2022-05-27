// Martin Duy Tat 27th May 2022
/**
 * EventDisplay is a class for storing and drawing all the objects that are in the event
 */

#ifndef EVENTDISPLAY
#define EVENTDISPLAY

#include<string>
#include<vector>
#include<utility>
#include<memory>
#include"TObject.h"

class EventDisplay {
 public:
  /**
   * Default constructor
   */
  EventDisplay() = default;
  /**
   * Draw the event display
   */
  void DrawEventDisplay(const std::string &Filename);
  /**
   * Add object to event display
   */
  void AddObject(std::unique_ptr<TObject> Object);
  /**
   * Add objects to event display
   */
  void AddObject(std::vector<std::unique_ptr<TObject>> Objects);
 private:
  /**
   * Vector storing all the objects to be drawn
   */
  std::vector<std::unique_ptr<TObject>> m_EventObjects;
};

#endif
