// Martin Duy Tat 27th May 2022
/**
 * EventDisplay is a class for storing and drawing all the objects that are in the event
 */

#ifndef EVENTDISPLAY
#define EVENTDISPLAY

#include<string>
#include<vector>
#include<utility>
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
  void AddObject(std::unique_ptr<TObject> Object, const std::string &Option = "");
  /**
   * Add objects to event display
   */
  void AddObject(std::vector<std::pair<std::unique_ptr<TObject>, std::string>> Objects);
 private:
  /**
   * Vector storing all the objects to be drawn
   */
  std::vector<std::pair<std::unique_ptr<TObject>, std::string>> m_EventObjects;
};

#endif
