// Martin Duy Tat 5th May 2022
/**
 * Settings is a class with only static member functions and the class contains all the settings and configurations
 * Each setting will belong to a set of settings with a given name
 * To obtain a setting, pass a string of the form "Name/Key"
 */

#ifndef SETTINGS
#define SETTINGS

#include<string>
#include<map>
#include<vector>

class Settings {
 public:
  /**
   * Delete constructor since class is "static"
   */
  Settings() = delete;
  /**
   * Add settings
   * @param Name Name of these settings
   * @param Filename Filename with settings
   */
  static void AddSettings(const std::string &Name, const std::string &Filename);
  /**
   * Get string setting
   */
  static std::string GetString(const std::string &Setting);
  /**
   * Get double setting
   */
  static double GetDouble(const std::string &Setting);
  /**
   * Get integer setting
   */
  static int GetInt(const std::string &Setting);
  /**
   * Get bool setting (if string is "true" this evaluates to true)
   */
  static bool GetBool(const std::string &Setting);
  /**
   * Get vector of int, comma separated in options file
   */
  static std::vector<int> GetIntVector(const std::string &Setting);
  /**
   * Check if setting exists
   */
  static bool Exists(const std::string &Setting);
 private:
  /**
   * Map where keys are name of settings file and values are maps with settings
   */
  static std::map<std::string, std::map<std::string, std::string>> m_Settings;
};

#endif
