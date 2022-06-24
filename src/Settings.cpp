// Martin Duy Tat 5th May 2022

#include<string>
#include<map>
#include<stdexcept>
#include<utility>
#include<fstream>
#include<sstream>
#include<vector>
#include<sstream>
#include<algorithm>
#include"Settings.h"

std::map<std::string, std::map<std::string, std::string>> Settings::m_Settings;

void Settings::AddSettings(const std::string &Name, const std::string &Filename) {
  if(m_Settings.find(Name) != m_Settings.end()) {
    throw std::runtime_error("Settings " + Name + " already exists");
  }
  std::map<std::string, std::string> NewSettings;
  std::ifstream File(Filename);
  std::string Line;
  while(std::getline(File, Line)) {
    Line = Line.substr(0, Line.find('#'));
    if(Line.empty()) {
      continue;
    }
    std::stringstream ss(Line);
    std::string Key, Value;
    ss >> Key >> Value;
    if(NewSettings.find(Key) != NewSettings.end()) {
      throw std::runtime_error("Key " + Key + " in settings " + Name + " already exists");
    }
    NewSettings.insert({Key, Value});
  }
  File.close();
  m_Settings.insert({Name, std::move(NewSettings)});
}

std::string Settings::GetString(const std::string &Setting) {
  const std::size_t SlashPos = Setting.find('/');
  std::string Name = Setting.substr(0, SlashPos);
  auto iter1 = m_Settings.find(Name);
  if(iter1 == m_Settings.end()) {
    throw std::runtime_error("Cannot find settings name " + Name);
  }
  std::string Key = Setting.substr(SlashPos + 1);
  auto iter2 = iter1->second.find(Key);
  if(iter2 == iter1->second.end()) {
    throw std::runtime_error("Cannot find key " + Key + " in settings " + Name);
  }
  return iter2->second;
}

double Settings::GetDouble(const std::string &Setting) {
  return std::stod(GetString(Setting));
}

int Settings::GetInt(const std::string &Setting) {
  return std::stoi(GetString(Setting));
}

bool Settings::GetBool(const std::string &Setting) {
  return GetString(Setting) == "true";
}

std::vector<int> Settings::GetIntVector(const std::string &Setting) {
  std::string CommaSeparatedList = GetString(Setting);
  std::replace(CommaSeparatedList.begin(), CommaSeparatedList.end(), ',', ' ');
  std::stringstream ss(CommaSeparatedList);
  std::vector<int> List;
  int Number;
  while(ss >> Number) {
    List.push_back(Number);
  }
  return List;
}

bool Settings::Exists(const std::string &Setting) {
  try {
    GetString(Setting);
  } catch(const std::runtime_error&) {
    return false;
  }
  return true;
}
