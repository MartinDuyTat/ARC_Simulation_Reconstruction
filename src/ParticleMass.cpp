// Martin Duy Tat 1st May 2022

#include<map>
#include"ParticleMass.h"

namespace ParticleMass {
  
  double GetMass(int PID) {
    static const std::map<int, double> Mass{
      {211, 0.139568},
      {321, 0.493677},
      {2212, 0.93827208816}
    };
    auto iter = Mass.find(PID);
    if(iter == Mass.end()) {
      return 0.0;
    } else {
      return iter->second;
    }
  }

}
