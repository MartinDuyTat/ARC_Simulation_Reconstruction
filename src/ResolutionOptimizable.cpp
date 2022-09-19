// Martin Duy Tat 28th July 2022

#include<vector>
#include<string>
#include"ResolutionOptimizable.h"
#include"Settings.h"
#include"RadiatorCell.h"
#include"ParticleTrack.h"
#include"ResolutionUtilities.h"

using Constraints = de::IOptimizable::Constraints;

ResolutionOptimizable::ResolutionOptimizable(RadiatorCell &radiatorCell,
					     const Tracks &Particles):
  IOptimizable(),
  m_RadiatorCell(&radiatorCell),
  m_Particles(&Particles),
  m_Seed(Settings::GetSizeT("General/Seed")) {
  for(std::size_t i = 0; i < m_Parameters.size(); i++) {
    std::string ParameterName = "Optimisation/" + std::string(m_Parameters[i]) + "_";
    if(Settings::GetBool(ParameterName + "IsFixed")) {
      double Value = Settings::GetDouble(ParameterName + "value");
      m_FixedParameters.insert({i, Value});
    } else {
      double Min = Settings::GetDouble(ParameterName + "min");
      double Max = Settings::GetDouble(ParameterName + "max");
      m_ParameterSpace.emplace_back(Min, Max, true);
    }
  }
}

double ResolutionOptimizable::EvaluateCost(std::vector<double> x) const {
  std::array<double, 5> xx;
  std::size_t j = 0;
  for(std::size_t i = 0; i < m_Parameters.size(); i++) {
    const auto iter = m_FixedParameters.find(i);
    if(iter == m_FixedParameters.end()) {
      xx[i] = x[j];
      j++;
    } else {
      xx[i] = iter->second;
    }
  }
  const double Resolution =
    ResolutionUtilities::fcn(xx[0], xx[1], xx[2], xx[3], xx[4],
			     *m_RadiatorCell, *m_Particles, m_Seed, true);
  return Resolution;
}

std::size_t ResolutionOptimizable::NumberOfParameters() const {
  return m_Parameters.size() - m_FixedParameters.size();
}

std::vector<Constraints> ResolutionOptimizable::GetConstraints() const {
  return m_ParameterSpace;
}

std::unordered_map<std::size_t, double>
ResolutionOptimizable::GetFixedParameters() const {
  return m_FixedParameters;
}
