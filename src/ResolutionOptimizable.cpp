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
					     const TracksPhotons &ParticlesPhotons):
  IOptimizable(),
  m_RadiatorCell(&radiatorCell),
  m_ParticlesPhotons(&ParticlesPhotons) {
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
  const double Resolution =
    ResolutionUtilities::fcn(x[0], x[1], x[2], *m_RadiatorCell, *m_ParticlesPhotons);
  return Resolution;
}

unsigned int ResolutionOptimizable::NumberOfParameters() const {
  return 3 - m_FixedParameters.size();
}

std::vector<Constraints> ResolutionOptimizable::GetConstraints() const {
  return m_ParameterSpace;
}

std::unordered_map<std::size_t, double>
ResolutionOptimizable::GetFixedParameters() const {
  return m_FixedParameters;
}
