// Martin Duy Tat 6th September 2022

#include"TRandom.h"
#include"TMath.h"
#include"Utilities.h"
#include"Settings.h"
#include"RadiatorCell.h"
#include"BarrelRadiatorCell.h"

namespace Utilities {

  Vector VectorFromSpherical(double R, double CosTheta, double Phi) {
    const double SinTheta = TMath::Sqrt(1.0 - CosTheta*CosTheta);
    const double CosPhi = TMath::Cos(Phi);
    const double SinPhi = TMath::Sin(Phi);
    return Vector{R*CosPhi*SinTheta, R*SinPhi*SinTheta, R*CosTheta};
  }

  Vector GenerateRandomBarrelTrackZRange(double z_min, double z_max) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double z = gRandom->Uniform(z_min, z_max);
    const double CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
    const double Phi = Settings::GetBool("Particle/RandomPhi")
      ? gRandom->Uniform(-TMath::Pi(), TMath::Pi())
      : gRandom->Uniform(Settings::GetDouble("Particle/Phi_min"),
			 Settings::GetDouble("Particle/Phi_max"));
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
    return Momentum;
  }

  Vector GenerateRandomBarrelTrack(double &CosTheta, double &Phi) {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius");
    const double z = gRandom->Uniform(Settings::GetDouble("Particle/z_min"),
				      Settings::GetDouble("Particle/z_max"));
    CosTheta = z/TMath::Sqrt(z*z + Radius*Radius);
    Phi = Settings::GetBool("Particle/RandomPhi")
      ? gRandom->Uniform(-TMath::Pi(), TMath::Pi())
      : gRandom->Uniform(Settings::GetDouble("Particle/Phi_min"),
			 Settings::GetDouble("Particle/Phi_max"));
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = VectorFromSpherical(MomentumMag, CosTheta, Phi);
    return Momentum;
  }

  Vector GenerateRandomEndCapTrack() {
    const double Radius = Settings::GetDouble("ARCGeometry/Radius") + 0.10;
    const double z = Settings::GetDouble("ARCGeometry/BarrelZ");
    auto GetUniformCircle = [Radius] () {
      double x = Radius, y = Radius;
      while(x*x + y*y > Radius*Radius) {
	x = gRandom->Uniform(-Radius, Radius);
	y = gRandom->Uniform(-Radius, Radius);
      }
      return std::make_pair(x, y);
    };
    auto [x, y] = GetUniformCircle();
    const double MomentumMag = Settings::GetDouble("Particle/Momentum");
    const Vector Momentum = Vector(x, y, z).Unit()*MomentumMag;
    return Momentum;
  }

  Vector SwapXZForEndCap(const RadiatorCell *radiatorCell, Vector v) {
    const auto radiatorCellCast =
      dynamic_cast<const BarrelRadiatorCell*>(radiatorCell);
    if(!radiatorCellCast) {
      const double Temp = v.X();
      v.SetX(v.Z());
      v.SetZ(Temp);
    }
    return v;
  };

}
