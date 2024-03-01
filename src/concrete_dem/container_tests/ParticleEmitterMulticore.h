// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mariusz Warzecha
// =============================================================================
//
// Header file with implementation of a particle emitter working with multi-core systems.
// The particle emitters available in chrono particlefactory seems to not work with objects of
// type ChSystemMulticore (the particles are not added to system). Its intended to work only
// with spheres.

#ifndef PARTICLEEMITTERMULTICORE_H
#define PARTICLEEMITTERMULTICORE_H

#include "chrono/core/ChDistribution.h"
#include "chrono/physics/ChMaterialSurfaceSMC.h"
#include "chrono/assets/ChSphereShape.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

namespace chrono {
namespace particlefactory{
  class ParticleEmitterMulticore {
  public:
    ParticleEmitterMulticore(){
      diameter = chrono_types::make_shared<ChConstantDistribution>(0.02);
      density = chrono_types::make_shared<ChConstantDistribution>(1000);
      material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
      material->SetYoungModulus(2.05e11);
      material->SetPoissonRatio(0.3);
      material->SetRestitution(0.5);
      material->SetFriction(0.2);
      sys = nullptr;
      vis = nullptr;
      mortar_layer = 0;
      particle_reservoir = 0;
      particles_per_second = 0;
      particle_positioner = chrono_types::make_shared<ChRandomParticlePositionRectangleOutlet>();
      particle_aligner = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
      particle_velocity = chrono_types::make_shared<ChRandomParticleVelocity>();
    }
    void SetSystem(ChSystemMulticore* msys){sys = msys;}
    void SetVisualisation(chrono::irrlicht::ChVisualSystemIrrlicht* mvis) { vis = mvis;}
    void SetMortarLayer(double h){mortar_layer = h;}
    void SetParticlePositioner(std::shared_ptr<ChRandomParticlePosition> mc) {
      particle_positioner = mc;
    }
    void SetParticleAligner(std::shared_ptr<ChRandomParticleAlignment> mc) {
      particle_aligner = mc;
    }
    void SetParticleVelocity(std::shared_ptr<ChRandomParticleVelocity> mc) {
      particle_velocity = mc;
    }
    void SetParticleReservoir(int mc) { particle_reservoir = mc;}
    void SetParticlesPerSecond(double mc) { particles_per_second = mc; }
    void EmitParticles(double mdt, ChFrameMoving<> pre_transform = ChFrameMoving<>()) {
      double particles_per_step = mdt * particles_per_second;
      double done_particles_per_step = 0;
      while (true) {
	if (particle_reservoir <= 0){
	  return;
	}
	if (done_particles_per_step > particles_per_step) {
	  return;
	}
	// Create the particle
	ChCoordsys<> mcoords;
	mcoords.pos = particle_positioner->RandomPosition();
	mcoords.rot = particle_aligner->RandomAlignment();
	ChCoordsys<> mcoords_abs;
	mcoords_abs = mcoords >> pre_transform.GetCoord();
	std::shared_ptr<ChBody> mbody = RandomGenerate(mcoords_abs);
	ChVector<> mv_loc = particle_velocity->RandomVelocity();
	ChVector<> mv_abs = pre_transform.TransformDirectionLocalToParent(mv_loc);
	mbody->SetPos_dt(mv_abs);
	sys->AddBody(mbody);
	vis->BindItem(mbody);
	particle_reservoir -= 1;
	done_particles_per_step += 1;
      }
    }
  private:
    std::shared_ptr<ChDistribution> diameter;  // generate random diameter
    std::shared_ptr<ChDistribution> density;  // generate random density
    std::shared_ptr<ChMaterialSurfaceSMC> material;  // set material, although not used by DFC
    std::shared_ptr<ChRandomParticlePosition> particle_positioner;
    std::shared_ptr<ChRandomParticleAlignment> particle_aligner;
    std::shared_ptr<ChRandomParticleVelocity> particle_velocity;
    int particle_reservoir;
    ChSystemMulticore* sys;
    chrono::irrlicht::ChVisualSystemIrrlicht* vis;
    double mortar_layer;  // thickness of mortar layer for visualisation
    double particles_per_second;

    std::shared_ptr<ChBody> RandomGenerate(ChCoordsys<> mcoords) {
      auto bsphere = std::shared_ptr<ChBody>(sys->NewBody());
      double mrad = 0.5 * diameter->GetRandom();
      double mdensity = density->GetRandom();
      double mass = (4.0 / 3.0) * CH_C_PI * pow(mrad, 3) * mdensity;
      bsphere->SetInertiaXX((2.0 / 5.0) * mass * pow(mrad, 3) * ChVector<>(1, 1, 1));
      bsphere->SetMass(mass);
      bsphere->GetCollisionModel()->ClearModel();
      utils::AddSphereGeometry(bsphere.get(), material, mrad);
      bsphere->SetCollide(true);
      bsphere->GetCollisionModel()->BuildModel();
    //bsphere->SetBodyFixed(false);
      auto sphere1 = chrono_types::make_shared<ChSphereShape>(mrad);
      sphere1->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
      sphere1->SetOpacity(0.4f);
      auto sphere2 = chrono_types::make_shared<ChSphereShape>(mrad - mortar_layer);
      sphere2->SetTexture(GetChronoDataFile("textures/rock.jpg"));
      auto bsphere_vis = chrono_types::make_shared<ChVisualModel>();
      bsphere_vis->AddShape(sphere1);
      bsphere_vis->AddShape(sphere2);
      bsphere->AddVisualModel(bsphere_vis);
      bsphere->SetCoord(mcoords);
      return bsphere;
    }
    
  };

}
}

#endif
