// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mariusz Warzecha
// =============================================================================
//
//  Application created to test DFC model implementation in chrono::multicore.
//  Implemented tests include three tests done in container: gravity-driven tests (1 and 2)
//  and hydrostatic confinement test.

#include <vector>
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono/collision/ChCollisionSystemBullet.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/particlefactory/ChParticleRemover.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/core/ChDistribution.h"
#include "ParticleEmitterMulticore.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::particlefactory;

chrono::collision::ChCollisionSystemType collision_type =
  chrono::collision::ChCollisionSystemType::BULLET;

std::vector<std::shared_ptr<ChBody>> create_container(ChSystemMulticoreSMC* sys,
						      double container_size){
  // container size given in [m], its positioning cannot be changed (bottom surface mid point
  // at the beginning of the global coordinate system
  // returned vector stores created objects 
  double thickness = 0.01;  // wall thickness (container will be build from boxes)
  double density = 1000;
  auto mat = chrono_types::make_shared<chrono::ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(2.05e11);
  mat->SetPoissonRatio(0.3);
  mat->SetRestitution(0.5);
  
  auto bottom_wall = std::shared_ptr<ChBody>(sys->NewBody());
  bottom_wall->SetMass(1);
  bottom_wall->SetInertiaXX(ChVector<>(1, 1, 1));
  bottom_wall->SetPos(ChVector<>(0, 0, -thickness/2));
  bottom_wall->SetBodyFixed(true);
  bottom_wall->SetCollide(true);
  bottom_wall->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
			ChVector<>(container_size, container_size, thickness));
  bottom_wall->GetCollisionModel()->BuildModel();

  auto side_wall_1 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_1->SetMass(1);
  side_wall_1->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_1->SetPos(ChVector<>(container_size/2 + thickness/2, 0, container_size/2));
  side_wall_1->SetBodyFixed(true);
  side_wall_1->SetCollide(true);
  side_wall_1->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(side_wall_1.get(), mat,
			ChVector<>(thickness, container_size, container_size));
  side_wall_1->GetCollisionModel()->BuildModel();

  auto side_wall_2 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_2->SetMass(1);
  side_wall_2->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_2->SetPos(ChVector<>(-container_size/2 - thickness/2, 0, container_size/2));
  side_wall_2->SetBodyFixed(true);
  side_wall_2->SetCollide(true);
  side_wall_2->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(side_wall_2.get(), mat,
			ChVector<>(thickness, container_size, container_size));
  side_wall_2->GetCollisionModel()->BuildModel();

  auto side_wall_3 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_3->SetMass(1);
  side_wall_3->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_3->SetPos(ChVector<>(0, -container_size/2 - thickness/2, container_size/2));
  side_wall_3->SetBodyFixed(true);
  side_wall_3->SetCollide(true);
  side_wall_3->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(side_wall_3.get(), mat,
			ChVector<>(container_size, thickness, container_size));
  side_wall_3->GetCollisionModel()->BuildModel();

  auto side_wall_4 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_4->SetMass(1);
  side_wall_4->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_4->SetPos(ChVector<>(0, container_size/2 + thickness/2, container_size/2));
  side_wall_4->SetBodyFixed(true);
  side_wall_4->SetCollide(true);
  side_wall_4->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(side_wall_4.get(), mat,
			ChVector<>(container_size, thickness, container_size));
  side_wall_4->GetCollisionModel()->BuildModel();
  
  sys->AddBody(bottom_wall);
  sys->AddBody(side_wall_1);
  sys->AddBody(side_wall_2);
  sys->AddBody(side_wall_3);
  sys->AddBody(side_wall_4);
  std::vector<std::shared_ptr<ChBody>> collection_of_walls;
  collection_of_walls.push_back(bottom_wall);
  collection_of_walls.push_back(side_wall_1);
  collection_of_walls.push_back(side_wall_2);
  collection_of_walls.push_back(side_wall_3);
  return collection_of_walls;
}

// define a class for concrete particle distribution
class DFCParticleDistr : public ChDistribution {
public:
  DFCParticleDistr(double min_size, double max_size, double mortar_layer, double mq)
    : d_0(max_size), d_a(min_size), h(mortar_layer), q(mq) {
  }
  virtual double GetRandom() override {
    double P = ChRandom();
    double aggregate_D = d_0 * pow((1 - P * (1 - pow(d_0, q)/pow(d_a, q))), (-1/q));
    return aggregate_D + 2 * h;
  }
private:
  double d_0;
  double d_a;
  double h;
  double q;
};

int main(int argc, char* argv[]) {
  GetLog() << "Test application for implementation of DFC model in chrono::multicore\n";
  GetLog() << "Based on open source library projectchrono.org\nChrono version: "
	   << CHRONO_VERSION << "\n\n\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  ChSystemMulticoreSMC sys;
  sys.Set_G_acc(ChVector<>(0, 0, -9.81));
  std::vector<std::shared_ptr<ChBody>> container_walls = create_container(&sys, 0.15);

  // concrete and DFC parameters (SI units)
  double containerVol = 0.15 * 0.15 * 0.15;
  double Vol = containerVol * 2;
  double h_layer = 4e-3;
  double minD = 0.005;  // minimum aggregate size
  double maxD = 0.01;  // maximum aggregate size
  double cement = 797;  // cement content kg/m^3
  double WtoC = 0.4;  // water-to-cement ratio
  double AtoC = 2.25;  // aggregate-to-cement ratio
  double rho_c = 3150;  // cement density kg/m^3
  double rho_w = 1000;  // water density kg/m^3
  double vair = 0;
  double nF = 0.5;  // Fuller exponent
  double va = 1 - cement / rho_c - (WtoC * cement) / rho_w - vair;  // aggregate volume fraction
  double va0 = (1 - pow((minD/maxD), nF)) * va;  // aggregate volume fraction with size bigger minD
  double Va0 = va0 * Vol;  // aggregate volume (with size bigger than minD)
  double rho_0 = cement * (1 + WtoC + AtoC);  // mixture density
  double targetmass = Va0 * rho_0;
  double targetVol = containerVol * 2.5;
  sys.GetSettings()->dfc_contact_param.E_Nm = 0.04e6;  // 2nd param --> 0.25 / 0.50 / 1 / 2 / 4 
  sys.GetSettings()->dfc_contact_param.E_Na = 100e6;
  sys.GetSettings()->dfc_contact_param.h = h_layer; // 1st param --> 0.75 / 1 / 1.10 / 1.20 / 1.40 
  sys.GetSettings()->dfc_contact_param.alfa_a = 0.25;
  sys.GetSettings()->dfc_contact_param.beta = 0.5;
  sys.GetSettings()->dfc_contact_param.sigma_t = 0.005e6;  // 3 rd param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.sigma_tau0 = 0.0005e6;  // 4th param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.eta_inf = 50;       /// 5th param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.kappa_0 = 100;
  sys.GetSettings()->dfc_contact_param.n = 1;
  sys.GetSettings()->dfc_contact_param.mi_a = 0.5;
  sys.GetSettings()->dfc_contact_param.E_Nm_s = 0.04e6;
  sys.GetSettings()->dfc_contact_param.E_Na_s = 100e6;
  sys.GetSettings()->dfc_contact_param.alfa_a_s = 0.25;
  sys.GetSettings()->dfc_contact_param.sigma_t_s = 0.005e6;
  sys.GetSettings()->dfc_contact_param.sigma_tau0_s = 0.0005e6;
  sys.GetSettings()->dfc_contact_param.eta_inf_s = 50;
  sys.GetSettings()->dfc_contact_param.mi_a_s = 0.5;
  sys.GetSettings()->dfc_contact_param.t = 4e-3;
  sys.GetSettings()->solver.contact_force_model = chrono::ChSystemSMC::ContactForceModel::DFC;

  // particle emitter
  ParticleEmitterMulticore emitter;
  emitter.SetSystem(&sys);
  emitter.SetParticlesPerSecond(10);
  emitter.SetParticleReservoir(9000);
  emitter.SetMortarLayer(h_layer);
  auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionRectangleOutlet>();
  emitter_positions->Outlet() = ChCoordsys<>(ChVector<>(0, 0, 0.2), QUNIT);
  emitter_positions->OutletWidth() = 0.135;
  emitter_positions->OutletHeight() = 0.135;
  emitter.SetParticlePositioner(emitter_positions);
  auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
  emitter.SetParticleAligner(emitter_rotations);
  auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityConstantDirection>();
  mvelo->SetDirection(-VECT_Z);
  mvelo->SetModulusDistribution(4);
  emitter.SetParticleVelocity(mvelo);	    

  std::shared_ptr<ChVisualSystem> vis;
  auto vis_irr = chrono_types::make_shared<ChVisualSystemIrrlicht>();
  vis_irr->AttachSystem(&sys);
  vis_irr->SetWindowSize(800, 600);
  vis_irr->SetWindowTitle("SMC callbacks");
  vis_irr->Initialize();
  vis_irr->AddLogo();
  vis_irr->AddSkyBox();
  vis_irr->AddCamera(ChVector<>(0.5, 0.5, 1));
  vis_irr->AddTypicalLights();
  vis = vis_irr;
  class MyCreatorForAll : public ChRandomShapeCreator::AddBodyCallback {
  public:
    virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
                           ChCoordsys<> mcoords,
                           ChRandomShapeCreator& mcreator) override {
      vis->BindItem(mbody);
      mbody->SetNoGyroTorque(true);
      sys->AddBody(mbody);
      //mbody->SetSleeping(true);
      std::cout << mbody->GetBodyFixed() << "\n";
    }
    ChVisualSystemIrrlicht* vis;
    ChSystemMulticoreSMC* sys;
    };
  emitter.SetVisualisation(vis_irr.get());
  double simulation_time = 0;
  double time_step = 1e-04;
  bool switch_val = true;
  while (vis->Run()) {
    vis->BeginScene();
    vis->Render();
    vis->RenderGrid(ChFrame<>(VNULL, Q_from_AngX(CH_C_PI_2)), 12, 0.5);
    //vis->RenderCOGFrames(1.0);
    emitter.EmitParticles(time_step);
    sys.Setup();
    auto body_list = sys.Get_bodylist();
    switch_val = true;
    sys.DoStepDynamics(time_step);
    for (auto body:body_list){
      std::cout << "Simulation time: " << simulation_time << "\n";
      std::cout << body->GetPos_dt() << "\n";
    }
    vis->EndScene();
    simulation_time += time_step;
  }
  return 0;
}
