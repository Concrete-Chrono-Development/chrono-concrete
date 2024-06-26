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
//  Application created to test DFC model implementation for simple uses cases
// (two particles, particle and container wall etc.). This is second approach.
//

// comment the line below two switch off irrlicht visualisation
#define IRR

#include <vector>
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono/collision/ChCollisionSystemBullet.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "MyContactReport.h"
#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;

#ifdef IRR
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif
#include "WriteParticlesVTK.h"

std::shared_ptr<ChBody> AddSphere(ChSystemMulticore& sys, ChVector<> pos, ChVector<> vel,
				  int id, bool fixed){
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);
  double radius = 0.008;   // in meters, I have doubt about used units
  double h = sys.GetSettings()->dfc_contact_param.h;
  double density = 797 * (1 + 0.4 + 2.25);
  double mass = ((4.0 / 3.0) * 3.1415 * pow(radius, 3)) * density;
  auto ball = std::shared_ptr<chrono::ChBody>(sys.NewBody());
  ball->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 3) * chrono::ChVector<>(1, 1, 1));
  ball->SetMass(mass);
  ball->SetIdentifier(id);
  ball->SetPos(pos);
  ball->SetPos_dt(vel);
  ball->SetWvel_par(ChVector<>(0, 0, 0));
  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get(), material, radius);
  ball->SetBodyFixed(fixed);
  ball->SetCollide(true);
  ball->GetCollisionModel()->BuildModel();
  #ifdef IRR
  auto sphere1 = chrono_types::make_shared<ChSphereShape>(radius);
  sphere1->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
  sphere1->SetOpacity(0.4f);
  auto sphere2 = chrono_types::make_shared<ChSphereShape>(radius - h);
  sphere2->SetTexture(GetChronoDataFile("textures/rock.jpg"));
  auto ball_vis = chrono_types::make_shared<ChVisualModel>();
  ball_vis->AddShape(sphere1);
  ball_vis->AddShape(sphere2);
  ball->AddVisualModel(ball_vis);
  #endif
  sys.AddBody(ball);
  return ball;
}

std::shared_ptr<ChBody> create_wall(ChSystemMulticoreSMC& sys,
				    double wall_size){
  // container size given in [m], its positioning cannot be changed (bottom surface mid point
  // at the beginning of the global coordinate system
  // returned vector stores created objects 
  double thickness = 0.01;  // wall thickness (container will be build from boxes)
  double density = 1000;
  auto mat = chrono_types::make_shared<chrono::ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(2.05e11);
  mat->SetPoissonRatio(0.3);
  mat->SetRestitution(0.5);
  mat->SetFriction(0.2);
  auto bottom_wall = std::shared_ptr<ChBody>(sys.NewBody());
  bottom_wall->SetMass(1);
  bottom_wall->SetIdentifier(-1);
  bottom_wall->SetInertiaXX(ChVector<>(1, 1, 1));
  bottom_wall->SetPos(ChVector<>(0, 0, -thickness/2));
  bottom_wall->SetBodyFixed(true);
  bottom_wall->SetCollide(true);
  bottom_wall->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
			ChVector<>(wall_size, wall_size, thickness));
  bottom_wall->GetCollisionModel()->BuildModel();
  sys.AddBody(bottom_wall);
  return bottom_wall;
}

void test_two_particles(ChSystemMulticoreSMC& sys){
  std::shared_ptr<ChBody> ball_1, ball_2;
  ball_1 = AddSphere(sys, ChVector<>(0, 0, 0), ChVector<>(0, 0, 0), 1, true);
  ball_2 = AddSphere(sys, ChVector<>(0.015, 0, 0), ChVector<>(0, 0, 0), 2, false);
  std::vector<vec2> pair_to_debug;
  pair_to_debug.push_back(vec2(ball_1->GetId(), ball_2->GetId()));
  sys.GetSettings()->dfc_contact_param.debug_contact_pairs = pair_to_debug;
}

void test_particle_against_wall(ChSystemMulticoreSMC& sys){
  std::shared_ptr<ChBody> wall, ball;
  wall = create_wall(sys, 0.1);
  ball = AddSphere(sys, ChVector<>(0, 0, 0.006), ChVector<>(0, 0, 0), 1, false);
  std::vector<vec2> pair_to_debug;
  pair_to_debug.push_back(vec2(wall->GetId(), ball->GetId()));
  sys.GetSettings()->dfc_contact_param.debug_contact_pairs = pair_to_debug;
}

int main(int argc, char* argv[]) {
  GetLog() << "Test application for implementation of DFC model in chrono::multicore\n";
  GetLog() << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  std::string out_dir = "OUT_VTK_screens_for_presentation";
  if (!filesystem::create_directory(filesystem::path(out_dir))) {
    std::cerr << "Error creating directory" << out_dir << std::endl;
    return 1;
  }
  ChSystemMulticoreSMC sys;
  sys.Set_G_acc(ChVector<>(0, 0, 0));
  // concrete and DFC parameters (SI units)
  double h_layer = 4e-3;
  sys.GetSettings()->dfc_contact_param.E_Nm = 0.04e6;
  sys.GetSettings()->dfc_contact_param.E_Na = 100e6;
  sys.GetSettings()->dfc_contact_param.h = h_layer;
  sys.GetSettings()->dfc_contact_param.alfa_a = 0.25;
  sys.GetSettings()->dfc_contact_param.beta = 0.5;
  sys.GetSettings()->dfc_contact_param.sigma_t = 0.005e6;
  sys.GetSettings()->dfc_contact_param.sigma_tau0 = 0.0005e6;
  sys.GetSettings()->dfc_contact_param.eta_inf = 50;
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
  sys.GetSettings()->dfc_contact_param.t = h_layer;
  sys.GetSettings()->dfc_contact_param.debug_verbose = true;
  sys.GetSettings()->solver.contact_force_model = chrono::ChSystemSMC::ContactForceModel::DFC;
  sys.GetSettings()->collision.narrowphase_algorithm = collision::ChNarrowphase::Algorithm::HYBRID;
  sys.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

  double simulation_time = 0;
  double time_step = 1e-05;
  bool continue_simulation = true;
  //  test_two_particles(sys);
  test_particle_against_wall(sys);
#ifdef IRR
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
  #endif
#ifdef IRR
  while (vis->Run()) {
    vis->BeginScene();
    vis->Render();
    vis->RenderGrid(ChFrame<>(VNULL, Q_from_AngX(CH_C_PI_2)), 12, 0.5);
    vis->EndScene();
#else
    while (continue_simulation) {
#endif
    sys.DoStepDynamics(time_step);
    simulation_time += time_step;
    if (simulation_time > 0.5)
#ifdef IRR
      break;
#else
    continue_simulation = false;
#endif
  }
  return 0;
}


