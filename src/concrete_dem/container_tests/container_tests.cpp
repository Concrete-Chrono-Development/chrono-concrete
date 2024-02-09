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

using namespace chrono;
using namespace chrono::irrlicht;

chrono::collision::ChCollisionSystemType collision_type =
  chrono::collision::ChCollisionSystemType::BULLET;

std::vector<std::shared_ptr<ChBody>> create_container(ChSystemMulticoreSMC* sys,
						      double container_size){
  // container size given in [m], its postionion cannot be changed (bottom surface mid point
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
  side_wall_1->SetPos(ChVector<>(0, container_size + thickness/2, 0));
  side_wall_1->SetBodyFixed(true);
  side_wall_1->SetCollide(true);
  side_wall_1->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
			ChVector<>(thickness, container_size, container_size));
  side_wall_1->GetCollisionModel()->BuildModel();

  auto side_wall_2 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_2->SetMass(1);
  side_wall_2->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_2->SetPos(ChVector<>(0, -container_size - thickness/2, 0));
  side_wall_2->SetBodyFixed(true);
  side_wall_2->SetCollide(true);
  side_wall_2->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
			ChVector<>(thickness, container_size, container_size));
  side_wall_2->GetCollisionModel()->BuildModel();

  auto side_wall_3 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_3->SetMass(1);
  side_wall_3->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_3->SetPos(ChVector<>(-container_size - thickness/2, 0, 0));
  side_wall_3->SetBodyFixed(true);
  side_wall_3->SetCollide(true);
  side_wall_3->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
			ChVector<>(container_size, thickness, container_size));
  side_wall_3->GetCollisionModel()->BuildModel();

  auto side_wall_4 = std::shared_ptr<ChBody>(sys->NewBody());
  side_wall_4->SetMass(1);
  side_wall_4->SetInertiaXX(ChVector<>(1, 1, 1));
  side_wall_4->SetPos(ChVector<>(container_size + thickness/2, 0, 0));
  side_wall_4->SetBodyFixed(true);
  side_wall_4->SetCollide(true);
  side_wall_4->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(bottom_wall.get(), mat,
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


void create_concrete(ChSystemMulticoreSMC* sys, double container_size,
		     double mortar_layer_thick){
  // function creates sphere bodies to represent concrete in a container
  // for now it fills container box with given size [m], rest of the parameters
  // are defined as local variables, except for thickness of mortat layer, which is
  // DFC parameter, so I don't want to have two paralell definitions
  
  //  return true
}

int main(int argc, char* argv[]) {
  GetLog() << "Test application for implementation of DFC model in chrono::multicore\n";
  GetLog() << "Based on open source library projectchrono.org\nChrono version: "
	   << CHRONO_VERSION << "\n\n\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  ChSystemMulticoreSMC sys;
  sys.Set_G_acc(ChVector<>(0, 0, -9.81));
  std::vector<std::shared_ptr<ChBody>> container_walls = create_container(&sys, 0.15);
  std::shared_ptr<ChVisualSystem> vis;
  auto vis_irr = chrono_types::make_shared<ChVisualSystemIrrlicht>();
  vis_irr->AttachSystem(&sys);
  vis_irr->SetWindowSize(800, 600);
  vis_irr->SetWindowTitle("SMC callbacks");
  vis_irr->Initialize();
  vis_irr->AddLogo();
  vis_irr->AddSkyBox();
  vis_irr->AddCamera(ChVector<>(0, 0.05, -0.1));
  vis_irr->AddTypicalLights();
  vis = vis_irr;
  double simulation_time = 0;
  double time_step = 1e-03;
  while (simulation_time <= 1) {
    vis->BeginScene();
    vis->Render();
    vis->RenderGrid(ChFrame<>(VNULL, Q_from_AngX(CH_C_PI_2)), 12, 0.5);
    vis->RenderCOGFrames(1.0);
    sys.DoStepDynamics(time_step);
    vis->EndScene();
    simulation_time += time_step;
  }
  return 0;
}
