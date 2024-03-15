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
#include "MyContactReport.h"

using namespace chrono;
using namespace chrono::irrlicht;
using namespace chrono::particlefactory;

std::shared_ptr<ChBody> AddSphere(ChSystemMulticore& sys, ChVector<> pos, ChVector<> vel, int id){
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);
  double radius = 0.008;   // in meters, I have doubt about used units
  double h = 4e-3;
  double density = 797 / (1 + 0.4 + 2.25);
  double mass = ((4.0 / 3.0) * 3.1415 * pow(radius, 3)) * 1000;
  auto ball = std::shared_ptr<chrono::ChBody>(sys.NewBody());
  ball->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 3) * chrono::ChVector<>(1, 1, 1));
  ball->SetMass(mass);
  ball->SetIdentifier(id);
  ball->SetPos(pos);
  ball->SetPos_dt(vel);
  //ball->SetWvel_par(ChVector<>(0, 0, 0));
  ball->GetCollisionModel()->ClearModel();
  utils::AddSphereGeometry(ball.get(), material, radius);
  ball->SetBodyFixed(false);
  ball->SetCollide(true);
  ball->GetCollisionModel()->BuildModel();
  auto sphere1 = chrono_types::make_shared<ChSphereShape>(radius);
  sphere1->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
  sphere1->SetOpacity(0.4f);
  auto sphere2 = chrono_types::make_shared<ChSphereShape>(radius - h);
  sphere2->SetTexture(GetChronoDataFile("textures/rock.jpg"));
  auto ball_vis = chrono_types::make_shared<ChVisualModel>();
  ball_vis->AddShape(sphere1);
  ball_vis->AddShape(sphere2);
  ball->AddVisualModel(ball_vis);
  sys.AddBody(ball);
  return ball;
}

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
  mat->SetFriction(0.2);
  
  auto bottom_wall = std::shared_ptr<ChBody>(sys->NewBody());
  bottom_wall->SetMass(1);
  bottom_wall->SetIdentifier(-1);
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
  side_wall_1->SetIdentifier(-2);
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
  side_wall_2->SetIdentifier(-3);
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
  side_wall_3->SetIdentifier(-4);
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
  side_wall_4->SetIdentifier(-5);
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
  //sys.SetCollisionSystemType(collision_type);
  sys.Set_G_acc(ChVector<>(0, 0, -9.81));

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
  real tolerance = 1e-3;
  uint max_iteration = 100;
  sys.GetSettings()->solver.max_iteration_bilateral = max_iteration;
  sys.GetSettings()->solver.tolerance = tolerance;
  sys.GetSettings()->collision.narrowphase_algorithm = collision::ChNarrowphase::Algorithm::HYBRID;
  sys.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
  std::vector<std::shared_ptr<ChBody>> container_walls = create_container(&sys, 0.15);
    
  // particle emitter
  ParticleEmitterMulticore emitter;
  emitter.SetSystem(&sys);
  emitter.SetParticlesPerSecond(1e6);
  emitter.SetParticleReservoir(9000);
  emitter.SetMortarLayer(h_layer);
  auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionRectangleOutlet>();
  emitter_positions->Outlet() = ChCoordsys<>(ChVector<>(0, 0, 0.01), QUNIT);
  emitter_positions->OutletWidth() = 0.135;
  emitter_positions->OutletHeight() = 0.135;
  emitter.SetParticlePositioner(emitter_positions);
  auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();
  emitter.SetParticleAligner(emitter_rotations);
  auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityConstantDirection>();
  mvelo->SetDirection(-VECT_Z);
  mvelo->SetModulusDistribution(0.004);
  emitter.SetParticleVelocity(mvelo);

  int ball_number = 0;
  std::vector<std::shared_ptr<ChBody>> created_balls;
  created_balls.push_back(AddSphere(sys, ChVector<>(0, 0,  0.008),
				    ChVector<>( 0, 0, 0), 1));
  created_balls.push_back(AddSphere(sys, ChVector<>(0, 0.012, 0.008), ChVector<>(0, 0, 0), 2));
  created_balls[0]->SetWvel_par(ChVector<>(0, 0, 100));
  //  created_balls.push_back(AddSphere(sys, ChVector<>(0.012, 0, 0.017), ChVector<>(0, 0, 0), 3));
  //  created_balls.push_back(AddSphere(sys, ChVector<>(0, 0.026, 0.017), ChVector<>(0, 0, 0), 4));
  //  created_balls.push_back(AddSphere(sys, ChVector<>(0.026, 0, 0.017), ChVector<>(0, 0, 0), 5));
  //  created_balls.push_back(AddSphere(sys, ChVector<>(0, 0.035, 0.017), ChVector<>(0, 0, 0), 6));
  ball_number += 2;
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
    }
    ChVisualSystemIrrlicht* vis;
    ChSystemMulticoreSMC* sys;
    };
  emitter.SetVisualisation(vis_irr.get());
  double simulation_time = 0;
  double time_step = 2e-05;
  bool switch_val = false;

  // text file for post-processing contact/collision data; first row stores column descriptors
  bool register_data = true;
  std::ofstream data_files[ball_number];
  if (register_data) {
    for (int i = 1; i <= ball_number; i++) {
      std::string file_name = "ball_" + std::to_string(i) + ".txt";
      const char* file_name_to_pass = file_name.c_str();  // convert type for object constructor
      data_files[i - 1].open(file_name_to_pass);
      data_files[i - 1] << "Sim time [s] ,"
			<< "Pos x [m], "
			<< "Pos y [m], "
			<< "Pos z [m], "
			<< "Vel x [m/s], "
			<< "Vel y [m/s], "
			<< "Vel z [m/s], "
			<< "Rot ang [-], "
			<< "Rot axis x, "
			<< "Rot axis y, "
			<< "Rot axis z, "
			<< "Ang_vel x [m/s], "
			<< "Ang_vel y [m/s], "
			<< "Ang_vel z [m/s], "
			<< "Force x [N], "
			<< "Force y [N], "
			<< "Force z [N], "
			<< "Torque x [Nm], "
			<< "Torque y [Nm], "
			<< "Torque z [Nm], "
			<< "Cont defor [m], "
			<< "Cont force x [N], "
			<< "Cont force res y [N], "
			<< "Cont force res z [N], "
			<< "Cont torque res x [Nm], "
			<< "Cont torque res y [Nm], "
			<< "Cont torque res z [Nm], " 
			<< "\n";
    }
  }
  while (vis->Run()) {
    vis->BeginScene();
    vis->Render();
    vis->RenderGrid(ChFrame<>(VNULL, Q_from_AngX(CH_C_PI_2)), 12, 0.5);
    //vis->RenderCOGFrames(1.0);
    if (switch_val) {
      emitter.EmitParticles(time_step);
      switch_val = false;
    }
    sys.DoStepDynamics(time_step);
    simulation_time += time_step;
    vis->EndScene();
    if (register_data) {
      std::shared_ptr<MyContactReport> contact_data_ptr = std::make_shared<MyContactReport>();
      sys.GetContactContainer()->ReportAllContacts(contact_data_ptr);
      for (int i = 0; i < ball_number; i++){
	data_files[i] << simulation_time << ", " ;
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetPos()[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetPos_dt()[j] << ", ";
	}
	data_files[i] << created_balls[i]->GetRotAngle() << ", ";
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetRotAxis()[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetWvel_par()[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetAppliedForce()[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  data_files[i] << created_balls[i]->GetAppliedTorque()[j] << ", ";
	}
	if (!contact_data_ptr->VectorOfCollisionData.empty()) {
	  bool added_contact_force = false;
	  for (auto e : contact_data_ptr->VectorOfCollisionData){
	    if (e.contactobjA == created_balls[i]->GetCollisionModel()->GetContactable() ||
		e.contactobjB == created_balls[i]->GetCollisionModel()->GetContactable()) {
	      data_files[i] << e.distance << ", ";
	      for (int j = 0; j < 3; j++) {
		data_files[i] << e.react_forces[j] << ", ";
	      }
	      for (int j = 0; j < 3; j++) {
		data_files[i] << e.react_torques[j] << ", ";
	      }
	      data_files[i] << "\n";
	      added_contact_force = true;
	    }
	  }
	  if (!added_contact_force) {
	    data_files[i] << "0, 0, 0, 0, 0, 0, 0\n";
	  }
	} else {
	  data_files[i] << "0, 0, 0, 0, 0, 0, 0\n";
	}
      }
    }
    
    
  }
  return 0;
}
