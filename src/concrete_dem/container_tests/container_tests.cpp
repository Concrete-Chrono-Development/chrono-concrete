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

// comment the line below two switch off irrlicht visualisation
//#define IRR

#include <vector>
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono/collision/ChCollisionSystemBullet.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/particlefactory/ChParticleRemover.h"
#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/core/ChDistribution.h"
#include "MyContactReport.h"
#include "chrono/core/ChMathematics.h"
#include "chrono_thirdparty/rapidjson/prettywriter.h"
#include "chrono_thirdparty/rapidjson/stringbuffer.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLinkMotorLinearForce.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

using namespace chrono;
using namespace chrono::particlefactory;

#ifdef IRR
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif
#include "WriteParticlesVTK.h"

double calc_aggVolFrac(std::vector<std::shared_ptr<ChBody>> &bodylist,
		       double hlayer, double specimenVol){
  double avfrac = 0;
  for (auto body : bodylist) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    double radius = body->GetCollisionModel()->GetShapeDimensions(0)[0] - hlayer;
    avfrac += 1.33333333333333333333333 * CH_C_PI * radius * radius * radius;
  }
  return avfrac / specimenVol;
}

void delete_particle(ChSystemMulticoreSMC& sys, double limit) {
  std::list<std::shared_ptr<ChBody>> to_delete;
  for (auto body : sys.Get_bodylist()) {
    if (body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<> pos = body->GetPos();
    if (pos.z() > limit || pos.z() < 0.0 || abs(pos.x()) > limit/2 || abs(pos.y()) > limit/2){
      to_delete.push_back(body);
    }
  }
  std::list<std::shared_ptr<ChBody>>::iterator ibody = to_delete.begin();
  while (ibody != to_delete.end()){
    sys.Remove((*ibody));
    ++ibody;
  }
}

std::vector<std::shared_ptr<ChBody>> list_particle_for_removal(ChSystemMulticoreSMC& sys, 
							      double limit) {
  std::vector<std::shared_ptr<ChBody>> to_delete;
  for (auto body : sys.Get_bodylist()) {
    if (body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<> pos = body->GetPos();
    if (pos.z() > limit || pos.z() < 0.0 || abs(pos.x()) > limit/2 || abs(pos.y()) > limit/2){
      to_delete.push_back(body);
    }
  }
  return to_delete;
}

// define a class for concrete particle distribution
class DFCParticleDistr : public ChDistribution {
public:
  DFCParticleDistr(double min_size, double max_size, double mortar_layer, double mq)
    : d_0(min_size), d_a(max_size), h(mortar_layer), q(mq) {
  }
  virtual double GetRandom() override {
    double P = ChRandom();
    double aggregate_D = d_0 * pow((1 - P * (1 - pow(d_0, q)/pow(d_a, q))), (-1/q));
    return aggregate_D + 2 * h;
  }
  double GetHLayer() {return h;}
  double GetMaxSize() {return d_a;}
private:
  double d_0;
  double d_a;
  double h;
  double q;
};

void AddSphereLayers(int layer_number, int box_number, double start_height,
		     ChSystemMulticoreSMC& sys, DFCParticleDistr d_param) {
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);
  double radius = 0;  0.5 * d_param.GetRandom();
  double box_size = d_param.GetMaxSize() + 0.001;  // max aggregate size + 1 mm for safety
  double density = 0.5* 797 * (1 + 0.4 + 2.25);
  double h = d_param.GetHLayer();
  double mass = 0;
  double shift_x[4] = {box_size/2, 0, -box_size/2, 0};
  double shift_y[4] = {0, box_size/2, 0, -box_size/2};
  double total_mass = 0;
  int l = -1;
  for (int k = 0; k < layer_number; ++k) {
    if (l < 4){
      ++l;
    }
    else {
      l = 0;
    }
    for (int j = 0; j < box_number; ++j) {
      for (int i = 0; i < box_number; ++i) {
	auto ball = std::shared_ptr<chrono::ChBody>(sys.NewBody());
	radius = 0.5 * d_param.GetRandom();
	double mass = ((4.0 / 3.0) * 3.1415 * pow(radius, 3)) * density;
	total_mass += mass;
	ball->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 3) * chrono::ChVector<>(1, 1, 1));
	ball->SetMass(mass);
	ball->SetPos(ChVector<>(i * box_size +  0.5 * box_size - 0.5 * box_size * box_number
				+ shift_x[l],
				j * box_size + 0.5 * box_size - 0.5 * box_size * box_number
				+ shift_y[l],
				start_height + k * box_size));
	ball->SetPos_dt(ChVector<>(0, 0, 0));
	ball->SetWvel_par(ChVector<>(0, 0, 0));
	ball->GetCollisionModel()->ClearModel();
	utils::AddSphereGeometry(ball.get(), material, radius);
	ball->SetBodyFixed(false);
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
      }
    }
  }
  GetLog() << "Total particle mass generated: " << total_mass << "\n";
  return;
}

std::shared_ptr<ChBody> AddSphere(ChSystemMulticore& sys, ChVector<> pos, ChVector<> vel, int id){
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);
  double radius = 0.008;   // in meters, I have doubt about used units
  double h = 4e-3;
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
  ball->SetBodyFixed(false);
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
  collection_of_walls.push_back(side_wall_4);
  return collection_of_walls;
}

// function to write reaction forces acting on container walls
void write_wall_forces(std::vector<std::shared_ptr<ChBody>> walls, std::string file_name,
		       double time){
  std::vector<std::string> wall_names {"Bottom wall: ", "Wall 1: ", "Wall 2: ",
				       "Wall 3: ", "Wall 4: "};
  std::ofstream file_to_write;
  file_to_write.open(file_name, std::ios_base::app);
  file_to_write << "Current time step: " << time << "\n";
  for (int i = 0; i < 5; ++i) {
    ChVector<> temp_force = walls[i]->GetContactForce();
    file_to_write << wall_names[i] << "x: " << temp_force.x() << " y: " << temp_force.y()
		  << " z: " << temp_force.z() << "\n";
  }
  file_to_write << "\n";
  file_to_write.close();
}

// function to export all contact forces acting on given container wall
void write_wall_contacts(ChSystemMulticore& sys, std::shared_ptr<ChBody> wall,
			 std::string file_name, double time){
  std::ofstream file_to_write;
  file_to_write.open(file_name, std::ios_base::app);
  file_to_write << "Current time step: " << time << "\n";
  std::shared_ptr<MyContactReport> contact_data_ptr = std::make_shared<MyContactReport>();
  sys.GetContactContainer()->ReportAllContacts(contact_data_ptr);
  if (!contact_data_ptr->VectorOfCollisionData.empty()) {
    for (auto e : contact_data_ptr->VectorOfCollisionData){
      if (e.contactobjA == wall->GetCollisionModel()->GetContactable()) {
	for (int j = 0; j < 3; j++) {
	  file_to_write << e.pA[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  file_to_write << e.react_forces[j] << ",";
	}
	file_to_write << e.distance << "\n";
      }
      if (e.contactobjB == wall->GetCollisionModel()->GetContactable()) {
	for (int j = 0; j < 3; j++) {
	  file_to_write << e.pB[j] << ", ";
	}
	for (int j = 0; j < 3; j++) {
	  file_to_write << e.react_forces[j] << ",";
	}
	file_to_write << e.distance << "\n";
      }
    }
  }
  else {
    file_to_write << "No contacts to report. Empty vector with contacts.\n";
    file_to_write.close();
    return;
  }
  file_to_write << "\n";
  file_to_write.close();
  return;
}

// function to write reaction forces acting on container walls
void write_cover_data(std::shared_ptr<ChBody> cover, std::string file_name,
		       double time){
  std::ofstream file_to_write;
  file_to_write.open(file_name, std::ios_base::app);
  file_to_write << "Current time step for cover data: " << time << "\n";
  ChVector<> temp_vector = cover->GetAppliedForce();
  file_to_write << "Force applied to cover " << "x= " << temp_vector.x() << " y= " << temp_vector.y()
		<< " z= " << temp_vector.z() << "\n";
  temp_vector = cover->Get_accumulated_force();
  file_to_write << "Accumulated force " << "x= " << temp_vector.x() << " y= " << temp_vector.y()
		<< " z= " << temp_vector.z() << "\n";
  temp_vector = cover->GetPos();
  file_to_write << "Position of cover " << "x= " << temp_vector.x() << " y= " << temp_vector.y()
		<< " z= " << temp_vector.z() << "\n";
  temp_vector = cover->GetPos_dt();
  file_to_write << "Velocity of cover " << "x= " << temp_vector.x() << " y= " << temp_vector.y()
		<< " z= " << temp_vector.z() << "\n";
  auto force_list = cover->GetForceList();
  file_to_write << "Size of force list is: " << force_list.size() << "\n";
  for (auto force : force_list){
    file_to_write << force << "\n";
  }
  file_to_write << "\n";
  file_to_write.close();
}


// function checking, if particles are inside the box
bool any_particle_inside(ChSystemMulticore& sys){
  auto body_list = sys.Get_bodylist();
  int particles_inside = 0;
  for (auto body : body_list) {
    if ((body->GetPos()).Length() < 0.15)
      ++particles_inside;
    if (particles_inside > 10)
      return true;
  }
   return false;
}

// find particle with highest coordinate of contact surface along z-direction, returns index
unsigned int find_particle_highest_z(ChSystemMulticoreSMC& sys) {
  double current_z_max = 0;
  unsigned int index = 0;
  auto body_list = sys.Get_bodylist();
  for (auto body : body_list) {
    if (body->GetIdentifier() > 0) {
      if (body->GetPos().z() + body->GetCollisionModel()->GetShapeDimensions(0)[0] > current_z_max ) {
	current_z_max = body->GetPos().z() + body->GetCollisionModel()->GetShapeDimensions(0)[0];
	index = body->GetId();
      }
    }
  }
  return index;
}

// function printing current energy status
ChVector<> print_energy_status(ChSystemMulticore& sys){
  float total_trans_kin_e = 0;
  float total_rot_kin_e = 0;
  auto body_list = sys.Get_bodylist();
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    auto mass = body->GetMass();
    auto inertia = body->GetInertiaXX();
    float trans_kin_e = 0;  // transnational kinetic energy
    float rot_kin_e = 0;  // rotational kinetic energy
    ChVector<float> vel_t = body->GetPos_dt();
    ChVector<float> vel_r = body->GetWvel_par();
    trans_kin_e = 0.5 * mass * (pow(vel_t.x(), 2) + pow(vel_t.y(), 2) + pow(vel_t.z(), 2));
    rot_kin_e = 0.5 * (inertia.x()*pow(vel_r.x(), 2) + inertia.y()*pow(vel_r.y(), 2)
		       + inertia.z()*pow(vel_r.z(), 2));
    total_trans_kin_e += trans_kin_e;
    total_rot_kin_e += rot_kin_e;
  }
  GetLog() << "Total transnational kinetic energy: " << total_trans_kin_e;
  GetLog() << " Total rotational kinetic energy: " << total_rot_kin_e << "\n";
  return ChVector<>(total_trans_kin_e, total_rot_kin_e, 0);
}

// function adding cover and 0.1 kPa pressure for formwork pressure
std::shared_ptr<ChBody> add_cover_formwork_pressure(ChSystemMulticore& sys, double container_size) {
  // create ChBody for cover
  double thickness = 0.01;  // wall thickness (same as for container walls)
  double density = 10; // it may impact the concrete load due to gravity, keep it low for now 
  auto mat = chrono_types::make_shared<chrono::ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(2.05e11);  // material parameters are not relevant, overloaded by DFC model
  mat->SetPoissonRatio(0.3);
  mat->SetRestitution(0.5);
  mat->SetFriction(0.2);
 
  auto cover = std::shared_ptr<ChBody>(sys.NewBody());
  cover->SetMass(2.25/9.81);  // use gravity load as pressure 
  cover->SetIdentifier(-6);
  cover->SetInertiaXX(ChVector<>(0.01, 0.01, 0.01));
  cover->SetPos(ChVector<>(0, 0, container_size + 1.32*thickness));
  cover->SetBodyFixed(false);
  cover->SetCollide(true);
  cover->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(cover.get(), mat,
			ChVector<>(container_size, container_size, thickness));
  cover->GetCollisionModel()->BuildModel();
  sys.AddBody(cover);

  // Create basement (copied from Bahar code, not sure if it is necessary, may use bottom plate)    
  auto mtruss = chrono_types::make_shared<ChBody>();
  mtruss->SetBodyFixed(true);
  mtruss->SetIdentifier(-7);
  sys.AddBody(mtruss);
  
  // constrain the cover (only movement in z axis is allowed)
  auto constr_cover = chrono_types::make_shared<ChLinkMateGeneric>(true, true, false,
  								   true, true, true);
  constr_cover->Initialize(cover, mtruss, false, cover->GetFrame_REF_to_abs(),
  			   mtruss->GetFrame_REF_to_abs());
  sys.Add(constr_cover);
  double pressure = 0.1 * 1000;  // 0.1 kPa expressed in Pa
  double pres_load = 0.15*0.15*pressure;
  // add force to compensate weight of cover
  //  auto load_container = chrono_types::make_shared<ChLoadContainer>();
  //  auto weight_compensation = chrono_types::make_shared<ChLoadBodyForce>(cover, ChVector<>(0, 0, 9.81),
  //									false, ChVector<>(0, 0, 0));
  //  auto pressure_load = chrono_types::make_shared<ChLoadBodyForce>(cover, ChVector<>(0, 0, -pres_load),
  //									false, ChVector<>(0, 0, 0));
  //  load_container->Add(weight_compensation);
  //  load_container->Add(pressure_load);
  //  sys.Add(load_container);

  // create pressure load
  //  auto motorZ = chrono_types::make_shared<ChLinkMotorLinearForce>();
  //  motorZ->Initialize(cover, mtruss, cover->GetFrame_REF_to_abs());
  //  auto loadZ = chrono_types::make_shared<ChFunction_Const>(-pressure * 0.15 * 0.15);
  //  motorZ->SetForceFunction(loadZ);
  //  sys.Add(motorZ);
  return cover;
}

int main(int argc, char* argv[]) {
  GetLog() << "Test application for implementation of DFC model in chrono::multicore\n";
  GetLog() << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  std::string out_dir = "OUT_VTK_AGV_Set_0_but_sigma_t_0_half_ro";
  if (!filesystem::create_directory(filesystem::path(out_dir))) {
    std::cerr << "Error creating directory" << out_dir << std::endl;
    return 1;
  }
  std::string terminal_log_file = out_dir + "/AGV_Set_0_terminal_log.txt";
  std::ofstream terminal_file(terminal_log_file);
  terminal_file << "Test application for implementation of DFC model in chrono::multicore\n";
  terminal_file << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";

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
  sys.GetSettings()->dfc_contact_param.sigma_t = 0;//0.005e6;  3 rd param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.sigma_tau0 = 0.0005e6;  // 4th param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.eta_inf = 50;       /// 5th param 0.25 / 0.50 / 1 / 2 / 4
  sys.GetSettings()->dfc_contact_param.kappa_0 = 100;
  sys.GetSettings()->dfc_contact_param.n = 1;
  sys.GetSettings()->dfc_contact_param.mi_a = 0.5;
  sys.GetSettings()->dfc_contact_param.E_Nm_s = 0.04e6;
  sys.GetSettings()->dfc_contact_param.E_Na_s = 100e6;
  sys.GetSettings()->dfc_contact_param.alfa_a_s = 0.25;
  sys.GetSettings()->dfc_contact_param.sigma_t_s = 0; //0.005e6;
  sys.GetSettings()->dfc_contact_param.sigma_tau0_s = 0.0005e6;
  sys.GetSettings()->dfc_contact_param.eta_inf_s = 50;
  sys.GetSettings()->dfc_contact_param.mi_a_s = 0.5;
  sys.GetSettings()->dfc_contact_param.t = 4e-3/2;
  sys.GetSettings()->dfc_contact_param.debug_verbose = false;
  sys.GetSettings()->solver.contact_force_model = chrono::ChSystemSMC::ContactForceModel::DFC;
  //real tolerance = 1e-3;
  //uint max_iteration = 100;
  //sys.GetSettings()->solver.max_iteration_bilateral = max_iteration;
  //sys.GetSettings()->solver.tolerance = tolerance;
  sys.GetSettings()->collision.narrowphase_algorithm = collision::ChNarrowphase::Algorithm::HYBRID;
  sys.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
  std::vector<std::shared_ptr<ChBody>> container_walls = create_container(&sys, 0.15);

  std::shared_ptr<ChBody> cover;
  //cover = add_cover_formwork_pressure(sys, 0.15);

  AddSphereLayers(36, 11, 0.007, sys, DFCParticleDistr(5e-3, 10e-3, h_layer, 2.5));
  AddSphereLayers(4, 15, 0.405, sys, DFCParticleDistr(5e-3, 10e-3, h_layer, 2.5));
  //  read_particles_VTK(sys, "OUT_VTK_AGV_Set_3/particle_time_steps_01000.vtk");
  //  read_particles_VTK_inside(sys, "OUT_VTK_DeComp_Set_0_Bahar_start/particle_time_steps_00500.vtk", 0.15);
  //  read_particles_VTK_Bahar_files(sys, "Bahar_start_system/particle_S1_coords.dat",
  //				 "Bahar_start_system/particle_S1_radius.dat", true, 1);
  //std::shared_ptr<ChBody> sphere_1, sphere_2;
  //sphere_1 = AddSphere(sys, ChVector<>(0, 0, 0.15), ChVector<>(0, 0, 0), 10);
  //sphere_2 = AddSphere(sys, ChVector<>(0.010, 0, 0.15), ChVector<>(0, 0, 0), 11);

  double simulation_time = 0;
  double time_step = 1e-06;
  // time interval for data storage expressed in simulation steps
  int save_step =  1e-3 / time_step;  
  bool switch_val = false;
  int saved_steps = 0;  // serves as index for file name generation
  bool register_data = true;
  int step_num = 0;
  bool continue_simulation = true;

  std::string reaction_forces_file = out_dir + "/wall_reaction_forces.txt";
  std::string contact_forces_file = out_dir + "/wall_contact_forces.txt";

  //std::vector<std::shared_ptr<ChBody>> list_of_particles_to_delete;
  //list_of_particles_to_delete = list_particle_for_removal(sys, 0.15);

  //  for (auto i : list_of_particles_to_delete){
  //GetLog() << "Before step dynamics" << "\n";
  //sys.DoStepDynamics(time_step);  // init system for correct delete
  //GetLog() << "Before energy calculation" << "\n";
  //ChVector<> temp_energy(print_energy_status(sys));
  //terminal_file << "Total transnational kinetic energy: " << temp_energy.x();
  //terminal_file << " Total rotational kinetic energy: " << temp_energy.y() << "\n";
  //GetLog() << "Particle for removal: " << i->GetId() << "\n"; 
  //sys.Remove(i);
  //sys.Update();
  //GetLog() << "Iteratively deleted particle: " << i->GetId() << "\n";
  // }
  //std::vector<vec2> pair_to_debug;
  //pair_to_debug.push_back(vec2(sphere_1->GetId(), sphere_2->GetId()));
  //sys.GetSettings()->dfc_contact_param.debug_contact_pairs = pair_to_debug;
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
    //if (simulation_time > 0.002 && !switch_val){
    //  cover = add_cover_formwork_pressure(sys, 0.15);
    //  switch_val = true;
    //}
    // if (cover)
    //	  write_cover_data(cover, reaction_forces_file, simulation_time);
    if (register_data) {
      if (std::fmod(step_num, save_step) == 0) {
	std::string file_name = out_dir + generate_file_name("/particle_time_steps", saved_steps);
	write_particles_VTK(sys, file_name);
	++saved_steps;
	write_wall_forces(container_walls, reaction_forces_file, simulation_time);
	if (cover)
	  write_cover_data(cover, reaction_forces_file, simulation_time);
	GetLog() << "Simulation is running. Current time step: " << simulation_time << "\n";
	terminal_file << "Simulation is running. Current time step: " << simulation_time << "\n";
      }
    }
    if (std::fmod(step_num, 1000) == 0){
      write_wall_contacts(sys, container_walls[1], contact_forces_file, simulation_time);
      ChVector<> temp_energy(print_energy_status(sys));
      terminal_file << "Total transnational kinetic energy: " << temp_energy.x();
      terminal_file << " Total rotational kinetic energy: " << temp_energy.y() << "\n";      
      if (!any_particle_inside(sys) && simulation_time > 0.35){
	continue_simulation = false;
	GetLog() << "Simulation stopped. No particles in container.";
	terminal_file << "Simulation stopped. No particles in container.";
      }
      if (simulation_time > 0.5)
#ifdef IRR
	break;
#else
	continue_simulation = false;
#endif
    }
    ++step_num;
  }

  std::vector<std::shared_ptr<ChBody>> body_list = sys.Get_bodylist();
  double aggVolFrac = calc_aggVolFrac(body_list, h_layer, containerVol);
  GetLog() << "Aggregate volume fraction is equal: " << aggVolFrac << "\n";
  terminal_file << "Aggregate volume fraction is equal: " << aggVolFrac << "\n";
  terminal_file.close();
  return 0;
}


