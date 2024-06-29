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
//  Implemented test includes mini slump test.
//  This simulation uses different units mm, MPa, N, s

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
// improve this include to remove full path
#include "/home/mariusz/PROJECT_CHRONO/SRC/src/concrete_dem/container_tests/MyContactReport.h"
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
// improve this include to remove full path
#include "/home/mariusz/PROJECT_CHRONO/SRC/src/concrete_dem/container_tests/WriteParticlesVTK.h"

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

std::shared_ptr<ChBody> create_floor_body(ChSystemMulticore& sys){
  auto mat = chrono_types::make_shared<chrono::ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(1e2);
  mat->SetPoissonRatio(0.3);
  mat->SetRestitution(0.0);
  mat->SetFriction(0.5);
  mat->SetAdhesion(0.0);
  auto floor_body = std::shared_ptr<ChBody>(sys.NewBody());
  floor_body->SetMass(1);
  floor_body->SetInertiaXX(ChVector<>(1, 1, 1));
  floor_body->SetPos(ChVector<>(0, -5, 0));
  floor_body->SetBodyFixed(true);
  floor_body->SetCollide(true);
  floor_body->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(floor_body.get(), mat, ChVector<>(400, 10, 400));
  floor_body->GetCollisionModel()->BuildModel();
  floor_body->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
  floor_body->GetVisualShape(0)->SetOpacity(0.2f);
  sys.AddBody(floor_body);
  return floor_body;
}

std::shared_ptr<ChBody> AddConicalContainer(ChSystemMulticoreSMC& sys,
					    std::shared_ptr<ChMaterialSurface> mat,
					    double density, double radius,
					    double height,
					    std::string current_dir) {    // create cylinder
  const std::string path = current_dir + "minislump-Cut4.obj";
  auto cylinder = std::shared_ptr<ChBody>(sys.NewBody());
  cylinder->SetMass(1);
  cylinder->SetInertiaXX(ChVector<>(1, 1, 1));
  cylinder->SetCollide(true);
  cylinder->SetBodyFixed(false);
  auto trimesh = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(path);
  cylinder->GetCollisionModel()->ClearModel();
  cylinder->GetCollisionModel()->AddTriangleMesh(mat, trimesh, true, true,
						 ChVector<>(0), ChMatrix33<>(1), 5);
  cylinder->GetCollisionModel()->BuildModel();
  auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
  trimesh_shape->SetMesh(trimesh);
  trimesh_shape->SetName("Container 1");
  cylinder->AddVisualShape(trimesh_shape);
  ChQuaternion<> q=Q_from_AngAxis(CH_C_PI_2, VECT_X); //mini slump
  cylinder->SetRot(q);
  std::cout << "COG: " << cylinder->GetFrame_COG_to_abs().GetPos() << std::endl;
  cylinder->SetPos(ChVector<>(0,59,0));
  sys.AddBody(cylinder);
  return cylinder;		
}

std::shared_ptr<ChBody> AddConicalContainer2(ChSystemMulticoreSMC& sys,
					     std::shared_ptr<ChMaterialSurface> mat,
					     double density, double radius,
					     double height,
					     std::string current_dir
					     ) {    // create cylinder 
  const std::string path = current_dir + "minislump-Cut4.obj";
  auto cylinder = std::shared_ptr<ChBody>(sys.NewBody());
  cylinder->SetMass(1);
  cylinder->SetInertiaXX(ChVector<>(1, 1, 1));
  cylinder->SetCollide(true);
  cylinder->SetBodyFixed(false);
  auto trimesh = geometry::ChTriangleMeshConnected::CreateFromWavefrontFile(path);
  cylinder->GetCollisionModel()->ClearModel();
  cylinder->GetCollisionModel()->AddTriangleMesh(mat, trimesh, true, true,
						 ChVector<>(0), ChMatrix33<>(1), 5);
  cylinder->GetCollisionModel()->BuildModel();
  auto trimesh_shape = chrono_types::make_shared<ChTriangleMeshShape>();
  trimesh_shape->SetMesh(trimesh);
  trimesh_shape->SetName("Container 1");
  cylinder->AddVisualShape(trimesh_shape);

  ChQuaternion<> q=Q_from_AngAxis(CH_C_PI_2, -VECT_X); //mini slump
  cylinder->SetRot(q);
  std::cout << "COG: " << cylinder->GetFrame_COG_to_abs().GetPos() << std::endl;
  cylinder->SetPos(ChVector<>(0,(50.8 - 22.3) + 50.8 + 3,0));
  cylinder->SetBodyFixed(true);
  sys.AddBody(cylinder); 
  return cylinder;		
}

std::shared_ptr<ChBody> create_obstacle_body(ChSystemMulticore& sys){
  auto mat = chrono_types::make_shared<chrono::ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(1e2);
  mat->SetPoissonRatio(0.3);
  mat->SetRestitution(0.0);
  mat->SetFriction(0.5);
  mat->SetAdhesion(0.0);
  auto obstacle_body = std::shared_ptr<ChBody>(sys.NewBody());
  obstacle_body->SetMass(1);
  obstacle_body->SetInertiaXX(ChVector<>(1, 1, 1));
  obstacle_body->SetPos(ChVector<>(0, 23, 0));
  obstacle_body->SetBodyFixed(true);
  obstacle_body->SetCollide(true);
  obstacle_body->GetCollisionModel()->ClearModel();
  utils::AddBoxGeometry(obstacle_body.get(), mat, ChVector<>(120, 3, 120));
  obstacle_body->GetCollisionModel()->BuildModel();
  obstacle_body->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
  obstacle_body->GetVisualShape(0)->SetOpacity(1.0f);
  obstacle_body->SetBodyFixed(true);
  obstacle_body->GetCollisionModel()->SetFamilyGroup(2);
  obstacle_body->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  sys.AddBody(obstacle_body);
  return obstacle_body;
}

void ReadDFCparticles(ChSystem& sys, std::string& data_path, std::string& file_name,
		      double rho, double h_layer
		      #ifdef IRR
		      ,std::shared_ptr<ChVisualSystemIrrlicht>& vis
		      #endif
		      ) {
  auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(1e2); //2e3
  mat->SetFriction(0.5);
  mat->SetRestitution(0.0f);
  mat->SetAdhesion(0);
  struct ParticleData {
    int n;		// no
    double x;   // x coordinate
    double y;	// y coordinate
    double z;	// z coordinate	
    double d;	// diameter 
  };	
  std::string nodesFilename=data_path+file_name+".dat";
  std::ifstream nodesFile(nodesFilename);   
  // ================================================================================
  // Particle Data File
  // ================================================================================
  //
  // Data Structure:
  // n x y z d
  //
  // ================================================================================
  std::vector<ParticleData> particles;
  ParticleData temp;
  if (nodesFile.is_open()) {  		
    std::string line; 
    while (std::getline(nodesFile, line)) {
      std::istringstream iss(line);
      if (!(iss >> temp.n >> temp.x >> temp.y >> temp.z >> temp.d )) {            
	continue;
      }
      particles.push_back(temp);
    }
  }	
  Quaternion q=Q_from_AngAxis(-CH_C_PI_2, VECT_X);
  ChVector<> t(0,0,0);
  chrono::Coordsys csys(t, q);
  for (auto particle:particles) {
    float x_pos=particle.x;
    float y_pos=particle.y;
    float z_pos=particle.z;
    ChVector<> pos={x_pos, y_pos, z_pos+3};            
    ChVector<> rotatedpos=csys.TransformLocalToParent(pos);
    std::cout<<"pos: "<<pos<<"\t";
    std::cout<<"rotated pos: "<<rotatedpos<<"\n";
    double radius=particle.d/2+h_layer;
    double mass = rho/2*4/3*pow(radius, 3)*CH_C_PI;		
    auto body = std::shared_ptr<ChBody>(sys.NewBody());;
    body->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector<>(1, 1, 1));
    body->SetMass(mass);
    body->SetPos(rotatedpos);
    body->GetCollisionModel()->ClearModel();
    body->GetCollisionModel()->AddSphere(mat, radius);
    body->GetCollisionModel()->BuildModel();
    body->SetCollide(true);
    #ifdef IRR
    auto sphereMor = chrono_types::make_shared<ChSphereShape>(radius);
    sphereMor->SetColor(ChColor(128.f/255, 128.f/255, 128.f/255));
    sphereMor->SetOpacity(0.25f);
    body->AddVisualShape(sphereMor);			
    auto sphereAgg = chrono_types::make_shared<ChSphereShape>(radius-h_layer);
    sphereAgg->SetColor(ChColor(5.f/255, 48.f/255, 173.f/255));	
    body->AddVisualShape(sphereAgg);
    vis->BindItem(body);
    #endif
    sys.AddBody(body);
    body->SetLimitSpeed(true);
    body->SetMaxSpeed(2000.);
    body->SetMaxWvel(2000.);       
  }
}



int main(int argc, char* argv[]) {
  GetLog() << "Test application for implementation of DFC model in chrono::multicore\n";
  GetLog() << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  std::string out_dir = "OUT_VTK_mini_slump_development";
  if (!filesystem::create_directory(filesystem::path(out_dir))) {
    std::cerr << "Error creating directory" << out_dir << std::endl;
    return 1;
  }
  std::string terminal_log_file = out_dir + "/Mini_slump_Set_xx_terminal_log.txt";
  std::ofstream terminal_file(terminal_log_file);
  terminal_file << "Test application for implementation of DFC model in chrono::multicore\n";
  terminal_file << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";
  std::string current_dir(argv[0]);
  int pos = current_dir.find_last_of("/\\");
  current_dir = current_dir.substr(0, pos+1);
  ChSystemMulticoreSMC sys;
  sys.Set_G_acc(ChVector<>(0, -9810, 0));

  // concrete and DFC parameters (SI units)
  double specimenVol=( CH_C_PI*pow(101.6/2,2)*50.8*2/3- CH_C_PI*pow(69.85/2,2)*50.8/3);
  double Vol=specimenVol;
  // Concrete properties
  double h_layer=2.0; //4;//3;
  double minD=2.0;
  double maxD=4.0;
  double cement=845; //620; //571; //827;
  double WtoC=0.21;  //0.3; //0.45;//0.55;
  double AtoC=1.2422;  //1.2422 : 0.9674
  double rho_c=3150;
  double rho_w=1000;
  double vair=0.03;
  double nF=0.5;	
  //
  double va=1.0-cement/rho_c-(WtoC*cement)/rho_w-vair;
  double va0=(1-pow((minD/maxD),nF))*va;	
  double Va0=va0*Vol;
  double rho=cement*(1.0+WtoC+AtoC)*1E-12;
  double targetmass=Va0*rho;
  double targetVol=specimenVol*2.5;
  //	
  std::cout<<"Volume: "<< Vol<<" targetmass: "<<targetmass<<"\n";	
  sys.GetSettings()->dfc_contact_param.E_Nm = 1.5e-2;
  sys.GetSettings()->dfc_contact_param.E_Na = 100;
  sys.GetSettings()->dfc_contact_param.h = h_layer;
  sys.GetSettings()->dfc_contact_param.alfa_a = 0.25;
  sys.GetSettings()->dfc_contact_param.beta = 0.5;
  sys.GetSettings()->dfc_contact_param.sigma_t = 1.0e-3;
  sys.GetSettings()->dfc_contact_param.sigma_tau0 = 2.0e-4;
  sys.GetSettings()->dfc_contact_param.eta_inf = 2.0e-5;
  sys.GetSettings()->dfc_contact_param.kappa_0 = 100;
  sys.GetSettings()->dfc_contact_param.n = 1;
  sys.GetSettings()->dfc_contact_param.mi_a = 0.5;
  sys.GetSettings()->dfc_contact_param.E_Nm_s = 1.5e-2;
  sys.GetSettings()->dfc_contact_param.E_Na_s = 100;
  sys.GetSettings()->dfc_contact_param.alfa_a_s = 0.25;
  sys.GetSettings()->dfc_contact_param.sigma_t_s = 1.0e-3;
  sys.GetSettings()->dfc_contact_param.sigma_tau0_s = 2.0e-4;
  sys.GetSettings()->dfc_contact_param.eta_inf_s = 50;
  sys.GetSettings()->dfc_contact_param.mi_a_s = 0.5;
  sys.GetSettings()->dfc_contact_param.t = h_layer/2;
  sys.GetSettings()->dfc_contact_param.debug_verbose = false;
  sys.GetSettings()->solver.contact_force_model = chrono::ChSystemSMC::ContactForceModel::DFC;
  sys.GetSettings()->collision.narrowphase_algorithm = collision::ChNarrowphase::Algorithm::HYBRID;

  auto floor_body = create_floor_body(sys);
  auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  mat->SetYoungModulus(1e8);
  mat->SetFriction(0.5);
  mat->SetRestitution(0.0f);
  mat->SetAdhesion(0.0);
  double cyl_radius=100;
  double cyl_height=150;
  double density=7.8E-9;
  /*
  auto cylinder_body = AddConicalContainer(sys, mat, density*1000,
					   cyl_radius, cyl_height, current_dir);
  cylinder_body->GetVisualShape(0)->SetOpacity(0.4f);
  cylinder_body->GetCollisionModel()->SetFamilyGroup(2);
  cylinder_body->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  cylinder_body->SetBodyFixed(true);
  
  auto cylinder_body2 = AddConicalContainer2(sys, mat, density,
					     cyl_radius, cyl_height, current_dir);
  cylinder_body2->GetVisualShape(0)->SetOpacity(0.4f);
  cylinder_body2->SetBodyFixed(true);
  cylinder_body2->SetCollide(true);
  auto obstacle_body = create_obstacle_body(sys);*/
  
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
  std::string data_path=current_dir;
  std::string file_name="DFCgeo000-data-particles-mini-2x4-dist";
  #ifdef IRR
  std::shared_ptr<ChVisualSystem> vis;
  auto vis_irr = chrono_types::make_shared<ChVisualSystemIrrlicht>();
  ReadDFCparticles(sys, data_path, file_name, rho, h_layer,  vis_irr);
  vis_irr->AttachSystem(&sys);
  vis_irr->SetWindowSize(800, 600);
  vis_irr->SetWindowTitle("SMC callbacks");
  vis_irr->Initialize();
  vis_irr->AddLogo();
  vis_irr->AddSkyBox();
  vis_irr->AddCamera(ChVector<>(200, 200, -300));
  vis_irr->AddTypicalLights();
  vis = vis_irr;
  #endif
  #ifndef IRR
  ReadDFCparticles(sys, data_path, file_name, rho, h_layer);
  GetLog() << "Particle Read";
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
    if (register_data) {
      if (std::fmod(step_num, save_step) == 0) {
	std::string file_name = out_dir + generate_file_name("/particle_time_steps", saved_steps);
	write_particles_VTK(sys, file_name);
	++saved_steps;
	GetLog() << "Simulation is running. Current time step: " << simulation_time << "\n";
	terminal_file << "Simulation is running. Current time step: " << simulation_time << "\n";
      }
    }
    if (std::fmod(step_num, 1000) == 0){
      ChVector<> temp_energy(print_energy_status(sys));
      terminal_file << "Total translation kinetic energy: " << temp_energy.x(); // 
      terminal_file << " Total rotational kinetic energy: " << temp_energy.y() << "\n";
    }
    if (simulation_time > 0.5)
#ifdef IRR
	break;
#else
      continue_simulation = false;
#endif
      ++step_num;
  }
  return 0;
}


