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
//  Implemented test includes BMC test and is based on Bahar code.
//  This simulation uses different units mm, MPa, N, s

// comment the line below two switch off irrlicht visualisation
#define IRR
#include <vector>
#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono/assets/ChVisualSystem.h"
#include "chrono/utils/ChUtilsCreators.h"
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
#include "chrono/physics/ChInertiaUtils.h"
#include "chrono/assets/ChVisualShapeTriangleMesh.h"
#include "chrono/assets/ChTexture.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/physics/ChBodyEasy.h"

// improve this include to remove full path
#include "/home/mariusz/PROJECT_CHRONO/SRC/src/concrete_dem/container_tests/WriteParticlesVTK.h"
// Use the main namespace of Chrono, and other chrono namespaces

#ifdef IRR
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"
using namespace chrono::irrlicht;
#endif
using namespace chrono;

// Define parameters
double density_steel=7.8E-9;
double bmcRod_AngVel = 300.;   //*CH_2PI;
double total_mass = 5.22;
std::string bmcRod_obj = "bmcRod.obj";
std::string bmcContainer_obj = "bmcContainer.obj";

// Initial Position and Velocity of the bmcRod and the container
ChVector3d bmcRod_IniPos(0.0, 0.16, 0.05);
ChVector3d bmcContainer_IniPos(0.0, 0.0, 0.0);
ChVector3d bmcRod_IniVel(0.0, 0.0, 0.0);
ChVector3d bmcContainer_IniVel(0.0, 0.0, 0.0);

// calculate kinetic energy of all moving bodies in the system
double calculateKE(ChSystemMulticoreSMC& sys){
  double KE=0;
  for (auto body:sys.GetBodies()){
    if (body->IsFixed() || body->GetCollisionModel()->GetShapeInstance(0).first->GetType()!=0 )
      continue;
    ChVector3d velo=body->GetFrameCOMToAbs().GetPosDt();
    ChVector3d Wvelo=body->GetFrameCOMToAbs().GetAngVelLocal();
    double velMag=velo.Length();
    double WvelMag=Wvelo.Length();
    KE+=0.5*body->GetMass()*velMag*velMag + 0.5* body->GetInertiaXX()[0]*WvelMag*WvelMag;
  }	
  return KE;
};

std::shared_ptr<ChBodyAuxRef> CreateBMCrod(ChSystemMulticoreSMC& sysMBS,
					   #ifdef IRR
					   std::shared_ptr<ChVisualSystemIrrlicht> vis,
					   #endif
					   std::shared_ptr<ChBody> ground,
					   std::string current_dir) {
    std::string bmcRod_obj = (current_dir+"/bmcRod.obj").c_str();
    std::cout << bmcRod_obj << "\n";
    // Common contact material
    auto cmaterial = chrono_types::make_shared<ChContactMaterialSMC>();
    cmaterial->SetYoungModulus(1e2);
    cmaterial->SetFriction(0.5f);
    cmaterial->SetRestitution(0.0f);
    cmaterial->SetAdhesion(0);
    
    // Create the BMC rod
    auto trimesh = chrono_types::make_shared<ChTriangleMeshConnected>();
    double scale_ratio = 1;//0.05555556;
    trimesh->LoadWavefrontMesh(bmcRod_obj, false, true);
    trimesh->Transform(ChVector3d(0, 0, 0), ChMatrix33<>(scale_ratio));  // scale to a different size
    trimesh->RepairDuplicateVertexes(1e-9);  // if meshes are not watertight

    // Compute mass inertia from mesh
    double mmass;
    double mdensity = density_steel;
    ChVector3d mcog;
    ChMatrix33<> minertia;
    trimesh->ComputeMassProperties(true, mmass, mcog, minertia);
    ChMatrix33<> principal_inertia_rot;
    ChVector3d principal_I;
    ChInertiaUtils::PrincipalInertia(minertia, principal_I, principal_inertia_rot);
    mcog = ChVector3d(0.0, 0.0, 0.0);
    std::cout << " minertia " << minertia << " p iner rot " << principal_inertia_rot << std::endl;
    
    // Set the abs orientation, position and velocity
    auto bmcRod = chrono_types::make_shared<ChBodyAuxRef>();
    ChQuaternion<> bmcRod_Rot = QUNIT;

    // Set the COG coordinates to barycenter, without displacing the REF reference.
    // Make the COG frame a principal frame.
    bmcRod->SetFrameCOMToRef(ChFrame<>(mcog, principal_inertia_rot));

    // Set inertia
    bmcRod->SetPosDt(bmcRod_IniVel);
    bmcRod->SetAngVelLocal(ChVector3d(0.0, 0.0, 0.0));  // set an initial angular velocity (rad/s)

    // Set the absolute position of the body:
    bmcRod->SetFrameRefToAbs(ChFrame<>(ChVector3d(bmcRod_IniPos), ChQuaternion<>(bmcRod_Rot)));
    sysMBS.AddBody(bmcRod);

    bmcRod->SetFixed(false);
    auto bmcRod_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(cmaterial,
										trimesh,
										false,
										false,
										0.0002);
    bmcRod->AddCollisionShape(bmcRod_shape);
    bmcRod->EnableCollision(false);
    #ifdef IRR
    auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    trimesh_shape->SetMesh(trimesh);
    trimesh_shape->SetMutable(false);
    trimesh_shape->SetColor(ChColor(0.2f, 0.32f, 0.48f));
    bmcRod->AddVisualShape(trimesh_shape, bmcRod->GetFrameRefToAbs());
    vis->BindItem(bmcRod);
    #endif
    
    //re-assign mass & inertia
    bmcRod->SetMass(mmass*density_steel);
    bmcRod->SetInertia(principal_I*density_steel);
    std::cout << " principal_I : " << bmcRod->GetInertia() << " mass : "
	      << bmcRod->GetMass() << bmcRod->GetPos() << std::endl;
    return bmcRod;
}

std::shared_ptr<ChBody> AddBMCContainer(ChSystemMulticoreSMC& sys,
					std::shared_ptr<ChContactMaterial> mat,
					double density,
					std::string current_dir) {
  std::string bmcContainer_obj = (current_dir+"/bmcContainer.obj").c_str();
  std::cout << bmcContainer_obj << "\n";
  auto container = chrono_types::make_shared<ChBodyEasyMesh>(bmcContainer_obj,
							     density,
							     true,  //auto. mass, COG pos., inertia
							     #ifdef IRR
							     true,  // attach visualization asset
							     #else
							     false,
							     #endif
							     true,  // enable the collision detection
							     mat,  // surface contact material
							     1.0  // radius of 'inflating' of mesh for
							          //(more robust collision detection)
							     );
  ChQuaternion<> q=QuatFromAngleAxis(CH_PI_2, VECT_Z);
  container->SetRot(q);
  //  container->SetPos(ChVector3d(0,0,0));	
  //container->GetVisualShape(0)->GetMaterial(0)->SetOpacity(0.5);
  //container->SetFixed(true);	
  sys.AddBody(container); 
  return container;		
}

std::shared_ptr<ChBody> AddBMCLid(ChSystemMulticoreSMC& sys,
				  std::shared_ptr<ChContactMaterial> mat,
				  double density,
				  std::string current_dir) {   
  std::string bmcLid_obj = (current_dir+"/bmc_lid.obj").c_str();
  std::cout << bmcLid_obj << "\n";
  auto lid = chrono_types::make_shared<ChBodyEasyMesh>(bmcLid_obj,
						       density,
						       true,  ////auto. mass, COG pos., inertia
                                                       #ifdef IRR
						       true,  // attach visualization asset
                                                       #else
						       false,
						       #endif
						       true,  // enable the collision detection
						       mat,  // surface contact material
						       1.0  // radius of 'inflating' of mesh for
						       //(more robust collision detection)
						       );	
  //ChQuaternion<> q=QuatFromAngleAxis(CH_PI_2, VECT_Z);
  //lid->SetRot(q);
  //lid->SetPos(ChVector3d(0,60,0));
  //container->GetVisualShape(0)->GetMaterial(0)->SetOpacity(0.5);
  lid->SetFixed(true);	
  sys.AddBody(lid); 		
  return lid;	
}

// Function to write particle positions and radii to a VTK file
void WriteParticlesVTK(ChSystem& sys, const std::string& filename) {
  // Get the number of particles
  auto body_list= sys.GetBodies();
  std::vector<std::shared_ptr<ChBody>> body_list_new;
  for (auto body:body_list){
    if (body->IsFixed() || body->GetCollisionModel()->GetShapeInstance(0).first->GetType()!=0)
      continue;
    body_list_new.push_back(body);
  }			
  int num_particles = body_list_new.size();	

  // Create the VTK file and write the header
  std::ofstream vtk_file(filename);
  vtk_file << "# vtk DataFile Version 3.0\n";
  vtk_file << "vtk output\n";
  vtk_file << "ASCII\n";
  vtk_file << "DATASET UNSTRUCTURED_GRID\n";

  // Write the particle positions
  vtk_file << "POINTS " << num_particles << " float\n";
  for (int i = 0; i < num_particles; i++) {
    ChVector3<float> pos = body_list_new[i]->GetPos();
    vtk_file << pos.x() << " " << pos.y() << " " << pos.z() << "\n";
  }

  // Write the particle IDs
  vtk_file << "\nCELLS " << num_particles << " " << num_particles * 2 << "\n";
  for (int i = 0; i < num_particles; i++) {
    vtk_file << "1 " << i << "\n";
  }

  // Write the cell types
  vtk_file << "\nCELL_TYPES " << num_particles << "\n";
  for (int i = 0; i < num_particles; i++) {
    vtk_file << "1\n";
  }

  // Write the particle radii
  vtk_file << "\nPOINT_DATA " << num_particles << "\n";
  vtk_file << "SCALARS radius float 1\n";
  vtk_file << "LOOKUP_TABLE default\n";
  for (int i = 0; i < num_particles; i++) {		
    auto shape=body_list_new[i]->GetCollisionModel()->GetShapeInstance(0).first;
    auto shape_sphere = std::static_pointer_cast<ChCollisionShapeSphere>(shape);             	
    double radius=shape_sphere->GetRadius();
    vtk_file << radius << "\n";
  }

  // Write the particle velocities
  vtk_file << "\nVECTORS velocity float\n";
  for (int i = 0; i < num_particles; i++) {
    ChVector3<float> vel = body_list_new[i]->GetPosDt();;
    vtk_file << vel.x() << " " << vel.y() << " " << vel.z() << "\n";
  }
	 
  // Write the particle angular velocities
  vtk_file << "\nVECTORS angular_velocity float\n";
  for (int i = 0; i < num_particles; i++) {
    ChVector3<float> w = body_list_new[i]->GetFrameCOMToAbs().GetAngVelLocal();
    vtk_file << w.x() << " " << w.y() << " " << w.z() << "\n";
  }

  // Close the file
  vtk_file.close();
}

// Function to save BMC to Paraview VTK files
void WritebmcSolidVTK(const std::string& filename,
		      ChTriangleMeshConnected& mesh,
		      const ChFrame<>& frame) {
  std::ofstream outf;
  outf.open(filename);
  outf << "# vtk DataFile Version 2.0" << std::endl;
  outf << "VTK from simulation" << std::endl;
  outf << "ASCII" << std::endl;
  outf << "DATASET UNSTRUCTURED_GRID" << std::endl;
  outf << "POINTS " << mesh.GetCoordsVertices().size() << " "
       << "float" << std::endl;
  for (auto& v : mesh.GetCoordsVertices()) {
    auto w = frame.TransformPointLocalToParent(v);
    outf << w.x() << " " << w.y() << " " << w.z() << std::endl;
  }
  auto nf = mesh.GetIndicesVertexes().size();
  outf << "CELLS " << nf << " " << 4 * nf << std::endl;
  for (auto& f : mesh.GetIndicesVertexes()) {
    outf << "3 " << f.x() << " " << f.y() << " " << f.z() << std::endl;
  }
  outf << "CELL_TYPES " << nf << std::endl;
  for (int i = 0; i < nf; i++) {
    outf << "5 " << std::endl;
  }
  outf.close();
}

void AddFallingItems(ChSystem& sys,
		     #ifdef IRR
		     std::shared_ptr<ChVisualSystemIrrlicht> vis,
		     #endif
		     std::string& data_path,
		     std::string& file_name,
		     double h_layer,
		     double rho) {
  auto mat = chrono_types::make_shared<ChContactMaterialSMC>();
  mat->SetYoungModulus(1e2); //2e3
  mat->SetFriction(0.5);
  mat->SetRestitution(0.0f);
  mat->SetAdhesion(0);
  
  std::vector<chrono::ChVector3<float>> body_points;
  std::vector<chrono::ChVector3<float>> body_points_new;
  std::string nodesFilename=data_path+file_name+"_coords.dat";
  std::ifstream nodesFile(nodesFilename);
  float shift_in_y_direction = 27;  // read particles were below container, shift by this value
  try {
    chrono::ChVector3d pos;
    double x, y, z;       	
    std::string line;
    
    while (std::getline(nodesFile, line)) {			
      if( line[0]!='#') {				
	std::istringstream sline(line);
	sline>> x >> y >> z;
	body_points.push_back(ChVector3d(x,y,z)); 
      }
    }
		
    } catch (std::exception myerr) {
    std::cerr << "ERROR opening nodes info file: "<< std::string(nodesFilename)
	      << myerr.what() << std::endl;        
    }
  // Read radius info:
  std::string radiusFilename=data_path+file_name+"_radius.dat";
  std::ifstream radiusFile(radiusFilename);   
  std::vector<double> body_radius;	
  try {
    double rad;       	
    std::string line;
    while (std::getline(radiusFile, line)) {			
      if( line[0]!='#') {				
	std::istringstream sline(line);
	sline>> rad;
	body_radius.push_back(rad); 
      }
    }		
  } catch (std::exception myerr) {
    std::cerr << "ERROR opening radius info file: "<< std::string(nodesFilename)
	      << myerr.what() << std::endl;        
  }
		
  auto numSpheres = body_points.size();
  body_points_new.resize(numSpheres);
  int ii=0;
  int i=0;
  for (i=0; i<numSpheres; i++) {
    float x_pos=body_points[i].x();
    float y_pos=body_points[i].y();
    float z_pos=body_points[i].z();
    double y_level=y_pos-0;
    if ((y_pos>60 || y_pos<-27) || (x_pos*x_pos+z_pos*z_pos)>30.*30.)
      continue; 
    body_points_new[ii]=body_points[i];
    ChVector3d pos={x_pos, y_pos+shift_in_y_direction, z_pos};            
    ii++;
    double radius=body_radius[i];
    double mass = rho/2*4/3*pow(radius, 3)*CH_PI;           
    auto body = chrono_types::make_shared<ChBody>();
    body->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector3d(1, 1, 1));
    body->SetMass(mass);
    body->SetPos(pos);	
    auto sphere_shape = chrono_types::make_shared<ChCollisionShapeSphere>(mat, radius);
    body->AddCollisionShape(sphere_shape);    
    body->EnableCollision(true);   
    #ifdef IRR
    auto sphereMor = chrono_types::make_shared<ChVisualShapeSphere>(radius);
    sphereMor->SetColor(ChColor(128.f/255, 128.f/255, 128.f/255));
    sphereMor->SetOpacity(0.25f);
    body->AddVisualShape(sphereMor);
    auto sphereAgg = chrono_types::make_shared<ChVisualShapeSphere>(radius-h_layer);
    sphereAgg->SetColor(ChColor(5.f/255, 48.f/255, 173.f/255));	
    body->AddVisualShape(sphereAgg);
    vis->BindItem(body);
    #endif
    sys.AddBody(body);
  }
}


int main(int argc, char* argv[]) {
  std::cout << "Test application for implementation of DFC model in chrono::multicore\n";
  std::cout << "Based on open source library projectchrono.org Chrono version: "
	    << CHRONO_VERSION << "\n";
  chrono::SetChronoDataPath(CHRONO_DATA_DIR);
  std::string out_dir = "OUT_VTK_bmc_first_test";
  if (!filesystem::create_directory(filesystem::path(out_dir))) {
    std::cerr << "Error creating directory" << out_dir << std::endl;
    return 1;
  }
  std::string terminal_log_file = out_dir + "/Bmc_first_test_terminal_log.txt";
  std::ofstream terminal_file(terminal_log_file);
  terminal_file << "Test application for implementation of DFC model in chrono::multicore\n";
  terminal_file << "Based on open source library projectchrono.org Chrono version: "
	   << CHRONO_VERSION << "\n";
  std::string current_dir(argv[0]);
  int pos = current_dir.find_last_of("/\\");
  current_dir = current_dir.substr(0, pos+1);
  current_dir += "Bahar_implementation";  // adjust the data storing folder
  ChSystemMulticoreSMC sys;
  sys.SetGravitationalAcceleration(ChVector3d(0, -9810, 0));
  sys.SetMaxPenetrationRecoverySpeed(100.);
  sys.SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

  #ifdef IRR
  // Create the Irrlicht visualization system
  auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
  vis->SetWindowSize(800, 600);
  vis->SetWindowTitle("BMC");
  vis->Initialize();
  vis->AddLogo();
  vis->AddSkyBox();
  vis->AddTypicalLights();
  vis->AddCamera(ChVector3d(0, 185, -120));
  #endif

  // Create the floor:
  auto floor_mat = chrono_types::make_shared<ChContactMaterialSMC>();
  floor_mat->SetYoungModulus(1e2);
  floor_mat->SetFriction(0.5f);
  floor_mat->SetRestitution(0.0f);
  floor_mat->SetAdhesion(0);
	
  auto floorBody = chrono_types::make_shared<ChBody>();
  floorBody->SetPos(ChVector3d(0, -5, 0));
  floorBody->SetFixed(true);
  sys.Add(floorBody);
	
  double density=7.8E-9;
  auto bmcContainer =AddBMCContainer(sys, floor_mat, density, current_dir);
  #ifdef IRR
  bmcContainer->GetVisualShape(0)->SetOpacity(0.25f);
  #endif
  bmcContainer->SetFixed(true);	
  bmcContainer->GetCollisionModel()->SetFamilyGroup(2);
  bmcContainer->GetCollisionModel()->DisallowCollisionsWith(1);
  bmcContainer->EnableCollision(true);
    
  auto bmcLid =AddBMCLid(sys, floor_mat, density, current_dir);
  bmcLid->GetCollisionModel()->SetFamilyGroup(2);
  bmcLid->GetCollisionModel()->DisallowCollisionsWith(1);
  bmcLid->EnableCollision(false);

  auto bmcRod =CreateBMCrod(sys,
                            #ifdef IRR
			    vis,
			    #endif
			    floorBody, current_dir);
  bmcRod->GetCollisionModel()->SetFamilyGroup(2);
  bmcRod->GetCollisionModel()->DisallowCollisionsWith(1);
  bmcRod->EnableCollision(true);
	
  auto rodRefBody = chrono_types::make_shared<ChBody>();
  rodRefBody->SetPos(bmcRod->GetPos());
  rodRefBody->SetFixed(true);
  sys.Add(rodRefBody);

  // concrete and DFC parameters (SI units)
  double specimenVol=40.*40.*250.;
  double Vol=specimenVol*2;
  // Concrete properties
  double h_layer=2.0; //4;//3;
  double minD=2.0;
  double maxD=4.0;
  double cement=845; //620; //571; //827;
  double WtoC=0.21;  //0.3; //0.45;//0.55;
  double AtoC=1.24;  //1.2422 : 0.9674
  double rho_c=3150;
  double rho_w=1000;
  double vair=0.03;
  double nF=0.5;	
  double va=1.0-cement/rho_c-(WtoC*cement)/rho_w-vair;
  double va0=(1-pow((minD/maxD),nF))*va;	
  double Va0=va0*Vol;
  double rho=cement*(1.0+WtoC+AtoC)*1E-12;
  double targetmass=Va0*rho;
  double targetVol=specimenVol*2.5;
  std::cout<<"Volume: "<< Vol<<" targetmass: "<<targetmass<<"\n";
  terminal_file << "Volume: "<< Vol<<" targetmass: "<<targetmass<<"\n";
  sys.GetSettings()->dfc_contact_param.E_Nm = 1.5e-2;
  sys.GetSettings()->dfc_contact_param.E_Na = 100*0.001;
  sys.GetSettings()->dfc_contact_param.h = h_layer;
  sys.GetSettings()->dfc_contact_param.alfa_a = 0.25;
  sys.GetSettings()->dfc_contact_param.beta = 0.5;
  sys.GetSettings()->dfc_contact_param.sigma_t = 4.0e-4;
  sys.GetSettings()->dfc_contact_param.sigma_tau0 = 2.5e-5;
  sys.GetSettings()->dfc_contact_param.eta_inf = 12.5e-6;
  sys.GetSettings()->dfc_contact_param.kappa_0 = 100;
  sys.GetSettings()->dfc_contact_param.n = 1;
  sys.GetSettings()->dfc_contact_param.mi_a = 0.5;
  sys.GetSettings()->dfc_contact_param.E_Nm_s = 1.5e-2;
  sys.GetSettings()->dfc_contact_param.E_Na_s = 100*0.001;
  sys.GetSettings()->dfc_contact_param.alfa_a_s = 0.25;
  sys.GetSettings()->dfc_contact_param.sigma_t_s = 4.0e-4;
  sys.GetSettings()->dfc_contact_param.sigma_tau0_s = 2.5e-5;
  sys.GetSettings()->dfc_contact_param.eta_inf_s = 12.5e-6;
  sys.GetSettings()->dfc_contact_param.mi_a_s = 0.5;
  sys.GetSettings()->dfc_contact_param.t = h_layer/2;
  sys.GetSettings()->dfc_contact_param.debug_verbose = false;
  sys.GetSettings()->solver.contact_force_model = chrono::ChSystemSMC::ContactForceModel::DFC;
  sys.GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;

  //  auto container = chrono_types::make_shared<MyContactContainer>();
  //  sys.SetContactContainer(container);	// I will need to check these lines
  std::string particle_coords = current_dir + "/bmc_particle_New4_coords.dat";
  std::string particle_radii = current_dir + "/bmc_particle_New4_radius.dat";
  std::string file_name = "/bmc_particle_New4";
  AddFallingItems(sys,
		  #ifdef IRR
		  vis,
		  #endif
		  current_dir,
		  file_name,
		  h_layer,
		  rho*0.7051*0.966);
  /* this part was replaced by AddFallingitems
  read_particles_VTK_Bahar_files(sys,
				 #ifdef IRR
				 vis,
				 #endif
				 particle_coords, particle_radii, false, 1);
  */
  // Create the linear motor
  //	
  auto motor = chrono_types::make_shared<ChLinkMotorRotationSpeed>();
  motor->SetName("engine_bmcRod_bmcShaft");
  motor->Initialize(bmcRod,
		    rodRefBody,
		    ChFrame<>(bmcRod->GetPos(), 
			      chrono::QuatFromAngleAxis(CH_PI / 2.0, VECT_X))); //Rotate on Y Axis
  auto f_sequence1 = chrono_types::make_shared<ChFunctionSequence>();
  f_sequence1->InsertFunct(chrono_types::make_shared<ChFunctionConst>(0), 0.5, 1, true);
  f_sequence1->InsertFunct(chrono_types::make_shared<ChFunctionRamp>(0,1.0472/60), 60., 1, true);
  motor->SetSpeedFunction(f_sequence1);
  sys.AddLink(motor);
  
  double simulation_time = 0;
  double time_step = 1e-06;
  // time interval for data storage expressed in simulation steps
  int save_step =  1e-2 / time_step;  
  bool switch_val = false;
  int saved_steps = 0;  // serves as index for file name generation
  bool register_data = true;
  int step_num = 0;
  bool continue_simulation = true;

  FILE *fptr;
  auto outfilename = out_dir + "/history.txt";
  fptr= fopen(outfilename.c_str(), "w");
  fprintf(fptr,"time\t");
  fprintf(fptr,"Rod x\t Rod y\t Rod z\t");
  fprintf(fptr,"Motor wx\t Motor wy\t Motor wz\t");
  fprintf(fptr,"Rod wx\t Rod wy\t Rod wz\t");
  fprintf(fptr,"Motor Tx\t Motor Ty\t Motor Tz\n");
  
  #ifdef IRR
  vis->AttachSystem(&sys);
  while (vis->Run()) {
    if (std::fmod(step_num, 10000) ==0) {
      vis->BeginScene();
      vis->Render();
      std::cout << "Check if I am updating the visualistation";
      std::string screen_file;
      screen_file += std::string("/home/mariusz/PROJECT_CHRONO/BUILD_DEBUG/src/") +
	std::string("concrete_dem/test_bmc/OUT_VTK_bmc_first_test/Pictures/screen_step") +
	std::to_string(step_num) + std::string(".png");
      std::cout << screen_file << std::endl;
      vis->WriteImageToFile(screen_file);
      vis->EndScene();
    }
  #else
  while (continue_simulation) {
  #endif
    sys.DoStepDynamics(time_step);
    simulation_time += time_step;
    if (std::fmod(step_num, 1000) == 0) {
      auto current_energy = calculateKE(sys);
      std::cout << "Current time step: " << simulation_time << "\n";
      std::cout << "System total kinetic energy: " << current_energy << "\n";
      std::cout << "bmcRod angular velocity: " << bmcRod->GetAngVelParent() << "\n";
      terminal_file << "Current time step: " << simulation_time << "\n";
      terminal_file << "System total kinetic energy: " << current_energy << "\n";
      terminal_file << "bmcRod angular velocity: " << bmcRod->GetAngVelParent() << "\n";
      ChVector3d w_pos = bmcRod->GetPos();
      ChVector3d w_vel = motor->GetMotorAngleDt();
      ChVector3d angvel = bmcRod->GetAngVelLocal();
      ChVector3d torque = motor->GetMotorTorque(); 
      fprintf(fptr, " %10.6f\t ", sys.GetChTime() );
      fprintf(fptr, " %10.6f %10.6f  %10.6f\t ", w_pos.x(), w_pos.y(), w_pos.z() );
      fprintf(fptr, " %10.6f %10.6f  %10.6f\t ", w_vel.x(), w_vel.y(), w_vel.z() );
      fprintf(fptr, " %10.6f %10.6f  %10.6f\t ", angvel.x(), angvel.y(), angvel.z() );
      fprintf(fptr, " %10.6f %10.6f  %10.6f\n ", torque.x(), torque.y(), torque.z() );
    }
    if (simulation_time > 60)
      #ifdef IRR
      break;
      #else
      continue_simulation = false;
      #endif
    ++step_num;
  }
  return 0;
}
