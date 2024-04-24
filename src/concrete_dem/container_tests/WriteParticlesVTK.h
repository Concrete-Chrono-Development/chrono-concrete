// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Bahar Ayhan, Mariusz Warzecha
// =============================================================================
//
// Header file with the implementation of VTK export function

#ifndef WRITEPARTICLESVTK_H
#define WRITEPARTICLESVTK_H

using namespace chrono;

void write_particles_VTK(ChSystemMulticoreSMC& sys, const std::string& filename) {
  // Get the number of particles
  auto body_list= sys.Get_bodylist();
  int num_particles = 0;
  
  // count spheres (for generality, as system may contain other shapes)
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ++num_particles;
  }
     
  // Create the VTK file and write the header
  std::ofstream vtk_file(filename);
  vtk_file << "# vtk DataFile Version 3.0\n";
  vtk_file << "vtk output\n";
  vtk_file << "ASCII\n";
  vtk_file << "DATASET UNSTRUCTURED_GRID\n";

  // Write the particle positions
  vtk_file << "POINTS " << num_particles << " float\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<float> pos = body->GetPos();
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

  // Write the particle radii, based on value stored in body shape
  vtk_file << "\nPOINT_DATA " << num_particles << "\n";
  vtk_file << "SCALARS radius float 1\n";
  vtk_file << "LOOKUP_TABLE default\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    vtk_file << body->GetCollisionModel()->GetShapeDimensions(0)[0] << "\n";
  }

  // Write the particle linear velocities
  vtk_file << "\nVECTORS l_velocity float\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<float> vel = body->GetPos_dt();
    vtk_file << vel.x() << " " << vel.y() << " " << vel.z() << "\n";
  }

  // Write the particle angular velocities
  vtk_file << "\nVECTORS ang_velocity float\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<float> vel = body->GetWvel_par();  // angular speed in parent coords
    vtk_file << vel.x() << " " << vel.y() << " " << vel.z() << "\n";
  }

  // Write the force acting on particle
  vtk_file << "\nVECTORS force float\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<float> force = body->GetAppliedForce();  // probably in parent (global) coords
    vtk_file << force.x() << " " << force.y() << " " << force.z() << "\n";
  }

  // Write the torque acting on particle
  vtk_file << "\nVECTORS torque float\n";
  for (auto body : body_list) {
    if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType() != 0)
      continue;
    ChVector<float> torque = body->GetAppliedTorque();  // probably in local coords
    vtk_file << torque.x() << " " << torque.y() << " " << torque.z() << "\n";
  }

  // Write the energy of particles
  vtk_file << "\nVECTORS energy float\n";
  float total_trans_kin_e = 0;
  float total_rot_kin_e = 0;
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
    vtk_file << trans_kin_e << " " << rot_kin_e << " " << trans_kin_e + rot_kin_e << "\n";
    total_trans_kin_e += trans_kin_e;
    total_rot_kin_e += rot_kin_e;
  }
  // Close the file
  vtk_file.close();
}

std::string generate_file_name(const std::string& base_name, int number){
  // convert number to string
  std::string num_str = std::to_string(number);

  // determine number of zeros to add
  int num_zeros = 5 - num_str.length();

  // add needed zeros
  std::string padded_num = std::string(num_zeros, '0') + num_str;

  return base_name + "_" + padded_num + ".vtk";
}

void read_particles_VTK(ChSystemMulticoreSMC& sys, const std::string& filename) {
  // restore system state saved in VTK file
  // open the file
  std::ifstream vtk_file(filename);
  if (!vtk_file.is_open()) {
    std::cerr << "Error: unable to open file " << filename << "\n";
  }

  // read number of particles
  int particle_number = 0;
  std::string line;
  while (std::getline(vtk_file, line)) {
    if (line.find("POINTS") != std::string::npos) {
      std::istringstream iss(line);
      std::string temp;
      while (!iss.eof()) {
	iss >> temp;  // read one word and move to next one
	if (std::stringstream(temp) >> particle_number)  // is true if temp is an int number
	  break;
      }
      break;
    }
  }

  // read particle positions
  std::vector<ChVector<>> particle_pos;
  std::string temp;
  float temp_number;
  for (int i = 0; i < particle_number; ++i){
    std::getline(vtk_file, line);
    std::istringstream iss(line);
    std::vector<float> coordinates;
    while (!iss.eof()) {
      iss >> temp;
      std::stringstream(temp) >> temp_number;
      coordinates.push_back(temp_number);
    }
    particle_pos.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
  }

  // read particle radius
  std::vector<float> particle_radiuses;
  while (std::getline(vtk_file, line)) {
    if (line.find("LOOKUP_TABLE") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	float temp_radius;
	iss >> temp_radius;
	particle_radiuses.push_back(temp_radius);
      }
      break;
    }
  }

  // read particle translation velocity
  std::vector<ChVector<>> particle_l_velocity;
  while (std::getline(vtk_file, line)) {
    if (line.find("l_velocity") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	std::vector<float> coordinates;
	while (!iss.eof()) {
	  iss >> temp;
	  std::stringstream(temp) >> temp_number;
	  coordinates.push_back(temp_number);
	}
	particle_l_velocity.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
      }
      break;
    }
  }

  // read particle angular velocity
  std::vector<ChVector<>> particle_ang_velocity;
  while (std::getline(vtk_file, line)) {
    if (line.find("ang_velocity") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	std::vector<float> coordinates;
	while (!iss.eof()) {
	  iss >> temp;
	  std::stringstream(temp) >> temp_number;
	  coordinates.push_back(temp_number);
	}
	particle_ang_velocity.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
      }
      break;
    }
  }
  vtk_file.close();

  // create particles
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);
  double density = 797 * (1 + 0.4 + 2.25);
  for (int i = 0; i < particle_number; ++i) {
    auto ball = std::shared_ptr<chrono::ChBody>(sys.NewBody());
    double mass = ((4.0 / 3.0) * 3.1415 * pow(particle_radiuses[i], 3)) * density;
    ball->SetInertiaXX((2.0 / 5.0)*mass * pow(particle_radiuses[i], 3) * 
		       chrono::ChVector<>(1, 1, 1));
    ball->SetMass(mass);
    ball->SetPos(particle_pos[i]);
    ball->SetPos_dt(particle_l_velocity[i]);
    ball->SetWvel_par(particle_ang_velocity[i]);
    ball->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(ball.get(), material, particle_radiuses[i]);
    ball->SetBodyFixed(false);
    ball->SetCollide(true);
    ball->GetCollisionModel()->BuildModel();
#ifdef IRR
    auto sphere1 = chrono_types::make_shared<ChSphereShape>(particle_radiuses[i]);
    sphere1->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
    sphere1->SetOpacity(0.4f);
    auto sphere2 = chrono_types::make_shared<ChSphereShape>(particle_radiuses[i] - 4e-3);
    sphere2->SetTexture(GetChronoDataFile("textures/rock.jpg"));
    auto ball_vis = chrono_types::make_shared<ChVisualModel>();
    ball_vis->AddShape(sphere1);
    ball_vis->AddShape(sphere2);
    ball->AddVisualModel(ball_vis);
#endif
    sys.AddBody(ball);
  }

}

void read_particles_VTK_inside(ChSystemMulticoreSMC& sys, const std::string& filename, double limit) {
  // restore system state saved in VTK file
  // open the file
  std::ifstream vtk_file(filename);
  if (!vtk_file.is_open()) {
    std::cerr << "Error: unable to open file " << filename << "\n";
  }

  // read number of particles
  int particle_number = 0;
  std::string line;
  while (std::getline(vtk_file, line)) {
    if (line.find("POINTS") != std::string::npos) {
      std::istringstream iss(line);
      std::string temp;
      while (!iss.eof()) {
	iss >> temp;  // read one word and move to next one
	if (std::stringstream(temp) >> particle_number)  // is true if temp is an int number
	  break;
      }
      break;
    }
  }

  // read particle positions
  std::vector<ChVector<>> particle_pos;
  std::string temp;
  float temp_number;
  for (int i = 0; i < particle_number; ++i){
    std::getline(vtk_file, line);
    std::istringstream iss(line);
    std::vector<float> coordinates;
    while (!iss.eof()) {
      iss >> temp;
      std::stringstream(temp) >> temp_number;
      coordinates.push_back(temp_number);
    }
    particle_pos.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
  }

  // read particle radius
  std::vector<float> particle_radiuses;
  while (std::getline(vtk_file, line)) {
    if (line.find("LOOKUP_TABLE") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	float temp_radius;
	iss >> temp_radius;
	particle_radiuses.push_back(temp_radius);
      }
      break;
    }
  }

  // read particle translation velocity
  std::vector<ChVector<>> particle_l_velocity;
  while (std::getline(vtk_file, line)) {
    if (line.find("l_velocity") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	std::vector<float> coordinates;
	while (!iss.eof()) {
	  iss >> temp;
	  std::stringstream(temp) >> temp_number;
	  coordinates.push_back(temp_number);
	}
	particle_l_velocity.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
      }
      break;
    }
  }

  // read particle angular velocity
  std::vector<ChVector<>> particle_ang_velocity;
  while (std::getline(vtk_file, line)) {
    if (line.find("ang_velocity") != std::string::npos) {
      for (int i = 0; i < particle_number; ++i) {
	std::getline(vtk_file, line);
	std::istringstream iss(line);
	std::vector<float> coordinates;
	while (!iss.eof()) {
	  iss >> temp;
	  std::stringstream(temp) >> temp_number;
	  coordinates.push_back(temp_number);
	}
	particle_ang_velocity.push_back(ChVector<>(coordinates[0], coordinates[1], coordinates[2]));
      }
      break;
    }
  }
  vtk_file.close();

  // create particles
  auto material = chrono_types::make_shared<ChMaterialSurfaceSMC>();
  material->SetYoungModulus(2.05e11);
  material->SetPoissonRatio(0.3);
  material->SetRestitution(0.5);
  material->SetFriction(0.2);

  // calculate updated density
  double density_old = 797 * (1 + 0.4 + 2.25);
  double V_l = 0;  // volume of particles in container
  for (int i =0; i < particle_number; ++i) {
    V_l += 0.75 * 3.141592653589793238462 * pow(particle_radiuses[i], 3);
  }
  double density_new = (density_old * 0.15 * 0.15 * 0.15) / V_l;
  GetLog() << "Recalculation of density. \n Old density: " << density_old << "\n";
  GetLog() << "Total volume of particles in container: " << V_l << "\n";
  GetLog() << "New density: " << density_new << " (should be smaller than old density) \n";

  for (int i = 0; i < particle_number; ++i) {
    if (particle_pos[i].z() > limit || particle_pos[i] < 0 || abs(particle_pos[i].x()) > limit/2 ||
	abs(particle_pos[i].y()) > limit/2)
      continue;
    auto ball = std::shared_ptr<chrono::ChBody>(sys.NewBody());
    double mass = ((4.0 / 3.0) * 3.141592653589793238462 * pow(particle_radiuses[i], 3)) * density_new;
    ball->SetInertiaXX((2.0 / 5.0)*mass * pow(particle_radiuses[i], 3) * 
		       chrono::ChVector<>(1, 1, 1));
    ball->SetMass(mass);
    ball->SetPos(particle_pos[i]);
    ball->SetPos_dt(particle_l_velocity[i]);
    ball->SetWvel_par(particle_ang_velocity[i]);
    ball->GetCollisionModel()->ClearModel();
    utils::AddSphereGeometry(ball.get(), material, particle_radiuses[i]);
    ball->SetBodyFixed(false);
    ball->SetLimitSpeed(true);
    ball->SetMaxSpeed(.2);
    ball->SetMaxWvel(.2);
    ball->SetCollide(true);
    ball->GetCollisionModel()->BuildModel();
#ifdef IRR
    auto sphere1 = chrono_types::make_shared<ChSphereShape>(particle_radiuses[i]);
    sphere1->SetTexture(GetChronoDataFile("textures/bluewhite.png"));
    sphere1->SetOpacity(0.4f);
    auto sphere2 = chrono_types::make_shared<ChSphereShape>(particle_radiuses[i] - 4e-3);
    sphere2->SetTexture(GetChronoDataFile("textures/rock.jpg"));
    auto ball_vis = chrono_types::make_shared<ChVisualModel>();
    ball_vis->AddShape(sphere1);
    ball_vis->AddShape(sphere2);
    ball->AddVisualModel(ball_vis);
#endif
    sys.AddBody(ball);
  }
}
#endif
