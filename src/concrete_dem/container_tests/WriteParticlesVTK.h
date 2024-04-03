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

#endif
