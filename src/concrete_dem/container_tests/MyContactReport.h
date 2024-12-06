// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Mariusz Warzecha
// =============================================================================
//
// Header file with the implementation of the OnReportContact() function stored in
// ReportContactCallback class, for each contact the data is stored in vector

#ifndef MYCONTACTREPORT_H
#define MYCONTACTREPORT_H

using namespace chrono;

class MyContactReport : public chrono::ChContactContainer::ReportContactCallback {
  struct CollisionData {
    ChVector3d pA; ///< contact pA
    ChVector3d pB; ///< contact pB
    ChMatrix33<> plane_coord;  ///< contact plane coordsystem (A column 'X' is contact normal)
    double distance;  ///< contact distance
    double eff_radius;  ///< effective radius of curvature at contact
    ChVector3d react_forces; ///< react.forces (if already computed). In coordsys 'plane_coord'
    ChVector3d react_torques;  ///< react.torques, if rolling friction (if already computed).
    ChContactable* contactobjA;  ///< model A (could be nullptr for some containers)
    ChContactable* contactobjB;  ///< model B (could be nullptr for some containers)
    CollisionData(
		  const ChVector3d& fpA,
		  const ChVector3d& fpB,
		  const ChMatrix33<>& fplane_coord,
		  const double& fdistance,
		  const double& feff_radius,
		  const ChVector3d& freact_forces,
		  const ChVector3d& freact_torques,
		  ChContactable* fcontactobjA,
		  ChContactable* fcontactobjB
        ) {
      pA = fpA;
      pB = fpB;
      plane_coord = fplane_coord;
      distance = fdistance;
      eff_radius = feff_radius;
      react_forces = freact_forces;
      react_torques = freact_torques;
      contactobjA = fcontactobjA;
      contactobjB = fcontactobjB;
    };
  };
  public:
  std::vector<CollisionData> VectorOfCollisionData;

  bool OnReportContact(
		       const ChVector3d& pA,
		       const ChVector3d& pB,
		       const ChMatrix33<>& plane_coord,
		       const double& distance,
		       const double& eff_radius,
		       const ChVector3d& react_forces,
		       const ChVector3d& react_torques,
		       ChContactable* contactobjA,
		       ChContactable* contactobjB
		       ) override {
    CollisionData TempCollisionData(pA, pB, plane_coord, distance, eff_radius, react_forces,
				    react_torques, contactobjA, contactobjB);
    VectorOfCollisionData.push_back(TempCollisionData);
    return 1;
  };
};

#endif
