// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Test mesh collision
//
// =============================================================================

#include "chrono/ChConfig.h"
#include "chrono/assets/ChVisualShapeBox.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChIterativeSolverMulticore.h"

#include "chrono_vsg/ChVisualSystemVSG.h"

using namespace chrono;
using namespace chrono::vsg3d;

// ====================================================================================

int main(int argc, char* argv[]) {
    // ---------------------
    // Simulation parameters
    // ---------------------

    double gravity = 9.81;    // gravitational acceleration
    double time_step = 1e-4;  // integration step size

    double tolerance = 0;
    double contact_recovery_speed = 10;
    double collision_envelope = .005;

    uint max_iteration_normal = 0;
    uint max_iteration_sliding = 0;
    uint max_iteration_spinning = 100;
    uint max_iteration_bilateral = 0;

    double mesh_swept_sphere_radius = 0.005;

    // ---------------------------
    // Contact material properties
    // ---------------------------

    ChContactMethod contact_method = ChContactMethod::NSC;

    float object_friction = 0.9f;
    float object_restitution = 0.0f;
    float object_young_modulus = 2e7f;
    float object_poisson_ratio = 0.3f;
    float object_kn = 2e5;
    float object_gn = 40;
    float object_kt = 2e5;
    float object_gt = 20;

    float ground_friction = 0.9f;
    float ground_restitution = 0.01f;
    float ground_young_modulus = 1e6f;
    float ground_poisson_ratio = 0.3f;
    float ground_kn = 2e5;
    float ground_gn = 40;
    float ground_kt = 2e5;
    float ground_gt = 20;

    // ---------------------------------
    // Parameters for the falling object
    // ---------------------------------

    ChVector3d pos(0.2, 0.55, 0.2);
    ChVector3d init_vel(0, 0, 0);
    ChVector3d init_omg(0, 0, 0);

    // ---------------------------------
    // Parameters for the containing bin
    // ---------------------------------

    double width = 4;
    double length = 2;
    double thickness = 0.2;

    // Create the system
    ChSystemMulticore* sys = nullptr;
    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto sysNSC = new ChSystemMulticoreNSC();
            sysNSC->ChangeSolverType(SolverType::APGD);
            sysNSC->GetSettings()->solver.solver_mode = SolverMode::SPINNING;
            sysNSC->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
            sysNSC->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
            sysNSC->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
            sysNSC->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
            sysNSC->GetSettings()->solver.alpha = 0;
            sysNSC->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
            sys = sysNSC;
            break;
        }
        case ChContactMethod::SMC: {
            auto sysSMC = new ChSystemMulticoreSMC();
            sysSMC->GetSettings()->solver.contact_force_model = ChSystemSMC::Hertz;
            sysSMC->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::OneStep;
            sys = sysSMC;
            break;
        }
    }

    sys->SetGravitationalAcceleration(ChVector3d(0, -gravity, 0));

    // Set associated collision detection system
    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    // Set number of threads
    sys->SetNumThreads(2);

    sys->GetSettings()->solver.use_full_inertia_tensor = false;
    sys->GetSettings()->solver.tolerance = tolerance;

    sys->GetSettings()->collision.collision_envelope = collision_envelope;
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::HYBRID;
    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // Rotation Z->Y (because meshes used here assume Z up)
    ChQuaternion<> z2y = QuatFromAngleX(-CH_PI_2);

    // Create the falling object
    auto object = chrono_types::make_shared<ChBody>();
    sys->AddBody(object);

    object->SetMass(200);
    object->SetInertiaXX(40.0 * ChVector3d(1, 1, 0.2));
    object->SetPos(pos);
    object->SetRot(z2y);
    object->SetPosDt(init_vel);
    object->SetAngVelParent(init_omg);
    object->EnableCollision(true);
    object->SetFixed(false);

    std::shared_ptr<ChContactMaterial> object_mat;
    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChContactMaterialNSC>();
            matNSC->SetFriction(object_friction);
            matNSC->SetRestitution(object_restitution);
            object_mat = matNSC;
            break;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChContactMaterialSMC>();
            matSMC->SetFriction(object_friction);
            matSMC->SetRestitution(object_restitution);
            matSMC->SetYoungModulus(object_young_modulus);
            matSMC->SetPoissonRatio(object_poisson_ratio);
            matSMC->SetKn(object_kn);
            matSMC->SetGn(object_gn);
            matSMC->SetKt(object_kt);
            matSMC->SetGt(object_gt);
            object_mat = matSMC;
            break;
        }
    }

    auto trimesh =
        ChTriangleMeshConnected::CreateFromWavefrontFile(GetChronoDataFile("vehicle/hmmwv/hmmwv_tire_coarse.obj"));

    auto object_ct_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(object_mat, trimesh, false, false,
                                                                                   mesh_swept_sphere_radius);
    object->AddCollisionShape(object_ct_shape);

    auto trimesh_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
    trimesh_shape->SetMesh(trimesh);
    trimesh_shape->SetName("tire");
    object->AddVisualShape(trimesh_shape);

    // Create ground body
    auto ground = chrono_types::make_shared<ChBody>();
    sys->AddBody(ground);

    ground->SetMass(1);
    ground->SetPos(ChVector3d(0, 0, 0));
    ground->SetRot(z2y);
    ground->EnableCollision(true);
    ground->SetFixed(true);

    std::shared_ptr<ChContactMaterial> ground_mat;
    switch (contact_method) {
        case ChContactMethod::NSC: {
            auto matNSC = chrono_types::make_shared<ChContactMaterialNSC>();
            matNSC->SetFriction(ground_friction);
            matNSC->SetRestitution(ground_restitution);
            ground_mat = matNSC;
            break;
        }
        case ChContactMethod::SMC: {
            auto matSMC = chrono_types::make_shared<ChContactMaterialSMC>();
            matSMC->SetFriction(ground_friction);
            matSMC->SetRestitution(ground_restitution);
            matSMC->SetYoungModulus(ground_young_modulus);
            matSMC->SetPoissonRatio(ground_poisson_ratio);
            matSMC->SetKn(ground_kn);
            matSMC->SetGn(ground_gn);
            matSMC->SetKt(ground_kt);
            matSMC->SetGt(ground_gt);
            ground_mat = matSMC;
            break;
        }
    }

    auto ground_ct_shape = chrono_types::make_shared<ChCollisionShapeBox>(ground_mat, width, length, thickness);
    ground->AddCollisionShape(ground_ct_shape, ChFrame<>(ChVector3d(0, 0, -thickness / 2), QUNIT));

    auto box = chrono_types::make_shared<ChVisualShapeBox>(width, length, thickness);
    ground->AddVisualShape(box, ChFrame<>(ChVector3d(0, 0, -thickness / 2)));

    // Create the visualization window
    auto vis = chrono_types::make_shared<ChVisualSystemVSG>();
    vis->AttachSystem(sys);
    vis->SetWindowTitle("Mesh-mesh test");
    vis->SetCameraVertical(CameraVerticalDir::Y);
    vis->AddCamera(ChVector3d(2, 1, 2), ChVector3d(0, 0, 0));
    vis->SetWindowSize(1280, 720);
    vis->SetBackgroundColor(ChColor(0.8f, 0.85f, 0.9f));
    vis->EnableSkyBox();
    vis->SetCameraAngleDeg(40.0);
    vis->SetLightIntensity(1.0f);
    vis->SetLightDirection(1.5 * CH_PI_2, CH_PI_4);
    vis->EnableShadows();
    vis->Initialize();

    // Simulation loop
    while (vis->Run()) {
        sys->DoStepDynamics(time_step);
        vis->Render();
    }

    delete sys;
    return 0;
}
