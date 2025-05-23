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
// Authors: Alessandro Tasora, Radu Serban
// =============================================================================

//// RADU:  this was always hard-coded for NSC.  Look into allowing either contact method

#ifndef CHCONVEYOR_H
#define CHCONVEYOR_H

#include <cmath>

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChLinkLock.h"

namespace chrono {

/// Class for conveyor belt.
/// A conveyor belt is approximated by a box collision shape, where the upper surface
/// has continuous motion in X direction. No cylindrical rounding is used at the ends.
class ChApi ChConveyor : public ChPhysicsItem {
  public:
    /// Build a conveyor belt, with motion along x axis.
    ChConveyor(double xlength = 1, double ythick = 0.1, double zwidth = 0.5);
    ChConveyor(const ChConveyor& other);
    ~ChConveyor();

    /// "Virtual" copy constructor (covariant return type).
    virtual ChConveyor* Clone() const override { return new ChConveyor(*this); }

    /// Set the pointer to the parent ChSystem().
    virtual void SetSystem(ChSystem* m_system) override;

    /// Set the speed of the conveyor belt (upper part, X direction).
    void SetConveyorSpeed(double mspeed) { conveyor_speed = mspeed; }

    /// Get the speed of the conveyor belt (upper part, X direction).
    double GetConveyorSpeed() { return conveyor_speed; }

    /// Access the internal body used as the truss of the moving belt.
    ChBody* GetTruss() const { return conveyor_truss; }

    /// Access the internal body used as the moving belt (a plate with const.vel.).
    ChBody* GetPlate() const { return conveyor_plate; }

    // Shortcuts for ChBody-like transformations.
    // These, and others, can also be done as my_conveyor->GetTruss()->Ch***(...).
    void SetFixed(bool mev) { GetTruss()->SetFixed(mev); }
    bool IsFixed() { return GetTruss()->IsFixed(); }

    const ChCoordsys<>& GetCoordsys() const { return GetTruss()->GetCoordsys(); }
    const ChVector3d& GetPos() const { return GetTruss()->GetPos(); }
    const ChQuaternion<>& GetRot() const { return GetTruss()->GetRot(); }
    void SetCoordsys(const ChCoordsys<>& mcoord) { return GetTruss()->SetCoordsys(mcoord); }
    void SetCoordsys(const ChVector3d& mv, const ChQuaternion<>& mq) { GetTruss()->SetCoordsys(mv, mq); }
    void SetRot(const ChQuaternion<>& mrot) { GetTruss()->SetRot(mrot); }
    void SetPos(const ChVector3d& mpos) { GetTruss()->SetPos(mpos); }

    /// Set the material surface properties by passing a ChContactMaterialNSC or ChContactMaterialSMC object.
    void SetMaterialSurface(std::shared_ptr<ChContactMaterial> mat) { conveyor_mat = mat; }

    /// Access the material surface properties of the conveyor belt (shortcut).
    std::shared_ptr<ChContactMaterial> GetMaterialSurface() const { return conveyor_mat; }

    /// Number of coordinates: this contains an auxiliary body, so it is 14 (with quaternions for rotations).
    virtual unsigned int GetNumCoordsPosLevel() override { return 7 + 7; }

    /// Number of coordinates of the particle cluster (for two bodies).
    virtual unsigned int GetNumCoordsVelLevel() override { return 6 + 6; }

    /// Get the number of scalar constraints. In this case, a lock constraint is embedded.
    virtual unsigned int GetNumConstraintsBilateral() override { return 6; }

    virtual bool IsCollisionEnabled() const override { return true; }
    virtual void AddCollisionModelsToSystem(ChCollisionSystem* coll_sys) const override;
    virtual void RemoveCollisionModelsFromSystem(ChCollisionSystem* coll_sys) const override;
    virtual void SyncCollisionModels() override;

    /// Method to allow serialization of transient data to archives.
    virtual void ArchiveOut(ChArchiveOut& archive_out) override;

    /// Method to allow deserialization of transient data from archives.
    virtual void ArchiveIn(ChArchiveIn& archive_in) override;

  private:
    double conveyor_speed;                            ///< speed of conveyor, along the X direction of the box.
    ChLinkLockLock* internal_link;                    ///< link between this body and conveyor plate
    ChBody* conveyor_truss;                           ///< used for the conveyor truss
    ChBody* conveyor_plate;                           ///< used for the conveyor plate
    std::shared_ptr<ChContactMaterial> conveyor_mat;  ///< surface contact material for the conveyor plate

    // Solver and integrator interface functions

    virtual void Update(double time, bool update_assets) override;

    virtual void IntStateGather(const unsigned int off_x,
                                ChState& x,
                                const unsigned int off_v,
                                ChStateDelta& v,
                                double& T) override;
    virtual void IntStateScatter(const unsigned int off_x,
                                 const ChState& x,
                                 const unsigned int off_v,
                                 const ChStateDelta& v,
                                 const double T,
                                 bool full_update) override;
    virtual void IntStateGatherAcceleration(const unsigned int off_a, ChStateDelta& a) override;
    virtual void IntStateScatterAcceleration(const unsigned int off_a, const ChStateDelta& a) override;
    virtual void IntStateGatherReactions(const unsigned int off_L, ChVectorDynamic<>& L) override;
    virtual void IntStateScatterReactions(const unsigned int off_L, const ChVectorDynamic<>& L) override;
    virtual void IntStateIncrement(const unsigned int off_x,
                                   ChState& x_new,
                                   const ChState& x,
                                   const unsigned int off_v,
                                   const ChStateDelta& Dv) override;
    virtual void IntStateGetIncrement(const unsigned int off_x,
                                      const ChState& x_new,
                                      const ChState& x,
                                      const unsigned int off_v,
                                      ChStateDelta& Dv) override;
    virtual void IntLoadResidual_F(const unsigned int off, ChVectorDynamic<>& R, const double c) override;
    virtual void IntLoadResidual_Mv(const unsigned int off,
                                    ChVectorDynamic<>& R,
                                    const ChVectorDynamic<>& w,
                                    const double c) override;
    virtual void IntLoadLumpedMass_Md(const unsigned int off,
                                      ChVectorDynamic<>& Md,
                                      double& err,
                                      const double c) override;
    virtual void IntLoadResidual_CqL(const unsigned int off_L,
                                     ChVectorDynamic<>& R,
                                     const ChVectorDynamic<>& L,
                                     const double c) override;
    virtual void IntLoadConstraint_C(const unsigned int off,
                                     ChVectorDynamic<>& Qc,
                                     const double c,
                                     bool do_clamp,
                                     double recovery_clamp) override;
    virtual void IntLoadConstraint_Ct(const unsigned int off, ChVectorDynamic<>& Qc, const double c) override;
    virtual void IntToDescriptor(const unsigned int off_v,
                                 const ChStateDelta& v,
                                 const ChVectorDynamic<>& R,
                                 const unsigned int off_L,
                                 const ChVectorDynamic<>& L,
                                 const ChVectorDynamic<>& Qc) override;
    virtual void IntFromDescriptor(const unsigned int off_v,
                                   ChStateDelta& v,
                                   const unsigned int off_L,
                                   ChVectorDynamic<>& L) override;

    virtual void InjectVariables(ChSystemDescriptor& descriptor) override;
    virtual void VariablesFbReset() override;
    virtual void VariablesFbLoadForces(double factor = 1) override;
    virtual void VariablesQbLoadSpeed() override;
    virtual void VariablesFbIncrementMq() override;
    virtual void VariablesQbSetSpeed(double step = 0) override;
    virtual void VariablesQbIncrementPosition(double step) override;

    virtual void InjectConstraints(ChSystemDescriptor& descriptor) override;
    virtual void ConstraintsBiReset() override;
    virtual void ConstraintsBiLoad_C(double factor = 1, double recovery_clamp = 0.1, bool do_clamp = false) override;
    virtual void ConstraintsBiLoad_Ct(double factor = 1) override;
    virtual void ConstraintsBiLoad_Qc(double factor = 1) override;
    virtual void LoadConstraintJacobians() override;
    virtual void ConstraintsFetch_react(double factor = 1) override;
};

CH_CLASS_VERSION(ChConveyor, 0)

}  // end namespace chrono

#endif
