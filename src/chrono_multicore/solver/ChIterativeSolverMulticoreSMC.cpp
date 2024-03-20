// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2016 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban, Arman Pazouki
// =============================================================================
//
// Implementation of methods specific to the multicore smooth-contact solver.
//
// These functions implement the basic time update for a multibody system using
// a penalty-based approach for including frictional contact. It is assumed that
// geometric contact information has been already computed and is available.
// The current algorithm is based on a semi-implicit Euler scheme and projection
// on the velocity manifold of the bilateral constraints.
//
// =============================================================================

//// RADU
//// When using the MultiStep tangential displacement mode, we need to calculate the
//// indeices of the shape in collision.  This conflicts with cases where contacts
//// are defined by a user custom callback as there are no shapes defined in that
//// case. Is there a solution?

#include <algorithm>
#include <stdexcept>

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChMaterialSurfaceSMC.h"
#include "chrono_multicore/solver/ChIterativeSolverMulticore.h"

#include <thrust/sort.h>

#if defined _WIN32
    #include <cstdint>
#endif

using namespace chrono;

// -----------------------------------------------------------------------------
// Main worker function for calculating contact forces. Calculates the contact
// force and torque for the contact pair identified by 'index' and stores them
// in the 'extended' output arrays. The calculated force and torque vectors are
// therefore duplicated in the output arrays, once for each body involved in the
// contact (with opposite signs for the two bodies).
// -----------------------------------------------------------------------------
void function_CalcContactForces(
    int index,                                            // index of this contact pair
    vec2* body_pairs,                                     // indices of the body pair in contact
    vec2* shape_pairs,                                    // indices of the shape pair in contact
    ChSystemSMC::ContactForceModel contact_model,         // contact force model
    ChSystemSMC::AdhesionForceModel adhesion_model,       // adhesion force model
    ChSystemSMC::TangentialDisplacementModel displ_mode,  // type of tangential displacement history
    bool use_mat_props,                                   // flag specifying how coefficients are obtained
    real char_vel,                                        // characteristic velocity (Hooke)
    real min_slip_vel,                                    // threshold tangential velocity
    real min_roll_vel,                                    // threshold rolling velocity
    real min_spin_vel,                                    // threshold spinning velocity
    real dT,                                              // integration time step
    real* body_mass,                                      // body masses (per body)
    real3* pos,                                           // body positions
    quaternion* rot,                                      // body orientations
    real* vel,                                            // body linear and angular velocities
    real3* friction,                                      // eff. coefficients of friction (per contact)
    real2* modulus,                                       // eff. elasticity and shear modulus (per contact)
    real3* adhesion,                                      // eff. adhesion paramters (per contact)
    real* cr,                                             // eff. coefficient of restitution (per contact)
    real4* smc_params,                                    // eff. SMC parameters k and g (per contact)
    real3* pt1,                                           // point on shape 1 (per contact)
    real3* pt2,                                           // point on shape 2 (per contact)
    real3* normal,                                        // contact normal (per contact)
    real* depth,                                          // penetration depth (per contact)
    real* eff_radius,                                     // effective contact radius (per contact)
    vec3* shear_neigh,                                    // neighbor list of contacting bodies and shapes (per body)
    char* shear_touch,                                    // flag if contact in neighbor list is persistent (per body)
    real3* shear_disp,                                    // accumulated shear displacement for each neighbor (per body)
    real* contact_relvel_init,                            // initial relative normal velocity per contact pair
    real* contact_duration,                               // duration of persistent contact between contact pairs
    int* ct_bid,                                          // [output] body IDs (two per contact)
    real3* ct_force,                                      // [output] body force (two per contact)
    real3* ct_torque                                      // [output] body torque (two per contact)
) {
    // Identify the two bodies in contact (global body IDs).
    int b1 = body_pairs[index].x;
    int b2 = body_pairs[index].y;

    // If the two contact shapes are actually separated, set zero forces and torques.
    if (depth[index] >= 0) {
        ct_bid[2 * index] = b1;
        ct_bid[2 * index + 1] = b2;
        ct_force[2 * index] = real3(0);
        ct_force[2 * index + 1] = real3(0);
        ct_torque[2 * index] = real3(0);
        ct_torque[2 * index + 1] = real3(0);
        return;
    }

    // Kinematic information
    // ---------------------

    // Express contact point locations in local frames
    //   s' = At * s = At * (rP - r)
    real3 pt1_loc = TransformParentToLocal(pos[b1], rot[b1], pt1[index]);
    real3 pt2_loc = TransformParentToLocal(pos[b2], rot[b2], pt2[index]);

    // Calculate velocities of the contact points (in global frame)
    //   vP = v + omg x s = v + A * (omg' x s')
    real3 v_body1 = real3(vel[b1 * 6 + 0], vel[b1 * 6 + 1], vel[b1 * 6 + 2]);
    real3 v_body2 = real3(vel[b2 * 6 + 0], vel[b2 * 6 + 1], vel[b2 * 6 + 2]);

    real3 o_body1 = real3(vel[b1 * 6 + 3], vel[b1 * 6 + 4], vel[b1 * 6 + 5]);
    real3 o_body2 = real3(vel[b2 * 6 + 3], vel[b2 * 6 + 4], vel[b2 * 6 + 5]);

    real3 vel1 = v_body1 + Rotate(Cross(o_body1, pt1_loc), rot[b1]);
    real3 vel2 = v_body2 + Rotate(Cross(o_body2, pt2_loc), rot[b2]);

    // Calculate relative velocity (in global frame)
    // Note that relvel_n_mag is a signed quantity, while relvel_t_mag is an
    // actual magnitude (always positive).
    real3 relvel = vel2 - vel1;
    real relvel_n_mag = Dot(relvel, normal[index]);
    real3 relvel_n = relvel_n_mag * normal[index];
    real3 relvel_t = relvel - relvel_n;
    real relvel_t_mag = Length(relvel_t);

    // Extract composite material properties
    // -------------------------------------

    real m_eff = body_mass[b1] * body_mass[b2] / (body_mass[b1] + body_mass[b2]);

    real mu_eff = friction[index].x;
    real muRoll_eff = friction[index].y;
    real muSpin_eff = friction[index].z;

    real E_eff = modulus[index].x;
    real G_eff = modulus[index].y;

    real adhesion_eff = adhesion[index].x;
    real adhesionMultDMT_eff = adhesion[index].y;
    real adhesionSPerko_eff = adhesion[index].z;

    real cr_eff = cr[index];

    real user_kn = smc_params[index].x;
    real user_kt = smc_params[index].y;
    real user_gn = smc_params[index].z;
    real user_gt = smc_params[index].w;

    // Contact force
    // -------------

    // All models use the following formulas for normal and tangential forces:
    //     Fn = kn * delta_n - gn * v_n
    //     Ft = kt * delta_t - gt * v_t
    // The stiffness and damping coefficients are obtained differently, based
    // on the force model and on how coefficients are specified.
    real kn = 0;
    real kt = 0;
    real gn = 0;
    real gt = 0;
    real kn_simple = 0;
    real gn_simple = 0;

    real t_contact = 0;
    real relvel_init = abs(relvel_n_mag);
    real delta_n = -depth[index];
    real3 delta_t = real3(0);

    int i;
    int contact_id = -1;
    int shear_body1 = -1;
    int shear_body2;
    int shear_shape1;
    int shear_shape2;
    bool newcontact = true;

    if (displ_mode == ChSystemSMC::TangentialDisplacementModel::OneStep) {
        delta_t = relvel_t * dT;
    } else if (displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep) {
        delta_t = relvel_t * dT;

        // Identify the two shapes in contact (global shape IDs).
        int s1 = shape_pairs[index].x;
        int s2 = shape_pairs[index].y;

        // Contact history information stored on the body with the smaller shape or else the body with larger index.
        // Currently, it is assumed that the smaller shape is on the body with larger ID. We call this body shear_body1.
        shear_body1 = std::max(b1, b2);
        shear_body2 = std::min(b1, b2);
        shear_shape1 = std::max(s1, s2);
        shear_shape2 = std::min(s1, s2);

        // Check if contact history already exists. If not, initialize new contact history.
        for (i = 0; i < max_shear; i++) {
            int ctIdUnrolled = max_shear * shear_body1 + i;
            if (shear_neigh[ctIdUnrolled].x == shear_body2 && shear_neigh[ctIdUnrolled].y == shear_shape1 &&
                shear_neigh[ctIdUnrolled].z == shear_shape2) {
                contact_duration[ctIdUnrolled] += dT;
                contact_id = i;
                newcontact = false;
                break;
            }
        }
        if (newcontact == true) {
            for (i = 0; i < max_shear; i++) {
                int ctIdUnrolled = max_shear * shear_body1 + i;
                if (shear_neigh[ctIdUnrolled].x == -1) {
                    contact_id = i;
                    shear_neigh[ctIdUnrolled].x = shear_body2;
                    shear_neigh[ctIdUnrolled].y = shear_shape1;
                    shear_neigh[ctIdUnrolled].z = shear_shape2;
                    shear_disp[ctIdUnrolled].x = 0;
                    shear_disp[ctIdUnrolled].y = 0;
                    shear_disp[ctIdUnrolled].z = 0;
                    contact_relvel_init[ctIdUnrolled] = relvel_init;
                    contact_duration[ctIdUnrolled] = 0;
                    break;
                }
            }
        }

        // Record that these two bodies are really in contact at this time.
        int ctSaveId = max_shear * shear_body1 + contact_id;
        shear_touch[ctSaveId] = true;

        // Increment stored contact history tangential (shear) displacement vector and project it onto the current
        // contact plane.
        if (shear_body1 == b1) {
            shear_disp[ctSaveId] += delta_t;
            shear_disp[ctSaveId] -= Dot(shear_disp[ctSaveId], normal[index]) * normal[index];
            delta_t = shear_disp[ctSaveId];
        } else {
            shear_disp[ctSaveId] -= delta_t;
            shear_disp[ctSaveId] -= Dot(shear_disp[ctSaveId], normal[index]) * normal[index];
            delta_t = -shear_disp[ctSaveId];
        }

        // Load the initial collision velocity and accumulated contact duration from the contact history.
        relvel_init = (contact_relvel_init[ctSaveId] < char_vel) ? char_vel : contact_relvel_init[ctSaveId];
        t_contact = contact_duration[ctSaveId];
    }

    auto eps = std::numeric_limits<double>::epsilon();

    switch (contact_model) {
        case ChSystemSMC::ContactForceModel::Hooke:
            if (use_mat_props) {
                real tmp_k = (16.0 / 15) * Sqrt(eff_radius[index]) * E_eff;
                char_vel = (displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep) ? relvel_init : char_vel;
                real v2 = char_vel * char_vel;
                real loge = (cr_eff < eps) ? Log(eps) : Log(cr_eff);
                loge = (cr_eff > 1 - eps) ? Log(1 - eps) : loge;
                real tmp_g = 1 + Pow(CH_C_PI / loge, 2);
                kn = tmp_k * Pow(m_eff * v2 / tmp_k, 1.0 / 5);
                kt = kn;
                gn = Sqrt(4 * m_eff * kn / tmp_g);
                gt = gn;
            } else {
                kn = user_kn;
                kt = user_kt;
                gn = m_eff * user_gn;
                gt = m_eff * user_gt;
            }

            kn_simple = kn;
            gn_simple = gn;

            break;

        case ChSystemSMC::ContactForceModel::Hertz:
            if (use_mat_props) {
                real sqrt_Rd = Sqrt(eff_radius[index] * delta_n);
                real Sn = 2 * E_eff * sqrt_Rd;
                real St = 8 * G_eff * sqrt_Rd;
                real loge = (cr_eff < eps) ? Log(eps) : Log(cr_eff);
                real beta = loge / Sqrt(loge * loge + CH_C_PI * CH_C_PI);
                kn = (2.0 / 3) * Sn;
                kt = St;
                gn = -2 * Sqrt(5.0 / 6) * beta * Sqrt(Sn * m_eff);
                gt = -2 * Sqrt(5.0 / 6) * beta * Sqrt(St * m_eff);
            } else {
                real tmp = eff_radius[index] * Sqrt(delta_n);
                kn = tmp * user_kn;
                kt = tmp * user_kt;
                gn = tmp * m_eff * user_gn;
                gt = tmp * m_eff * user_gt;
            }

            kn_simple = kn / Sqrt(delta_n);
            gn_simple = gn / Pow(delta_n, 1.0 / 4.0);

            break;

        case ChSystemSMC::Flores:
            if (use_mat_props) {
                real sqrt_Rd = Sqrt(eff_radius[index] * delta_n);
                real Sn = 2 * E_eff * sqrt_Rd;
                real St = 8 * G_eff * sqrt_Rd;
                cr_eff = (cr_eff < 0.01) ? 0.01 : cr_eff;
                cr_eff = (cr_eff > 1.0 - eps) ? 1.0 - eps : cr_eff;
                real loge = Log(cr_eff);
                real beta = loge / Sqrt(loge * loge + CH_C_PI * CH_C_PI);
                char_vel = (displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep) ? relvel_init : char_vel;
                kn = (2.0 / 3.0) * Sn;
                kt = (2.0 / 3.0) * St;
                gn = 8.0 * (1.0 - cr_eff) * kn * delta_n / (5.0 * cr_eff * char_vel);
                gt = -2 * Sqrt(5.0 / 6) * beta * Sqrt(St * m_eff);  // Need to multiply St by 2/3 here as well ?
            } else {
                real tmp = eff_radius[index] * Sqrt(delta_n);
                kn = tmp * user_kn;
                kt = tmp * user_kt;
                gn = tmp * m_eff * user_gn * delta_n;
                gt = tmp * m_eff * user_gt;
            }

            kn_simple = kn / Sqrt(delta_n);
            gn_simple = gn / Pow(delta_n, 3.0 / 2.0);

            break;

        case ChSystemSMC::ContactForceModel::PlainCoulomb:
            if (use_mat_props) {
                real sqrt_Rd = Sqrt(delta_n);
                real Sn = 2 * E_eff * sqrt_Rd;
                real loge = (cr_eff < eps) ? Log(eps) : Log(cr_eff);
                real beta = loge / Sqrt(loge * loge + CH_C_PI * CH_C_PI);
                kn = (2.0 / 3) * Sn;
                gn = -2 * Sqrt(5.0 / 6) * beta * Sqrt(Sn * m_eff);
            } else {
                real tmp = Sqrt(delta_n);
                kn = tmp * user_kn;
                gn = tmp * user_gn;
            }

            kn_simple = kn / Sqrt(delta_n);
            gn_simple = gn / Pow(delta_n, 1.0 / 4.0);

            kt = 0;
            gt = 0;

            {
                real forceN_mag = kn * delta_n - gn * relvel_n_mag;
                real forceT_mag = mu_eff * Tanh(5.0 * relvel_t_mag) * forceN_mag;

                // Accumulate normal and tangential forces
                real3 force = forceN_mag * normal[index];
                if (relvel_t_mag >= min_slip_vel)
                    force -= (forceT_mag / relvel_t_mag) * relvel_t;

                // Convert force into the local body frames and calculate induced torques
                real3 torque1_loc = Cross(pt1_loc, RotateT(force, rot[b1]));
                real3 torque2_loc = Cross(pt2_loc, RotateT(force, rot[b2]));

                // If the duration of the current contact is less than the durration of a typical collision,
                // do not apply friction. Rolling and spinning friction should only be applied to persistant contacts
                // Rolling and spinning friction are applied right away for critically damped or over-damped systems
                real d_coeff = gn_simple / (2.0 * m_eff * Sqrt(kn_simple / m_eff));
                if (d_coeff < 1.0) {
                    real t_collision = CH_C_PI * Sqrt(m_eff / (kn_simple * (1 - d_coeff * d_coeff)));
                    if (t_contact <= t_collision) {
                        muRoll_eff = 0.0;
                        muSpin_eff = 0.0;
                    }
                }

                // Compute some additional vales needed for the rolling and spinning friction calculations
                real3 v_rot = Rotate(Cross(o_body2, pt2_loc), rot[b2]) - Rotate(Cross(o_body1, pt1_loc), rot[b1]);
                real3 rel_o = Rotate(o_body2, rot[b2]) - Rotate(o_body1, rot[b1]);

                // Calculate rolling friction torque as M_roll = mu_r * R * (F_N x v_rot) / |v_rot| (Schwartz et al.
                // 2012)
                real3 m_roll1 = real3(0);
                real3 m_roll2 = real3(0);

                if (Length(v_rot) > min_roll_vel && muRoll_eff > eps) {
                    m_roll1 = muRoll_eff * Cross(forceN_mag * pt1_loc, RotateT(v_rot, rot[b1])) / Length(v_rot);
                    m_roll2 = muRoll_eff * Cross(forceN_mag * pt2_loc, RotateT(v_rot, rot[b2])) / Length(v_rot);
                }

                // Calculate spinning friction torque as M_spin = -mu_t * r_c * ((w_n - w_p) . F_n / |w_n - w_p|) * n
                // r_c is the radius of the circle resulting from the intersecting body surfaces (Schwartz et al. 2012)
                //
                // TODO: The spinning moment calculation is only valid for sphere-sphere collisions because of the
                // r1 and r2 terms. In order for the calculation to be valid for sphere-wall collisions, the wall
                // must be ~100x particle diameters in thickness
                real3 m_spin1 = real3(0);
                real3 m_spin2 = real3(0);

                if (Length(rel_o) > min_spin_vel && muSpin_eff > eps) {
                    real r1 = Length(pt1_loc);
                    real r2 = Length(pt2_loc);
                    real xc = (r1 * r1 - r2 * r2) / (2 * (r1 + r2 - delta_n)) + 0.5 * (r1 + r2 - delta_n);
                    real rc = r1 * r1 - xc * xc;
                    rc = (rc < eps) ? eps : Sqrt(rc);

                    m_spin1 = muSpin_eff * rc *
                              RotateT(Dot(rel_o, forceN_mag * normal[index]) * normal[index], rot[b1]) / Length(rel_o);
                    m_spin2 = muSpin_eff * rc *
                              RotateT(Dot(rel_o, forceN_mag * normal[index]) * normal[index], rot[b2]) / Length(rel_o);
                }

                // Account for adhesion
                switch (adhesion_model) {
                    case ChSystemSMC::AdhesionForceModel::Constant:
                        force -= adhesion_eff * normal[index];
                        break;
                    case ChSystemSMC::AdhesionForceModel::DMT:
                        force -= adhesionMultDMT_eff * Sqrt(eff_radius[index]) * normal[index];
                        break;
                    case ChSystemSMC::AdhesionForceModel::Perko:
                        force -= adhesionSPerko_eff * eff_radius[index] * normal[index];
                        break;
                }

                ct_bid[2 * index] = b1;
                ct_bid[2 * index + 1] = b2;
                ct_force[2 * index] = -force;
                ct_force[2 * index + 1] = force;
                ct_torque[2 * index] = -torque1_loc + m_roll1 + m_spin1;
                ct_torque[2 * index + 1] = torque2_loc - m_roll2 - m_spin2;
            }

            return;
    }

    // Calculate the the normal and tangential contact forces.
    // The normal force is a magnitude, and it will be applied along the contact
    // normal direction (negative relative to body1 & positive relative to body2).
    // The tangential force is a vector with two parts: one depends on the stored
    // contact history tangential (or shear) displacement vector delta_t, and the
    // other depends on the current relative velocity vector (for viscous damping).
    real forceN_mag = kn * delta_n - gn * relvel_n_mag;
    real3 forceT_stiff = kt * delta_t;
    real3 forceT_damp = gt * relvel_t;

    // Apply Coulomb friction law.
    // We must enforce force_T_mag <= mu_eff * |forceN_mag|.
    // If force_T_mag > mu_eff * |forceN_mag| and there is shear displacement
    // due to contact history, then the shear displacement is scaled so that
    // the tangential force will be correct if force_T_mag subsequently drops
    // below the Coulomb limit.
    //
    // TODO: This implementation currently assumes that mu_slip and mu_k are equal
    real3 forceT = forceT_stiff + forceT_damp;
    real forceT_mag = Length(forceT);
    real delta_t_mag = Length(delta_t);
    real forceT_slide = mu_eff * Abs(forceN_mag);
    if (forceT_mag > forceT_slide) {
        if (delta_t_mag > eps) {
            real ratio = forceT_slide / forceT_mag;
            forceT *= ratio;
            if (displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep) {
                delta_t = (forceT - forceT_damp) / kt;
                if (shear_body1 == b1) {
                    shear_disp[max_shear * shear_body1 + contact_id] = delta_t;
                } else {
                    shear_disp[max_shear * shear_body1 + contact_id] = -delta_t;
                }
            }
        } else {
            forceT = real3(0);
        }
    }

    // Accumulate normal and tangential forces
    real3 force = forceN_mag * normal[index] - forceT;

    // Body forces (in global frame) & torques (in local frame)
    // --------------------------------------------------------

    // Convert force into the local body frames and calculate induced torques
    //    n' = s' x F' = s' x (A*F)
    real3 torque1_loc = Cross(pt1_loc, RotateT(force, rot[b1]));
    real3 torque2_loc = Cross(pt2_loc, RotateT(force, rot[b2]));

    // If the duration of the current contact is less than the durration of a typical collision,
    // do not apply friction. Rolling and spinning friction should only be applied to persistant contacts.
    // Rolling and spinning friction are applied right away for critically damped or over-damped systems.
    real d_coeff = gn_simple / (2.0 * m_eff * Sqrt(kn_simple / m_eff));
    if (d_coeff < 1.0) {
        real t_collision = CH_C_PI * Sqrt(m_eff / (kn_simple * (1 - d_coeff * d_coeff)));
        if (t_contact <= t_collision) {
            muRoll_eff = 0.0;
            muSpin_eff = 0.0;
        }
    }

    // Compute some additional vales needed for the rolling and spinning friction calculations
    real3 v_rot = Rotate(Cross(o_body2, pt2_loc), rot[b2]) - Rotate(Cross(o_body1, pt1_loc), rot[b1]);
    real3 rel_o = Rotate(o_body2, rot[b2]) - Rotate(o_body1, rot[b1]);

    // Calculate rolling friction torque as M_roll = mu_r * R * (F_N x v_rot) / |v_rot| (Schwartz et al. 2012)
    real3 m_roll1 = real3(0);
    real3 m_roll2 = real3(0);

    if (Length(v_rot) > min_roll_vel && muRoll_eff > eps) {
        m_roll1 = muRoll_eff * Cross(forceN_mag * pt1_loc, RotateT(v_rot, rot[b1])) / Length(v_rot);
        m_roll2 = muRoll_eff * Cross(forceN_mag * pt2_loc, RotateT(v_rot, rot[b2])) / Length(v_rot);
    }

    // Calculate spinning friction torque as M_spin = -mu_t * r_c * ((w_n - w_p) . F_n / |w_n - w_p|) * n
    // r_c is the radius of the circle resulting from the intersecting body surfaces (Schwartz et al. 2012)
    //
    // TODO: The spinning moment calculation is only valid for sphere-sphere collisions because of the
    // r1 and r2 terms. In order for the calculation to be valid for sphere-wall collisions, the wall
    // must be ~100x particle diameters in thickness
    real3 m_spin1 = real3(0);
    real3 m_spin2 = real3(0);

    if (Length(rel_o) > min_spin_vel && muSpin_eff > eps) {
        real r1 = Length(pt1_loc);
        real r2 = Length(pt2_loc);
        real xc = (r1 * r1 - r2 * r2) / (2 * (r1 + r2 - delta_n)) + 0.5 * (r1 + r2 - delta_n);
        real rc = r1 * r1 - xc * xc;
        rc = (rc < eps) ? eps : Sqrt(rc);

        m_spin1 =
            muSpin_eff * rc * RotateT(Dot(rel_o, forceN_mag * normal[index]) * normal[index], rot[b1]) / Length(rel_o);
        m_spin2 =
            muSpin_eff * rc * RotateT(Dot(rel_o, forceN_mag * normal[index]) * normal[index], rot[b2]) / Length(rel_o);
    }

    // Account for adhesion
    switch (adhesion_model) {
        case ChSystemSMC::AdhesionForceModel::Constant:
            force -= adhesion_eff * normal[index];
            break;
        case ChSystemSMC::AdhesionForceModel::DMT:
            force -= adhesionMultDMT_eff * Sqrt(eff_radius[index]) * normal[index];
            break;
        case ChSystemSMC::AdhesionForceModel::Perko:
            force -= adhesionSPerko_eff * eff_radius[index] * normal[index];
            break;
    }

    // Store body forces and torques, duplicated for the two bodies.
    ct_bid[2 * index] = b1;
    ct_bid[2 * index + 1] = b2;
    ct_force[2 * index] = -force;
    ct_force[2 * index + 1] = force;
    ct_torque[2 * index] = -torque1_loc + m_roll1 + m_spin1;
    ct_torque[2 * index + 1] = torque2_loc - m_roll2 - m_spin2;
}

// -----------------------------------------------------------------------------
// Main worker function for calculating contact forces. Calculates the contact
// force and torque for the contact pair identified by 'index' and stores them
// in the 'extended' output arrays. The calculated force and torque vectors are
// therefore duplicated in the output arrays, once for each body involved in the
// contact (with opposite signs for the two bodies).
// -----------------------------------------------------------------------------
void function_CalcDFCForces(int index,               // index of this contact pair
                            vec2* body_pairs,        // indices of the body pair in contact
                            vec2* shape_pairs,       // indices of the shape pair in contact
                            real dT,                 // integration time step
                            real3* pos,              // body positions
                            quaternion* rot,         // body orientations
                            real* vel,               // body linear and angular velocities
                            real3* pt1,              // point on shape 1 (per contact)
                            real3* pt2,              // point on shape 2 (per contact)
                            real3* normal,           // contact normal (per contact)
                            real* depth,             // penetration depth (per contact)
                            real3* radiuses,          // radiuses of contacting spheres (or -1 if other shape)
                            vec3* cont_neigh,        // neighbor list of contacting bodies and shapes (per body)
                            char* cont_touch,        // flag if contact in neighbor list is persistent (per body)
                            real3* DFC_stress,       // accumulated shear displacement for each neighbor (per body)
                            real* contact_duration,  // duration of persistent contact between contact pairs
                            dfc_parameters& param,   // DFC parameters stored in one object
                            int* ct_bid,             // [output] body IDs (two per contact)
                            real3* ct_force,         // [output] body force (two per contact)
                            real3* ct_torque         // [output] body torque (two per contact)
) {
    // Identify the two bodies in contact (global body IDs).
    int b1 = body_pairs[index].x;
    int b2 = body_pairs[index].y;

    // If the two contact shapes are actually separated, set zero forces and torques.
    if (depth[index] >= 0) {
        ct_bid[2 * index] = b1;
        ct_bid[2 * index + 1] = b2;
        ct_force[2 * index] = real3(0);
        ct_force[2 * index + 1] = real3(0);
        ct_torque[2 * index] = real3(0);
        ct_torque[2 * index + 1] = real3(0);
        return;
    }

    // get proper model parameters
    /// Mortar to mortar and mortar to aggregate stiffness
    real E_Nm;
    /// Aggregate to aggregate stiffness
    real E_Na;
    /// Thickness of mortar layer around an aggregate
    real h;
    /// Normal-shear coupling parameter inside concrete
    real alfa_a;
    /// Parameter governing viscous behaviour in normal direction
    real beta;
    /// Tensile strength of mortar
    real sigma_t;
    /// Mortar shear yield stress
    real sigma_tau0;
    /// Mortar plastic viscosity
    real eta_inf;
    /// Peanalty constant
    real kappa_0;
    /// Constant defining flow (n = 1 -> Newtonian, n > 1 -> shear-thickening,
    /// n < 1 shear-thinning)
    real n;
    /// Aggregate to aggregate friction coefficient
    real mi_a;
    /// thickness of mortar layer on surfaces
    real t;
    if (radiuses[index].x != -1 && radiuses[index].y != -1) {   // both contacting shapes are spheres
        E_Nm = param.E_Nm;
        E_Na = param.E_Na;
        h = param.h;
        alfa_a = param.alfa_a;
        beta = param.beta;
        sigma_t = param.sigma_t;
        sigma_tau0 = param.sigma_tau0;
        eta_inf = param.eta_inf;
        kappa_0 = param.kappa_0;
        n = param.n;
        mi_a = param.mi_a;
        t = 0;
    } else {  // any of the contacting shapes is not a sphere
        E_Nm = param.E_Nm_s;
        E_Na = param.E_Na_s;
        h = param.h;
        alfa_a = param.alfa_a_s;
        beta = param.beta;
        sigma_t = param.sigma_t_s;
        sigma_tau0 = param.sigma_tau0_s;
        eta_inf = param.eta_inf_s;
        kappa_0 = param.kappa_0;
        n = param.n;
        mi_a = param.mi_a_s;
        t = param.t;
    }
    // Identify the two shapes in contact (global shape IDs).
    int s1 = shape_pairs[index].x;
    int s2 = shape_pairs[index].y;

    // Define body I and J
    int body_I = std::max(b1, b2);
    int body_J = std::min(b1, b2);
    int shape_body_I = std::max(s1, s2);
    int shape_body_J = std::min(s1, s2);
    real R_I, R_J;
    real3 pt_I, pt_J;
    if (body_I == b1) {
        R_I = radiuses[index].x;
        R_J = radiuses[index].y;
        pt_I = pt1[index];
        pt_J = pt2[index];
    } else {
        R_I = radiuses[index].y;
        R_J = radiuses[index].x;
        pt_I = pt2[index];
        pt_J = pt1[index];
    }

    // Kinematic information
    // ---------------------
    // Calculate distance between bodies
    real l_IJ;
    if (R_I != -1 && R_J != -1) {
        l_IJ = Length(pos[body_J] - pos[body_I]);
    } else if (R_I == -1) {
        l_IJ = R_J + depth[index] + t;  // distance of the particle center to surface (depth is minus if penetration occurs)
    } else {
        l_IJ = R_I + depth[index] + t;
    }

    // Calculate versor of normal direction, pointing from I to J
    real3 e_IJ_N_vec;
    if (R_I != -1 && R_J != -1) {
      e_IJ_N_vec = (pos[body_J] - pos[body_I]) / l_IJ;
    } else {
      e_IJ_N_vec = -normal[index];
    }
    // Calculate distance from center o body I to contact surface
    real a_I;
    if (R_I != -1 && R_J != -1) {
        a_I = (Pow(R_I, 2) - Pow(R_J, 2) + Pow(l_IJ, 2)) / (2*l_IJ);
    } else {
        a_I = l_IJ - t;
    }

    // Calculate squared radius of contact surface 
    real H_IJ = Pow(R_I, 2) - Pow(a_I, 2);
    // Calculate area of contact surface H_IJ is already squared
    real A_IJ = H_IJ * CH_C_PI;
    // Create vectors to express location of the contact area
    real3 a_I_vec, a_J_vec, a_I_vec_loc, a_J_vec_loc;
    if (R_I != -1 && R_J != -1) {
      a_I_vec = a_I * e_IJ_N_vec;
      a_I_vec_loc = RotateT(a_I_vec, rot[body_I]);  // transform a_I_vec from global to local
      a_J_vec = -(l_IJ - a_I) * e_IJ_N_vec;
      a_J_vec_loc = RotateT(a_J_vec, rot[body_J]);
    } else if (R_I == -1) {
      a_I_vec_loc = TransformParentToLocal(pos[body_I], rot[body_I], pt_I);
      a_I_vec = Rotate(a_I_vec_loc, rot[body_I]);  // transform a_I_vec from local to global
      a_J_vec = -(l_IJ - a_I) * e_IJ_N_vec;
      a_J_vec_loc = RotateT(a_J_vec, rot[body_J]);
    } else {
      a_I_vec = a_I * e_IJ_N_vec;
      a_I_vec_loc = RotateT(a_I_vec, rot[body_I]);
      a_J_vec_loc = TransformParentToLocal(pos[body_J], rot[body_J], pt_J);
      a_J_vec = Rotate(a_J_vec_loc, rot[body_J]);
    }

    // Calculate velocities of the contact points (in global frame)
    //   vP = v + omg x s = v + A * (omg' x s')
    real3 v_body_I = real3(vel[body_I * 6 + 0], vel[body_I * 6 + 1], vel[body_I * 6 + 2]);
    real3 v_body_J = real3(vel[body_J * 6 + 0], vel[body_J * 6 + 1], vel[body_J * 6 + 2]);

    real3 o_body_I = real3(vel[body_I * 6 + 3], vel[body_I * 6 + 4], vel[body_I * 6 + 5]);
    real3 o_body_J = real3(vel[body_J * 6 + 3], vel[body_J * 6 + 4], vel[body_J * 6 + 5]);
    // Calculate relative velocity vectors
    real3 vel_I, vel_J;
    // vectors o_body_I and o_body_J are given in local coordinate system
    vel_I = v_body_I + Rotate(Cross(o_body_I, a_I_vec_loc), rot[body_I]);
    vel_J = v_body_J + Rotate(Cross(o_body_J, a_J_vec_loc), rot[body_J]);
    
    real3 u_IJ_dt_vec =  vel_J - vel_I;
    real3 u_IJ_ML_dt_vec = u_IJ_dt_vec - Dot(u_IJ_dt_vec, e_IJ_N_vec) * e_IJ_N_vec;
    real3 e_IJ_ML_vec;
    real3 u_IJ_ML_dt_vec_norm;
    if (Length(u_IJ_ML_dt_vec) != 0) {
      u_IJ_ML_dt_vec_norm = u_IJ_ML_dt_vec / Length(u_IJ_ML_dt_vec);
    } else {
      u_IJ_ML_dt_vec_norm = real3(0, 0, 0);
    }

    // Calculate versor of tangent direction
    real alfa = acos(Dot(e_IJ_N_vec, real3(1, 0, 0)));  // angle between normal unit vector and global X-axis
    if (alfa != 0) {
        real3 rotation_axis = Cross(e_IJ_N_vec, real3(1, 0, 0));
        if (Length(rotation_axis) != 0){    // normalize rotation axis
            rotation_axis = rotation_axis/Length(rotation_axis);
        } else {
            rotation_axis = real3(0, 0, 0);
        }
        chrono::ChQuaternion<real> transformation_quaternion;
        transformation_quaternion.Q_from_AngAxis(alfa,
                                                 chrono::ChVector(rotation_axis.x, rotation_axis.y, rotation_axis.z));
        auto global_normal_direction = transformation_quaternion.Rotate(chrono::ChVector(1, 0, 0));
        if (global_normal_direction.x() != e_IJ_N_vec.x || global_normal_direction.y() != e_IJ_N_vec.y ||
            global_normal_direction.z() != e_IJ_N_vec.z) {
            // the difference between calculated global normal direction and vector e_IJ_N_vec results from wrong sign
            transformation_quaternion.Q_from_AngAxis(-alfa, chrono::ChVector(rotation_axis.x, rotation_axis.y, rotation_axis.z));
            global_normal_direction = transformation_quaternion.Rotate(chrono::ChVector(1, 0, 0));
        }
        auto global_tangent_direction1 = transformation_quaternion.Rotate(chrono::ChVector(0, 1, 0));
        auto global_tangent_direction2 = transformation_quaternion.Rotate(chrono::ChVector(0, 0, 1));
	real dir1_val = Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction1.x(),
						  global_tangent_direction1.y(),
						  global_tangent_direction1.z()));
	real dir2_val = Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction2.x(),
						  global_tangent_direction2.y(),
						  global_tangent_direction2.z()));
	real3 dir1_vec = real3(global_tangent_direction1.x(), global_tangent_direction1.y(),
			       global_tangent_direction1.z()) * dir1_val;
	real3 dir2_vec = real3(global_tangent_direction2.x(), global_tangent_direction2.y(),
			       global_tangent_direction2.z()) * dir2_val;
	auto temp_direction = dir1_vec + dir2_vec;
	if (Length(temp_direction) != 0) {
	  e_IJ_ML_vec = temp_direction/Length(temp_direction);
	} else {
	  e_IJ_ML_vec = real3(0, 0, 0);
	}
	/*
        if (Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction1.x(), global_tangent_direction1.y(), global_tangent_direction1.z())) == 0) {
            e_IJ_ML_vec = real3(global_tangent_direction2.x(), global_tangent_direction2.y(), global_tangent_direction2.z());
        }
        else if(Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction2.x(), global_tangent_direction2.y(), global_tangent_direction2.z())) == 0){ 
            e_IJ_ML_vec = real3(global_tangent_direction1.x(), global_tangent_direction1.y(), global_tangent_direction1.z());
        } else {
            std::cout << "Definition of tangent direction failed. Check the problem.";
	    }*/
    } else {
      auto global_tangent_direction1 = chrono::ChVector(0, 1, 0);
      auto global_tangent_direction2 = chrono::ChVector(0, 0, 1);
      real dir1_val = Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction1.x(),
						global_tangent_direction1.y(),
						global_tangent_direction1.z()));
      real dir2_val = Dot(u_IJ_ML_dt_vec, real3(global_tangent_direction2.x(),
						global_tangent_direction2.y(),
						global_tangent_direction2.z()));
      real3 dir1_vec = real3(global_tangent_direction1.x(), global_tangent_direction1.y(),
			     global_tangent_direction1.z()) * dir1_val;
      real3 dir2_vec = real3(global_tangent_direction2.x(), global_tangent_direction2.y(),
			     global_tangent_direction2.z()) * dir2_val;
      auto temp_direction = dir1_vec + dir2_vec;
      if (Length(temp_direction) != 0) {
	e_IJ_ML_vec = temp_direction/Length(temp_direction);
      } else {
	e_IJ_ML_vec = real3(0, 0, 0);
      }
    }
    
    // Calculate strain rates
    real epsilon_IJ_N_dt = Dot(u_IJ_dt_vec, e_IJ_N_vec) / l_IJ;
    real epsilon_IJ_ML_dt = Dot(u_IJ_ML_dt_vec, e_IJ_ML_vec) / l_IJ;  // previous code: Length(u_IJ_ML_dt_vec)

    // Calculate gamma0 prim
    real gamma_0_dt = sigma_tau0 / (kappa_0 * eta_inf);

    // Calculate gamma prim
    real gamma_dt = Sqrt(beta * Pow(epsilon_IJ_N_dt, 2) + Pow(epsilon_IJ_ML_dt, 2));

    // Calculate eta_gamma prim
    real eta_gamma_dt;
    if (gamma_dt <= gamma_0_dt) {
        eta_gamma_dt = kappa_0 * eta_inf;
    }
    else {
        eta_gamma_dt = eta_inf * Pow(Abs(gamma_dt), n-1);
    }

    // Calculate viscous stresses
    real sigma_N_tau = beta * eta_gamma_dt * epsilon_IJ_N_dt;
    //real sigma_ML_tau = eta_gamma_dt * epsilon_IJ_ML_dt;
    real sigma_ML_tau = eta_gamma_dt * (Length(u_IJ_ML_dt_vec) / l_IJ);  // always positive

    // Calculate epsilon_N
    real epsilon_N;
    if (R_I != -1 && R_J != -1) {
        epsilon_N = log(l_IJ / (R_I + R_J - h));
    } else if (R_I == -1) {
        epsilon_N = log(l_IJ/ (R_J + t - h/2));
    } else {
        epsilon_N = log(l_IJ/ (R_I + t - h/2));
    }

    // Calculate epsilon_a
    real epsilon_a;
    if (R_I != -1 && R_J != -1) {
      epsilon_a = -log(1 + h / (R_I + R_J - 2*h));
    } else if (R_I == -1) {
      epsilon_a = -log(1 + h / (R_J - h + t));
    } else {
      epsilon_a = -log(1 + h / (R_I - h + t));
    }

    // Check if contact history already exists. If not, initialize new contact history.
    int i;
    int contact_id = -1;
    bool newcontact = true;
    for (i = 0; i < max_shear; i++) {
        int ctIdUnrolled = max_shear * body_I + i;
        if (cont_neigh[ctIdUnrolled].x == body_J && cont_neigh[ctIdUnrolled].y == shape_body_I &&
            cont_neigh[ctIdUnrolled].z == shape_body_J) {
            contact_id = i;
            newcontact = false;
            break;
        }
    }
    if (newcontact == true) {
        for (i = 0; i < max_shear; i++) {
            int ctIdUnrolled = max_shear * body_I + i;
            if (cont_neigh[ctIdUnrolled].x == -1) {
                contact_id = i;
                cont_neigh[ctIdUnrolled].x = body_J;
                cont_neigh[ctIdUnrolled].y = shape_body_I;
                cont_neigh[ctIdUnrolled].z = shape_body_J;
		if (epsilon_N > epsilon_a) {
		  DFC_stress[ctIdUnrolled].x = epsilon_N * E_Nm;  // for mortar initial contact
		} else {
		  DFC_stress[ctIdUnrolled].x = epsilon_N * E_Na; // for aggregate initial contact
		}
                DFC_stress[ctIdUnrolled].y = 0;
                DFC_stress[ctIdUnrolled].z = 0;
                break;
                }
            }
        }
    
    // calculate stiffness stresses
    int ctSaveId = max_shear * body_I + contact_id;
    real input_sigma_N_s = DFC_stress[ctSaveId].x;
    real input_sigma_ML_s = DFC_stress[ctSaveId].y;
    real sigma_N_s;
    real sigma_ML_s;
    real delta_sigma_N_s;
    real delta_sigma_ML_s;
    if (epsilon_N >= 0) {
        if (input_sigma_N_s <= sigma_t) {
            sigma_N_s = input_sigma_N_s;
            sigma_ML_s = 0;
        }
        else {
            sigma_N_s = sigma_t;
            sigma_ML_s = 0;
        }
        delta_sigma_N_s = E_Nm * epsilon_IJ_N_dt *dT;
        delta_sigma_ML_s = 0;
    }
    else {
        if (epsilon_N > epsilon_a){
            sigma_N_s = input_sigma_N_s;
            sigma_ML_s = 0;
            delta_sigma_N_s = E_Nm * epsilon_IJ_N_dt * dT;
            delta_sigma_ML_s = 0;
        } else {
            sigma_N_s = input_sigma_N_s;
            if (abs(input_sigma_ML_s) <= abs(mi_a * sigma_N_s)) {
                sigma_ML_s = input_sigma_ML_s;
            }
            else {
                sigma_ML_s = mi_a * abs(sigma_N_s) * Sign(input_sigma_ML_s);
            }
            delta_sigma_N_s = E_Na * epsilon_IJ_N_dt * dT;
            delta_sigma_ML_s = alfa_a * E_Na * epsilon_IJ_ML_dt *dT;
            // Record that these two bodies are in aggregate to aggregate contact
            cont_touch[ctSaveId] = true;
        }
    }
    
    // Calculate contact force
    real3 contact_force = (sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
      (sigma_ML_s) * A_IJ * e_IJ_ML_vec + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm;

    // Convert force into the local body frames and calculate induced torques
    //    n' = s' x F' = s' x (A*F)
    real3 contact_torque_I;
    real3 contact_torque_J;
    contact_torque_I = Cross(a_I_vec_loc, RotateT(contact_force, rot[body_I]));
    contact_torque_J = Cross(a_J_vec_loc, RotateT(contact_force, rot[body_J]));
			     
    /*
    if (R_I != -1 and R_J != -1){
      contact_torque_I = Cross(RotateT(a_I_vec, rot[body_I]),
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_I]));
      contact_torque_J = Cross(RotateT(a_J_vec, rot[body_J]),
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_J]));
    } else if (R_I == -1) {
      contact_torque_I = Cross(a_I_vec,
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_I]));
      contact_torque_J = Cross(RotateT(a_J_vec, rot[body_J]),
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_J]));
    } else {
      contact_torque_I = Cross(RotateT(a_I_vec, rot[body_I]),
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_I]));
      contact_torque_J = Cross(a_J_vec,
			       RotateT((sigma_N_s + sigma_N_tau) * A_IJ * e_IJ_N_vec +
				       sigma_ML_s * A_IJ * e_IJ_ML_vec
				       + sigma_ML_tau * A_IJ * u_IJ_ML_dt_vec_norm, rot[body_J]));
				       } */
    // Increment stored stiffness stresses and return contact force and torque
    DFC_stress[ctSaveId].x += delta_sigma_N_s;
    DFC_stress[ctSaveId].y += delta_sigma_ML_s;
    
    // Store body forces and torques, duplicated for the two bodies.
    ct_bid[2 * index] = body_I;
    ct_bid[2 * index + 1] = body_J;
    ct_force[2 * index] = contact_force;
    ct_force[2 * index + 1] = -contact_force;
    ct_torque[2 * index] = contact_torque_I;
    ct_torque[2 * index + 1] = -contact_torque_J;

    
    // only for validation purposes
    std::ofstream DFC_strain_stress;
    DFC_strain_stress.open("DFC_strain_stress.txt", std::ios_base::app);
    DFC_strain_stress << depth[index] << ", " << epsilon_N << ", " << sigma_N_s
		      << ", " << sigma_N_tau << ", " << epsilon_IJ_ML_dt << ", "
                      << sigma_ML_s << ", " << sigma_ML_tau << ", " << DFC_stress[ctSaveId].z
		      << ", " << epsilon_IJ_N_dt << " ," << contact_torque_I.x << " ,"
		      << contact_torque_I.y << " ," << contact_torque_I.z << "\n";
    DFC_strain_stress.close();
    DFC_stress[ctSaveId].z += epsilon_IJ_N_dt * dT;
    
    return;
}

// -----------------------------------------------------------------------------
// Calculate contact forces and torques for all contact pairs.
// -----------------------------------------------------------------------------

void ChIterativeSolverMulticoreSMC::host_CalcContactForces(custom_vector<int>& ct_bid,
                                                           custom_vector<real3>& ct_force,
                                                           custom_vector<real3>& ct_torque,
                                                           custom_vector<vec2>& shape_pairs,
                                                           custom_vector<char>& shear_touch,
                                                           custom_vector<real3>& shape_radiuses) {
#pragma omp parallel for
    for (int index = 0; index < (signed)data_manager->cd_data->num_rigid_contacts; index++) {
        if (data_manager->settings.solver.contact_force_model == ChSystemSMC::ContactForceModel::DFC) {
            function_CalcDFCForces(
                index,                                             // index of this contact pair
                data_manager->cd_data->bids_rigid_rigid.data(),    // indices of the body pair in contact
                shape_pairs.data(),                                 // indices of the shape pair in contact
                data_manager->settings.step_size,                  // integration time step
                data_manager->host_data.pos_rigid.data(),          // body positions
                data_manager->host_data.rot_rigid.data(),          // body orientations
                data_manager->host_data.v.data(),                  // body linear and angular velocities
                data_manager->cd_data->cpta_rigid_rigid.data(),    // point on shape 1 (per contact)
                data_manager->cd_data->cptb_rigid_rigid.data(),    // point on shape 2 (per contact)
                data_manager->cd_data->norm_rigid_rigid.data(),    // contact normal (per contact)
                data_manager->cd_data->dpth_rigid_rigid.data(),    // penetration depth (per contact)
                shape_radiuses.data(),                             // radiuses of contacting spheres (-1 if other shape)
                data_manager->host_data.shear_neigh.data(),         // neighbor list of contacting bodies and shapes (per body)
                shear_touch.data(),                                 // flag if contact in neighbor list is persistent (per body)
                data_manager->host_data.shear_disp.data(),          // accumulated stresses (variable name left from other contact force model)
                data_manager->host_data.contact_duration.data(),  // duration of persistent contact between contact pairs
                data_manager->settings.dfc_contact_param,         // DFC contact parameters
                ct_bid.data(),                                      // [output] body IDs (two per contact)
                ct_force.data(),                                    // [output] body force (two per contact)
                ct_torque.data()                                    // [output] body torque (two per contact)
            );
        } else {
            function_CalcContactForces(
                index,                                                  // index of this contact pair
                data_manager->cd_data->bids_rigid_rigid.data(),         // indices of the body pair in contact
                shape_pairs.data(),                                     // indices of the shape pair in contact
                data_manager->settings.solver.contact_force_model,      // contact force model
                data_manager->settings.solver.adhesion_force_model,     // adhesion force model
                data_manager->settings.solver.tangential_displ_mode,    // type of tangential displacement history
                data_manager->settings.solver.use_material_properties,  // flag specifying how coefficients are obtained
                data_manager->settings.solver.characteristic_vel,       // characteristic velocity (Hooke)
                data_manager->settings.solver.min_slip_vel,             // threshold tangential velocity
                data_manager->settings.solver.min_roll_vel,             // threshold rolling velocity
                data_manager->settings.solver.min_spin_vel,             // threshold spinning velocity
                data_manager->settings.step_size,                       // integration time step
                data_manager->host_data.mass_rigid.data(),              // body masses
                data_manager->host_data.pos_rigid.data(),               // body positions
                data_manager->host_data.rot_rigid.data(),               // body orientations
                data_manager->host_data.v.data(),                       // body linear and angular velocities
                data_manager->host_data.fric_rigid_rigid.data(),        // eff. coefficients of friction (per contact)
                data_manager->host_data.modulus_rigid_rigid.data(),   // eff. elasticity and shear modulus (per contact)
                data_manager->host_data.adhesion_rigid_rigid.data(),  // eff. adhesion paramters (per contact)
                data_manager->host_data.cr_rigid_rigid.data(),        // eff. coefficient of restitution (per contact)
                data_manager->host_data.smc_rigid_rigid.data(),       // eff. SMC parameters k and g (per contact)
                data_manager->cd_data->cpta_rigid_rigid.data(),       // point on shape 1 (per contact)
                data_manager->cd_data->cptb_rigid_rigid.data(),       // point on shape 2 (per contact)
                data_manager->cd_data->norm_rigid_rigid.data(),       // contact normal (per contact)
                data_manager->cd_data->dpth_rigid_rigid.data(),       // penetration depth (per contact)
                data_manager->cd_data->erad_rigid_rigid.data(),       // effective contact radius (per contact)
                data_manager->host_data.shear_neigh.data(),  // neighbor list of contacting bodies and shapes (per body)
                shear_touch.data(),  // flag if contact in neighbor list is persistent (per body)
                data_manager->host_data.shear_disp.data(),  // accumulated shear displacement for each neighbor (per body)
                data_manager->host_data.contact_relvel_init.data(),  // initial relative normal velocity per contact pair
                data_manager->host_data.contact_duration.data(),      // duration of persistent contact between contact pairs
                ct_bid.data(),    // [output] body IDs (two per contact)
                ct_force.data(),  // [output] body force (two per contact)
                ct_torque.data()  // [output] body torque (two per contact)
            );
        }
    }
}

// -----------------------------------------------------------------------------
// Include contact impulses (linear and rotational) for all bodies that are
// involved in at least one contact. For each such body, the corresponding
// entries in the arrays 'ct_body_force' and 'ct_body_torque' contain the
// cummulative force and torque, respectively, over all contacts involving that
// body.
// -----------------------------------------------------------------------------
void ChIterativeSolverMulticoreSMC::host_AddContactForces(uint ct_body_count, const custom_vector<int>& ct_body_id) {
    const custom_vector<real3>& ct_body_force = data_manager->host_data.ct_body_force;
    const custom_vector<real3>& ct_body_torque = data_manager->host_data.ct_body_torque;

#pragma omp parallel for
    for (int index = 0; index < (signed)ct_body_count; index++) {
        real3 contact_force = data_manager->settings.step_size * ct_body_force[index];
        real3 contact_torque = data_manager->settings.step_size * ct_body_torque[index];
        data_manager->host_data.hf[ct_body_id[index] * 6 + 0] += contact_force.x;
        data_manager->host_data.hf[ct_body_id[index] * 6 + 1] += contact_force.y;
        data_manager->host_data.hf[ct_body_id[index] * 6 + 2] += contact_force.z;
        data_manager->host_data.hf[ct_body_id[index] * 6 + 3] += contact_torque.x;
        data_manager->host_data.hf[ct_body_id[index] * 6 + 4] += contact_torque.y;
        data_manager->host_data.hf[ct_body_id[index] * 6 + 5] += contact_torque.z;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void ChIterativeSolverMulticoreSMC::host_SetContactForcesMap(uint ct_body_count, const custom_vector<int>& ct_body_id) {
    custom_vector<int>& ct_body_map = data_manager->host_data.ct_body_map;

#pragma omp parallel for
    for (int index = 0; index < (signed)ct_body_count; index++) {
        ct_body_map[ct_body_id[index]] = index;
    }
}

// Binary operation for adding two-object tuples
struct sum_tuples {
    thrust::tuple<real3, real3> operator()(const thrust::tuple<real3, real3>& a,
                                           const thrust::tuple<real3, real3>& b) const {
        return thrust::tuple<real3, real3>(thrust::get<0>(a) + thrust::get<0>(b),
                                           thrust::get<1>(a) + thrust::get<1>(b));
    }
};

// -----------------------------------------------------------------------------
// Process contact information reported by the narrowphase collision detection,
// generate contact forces, and update the (linear and rotational) impulses for
// all bodies involved in at least one contact.
// -----------------------------------------------------------------------------
void ChIterativeSolverMulticoreSMC::ProcessContacts() {
    const auto num_rigid_contacts = data_manager->cd_data->num_rigid_contacts;

    // 1. Calculate contact forces and torques - per contact basis
    //    For each pair of contact shapes that overlap, we calculate and store the
    //    IDs of the two corresponding bodies and the resulting contact forces and
    //    torques on the two bodies.

    custom_vector<int> ct_bid(2 * num_rigid_contacts);
    custom_vector<real3> ct_force(2 * num_rigid_contacts);
    custom_vector<real3> ct_torque(2 * num_rigid_contacts);

    // Set up additional vectors for multi-step tangential model
    custom_vector<vec2> shape_pairs;
    custom_vector<real3> shape_radiuses;
    custom_vector<char> shear_touch;
    if (data_manager->settings.solver.tangential_displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep ||
        data_manager->settings.solver.contact_force_model == ChSystemSMC::ContactForceModel::DFC) {
        shape_pairs.resize(num_rigid_contacts);
        shape_radiuses.resize(num_rigid_contacts);
        shear_touch.resize(max_shear * data_manager->num_rigid_bodies);
        Thrust_Fill(shear_touch, false);
#pragma omp parallel for
        for (int i = 0; i < (signed)num_rigid_contacts; i++) {
            vec2 pair = I2(int(data_manager->cd_data->contact_shapeIDs[i] >> 32),
                           int(data_manager->cd_data->contact_shapeIDs[i] & 0xffffffff));
            shape_pairs[i] = pair;
            double temp_R1;
            double temp_R2;
	    int body1 = data_manager->cd_data->bids_rigid_rigid[i].x;
	    int body2 = data_manager->cd_data->bids_rigid_rigid[i].y;
            if (data_manager->cd_data->shape_data.typ_rigid[pair.x] == 0) {
	      temp_R1 = (*data_manager->body_list)[body1]->
		GetCollisionModel()->GetShapeDimensions(0)[0];
            } else {
	      temp_R1 = -1;
            }
            if (data_manager->cd_data->shape_data.typ_rigid[pair.y] == 0) {
	      temp_R2 = (*data_manager->body_list)[body2]->
		GetCollisionModel()->GetShapeDimensions(0)[0];
            } else {
	      temp_R2 = -1;
            }
            shape_radiuses[i] = (real3(temp_R1, temp_R2, 0));
        }
    }

    host_CalcContactForces(ct_bid, ct_force, ct_torque, shape_pairs, shear_touch, shape_radiuses);

    data_manager->host_data.ct_force.resize(2 * num_rigid_contacts);
    data_manager->host_data.ct_torque.resize(2 * num_rigid_contacts);
    thrust::copy(THRUST_PAR ct_force.begin(), ct_force.end(), data_manager->host_data.ct_force.begin());
    thrust::copy(THRUST_PAR ct_torque.begin(), ct_torque.end(), data_manager->host_data.ct_torque.begin());

    if (data_manager->settings.solver.tangential_displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep &&
        data_manager->settings.solver.contact_force_model != ChSystemSMC::ContactForceModel::DFC) {
#pragma omp parallel for
        for (int index = 0; index < (signed)data_manager->num_rigid_bodies; index++) {
            for (int i = 0; i < max_shear; i++) {
                if (shear_touch[max_shear * index + i] == false)
                    data_manager->host_data.shear_neigh[max_shear * index + i].x = -1;
            }
        }
    }

    if (data_manager->settings.solver.contact_force_model == ChSystemSMC::ContactForceModel::DFC) {
#pragma omp parallel for
        for (int index = 0; index < (signed)data_manager->num_rigid_bodies; index++) {
            for (int i = 0; i < max_shear; i++) {
                if (shear_touch[max_shear * index + i] == false)
                    data_manager->host_data.shear_disp[max_shear * index + i].y = 0;
            }
        }
    }

    // 2. Calculate contact forces and torques - per body basis
    //    Accumulate the contact forces and torques for all bodies that are
    //    involved in at least one contact, by reducing the contact forces and
    //    torques from all contacts these bodies are involved in. The number of
    //    bodies that experience at least one contact is 'ct_body_count'.
    thrust::sort_by_key(THRUST_PAR ct_bid.begin(), ct_bid.end(),
                        thrust::make_zip_iterator(thrust::make_tuple(ct_force.begin(), ct_torque.begin())));

    custom_vector<int> ct_body_id(data_manager->num_rigid_bodies);
    custom_vector<real3>& ct_body_force = data_manager->host_data.ct_body_force;
    custom_vector<real3>& ct_body_torque = data_manager->host_data.ct_body_torque;

    ct_body_force.resize(data_manager->num_rigid_bodies);
    ct_body_torque.resize(data_manager->num_rigid_bodies);

    // Reduce contact forces from all contacts and count bodies currently involved
    // in contact. We do this simultaneously for contact forces and torques, using
    // zip iterators.
    auto end_range = thrust::reduce_by_key(
        THRUST_PAR ct_bid.begin(), ct_bid.end(),
        thrust::make_zip_iterator(thrust::make_tuple(ct_force.begin(), ct_torque.begin())), ct_body_id.begin(),
        thrust::make_zip_iterator(thrust::make_tuple(ct_body_force.begin(), ct_body_torque.begin())),
#if defined _WIN32
        thrust::equal_to<int64_t>(), sum_tuples()  // Windows compilers require an explicit-width type
#else
        thrust::equal_to<int>(), sum_tuples()
#endif
    );

    uint ct_body_count = (uint)(end_range.first - ct_body_id.begin());

    ct_body_force.resize(ct_body_count);
    ct_body_torque.resize(ct_body_count);

    // 3. Add contact forces and torques to existing forces (impulses):
    //    For all bodies involved in a contact, update the body forces and torques
    //    (scaled by the integration time step).
    host_AddContactForces(ct_body_count, ct_body_id);

    // 4. Set up map from all bodies in the system to bodies involved in a contact.
    host_SetContactForcesMap(ct_body_count, ct_body_id);
}

void ChIterativeSolverMulticoreSMC::ComputeD() {
    uint num_constraints = data_manager->num_constraints;
    if (num_constraints <= 0) {
        return;
    }

    uint num_dof = data_manager->num_dof;
    uint nnz_bilaterals = data_manager->nnz_bilaterals;

    CompressedMatrix<real>& D_T = data_manager->host_data.D_T;
    if (D_T.capacity() > 0) {
        clear(D_T);
    }

    D_T.reserve(nnz_bilaterals);
    D_T.resize(num_constraints, num_dof, false);

    data_manager->bilateral->GenerateSparsity();
    data_manager->bilateral->Build_D();

    data_manager->host_data.D = trans(D_T);
    data_manager->host_data.M_invD = data_manager->host_data.M_inv * data_manager->host_data.D;
}

void ChIterativeSolverMulticoreSMC::ComputeE() {
    if (data_manager->num_constraints <= 0) {
        return;
    }

    data_manager->host_data.E.resize(data_manager->num_constraints);
    reset(data_manager->host_data.E);

    data_manager->bilateral->Build_E();
}

void ChIterativeSolverMulticoreSMC::ComputeR() {
    if (data_manager->num_constraints <= 0) {
        return;
    }

    data_manager->host_data.b.resize(data_manager->num_constraints);
    reset(data_manager->host_data.b);
    data_manager->bilateral->Build_b();

    data_manager->host_data.R_full =
        -data_manager->host_data.b - data_manager->host_data.D_T * data_manager->host_data.M_invk;
}

// -----------------------------------------------------------------------------
// This is the main function for advancing the system state in time. On entry,
// geometric contact information is available as calculated by the narrowphase
// collision detection. This function calculates contact forces, updates the
// generalized velocities, then enforces the velocity-level constraints for any
// bilateral (joint) constraints present in the system.
// -----------------------------------------------------------------------------
void ChIterativeSolverMulticoreSMC::RunTimeStep() {
    // This is the total number of constraints, note that there are no contacts
    data_manager->num_constraints = data_manager->num_bilaterals;
    data_manager->num_unilaterals = 0;

    // Calculate contact forces (impulses) and append them to the body forces
    data_manager->host_data.ct_body_map.resize(data_manager->num_rigid_bodies);
    Thrust_Fill(data_manager->host_data.ct_body_map, -1);

    if (data_manager->cd_data->num_rigid_contacts > 0) {
        data_manager->system_timer.start("ChIterativeSolverMulticoreSMC_ProcessContact");
        ProcessContacts();
        data_manager->system_timer.stop("ChIterativeSolverMulticoreSMC_ProcessContact");
    }

    // Generate the mass matrix and compute M_inv_k
    ComputeInvMassMatrix();

    // If there are (bilateral) constraints, calculate Lagrange multipliers.
    if (data_manager->num_constraints != 0) {
        data_manager->system_timer.start("ChIterativeSolverMulticore_Setup");

        data_manager->bilateral->Setup(data_manager);

        solver->current_iteration = 0;
        data_manager->measures.solver.total_iteration = 0;
        data_manager->measures.solver.maxd_hist.clear();            ////
        data_manager->measures.solver.maxdeltalambda_hist.clear();  ////  currently not used

        solver->Setup(data_manager);

        data_manager->system_timer.stop("ChIterativeSolverMulticore_Setup");

        // Set the initial guess for the iterative solver to zero.
        data_manager->host_data.gamma.resize(data_manager->num_constraints);
        data_manager->host_data.gamma.reset();

        // Compute the jacobian matrix, the compliance matrix and the right hand side
        data_manager->system_timer.start("ChIterativeSolverMulticore_Matrices");
        ComputeD();
        ComputeE();
        ComputeR();
        data_manager->system_timer.stop("ChIterativeSolverMulticore_Matrices");

        ShurProductBilateral.Setup(data_manager);

        bilateral_solver->Setup(data_manager);

        // Solve for the Lagrange multipliers associated with bilateral constraints.
        PerformStabilization();
    }

    // Update velocity (linear and angular)
    ComputeImpulses();

    for (int i = 0; i < data_manager->measures.solver.maxd_hist.size(); i++) {
        AtIterationEnd(data_manager->measures.solver.maxd_hist[i], data_manager->measures.solver.maxdeltalambda_hist[i],
                       i);
    }
    m_iterations = (int)data_manager->measures.solver.maxd_hist.size();
}

void ChIterativeSolverMulticoreSMC::ComputeImpulses() {
    DynamicVector<real>& v = data_manager->host_data.v;
    const DynamicVector<real>& M_invk = data_manager->host_data.M_invk;
    const DynamicVector<real>& gamma = data_manager->host_data.gamma;

    uint num_unilaterals = data_manager->num_unilaterals;
    uint num_bilaterals = data_manager->num_bilaterals;

    if (data_manager->num_constraints > 0) {
        ConstSubVectorType gamma_b = blaze::subvector(gamma, num_unilaterals, num_bilaterals);
        v = M_invk + data_manager->host_data.M_invD * gamma_b;
    } else {
        v = M_invk;
    }
}
