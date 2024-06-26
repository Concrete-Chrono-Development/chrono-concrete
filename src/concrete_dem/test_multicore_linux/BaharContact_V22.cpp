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
// Authors: Bahar Ayhan
// =============================================================================
//
//  Contact force model for discrete fresh concrete
//  Details about the model can be found from the paper of "Ramyar, Elham, and Gianluca Cusatis"
//  Ramyar, Elham, and Gianluca Cusatis. "Discrete Fresh Concrete Model for Simulation of Ordinary, 
//                                        Self-Consolidating, and Printable Concrete Flow." 
//											Journal of Engineering Mechanics 148.2 (2022): 04021142.
//
//			    Material Parameters
//		 	  	float ENm		: mortar to mortar and mortar to aggregate stiffness
//				float ENa		: aggregate to aggregate stiffness
//				float h			: thickness off mortar layer around an aggregate 		
//				float alpha		: parameter for compressive contact  
//				float beta		: parameter for tensile contact
//				float np 		: np=1 newtonian fluid, np<1 shear-thinning and np>1 shear-thickening 
//				float sgmTmax 	: tensile strength of mortar
//				float sgmTau0 	: shear yield stress
//				float kappa0  	: is a constant 
//				float eta_inf 	: mortar plastic viscosity
// =============================================================================

#include "ChSystemSMC.h"
#include "ChContactSMC.h"
#include "ChContactContainerSMC.h"

#include "chrono/core/ChDistribution.h"

#include "chrono/collision/ChCollisionSystem.h"
#include "chrono/collision/ChCollisionSystemBullet.h" 
 

#include "chrono_thirdparty/rapidjson/prettywriter.h"
#include "chrono_thirdparty/rapidjson/stringbuffer.h"
#include "chrono_thirdparty/filesystem/path.h"



using namespace chrono;



// -----------------------------------------------------------------------------
// Class for overriding the default SMC contact force calculation
// -----------------------------------------------------------------------------
class MyContactForce : public ChSystemSMC::ChContactForceSMC {
  public:
    // Demonstration only.	
	
    virtual ChVector<> CalculateForce(
        const ChSystemSMC& sys,             ///< containing sys
        const ChVector<>& normal_dir,       ///< normal contact direction (expressed in global frame)
        const ChVector<>& p1,               ///< most penetrated point on obj1 (expressed in global frame)
        const ChVector<>& p2,               ///< most penetrated point on obj2 (expressed in global frame)
        const ChVector<>& vel1,             ///< velocity of contact point on obj1 (expressed in global frame)
        const ChVector<>& vel2,             ///< velocity of contact point on obj2 (expressed in global frame)
        const ChMaterialCompositeSMC& mat,  ///< composite material for contact pair
        double delta,                       ///< overlap in normal direction
        double eff_radius,                  ///< effective radius of curvature at contact
        double mass1,                       ///< mass of obj1
        double mass2,                        ///< mass of obj2
		ChContactable* objA,
		ChContactable* objB		
    )  const override {
		
		//
		// Get material properties
		//
		float mortar_layer=this->material->Get_mortar_h();
		float ENm=this->material->Get_ENm();
		float ENa=this->material->Get_ENa();			
		float alpha=this->material->Get_alpha();
		float beta=this->material->Get_beta();
		float np=this->material->Get_np();
		float sgmTmax=this->material->Get_sgmTmax();
		float sgmTau0=this->material->Get_sgmTau0();
		float kappa0=this->material->Get_kappa0();
		float eta_inf=this->material->Get_eta_inf();
		//printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f10 \t %f \t %f \n", mortar_layer, ENm, ENa, alpha, beta, 
		//	np, sgmTmax, sgmTau0, kappa0, eta_inf);
		
		//
		float h=mortar_layer;
		float R1=5.0;
		float R2=5.0; 
		//Material parameters of viscous part which should be defined inside contactinfo
		double eta0=kappa0*eta_inf;
		double Deps0=sgmTau0/eta0;
		
		
		//
    	//Get state variables form the map defined globally
    	//    	
		auto bodyA=dynamic_cast<ChBody*>(objA);
		auto bodyB=dynamic_cast<ChBody*>(objB);
		
		
		std::vector<unsigned int> ID={bodyA->GetId(),bodyB->GetId()};
		std::string mykey;
		
		if (ID[0]>ID[1]){			
		       mykey=std::to_string(ID[0])+"_"+std::to_string(ID[1]);		      
		}else{
		       mykey=std::to_string(ID[1])+"_"+std::to_string(ID[0]);		       
		}
		//
    	// Get state variables relevant to this contact	
		//
    	ChVector<> statevar=map_contact_info[mykey].strain;
    	
		//
		// initialize force values
		//
		ChVector<> force=(0,0,0);		
		//
		// If delta>h, objects are saparated
		//	
		if (delta<0)
			return force;
		//
		// Get local frame
		//
		//
		ChVector<> Vx, Vy, Vz;
		XdirToDxDyDz(normal_dir, VECT_Y, Vx, Vy, Vz);
        //std::cout<<" Vy: "<<Vy<<" Vz: "<<Vz<<"\t";	
				
		//
		// Get current_time
		//
		auto current_time=sys.GetChTime();
		//
		// Get the dimension of the object
		//
		if (bodyA->GetCollisionModel()->GetShape(0)->GetType()==0){        
			R1=bodyA->GetCollisionModel()->GetShapeDimensions(0)[0];
			
		}
		if (bodyB->GetCollisionModel()->GetShape(0)->GetType()==0){        
			R2=bodyB->GetCollisionModel()->GetShapeDimensions(0)[0];			
		}
		
		if (bodyA->GetCollisionModel()->GetShape(0)->GetType()!=0 || bodyB->GetCollisionModel()->GetShape(0)->GetType()!=0){ 
			h=mortar_layer/2;
		}
		//
		// Calculate contact plane rotation
		//		
		auto ref_rot=map_contact_info[mykey].mquaternion;
		ChMatrix33<> A0(ref_rot);	
	    ChMatrix33<> Aabs;
	    ChQuaternion<> abs_rot;
	    //
	    ChVector<> mXele = normal_dir;
	    ChVector<> myele = (bodyA->GetFrame_REF_to_abs().GetA().Get_A_Yaxis() + 
							bodyB->GetFrame_REF_to_abs().GetA().Get_A_Yaxis()).GetNormalized();
	    Aabs.Set_A_Xdir(mXele, myele);
	    abs_rot = Aabs.Get_A_quaternion();			
		//
		ChQuaternion<> q_delta=(abs_rot %  ref_rot.GetConjugate());
		//
		// update quaternion
		//
		//std::cout<<"current_time: "<<current_time<<" ref_rot : "<<ref_rot<< " abs_rot : "<<abs_rot<<"\t"<<normal_dir<<std::endl;
		//std::cout<<"q_delta  "<<q_delta<<std::endl;
		map_contact_info[mykey].mquaternion=abs_rot;
		//
		// Calculate contact area
		//		
		
		// center to center distance between 2 objects	
		double L0ij=R1+R2-h;
		double Lij=abs(R1+R2-delta);	
		double Rmin=std::min(R1, R2);	
		double Rmax=std::min(R1, R2);
		/*if(Lij<0.01*(Rmin-h))
			Lij=0.01*(Rmin-h);*/
		//double ai=(Rmax*Rmax-Rmin*Rmin+Lij*Lij)/Lij/2;
		double ai=abs(Rmin-delta/2);			
		double radius2=Rmin*Rmin-ai*ai;
		double contact_area=CH_C_PI*radius2;
		
		//
		// modify penetration according to initial overlap
		//
		double delta_new=delta-h;
		//
		//	
        // Relative velocities at contact at midplane (p0=(p1+p2)/2) 	
		//
		//			
		ChVector<> p0=(p1+p2)/2;		
		/*
		ChVector<> vcA=bodyA->GetFrame_COG_to_abs().GetPos_dt();
		ChVector<> vcB=bodyB->GetFrame_COG_to_abs().GetPos_dt();		
		//
		ChVector<> wcA=bodyA->GetFrame_COG_to_abs().GetWvel_loc();
		ChVector<> wcB=bodyB->GetFrame_COG_to_abs().GetWvel_loc();
		//
		ChVector<> velA=vcA+Vcross( bodyA->GetPos()-p0, wcA);
		ChVector<> velB=vcB+Vcross( bodyB->GetPos()-p0, wcB);
		*/				
		ChVector<> m_p1_loc = bodyA->Point_World2Body(p0);
		ChVector<> velA=bodyA->PointSpeedLocalToParent(m_p1_loc);
		
		ChVector<> m_p2_loc = bodyB->Point_World2Body(p0);
		ChVector<> velB=bodyB->PointSpeedLocalToParent(m_p2_loc);
		/*
		printf("  %10.6f  %10.6f  %10.6f\n ", vcA.x(), vcA.y(), vcA.z());
		printf("  %10.6f  %10.6f  %10.6f\n ", vcB.x(), vcB.y(), vcB.z());
		printf("  %10.6f  %10.6f  %10.6f\n ", wcA.x(), wcA.y(), wcA.z());
		printf("  %10.6f  %10.6f  %10.6f\n ", wcB.x(), wcB.y(), wcB.z());		
				
		printf("VA:  %10.6f  %10.6f  %10.6f\n ", velA.x(), velA.y(), velA.z());
		printf("VB:  %10.6f  %10.6f  %10.6f\n ", velB.x(), velB.y(), velB.z());	
		
		printf("VA2:  %10.6f  %10.6f  %10.6f\n ", velA2.x(), velA2.y(), velA2.z());
		printf("VB2:  %10.6f  %10.6f  %10.6f\n ", velB2.x(), velB2.y(), velB2.z());	
		
		printf("V1:  %10.6f  %10.6f  %10.6f\n ", vel1.x(), vel1.y(), vel1.z());
		printf("V2:  %10.6f  %10.6f  %10.6f\n ", vel2.x(), vel2.y(), vel2.z());
		*/
		//
        ChVector<> relvel = velB - velA;
        double relvel_n_mag = relvel.Dot(normal_dir);
        ChVector<> relvel_n = relvel_n_mag * normal_dir;
        ChVector<> relvel_t = relvel - relvel_n;
        double relvel_t_mag = relvel_t.Length();		
		//
		// Calculate displacement increment in normal and tangential direction
		//
		double dT = sys.GetStep();
		double delta_t = relvel_t_mag * dT;
		double delta_n = relvel_n_mag * dT;			
        ChVector<> v_delta_t = relvel_t * dT;
		//
		//  Calculate the strain increment in each local direction
		//
		double depsN=delta_n/Lij;
		double delta_M = relvel.Dot(Vy)*dT; //v_delta_t.Dot(Vy);
		double delta_L = relvel.Dot(Vz)*dT; //v_delta_t.Dot(Vz);
		double depsM=delta_M/Lij;
		double depsL=delta_L/Lij;		
		// 
		//
		//
		double epsA=log(1-h/(L0ij));			
		double epsN=log(Lij/L0ij);		
		//double epsN=statevar[0]+depsN;			
		double epsM=statevar[1]+depsM;
		double epsL=statevar[2]+depsL;
		//
		double epsT=pow((epsM*epsM+epsL*epsL),0.5);
		double epsQ=pow(epsN*epsN+alpha*epsT*epsT,0.5);
		//
		//
		//
		statevar[0]=epsN; statevar[1]=epsM;	 statevar[2]=epsL;		
		map_contact_info[mykey].strain=statevar;
		map_contact_info[mykey].step_time=current_time;
		//
		//
		//
		if (epsN<0){
			//
			// Compressive contact;
			//
			double sgmN=0;
			double sgmM=0;
			double sgmL=0;	
			double sgmT=0;			
			
			
			double stot;
			if (epsN>=epsA){				
				//stot=epsQ*ENm;	
				sgmN=epsN*ENm;
				map_contact_info[mykey].strain[1]=0; map_contact_info[mykey].strain[2]=0;
			}else{
				sgmN=(epsA)*ENm+(epsN-epsA)*ENa;
				sgmM=alpha*ENa*epsM;
				sgmM = sgn(epsM)*std::min<double>(abs(sgmM), mat.mu_eff * std::abs(sgmN))*0;
				sgmL=alpha*ENa*epsL;
				sgmL = sgn(epsL)* std::min<double>(abs(sgmL), mat.mu_eff * std::abs(sgmN))*0;
				//std::cout<<"delta_new: "<<delta_new <<" epsQ: "<<epsQ<<"  epsA: "<<epsA<<"  stot "<<stot<<"\n";
				//exit(0);
				/*
				stot=((epsA)*ENm+(epsQ-epsA)*ENa);
				
				if (epsQ!=0){
					sgmN=stot*epsN/epsQ;
					sgmT=alpha*stot*epsT/epsQ;
					sgmT = std::min<double>(sgmT, mat.mu_eff * std::abs(sgmN));
					if (epsT!=0){
						sgmM=sgmT*epsM/epsT;
						sgmL=sgmT*epsL/epsT;
					}
					//std::cout<<"sgmN: "<<sgmN<<"  sgmT: "<<sgmT<<"  sgmM: "<<sgmM<<"  sgmL: "<<sgmL<<std::endl;
				}
				*/
			}		
			
						
			//std::cout<<"delta_new: "<<delta_new<<" epsN: "<<epsN<< " epsQ: "<<epsQ<<" stot: "<<stot<<" sgmN: "<<sgmN<<std::endl;
			//////////////////////////////////////////////////////////
			//
			// Viscous stresses
			//		
			//////////////////////////////////////////////////////////
			//
			//depsN=relvel_n_mag/Lij*dT;
			//depsT=relvel_t_mag/Lij*dT;
			// calculate strain rates
			double vdepsN=depsN/dT;
			double vdepsM=depsM/dT;
			double vdepsL=depsL/dT;
			double vDeps=pow((beta*vdepsN*vdepsN+vdepsM*vdepsM+vdepsL*vdepsL),0.5);
			//
			double eta;
			if(vDeps<=Deps0) {
				eta=eta0;
			}else{
				eta=eta_inf*pow(abs(vDeps),np-1.);
			}			
			
			double sgmN_vis=beta*eta*vdepsN;
			double sgmM_vis=eta*vdepsM;
			double sgmL_vis=eta*vdepsL;
			double sgmT_vis=pow(sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis,0.5);
			//
			//	
			//
			double forceN=contact_area*(sgmN+sgmN_vis);	
			double forceM=contact_area*(sgmM+sgmM_vis);
			double forceL=contact_area*(sgmL+sgmL_vis);
			//std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
			//force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
			force=forceN*normal_dir+forceM*Vy+forceL*Vz;
			//std::cout<<" force "<<force<<std::endl;
			/*
			force = -forceN * normal_dir;		
			if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
				force -= (forceT / relvel_t_mag) * relvel_t;			 
            */			
			return -force;
			
		
		} else{
			//
			// Tensile contact
			//
						
			//////////////////////////////////////////////////////////
			//
			// Calculate Stress from material stiffness
			//
			//////////////////////////////////////////////////////////
			double sgmN=ENm*epsN;
			double sgmM=0;
			double sgmL=0;
			if (sgmN>sgmTmax)
				sgmN=sgmTmax;
			//
			//////////////////////////////////////////////////////////
			//
			// Viscous stresses
			//		
			//////////////////////////////////////////////////////////
			//			
			// calculate strain rates
			double vdepsN=depsN/dT;
			double vdepsM=depsM/dT;
			double vdepsL=depsL/dT;
			double vDeps=pow((beta*vdepsN*vdepsN+vdepsM*vdepsM+vdepsL*vdepsL),0.5);
			//
			double eta;
			if(vDeps<=Deps0) {
				eta=eta0;
			}else{
				eta=eta_inf*pow(abs(vDeps),np-1.);
			}			
			//std::cout<<"relvel"<<relvel<<" Deps0: "<<Deps0<<" vDeps"<<vDeps<<" eta: "<<eta<<std::endl;
			double sgmN_vis=beta*eta*vdepsN;
			double sgmM_vis=eta*vdepsM;
			double sgmL_vis=eta*vdepsL;
			//double sgmT=pow(sgmM*sgmM+sgmL*sgmL,0.5);
			//////////////////////////////////////////////////////////
			//
			// Combine Viscous stresses and stiffness stresses and calculate forces
			//		
			//////////////////////////////////////////////////////////				
			//exit(0);
			double forceN=contact_area*(sgmN+sgmN_vis);	
			double forceM=contact_area*(sgmM+sgmM_vis);
			double forceL=contact_area*(sgmL+sgmL_vis);
			//std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
			//double forceT = sgmT * contact_area;
			//force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
			force=forceN*normal_dir+forceM*Vy+forceL*Vz;
			/*			
			force = -forceN * normal_dir;
			if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
				force -= (forceT / relvel_t_mag) * relvel_t;			
			*/
			return -force;
		}
        
    }
	
	std::shared_ptr<ChMaterialFCM> Get_Material() const { return material; }  
    void Set_Material(std::shared_ptr<ChMaterialFCM> mmat) { material=mmat; }
    
	public:	
	
	std::shared_ptr<ChMaterialFCM> material;
	
};
