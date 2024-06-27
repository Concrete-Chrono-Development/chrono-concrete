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




class MyStateVar {
    public:
    
    MyStateVar() {};
    MyStateVar(double mradA, double mradB, ChVector<> mstrain, double mlambda=0)
    : radA(mradA), radB(mradB), strain(mstrain), mquaternion(QUNIT), lambda(mlambda)
    {};
    virtual ~MyStateVar() {}
    
    public:
	double step_time=0;
    float radA=1;
    float radB=1;
	float lambda=0;
    ChVector<> strain={0,0,0};
	ChQuaternion<> mquaternion;
    
};

std::map<std::string, MyStateVar> map_contact_info;




template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}


template<typename T>
T rungeKutta4(std::function<T(T, T, T, T, T, T)> f, T y0, T t0, T beta, T gammadot, T m, T Tcr, T h, int steps) {    
    T t = t0;
    T y = y0;
	h=h/steps;	
    for (int i = 0; i < steps; ++i) {
        T f1 = h * f(t, y, beta, gammadot, m, Tcr);
        T f2 = h * f(t + h / 2, y + f1 / 2, beta, gammadot, m, Tcr);
        T f3 = h * f(t + h / 2, y + f2 / 2, beta, gammadot, m, Tcr);
        T f4 = h * f(t + h, y + f3, beta, gammadot, m, Tcr);

        y = y + (f1 + 2 * f2 + 2 * f3 + f4) / 6;
        t = t + h;
    }   
    return y;
}


double fode(double t, double y) {
    return t * t - y; // Example ODE: y' = t^2 - y
}

double dfloc(double t, double y, double beta, double gammadot,  double m=0, double T=300.){ 
    return (1./T/pow(y,m) -beta*gammadot*y);
}



/// Class for fuller curve distribution of concrete aggregates
class ChApi ChConcreteDistribution : public ChDistribution {
  public:
    ChConcreteDistribution(double mminD, double mmaxD, double mcement, double mWtoC, 
						   double mAtoC, double mmortar_layer) : minD(mminD), maxD(mmaxD), cement(mcement), WtoC(mWtoC), AtoC(mAtoC), mortar_layer(mmortar_layer) {}

    /// Compute a random value whose probability is defined by the distribution,
    /// that is a value between min and max.
    virtual double GetRandom() override {
		//double va=1+c/rho_c+(WC*c)/rho_w+vair;
		//double va0=(1-(d0/da)**nF)*va;
		double q=3.-nF;
		double di=minD*pow( (1-::chrono::ChRandom() * ( 1.- pow((minD/maxD),q)) ), (-1./q) );	
		//double di=minD*pow( ( 1-::chrono::ChRandom() * ( 1- pow((0.5),q)) ), (-1./q) );
		return (di+2.*mortar_layer);
	};

  private:
    double minD;
    double maxD;
	double cement;
	double WtoC;
	double AtoC;
	double rho_c=3150;
	double rho_w=1000;
	double vair=0.03;
	double nF=0.5;	
	double mortar_layer=1.0;
};




class ChMaterialFCM {

  public:
    /// Construct a material.
    ChMaterialFCM(float mENm,
		float mENa,
		float mmortar_h,		
		float malpha,
		float mbeta,
		float mnp,
		float msgmTmax,
		float msgmTau0,
		float mkappa0,
		float meta_inf)
		: ENm(mENm),  ENa(mENa), mortar_h(mmortar_h),	alpha(malpha), beta(mbeta), np(mnp), 
		sgmTmax(msgmTmax), sgmTau0(msgmTau0), kappa0(mkappa0), eta_inf(meta_inf)
		 { };
		 

    ChMaterialFCM(){};

    // Destructor declared:
    ~ChMaterialFCM(){};



    /// Getter and Setter for material parameters.
    float Get_ENm() const { return ENm; }  
    void Set_ENm(float mENm) { ENm=mENm; }
    //
    float Get_ENa() const { return ENa; }  
    void Set_ENa(float mENa) { ENa=mENa; }
	//
	float Get_mortar_h() const { return mortar_h; }  
    void Set_mortar_h(float mmortar_h) { mortar_h=mmortar_h; }
	//
	float Get_alpha() const { return alpha; }  
    void Set_alpha(float malpha) { alpha=malpha; }
	//
	float Get_beta() const { return beta; }  
    void Set_beta(float mbeta) { beta=mbeta; }
	//
	float Get_np() const { return np; }  
    void Set_np(float mnp) { np=mnp; }
	//
	float Get_sgmTmax() const { return sgmTmax; }  
    void Set_sgmTmax(float msgmTmax) { sgmTmax=msgmTmax; }
	//
	float Get_sgmTau0() const { return sgmTau0; }  
    void Set_sgmTau0(float msgmTau0) { sgmTau0=msgmTau0; }
	//
	float Get_kappa0() const { return kappa0; }  
    void Set_kappa0(float mkappa0) { kappa0=mkappa0; }
	//
	float Get_eta_inf() const { return eta_inf; }  
    void Set_eta_inf(float meta_inf) { eta_inf=meta_inf; }
	
	float Get_flocbeta() const { return flocbeta;}  
    void Set_flocbeta(float mflocbeta) { flocbeta=mflocbeta;}
	
	float Get_flocm() const { return flocm; }  
    void Set_flocm(float mflocm) { flocm=mflocm; }
	
	float Get_flocTcr() const { return flocTcr; }  
    void Set_flocTcr(float mflocTcr) { flocTcr=mflocTcr; }
	
	float Get_FlocOnFlag() const { return FlocOnFlag; }  
    void Set_FlocOnFlag(float mFlocOnFlag) { FlocOnFlag=mFlocOnFlag; }
	
	double CalculateEffDensity(ChSystem& sys, double specVol, double rho_0){		
		double totVol=0;
		for (auto body:sys.Get_bodylist()){
			if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0 )
				continue;
			double radius=body->GetCollisionModel()->GetShapeDimensions(0)[0];
			totVol+=1.3333333333333333333*CH_C_PI*radius*radius*radius;
		}
		return rho_0*specVol/totVol;
	}
	
	void ModifyDensity(ChSystem& sys, double rho){
		double radius;
		double vol;
		double mass;
		for (auto body:sys.Get_bodylist()){
			if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0 )
				continue;
			radius=body->GetCollisionModel()->GetShapeDimensions(0)[0];
			mass=1.3333333333333333333*CH_C_PI*radius*radius*radius*rho;
			body->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector<>(1, 1, 1));
			body->SetMass(mass);
			body->SetDensity(rho);
		}		
	}
	
  public:
    float ENm=4.0E-2;
	float ENa=100;
	float mortar_h=3.0;		
	float alpha=0.25;
	float beta=0.5;
	float np=1.0;
	float sgmTmax=9.E-3;
	float sgmTau0=5.E-4;
	float kappa0=100.;
	float eta_inf=10.E-6;
	float flocbeta=1.0e-2;	
	float flocm=1.0;
	float flocTcr=120;
	bool FlocOnFlag = false;
};



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
		float flocbeta=this->material->Get_flocbeta();	
		float flocm=this->material->Get_flocm();
		float flocTcr=this->material->Get_flocTcr();
		//printf("%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f10 \t %f \t %f \n", mortar_layer, ENm, ENa, alpha, beta, 
		//	np, sgmTmax, sgmTau0, kappa0, eta_inf);
		
		//
		float h=mortar_layer;
		float R1=5000.0;
		float R2=5000.0; 
		double Rmin, Rmax, delta_rmin;
		double Lij, L0ij;
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
    	double lambda0=map_contact_info[mykey].lambda;
		//
		// initialize force values
		//
		ChVector<> force=(0,0,0);		
		//
		// If delta>h, objects are saparated
		//	
		if (delta<0){
			map_contact_info.erase(mykey);
			return force;			
		}
		
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
			// contact of DFC sphere with outer wall surface
			h=mortar_layer/2;
			/*ENm=ENm*250;
			sgmTmax=0;
			sgmTau0=0;*/
			Rmin=std::min(R1, R2);	
			Rmax=std::max(R1, R2);
			L0ij=Rmin-h;
			Lij=abs(Rmin-delta);	
			delta_rmin=delta;
			
		}else{
			// Contact of two DFC spheres
			// center to center distance between 2 objects	
			L0ij=R1+R2-h;
			Lij=abs(R1+R2-delta);	
			Rmin=std::min(R1, R2);	
			Rmax=std::max(R1, R2);
			delta_rmin=delta/2*(2.0*Rmax-delta)/Lij;
		}
		//
		// Calculate contact plane rotation
		//		
		/*auto ref_rot=map_contact_info[mykey].mquaternion;
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
		map_contact_info[mykey].mquaternion=abs_rot;*/
		//
		// Calculate contact area
		//		
		// Output the solution
		 /*{
		double y0 = 1; // Initial condition: y(0) = 1
		double t0 = 0; // Initial value of independent variable
		double h = 0.1; // Step size
		int NsubIter = 100; // Number of steps
		double flocbeta=0.01;
		double flocTcr=300;
		double gammadot=3;
		double flocm=1;
		
		std::cout << "Solution using Runge-Kutta 4th order method:" << std::endl;
		for (int i = 0; i < 10; ++i) {
			double solution = rungeKutta4<double>(dfloc, y0, t0, flocbeta, gammadot, flocm, flocTcr, h, NsubIter);
			std::cout << "t = " << t0 + h << ", y = " << solution << std::endl;
			t0=t0+h;
			y0=solution;
		}
		
		exit(9);
		 };*/
		 
		double ai=abs(Rmin-delta_rmin);			
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
			if(this->material->FlocOnFlag){	
				double lambda = rungeKutta4<double>(dfloc, lambda0, current_time, flocbeta, vDeps, flocm, flocTcr, dT, 1);
				map_contact_info[mykey].lambda=lambda;
				lambda=0;
				if (vDeps>0){				
					eta=(1.- exp(-vDeps/Deps0) ) * ( eta_inf + sgmTau0 * (1.0+lambda)/vDeps );				
				}else{
					eta= eta0*(1.0+lambda); //sgmTau0 * (1.0+lambda)/Deps0;
				}
			}else{
				if(vDeps<=Deps0) {
					eta=eta0;
				}else{
					eta=eta_inf*pow(abs(vDeps),np-1.)+sgmTau0*(1.0+0)/vDeps;
				}	
				
			}
			
			
			double sgmN_vis=beta*eta*vdepsN;
			double sgmM_vis=eta*vdepsM;
			double sgmL_vis=eta*vdepsL;
			//double sgmT_vis=pow(sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis,0.5);
			//std::cout<<"t: "<<current_time<<"\tgammadot: "<<vDeps<<"\t eta: "<<eta<<"\t"<<"sgm_vis: "<<std::sqrt(sgmN_vis*sgmN_vis+sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis)<< " flocm "<< flocm << "\n";
			//
			//FILE *fptr_floc ;
			//auto outfilename_floc="/home/erol/BAHAR/FLOC/case1_floc_pos_ux/2Spheres_Case1/ReportFloc.txt";
			//fptr_floc= fopen(outfilename_floc, "w");
			//fprintf(fptr_floc,"time: %f \t", current_time);
    			//fprintf(fptr_floc,"gammadot: %f \t", vDeps);
    			//fprintf(fptr_floc,"eta: %f \t", eta);
    			//fprintf(fptr_floc,"sigma_vis: %f \n", std::sqrt(sgmN_vis*sgmN_vis+sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis));
    			//fclose(fptr_floc);
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
			if(this->material->FlocOnFlag){	
				double lambda = rungeKutta4<double>(dfloc, lambda0, current_time, flocbeta, vDeps, flocm, flocTcr, dT, 1);
				map_contact_info[mykey].lambda=lambda;
				lambda=0;
				if (vDeps>0){				
					eta=(1.- exp(-vDeps/Deps0) ) * ( eta_inf + sgmTau0 * (1.0+lambda)/vDeps );				
				}else{
					eta= eta0*(1.0+lambda); //sgmTau0 * (1.0+lambda)/Deps0;
				}
			}else{
				if(vDeps<=Deps0) {
					eta=eta0;
				}else{
					eta=eta_inf*pow(abs(vDeps),np-1.)+sgmTau0*(1.0+0)/vDeps;
				}	
				
			}
			
			//std::cout<<"relvel"<<relvel<<" Deps0: "<<Deps0<<" vDeps"<<vDeps<<" eta: "<<eta<<std::endl;
			double sgmN_vis=beta*eta*vdepsN;
			double sgmM_vis=eta*vdepsM;
			double sgmL_vis=eta*vdepsL;
			//std::cout<<"t: "<<current_time<<"\tgammadot: "<<vDeps<<"\t eta: "<<eta<<"\t"<<"sgm_vis: "<<std::sqrt(sgmN_vis*sgmN_vis+sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis)<< " flocm "<< flocm << "\n";
			//
			//FILE *fptr_floc ;
			//auto outfilename_floc="/home/erol/BAHAR/FLOC/case1_floc_pos_ux/2Spheres_Case1/ReportFloc.txt";
			//fptr_floc= fopen(outfilename_floc, "w");
			//fprintf(fptr_floc,"time: %f \t", current_time);
    			//fprintf(fptr_floc,"gammadot: %f \t", vDeps);
    			//fprintf(fptr_floc,"eta: %f \t", eta);
    			//fprintf(fptr_floc,"sigma_vis: %f \n", std::sqrt(sgmN_vis*sgmN_vis+sgmM_vis*sgmM_vis+sgmL_vis*sgmL_vis));
    			//fclose(fptr_floc);
    			//
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


class MyContactContainer : public ChContactContainerSMC {
  public:
    MyContactContainer() {}
    ~MyContactContainer() {}
    // Traverse the list contactlist_6_6
    
    
    void IterOnContactList(double& current_time, double& IE) {
		int num_contact = 0;
		IE=0;
		std::map<std::string, MyStateVar> updated_map_contact_info;
        auto iter = contactlist_6_6.begin();        
        while (iter != contactlist_6_6.end()) {
            ChVector<> p1 = (*iter)->GetContactP1();
            ChVector<> p2 = (*iter)->GetContactP2();
			ChVector<> force = (*iter)->GetContactForce();
			ChVector<> dp=p1-p2;
			IE=IE+std::abs(force.x()*dp.x()+force.y()*dp.y()+force.z()*dp.z());
            double CD = (*iter)->GetContactDistance();
            ChBody* bodyA = dynamic_cast<ChBody*>((*iter)->GetObjA());
			ChBody* bodyB = dynamic_cast<ChBody*>((*iter)->GetObjB());
	        ///
	        ///
	        ///
            //std::cout<<"objA: "<<bodyA->GetId();
			//std::cout<<"\t objB: "<< bodyB->GetId()<<std::endl;                   
			std::vector<unsigned int> ID={bodyA->GetId(),bodyB->GetId()};
			std::string mykey;
		    ///
			///
			if (ID[0]>ID[1]){
				
				   mykey=std::to_string(ID[0])+"_"+std::to_string(ID[1]);
				   //std::cout<<mykey<<std::endl;
			}else{
				   mykey=std::to_string(ID[1])+"_"+std::to_string(ID[0]);
				   //std::cout<<mykey<<std::endl;
			}
			
			updated_map_contact_info[mykey]=map_contact_info[mykey];			
			
            ///
			///
            num_contact++;
            ++iter;
        }
		
		map_contact_info=updated_map_contact_info;	
        
    }



	int TrackHardContact(double& h_layer) {
		int num_hard_contact = 0;
		double max_CD=0;
        auto iter = contactlist_6_6.begin();        
        while (iter != contactlist_6_6.end()) {           
            double CD = (*iter)->GetContactDistance(); 
			if (CD<max_CD)
				max_CD=CD;
            if (CD<=-2*h_layer) 
				num_hard_contact++;
            ++iter;
        }	
		//std::cout<<"  max_CD "<<max_CD<<"\n";
        return num_hard_contact;
		
    }	
 
    
};


// -----------------------------------------------------------------------------
// Callback class for contact reporting
// -----------------------------------------------------------------------------
class ContactReporter : public ChContactContainer::ReportContactCallback {
  public:
    ContactReporter(FILE* fptr, std::shared_ptr<ChBody> objA, std::shared_ptr<ChBody> objB) : m_objA(objA), m_objB(objB), m_fptr(fptr) {}

  private:
    virtual bool OnReportContact(const ChVector<>& pA,
                                 const ChVector<>& pB,
                                 const ChMatrix33<>& plane_coord,
                                 const double& distance,
                                 const double& eff_radius,
                                 const ChVector<>& cforce,
                                 const ChVector<>& ctorque,
                                 ChContactable* modA,
                                 ChContactable* modB) override {
		//
		//
		auto bodyA=dynamic_cast<ChBody*>(modA);
		auto bodyB=dynamic_cast<ChBody*>(modB);
		ChVector<> pcA=bodyA->GetFrame_COG_to_abs().GetPos();
		ChVector<> pcB=bodyB->GetFrame_COG_to_abs().GetPos();
		
		//fprintf(m_fptr,"center of sphere A: %10.6f  %10.6f  %10.6f", pcA.x(), pcA.y(), pcA.z());
		//fprintf(m_fptr,"  center of sphere B: %10.6f  %10.6f  %10.6f", pcB.x(), pcB.y(), pcB.z());
		fprintf(m_fptr," %10.6f  %10.6f  %10.6f ", pcA.x(), pcA.y(), pcA.z());
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", pcB.x(), pcB.y(), pcB.z());
		
		ChVector<> pcA0=bodyA->GetFrame_REF_to_abs().GetPos();
		ChVector<> pcB0=bodyB->GetFrame_REF_to_abs().GetPos();
		
		ChVector<> vcA=bodyA->GetFrame_COG_to_abs().GetPos_dt();
		ChVector<> vcB=bodyB->GetFrame_COG_to_abs().GetPos_dt();
		//fprintf(m_fptr,"  vel of sphere A: %10.6f  %10.6f  %10.6f", vcA.x(), vcA.y(), vcA.z());
		//fprintf(m_fptr,"  vel of sphere B: %10.6f  %10.6f  %10.6f", vcB.x(), vcB.y(), vcB.z());
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", vcA.x(), vcA.y(), vcA.z());
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", vcB.x(), vcB.y(), vcB.z());
		
		ChVector<> wcA=bodyA->GetFrame_COG_to_abs().GetWvel_loc();
		ChVector<> wcB=bodyB->GetFrame_COG_to_abs().GetWvel_loc();		
		//fprintf(m_fptr,"  angular vel of sphere A: %10.6f  %10.6f  %10.6f", wcA.x(), wcA.y(), wcA.z());
		//fprintf(m_fptr,"  angular vel of sphere B: %10.6f  %10.6f  %10.6f", wcB.x(), wcB.y(), wcB.z());
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", wcA.x(), wcA.y(), wcA.z());
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", wcB.x(), wcB.y(), wcB.z());
		//
		//
		
        // Check if contact involves objA
        if (modA == m_objA.get()) {
            //fprintf(m_fptr,"  A contact on sphere 1 at pos: %10.6f  %10.6f  %10.6f", pA.x(), pA.y(), pA.z());
			fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", pA.x(), pA.y(), pA.z());
        } else if (modB == m_objA.get()) {
            //fprintf(m_fptr,"  B contact on sphere 1 at pos: %10.6f  %10.6f  %10.6f", pB.x(), pB.y(), pB.z());
			fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", pB.x(), pB.y(), pB.z());
        }

        // Check if contact involves objB
        if (modA == m_objB.get()) {
            //fprintf(m_fptr,"  A contact on sphere 2 at pos: %10.6f  %10.6f  %10.6f", pA.x(), pA.y(), pA.z());
			fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", pA.x(), pA.y(), pA.z());
        } else if (modB == m_objB.get()) {
            //fprintf(m_fptr,"  B contact on sphere 2 at pos: %10.6f  %10.6f  %10.6f", pB.x(), pB.y(), pB.z());
			fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", pB.x(), pB.y(), pB.z());
        }
				
		
        const ChVector<>& nrm = plane_coord.Get_A_Xaxis();
		const ChVector<>& nrt1 = plane_coord.Get_A_Yaxis();
		const ChVector<>& nrt2 = plane_coord.Get_A_Zaxis();
        //fprintf(m_fptr,"  nrm: %10.6f, %10.6f  %10.6f", nrm.x(), nrm.y(), nrm.z());
        //fprintf(m_fptr,"  frc: %12.8e  %12.8e  %12.8e", cforce.x(), cforce.y(), cforce.z());		
        ////printf("  trq: %7.3f, %7.3f  %7.3f", ctorque.x(), ctorque.y(), ctorque.z());
        //fprintf(m_fptr,"  penetration: %12.8e   eff. radius: %7.3f\n", distance, eff_radius);
		
		fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", nrm.x(), nrm.y(), nrm.z());
		//fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", nrt1.x(), nrt1.y(), nrt1.z());
		//fprintf(m_fptr,"  %10.6f  %10.6f  %10.6f ", nrt2.x(), nrt2.y(), nrt2.z());
        fprintf(m_fptr,"  %12.8e  %12.8e  %12.8e ", cforce.x(), cforce.y(), cforce.z());        
        fprintf(m_fptr,"  %12.8e   %7.3f\n", distance, eff_radius);

        return true;
    }

    std::shared_ptr<ChBody> m_objA;
    std::shared_ptr<ChBody> m_objB;
	FILE *m_fptr;
};



 
 // Function to write particle positions and radii to a VTK file
void WriteParticlesVTK(ChSystem& sys, const std::string& filename, double h_layer, bool isAgregateOut) {
     // Get the number of particles
     auto body_list= sys.Get_bodylist();
     
	 std::vector<std::shared_ptr<ChBody>> body_list_new;
	 for (auto body:body_list){
		 if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0)
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
         ChVector<float> pos = body_list_new[i]->GetPos();
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
		if (isAgregateOut)
			vtk_file << (body_list_new[i]->GetCollisionModel()->GetShapeDimensions(0)[0]-h_layer) << "\n";
		else
			vtk_file << body_list_new[i]->GetCollisionModel()->GetShapeDimensions(0)[0] << "\n";
     }

     // Write the particle velocities
     vtk_file << "\nVECTORS velocity float\n";
     for (int i = 0; i < num_particles; i++) {
         ChVector<float> vel = body_list_new[i]->GetPos_dt();;
         vtk_file << vel.x() << " " << vel.y() << " " << vel.z() << "\n";
     }
	 
	 // Write the particle angular velocities
     vtk_file << "\nVECTORS angular_velocity float\n";
     for (int i = 0; i < num_particles; i++) {
         ChVector<float> w = body_list_new[i]->GetFrame_COG_to_abs().GetWvel_loc();
         vtk_file << w.x() << " " << w.y() << " " << w.z() << "\n";
     }

     // Close the file
     vtk_file.close();
}
