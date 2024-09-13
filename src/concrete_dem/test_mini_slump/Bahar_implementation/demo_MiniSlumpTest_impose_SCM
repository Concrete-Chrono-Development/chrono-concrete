
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
// Demo code about
// - UHPC mini slump test
// - using DEM 
//
// =============================================================================

#include "mycontact_V2_2_floc.cpp"

#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLinkMotionImposed.h"

#include "ChParticleEmitter_SMC.h"
#include "ChRandomShapeCreator_SMC.h"
//#include "chrono/particlefactory/ChParticleEmitter.h"
#include "chrono/particlefactory/ChParticleRemover.h"
#include "chrono/assets/ChTexture.h"

#include "chrono_thirdparty/rapidjson/prettywriter.h"
#include "chrono_thirdparty/rapidjson/stringbuffer.h"
#include "chrono_thirdparty/filesystem/path.h"

#include "chrono/motion_functions/ChFunctionPosition_XYZfunctions.h"
#include "chrono/motion_functions/ChFunctionRotation_axis.h"

#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

// Use the main namespace of Chrono, and other chrono namespaces
using namespace chrono;
using namespace chrono::particlefactory;
using namespace chrono::irrlicht;


double calc_aggVolFrac(std::vector<std::shared_ptr<ChBody>> &bodylist, double hlayer, double specimenVol){
	double avfrac=0;
	for (auto body:bodylist){
		if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0 )
			continue;
		double radius=body->GetCollisionModel()->GetShapeDimensions(0)[0]-hlayer;
		avfrac+=1.3333333333333333333*CH_C_PI*radius*radius*radius;
	}
	return avfrac/specimenVol;
}


/// Event trigger for particles inside a Cyl volume.
/// The event is triggered every time that the particle processor is run, not only the 1st time that particle enters the volume.
/// Only the center of gravity of particle is taken into consideration, regardless of its size.
class ChParticleEventTriggerCyl : public ChParticleEventTrigger {
  public:
    ChParticleEventTriggerCyl() { invert_volume = false; }

    /// This function triggers the a particle event according to the fact the particle is inside a Cyl.
    /// If SetTriggerOutside(true), viceversa triggers event outside the Cyl.
    virtual bool TriggerEvent(std::shared_ptr<ChBody> mbody, ChSystem& msystem) {
        ChVector<> particle_pos = mbody->GetPos();

        ChVector<> pos = m_frame.TransformPointParentToLocal(particle_pos);
		
		ChVector<> particle_vel = mbody->GetPos_dt();
		double vel_mag=particle_vel.Length();
		if (vel_mag>9000)
			mbody->SetPos_dt(particle_vel*9000/vel_mag);

        if ( !mbody->GetBodyFixed() && ( (fabs(pos.y()) < m_Cyl.h) && (  pow(pos.x()*pos.x()+pos.z()*pos.z(),0.5) < (m_Cyl.r/2+fabs(300.-pos.y())/300.*m_Cyl.r/2) )  ) ^ invert_volume)
            return true;
        else
            return false;
    }

    void SetTriggerOutside(bool minvert) { invert_volume = minvert; }

    geometry::ChCylinder m_Cyl;  ///< Cylinder volume
    ChFrame<> m_frame;      	///< Cylinder position and orientation

  protected:
    bool invert_volume;
};


class ChParticleRemoverCyl : public ChParticleProcessor {
  public:
    ChParticleRemoverCyl() {
        this->SetEventTrigger(std::shared_ptr<ChParticleEventTriggerCyl>(new ChParticleEventTriggerCyl));
        this->SetParticleEventProcessor(std::shared_ptr<ChParticleProcessEventRemove>(new ChParticleProcessEventRemove));
    }

    /// Set the dimensions and position of the trigger Cyl.
    void SetCyl(const double& radius, const double& height, const ChFrame<>& frame) {
        auto trigCyl = std::dynamic_pointer_cast<ChParticleEventTriggerCyl>(trigger);
        if (!trigCyl)
            throw ChException("ChParticleRemoverCyl had trigger replaced to non-Cyl type");
        trigCyl->m_Cyl.r=radius;
		trigCyl->m_Cyl.h=height;
        trigCyl->m_frame = frame;
    }

    geometry::ChCylinder& GetCyl() {
        auto trigCyl = std::dynamic_pointer_cast<ChParticleEventTriggerCyl>(trigger);
        if (!trigCyl)
            throw ChException("ChParticleRemoverCyl had trigger replaced to non-Cyl type");
        return trigCyl->m_Cyl;
    }

    /// Toggle inside/outside trigger.
    void SetRemoveOutside(bool invert) {
        auto trigCyl = std::dynamic_pointer_cast<ChParticleEventTriggerCyl>(trigger);
        if (!trigCyl)
            throw ChException("ChParticleRemoverCyl had trigger replaced to non-Cyl type");
        trigCyl->SetTriggerOutside(invert);
    }
};


class ChParticleEventTriggerBoxN : public ChParticleEventTrigger {
  public:
    ChParticleEventTriggerBoxN() { invert_volume = false; }

    /// This function triggers the a particle event according to the fact the the particle is inside a box.
    /// If SetTriggerOutside(true), viceversa triggers event outside the box.
    virtual bool TriggerEvent(std::shared_ptr<ChBody> mbody, ChSystem& msystem) {
        ChVector<> particle_pos = mbody->GetPos();
		ChVector<> particle_vel = mbody->GetPos_dt();
		double vel_mag=particle_vel.Length();
		if (vel_mag>9000)
			mbody->SetPos_dt(particle_vel*9000/vel_mag);

        ChVector<> pos = m_frame.TransformPointParentToLocal(particle_pos);

        if ( !mbody->GetBodyFixed()  && ( (fabs(pos.x()) < m_box.hlen.x()) && (fabs(pos.y()) < m_box.hlen.y()) && (fabs(pos.z()) < m_box.hlen.z())) ^
            invert_volume)
            return true;
        else
            return false;
    }

    void SetTriggerOutside(bool minvert) { invert_volume = minvert; }

    geometry::ChBox m_box;  ///< box volume
    ChFrame<> m_frame;      ///< box position and orientation

  protected:
    bool invert_volume;
};


class ChParticleRemoverBoxN : public ChParticleProcessor {
  public:
    ChParticleRemoverBoxN() {
        this->SetEventTrigger(std::shared_ptr<ChParticleEventTriggerBoxN>(new ChParticleEventTriggerBoxN));
        this->SetParticleEventProcessor(std::shared_ptr<ChParticleProcessEventRemove>(new ChParticleProcessEventRemove));
    }

    /// Set the dimensions and position of the trigger box.
    void SetBox(const ChVector<>& lengths, const ChFrame<>& frame) {
        auto trigbox = std::dynamic_pointer_cast<ChParticleEventTriggerBoxN>(trigger);
        if (!trigbox)
            throw ChException("ChParticleRemoverBox had trigger replaced to non-box type");
        trigbox->m_box.SetLengths(lengths);
        trigbox->m_frame = frame;
    }

    geometry::ChBox& GetBox() {
        auto trigbox = std::dynamic_pointer_cast<ChParticleEventTriggerBoxN>(trigger);
        if (!trigbox)
            throw ChException("ChParticleRemoverBox had trigger replaced to non-box type");
        return trigbox->m_box;
    }

    /// Toggle inside/outside trigger.
    void SetRemoveOutside(bool invert) {
        auto trigbox = std::dynamic_pointer_cast<ChParticleEventTriggerBoxN>(trigger);
        if (!trigbox)
            throw ChException("ChParticleRemoverBox had trigger replaced to non-box type");
        trigbox->SetTriggerOutside(invert);
    }
};


double calculateKE(ChSystem& sys){
	double KE=0;
	for (auto body:sys.Get_bodylist()){
		if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0 )
			continue;
		ChVector<> velo=body->GetFrame_COG_to_abs().GetPos_dt();
		double velMag=velo.Length();
		KE+=body->GetMass()*velMag*velMag/2.0;			
	}
	
	return KE;
};


std::shared_ptr<ChBody> AddConicalContainer(ChSystem& sys, std::shared_ptr<ChMaterialSurface> mat, double density, double radius, double height, std::string current_dir) {
    // create cylinder 
   
	auto cylinder = chrono_types::make_shared<ChBodyEasyMesh>(
        //(current_dir+"/coni_100_50_300.obj").c_str(),  // file name for OBJ Wavefront mesh   
        (current_dir+"/minislump-Cut4.obj").c_str(),  // file name for OBJ Wavefront mesh   **mini_slump2-Cut**
        density,                                         // density of the body
        true,                                         // automatic evaluation of mass, COG position, inertia tensor
        true,                                         // attach visualization asset
        true,                                         // enable the collision detection
        mat,                                     // surface contact material
        5  // radius of 'inflating' of mesh (for more robust collision detection)
    );	
	//
	//ChQuaternion<> q=Q_from_AngAxis(CH_C_PI_2, VECT_Z);	
	ChQuaternion<> q=Q_from_AngAxis(CH_C_PI_2, -VECT_X); //mini slump
	cylinder->SetRot(q);
	std::cout << "COG: " << cylinder->GetFrame_COG_to_abs().GetPos() << std::endl;
	cylinder->SetPos(ChVector<>(0,27.0167,0));
	//cylinder->GetVisualShape(0)->GetMaterial(0)->SetOpacity(0.5);
	//cylinder->SetBodyFixed(true);	
    sys.AddBody(cylinder); 
    return cylinder;	
	
}

std::shared_ptr<ChBody> AddConicalContainer2(ChSystem& sys, std::shared_ptr<ChMaterialSurface> mat, double density, double radius, double height, std::string current_dir) {
    // create cylinder 
   
	auto cylinder = chrono_types::make_shared<ChBodyEasyMesh>(
        (current_dir+"/minislump-Cut4.obj").c_str(),  // file name for OBJ Wavefront mesh
        density,                                         // density of the body
        true,                                         // automatic evaluation of mass, COG position, inertia tensor
        true,                                         // attach visualization asset
        true,                                         // enable the collision detection
        mat,                                     // surface contact material
        5  // radius of 'inflating' of mesh (for more robust collision detection)
    );	
	//	
	ChQuaternion<> q=Q_from_AngAxis(CH_C_PI_2*3, -VECT_X);
	cylinder->SetRot(q);	
	cylinder->SetPos(ChVector<>(0,(50.8-22.3)+50.8+3,0));
	cylinder->SetBodyFixed(true);	
    sys.AddBody(cylinder); 
    return cylinder;	
	
}


void ReadDFCparticles(ChSystem& sys, std::shared_ptr<ChVisualSystemIrrlicht>& vis, std::string& data_path, std::string& file_name, double rho, double h_layer) {
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
    //ChVector<> rotatedpos = rot.Rotate(pos);
		//ChVector<> rotatedpos=ChTransform<>::TransformParentToLocal(pos, t, q);
		ChVector<> rotatedpos=csys.TransformLocalToParent(pos);
		std::cout<<"pos: "<<pos<<"\t";
		std::cout<<"rotated pos: "<<rotatedpos<<"\n";
		
		double radius=particle.d/2+h_layer;
		double mass = rho/2*4/3*pow(radius, 3)*CH_C_PI;
					
		auto body = chrono_types::make_shared<ChBody>();
		body->SetInertiaXX((2.0 / 5.0) * mass * pow(radius, 2) * ChVector<>(1, 1, 1));
		body->SetMass(mass);
		//body->SetPos(ChVector<>(11 * ix , 50+11*iy, 11 * iz));
		body->SetPos(rotatedpos);
		body->GetCollisionModel()->ClearModel();
		body->GetCollisionModel()->AddSphere(mat, radius);
		body->GetCollisionModel()->BuildModel();
		body->SetCollide(true); 
		
		auto sphereMor = chrono_types::make_shared<ChSphereShape>(radius);
		//sphereMor->SetTexture(GetChronoDataFile("textures/blue.png"));
		sphereMor->SetColor(ChColor(128.f/255, 128.f/255, 128.f/255));
		sphereMor->SetOpacity(0.25f);
		body->AddVisualShape(sphereMor);
					
		auto sphereAgg = chrono_types::make_shared<ChSphereShape>(radius-h_layer);
		//sphereAgg->SetTexture(GetChronoDataFile("textures/pink.png"));
		sphereAgg->SetColor(ChColor(5.f/255, 48.f/255, 173.f/255));	
		body->AddVisualShape(sphereAgg);
		vis->BindItem(body);	
		sys.AddBody(body);
		
		body->SetLimitSpeed(true);
		body->SetMaxSpeed(2000.);
		body->SetMaxWvel(2000.);
		//body->SetBodyFixed(true);
		//sys.AddBatch(body);
    
       
    }
}


// =================================================================================================

int main(int argc, char* argv[]) {
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";
	SetChronoDataPath(CHRONO_DATA_DIR);
     ///////////////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// Get the current diectory
    ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
	//
    std::string current_dir(argv[0]);
    int pos = current_dir.find_last_of("/\\");
    current_dir=current_dir.substr(0, pos-5); 
	//
    // Create a Chrono physical system
    ChSystemSMC sys;
	  sys.Set_G_acc(ChVector<>(0,-9810,0));
	  sys.SetNumThreads(16,16,1);
    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("Particle emitter processor Motor MINI");
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox("/chrono-concrete/data/skybox/");
    vis->AddLight(ChVector<>(200,200,-250), 20, ChColor(255,255,255));
    //vis->AddTypicalLights();
    vis->AddCamera(ChVector<>(200, 200, -300));

    // Create the floor body
    auto floor_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    floor_mat->SetYoungModulus(1e2); //2e3
    floor_mat->SetFriction(0.5);
    floor_mat->SetRestitution(0.0f);
    floor_mat->SetAdhesion(0.0);
  	//floor_mat->SetCohesion(0);
    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(400, 10, 400, 2.4e-9, true, true, floor_mat);
    floorBody->SetPos(ChVector<>(0, -5, 0));
    floorBody->SetBodyFixed(true);
    floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    floorBody->GetVisualShape(0)->SetOpacity(0.2f);
    sys.Add(floorBody);
	
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat->SetYoungModulus(1e2); //2e3
    mat->SetFriction(0.5);
    mat->SetRestitution(0.0f);
    mat->SetAdhesion(0.0);
    //mat->SetCohesion(0);	
	  double cyl_radius=100;
   	double cyl_height=150;
	  double density=7.8E-9;
	//
	auto cylinderBody =AddConicalContainer(sys, mat, density*1000, cyl_radius, cyl_height, current_dir);
	cylinderBody->GetVisualShape(0)->SetOpacity(0.4f);
	cylinderBody->GetCollisionModel()->SetFamilyGroup(2);
  cylinderBody->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
  //cylinderBody->SetBodyFixed(true);
	//
	auto cylinderBody2 =AddConicalContainer2(sys, mat, density, cyl_radius, cyl_height, current_dir);
	cylinderBody2->GetVisualShape(0)->SetOpacity(0.4f);
	cylinderBody2->SetBodyFixed(true);
	cylinderBody2->SetCollide(true);	
	//
	//
	// Create the obstacle
    auto obstacle_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    obstacle_mat->SetYoungModulus(1e2); //2e3
    obstacle_mat->SetFriction(0.5);
    obstacle_mat->SetRestitution(0.0f);
    obstacle_mat->SetAdhesion(0.0);
	//obstacle_mat->SetCohesion(0);
    auto obstacleBody = chrono_types::make_shared<ChBodyEasyBox>(120., 3., 120., 2.4e-9, true, true, obstacle_mat);
    obstacleBody->SetPos(ChVector<>(0, 51.0, 0));
    obstacleBody->SetBodyFixed(true);
    obstacleBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/concrete.jpg"));
    obstacleBody->GetVisualShape(0)->SetOpacity(1.0f);
    //obstacleBody->GetVisualShape(0)->SetColor(ChColor(1.0f, 1.0f, (float)ChRandom()));
    sys.Add(obstacleBody);
    obstacleBody->GetCollisionModel()->SetFamilyGroup(2);
    obstacleBody->GetCollisionModel()->SetFamilyMaskNoCollisionWithFamily(1);
	//
	
	//volume of specimen
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
	//
	//
	// Define material parameters of fresh concrete
	//
	float mortar_layer=h_layer;  /// 1st param --> 0.75 / 1 / 1.10 / 1.20 / 1.40 
	float ENm=1.5E-2;     /// 2nd param --> 0.25 / 0.50 / 1 / 2 / 4
	float ENa=100;
	float h=mortar_layer;		
	float alpha=0.25;
	float beta=0.5;
	float np=1.0;
	float sgmTmax=1.0E-3;  ///3rd param 0.25 / 0.50 / 1 / 2 / 4 
	float sgmTau0=2.0E-4;  ///4th param 0.25 / 0.50 / 1 / 2 / 4  
	float kappa0=100.;
	float eta_inf=2.0E-5;/// 5th param 0.25 / 0.50 / 1 / 2 / 4 
	//
	//
  //
 //	Change default contact force algorithm to fresh concrete ones
	//
  auto Contact_force_algorithm = chrono_types::make_unique<MyContactForce>();
	auto materialFCM=chrono_types::make_shared<ChMaterialFCM>(ENm, ENa, mortar_layer, 
									alpha, beta, np, sgmTmax, sgmTau0, kappa0, eta_inf);
	Contact_force_algorithm->Set_Material(materialFCM);
  sys.SetContactForceAlgorithm(std::move(Contact_force_algorithm));
  //
    auto container = chrono_types::make_shared<MyContactContainer>();
    sys.SetContactContainer(container);	
  //
  //
	// Read and load predefined particle list
	//
	std::string data_path=current_dir;
	std::string file_name="DFCgeo000-data-particles-mini-2x4-dist"; 

	ReadDFCparticles(sys, vis, data_path, file_name, rho, mortar_layer);
	// 
  // Create an emitter:
    ChParticleEmitter_SMC emitter;

    // Set the flow rate [particles/s]:
    emitter.ParticlesPerSecond() = 30000;


    // Optional: limit the total n. of particles that can be generated
    emitter.SetUseParticleReservoir(true);
    emitter.ParticleReservoirAmount() = Va0/(4/3*CH_C_PI*pow(minD/2,3))*2; //was *2


	
	emitter.SetUseAgregateVolReservoir(true);
	emitter.SetMortarLayerThickness(mortar_layer);
  emitter.AggregateVolReservoirAmount() = targetVol;


    // ---Initialize the randomizer for positions
    auto emitter_positions = chrono_types::make_shared<ChRandomParticlePositionRectangleOutlet>();
    emitter_positions->Outlet() =
        ChCoordsys<>(ChVector<>(0, 105, 0), Q_from_AngAxis(CH_C_PI_2, VECT_X));  // center and alignment of the outlet
    emitter_positions->OutletWidth() = 90;
    emitter_positions->OutletHeight() = 90;

    emitter.SetParticlePositioner(emitter_positions);

    // ---Initialize the randomizer for alignments
    auto emitter_rotations = chrono_types::make_shared<ChRandomParticleAlignmentUniform>();

    emitter.SetParticleAligner(emitter_rotations);

    // ---Initialize the randomizer for velocities, with statistical distribution
    auto mvelo = chrono_types::make_shared<ChRandomParticleVelocityConstantDirection>();
    mvelo->SetDirection(-VECT_Y);
    mvelo->SetModulusDistribution(4000.0);

    emitter.SetParticleVelocity(mvelo);

	
	auto mcreator_spheres = chrono_types::make_shared<ChRandomShapeCreatorSpheres>();
	mcreator_spheres->SetDiameterDistribution(chrono_types::make_shared<ChConcreteDistribution>(minD, maxD, cement, AtoC, WtoC, mortar_layer));
	mcreator_spheres->SetDensityDistribution(chrono_types::make_shared<ChConstantDistribution>(rho));
    // Optional: define a callback to be exectuted at each creation of a box particle:
    class MyCreator_spheres : public ChRandomShapeCreator::AddBodyCallback {
        // Here do custom stuff on the just-created particle:
      public:
        virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
                               ChCoordsys<> mcoords,
                               ChRandomShapeCreator& mcreator) override {
            mbody->GetVisualShape(0)->SetColor(ChColor(1.0f, 1.0f, (float)ChRandom()));
        }
    };
    auto callback_spheres = chrono_types::make_shared<MyCreator_spheres>();
    mcreator_spheres->RegisterAddBodyCallback(callback_spheres);
	//
	// Create material for sphere body
	//
    auto sphere_mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
	sphere_mat->SetYoungModulus(1e2); //2e3
    sphere_mat->SetFriction(0.5);
    sphere_mat->SetRestitution(0.0f);
    sphere_mat->SetAdhesion(0.2);
	//
	mcreator_spheres->SetMaterial(sphere_mat);
    // Finally, tell to the emitter that it must use the 'mixer' above:
    emitter.SetParticleCreator(mcreator_spheres);
	
    // --- Optional: what to do by default on ALL newly created particles?
    //     A callback executed at each particle creation can be attached to the emitter.
    //     For example, we need that new particles will be bound to Irrlicht visualization:

    // a- define a class that implement your custom OnAddBody method...
    class MyCreatorForAll : public ChRandomShapeCreator::AddBodyCallback {
      public:
        virtual void OnAddBody(std::shared_ptr<ChBody> mbody,
                               ChCoordsys<> mcoords,
                               ChRandomShapeCreator& mcreator) override {
            // Enable Irrlicht visualization for all particles
            vis->BindItem(mbody);
        //     
	mbody->SetLimitSpeed(true);
	mbody->SetMaxSpeed(4000.);
	mbody->SetMaxWvel(4000.);
	//
            // Other stuff, ex. disable gyroscopic forces for increased integrator stabilty
            mbody->SetNoGyroTorque(true);
        }
        ChVisualSystemIrrlicht* vis;
    };
    // b- create the callback object...
    auto mcreation_callback = chrono_types::make_shared<MyCreatorForAll>();
    // c- set callback own data that he might need...
    mcreation_callback->vis = vis.get();
    // d- attach the callback to the emitter!
    emitter.RegisterAddBodyCallback(mcreation_callback);
    
	//
	
	double c_radius=51;
	double c_height=120;
	ChParticleRemoverCyl remover;
	remover.SetRemoveOutside(true);
  remover.SetCyl(c_radius, c_height, ChFrame<>(ChVector<>(0,0,0), QUNIT)); 
	
    // Bind all existing visual shapes to the visualization system
    //vis->AttachSystem(&sys);

    // Modify some setting of the physical system for the simulation, if you want
    //sys.SetSolverType(ChSolver::Type::PSOR);
    //sys.SetSolverMaxIterations(40);

    // Simulation loop
	int step_num=0;
	char filename[200];
	const std::string out_dir = current_dir + "/OUT_MiniMotor_h2_E1p5e-2_S1e-3_T2e-4_V2e-5/";
	// Create oputput directories
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        std::cerr << "Error creating directory " << out_dir << std::endl;
        return 1;
    }     
    double timestep = 5.0e-6;
    double current_time;
    
    auto body_list = sys.Get_bodylist();
    	
	
	// IMPOSE REMOVE CONE
	//cylinderBody->SetBodyFixed(false);
	auto f_xyz = chrono_types::make_shared<ChFunctionPosition_XYZfunctions>();
	//f_xyz->SetFunctionY(chrono_types::make_shared<ChFunction_Ramp>(0, -7));
	auto my_motion_function0 = chrono_types::make_shared<ChFunction_Const>(0.0);
    	auto my_motion_function1 = chrono_types::make_shared<ChFunction_Ramp>(0, 200.0);
    	auto my_motion_function2 = chrono_types::make_shared<ChFunction_Const>(200.0);
    	auto f_sequence = chrono_types::make_shared<ChFunction_Sequence>();
    	f_sequence->InsertFunct(my_motion_function0, 1.0, 0, true);
    	f_sequence->InsertFunct(my_motion_function1, 0.5, 0, true);
    	f_sequence->InsertFunct(my_motion_function2, 5.5, 0, true);
    	f_xyz->SetFunctionY(f_sequence);
    	//
    	auto impose = chrono_types::make_shared<ChLinkMotionImposed>();
    	sys.Add(impose);
    	impose->Initialize(cylinderBody, floorBody, ChFrame<>(cylinderBody->GetPos()));
    	impose->SetPositionFunction(f_xyz);	
	//
	// 2. phase: remove particles above the prism height
	//
    std::cout<<"\n phase: remove particles above the cone height\n"<<std::endl;	
//////////////////////////////////////////////////////////
	
	
    //ChParticleRemoverBox remover;
    //remover.SetRemoveOutside(true);
    //remover.SetBox(ChVector<>(210, 1000, 210), ChFrame<>(ChVector<>(0,0,0), QUNIT));

	//original back	
	//materialFCM->Set_ENm(ENm);
 	//materialFCM->Set_sgmTmax(sgmTmax);
	//materialFCM->Set_eta_inf(eta_inf);
	//materialFCM->Set_sgmTau0(sgmTau0); 
		
	c_radius=50.8-0.08;
	c_height=50.0*2;	
	remover.SetRemoveOutside(true);
        remover.SetCyl(c_radius, c_height, ChFrame<>(ChVector<>(0,0,0), QUNIT));

	
	auto t1=sys.GetChTime();
	int kk=0;

	timestep = 1.0e-8;
    while ( sys.GetChTime()<t1+0.50) {
		//vis->Run() & 
		//vis->BeginScene(true, true, ChColor(185./255,202./255,222./255));
		//vis->Render();  
		//vis->EndScene();
		//		
        // Continuosly create particle flow:
        //emitter.EmitParticles(sys, timestep);
		remover.ProcessParticles(sys);
		
		if (step_num%20==0){
		timestep=std::min(timestep*=1.05,5.0e-6);
		std::cout<<"timestep: "<<timestep<<std::endl;
		}
		
		
		// Continuosly check if some particle must be removed:

		//
        sys.DoStepDynamics(timestep);
					
		 if (std::fmod(step_num, 800) ==0) {
			double max_y=std::max(100-kk*5,100);
			//remover.SetBox(ChVector<>(210, max_y, 210), ChFrame<>(ChVector<>(0,0,0), QUNIT));
			remover.SetCyl(c_radius, max_y, ChFrame<>(ChVector<>(0,0,0), QUNIT));
		    //
		    sprintf(filename, "%s/slump_va_A%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, true);
		    sprintf(filename, "%s/slump_va_M%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, false);		    
		    std::cout<< "filename: "<<filename<<"\t";
			//
			double rho_new=materialFCM->CalculateEffDensity(sys, specimenVol, rho);
			if (rho_new<rho/4.)
				rho_new=rho/4.;
			std::cout<<"max_y: "<<max_y<<" rho_new: "<<rho_new<<"\t";
			materialFCM->ModifyDensity(sys, rho_new);
			//
			std::vector<std::shared_ptr<ChBody>> bodylist=sys.Get_bodylist();
			double aggVolFrac=calc_aggVolFrac(bodylist, mortar_layer, specimenVol);
			double current_time= sys.GetChTime();
			double IE;
			container->IterOnContactList(current_time,IE);
			double KE=calculateKE(sys);
			std::cout<<"time: "<< current_time << " va0: " << va0 <<  " aggVolFrac: "<< aggVolFrac << " IE: "<< IE << " KE: "<< KE << "\n";
			
			/*
			if(KE<1.0){
                        remover.ProcessParticles(sys);
                        }
                        */
			kk++;
		 }

		 step_num++;
		 
		 }
		 
	for (auto body:sys.Get_bodylist()){
		if (body->GetBodyFixed() || body->GetCollisionModel()->GetShape(0)->GetType()!=0 )
			continue;
		
		body->SetLimitSpeed(false);
		body->SetMaxSpeed(10000.);
		body->SetMaxWvel(10000.);				
	}
		 
	// WAIT	 
	std::cout<<"\n  phase: waiting  \n"<<std::endl;
	

	auto t2=sys.GetChTime();
	timestep = 5.0e-6;
	
	c_height=50.8-0.08;
	remover.SetCyl(c_radius, c_height, ChFrame<>(ChVector<>(0,0,0), QUNIT));
	//
    while ( sys.GetChTime()<t2+0.5) {
		//vis->Run() &
		//vis->BeginScene(true, true, ChColor(185./255,202./255,222./255));
		//vis->Render();  
		//vis->EndScene(); 
		//
        // Continuosly create particle flow:
        //emitter.EmitParticles(sys, timestep);
		
		// Continuosly check if some particle must be removed:
        //remover.ProcessParticles(sys);

        sys.DoStepDynamics(timestep);
		
		
		 if (std::fmod(step_num, 800) ==0) {
		    remover.ProcessParticles(sys);
		    sprintf(filename, "%s/slump_va_A%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, true);
		    sprintf(filename, "%s/slump_va_M%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, false);		    
		    std::cout<< "filename: "<<filename<<"\t";
			//
			std::cout<<"max_y: "<<"max_y_value "<<" rho_new: "<<" rho_value "<<"\t";
			//
			double current_time= sys.GetChTime();
			std::vector<std::shared_ptr<ChBody>> bodylist=sys.Get_bodylist();
			double aggVolFrac=calc_aggVolFrac(bodylist, mortar_layer, specimenVol);
			double IE;
			container->IterOnContactList(current_time,IE);
			double KE=calculateKE(sys);
			std::cout<<"time: "<< current_time << " va0: " << va0 <<  " aggVolFrac: "<< aggVolFrac << " IE: "<< IE << " KE: "<< KE << "\n";
	     
	     		 if(KE<1.0e-4){
                           break;
                        }
	     
	     
		 }
		 
                      
		 step_num++;
    }
		 
	//
	// 4. phase: remove containers and allow flow
	//
	std::cout<<"\ phase: flow  \n"<<std::endl;
	//
	//cylinderBody->SetCollide(true);
	cylinderBody2->SetCollide(false);
	//
	// 
	timestep=1e-8;
	//
	auto t3=sys.GetChTime();
    	while ( sys.GetChTime()<t3+6.0) {
		//vis->Run() & 
		//vis->BeginScene(true, true, ChColor(185./255,202./255,222./255));
		//vis->Render();  
		//vis->EndScene(); 
		//
		sys.DoStepDynamics(timestep);
       		if (step_num%10)
			timestep=std::min(timestep*=1.01, 5.0e-6);
		// 
                        
		if (std::fmod(step_num, 800) ==0) {
		    sprintf(filename, "%s/slump_va_A%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, true);
		    sprintf(filename, "%s/slump_va_M%06d.vtk", out_dir.c_str(), step_num);
		    WriteParticlesVTK(sys, filename, mortar_layer, false);		    
		    std::cout<< "filename: "<<filename<<"\t";

			//
			std::cout<<"max_y: "<<"max_y_value "<<" rho_new: "<<" rho_value "<<"\t";		 
		     //
			double current_time= sys.GetChTime();
			std::vector<std::shared_ptr<ChBody>> bodylist=sys.Get_bodylist();
			double aggVolFrac=calc_aggVolFrac(bodylist, mortar_layer, specimenVol);
			double IE;
			container->IterOnContactList(current_time,IE);
			double KE=calculateKE(sys);
			std::cout<<"time: "<< current_time << " va0: " << va0 <<  " aggVolFrac: "<< aggVolFrac << " IE: "<< IE << " KE: "<< KE << "\n";
			
			
			if(sys.GetChTime()>1.50 & obstacleBody->GetCollide()) {
                		cylinderBody->SetCollide(false);
                		obstacleBody->SetCollide(false);
                        }
		     
		 }
		 
		 step_num++;
    } 

    return 0;
}
