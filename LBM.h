#include<iostream>
#include<cmath>
using namespace std;



// -----------------------------------------DEFINITIONS OF SIMULATION PARAMETERS ----------------------------------
#ifndef DEFINES_H
	#define DEFINES_H
		
	//FLOW PARAMETERS
	#define rho0 1.0
	#define Re 100.0
	
	
	//SIMULATION DOMAIN 
	#define BCXM -2.5
	#define BCXP 7.5
	#define BCYM -2.5
	#define BCYP 2.5

	//SIMULATION PARAMETERS
	static const int Nx = 800;
	static const int Ny = 800;
	static const int Nz = 1;

	#define dt 1.0
	#define dx 1.0
	#define dy 1.0
	#define dz 1.0
	static const double cs2 = pow(dx,2)/(3*pow(dt,2));
	#define p0 rho0*cs2

	//static const double nu = 0.1*40/Re; //18.75/(75*5)*2.0/3.0;  //problem here, no need to change the computational method for lift and drag (probably)
	//static const double tau = 0.5*dt+nu/cs2;

	// or
	static const double tau = 0.545;//3.5;//0.95;
	static const double nu = cs2*(tau - dt/2);
	
	//FLOW CASES (Basic Bounce Back)
	//#define Couette_2D
	//#define TaylorGreen_2D
	//#define Cylinder_Basic
	//#define Poiseuille_2D
	
	//Special Case, needed for some occasion of outputing pressure
	//#define Special_Case

	//Flow Cases (IBM)
	#define IBM
	//#define Cylinder_IBM
	#define Cylinder_Translationally_Oscillating
	//#define Cylinder_Rotationally_Oscillating
	//#define Accelerated_Normal_Flat_Plate
	
	//For Couette flow
	#if defined Couette_2D 
		#define u_wall 0.1
	#else
		#define u_wall 0
	#endif

	//Parameters only for Taylor-Green vortex : 
	#if defined TaylorGreen_2D
		//Variables : 
		static const double lx = (Nx-1)*dx, ly = (Ny-1)*dy; 
		static const double PI = 3.14159265;
		static const double kx = 2*PI/(lx+dx), ky = 2*PI/(ly+dy);
		static const double td = 1/(nu*(kx*kx + ky*ky)); //decay time;
		#define u_max  0.03*dx/dt
	#endif

	//Parameters for Poiseuille flow : Pressure gradient 
	#if defined Poiseuille_2D
		static const double dp_dx = -0.0001/5; //Pressure gradient, negative for the flow to accelerate
	#endif
	
	//Parameters for Cylinder IBM : Stencils, stiffness, number of nodes (Explicit Direct Forcing)
	#if defined Cylinder_IBM
		//#define EXPLICIT_DIRECT_FORCING
		#define MULTI_DIRECT_FORCING
		
		//Properties of IBM in general
		static const int stencil = 4; //Stencils can be 2 3 and 4
		static const int number_nodes = 204; //Number of marker at the cylinder

		#if defined EXPLICIT_DIRECT_FORCING
			static const double stiffness = 2; //Stiffness of the marker for forcing
		#elif defined MULTI_DIRECT_FORCING
			static const int m_max = 0;
		#endif
		
	#endif



	//Moving boundary Flow
	//-----------------------------------------------------------------------------------------	
	//For the case of translationally oscillating cylinder 
	#if defined Cylinder_Translationally_Oscillating
		#define MULTI_DIRECT_FORCING

		//Properties of IBM in general
		static const int stencil = 4; //Stencils can be 2 3 and 4
		static const int number_nodes = 204; //Number of marker at the cylinder

		#if defined MULTI_DIRECT_FORCING
			static const int m_max = 5;
		#endif

	#endif

	//For the case of rotationally oscillating cylinder
	#if defined Cylinder_Rotationally_Oscillating
		#define MULTI_DIRECT_FORCING

		//Properties of IBM in general
		static const int stencil = 4; //Stencils can be 2 3 and 4
		static const int number_nodes = 204; //Number of marker at the cylinder

		#if defined MULTI_DIRECT_FORCING
			static const int m_max = 5;
		#endif
	#endif

	//OUTPUT VALUE
	#if defined TaylorGreen_2D 
		#define T_OUT td*5
		#define Nt 10
	#else
		#define T_OUT 50000
		#define Nt 50
	#endif
	

	//VELOCITY SET 
	#define D2Q9
	
	//BOUNDARY CONDITION
	#define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_P 3 // periodic boundary
	#define TYPE_L 4 // Slip boundary
	
#endif
// ----------------------------------------------------------------------------------------------------------------------






//-------------------------------DEFINITION FOR LBM SIMULATION----------------------------------------------------------
#ifndef LBM_H
#define LBM_H
	//Velocity set 
	#ifdef D2Q9 
		const double ndim = 2;
		const double npop = 9;
		//Velocity set 
		const int cx[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
		const int cy[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
		const int cz[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
		//weight function 
		const double w[9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
		//opposite for D2Q9 velocity set : 
		const int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //For No slip boundary condition
		////////////////////////{0, 1, 2, 3, 4, 5, 6, 7, 8};
		const int slip[9] = {0, 1, 4, 3, 2, 7, 8, 5, 6}; //For slip boundary condition
	#endif
	
	//lattice definition 
	class LATTICE{
		public : 
			short type = TYPE_F; //details in setup
			double density = 1.0;
			double ux = 0.0;
			double uy = 0.0;
			double uz = 0.0;
			double f[9], f_star[9];
			double F[3] = {0,0,0}; //External force, forceless if Couette flow, Taylor-Green and Cylinder Bounce Back
			bool inside = 0;
	};
	
	//Marker definition, will be used in the form of array
	class MARKER {
		public : 
			//Variables : 
			double x = 0, y = 0, z = 1; //Current Coordinate, default to be 2D
			double x_ref = 0, y_ref = 0, z_ref = 1; //Reference coordinate, default to be 2D
			double rdot[3] = {0,0,0}; //velocity of marker, default to be not moving
			double F[3] = {0,0,0}; //Lagrangian force, default to be forceless

			#if defined MULTI_DIRECT_FORCING
				double u_boundary[3] = {0,0,0}; //Desired Boundary velocity 
			#endif

			
	};

	//LBM class definition 
	class LBM {
		private : 
			int Nx, Ny, Nz = 1;
			double tau;
		public : 
			LATTICE ***fluid; //Definition of lattice for fluid (IBM), accessible to all
			

			LBM(int Nx, int Ny, int Nz, double tau); //Parametric constructor for Lattice
			
			//Mutator
			void Init();
			

			void Collision(); 
			void Streaming(); 
			void MacroProp();
			
			//Accessor
			int getNx() const;
			int getNy() const;
			int getNz() const;
			double getTAU() const;
	};
	
#endif
// ----------------------------------------------------------------------------------------------------------------------





//-----------------------------------STEPS IMMERSED BOUNDARY METHODS ---------------------------------------------
#ifndef IBM_H
#define IBM_H
	//steps for Direct Forcing - Immersed Boundary Method (DF-IB-LBM) (RIGID BOUNDARY): 
	//0. Compute the new position of the center 
	//1. Compute Lagrangian Force from current boundary configuration (Varies for rigid and deformable)
	//2. Spread Lagrangian force to lattice via Eulerian force density (remember to use the stencil of Dirac delta kernel function)
	//3. Compute uncorrected velocity
	//4. LBM (Collision, streaming) with forcing
	//5. compute macroscopic properties.
	//6. Interpolate fluid velocity at Lagrangian node to obtain the speed of marker
	//7. Update marker configuration (relative to COM and COM motion)
	//Step 6 and 7 doesn't need to be done for rigid boundary

	double kernel(double distance, double delta);

	//Step 0 : Compute new position
	#if defined Cylinder_Translationally_Oscillating 
		void Update_Position_COM(vector<MARKER>&, double&, double*, double*); //For COM motion with externally specified motion
	#elif defined Cylinder_Rotationally_Oscillating
		void Update_Position_COM(vector<MARKER>&, double&, double*, double*, double&, double&, double&);
	#endif
	

	//Step 1 : Compute Lagrangian Force
	void Compute_Lagrangian(vector<MARKER>&, LBM, LBM&);
	
	//Step 2 : Spread Force to Eulerian Force
	void Spread_Force(LBM&, vector<MARKER>&);
	
	//Step 3 : Compute uncorrected velocity 
	void Compute_Uncorrected_Velocity(LBM&);

	//Step 4 : LBM 					
	//Step 5 : Macroprop
	void MacroProp_IBM(LBM&);

	//Step 6 : Interpolate fluid velocity to node position
	void Interpolate(LBM&, vector<MARKER>&);

	//Step 7 : New boundary configuration (relative to center of mass)
	void Update_Position(vector<MARKER>&); //For relative motion
	
	

#endif
//---------------------------------------------------------------------------------------------------------------





//DONE FOR IBM
// --------------------------------------SETUP OF THE PROBLEM -----------------------------------------------------------
#ifndef SETUP_H
	#define SETUP_H
	
	#if defined Couette_2D
		LBM main_setup_Couette();
	#elif defined TaylorGreen_2D
		LBM main_setup_Taylor_Green_2D();
	#elif defined Cylinder_Basic
		static const double D = 50;
		static const double u_max = Re*nu/D;
		static double center_x = 2*D;//dx*(Nx/5.0);
		static double center_y = 2*D;//dy*(Ny/2.0);
		LBM main_setup_Cylinder();
	#elif defined Poiseuille_2D
		LBM main_setup_Poiseuille_2D();
	#elif defined Cylinder_IBM
		static const double D = 40;
		static const double u_max = Re*nu/D;
		static double center_x = dx*(Nx/5.0);
		static double center_y = dy*(Ny/2.0);
		LBM main_setup_Cylinder_IBM(); //Fluid
	#elif defined Cylinder_Translationally_Oscillating
		static const double D = 50;
		static const double u_max = Re*nu/D;
		static double center_x = Nx/2.0;
		static double center_y = Ny/2.0;

		//Parameter for this case
		static double KC = 5; //Keulegan Carpenter number 
		static double T = KC*D/u_max; //Period of oscillation

		//For calculation of internal mass effect
		static double rho_f = 1.0;
		static double rho_b = 2.8;
		static double M = rho_b*3.14159265*0.25*pow(D,2);

		//function for this case
		LBM main_setup_Cylinder_Translation_Oscillate(); //Fluid
	
	#elif defined Cylinder_Rotationally_Oscillating
		static const double D = 50.0;
		static const double PI = 3.14159265;

		//Parameters used in simulation 
		static const double Re_r = 1.0; //Rotation Reynolds number
		static const double dtheta = PI/4.0; //Amplitude for angle
		static const double u_max = nu*Re_r/D; //Maximum velocity
		static const double Omega_max = 2*u_max/D; //Maximum rotational speed
		static const double T = 2*PI*dtheta/Omega_max; //Period of oscillation

		//Position of cylinder
		static double center_x = Nx/2.0; 
		static double center_y = Ny/2.0;

		//For calculation of internal mass effect
		static double rho_f = 1.0;
		static double rho_b = 2.8;
		static double M = rho_b*3.14159265*0.25*pow(D,2); //Total mass
		static double I_B = 0.5*M*pow(D/2.0,2.0); //moment of inertia 

		//function for this case
		LBM main_setup_Cylinder_Rotation_Oscillate(); //Fluid
	#endif

#endif


//-----------------------------------------------------------------------------------------------------------------------










//DONE FOR GEOMETRY DEFINITION IN IBM
// --------------------------------------DEFINITION OF PROBLEM GEOMETRY -------------------------------------------------
/* Considered Problem : 
1. Couette Flow :
2. Taylor Green vortex : 
3. Flow over cylinder : 
4. Poiseuille flow : 
5. IBM
*/



#ifndef GEOM_H
	#define GEOM_H
	
	#if defined Cylinder_Basic 
		void Generate_Cylinder(LBM& lb, int& count);
	#elif defined Cylinder_IBM 
		void main_setup_marker(LBM& lb, vector<MARKER>& marker); //Marker for cylinder and update whether it is inside or outside the boundary
	#elif defined Cylinder_Translationally_Oscillating
		void main_setup_marker_Translation(LBM& lb, vector<MARKER>& marker); //Marker only for cylinder
	#elif defined Cylinder_Rotationally_Oscillating
	 	void main_setup_marker_Rotation(LBM& lb, vector<MARKER>& marker); //Marker only for cylinder
	#endif
	
#endif

// ----------------------------------------------------------------------------------------------------------------------









// -------------------------------------------OUTPUT HEADER--------------------------------------------------------------
#ifndef OUTPUT_H
	#define OUTPUT_H
	void print_Logo();
	void OutputVTK(int&, LBM&); //For fluid
	
	#if defined Special_Case
		void write_fluid_profile(int, LBM&);
	#endif

	#if defined IBM
		void OutputVTK_Marker(int&, vector<MARKER>&); // For boundary marker
	#endif

	//Analytical result
	#if defined Couette_2D
		//double Couette_Analytic(int step) ;
		double Error_Calculation(LBM &lb, int step);
	#elif defined TaylorGreen_2D //need for checking : velocity, pressure, sigma_xx and sigma_xy
		//Analytic result
		double TaylorGreen_Analytic_ux(double t, double x, double y);
		double TaylorGreen_Analytic_uy(double t, double x, double y);
		double TaylorGreen_Analytic_p(double t, double x, double y);
		double TaylorGreen_Analytic_sxx(double t, double x, double y);
		double TaylorGreen_Analytic_sxy(double t, double x, double y);

		//Error calculation 
		double Error_Calculation_TG_u(LBM &lb, int step);
		double Error_Calculation_TG_p(LBM &lb, int step);
		double Error_Calculation_TG_sxx(LBM &lb, int step);
		double Error_Calculation_TG_sxy(LBM &lb, int step);
	#elif defined Poiseuille_2D
		double Poiseuille_Analytic (double t, double y);
		double Error_Calculation_Poiseuille(LBM &lb, int step);
	#endif

	//Momentum Exchange Algorithm for stationary boundary : 
	#if defined Cylinder_Basic or defined Cylinder_IBM 
		double MEA_CL(LBM &lb);
		double MEA_CD(LBM &lb);
	#endif

	//Lift and Drag coefficient for moving boundary :
	#if defined Cylinder_Translationally_Oscillating
		double Calc_CL(LBM &lb, double&, double&);
		double Calc_CD(LBM &lb, double&, double&);
	#elif defined Cylinder_Rotationally_Oscillating
		double Calc_CT(LBM &lb, double&, double&);
	#endif


	#if defined TaylorGreen_2D 
		void OutputCSV_TG_V(int &nout, double *error); //velocity
		void OutputCSV_TG_p(int &nout, double *error); //pressure
		void OutputCSV_TG_sxx(int &nout, double *error); //sigma xx
		void OutputCSV_TG_sxy(int &nout, double *error); //sigma xy
	#elif defined Cylinder_Basic or defined Cylinder_IBM or defined Cylinder_Translationally_Oscillating
		void OutputCSV_Cylinder_CL(int nout, vector<double>& error); //Output CL
		void OutputCSV_Cylinder_CD(int nout, vector<double>& error); //Output CD
	#elif defined Cylinder_Rotationally_Oscillating
		void OutputCSV_Cylinder_CT(int nout, vector<double>& error); //Output CT
	#else
		void OutputCSV(int &nout, double *error); //only for velocity
	#endif
	

#endif

//-----------------------------------------------------------------------------------------------------------------------