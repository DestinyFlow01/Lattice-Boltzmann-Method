#include<iostream>
#include<cmath>
using namespace std;



// -----------------------------------------DEFINITIONS OF SIMULATION PARAMETERS ----------------------------------
#ifndef DEFINES_H
	#define DEFINES_H
		
	//FLOW PARAMETERS
	#define rho0 1.0
	#define Re 300.0
	#define p0 0.0
	
	//SIMULATION DOMAIN 
	#define BCXM -2.5
	#define BCXP 7.5
	#define BCYM -2.5
	#define BCYP 2.5

	//SIMULATION PARAMETERS
	static const int Nx = 200;
	static const int Ny = 100;
	static const int Nz = 1;
	static const double tau = 1;
	#define dt 1.0
	#define dx 1.0
	#define dy 1.0
	#define dz 1.0
	static const double cs2 = pow(dx,2)/(3*pow(dt,2));
	
	
	
	
	//OUTPUT VALUE
	#define T_OUT 1000
	#define Nt 100
	
	//FLOW CASES
	//#define Couette_2D
	//#define TaylorGreen_2D
	#define Cylinder_Basic
	//#define Poiseuille_2D
	//#define Cylinder_IBM
	//#define Tandem_2_Cylinder
	
	//For Couette flow
	#if defined Couette_2D 
		#define u_wall 0.4
	#else
		#define u_wall 0
	#endif


	//VELOCITY SET 
	#define D2Q9
	
	//BOUNDARY CONDITION
	#define TYPE_F 0 // fluid domain
    #define TYPE_S 1 // (stationary or moving) solid boundary
    #define TYPE_E 2 // equilibrium boundary (inflow/outflow)
    #define TYPE_P 3 // periodic boundary
	
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
		const double cx[9] = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0};
		const double cy[9] = {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0};
		const double cz[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
		//weight function 
		const double w[9] = {4.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/9, 1.0/36, 1.0/36, 1.0/36, 1.0/36};
		//opposite for D2Q9 velocity set : 
		const int opposite[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
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
	};
	
	//LBM class definition 
	class LBM {
		private : 
			int Nx, Ny, Nz = 1;
			double tau;
		public : 
			LATTICE ***fluid; //Definition of lattice, accessible to all
			LBM(int Nx, int Ny, int Nz, double tau); //Parametric constructor
			
			//Mutator
			void Init();
			LATTICE*** Collision(); 
			LATTICE*** Streaming(); 
			LATTICE*** MacroProp();
			
			//Accessor
			int getNx() const;
			int getNy() const;
			int getNz() const;
			double getTAU() const;
	};
	
#endif
// ----------------------------------------------------------------------------------------------------------------------








// --------------------------------------SETUP OF THE PROBLEM -----------------------------------------------------------
#ifndef SETUP_H
	#define SETUP_H
	
	#if defined Couette_2D
		LBM main_setup_Couette();
		
	#elif defined TaylorGreen_2D
		LBM main_setup_Taylor_Green_2D();
	#elif defined Cylinder_Basic || defined Cylinder_IBM
		LBM main_setup_Cylinder();
	#elif defined Poiseuille_2D
		LBM main_setup_Poiseuille_2D();
	#elif defined Tandem_2_Cylinder
	
	#endif

#endif


//-----------------------------------------------------------------------------------------------------------------------











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
	
	#if defined Cylinder_Basic || defined Cylinder_IBM 
		void Generate_Cylinder(LBM& lb, int D, int& count);
	#elif defined Tandem_2_Cylinder 
		
	#endif
	
#endif

// ----------------------------------------------------------------------------------------------------------------------









// -------------------------------------------OUTPUT HEADER--------------------------------------------------------------
#ifndef OUTPUT_H
	#define OUTPUT_H
	void print_Logo();
	void OutputVTK(int&, LBM&);
#endif

//-----------------------------------------------------------------------------------------------------------------------