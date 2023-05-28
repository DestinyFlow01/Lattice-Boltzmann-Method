#include<iostream>
#include<ctime>
#include<vector>
//#include "OUTPUT.h"
#include "LBM.h"
using namespace std;



int main() {
	print_Logo();
	//ofstream ofs;
	
	clock_t c_start = clock();
	//Defining variable for simulation : 
	#if defined Couette_2D
		LBM lb = main_setup_Couette();
	#elif defined Cylinder_Basic
		LBM lb = main_setup_Cylinder();
	#elif defined TaylorGreen_2D
		LBM lb= main_setup_Taylor_Green_2D();
	#elif defined Poiseuille_2D 
		LBM lb = main_setup_Poiseuille_2D();
	#elif defined Cylinder_IBM	
		LBM lb = main_setup_Cylinder_IBM();
		vector<MARKER> marker(number_nodes + 1);
		main_setup_marker(lb, marker);
	#elif defined Cylinder_Translationally_Oscillating
		LBM lb = main_setup_Cylinder_Translation_Oscillate();
		vector<MARKER> marker(number_nodes + 1);
		main_setup_marker_Translation(lb, marker);

		//center of mass velocity : 
		double COM_Velocity[3] = {0,0,0}; //New center of mass velocity
		double COM_Velocity_Old[3] = {0,0,0}; //Old center of mass velocity, will be changed in update position program
	#elif defined Cylinder_Rotationally_Oscillating
		LBM lb = main_setup_Cylinder_Rotation_Oscillate();
		vector<MARKER> marker(number_nodes + 1);
		main_setup_marker_Rotation(lb, marker);

		//Angular velocity : 
		double Omega[3] = {0,0,0}; //New angular velocity
		double Omega_Old[3] = {0,0,0}; //Old angular velocity, will be changed in update position program

		//center of mass velocity : 
		double COM_Velocity[3] = {0,0,0}; //New center of mass velocity
		double COM_Velocity_Old[3] = {0,0,0}; //Old center of mass velocity, will be changed in update position program

		//Angle 
		double theta = 0;
	#endif
	

	lb.Init();
	
	//Preliminary Result : 
	int step = 0, n = 0;
	#if defined Special_Case
		write_fluid_profile(step*Nt, lb);
	#endif

	OutputVTK(step,lb); //Fluid
	#if defined IBM
		cout<<"Number of nodes = "<<number_nodes<<endl;
		OutputVTK_Marker(step,marker); //Marker
		cout<<"Stencil = "<<stencil<<endl;
	#endif
	
	int _NUMBER = floor(T_OUT/Nt);

	#if defined TaylorGreen_2D 
		double error_v[_NUMBER]; // velocity
		double error_p[_NUMBER]; // pressure
		double error_sigma_xx[_NUMBER]; // Sigma xx
		double error_sigma_xy[_NUMBER]; // Sigma xy
	#elif defined Cylinder_Basic or defined Cylinder_IBM or defined Cylinder_Translationally_Oscillating
		vector <double> CL; //Lift Coefficient
		vector <double> CD; //Drag Coefficient
		double total_data = 0; //Total data for the simulation right now
	#elif defined Cylinder_Rotationally_Oscillating
		vector <double> CT; //Torque coefficient
		double total_data = 0; //Total data for the simulation right now
	#else
		double error[_NUMBER]; // only for velocity
	#endif
	
	
	//Iterations 
	int numstep = T_OUT/dt;
	cout<<"Total number of step = "<<numstep<<endl;
	for(step = 1; step<=numstep; step++) {
		#if defined IBM
			//Step 0 : Change the center position
			double t = step*dt;
			#if defined Cylinder_Translationally_Oscillating
				Update_Position_COM(marker, t, COM_Velocity_Old, COM_Velocity);
			#elif defined Cylinder_Rotationally_Oscillating
				Update_Position_COM(marker, t, COM_Velocity_Old, COM_Velocity, Omega_Old[2], Omega[2], theta);
			#endif
			
			//Step 1 : Compute Lagrangian force
			Compute_Lagrangian(marker, lb, lb);

			//Step 2 : Spread Lagrangian force to Eulerian force density
			Spread_Force(lb,marker);

			//Step 3 : Compute Uncorrected Velocity
			Compute_Uncorrected_Velocity(lb);
		#endif
		

		//Step 4 : LBM
		//Collision
		lb.Collision();
		
		//Streaming and BC (Since the initialization has been done)
		lb.Streaming();
		
		
		//Step 5 : Update velocity
		#if defined IBM
			MacroProp_IBM(lb);
		#else
			lb.MacroProp();
		#endif
		
		
		
		//Outputing file : 
		if(step%Nt == 0) {
			cout<<lb.fluid[Nx-1][Ny/2][1].density<<endl;
			cout<<lb.fluid[10][Ny/2][1].ux<<endl;
			cout<<lb.fluid[Nx-1][Ny/2][1].f[2]<<endl;
			cout<<lb.fluid[10][Ny/2][1].F[0]<<endl;
			cout<<"Near inlet :  "<<lb.fluid[int(center_x-D/2)][int(center_y)][1].F[0]<<" "<<lb.fluid[int(center_x-D/2)][int(center_y)][1].F[1]<<endl;
			cout<<"Far from inlet :  "<<lb.fluid[int(center_x+D/2)][int(center_y)][1].F[0]<<" "<<lb.fluid[int(center_x+D/2)][int(center_y)][1].F[1]<<endl;
			cout<<"Top of cylinder :  "<<lb.fluid[int(center_x)][int(center_y + D/2)][1].F[0]<<" "<<lb.fluid[int(center_x)][int(center_y + D/2)][1].F[1]<<endl;
			cout<<"Bottom of cylinder :  "<<lb.fluid[int(center_x)][int(center_y - D/2)][1].F[0]<<" "<<lb.fluid[int(center_x)][int(center_y - D/2)][1].F[1]<<endl;

			n++;
			OutputVTK(n,lb); // For fluid
			
			#if defined IBM
				OutputVTK_Marker(n,marker); //Marker
				cout<<marker[number_nodes - 10].rdot[1]<<endl;
			#endif
			
			#if defined Special_Case
				write_fluid_profile(step, lb);
			#endif
			cout<<"Step = "<<step<<endl;
			
			//Additional steps :
			#if defined Couette_2D
				error[step/Nt] = Error_Calculation(lb, step);
				cout<<"L2 error = "<<error[step/Nt]<<endl;

			#elif defined TaylorGreen_2D 
				error_v[step/Nt] = Error_Calculation_TG_u(lb, step);
				cout<<"L2 error for velocity = "<<error_v[step/Nt]<<endl;

				error_p[step/Nt] = Error_Calculation_TG_p(lb, step);
				cout<<"L2 error for pressure = "<<error_p[step/Nt]<<endl;

				error_sigma_xx[step/Nt] = Error_Calculation_TG_sxx(lb, step);
				cout<<"L2 error for sigma_xx = "<<error_sigma_xx[step/Nt]<<endl;

				error_sigma_xy[step/Nt] = Error_Calculation_TG_sxy(lb, step);
				cout<<"L2 error for sigma_xy = "<<error_sigma_xy[step/Nt]<<endl;

			#elif defined Cylinder_Basic || defined Cylinder_IBM
				total_data++;
				
				CL.push_back(MEA_CL(lb));
				cout<<"Lift coefficient at time step "<<step<<" = "<<CL[total_data-1]<<endl;

				CD.push_back(MEA_CD(lb));
				cout<<"Drag coefficient at time step "<<step<<" = "<<CD[total_data-1]<<endl;

				//Output partial solution 
				OutputCSV_Cylinder_CL(total_data, CL); //Output CL
				OutputCSV_Cylinder_CD(total_data, CD); //Output CD
				
			#elif defined Poiseuille_2D 
				error[step/Nt] = Error_Calculation_Poiseuille(lb, step);
				cout<<"L2 error = "<<error[step/Nt]<<endl;
			#elif defined Cylinder_Translationally_Oscillating
				total_data++;
				
				CL.push_back( Calc_CL(lb, COM_Velocity_Old[1], COM_Velocity[1]) );
				cout<<"Lift coefficient at time step "<<step<<" = "<<CL[total_data-1]<<endl;

				CD.push_back( Calc_CD(lb, COM_Velocity_Old[0], COM_Velocity[0]) );
				cout<<"Drag coefficient at time step "<<step<<" = "<<CD[total_data-1]<<endl;

				//Output partial solution 
				OutputCSV_Cylinder_CL(total_data, CL); //Output CL
				OutputCSV_Cylinder_CD(total_data, CD); //Output CD
			#elif defined Cylinder_Rotationally_Oscillating
				total_data++;

				CT.push_back( Calc_CT(lb, Omega[2], Omega_Old[2]) );
				cout<<"Torque coefficient at time step "<<step<<" = "<<CT[total_data-1]<<endl;
				cout<<"Theta = "<<theta*360/(2*PI)<<" with maximum angle = "<<dtheta<<endl;

				//Output partial solution 
				OutputCSV_Cylinder_CT(total_data, CT); //Output CT
			#endif
			
			
			

		}

		#if defined IBM
			//Step 6 : Interpolate velocity
			Interpolate(lb, marker);

			//Step 7 : Update node position 
			//Update_Position(marker); //relative motion with respect to center of mass

		#endif
		
		
		
	}
	
	//--------------------------------------------------MAKING CSV ------------------------------------------
	#if defined Couette_2D or defined Poiseuille_2D
		OutputCSV(_NUMBER, error);
	#elif defined TaylorGreen_2D 
		OutputCSV_TG_V(_NUMBER, error_v); //velocity 
		OutputCSV_TG_p(_NUMBER, error_p); //pressure
		OutputCSV_TG_sxx(_NUMBER, error_sigma_xx); //sigma xx
		OutputCSV_TG_sxy(_NUMBER, error_sigma_xy); //sigma xy
	#elif defined Cylinder_Basic or defined Cylinder_IBM or defined Cylinder_Translationally_Oscillating
		OutputCSV_Cylinder_CL(_NUMBER, CL); //Output CL
		OutputCSV_Cylinder_CD(_NUMBER, CD); //Output CD
	#elif defined Cylinder_Rotationally_Oscillating
		OutputCSV_Cylinder_CT(_NUMBER, CT); //Output CT
	#endif
	
	
	//-------------------------------------------------------------------------------------------------
	clock_t c_end = clock();
	long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
}