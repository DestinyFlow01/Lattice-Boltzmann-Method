//Definition of member function in LBM.h
#include<iostream>
#include<stdio.h>
#include<cmath>
#include<omp.h>
#include<vector>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "LBM.h"
using namespace std;

LBM::LBM(int Nx, int Ny, int Nz, double tau) : Nx(Nx), Ny(Ny), Nz(Nz), tau(tau) {
	double omega = dt/tau;
	//Memory allocation for LATTICE **fluid
	fluid = new LATTICE **[Nx+2];
	
	for (int i = 0; i<Nx+2; i++) {
		fluid[i] = new LATTICE *[Ny+2];
		for(int j = 0; j<Ny+2; j++) {
			fluid[i][j] = new LATTICE [Nz+2];
		}
	}
}


void LBM::Init() { //Initialization using equilibrium 
	#pragma omp parallel for 
	
	
	for (int i = 0; i<=Nx+1; i++) { //x direction
		for (int j = 0; j<=Ny+1; j++) {  //y direction
			for(int k = 0; k<=Nz+1; k++) { // z direction
				if(fluid[i][j][k].type == TYPE_F || fluid[i][j][k].type == TYPE_E || fluid[i][j][k].type == TYPE_L)  {
					double u_dot_u = pow(fluid[i][j][k].ux,2) + pow(fluid[i][j][k].uy,2) + pow(fluid[i][j][k].uz,2);
					//cout<<u_dot_u<<endl;
					for (int l = 0; l<npop; l++) {
						double c_dot_u = fluid[i][j][k].ux*cx[l] + fluid[i][j][k].uy*cy[l] + fluid[i][j][k].uz*cz[l];
						double f_eq = fluid[i][j][k].density*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));
						fluid[i][j][k].f[l] = f_eq;
						
						fluid[i][j][k].f_star[l] = f_eq;
					}
				}

				#if defined Poiseuille_2D //Force density = pressure gradient -dp/dx, initial force profile
					fluid[i][j][k].F[0] = -dp_dx; //x direction force Fx, totallny will be positive to accelerate the flow
				#endif
			}
		}
	}
	
}

//Terminating in Collision function
void LBM::Collision() {
	double omega = dt/tau;
	#pragma omp parallel for 
	
	for(int i = 0; i<=Nx+1; i++) {
		for (int j = 0; j<=Ny+1; j++) {
			for(int k = 0; k<=Nz+1; k++) {
				if(fluid[i][j][k].type == TYPE_F or fluid[i][j][k].type == TYPE_E) {
					//Macroscopic property of the fluid 
					
					/*
					double rho = 0, rho_ux = 0, rho_uy = 0, rho_uz = 0;
					for(int l = 0; l<npop; l++) {
						rho += fluid[i][j][k].f[l];
						rho_ux += fluid[i][j][k].f[l]*cx[l];
						rho_uy += fluid[i][j][k].f[l]*cy[l];
						rho_uz += fluid[i][j][k].f[l]*cz[l];
					}
					
					
					double ux = rho_ux/rho;
					double uy = rho_uy/rho;
					double uz = rho_uz/rho;
					*/
					
					double rho = fluid[i][j][k].density;
					double ux = fluid[i][j][k].ux;
					double uy = fluid[i][j][k].uy;
					double uz = fluid[i][j][k].uz;
					double Fx = fluid[i][j][k].F[0];
					double Fy = fluid[i][j][k].F[1];
					double Fz = fluid[i][j][k].F[2];
					
					
					if(j == Ny) {
						//cout<<"ux = "<<ux<<endl;
					}
					

					//Equilibrium distribution 
					double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
					double F_dot_u = Fx*ux + Fy*uy + Fz*uz;
					for(int l = 0; l<npop; l++) {
						//Preparation for collision term
						double c_dot_u = ux*cx[l] + uy*cy[l] + uz*cz[l];
						double f_eq = rho*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	

						//Preparation for source term
						double F_dot_c = Fx*cx[l] + Fy*cy[l] + Fz*cz[l];
						double F = w[l]*(F_dot_c/cs2 + (F_dot_c)*(c_dot_u)/(cs2*cs2) - F_dot_u/cs2 );
						double S = (1-dt/(2*tau))*F;

						//Collision process
						fluid[i][j][k].f_star[l] = fluid[i][j][k].f[l]*(1-omega) + omega*f_eq + S*dt;
						if(f_eq <0) {
							cout<<"Terminating in Collision function"<<endl;
							terminate();
						}
					}
				}
			}
		}
	}
	//cout<<fluid[Nx/2][Ny/2][1].f[2]<<endl;
}

//will be updated for IBM
void LBM::Streaming() {
    #pragma omp parallel for
    for(int i=0; i<=Nx+1; ++i)
    {
        for(int j=0; j<=Ny+1; ++j)
        {
            for(int k = 0; k<=Nz+1; ++k)
            {
				
				if(fluid[i][j][k].type==TYPE_F) {
					int i_nb, j_nb, k_nb;
					for (int l=0; l < npop; ++l) {
						i_nb = (i - cx[l] + (Nx+2))%(Nx+2) ;
						j_nb = (j - cy[l] + (Ny+2))%(Ny+2) ;
						k_nb = (k - cz[l] + (Nz+2))%(Nz+2) ;
						//---- Solid Boundary Condition ----------------------
						if(fluid[i_nb][j_nb][k_nb].type==TYPE_S) {
							#if defined Couette_2D 
								if(j == Ny) {
									double uwall_dot_c = u_wall*cx[l];
									fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]] + 2*w[l]*fluid[i][j][k].density*uwall_dot_c/cs2;
								}
								else {
									fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]];
								}

							#elif defined Cylinder_Basic
								//Interpolated bounce back only on cylinder

								//distance at i,j,k (solid)
								double rx1 = dx*i_nb-dx*0.5;
								double ry1 = dy*j_nb-dy*0.5;
								double distance1 = sqrt(pow(rx1-center_x,2.0) + pow(ry1 - center_y,2.0));
								//cout<<"rx1 = "<<rx1<<" ry1 = "<<ry1<<" Awal : \n";
								if(distance1<=0.5*D) {
									//cout<<"Masuk : \n";
									//distance at i,j,k (fluid)
									double rx2 = dx*i-dx*0.5;
									double ry2 = dy*j-dy*0.5;
									double distance2 = sqrt(pow(rx2-center_x,2.0) + pow(ry2 - center_y,2.0));

									//Calculating q :
									//double q = abs(distance2 - 0.5*D)/(sqrt(pow(cx[l]*dx,2) + pow(cy[l]*dy,2))*dt);
									double q = 0.5;
									//cout<<"q = "<<q<<endl;

									//Boundary condition
									if(q <= 0.5) {
										int i_f = (i + cx[l] + (Nx+2))%(Nx+2) ;
										int j_f = (j + cy[l] + (Ny+2))%(Ny+2) ;
										int k_f = (k + cz[l] + (Nz+2))%(Nz+2) ;
										fluid[i][j][k].f[l] = 2*q*fluid[i][j][k].f_star[opposite[l]] + (1-2*q)*fluid[i_f][j_f][k_f].f_star[opposite[l]];
									}
										
									else if (q<=1) 
										fluid[i][j][k].f[l] = 1.0/(2.0*q)*fluid[i][j][k].f_star[opposite[l]] + (2.0*q-1.0)/(2*q)*fluid[i][j][k].f_star[l];
									else {
										cout<<"q = "<<q<<", q >= 1"<<endl; terminate();
									}
								}
								else {
									fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]];
								}
								

								
									


							#else
								fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]];
							#endif

						}
						
						
						
						//---- Inlet/Outlet Boundary Condition
						else if (fluid[i_nb][j_nb][k_nb].type==TYPE_E)
						{
							fluid[i_nb][j_nb][k_nb].f[l] = fluid[i][j][k].f_star[l];
						}
						
						//Slip boundary condition
						else if(fluid[i_nb][j_nb][k_nb].type==TYPE_L) {
							fluid[i_nb][j_nb][k_nb].f[l] = fluid[i][j][k].f_star[slip[l]];
						}
						else //---- Periodic Boundary Condition --------------------
						{
							//if(i!= Nx+1 or i!=0) {
								i_nb = ((i_nb - 1 + (Nx)) % (Nx)) + 1;
								j_nb = ((j_nb - 1 + (Ny)) % (Ny)) + 1;
								k_nb = ((k_nb - 1 + (Nz)) % (Nz)) + 1;
								fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].f_star[l];
							//}
						}
						
					}
				}
            }
        }
    } 
}


void LBM::MacroProp() {
	#pragma omp parallel for 
	
	for (int i = 1; i<=Nx; i++) {
		for(int j = 1; j<=Ny; j++) {
			for(int k = 1; k<=Nz; k++) {
				
				if(fluid[i][j][k].type == TYPE_F) {
					double rho = 0, rho_ux = 0, rho_uy = 0, rho_uz = 0;
					for(int l = 0; l<npop; l++) {
						rho += fluid[i][j][k].f[l];
						rho_ux += fluid[i][j][k].f[l]*cx[l];
						rho_uy += fluid[i][j][k].f[l]*cy[l];
						rho_uz += fluid[i][j][k].f[l]*cz[l];
					}
					
					fluid[i][j][k].density = rho;
					fluid[i][j][k].ux = (rho_ux + 0.5*dt*fluid[i][j][k].F[0])/rho;
					fluid[i][j][k].uy = (rho_uy + 0.5*dt*fluid[i][j][k].F[1])/rho;
					fluid[i][j][k].uz = (rho_uz + 0.5*dt*fluid[i][j][k].F[2])/rho;
					//cout<<"ux = "<<fluid[i][j][k].ux<<endl;
					
				}
				else {
					fluid[i][j][k].density = 0;
					fluid[i][j][k].ux = 0;
					fluid[i][j][k].uy = 0;
					fluid[i][j][k].uz = 0;
				}
				if(fluid[i][j][k].density <0 or fluid[i][j][k].density >100 ) {
							cout<<"Terminating in MacroProp function"<<endl;
							terminate();
						}
			}
		}
	}
	//cout<<"ux = "<<fluid[Nx/2][Ny/2][1].ux<<endl;
}

int LBM::getNx() const {
	return Nx;
}

int LBM::getNy() const {
	return Ny;
}

int LBM::getNz() const {
	return Nz;
}

double LBM::getTAU() const{
	return tau;
}


//------------------------------------------STEPS FOR IMMERSED BOUNDARY METHOD---------------------------------
//steps for Multi Direct Forcing - Immersed Boundary Method (MDF-IB-LBM) (RIGID BOUNDARY): 
	//1. Compute Lagrangian Force from current boundary configuration (Varies for rigid and deformable)
	//2. Spread Lagrangian force to lattice via Eulerian force density (remember to use the stencil of Dirac delta kernel function) (DONE)
	//3. Compute uncorrected velocity (DONE)
	//4. LBM (Collision, streaming) with forcing (DONE IN PREVIOUS PART)
	//5. compute macroscopic properties. (DONE IN PREVIOUS PART)
	//6. Interpolate fluid velocity at Lagrangian node to obtain the speed of marker (DONE)
	//7. Update marker configuration
	//Step 6 and 7 doesn't need to be done for rigid boundary

#if defined IBM
	
	double kernel(double x, double delta) { //x = distance
		if(stencil == 2) {
			return 1-abs(x/delta);
		}
		
		else if (stencil == 3) {
			if(x <= 0.5*delta) return ( 1 + sqrt(1-3*pow(x/delta,2)))/3.0;
			else if(x<=1.5*delta) return ( 5 - 3*abs(x/delta) - sqrt(-3*pow(x/delta,2) + 6*abs(x/delta) - 2) )/6.0 ;
		}
		
		else if (stencil == 4) {
			if(x <= delta) return ( 3 - 2*abs(x/delta) + sqrt(1+4*abs(x/delta) - 4*pow(x/delta,2)) )/8.0;
			else if(x<=2*delta) return ( 5 - 2*abs(x/delta) - sqrt(-7 + 12*abs(x/delta) - 4*pow(x/delta,2)))/8.0;
		}
	}

	//Step 0 : Compute the new position of the center 
	#if defined Cylinder_Translationally_Oscillating
		void Update_Position_COM(vector<MARKER>& marker, double& t, double* COM_Velocity_Old, double* COM_Velocity) {
		//changing position from new to old velocity
		COM_Velocity_Old[0] = COM_Velocity[0];
		COM_Velocity_Old[1] = COM_Velocity[1];
		COM_Velocity_Old[2] = COM_Velocity[2];

		//Changing new velocity part (Externally prescribed motion) : 
		#if defined Cylinder_Translationally_Oscillating
			const double PI = 3.14159265;
			COM_Velocity[0] = -u_max*cos(2*PI*t/T);
			COM_Velocity[1] = 0; COM_Velocity[2] = 0;
		#endif
		

		//boundary velocity
		//Translation effect
		
		for (int n = 0; n<number_nodes; n++) {
			marker[n].u_boundary[0] = COM_Velocity[0]; //x component
			marker[n].u_boundary[1] = COM_Velocity[1]; //y component
			marker[n].u_boundary[2] = COM_Velocity[2]; //z component

		}
		

		//Rotation effect

		//changing the position of the center of mass
		center_x += COM_Velocity[0]*dt; center_y += COM_Velocity[1]*dt;

		//generating new coordinate of marker
		double radius = D/2.0;
		double x,y;
		//#pragma omp parallel for 
		
		#if defined Cylinder_Translationally_Oscillating
			for (int i = 0; i<number_nodes; i++) {
				x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
				marker[i].x = x; 
				marker[i].y = y;
			}
		#endif


		}
	#elif defined Cylinder_Rotationally_Oscillating
		void Update_Position_COM(vector<MARKER>& marker, double& t, double* COM_Velocity_Old, double* COM_Velocity, double& Omega_Old, double& Omega, double& theta) { //generalize only for 2D
			//changing position from new to old velocity
			COM_Velocity_Old[0] = COM_Velocity[0];
			COM_Velocity_Old[1] = COM_Velocity[1];
			COM_Velocity_Old[2] = COM_Velocity[2];

			Omega_Old = Omega;

			//Changing new velocity part (Externally prescribed motion) : 
			COM_Velocity[0] = 0; COM_Velocity[1] = 0; COM_Velocity[2] = 0;
			Omega = Omega_max*sin(2*PI*t/T);

			//Translation effect
			for (int n = 0; n<number_nodes; n++) {
				marker[n].u_boundary[0] = COM_Velocity[0]; //x component
				marker[n].u_boundary[1] = COM_Velocity[1]; //y component
				marker[n].u_boundary[2] = COM_Velocity[2]; //z component
			}

			//Rotation effect
			for (int n = 0; n<number_nodes; n++) {
				marker[n].u_boundary[0] += -Omega*(marker[n].y - center_y);
				marker[n].u_boundary[1] += Omega*(marker[n].x - center_x);
			}


			//changing the position of the center of mass
			center_x += COM_Velocity[0]*dt; center_y += COM_Velocity[1]*dt;

			//changing the angle of the cylinder
			theta += Omega*dt;

			//generating new coordinate of marker
			double radius = D/2.0;
			double x,y;
			
			for (int i = 0; i<number_nodes; i++) {
				x = radius*cos(2*PI/number_nodes*i + theta) + center_x  +0.5*dx; y = radius*sin(2*PI/number_nodes*i + theta) + center_y + 0.5*dy;
				marker[i].x = x; 
				marker[i].y = y;
			}
		}
	#endif
	

	//Step 1 : Compute Lagrangian Force
	void Compute_Lagrangian(vector<MARKER>& marker, LBM lb_copy, LBM& lb) {
		#if defined EXPLICIT_DIRECT_FORCING
			#pragma omp parallel for
			for (int n = 0; n<number_nodes; n++) {
				//Delta r
				double delta_r_x = marker[n].x - marker[n].x_ref;
				double delta_r_y = marker[n].y - marker[n].y_ref;
				double delta_r_z = marker[n].z - marker[n].z_ref;

				double delta_r = sqrt(pow(delta_r_x,2) + pow(delta_r_y,2) + pow(delta_r_z,2));
				

				//Force : 
				if(delta_r > 0) {
					marker[n].F[0] = -stiffness*delta_r_x/delta_r;
					marker[n].F[1] = -stiffness*delta_r_y/delta_r;
					marker[n].F[2] = -stiffness*delta_r_z/delta_r;
				}
				else {
					marker[n].F[0] = 0;
					marker[n].F[1] = 0;
					marker[n].F[2] = 0;
				}

			}
		#elif defined MULTI_DIRECT_FORCING
			//To calculate initial volume force at boundary point
			//Calculate Uncorrected Velocity
			Compute_Uncorrected_Velocity(lb_copy);

			//Interpolate to boundary velocity
			Interpolate(lb_copy, marker);

			//Initial volume force for l = 0
			for (int n = 0; n<number_nodes; n++) {
				marker[n].F[0] = (marker[n].u_boundary[0] - marker[n].rdot[0])/dx;
				marker[n].F[1] = (marker[n].u_boundary[1] - marker[n].rdot[1])/dx;
				marker[n].F[2] = (marker[n].u_boundary[2] - marker[n].rdot[2])/dx;
			}
					
			//Iteration step 
			for (int l = 1; l<=m_max; l++) {
				if(l != 1) {
					for (int n = 0; n<number_nodes; n++) {
						marker[n].F[0] += (marker[n].u_boundary[0] - marker[n].rdot[0])/dx;
						marker[n].F[1] += (marker[n].u_boundary[1] - marker[n].rdot[1])/dx;
						marker[n].F[2] += (marker[n].u_boundary[2] - marker[n].rdot[2])/dx;
					}
					
				}

				//Step 1 : Compute volume force at the lattice points during lth iteration
				Spread_Force(lb_copy, marker);

				//Step 2 : Correct Flow velocity
				MacroProp_IBM(lb_copy);

				//Step 3 : Interpolate flow velocity
				Interpolate(lb_copy, marker);
			}

			//copying to the real lb
			for (int i = 1; i<=Nx; i++) {
				for (int j = 1; j<=Ny; j++) {
					for (int k = 1; k<=Nz; k++) {
						lb.fluid[i][j][k].F[0] = lb_copy.fluid[i][j][k].F[0];
						lb.fluid[i][j][k].F[1] = lb_copy.fluid[i][j][k].F[1];
						lb.fluid[i][j][k].F[2] = lb_copy.fluid[i][j][k].F[2];
					}
				}
			}

		#endif
		
	}
	
	//Step 2 : Spread Force to Eulerian Force
	void Spread_Force(LBM& lb, vector<MARKER>& marker) {
		#pragma omp parallel for 
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					double x_lattice = i + 0.5, y_lattice = j + 0.5;
					//if(lb.fluid[i][j][k].inside == 0) {
						//Initialization
						lb.fluid[i][j][k].F[0] = 0; lb.fluid[i][j][k].F[1] = 0; lb.fluid[i][j][k].F[2] = 0;

						//Calculation
						for (int n = 0; n<number_nodes; n++) {
							//Distance 
							double x_marker = marker[n].x, y_marker = marker[n].y;
							double delta_x = x_lattice - x_marker; 
							double delta_y = y_lattice - y_marker;

							//Interpolation
							if(stencil == 2) {
								if(abs(delta_x) <= dx and abs(delta_y) <= dy) {
									lb.fluid[i][j][k].F[0] += marker[n].F[0]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[1] += marker[n].F[1]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[2] += marker[n].F[2]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
								}

								
							}
							else if (stencil == 3) {
								if(abs(delta_x) <= 1.5*dx and abs(delta_y) <= 1.5*dy) {
									lb.fluid[i][j][k].F[0] += marker[n].F[0]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[1] += marker[n].F[1]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[2] += marker[n].F[2]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
								}

							}
							else if (stencil == 4) {
								if(abs(delta_x) <= 2*dx and abs(delta_y) <= 2*dy) {
									lb.fluid[i][j][k].F[0] += marker[n].F[0]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[1] += marker[n].F[1]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									lb.fluid[i][j][k].F[2] += marker[n].F[2]*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
								}
							}
						}
					//}
					
				}
			}
		}
	}

	//Step 3 : Compute uncorrected velocity 
	void Compute_Uncorrected_Velocity(LBM& lb) {
		#pragma omp parallel for 
	
		for (int i = 1; i<=Nx; i++) {
			for(int j = 1; j<=Ny; j++) {
				for(int k = 1; k<=Nz; k++) {
					
					//if(lb.fluid[i][j][k].inside == 0) {
						double rho = 0, rho_ux = 0, rho_uy = 0, rho_uz = 0;
						for(int l = 0; l<npop; l++) {
							rho += lb.fluid[i][j][k].f[l];
							rho_ux += lb.fluid[i][j][k].f[l]*cx[l];
							rho_uy += lb.fluid[i][j][k].f[l]*cy[l];
							rho_uz += lb.fluid[i][j][k].f[l]*cz[l];
						}
						
						lb.fluid[i][j][k].density = rho;
						lb.fluid[i][j][k].ux = (rho_ux)/rho;
						lb.fluid[i][j][k].uy = (rho_uy)/rho;
						lb.fluid[i][j][k].uz = (rho_uz)/rho;
						
						
					//}
					
					
				}
			}
		}
	}

	//Step 4 : LBM 					
	//Step 5 : Macroprop
	void MacroProp_IBM(LBM& lb) {
		#pragma omp parallel for 

		for(int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					lb.fluid[i][j][k].ux += lb.fluid[i][j][k].F[0]*dt/(2*lb.fluid[i][j][k].density);
					lb.fluid[i][j][k].uy += lb.fluid[i][j][k].F[1]*dt/(2*lb.fluid[i][j][k].density);
					lb.fluid[i][j][k].uz += lb.fluid[i][j][k].F[2]*dt/(2*lb.fluid[i][j][k].density);
				}
			}
		}
	}

	//Step 6 : Interpolate fluid velocity to node position
	void Interpolate(LBM& lb, vector<MARKER>& marker) {
		#pragma omp parallel for
		for (int n = 0; n<number_nodes; n++) {
			//Initialization
			marker[n].rdot[0] = 0; marker[n].rdot[1] = 0; marker[n].rdot[2] = 0;
			double x_marker = marker[n].x; double y_marker = marker[n].y; 

			//Calculation
			for(int i = 1; i<=Nx; i++) {
				for (int j = 1; j<=Ny; j++) {
					for (int k = 1; k<=Nz; k++) {
						
						//if(lb.fluid[i][j][k].inside == 0) {
							//Distance 
							double x_lattice = i + 0.5, y_lattice = j + 0.5;
							double delta_x = x_lattice - x_marker; 
							double delta_y = y_lattice - y_marker;

							//Interpolation
							if(stencil == 2) {
								if(abs(delta_x) <= dx and abs(delta_y) <= dy) {
									marker[n].rdot[0] += pow(dx,3)*lb.fluid[i][j][k].ux*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[1] += pow(dx,3)*lb.fluid[i][j][k].uy*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[2] += pow(dx,3)*lb.fluid[i][j][k].uz*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
								}

								
							}
							else if (stencil == 3) {
								if(abs(delta_x) <= 1.5*dx and abs(delta_y) <= 1.5*dy) {
									marker[n].rdot[0] += pow(dx,3)*lb.fluid[i][j][k].ux*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[1] += pow(dx,3)*lb.fluid[i][j][k].uy*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[2] += pow(dx,3)*lb.fluid[i][j][k].uz*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
								}

							}
							else if (stencil == 4) {
								if(abs(delta_x) <= 2*dx and abs(delta_y) <= 2*dy) {
									//cout<<"Delta x = "<<delta_x<<", Delta y = "<<delta_y<<endl;
									marker[n].rdot[0] += pow(dx,3)*lb.fluid[i][j][k].ux*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[1] += pow(dx,3)*lb.fluid[i][j][k].uy*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									marker[n].rdot[2] += pow(dx,3)*lb.fluid[i][j][k].uz*kernel(abs(delta_x), dx)*kernel(abs(delta_y), dy)/pow(dx,ndim);
									
								}
							}

						//}
						
					}
				}
			}
			//cout<<"Marker r_dot "<<marker[n].rdot[0]<<endl;
		}
	}

	//Step 7 : New boundary configuration
	//Need to be updated to account for deformable boundary and the code to determine whether the point is inside or outside the boundary
	void Update_Position(vector<MARKER>& marker) {
		#pragma omp parallel for
		for(int n = 0; n<= number_nodes; n++) {
			marker[n].x += marker[n].rdot[0]*dt;
			marker[n].y += marker[n].rdot[1]*dt;
			marker[n].z += marker[n].rdot[2]*dt;
		}
	}

#endif
//-------------------------------------------------------------------------------------------------------------




//------------------------------------------------SETUP FOR SIMULATION ----------------------------------------


#if defined Cylinder_Basic
	
	//Something wrong with the generate cylinder function 
	void Generate_Cylinder(LBM& lb, int& count) {
		double radius = D/2.0;

		#pragma omp parallel for 
		for(int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 0; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0)) <= radius) {
						//cout<<"Type before= "<<lb.fluid[i][j][k].type<<endl;
						lb.fluid[i][j][k].type = TYPE_S;
						//cout<<"Type after= "<<lb.fluid[i][j][k].type<<endl;
						count++;
					}

				}
			}	
		}
		cout<<"count = "<<count<<"and Nx*Ny = "<<Nx*Ny<<endl;
	}

	LBM main_setup_Cylinder() { //Cylinder 2D flow Bounce Back 
		LBM lb(Nx, Ny, Nz, tau);
		
		int count = 0;
		
		cout<<"Reynolds number = "<<Re<<", nu = "<<nu<<", u_max = "<<u_max<<endl;

		//generating cylinder and its solid content 
		Generate_Cylinder(lb, count);
		cout<<"Count = "<<count<<endl;
		//cout<<"Type at center = "<<lb.fluid[int (dx*(Nx/5))][int (dy*(Ny/2))][0].type<<endl;
		cout<<"Cylinder generated\n";
		
		count = 0;
		
		#pragma omp parallel for 
		//applying BC on all boundaries : 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //Slip boundary
					if(i == 0 or i == Nx+1) lb.fluid[i][j][k].type = TYPE_E; //Inlet outlet boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_E || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 4*u_max*j*(Ny+1-j)/pow(Ny+1,2);//u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
						//count++;
					}
					
					
				}
			}
		}
		//cout<<"Count = "<<count<<endl;
		
		return lb;
	}

#elif defined Couette_2D 

	LBM main_setup_Couette() {
		LBM lb(Nx, Ny, Nz, tau);
		
		double nu = cs2*(tau - 0.5*dt);
		cout<<"nu = "<<nu<<", u_wall = "<<u_wall<<", tau = "<<tau<<endl;

		#pragma omp parallel for 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
				

					if(i == 0) lb.fluid[i][j][k].type = TYPE_F; //Inlet and outlet
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_P; //Periodic boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //solid boundary
					if( Nz!=1 and (k == 0 or  k == Nz-1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F ) {
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
				}
			}
		}
		
		return lb;

	}

#elif defined TaylorGreen_2D

	LBM main_setup_Taylor_Green_2D() {
		LBM lb(Nx, Ny, Nz, tau);
		cout<<"Decay time : "<<td<<endl;


		#pragma omp parallel for 
		for (int i = 0; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++ ){
				for(int k = 0; k<=Nz+1; k++) {
					
					if(i == 0 or i == Nx+1 or j == 0 or j == Ny+1 or (Nz!=1 and (k == 0 or  k == Nz-1)))  {
						lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					else {
						//initialization of properties : 
						double t = 0, x = i*dx, y = j*dy;
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = TaylorGreen_Analytic_ux(0,x,y);
						lb.fluid[i][j][k].uy = TaylorGreen_Analytic_uy(0,x,y);
						lb.fluid[i][j][k].uz = 0.0;
					}
					
					
				}
			}
		}
		return lb;
	}

#elif defined Poiseuille_2D

	LBM main_setup_Poiseuille_2D() {
		LBM lb(Nx, Ny, Nz, tau);
		
		for(int i = 0; i<=Nx+1; i++) {
			for(int j = 0; j<=Ny+1; j++) {
				for(int k = 0; k<Nz+1; k++) {
					if(i == 0) lb.fluid[i][j][k].type = TYPE_F; //Inlet and outlet
					if(i == Nx+1) lb.fluid[i][j][k].type = TYPE_P; //Periodic boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //solid boundary
					if( Nz!=1 and (k == 0 or  k == Nz-1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F ) {
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
				}
			}
		}
		return lb;
	}

#elif defined Cylinder_IBM
	LBM main_setup_Cylinder_IBM() {
		LBM lb(Nx, Ny, Nz, tau);
		
		cout<<"Reynolds number = "<<Re<<", nu = "<<nu<<", u_max = "<<u_max<<endl;
		
		//applying BC on all boundaries : 
		#pragma omp parallel for 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_L; //Slip boundary
					if(i == 0 or i == Nx+1) lb.fluid[i][j][k].type = TYPE_E; //Inlet outlet boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					
					
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_E || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
					//if(i == 0 or i == 1) {
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;//6*u_max*(j)*(Ny-j)/pow(Ny,2); //u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
						
				}
			}
		}
		return lb;
	}
	
	void main_setup_marker(LBM &lb, vector<MARKER>& marker) { //Generating cylinder
		//Memory allocation
		//MARKER *marker;
		//marker = new MARKER [number_nodes];
		
		//Setup for marker 
		double radius = D/2.0;
		const double PI = 3.14159265;
		double x,y;
		
		for (int i = 0; i<number_nodes; i++) {
			x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].x = x; 
			marker[i].y = y;
			marker[i].x_ref = x; 
			marker[i].y_ref = y; 
		}
		
		//Setup for lattice
		#pragma omp parallel for 
		
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0)) <= radius) {
						lb.fluid[i][j][k].inside = 1;
					}
				}
			}
		}
		
	}

#elif defined Cylinder_Translationally_Oscillating
	LBM main_setup_Cylinder_Translation_Oscillate() {
		LBM lb(Nx, Ny, Nz, tau);
		
		cout<<"Reynolds number = "<<Re<<", nu = "<<nu<<", u_max = "<<u_max<<endl;
		cout<<"Oscillation period = "<<T<<endl;
		
		//applying BC on all boundaries : 
		#pragma omp parallel for 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_L; //Slip boundary
					if(i == 0 or i == Nx+1) lb.fluid[i][j][k].type = TYPE_E; //Inlet outlet boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_E || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
						
				}
			}
		}
		return lb;
	}

	void main_setup_marker_Translation(LBM &lb, vector<MARKER>& marker) {
		//Setup for marker 
		double radius = D/2.0;
		const double PI = 3.14159265;
		double x,y;
		
		for (int i = 0; i<number_nodes; i++) {
			x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].x = x; 
			marker[i].y = y;
			marker[i].x_ref = x; 
			marker[i].y_ref = y; 
		}
	}
#elif defined Cylinder_Rotationally_Oscillating
	LBM main_setup_Cylinder_Rotation_Oscillate() {
		LBM lb(Nx, Ny, Nz, tau);
		
		cout<<"Reynolds number = "<<Re_r<<", nu = "<<nu<<", u_max = "<<u_max<<endl;
		cout<<"Oscillation period = "<<T<<endl;
		
		//applying BC on all boundaries : 
		#pragma omp parallel for 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_L; //Slip boundary
					if(i == 0 or i == Nx+1) lb.fluid[i][j][k].type = TYPE_E; //Inlet outlet boundary
					if( Nz!=1 and (k == 0 or  k == Nz+1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_E || lb.fluid[i][j][k].type == TYPE_L) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = 0.0;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
					}
					
						
				}
			}
		}
		return lb;
	}
	
	void main_setup_marker_Rotation(LBM& lb, vector<MARKER>& marker) {
		//Setup for marker 
		double radius = D/2.0;
		const double PI = 3.14159265;
		double x,y;
		
		for (int i = 0; i<number_nodes; i++) {
			x = radius*cos(2*PI/number_nodes*i) + center_x  +0.5*dx; y = radius*sin(2*PI/number_nodes*i) + center_y + 0.5*dy;
			marker[i].x = x; 
			marker[i].y = y;
			marker[i].x_ref = x; 
			marker[i].y_ref = y; 
		}
	}
#endif

//---------------------------------------------------------------------------------------------------------------
//usual function definition
void print_Logo() {
	cout<<"Flow Diagnostic LBM Club .... Presents \nThe Computational Method : \n";
	cout <<R"(
		    ____                ___________          __        __ 
		   |    |               |          \        |  \      /  |
		   |    |               |     __    \       |   \    /   |
		   |    |               |    |  \    |      |    \  /    |
		   |    |               |    |__/    |      |     \/     |
		   |	|               |           /       |  |\    /|  |
		   |	|               |     __    \       |  | \  / |  |
		   |	|               |    |  \    |      |  |  \/  |  |
		   |	|________       |    |__/    |      |  |      |  |
		   |             |      |            /      |  |      |  |
		   |_____________|      |___________/       |__|      |__|
	)";
	cout<<endl;
}

//Output fluid
void OutputVTK(int &nout, LBM &lb) {
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	// short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	sprintf(filename,"./case 1 field%06d.vtr",nout);
	fp=fopen(filename,"wb");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,Nz);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx+2,Ny+2,Nz);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(Nz);
	#if defined Special_Case
		fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Pressure\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(Nz);
	#endif
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx+2)*(Ny+2)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+2)*(Ny+2)*(Nz);
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+3);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+3);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nz+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

    // Density (cell)
	array_size=1*4*(Nx+2)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				val32=(float)lb.fluid[i][j][k].density; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}

	#if defined Special_Case
		// Pressure (cell)
		array_size=1*4*(Nx+2)*(Ny+2)*(Nz);
		fwrite(&array_size,sizeof(int),1,fp);
		for(k=1;k<=Nz;k++){
			for(j=0;j<Ny+2;j++){
				for(i=0;i<Nx+2;i++){
					val32=(float)lb.fluid[i][j][k].density*cs2; fwrite(&val32,sizeof(float),1,fp);
				}
			}
		}
	#endif

    // Velocity (cell)
	array_size=3*4*(Nx+2)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				val32=(float)lb.fluid[i][j][k].ux; fwrite(&val32,sizeof(float),1,fp);
				val32=(float)lb.fluid[i][j][k].uy; fwrite(&val32,sizeof(float),1,fp);
				val32=0.0; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}

    // CellType (cell)
	array_size=1*4*(Nx+2)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<=Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx+2;i++){
				val32=(int)lb.fluid[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
			}
		}
	}
	// Coordinates (vertices)
	array_size=1*4*(Nx+3);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<Nx+3;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Ny+3);
	fwrite(&array_size,sizeof(int),1,fp);
	for(j=0;j<Ny+3;j++){ val32=(float)(j*dy); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Nz+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=1;k<Nz+1;k++){ val32=(float)(k*dz); fwrite(&val32,sizeof(float),1,fp); }

	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}

#if defined Special_Case
	void write_fluid_profile(int time, LBM& lb) {

		// Create filename.
		stringstream output_filename;
		output_filename << "fluid_t" << time << ".dat";
		ofstream output_file;
		output_file.open(output_filename.str().c_str());
		
		// Write header.
		output_file << "X Y Z density pressure vel_x vel_y vel_z\n";
		
		// Write data.
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					if(lb.fluid[i][j][k].type == TYPE_F) {
						output_file << i << " " << j << " " << k <<" " << lb.fluid[i][j][k].density <<" " << lb.fluid[i][j][k].density * cs2 << " " << lb.fluid[i][j][k].ux <<" "<<lb.fluid[i][j][k].uy << " " <<lb.fluid[i][j][k].uz << "\n";
					}
					
				}
			}
		}


		output_file.close();

		
	}
#endif

#if defined IBM
	void OutputVTK_Marker(int& time, vector<MARKER>& marker) {
		//Create filename 
		stringstream output_filename;
		output_filename << "case 2 wall_t"<< time << ".vtk";
		ofstream output_file;
		output_file.open(output_filename.str().c_str());
		
		// Write VTK header
	    output_file << "# vtk DataFile Version 3.0\n";
	    output_file << "wall_state\n";
	    output_file << "ASCII\n";
	    output_file << "DATASET POLYDATA\n";
	    
	    //Write node positions
	    output_file<< "POINTS " << number_nodes << " float\n";
	    

	    for (int n = 0; n<number_nodes; n++) {
	    	output_file << marker[n].x << " " << marker[n].y<< " " <<marker[n].z << " \n";
		}

	// Write lines between neighboring nodes
	/*
	output_file << "LINES " << number_nodes - 2 << " " << 3 * (number_nodes - 2) << "\n";
	for (int i = 0; i<number_nodes; i++) {
		output_file << "2 " << i%number_nodes << " " << (i + 1) << "\n";
	}
	*/

	output_file << "LINES " << number_nodes << " " << 3 * (number_nodes) << "\n";
      for(int n = 0; n < number_nodes; ++n) {
        output_file << "2 " << n%(number_nodes) << " " << (n + 1)%(number_nodes) << "\n";
      }
      
	
    // Write vertices
    output_file << "VERTICES 1 " << number_nodes + 1 << "\n";
    output_file << number_nodes << " ";

    for(int n = 0; n < number_nodes; n++) {
      output_file << n << " ";
    }

		output_file.close();
	}
	
#endif
//-------------------------------------------------ANALYTICAL RESULTS---------------------------------------------
#if defined Couette_2D
	
	double Error_Calculation(LBM &lb, int step) {
		//----------------------------------------------Analytic calculation------------------------------------
		double u_analytic[Ny+1]; //only at one time 
		double t = step*dt;
		double sum;
		double h = Ny*dy;
		const double PI = 3.14159265;

		for(int i = 0; i<=Ny; i++) {
			sum = 0;
			for(int n = 1; n<=100; n++) {
				sum += pow(-1,n)/n*exp(-nu*pow(n*PI/h,2)*t)*sin(n*PI*i*dy/h);
			}
			u_analytic[i] = u_wall*i*dy/h + 2/PI*sum;
		}

		//------------------------------------------------------------------------------------------------
		
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double u = 0;
		
		for (int j = 0; j<=Ny; j++) {
			u = lb.fluid[Nx][j][Nz].ux;
			sum_error_num += pow( (u - u_analytic[j]) ,2);
			sum_error_denom += pow( (u_analytic[j]) ,2);

			if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
			else sum_error = sqrt(sum_error_num/sum_error_denom);
		}
		return sum_error;
	}

#elif defined TaylorGreen_2D
	//---------------------------------------------Analytical Calculation-----------------------------------------
	double TaylorGreen_Analytic_ux(double t, double x, double y) {
		double ux = u_max*(-sqrt(ky/kx)*cos(kx*x)*sin(ky*y))*exp(-t/td);
		return ux;
	}
	
	double TaylorGreen_Analytic_uy(double t, double x, double y) {
		double uy = u_max*(sqrt(kx/ky)*sin(kx*x)*cos(ky*y))*exp(-t/td);
		return uy;
	}

	double TaylorGreen_Analytic_p(double t, double x, double y) {
		double p = p0 - rho0*pow(u_max,2)/4*(ky/kx*cos(2*kx*x) + kx/ky*cos(2*ky*y))*exp(-2*t/td);
		return p;
	}

	double TaylorGreen_Analytic_sxx(double t, double x, double y) {
		double sigma_xx = 2*rho0*nu*u_max*sqrt(kx*ky)*sin(kx*x)*sin(ky*y)*exp(-t/td);
		return sigma_xx;
	}

	double TaylorGreen_Analytic_sxy(double t, double x, double y) {
		double sigma_xy = rho0*nu*u_max*(sqrt(pow(kx,3)/ky) - sqrt(pow(ky,3)/kx) )*cos(kx*x)*cos(ky*y)*exp(-t/td);
		return sigma_xy;
	}
	//-------------------------------------------------------------------------------------------------------------


	//----------------------------------------------------Error calculation----------------------------------------
	double Error_Calculation_TG_u(LBM &lb, int step) { //velocity
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double u = 0, v = 0, p = 0;
		double u_mag = 0;
		double x, y, u_analytic, v_analytic, u_mag_analytic;
		double t = step*dt;

		
		for (int i = 1; i<=Nx; i++) {
			for (int j =1; j<=Ny; j++) {
				//Numerical result
				u = lb.fluid[i][j][Nz].ux; v = lb.fluid[i][j][Nz].uy;
				u_mag = sqrt(pow(u,2) + pow(v,2));

				//Analytical result : 
				x = i*dx, y = j*dy;
				u_analytic = TaylorGreen_Analytic_ux(t,x,y);
				v_analytic = TaylorGreen_Analytic_uy(t,x,y);
				u_mag_analytic = sqrt(pow(u_analytic,2) + pow(v_analytic,2));

				//error calculation
				sum_error_num += pow( (u_mag - u_mag_analytic) ,2);
				sum_error_denom += pow( (u_mag_analytic) ,2);

				if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
				else sum_error = sqrt(sum_error_num/sum_error_denom);
			}
		}
		return sum_error;
	}

	double Error_Calculation_TG_p(LBM &lb, int step) { //pressure
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double p = 0;
		double x, y, p_analytic;
		double t = step*dt;

		
		for (int i = 1; i<=Nx; i++) {
			for (int j =1; j<=Ny; j++) {
				//Numerical result
				p = lb.fluid[i][j][Nz].density*cs2;

				//Analytical result : 
				x = i*dx, y = j*dy;
				p_analytic = TaylorGreen_Analytic_p(t,x,y);
				
				//error calculation
				sum_error_num += pow( (p - p_analytic) ,2);
				sum_error_denom += pow( (p_analytic) ,2);

				if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
				else sum_error = sqrt(sum_error_num/sum_error_denom);
			}
		}
		return sum_error;
	}

	double Error_Calculation_TG_sxx(LBM &lb, int step) { //pressure
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double sxx;
		double x, y, sxx_analytic;
		double t = step*dt;

		
		for (int i = 1; i<=Nx; i++) {
			for (int j =1; j<=Ny; j++) {
				//Numerical result
				sxx = 0;
				if(lb.fluid[i][j][Nz].type == TYPE_F)  {
					double u_dot_u = pow(lb.fluid[i][j][Nz].ux,2) + pow(lb.fluid[i][j][Nz].uy,2) + pow(lb.fluid[i][j][Nz].uz,2);
					//cout<<u_dot_u<<endl;
					for (int l = 0; l<npop; l++) {
						double c_dot_u = lb.fluid[i][j][Nz].ux*cx[l] + lb.fluid[i][j][Nz].uy*cy[l] + lb.fluid[i][j][Nz].uz*cz[l];
						double f_eq = lb.fluid[i][j][Nz].density*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));
						double f_neq = lb.fluid[i][j][Nz].f[l] - f_eq;
						sxx += cx[l]*cx[l]*f_neq;
					}
				}
				sxx = -(1-dt/(2*tau)) *sxx;

				//Analytical result : 
				x = i*dx, y = j*dy;
				sxx_analytic = TaylorGreen_Analytic_sxx(t,x,y);
				
				//error calculation
				sum_error_num += pow( (sxx - sxx_analytic) ,2);
				sum_error_denom += pow( (sxx_analytic) ,2);

				if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
				else sum_error = sqrt(sum_error_num/sum_error_denom);
			}
		}
		return sum_error;
	}

	double Error_Calculation_TG_sxy(LBM &lb, int step) { //pressure
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double sxy;
		double x, y, sxy_analytic;
		double t = step*dt;

		
		for (int i = 1; i<=Nx; i++) {
			for (int j =1; j<=Ny; j++) {
				//Numerical result
				sxy = 0;
				if(lb.fluid[i][j][Nz].type == TYPE_F)  {
					double u_dot_u = pow(lb.fluid[i][j][Nz].ux,2) + pow(lb.fluid[i][j][Nz].uy,2) + pow(lb.fluid[i][j][Nz].uz,2);
					//cout<<u_dot_u<<endl;
					for (int l = 0; l<npop; l++) {
						double c_dot_u = lb.fluid[i][j][Nz].ux*cx[l] + lb.fluid[i][j][Nz].uy*cy[l] + lb.fluid[i][j][Nz].uz*cz[l];
						double f_eq = lb.fluid[i][j][Nz].density*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));
						double f_neq = lb.fluid[i][j][Nz].f[l] - f_eq;
						sxy += cx[l]*cy[l]*f_neq;
					}
				}
				sxy = -(1-dt/(2*tau)) *sxy;

				//Analytical result : 
				x = i*dx, y = j*dy;
				sxy_analytic = TaylorGreen_Analytic_sxy(t,x,y);
				
				//error calculation
				sum_error_num += pow( (sxy - sxy_analytic) ,2);
				sum_error_denom += pow( (sxy_analytic) ,2);

				if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
				else sum_error = sqrt(sum_error_num/sum_error_denom);
			}
		}
		return sum_error;
	}
	//-------------------------------------------------------------------------------------------------------------

#elif defined Poiseuille_2D 
	double Poiseuille_Analytic (double t, double y) {
			double mu = nu*rho0;
			double h = Ny*dy;
			const double PI = 3.14159265;
			double sum = 0;
			

			for(int n = 1; n<=1000; n++) {
				sum+= 1/pow(n,3) * (1-pow(-1,n))*exp(-nu*pow(n*PI/h,2)*t)*sin(n*PI*y/h);
			}
		double velocity = 1/mu*dp_dx*(0.5*y*(y-h) -2*pow(h,2)/pow(PI,3) *sum);
		return velocity;
	}

	double Error_Calculation_Poiseuille(LBM &lb, int step) {
		//-----------------------------------------Analytic Calculation -------------------------------------
		double u_analytic[Ny+1];
		double t = step*dt;
		double y;

		for(int i = 0; i<=Ny; i++) { 
			y = i*dy;
			u_analytic[i] = Poiseuille_Analytic(t,y);
		}
		
		//---------------------------------------------------------------------------------------------------
		double sum_error_num = 0;
		double sum_error_denom = 0;
		double sum_error = 0;
		double u = 0;
		
		for (int j = 0; j<=Ny; j++) {
			u = lb.fluid[Nx][j][Nz].ux;
			sum_error_num += pow( (u - u_analytic[j]) ,2);
			sum_error_denom += pow( (u_analytic[j]) ,2);

			if(sum_error_denom == 0) sum_error = sqrt(sum_error_num);
			else sum_error = sqrt(sum_error_num/sum_error_denom);
		}
		return sum_error;
		
	}
#endif

//----------------------------------------------------------------------------------------------------------------

//---------------------------------------------------ONLY FOR CYLINDER--------------------------------------------
//Momentum exchange Algorithm : to calculate CL and CD
#if defined Cylinder_Basic or defined Cylinder_IBM
	double MEA_CL(LBM &lb) {
		double Py = 0; //Total momentum exchange in y direction
		//#pragma omp parallel for
		for(int i=1; i<=Nx; ++i) {
			for(int j=1; j<=Ny; ++j) {
				for(int k = 1; k<=Nz; ++k) {

					#if defined Cylinder_IBM
						
						Py += - lb.fluid[i][j][k].F[1]*pow(dx,ndim);
						if(lb.fluid[i][j][k].inside==0) {
							int i_nb, j_nb, k_nb;
							for (int l=0; l < npop; ++l) {
								i_nb = i + cx[l];
								j_nb = j + cy[l];
								k_nb = k + cz[l]; 
								//---- Solid Boundary Condition ----------------------
								
								if(lb.fluid[i_nb][j_nb][k_nb].inside == 1)
								{
									//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] + lb.fluid[i][j][k].f[l]*cy[l]; //for Re = 300
									//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cy[l];										
									//Py += (lb.fluid[i][j][k].f[opposite[l]] + lb.fluid[i][j][k].f_star[l])*cy[l];
									//Py += 2*lb.fluid[i][j][k].f[l]*cy[l];
									//Py += - lb.fluid[i][j][k].F[1]*pow(dx,ndim);
									
								}
							}
						}
						
					#else
						if(lb.fluid[i][j][k].type==TYPE_F) {
							int i_nb, j_nb, k_nb;
							for (int l=0; l < npop; ++l)
							{
								i_nb = i - cx[l];
								j_nb = j - cy[l];
								k_nb = k - cz[l]; 
								//---- Solid Boundary Condition ----------------------
								
								if(lb.fluid[i_nb][j_nb][k_nb].type==TYPE_S)
								{
									//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] + lb.fluid[i][j][k].f[l]*cy[l];
									//Py += lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cy[l]; //usually used
									Py += (lb.fluid[i][j][k].f_star[opposite[l]] + lb.fluid[i][j][k].f[l])*cy[opposite[l]];

									/*
									double term1 = lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] + lb.fluid[i][j][k].f[l]*cy[l];
									double term2 = lb.fluid[i][j][k].f[opposite[l]]*cy[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cy[l];
									Py +=0.5*(term1 + term2);
									*/
								}
							}
						}
					#endif
					
				}
			}
		} 
		
		#if defined Cylinder_IBM
			const double PI = 3.14159265;
			//Py *= PI*D/number_nodes;
		#else
			Py *= pow(dx,3)/dt;
		#endif

		double CL = Py/(0.5*rho0*pow(u_max,2)*D); //Here, u_max is the maximum velocity located at the center of the inlet
		return CL;
	}
	
	
	double MEA_CD(LBM &lb) {
		double Px = 0; //Total momentum exchange in x direction
		//#pragma omp parallel for 
		for(int i=1; i<=Nx; ++i) {
			for(int j=1; j<=Ny; ++j) {
				for(int k = 1; k<=Nz; ++k) {
					
					#if defined Cylinder_IBM
						Px += - lb.fluid[i][j][k].F[0]*pow(dx,ndim);
						if(lb.fluid[i][j][k].inside == 0) {
							int i_nb, j_nb, k_nb;
							for (int l=0; l < npop; ++l)
							{
								i_nb = i + cx[l];
								j_nb = j + cy[l];
								k_nb = k + cz[l]; 
								//---- Solid Boundary Condition ----------------------
								
								if(lb.fluid[i_nb][j_nb][k_nb].inside == 1)
								{
									//Px += -(lb.fluid[i][j][k].f[opposite[l]]+ lb.fluid[i][j][k].f[l])*cx[l];
									//Px += lb.fluid[i][j][k].f[opposite[l]]*cx[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cx[l]; //usually used
									//Px += 2*lb.fluid[i][j][k].f_star[l]*cx[l];
									//Px += (lb.fluid[i][j][k].f[opposite[l]] + lb.fluid[i][j][k].f_star[l])*cx[l];
									//Px += - lb.fluid[i][j][k].F[0]*pow(dx,ndim);
								}
							}
						}
						
					#else
						if(lb.fluid[i][j][k].type==TYPE_F) {
							int i_nb, j_nb, k_nb;
							for (int l=0; l < npop; ++l)
							{
								i_nb = i - cx[l];
								j_nb = j - cy[l];
								k_nb = k - cz[l]; 
								//---- Solid Boundary Condition ----------------------
								
								if(lb.fluid[i_nb][j_nb][k_nb].type==TYPE_S)
								{
									//Px += lb.fluid[i][j][k].f[opposite[l]]*cx[opposite[l]] - lb.fluid[i][j][k].f_star[l]*cx[l]; //usually used
									Px += (lb.fluid[i][j][k].f_star[opposite[l]] + lb.fluid[i][j][k].f[l])*cx[opposite[l]];
								}
							}
						}
					#endif
					
					
				}
			}
		} 
		#if defined Cylinder_IBM
			const double PI = 3.14159265;
			//Px *= PI*D/number_nodes;
		#else
			Px *= pow(dx,3)/dt;
		#endif
		
		double CD = Px/(0.5*rho0*pow(u_max,2)*D);
		return CD;
	}
#endif

//Volume force method
#if defined Cylinder_Translationally_Oscillating
	double Calc_CL(LBM& lb, double& COM_Velocity_old_Y, double& COM_Velocity_new_Y) {
		double Py = 0;

		//Total force calculation
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					Py += - lb.fluid[i][j][k].F[1]*pow(dx,ndim);
				}
			}
		}

		//Internal force calculation (Feng's rigid body approximation)
		Py += rho_f/rho_b*M*(COM_Velocity_new_Y - COM_Velocity_old_Y)/dt;

		//calculate lift coefficient
		double CL = Py/(0.5*rho0*pow(u_max,2)*D); 
		return CL;

	}

	double Calc_CD(LBM& lb, double& COM_Velocity_old_X, double& COM_Velocity_new_X) {
		double Px = 0;

		//Total force calculation
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					Px += - lb.fluid[i][j][k].F[0]*pow(dx,ndim);
				}
			}
		}

		//Internal force calculation (Feng's rigid body approximation)
		Px += rho_f/rho_b*M*(COM_Velocity_new_X - COM_Velocity_old_X)/dt;

		//calculate lift coefficient
		double CD = Px/(0.5*rho0*pow(u_max,2)*D); 
		return CD;
	}

#elif defined Cylinder_Rotationally_Oscillating
	double Calc_CT(LBM &lb, double& Omega_Old, double& Omega) {
		double Tz = 0;

		//Total torque calculation 
		for (int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					Tz += -(lb.fluid[i][j][k].F[1]*(i - center_x) - lb.fluid[i][j][k].F[0]*(j - center_y) )*pow(dx, ndim);
				}
			}
		}

		//Internal torque calculation (Feng's rigid body approximation)
		//Tz += rho_f/rho_b*I_B*(Omega - Omega_Old)/dt;

		//Torque coefficient calculation
		double CT = Tz/(0.5*pow(u_max,2.0)*pow(D,2.0));
		return CT;
	}
#endif

//----------------------------------------------------------------------------------------------------------------

//--------------------------------------------------OUTPUTING IN CSV----------------------------------------------
#if defined TaylorGreen_2D
	void OutputCSV_TG_V(int &nout, double *error) { //only for velocity
	
		ofstream ofs;
		ofs.open("error for velocity.csv");
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << i * dt * Nt/td<< "," << *(error+i)<<"\n";
		}
		ofs.close();
	}

	void OutputCSV_TG_p(int &nout, double *error) { //only for pressure	
		ofstream ofs;
		ofs.open("error for pressure.csv");
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << i * dt * Nt/td << "," << *(error+i)<<"\n";
		}
		ofs.close();
	}

	void OutputCSV_TG_sxx(int &nout, double *error) { //only for sigma_xx
	
		ofstream ofs;
		ofs.open("error for sigma xx.csv");
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << i * dt * Nt/td << "," << *(error+i)<<"\n";
		}
		ofs.close();
	}

	void OutputCSV_TG_sxy(int &nout, double *error) { //only for sigma_xy
	
		ofstream ofs;
		ofs.open("error for sigma xy.csv");
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << i * dt * Nt/td << "," << *(error+i)<<"\n";
		}
		ofs.close();
	}
#elif defined Cylinder_Basic or defined Cylinder_IBM or defined Cylinder_Translationally_Oscillating
	void OutputCSV_Cylinder_CL(int nout, vector<double>& error) { //Output CL
		ofstream ofs;
		ofs.open("case 1 error for Cl.csv");
		ofs << "t*,Cl\n";
		for (int i = 0; i<nout; i++) {
			#if defined Cylinder_Translationally_Oscillating
				ofs << (i+1)*Nt*dt/T << "," << error[i]<<"\n";
			#else
				ofs << (i+1)*Nt * dt *u_max/(D) << "," << error[i]<<"\n";
			#endif
		}
		ofs.close();
	} 

	void OutputCSV_Cylinder_CD(int nout, vector<double>& error) { //Output CD
		ofstream ofs;
		ofs.open("case 1 error for Cd.csv");
		ofs << "t*,Cd\n";
		for (int i = 0; i<nout; i++) {
			#if defined Cylinder_Translationally_Oscillating
				ofs << (i+1)*Nt*dt/T << "," << error[i]<<"\n";
			#else
				ofs << (i+1)*Nt * dt *u_max/(D) << "," << error[i]<<"\n";
			#endif
			
		}
		ofs.close();
	} 
#elif defined Cylinder_Rotationally_Oscillating
	void OutputCSV_Cylinder_CT(int nout, vector<double>& error) { //Output CT
		ofstream ofs;
		ofs.open("case 2 error for CT.csv");
		ofs << "t*,CT\n";
		for (int i = 0; i<nout; i++) {
			#if defined Cylinder_Rotationally_Oscillating
				ofs << (i+1)*Nt*dt/T << "," << error[i]<<"\n";
			#endif
		}
		ofs.close();
	} 
#else
	void OutputCSV(int &nout, double *error) { //only for velocity
	
		ofstream ofs;
		ofs.open("error.csv");
		ofs << "time,error\n";
		for (int i = 0; i<nout; i++) {
			ofs << i * dt << "," << *(error+i)<<"\n";
		}
		ofs.close();
	}
#endif
//----------------------------------------------------------------------------------------------------------------
