//Definition of member function in LBM.h
#include<iostream>
#include<cmath>
#include<omp.h>
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
	
	
	for (int i = 1; i<=Nx; i++) { //x direction
		for (int j = 1; j<=Ny; j++) {  //y direction
			for(int k = 1; k<=Nz; k++) { // z direction
				if(fluid[i][j][k].type == TYPE_F) {
					double u_dot_u = pow(fluid[i][j][k].ux,2) + pow(fluid[i][j][k].uy,2) + pow(fluid[i][j][k].uz,2);
					
					for (int l = 0; l<npop; l++) {
						double c_dot_u = fluid[i][j][k].ux*cx[l] + fluid[i][j][k].uy*cy[l] + fluid[i][j][k].uz*cz[l];
						double f_eq = fluid[i][j][k].density*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));
						fluid[i][j][k].f[l] = f_eq;
						fluid[i][j][k].f_star[l] = f_eq;
					}
				}
			}
		}
	}
}

LATTICE*** LBM::Collision() {
	double omega = dt/tau;
	#pragma omp parallel for 
	
	for(int i = 1; i<=Nx; i++) {
		for (int j = 1; j<=Ny; j++) {
			for(int k = 1; k<=Nz; k++) {
				if(fluid[i][j][k].type == TYPE_F) {
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
					
					
					if(j == Ny) {
						//cout<<"ux = "<<ux<<endl;
					}
					

					//Equilibrium distribution 
					double u_dot_u = pow(ux,2) + pow(uy,2) + pow(uz,2);
					for(int l = 0; l<npop; l++) {
						double c_dot_u = ux*cx[l] + uy*cy[l] + uz*cz[l];
						double f_eq = rho*w[l]*(1 + c_dot_u/cs2 + pow(c_dot_u,2)/(2*cs2*cs2) - u_dot_u/(2*cs2));	
						fluid[i][j][k].f_star[l] = fluid[i][j][k].f[l]*(1-omega) + omega*f_eq;
					}
				}
			}
		}
	}
	return fluid;
}




LATTICE*** LBM::Streaming() {
	#pragma omp parallel for 
	

	for(int i = 1; i<=Nx; i++) {
		for (int j = 1; j<=Ny; j++) {
			for(int k = 1; k<=Nz; k++) {
				if(fluid[i][j][k].type == TYPE_F) {
					int i_nb, j_nb, k_nb;
                    for (int l=0; l < npop; ++l)
                    {
                        i_nb = i - cx[l];
                        j_nb = j - cy[l];
                        k_nb = k - cz[l];
                        //---- Solid Boundary Condition ----------------------
                        if(fluid[i_nb][j_nb][k_nb].type==TYPE_S)
                        {
							if(j_nb == Ny+1) {
								double uwall_dot_c = u_wall*cx[l];
								
                            	fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]] + 2*fluid[i][j][k].density*w[l]*(uwall_dot_c)/(cs2);
							}
							else {
								fluid[i][j][k].f[l] = fluid[i][j][k].f_star[opposite[l]];
							}
							
                        }
                        //---- Inlet/Outlet Boundary Condition
                        else if (fluid[i_nb][j_nb][k_nb].type==TYPE_E)
                        {
                            fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].f_star[l];
                        }
                        else //---- Periodic Boundary Condition and usual Streaming--------------------
                        {
                            i_nb = ((i_nb + (Nx)) % (Nx)) + 1;
                            j_nb = ((j_nb + (Ny)) % (Ny)) + 1;
                            k_nb = ((k_nb + (Nz)) % (Nz)) + 1;
                            fluid[i][j][k].f[l] = fluid[i_nb][j_nb][k_nb].f_star[l];
                        }
                    }
				}
			}
		}
	}
	return fluid;
}

LATTICE*** LBM::MacroProp() {
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
					fluid[i][j][k].ux = rho_ux/rho;
					fluid[i][j][k].uy = rho_uy/rho;
					fluid[i][j][k].uz = rho_uz/rho;
					
				}
				else {
					fluid[i][j][k].density = 1;
					fluid[i][j][k].ux = 0;
					fluid[i][j][k].uy = 0;
					fluid[i][j][k].uz = 0;
				}
			}
		}
	}
	return fluid;
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

//------------------------------------------------SETUP FOR SIMULATION ----------------------------------------


#if defined Cylinder_Basic || defined Cylinder_IBM
	
	//Something wrong with the generate cylinder function 
	void Generate_Cylinder(LBM& lb, double D, int& count) {
		double radius = D/2.0;
		double center_x = dx*(Nx/5.0);
		double center_y = dy*(Ny/2.0);

		#pragma omp parallel for 
		for(int i = 1; i<=Nx; i++) {
			for (int j = 1; j<=Ny; j++) {
				for (int k = 1; k<=Nz; k++) {
					
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					if(sqrt(pow(center_x-rx,2.0) + pow(center_y-ry,2.0)) <= radius) {
						lb.fluid[i][j][k].type = TYPE_S;
						count++;
					}

				}
			}	
		}
		cout<<"count = "<<count<<"and Nx*Ny = "<<Nx*Ny<<endl;
	}

	LBM main_setup_Cylinder() { //Cylinder 2D flow Bounce Back 
		LBM lb(Nx, Ny, Nz, tau);
		
		double D = Ny/5;
		double nu = cs2*(tau - 0.5);
		double u_max = Re*nu/D;
		int count = 0;
		
		cout<<"Reynolds number = "<<Re<<", nu = "<<nu<<", u_max = "<<u_max<<endl;

		//generating cylinder and its solid content 
		Generate_Cylinder(lb,D, count);
		cout<<"Count = "<<count<<endl;
		cout<<"Cylinder generated\n";
		

		#pragma omp parallel for 
		//applying BC on all boundaries : 
		count = 0;
		
		double radius = D/2;
		double center_x = dx*(Nx/5);
		double center_y = dy*(Ny/2);
		
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
				
					//-----------------------------------Constructing cylinder----------------------------------
					double rx = dx*i-dx*0.5;
					double ry = dy*j-dy*0.5;
					if(sqrt(pow(center_x-rx,2) + pow(center_y-ry,2)) <= radius) {
						lb.fluid[i][j][k].type = TYPE_S;
						count++;
					}
					
					//-----------------------------------------------------------------------------------

					if(i == 0) lb.fluid[i][j][k].type = TYPE_E; //inlet outlet boundary
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //solid boundary
					if( Nz!=1 and (k == 0 or  k == Nz-1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
					if(lb.fluid[i][j][k].type == TYPE_F || lb.fluid[i][j][k].type == TYPE_E) { //Fluid domain
						lb.fluid[i][j][k].density = 1.0;
						lb.fluid[i][j][k].ux = u_max;
						lb.fluid[i][j][k].uy = 0.0;
						lb.fluid[i][j][k].uz = 0.0;
						count++;
					}
					
					
				}
			}
		}
		cout<<"Count = "<<count<<endl;
		
		return lb;
	}

#elif defined Couette_2D 

	LBM main_setup_Couette() {
		LBM lb(Nx, Ny, Nz, tau);
		
		double nu = cs2*(tau - 0.5);
		cout<<"nu = "<<nu<<", u_wall = "<<u_wall<<endl;

		#pragma omp parallel for 
		for (int i = 0 ; i<=Nx+1; i++) {
			for (int j = 0; j<=Ny+1; j++) {
				for (int k = 0; k<=Nz+1; k++) {
				

					if(i == 0 or i == Nx+1) lb.fluid[i][j][k].type = TYPE_F; //periodic boundary for inlet and outlet
					if(j == 0 or j == Ny+1) lb.fluid[i][j][k].type = TYPE_S; //solid boundary
					if( Nz!=1 and (k == 0 or  k == Nz-1) ) lb.fluid[i][j][k].type = TYPE_P; //periodic boundary
				}
			}
		}
		
		return lb;

	}

#endif

//---------------------------------------------------------------------------------------------------------------
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

void OutputVTK(int &nout, LBM &lb) {
	int		i,j,k;
	char	filename[128];
	FILE	*fp;
	unsigned int array_size;
	unsigned long int offset=0;
	// short  num16; // Int16 2byte
	float  val32; // Float32 4byte

	sprintf(filename,"./field%06d.vtr",nout);
	fp=fopen(filename,"wb");

	fprintf(fp,"<?xml version=\"1.0\"?>\n");
	fprintf(fp,"<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
	fprintf(fp,"  <RectilinearGrid WholeExtent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny+2,Nz);
	fprintf(fp,"  <Piece Extent=\"0 %d 0 %d 0 %d\">\n",Nx,Ny+2,Nz);
	fprintf(fp,"    <PointData>\n");
	fprintf(fp,"    </PointData>\n");
	fprintf(fp,"    <CellData>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Density\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny+2)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+3*4*(Nx)*(Ny+2)*(Nz);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CellType\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx)*(Ny+2)*(Nz);
	fprintf(fp,"    </CellData>\n");
	fprintf(fp,"    <Coordinates>\n");
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateX\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nx+1);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateY\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Ny+3);
	fprintf(fp,"      <DataArray type=\"Float32\" Name=\"CoordinateZ\" format=\"appended\" offset=\"%ld\" />\n",offset); offset+=4+1*4*(Nz+1);
	fprintf(fp,"    </Coordinates>\n");
	fprintf(fp,"  </Piece>\n");
	fprintf(fp,"  </RectilinearGrid>\n");
	fprintf(fp,"  <AppendedData encoding=\"raw\">");
	fprintf(fp,"_");

    // Density (cell)
	array_size=1*4*(Nx)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx;i++){
				val32=(float)lb.fluid[i][j][k].density; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}
    // Velocity (cell)
	array_size=3*4*(Nx)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx;i++){
				val32=(float)lb.fluid[i][j][k].ux; fwrite(&val32,sizeof(float),1,fp);
				val32=(float)lb.fluid[i][j][k].uy; fwrite(&val32,sizeof(float),1,fp);
				val32=0.0; fwrite(&val32,sizeof(float),1,fp);
			}
		}
	}

    // CellType (cell)
	array_size=1*4*(Nx)*(Ny+2)*(Nz);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz;k++){
		for(j=0;j<Ny+2;j++){
			for(i=0;i<Nx;i++){
				val32=(int)lb.fluid[i][j][k].type; fwrite(&val32,sizeof(int),1,fp);
			}
		}
	}
	// Coordinates (vertices)
	array_size=1*4*(Nx+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(i=0;i<Nx+1;i++){ val32=(float)(i*dx); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Ny+3);
	fwrite(&array_size,sizeof(int),1,fp);
	for(j=0;j<Ny+3;j++){ val32=(float)(j*dy); fwrite(&val32,sizeof(float),1,fp); }

	array_size=1*4*(Nz+1);
	fwrite(&array_size,sizeof(int),1,fp);
	for(k=0;k<Nz+1;k++){ val32=(float)(k*dz); fwrite(&val32,sizeof(float),1,fp); }

	fprintf(fp,"  </AppendedData>\n");
	fprintf(fp,"</VTKFile>\n");

	fclose(fp);
}
