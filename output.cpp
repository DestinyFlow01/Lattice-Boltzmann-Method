//Definition of another function other than class member function 
#include<iostream>
#include<stdio.h>
#include "LBM.h"
#include "OUTPUT.h"
using namespace std;


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