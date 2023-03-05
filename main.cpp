#include<iostream>
//#include "OUTPUT.h"
#include "LBM.h"
using namespace std;



int main() {
	print_Logo();
	
	//Defining variable for simulation : 
	LBM lb = main_setup_Cylinder();
	//LBM lb = main_setup_Cylinder();
	lb.Init();
	cout<<lb.fluid[Nx/2][40][0].f[2]<<endl;
	//Preliminary Result : 
	int step = 0, n = 0;
	OutputVTK(step,lb);
	
	
	//Iterations 
	
	for(step = 1; step<=T_OUT; step++) {
		//Collision
		lb.Collision();
		
		//Streaming and BC
		lb.Streaming();
		lb.MacroProp();
		
		//Outputing file : 
		if(step%Nt == 0) {
			cout<<lb.fluid[Nx/2][Ny/2][0].ux<<endl;
			n++;
			OutputVTK(n,lb);
			cout<<"Step = "<<step<<endl;
		}
	}
	
}