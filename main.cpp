#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <Eigen/Dense>
#include "Mesh.hpp"
#include "Geometry.hpp"
#include "FiniteElementP1.hpp"
#include "FiniteElementP2.hpp"

using namespace std;

double f2(double x, double y){
	return(pow(3-x,2)+4*pow(0.5-y,2)-x*y);
}

VectorXd Full_Ones_Vect(int size){
	VectorXd V(size);
	for(int k=0; k<size; k++){
		V(k) = 0.0;
	}
	return(V);
}



int main1(){
	Mesh M0("FreeFem/Th.msh");
	M0.show_parameters();
	VectorXd Second_member = Full_Ones_Vect(M0.get_nv());
	FiniteElementP1 FEP1(M0, Second_member);
	FEP1.Compute_Second_Member();
	FEP1.Compute_Rigidity_Matrix();
	FEP1.Direct_Method_Solve_Systeme("LU");
	cout<<"error:"<<FEP1.Compute_l2_error()<<endl;
	FEP1.Export_Solution("Solution.data");
	//FEP1.Display_Linear_Syst();
	FEP1.Display_Solution();

	return(0);
}


int main(){
	Mesh M0("FreeFem/Th.msh");
	//M0.build_half_bridges();
	M0.show_parameters();
	VectorXd Second_member = Full_Ones_Vect(M0.get_nv());
	FiniteElementP2 FEP2(M0, Second_member);
	FEP2.Compute_Rigidity_Matrix();
	FEP2.Direct_Method_Solve_Systeme("LU");
	FEP2.Export_Solution("SolutionP22.data");



	return(0);
}
