#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <Eigen/Dense>
#include "Mesh.hpp"
#include "Geometry.hpp"
#include "FiniteElementP1.hpp"
#include "FiniteElementP2.hpp"
#include "Stokes.hpp"
#include "Quadratures.hpp"

using namespace std;

double f2(double x, double y){
	return(1.0);
}

VectorXd Full_Ones_Vect(int size){
	VectorXd V(size);
	for(int k=0; k<size; k++){
		V(k) = 0.0;
	}
	return(V);
}



int  main(int argc, const char** argv){
	if (strcmp( argv[1], "P1") == 0){
		Mesh M0("FreeFem/Square.msh");
		cout<<"Mesh chargée : Done"<<endl;
		M0.show_parameters();
		VectorXd Second_member = Full_Ones_Vect(M0.get_nv());
		FiniteElementP1 FEP1(M0, Second_member);
		FEP1.Compute_Second_Member();
		FEP1.Compute_Stiffness_Matrix();
		FEP1.Direct_Method_Solve_Systeme("LU");
		cout<<"error:"<<FEP1.Compute_l2_error()<<endl;
		FEP1.Export_Solution("SolutionP1.csv");
		FEP1.Display_Linear_Syst();
		FEP1.Display_Solution();
	}
	else if (strcmp( argv[1], "P2") == 0){
		Mesh M0("FreeFem/Square.msh");
		cout<<"Mesh chargée : Done"<<endl;
		M0.show_parameters();
		VectorXd Second_member = Full_Ones_Vect(M0.get_nv());
		FiniteElementP2 FEP1(M0, Second_member);
		FEP1.Compute_Second_Member();
		FEP1.Compute_Stiffness_Matrix();
		FEP1.Direct_Method_Solve_Systeme("LU");
		cout<<"error:"<<FEP1.Compute_l2_error()<<endl;
		FEP1.Export_Solution("SolutionP2.csv");
		FEP1.Display_Linear_Syst();
		FEP1.Display_Solution();
	}

	return(0);
}


// int main1(){
// 	Mesh M0("FreeFem/Circle.msh");
// 	//M0.build_half_bridges();
// 	M0.show_parameters();
// 	VectorXd Second_member = Full_Ones_Vect(M0.get_nv());
// 	FiniteElementP2 FEP2(M0, Second_member);
// 	FEP2.Compute_Second_Member();
// 	FEP2.Compute_Rigidity_Matrix();
// 	FEP2.Direct_Method_Solve_Systeme("ITERATIVE");
// 	FEP2.Export_Solution("SolutionP2.data");
// 	//FEP2.Display_Solution();
// 	cout <<FEP2.Compute_l2_error()<<endl;
//
// 	return(0);
// }
