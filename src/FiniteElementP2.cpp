#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FiniteElementP2.hpp"


using namespace Eigen;
using namespace std;
#define assertm(exp, msg) assert(((void)msg, exp))

FiniteElementP2::FiniteElementP2(Mesh & mesh0, VectorXd B0):mesh(mesh0), B(B0){
	mesh.build_half_bridges();
	const int d = mesh.get_nv();
	A = MatrixXd::Constant(d, d, 0);
	M = MatrixXd::Constant(d, d, 0);
	B = VectorXd::Constant(d, 1.0);
	U.resize(d);
}

Matrix2d FiniteElementP2::Compute_Jacobi_Matrix(const Triangle &T){
	vector<Vertice> T_vertices = T.get_vertices();		//3 sommets du triangle T
	double J11 = T_vertices[1].x() - T_vertices[0].x();		//xn2-xn1
	double J12 = T_vertices[2].x() - T_vertices[0].x();		//xn3-xn1
	double J21 = T_vertices[1].y() - T_vertices[0].y();		//yn2-yn1
	double J22 = T_vertices[2].y() - T_vertices[0].y();		//yn3-yn1
	Matrix2d J;
	J << J11, J12, J21, J22;
	return(J);
}

Matrix2d FiniteElementP2::Compute_Inverse_Jacobi_Matrix(const Triangle &T){
	Matrix2d J = Compute_Jacobi_Matrix(T);
	double area_K = T.area();
	Matrix2d J_inv;
	J_inv << J(1,1), -J(0,1),  -J(1,0), J(0,0);
	J_inv = (1.0/(2*area_K)) * J_inv;
	return(J_inv);
}

Vector2d FiniteElementP2::Compute_Grad_Triangle_Ref(int local_vertice_index){
	assertm(local_vertice_index<=3, "local vertice index over 3");
	Vector2d Grad;
	switch(local_vertice_index){
		case 0: Grad << -1, -1;
				break;
		case 1: Grad << 1, 0;
				break;
		case 2: Grad << 0, 1;
				break;
	}
	return(Grad);
}


Matrix<double, 6, 6> FiniteElementP2::Compute_Elementary_Matrix(const Triangle &T){
	Matrix<double, 6, 6> A_K = Matrix<double, 6, 6>::Zero();
	// double area_K = T.area();
	// double area_K0 = 0.5;
	double det_JK = Compute_Jacobi_Matrix(T).determinant();
	Matrix2d J_inv = Compute_Inverse_Jacobi_Matrix(T);
	Matrix2d J_prod = J_inv * J_inv.transpose();
	for(int sommet_p=0; sommet_p<3; sommet_p++){
		for(int sommet_q=0; sommet_q<3; sommet_q++){
			Vector2d Grad_p = 	Compute_Grad_Triangle_Ref(sommet_p);
			Vector2d Grad_q = Compute_Grad_Triangle_Ref(sommet_q);
			Vector2d V1 = J_prod * Grad_p;
			Vector2d V2 =  Grad_q;
			double alpha_k = det_JK * (1.0/3) * ( 2*int(sommet_p==sommet_q) - 0.5 );
			//cout<<sommet_p<<"   "<<sommet_q<<"   "<<alpha_k<<endl;
			//cout<<sommet_p<<"   "<<sommet_q<<"   "<<V1.dot(V2)<<endl;
		double A_pq_K = alpha_k * V1.dot(V2);
		A_K(sommet_p, sommet_q) += A_pq_K;
		}
	}

	vector<int> Vertice_4({0, 1});		//Meaning : Le sommet 4 est le milei du coté [0, 1]
	vector<int> Vertice_5({1, 2});
	vector<int> Vertice_6({0, 2});

	map<int, vector<int>> Midpoint_to_extremities;		//associe un milieu de segment avec les extremité du segment
	Midpoint_to_extremities[3] = Vertice_4;
	Midpoint_to_extremities[4] = Vertice_5;
	Midpoint_to_extremities[5] = Vertice_6;

	for(int sommet_i=0; sommet_i<3; sommet_i++){
		for(int sommet_k=3; sommet_k<6; sommet_k++){
			int sommet_p = Midpoint_to_extremities[sommet_k][0];
			int sommet_q = Midpoint_to_extremities[sommet_k][1];
			Vector2d Grad_i = 	Compute_Grad_Triangle_Ref(sommet_i);
			Vector2d Grad_p = 	Compute_Grad_Triangle_Ref(sommet_p);
			Vector2d Grad_q = 	Compute_Grad_Triangle_Ref(sommet_q);
			Vector2d V1 = J_prod * Grad_i;
			double alpha_k_1 = det_JK *( (2.0/3)*int(sommet_i==sommet_p) + 0.0);
			double alpha_k_2 = det_JK *( (2.0/3)*int(sommet_i==sommet_q) + 0.0);
			double A_pq_K = alpha_k_1 * V1.dot(Grad_q) + alpha_k_2 * V1.dot(Grad_p);
			A_K(sommet_i, sommet_k) += A_pq_K;
			A_K(sommet_k, sommet_i) += A_pq_K;
			cout<<int(sommet_i==sommet_p)<<endl;
		}
	}
	for(int sommet_p=3; sommet_p<6; sommet_p++){
		for(int sommet_q=3; sommet_q<6; sommet_q++){
			int sommet_i = Midpoint_to_extremities[sommet_p][0];
			int sommet_j = Midpoint_to_extremities[sommet_p][1];
			int sommet_k = Midpoint_to_extremities[sommet_q][0];
			int sommet_l = Midpoint_to_extremities[sommet_q][1];
			Vector2d Grad_i = 	Compute_Grad_Triangle_Ref(sommet_i);
			Vector2d Grad_j = 	Compute_Grad_Triangle_Ref(sommet_j);
			Vector2d Grad_k = 	Compute_Grad_Triangle_Ref(sommet_k);
			Vector2d Grad_l = 	Compute_Grad_Triangle_Ref(sommet_l);
			Vector2d V1 = J_prod * Grad_j;
			Vector2d V2 = J_prod * Grad_i;
			double alpha_k_1 = det_JK * (2.0/3)*( 1 + int(sommet_i==sommet_k) );
			double alpha_k_2 = det_JK * (2.0/3)*( 1 + int(sommet_i==sommet_l) );
			double alpha_k_3 = det_JK * (2.0/3)*( 1 + int(sommet_j==sommet_k) );
			double alpha_k_4 = det_JK * (2.0/3)*( 1 + int(sommet_j==sommet_l) );
			double A_pq_K = alpha_k_1 * V1.dot(Grad_l) + alpha_k_2 * V1.dot(Grad_k)
							+ alpha_k_3 * V2.dot(Grad_l) + alpha_k_4 * V2.dot(Grad_k);
			cout<<"ici<<  "<<A_pq_K<<endl;
			A_K(sommet_p, sommet_q) += A_pq_K;

		}
	}


	cout<<A_K<<endl<<endl;
	return(A_K);
}

bool FiniteElementP2::check_sizes(){
	int d = mesh.get_nv();
	bool size_ok;
	size_ok = ((A.rows()==d) && (A.rows()==d) && (B.size()==d) && (U.size()==d));
	return(size_ok);
}

void FiniteElementP2::Compute_Rigidity_Matrix(){
	//assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix<double, 6, 6> Elementary_matrix_K = Compute_Elementary_Matrix(T_k);		//calcul de la matrice élémentaire (3x3)

		for(int p=0; p<6; p++){
			for(int q=0; q<6; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k
				if(mesh.get_vertice(i).get_label()==1)  {A(i,j) = int(i==j); B(i) = 0;}
				else{A(i,j) = A(i,j) + Elementary_matrix_K(p,q);}		//On update le coefficients i, j de la matrice de rigidité

			}
		}
	}
}

inline double kron(int i, int k){
	double delta;
	int d = (i==k);
	delta = d;
	return(delta);
}

Matrix<double, 6, 6> FiniteElementP2::Compute_Elementary_Mass_Matrix(const Triangle &T){
	double det_JK = Compute_Jacobi_Matrix(T).determinant();
	Matrix<double, 6, 6> Elementary_mass_matrix_K = Matrix<double, 6, 6>::Zero();
	for(int sommet_p=0; sommet_p<3; sommet_p++){
		for(int sommet_q=0; sommet_q<3; sommet_q++){
			double M_p_q = 0.5 * det_JK * (1.0/60.0)*(13 - 11 * kron(sommet_p, sommet_q));
			Elementary_mass_matrix_K(sommet_p, sommet_q) += M_p_q;
		}
	}
	return(Elementary_mass_matrix_K);
}


void FiniteElementP2::Compute_Mass_Matrix(){
	//assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix<double, 6, 6> Elementary_mass_matrix_K = Compute_Elementary_Mass_Matrix(T_k);		//calcul de la matrice élémentaire (3x3)

		for(int p=0; p<6; p++){
			for(int q=0; q<6; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k
				if(mesh.get_vertice(i).get_label()==1)  {A(i,j) = int(i==j); B(i) = 0;}
				else{M(i,j) = M(i,j) + Elementary_mass_matrix_K(p,q);}		//On update le coefficients i, j de la matrice de rigidité

			}
		}
	}
	//cout<<A<<endl;
}


void FiniteElementP2::Direct_Method_Solve_Systeme(string Solver_type){
	map<string, int> map_solver = {{"LU", 1}, {"CHOLESKY",2}};
	assertm( 1, "Wrong Direct Method Solver Name" );
	int solver_index = map_solver[Solver_type];
	switch(solver_index){
		case 1:{ U = A.fullPivLu().solve(B);
				 break; }
		case 2 :{ LLT<MatrixXd> Decomposition(A);
				  Decomposition.solve(B);
				  break; }
				}
	A_is_modified = true;
}

double FiniteElementP2::Compute_l2_error(){
	double error_l2 = (A*U - B).norm();		//.norm() is an Eigen method that compute l2 norm of an array
	return(error_l2);
}

void FiniteElementP2::Display_Linear_Syst(){
	cout<<"Nombre de sommets : "<<mesh.get_nv()<<endl;
	cout<<"Nombre de triangle : "<<mesh.get_nt()<<endl;
	cout<<" Matrice de rigidité : "<<endl;
	cout<<A<<endl;
	cout<<"Second Membre : "<<endl;
	cout<<B<<endl;
	cout<<"Solution : "<<endl;
	cout<<U<<endl;
}

void FiniteElementP2::Export_Solution(string filename){
	ofstream F(filename);
	int Nb_edge = mesh.get_nae();
	for(int k=0; k<Nb_edge; k++){
		Edge E_k = mesh.get_edge(k);
		vector<Vertice> Vk = E_k.get_vertices();
		for(int p=0; p<2; p++){
			int idx = Vk[p].get_index();
			F << Vk[p].x() <<"\t"<< Vk[p].y()<<"\t" <<U[idx]<<endl;
		}
	F<<endl;
	}
}

void FiniteElementP2::Compute_second_member(){


	}
