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
	mesh.show_parameters();
	const int d = mesh.get_nv();
	int Number_vertices = mesh.get_nv();
	A = SparseMatrix<double>(Number_vertices, Number_vertices);
	A.reserve(VectorXi::Constant(Number_vertices, 12));
	// std::cout<<A<<std::endl;
	M = MatrixXd::Constant(d, d, 0);
	B = VectorXd::Constant(d, 0);
	U = VectorXd::Constant(d, 0);

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


Matrix<double, 6, 6> FiniteElementP2::Compute_Elementary_Stiffness_Matrix(const Triangle &T){
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
			A_K(sommet_p, sommet_q) += A_pq_K;
		}
	}
	return(A_K);

}

Matrix<double, 6, 6> FiniteElementP2::Compute_Elementary_Stiffness_Matrix_test(const Triangle &T){
	Matrix<double, 6, 6> A_K = Matrix<double, 6, 6>::Zero();
	// double area_K = T.area();
	// double area_K0 = 0.5;
	double det_JK = Compute_Jacobi_Matrix(T).determinant();
	Matrix2d J_inv = Compute_Inverse_Jacobi_Matrix(T);
	Matrix2d J_prod = J_inv * J_inv.transpose();
	for(int sommet_p=0; sommet_p<6; sommet_p++){
		for(int sommet_q=0; sommet_q<6; sommet_q++){
			A_K(sommet_p, sommet_q) = Eval_Quadrature_Stiff_P2(T, sommet_q, sommet_p);
		}
	}
	return(A_K);
}

bool FiniteElementP2::check_sizes(){
	int d = mesh.get_nv();
	bool size_ok;
	size_ok = ((A.rows()==d) && (A.rows()==d) && (B.size()==d) && (U.size()==d));
	return(size_ok);
}

void FiniteElementP2::Compute_Stiffness_Matrix(){
	//assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix<double, 6, 6> Elementary_matrix_K = Compute_Elementary_Stiffness_Matrix(T_k);		//calcul de la matrice élémentaire (3x3)
		Matrix<double, 6, 6> Elementary_matrix_K_test = Compute_Elementary_Stiffness_Matrix_test(T_k);		//calcul de la matrice élémentaire (3x3)

		std::cout<<Elementary_matrix_K<<std::endl<<std::endl;
		std::cout<<Elementary_matrix_K_test<<std::endl<<std::endl<<std::endl;
		for(int p=0; p<6; p++){
			for(int q=0; q<6; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k
				if(mesh.get_vertice(i).get_label()!=0)  {A.coeffRef(i, j) = int(i==j); B(i) = 0;}
				else{A.coeffRef(i, j) += Elementary_matrix_K_test(q,p);}		//On update le coefficients i, j de la matrice de rigidité
			}
		}
	}
	A.makeCompressed();
}


inline double kron(int i, int k){
	double delta;
	int d = (i==k);
	delta = d;
	return(delta);
}

//determine l'intégrande entre l⁴, l²*l², l³*l, l*l²*l
inline int Compute_Integral_type(int i, int j, int k ,int l){
	int integrande_type = -1;
	if(i==j && j==k && k==l){
		integrande_type = 1;
	}
	else if( (i==j && j==k && k!=l) || (i==j && j==l && l!=k)
			||(i==k && k==l && l!=j) ||(j==k && k==l && l!=i) ){
		integrande_type = 2;

	}
	else if((i==j && k==l && j!=k) || (i==k && j==l && l!=k) || (i==l && k==j && l!=k)){
		integrande_type = 3;

	}
	else if((i==j && i!=l && i!=k) || (i==k && i!=l && i!=j) || (i==l && i!=j && i!=k)
			|| (j==k && i!=j && k!=l)  || (j==l && j!=k && i!=l)){
		integrande_type = 4;
	}
	if(integrande_type==-1){cout<<"Error no integrande type match for (i, j, k, l) = "<<i<<" "<<j<<" "<<k<<" "<<l<<endl;}
	return(integrande_type);
}


//Calcul de l'intégral des coordonée bary_centriques sur le triangle de ref  pour p, q des milieux d'arretes
inline double Compute_Integral(int sommet_p, int sommet_q){
	double coeff = 0;

	vector<int> Vertice_4({0, 1});		//Meaning : Le sommet 4 est le milei du coté [0, 1]
	vector<int> Vertice_5({1, 2});
	vector<int> Vertice_6({0, 2});

	map<int, vector<int>> Midpoint_to_extremities;		//associe un milieu de segment avec les extremité du segment
	Midpoint_to_extremities[3] = Vertice_4;
	Midpoint_to_extremities[4] = Vertice_5;
	Midpoint_to_extremities[5] = Vertice_6;
	int sommet_i = Midpoint_to_extremities[sommet_p][0];
	int sommet_j = Midpoint_to_extremities[sommet_p][1];
	int sommet_k = Midpoint_to_extremities[sommet_q][0];
	int sommet_l = Midpoint_to_extremities[sommet_q][1];

	int type_integrande = Compute_Integral_type(sommet_i, sommet_j, sommet_k, sommet_l);

	if (type_integrande==3){coeff = 4.0/45;};
	if (type_integrande==4){coeff = 2.0/45;};

	return(coeff);
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
	for(int sommet_p=3; sommet_p<6; sommet_p++){
		for(int sommet_q=3; sommet_q<6; sommet_q++){
			double M_p_q = det_JK * Compute_Integral(sommet_p, sommet_q);
			Elementary_mass_matrix_K(sommet_p, sommet_q) += M_p_q;
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
		for(int sommet_p=3; sommet_p<6; sommet_p++){
			double M_p_q = 0;
			int sommet_l = Midpoint_to_extremities[sommet_p][0];
			int sommet_k = Midpoint_to_extremities[sommet_p][1];
			if (sommet_i == sommet_l || sommet_i==sommet_k){
				M_p_q = det_JK * 1 ;
			}
			else if(sommet_i!=sommet_l && sommet_i!=sommet_k){
				M_p_q = det_JK * (29.0/90);
			}
			else{cout<<"ERROR CASE IN MATRIX ELEMENTARY MASS"<<endl;}
			Elementary_mass_matrix_K(sommet_i, sommet_p) = M_p_q;
			Elementary_mass_matrix_K(sommet_p, sommet_i) = M_p_q;
		}
	}
	return(Elementary_mass_matrix_K);
}

Matrix<double, 6, 6> FiniteElementP2::Compute_Elementary_Mass_Matrix_test(const Triangle &T){
	Matrix<double, 6, 6> Elementary_mass_matrix_K(6,6);
	Elementary_mass_matrix_K << 1./60, -1./360, -1./360, 0, -1./90, 0,
	 							-1./360, 1./60, -1./360, 0, 0, -1./90,
	  							-1./360, -1./360, 1./60, -1./90, 0, 0,
	   							0, 0, -1./90, 4./45, 2./45, 2./45,
								-1./90, 0, 0, 2./45, 4./45, 2./45,
								0, -1./90, 0, 2./45, 2./45, 4./45;
	Elementary_mass_matrix_K = T.area() * Elementary_mass_matrix_K;
	return(Elementary_mass_matrix_K);
}


void FiniteElementP2::Compute_Mass_Matrix(){
	//assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix<double, 6, 6> Elementary_mass_matrix_K = Compute_Elementary_Mass_Matrix_test(T_k);		//calcul de la matrice élémentaire (3x3)

		for(int p=0; p<6; p++){
			for(int q=0; q<6; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k

				if(mesh.get_vertice(i).get_label()!=0)  {M(i,j) = int(i==j); B(i) = 0;}
				else{M(i,j) = M(i,j) + Elementary_mass_matrix_K(p,q);}		//On update le coefficients i, j de la matrice de rigidité
			}
		}
	}
}

Vector6d FiniteElementP2::Compute_Elementary_Second_Member(const Triangle &T){
	Matrix<double, 6, 6> Mass_Matrix_K = Compute_Elementary_Mass_Matrix_test(T);
	Vector6d F;
	F << 1, 1, 1, 1, 1, 1;
	Vector6d B_K = Mass_Matrix_K * F;
	return(B_K);
}


// Vector6d  FiniteElementP2::Compute_Elementary_Second_Member_Quadrature(const Triangle &T){
// 	double area_K = T.area();
// 	Vector6d F0, F;
// 	F0 << 1, 1, 1, 1, 1, 1;
// 	for(int i=0; i<6; i++){
// 		F(i) = Eval_Quadrature_Mass(T, i);
// 	}
// 	return(F);
// }


void FiniteElementP2::Compute_Second_Member(){
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){
		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Vector6d Elementary_Vector_K = Compute_Elementary_Second_Member(T_k);		//calcul de la matrice élémentaire (3x3)
		for(int sommet_p=0; sommet_p<6; sommet_p++){
			int i = global_index_vertices[sommet_p];
			B(i) += Elementary_Vector_K(sommet_p);	//On update le coefficients i
		}
	}

}


void FiniteElementP2::Direct_Method_Solve_Systeme(string Solver_type){
	map<string, int> map_solver = {{"LU", 1}, {"CHOLESKY",2}, {"ITERATIVE", 3}};
	assertm( 1, "Wrong Direct Method Solver Name" );
	int solver_index = map_solver[Solver_type];
	switch(solver_index){
		case 1:{ SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
				 solver.analyzePattern(A);
				 solver.factorize(A);
				 U = solver.solve(B);
				 break; }
		case 2 :{ LLT<MatrixXd> Decomposition(A);
				  Decomposition.solve(B);
				  break; }
		case 3 :{ BiCGSTAB<SparseMatrix<double> > solver;
		  				  solver.compute(A);
						  U = solver.solve(B);
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
	int Nb_vertice = mesh.get_nv();
	F<<'X'<<','<<'Y'<<','<<'U'<<endl;
	for(int k=0; k<Nb_vertice; k++){
		Vertice V_k = mesh.get_vertice(k);
			int idx = V_k.get_index();
			F << V_k.x() <<","<< V_k.y()<<"," <<U[idx]<<endl;
		}
	}
