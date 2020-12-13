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
#include <Eigen/Sparse>

#include "FiniteElementP1.hpp"


using namespace Eigen;
using namespace std;


#define assertm(exp, msg) assert(((void)msg, exp))

double func_force(double x, double y){
	double norm = sqrt(x*x + y*y);
	double f = 0;
	if (norm < 0.15)
		f = -100;
	else
		f = 0;
	return(10);
}


FiniteElementP1::FiniteElementP1(const Mesh &mesh0, VectorXd B0):mesh(mesh0){
	int Number_vertices = mesh.get_nv();
	A = SparseMatrix<double>(Number_vertices, Number_vertices);
	A.reserve(VectorXi::Constant(Number_vertices, 12));
	M.resize(Number_vertices, Number_vertices);
	B.resize(Number_vertices);
	U.resize(Number_vertices);
}

Matrix2d FiniteElementP1::Compute_Jacobi_Matrix(const Triangle &T){
	vector<Vertice> T_vertices = T.get_vertices();		//3 sommets du triangle T
	double J11 = T_vertices[1].x() - T_vertices[0].x();		//xn2-xn1
	double J12 = T_vertices[2].x() - T_vertices[0].x();		//xn3-xn1
	double J21 = T_vertices[1].y() - T_vertices[0].y();		//yn2-yn1
	double J22 = T_vertices[2].y() - T_vertices[0].y();		//yn3-yn1
	Matrix2d J;
	J << J11, J12, J21, J22;
	return(J);
}

Matrix2d FiniteElementP1::Compute_Inverse_Jacobi_Matrix(const Triangle &T){
	Matrix2d J = Compute_Jacobi_Matrix(T);
	double area_K = T.area();
	Matrix2d J_inv;
	J_inv << J(1,1), -J(0,1),  -J(1,0), J(0,0);
	J_inv = (1.0/(2*area_K)) * J_inv;
	return(J_inv);
}

Vector2d FiniteElementP1::Compute_Grad_Triangle_Ref(int local_vertice_index){
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
	return(Grad.transpose());
}

bool FiniteElementP1::check_sizes(){
	int Number_vertices = mesh.get_nv();
	bool size_ok;
	size_ok = ((A.rows()==Number_vertices) && (A.rows()==Number_vertices) && (B.size()==Number_vertices) && (U.size()==Number_vertices));
	return(size_ok);
}

Matrix3d FiniteElementP1::Compute_Elementary_Stiffness_Matrix(const Triangle &T){
	Matrix3d A_K = Matrix3d::Constant(3,3,0);
	double area_K = T.area();
	vector<Vertice> T_vertices = T.get_vertices();		//3 sommets du triangle T
	double J11 = T_vertices[1].x() - T_vertices[0].x();		//xn2-xn1
	double J12 = T_vertices[2].x() - T_vertices[0].x();		//xn3-xn1
	double J21 = T_vertices[1].y() - T_vertices[0].y();		//yn2-yn1
	double J22 = T_vertices[2].y() - T_vertices[0].y();		//yn3-yn1
	double A22 = (1.0/(4*area_K)) * ( pow((J12),2) + pow((J22),2) );
	double A33 = (1.0/(4*area_K)) * ( pow((J11),2) + pow((J21),2) );
	double A23 = (1.0/(4*area_K)) * ( J12 * (-J11) + (J22) * (-J21) );

	A_K << A22+A33+2*A23, -A23-A22, -A23-A33,
				 -A23-A22,       A22,      A23,
				 -A23-A33,       A23, 		 A33;

	return(A_K);
}

Matrix3d FiniteElementP1::Compute_Elementary_Stiffness_Matrix_test(const Triangle &T){
	Matrix3d A_K = Matrix3d::Constant(3,3,0);
	double area_K = T.area();
	vector<Vertice> T_vertices = T.get_vertices();		//3 sommets du triangle T
	for(int p=0; p<3; p++){
		for(int q=0; q<3; q++){
			A_K(p,q) = Eval_Quadrature_Stiff_P1(T, p, q);
		}
	}
	return(A_K);
}


Matrix3d FiniteElementP1::Compute_Elementary_Mass_Matrix(const Triangle &T){
	Matrix3d Mass_Matrix(3,3);
	Mass_Matrix << 2, 1, 1,
	 			   1, 2, 1,
				   1, 1, 2;
	double a_K = T.area() / 12.0;
	Mass_Matrix = a_K * Mass_Matrix;
	return(Mass_Matrix);
}


void FiniteElementP1::Compute_Stiffness_Matrix(){
	//assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix3d Elementary_matrix_K = Compute_Elementary_Stiffness_Matrix(T_k);		//calcul de la matrice élémentaire (3x3)
		Matrix3d Elementary_matrix_K_test = Compute_Elementary_Stiffness_Matrix_test(T_k);		//calcul de la matrice élémentaire (3x3)
		std::cout<<Elementary_matrix_K<<std::endl<<std::endl;
		std::cout<<Elementary_matrix_K_test<<std::endl<<std::endl<<std::endl;

		for(int p=0; p<3; p++){
			for(int q=0; q<3; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k
				// std::cout<<p<<"\t"<<i<<"\t"<<q<<"\t"<<j<<"\t"<<std::endl;

				if(mesh.get_vertice(i).get_label()!=0)  {A.coeffRef(i, j) = int(i==j); B(i) = 0;}
				else{A.coeffRef(i, j) += Elementary_matrix_K_test(p,q);}		//On update le coefficients i, j de la matrice de rigidité
			}
		}
	}
	A.makeCompressed();
}

void FiniteElementP1::Compute_Mass_Matrix(){
	assertm(check_sizes(), "Sizes error :sizes doesn't match");
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Matrix3d Elementary_matrix_K = Compute_Elementary_Mass_Matrix(T_k);		//calcul de la matrice élémentaire (3x3)

		for(int p=0; p<3; p++){
			for(int q=0; q<3; q++){
				int i = global_index_vertices[p];		//On récupère les indices globaux des sommets du triangle T_k
				int j = global_index_vertices[q];		//On récupère les indices globaux des sommets du triangle T_k
				M(i,j) = M(i,j) + Elementary_matrix_K(p,q);	//On update le coefficients i, j de la matrice de rigidité

			}
		}
	}
}

Vector3d FiniteElementP1::Compute_Elementary_Second_Member(const Triangle &T){
	Vector3d F;
	Matrix3d Elementary_mass_matrix_K = Compute_Elementary_Mass_Matrix(T);		//calcul de la matrice élémentaire (3x3)
	vector<Vertice> T_vertices = T.get_vertices();		//3 sommets du triangle T
	double x1 = T_vertices[0].x();
	double x2 = T_vertices[1].x();
	double x3 = T_vertices[2].x();
	double y1 = T_vertices[0].y();
	double y2 = T_vertices[1].y();
	double y3 = T_vertices[2].y();
	double f_vertex_1 = func_force(x1, y1);
	double f_vertex_2 = func_force(x2, y2);
	double f_vertex_3 = func_force(x3, y3);
	F << f_vertex_1, f_vertex_2, f_vertex_3;
	F = 2*(T.area()/3.0) * F;
	return(F);
}


void FiniteElementP1::Compute_Second_Member(){
	std::cout<<"COMPUTE SECONDE MEMBER"<<std::endl;
	unsigned int N_T = mesh.get_nt();
	for(unsigned int k = 0; k<N_T; k++){

		Triangle T_k = mesh.get_triangle(k);
		vector<int> global_index_vertices = T_k.get_vertices_index();
		Vector3d Elementary_Vector_K = Compute_Elementary_Second_Member(T_k);		//calcul de la matrice élémentaire (3x3)
		std::cout<<Elementary_Vector_K<<std::endl;
		for(int sommet_p=0; sommet_p<3; sommet_p++){
				int i = global_index_vertices[sommet_p];		//On récupère les indices globaux des sommets du triangle T_k
				B(i) = B(i) + Elementary_Vector_K(sommet_p);	//On update le coefficients i, j de la matrice de rigidité
		}
	}
}

void FiniteElementP1::Direct_Method_Solve_Systeme(string Solver_type){
	map<string, int> map_solver = {{"LU", 1}, {"CHOLESKY",2}, {"ITERATIVE", 3}};
	assertm( 1, "Wrong Direct Method Solver Name" );
	int solver_index = map_solver[Solver_type];
	switch(solver_index){
		case 1:{ SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
				 solver.analyzePattern(this->A);
				 solver.factorize(this->A);
				 U = solver.solve(this->B);
				 break; }
		case 2 :{ LLT<MatrixXd> Decomposition(A);
				  Decomposition.solve(this->B);
				  break; }
		case 3 :{ BiCGSTAB<SparseMatrix<double> > solver;
		  				  solver.compute(this->A);
						  U = solver.solve(this->B);
		  				  break; }
				}
	A_is_modified = true;
}

double FiniteElementP1::Compute_l2_error(){
	double error_l2 = (A*U - B).norm();		//.norm() is an Eigen method that compute l2 norm of an array
	return(error_l2);
}

void FiniteElementP1::Display_Linear_Syst(){
	cout<<"Nombre de sommets : "<<mesh.get_nv()<<endl;
	cout<<"Nombre de triangle : "<<mesh.get_nt()<<endl;
	cout<<" Matrice de rigidité : "<<endl;
	cout<<A<<endl;
	cout<<"Second Membre : "<<endl;
	cout<<B<<endl;
	cout<<"Solution : "<<endl;
	cout<<U<<endl;
}

void FiniteElementP1::Export_Solution(string filename){
	ofstream F(filename);
	int Nb_edge = mesh.get_nae();
	F<<'X'<<','<<'Y'<<','<<'U'<<endl;
	for(int k=0; k<Nb_edge; k++){
		Edge E_k = mesh.get_edge(k);
		vector<Vertice> Vk = E_k.get_vertices();
		for(int p=0; p<2; p++){
			int idx = Vk[p].get_index();
			F << Vk[p].x() <<","<< Vk[p].y()<<"," <<U[idx]<<endl;
		}
	F<<endl;
	}
}
