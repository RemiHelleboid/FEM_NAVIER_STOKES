#pragma once
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include "Mesh.hpp"

using namespace Eigen;
using namespace std;

class FiniteElementP2{
	private:
		Mesh mesh;	//Maillage
		MatrixXd A;			//Matrice de rigidité
        MatrixXd M;			//Matrice de rigidité
		VectorXd B;			//second membre
		VectorXd U;			//Vecteur solution
		bool A_is_modified = false;
    public:
    //Constructeurs
		FiniteElementP2();
		FiniteElementP2(Mesh &mesh0, VectorXd B0);
		FiniteElementP2(const FiniteElementP2 &FE):mesh(FE.mesh), A(FE.A),M(FE.M), B(FE.B), U(FE.U), A_is_modified(FE.A_is_modified){};
		FiniteElementP2 & operator=(const FiniteElementP2 &) = default;

	//Méthodes Transformation triangles et interpollation
		Matrix2d Compute_Jacobi_Matrix(const Triangle &T);		//Compute la matrice J de la transformation de T vers l'éléments de référence
		Matrix2d Compute_Inverse_Jacobi_Matrix(const Triangle &T);		//Compute la matrice inverse J^-1 de la Jacobienne
		Vector2d Compute_Grad_Triangle_Ref(int i);		//Gradient de la fonction de base associé au sommet i du triangle de ref

	//Méthodes : Création et résolution du système linéaire
		Matrix<double, 6, 6> Compute_Elementary_Matrix(const Triangle &T);
		Matrix<double, 6, 6> Compute_Elementary_Mass_Matrix(const Triangle &T);
		Matrix3d Compute_Elementary_Matrix_Formal(Triangle T);
		void Compute_Rigidity_Matrix();
		void Compute_Mass_Matrix();
		void Compute_Dirichlet_Bound_Condition();	//To do or not ?
		void Compute_second_member();		//If needed
		void Direct_Method_Solve_Systeme(string Solver_type);
		double Compute_l2_error();
		bool check_sizes();
		void Display_Linear_Syst();
		void Display_Solution(){cout<<"Solution : \n"<<U<<endl;cout<<"MAX : \n"<<U.maxCoeff()<<endl;};
		void Export_Solution(string filename);
	};
