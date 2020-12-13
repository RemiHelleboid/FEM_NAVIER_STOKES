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
#include "Geometry.hpp"

using namespace Eigen;
using namespace std;
typedef Matrix<double, 6, 1> Vector6d;

Vector3d Eval_Local_Shape_Func_Triangle_Ref_P1(double x, double y);
vector<Vector2d> Eval_Local_Shape_Grad_Func_Triangle_Ref_P1(const Triangle &T, double x, double y);
Vector6d Eval_Local_Shape_Func_Triangle_P2(double x, double y);
std::vector<Vector2d> Eval_Grad_Local_Shape_Func_Triangle_P2(const Triangle &T, double x, double y);
std::vector<Vector2d> Eval_Grad_Local_Shape_Func_Triangle_REF_P2(const Triangle &T, double lambda1, double lambda2, double lambda3);

double Eval_Quadrature_Stiff_P1(const Triangle &T, int i, int j);
double Eval_Quadrature_Stiff_P2(const Triangle &T, int i, int j);
double Eval_Quadrature_Mass_P2(const Triangle &T, int i, int j);
Matrix2d Compute_Jacobi_Matrix(const Triangle &T);
Matrix2d Compute_Inverse_Jacobi_Matrix(const Triangle &T);
