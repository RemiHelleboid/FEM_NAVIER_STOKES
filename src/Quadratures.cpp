#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include <Eigen/Dense>
#include "Quadratures.hpp"

using namespace Eigen;
using namespace std;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 3, 1> Vector3d;


Vector3d Eval_Local_Shape_Func_Triangle_Ref_P1(double x, double y){
    const double lambda1 = 1 - x - y;
    const double lambda2 = x;
    const double lambda3 = y;
    Vector3d Value(3);
    Value[0] = lambda1;
    Value[1] = lambda2;
    Value[2] = lambda3;
    return(Value);
}
vector<Vector2d> Eval_Local_Shape_Grad_Func_Triangle_Ref_P1(const Triangle &T, double x, double y){
    vector<Vector2d> Value(3);
    Value[0] = {-1, -1};
    Value[1] = {1, 0};
    Value[2] = {0, 1};
    return(Value);
}

double Eval_Quadrature_Stiff_P1(const Triangle &T, int i, int j){
    double sum = 0;
    vector<Vector2d> V_grad_xi_1 = Eval_Local_Shape_Grad_Func_Triangle_Ref_P1(T, 0.5, 0.5);
    vector<Vector2d> V_grad_xi_2 = Eval_Local_Shape_Grad_Func_Triangle_Ref_P1(T, 0.5, 0.0);
    vector<Vector2d> V_grad_xi_3 = Eval_Local_Shape_Grad_Func_Triangle_Ref_P1(T, 0, 0.5);
    sum = T.GradLambda(i).dot(T.GradLambda(j));
    double integral_K = sum * T.get_Jacobian_det() * (1.0/8);
    return(integral_K);
}




//##################################################################################################




Vector6d Eval_Local_Shape_Func_Triangle_Ref_P2(double x, double y){
    const double lambda1 = 1 - x - y;
    const double lambda2 = x;
    const double lambda3 = y;
    Vector6d Value(6);

    Value[0] = lambda1 * (2 * lambda1 - 1);
    Value[1] = lambda2 * (2 * lambda2 - 1);
    Value[2] = lambda3 * (2 * lambda3 - 1);

    Value[3] = 4 * lambda1 * lambda2;
    Value[4] = 4 * lambda2 * lambda3;
    Value[5] = 4 * lambda3 * lambda1;

    return(Value);
}
std::vector<Vector2d> Eval_Grad_Local_Shape_Func_Triangle_P2(const Triangle &T, double x, double y){
    Matrix2d J_inv = T.get_Jacobian_inv_mat();
    const double lambda1 =  T.computeBarryCoord(x ,y)[0];
    const double lambda2 =  T.computeBarryCoord(x ,y)[1];
    const double lambda3 =  T.computeBarryCoord(x ,y)[2];

    std::vector<Vector2d> Value(6);
    Value[0] =( lambda1 * (2 * T.GradLambda(0)) + (2 * lambda1 - 1) * T.GradLambda(0) );
    Value[1] = ( lambda2 * (2 * T.GradLambda(1)) + (2 * lambda2 - 1) * T.GradLambda(1) );
    Value[2] = ( lambda3 * (2 * T.GradLambda(2)) + (2 * lambda3 - 1) * T.GradLambda(2) );

    Value[3] = 4 * T.GradLambda(0) * lambda2 + 4 * lambda1 * T.GradLambda(1);
    Value[4] = 4 * T.GradLambda(1) * lambda3 + 4 * lambda2 * T.GradLambda(2);
    Value[5] = 4 * T.GradLambda(0) * lambda3 + 4 * lambda1 * T.GradLambda(2);

    return(Value);
}

std::vector<Vector2d> Eval_Grad_Local_Shape_Func_Triangle_REF_P2(const Triangle &T, double lambda1, double lambda2){

    double lambda3 = 1 - lambda1 - lambda2;
    // std::cout<<lambda1<<"   "<<lambda2<<"   "<<lambda3<<std::endl;
    std::vector<Vector2d> Value(6);

    Vector2d GRAD0 = {-1, -1};
    Vector2d GRAD1 = {1, 0};
    Vector2d GRAD2 = {0, 1};

    Value[0] = lambda1 * (2 * GRAD0) + (2 * lambda1 - 1) * GRAD0;
    Value[1] = lambda2 * (2 * GRAD1) + (2 * lambda2 - 1) * GRAD1;
    Value[2] = lambda3 * (2 * GRAD2) + (2 * lambda3 - 1) * GRAD2;

    Value[3] = 4 * GRAD0 * lambda2 + 4 * lambda1 * GRAD1;
    Value[4] = 4 * GRAD1 * lambda3 + 4 * lambda2 * GRAD2;
    Value[5] = 4 * GRAD0 * lambda3 + 4 * lambda1 * GRAD2;

    std::cout<<"VALEUR "<<Value[1]<<std::endl;

    return(Value);
}

double Eval_Quadrature_Stiff_P2(const Triangle &T, int i, int j){
    vector<Vertice> V_half = T.get_vertices_middle();
    Vertice M12 = V_half[0];
    Vertice M23 = V_half[1];
    Vertice M31 = V_half[2];
    double sum = 0;
    double l1 = 1.0/3;
    double l2 = 1.0/3;
    sum = T.GradLambdaP2(i, l1, l2).dot(T.GradLambdaP2(j, 1.0/3, 1.0/3)) + T.GradLambdaP2(i, 1.0/3, 2.0/3).dot(T.GradLambdaP2(j, 1.0/3, 2.0/3)) + T.GradLambdaP2(i, l1, l2).dot(T.GradLambdaP2(j, 2.0/3, 2.0/3));
    double integral_K = 0.33 * sum * T.get_Jacobian_det() * (1.0/8);
    return(integral_K);
}

double Eval_Quadrature_Mass_P2(Triangle &T, int i, int j){
    vector<Vertice> V_half = T.get_vertices_middle();
    Vertice M12 = V_half[0];
    Vertice M23 = V_half[1];
    Vertice M31 = V_half[2];
    double sum = 0;
    double l1 = 1.0/3;
    double l2 = 1.0/3;
    sum = T.GradLambdaP2(i, l1, l2).dot(T.GradLambdaP2(j, l1, l2));
    double integral_K = sum * T.get_Jacobian_det() * (1.0/8);
    return(integral_K);
}
