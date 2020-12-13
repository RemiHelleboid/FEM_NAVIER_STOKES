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
#include "Quadra.hpp"

using namespace std;


double f(double x, double y){
    return(x+y);
}

int main(){
    Mesh M0("FreeFem/Square.msh");
	cout<<"Mesh chargée : Done"<<endl;
	M0.show_parameters();
    Triangle T1 = M0.get_triangle(0);
    cout<<"Area of triangle "<<T1.area()<<endl;
    double Int = Quad_Triangle_Ref(T1, f);
    cout<<"Intégrale = "<<Int<<endl;
}
