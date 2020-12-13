#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>
#include "Geometry.hpp"

using namespace std;

Edge::Edge(int l, vector<Vertice> V0):label(l){
	V.resize(2);
	int i = V0[0].get_index();				//On veut toujours avoir les indices croissants dans les edges pour éviter le doublons
	int j = V0[1].get_index();
	if(i<=j || 1==1){
		V = V0;
	}
	else{
		V[0] = V0[1];
		V[1] = V0[0];
	}
}

bool Edge::operator==(const Edge & E0)const{
	bool b;
	b = (E0.get_vertices()==get_vertices());
	return(b);
}

vector<double> Edge::vect(){
	vector<double> vect(2);
	vector<double> C0, C1;
	C0 = V[0].get_coord();
	C1 = V[1].get_coord();
	vect[0] = C1[0] - C0[0];
	vect[1] = C1[1] - C0[1];
	return(vect);
}

vector<double> Edge::perp(){
	vector<double> vec = vect();
	vector<double> perp(2);
	perp[0] = -vec[1];
	perp[1] = vec[0];
	return(perp);
}

Vertice Edge::Compute_middle(){
	double x_half =  0.5 * ( V[0].x() + V[1].x() );
	double y_half =  0.5 * ( V[0].y() + V[1].y() );
	int index_half = -1;		//Il faudra changer l'index si on ajoute le sommet a la liste des sommets du maillage
	int label_half = label;
	Vertice V_half(index_half, label_half, x_half, y_half);
	return(V_half);
}

double Triangle::area()const{
	vector<double> V1(V[0].get_coord()), V2(V[1].get_coord()), V3(V[2].get_coord());
	double S = 0.5 * abs( (V2[0] - V1[0]) * (V3[1] - V1[1]) - (V3[0] - V1[0]) * (V2[1] - V1[1]) ) ;
	return(S);
}

double Triangle::signed_area()const{
	vector<double> V1(V[0].get_coord()), V2(V[1].get_coord()), V3(V[2].get_coord());
	double S = 0.5 * abs((V2[0] - V1[0]) * (V3[1] - V1[1]) - (V3[0] - V1[0]) * (V2[1] - V1[1])) ;
	return(S);
}

Edge Triangle::edge_opp(int i)const{
	assert(i>=0 && i<=2);
	int j = (i+2)%3;
	int k = (i+1)%3;
	vector<Vertice> V_opp(2);
	V_opp = {V[j], V[k]};
	Edge e(0, V_opp);
	return(e);
}

Vector2d Triangle::GradLambda(int i)const{
	Edge E = edge_opp(i);
	Vector2d G(2);
	G[0] =  1.0/(signed_area()) * E.perp()[0];
	G[1] =  1.0/(signed_area()) * E.perp()[1];
	return(G);
}

Vector2d Triangle::GradLambdaP2(int i, double lambda1, double lambda2)const{
	double lambda3 = 1 - lambda1 - lambda2;
	Vector2d G;
	if (i==0){
		G = 4 * lambda1 * GradLambda(0);
	}
	if (i==1){
		G = 4 * lambda2 * GradLambda(1);
	}
	if (i==2){
		G = 4 * lambda3 * GradLambda(2);
	}
	else if (i==3){
		G = 4 * (GradLambda(0) * lambda1 + GradLambda(1) * lambda2);
	}
	else if (i==4){
		G = 4 * (GradLambda(1) * lambda2 + GradLambda(2) * lambda3);
	}
	else if (i==5){
		G = 4 * (GradLambda(2) * lambda3 + GradLambda(0) * lambda1);
	}

	return(G);
}       // H(i) dans le cours ???



vector<int> Triangle::get_vertices_index(){
	vector<int> vect_index;
	for(unsigned int k=0; k<V.size(); k++){
		vect_index.push_back(V[k].get_index());
	}
	for(unsigned int k=0; k<V_half.size(); k++){
		vect_index.push_back(V_half[k].get_index());
	}
	return(vect_index);
}

void Triangle::compute_middle(){
	if(V_half.size()!=3){
		for(int k=0; k<3; k++){
			Vertice M_half = edge_opp(k).Compute_middle();
			Add_Half_Vertices(M_half);
		}
	}
}

vector<double> Triangle::computeBarryCenter()const{
	vector<Vertice> T_vertices = get_vertices();
	double x_barry = (1.0/3.0) * (T_vertices[0].x() + T_vertices[1].x() + T_vertices[2].x());
	double y_barry = (1.0/3.0) * (T_vertices[0].y() + T_vertices[1].y() + T_vertices[2].y());
	vector<double> Barry_coord = {x_barry, y_barry};
	return(Barry_coord);
}


vector<double> Triangle::computeBarryCoord(double x, double y)const{
	double lambda1, lambda2, lambda3;
	vector<Vertice> T_vertices = get_vertices();
	double x1 = T_vertices[0].x();		//xn2-xn1
	double x2 = T_vertices[1].x();		//xn2-xn1
	double x3 = T_vertices[2].x();		//xn2-xn1
	double y1 = T_vertices[0].y();
	double y2 = T_vertices[1].y();
	double y3 = T_vertices[2].y();
	double detT = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);
	lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
	lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
	lambda3 = 1 - lambda1 - lambda2;
	vector<double> CoordBarry = {lambda1, lambda2, lambda3};
	return(CoordBarry);
}

vector<double> Triangle::computeBarryCoord(Vertice V)const{
	double x = V.x();
	double y = V.y();
	vector<double> CoordBarry = this->computeBarryCoord(x, y);
	return(CoordBarry);
}

Matrix2d Triangle::get_Jacobian_mat()const{
	vector<Vertice> T_vertices = get_vertices();		//3 sommets du triangle T
	double J11 = T_vertices[1].x() - T_vertices[0].x();		//xn2-xn1
	double J12 = T_vertices[2].x() - T_vertices[0].x();		//xn3-xn1
	double J21 = T_vertices[1].y() - T_vertices[0].y();		//yn2-yn1
	double J22 = T_vertices[2].y() - T_vertices[0].y();		//yn3-yn1
	Matrix2d J;
	J << J11, J12, J21, J22;
	return(J);
}

Matrix2d Triangle::get_Jacobian_inv_mat()const{
	Matrix2d J = get_Jacobian_mat();
	Matrix2d J_inv;
	J_inv << J(1,1), -J(0,1),  -J(1,0), J(0,0);
	J_inv = (1.0/(J.determinant())) * J_inv;
	return(J_inv);
}

double Triangle::get_Jacobian_det()const{
	Matrix2d J = get_Jacobian_mat();
	double det = J.determinant();
	return(det);
}


void Triangle::display()const{
	cout<<"Triangle coordonnée"<<endl;
	cout<<V[0].x()<<"	"<<V[0].y()<<endl;
	cout<<V[1].x()<<"	"<<V[1].y()<<endl;
	cout<<V[2].x()<<"	"<<V[2].y()<<endl;
}
