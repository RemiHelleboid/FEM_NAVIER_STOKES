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
	if(i<=j){
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
	double S = 0.5 * abs( (V2[0] - V1[0])*(V3[1] - V1[1]) - (V3[0]-V1[0])*(V2[1]-V1[1]) ) ;
	return(S);
}

Edge Triangle::edge_opp(int i){
	assert(i>=0 && i<=3);
	int j = (i+1)%3;
	int k = (i+2)%3;
	vector<Vertice> V_opp{V[j], V[k]};
	Edge e(0, V_opp);
	return(e);
}

vector<double> Triangle::GradLambda(int i){
	Edge E = edge_opp(i);
	vector<double> G(2);
	G[0] =  1.0/(2*area())*E.perp()[0];
	G[1] =  1.0/(2*area())*E.perp()[1];
	return(G);
}

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

void Triangle::Display(){
	cout<<"Triangle coordonnée"<<endl;
	cout<<V[0].x()<<"	"<<V[0].y()<<endl;
	cout<<V[1].x()<<"	"<<V[1].y()<<endl;
	cout<<V[2].x()<<"	"<<V[2].y()<<endl;
}
