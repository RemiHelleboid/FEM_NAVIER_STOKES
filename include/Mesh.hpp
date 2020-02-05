#pragma once
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <functional>
#include <cmath>
#include "Geometry.hpp"

using namespace std;

typedef vector<double> vect;


class Mesh{
	private:
		unsigned int n_v; //Number of vertices
		unsigned int n_t; //Number of triangles
		unsigned int n_s; //Number of bound edges
		unsigned int n_e; //Number of edges
		vector<Vertice> Vertices;
		vector<Triangle> Triangles;
		vector<Edge> BoundEdges;
		vector<Edge> AllEdges;

	public:
		Mesh(string filename);
		void show_parameters();
		void display();
		void build_edges();
		void build_half_bridges();
		int get_nt()const{return(n_t);};
		int get_nv()const{return(n_v);};
		int get_nbe()const{return(n_s);};
		int get_nae()const{return(n_e);};
		Triangle get_triangle(int k)const{return(Triangles[k]);};	//retourne le k-ieme triangle
		void set_triangle(int k, Triangle T_new){Triangles[k] = T_new;};	//retourne le k-ieme triangle
		Edge get_edge(int k)const{return(AllEdges[k]);};	//retourne le k-ieme triangle
		Vertice get_vertice(int k)const{return(Vertices[k]);};
};
