#pragma once
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

typedef vector<double> vect;
typedef vector<int> vectint;


class Vertice{
    private:
        int index;	//numéro global du sommet entre 1 et n_v !!!!!
        int label;	//label du sommet
        vector<double> coordinates;   //(x, y) coord
    public:
//Constructeurs
        Vertice(){label=0; index = 0; coordinates.resize(2); coordinates[0]=0, coordinates[1]=0;};
        Vertice(int index0, int label0, double x0, double y0){coordinates.resize(2); coordinates[0]=x0, coordinates[1]=y0; index=index0; label=label0;};
        Vertice(const Vertice &V):index(V.index), label(V.label), coordinates(V.coordinates){};
//Methods
        vector<double> get_coord()const{return(coordinates);}
        double x()const{return(coordinates[0]);}
        double y()const{return(coordinates[1]);}
        int get_index()const{return(index);}
        int get_label()const{return(label);}
        void set_label(int new_label){label = new_label;}
        void set_index(int new_index){index = new_index;}
        Vertice & operator=(const Vertice &) = default;
        bool operator==(const Vertice & V0)const{return(V0.coordinates==coordinates);};
};

class Edge{
    private:
        int label;
        vector<Vertice> V;
    public:
        Edge(){label = 0; V = vector<Vertice>(2);}
        Edge(int l, vector<Vertice> V0);
        Edge(const Edge &E):label(E.label), V(E.V){};
        Edge & operator=(const Edge &) = default;
        bool operator==(const Edge &)const;
        vector<double> vect();
        vector<double> perp();
        vector<Vertice> get_vertices()const{return(V);};
        int get_label()const{return(label);};
        Vertice Compute_middle();		//retourne le milieu du segment
};

class Triangle{
    private:
        int label;
        vector<Vertice> V;
        vector<Vertice>	V_half;
    public:
        Triangle(){label=0; V = vector<Vertice>(3);};
        Triangle(int l, vector<Vertice> V0):label(l), V(V0){};
        Triangle(const Triangle &T):label(T.label), V(T.V), V_half(T.V_half){};
        Triangle& operator=(const Triangle &) = default;
//        void Add_Half_Vertices(vector<Vertice> V_half_add){V_half = V_half_add;};
        void Add_Half_Vertices(Vertice V_half_add){V_half.push_back(V_half_add);};
        double area()const;
        Edge edge_opp(int i);
        vector<double> GradLambda(int i);       // H(i) dans le cours ???
        vector<int> get_vertices_index();		//renvoie les numéros i,j,k,ij, ik, jk  global des sommets du triangle
        vector<Vertice> get_vertices()const{return(V);}
        void Display();
};
