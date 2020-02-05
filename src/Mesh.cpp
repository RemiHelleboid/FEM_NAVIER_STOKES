#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cassert>
#include <vector>
#include <functional>
#include <cmath>
#include "Mesh.hpp"

using namespace std;


Mesh::Mesh(string filename){
	ifstream meshfile(filename);
	assert(meshfile.is_open() == 1);
	meshfile >> n_v;
	meshfile >> n_t;
	meshfile >> n_s;

	Vertices.resize(n_v);
	Triangles.resize(n_t);
	BoundEdges.resize(n_s);

	for(unsigned int l=0; l<n_v; l++){
		vector<double> V_buf(3);
		meshfile >> V_buf[0];	//abcisse
		meshfile >> V_buf[1];	//ordonnée
		meshfile >> V_buf[2];	//boundary label
		int global_index = l;
		Vertice V(global_index, V_buf[2], V_buf[0], V_buf[1]);
		Vertices[l] = V;
	}

	for(unsigned int l=0; l<n_t; l++){
		vector<int> V_index(3);
		int label;
		meshfile >> V_index[0];	//index of the first vertice of the triangle
		meshfile >> V_index[1];	//second one
		meshfile >> V_index[2];	//third one
		meshfile >> label;	//region label
		vector<Vertice> T_vertices({Vertices[V_index[0]-1], Vertices[V_index[1]-1], Vertices[V_index[2]-1]});
		Triangle T(label, T_vertices);
		Triangles[l] = T;
	}

	for(unsigned int l=0; l<n_s; l++){
		vector<int> E_index(3);
		int label;
		meshfile >> E_index[0];	//index of the first vertice of the edge on a bound
		meshfile >> E_index[1];	//second one
		meshfile >> label;	//boundary label
		vector<Vertice> E_vertices({Vertices[E_index[0]-1], Vertices[E_index[1]-1]});
		Edge E(label, E_vertices);
		BoundEdges[l] = E;
	}
	meshfile.close();
	build_edges();
}

void Mesh::show_parameters(){
	cout<<"Number of triangles : "<<n_t<<endl;
	cout<<"Number of vertices : "<<n_v<<endl;
	cout<<"Number of edges on bound : "<<n_s<<endl;
	cout<<"Number of all edges: "<<n_e<<endl;
}

void Mesh::display(){
		for(unsigned int l=0; l<n_v; l++){
			cout << Vertices[l].x()<<" ";	//abcisse
			cout << Vertices[l].y()<<" ";	//ordonnée
			cout << Vertices[l].get_label()<<endl;	//boundary label
		}
}
		//for(unsigned int l=0; l<n_t; l++){
			//cout << Triangles[l][0]<<" ";	//index of the first vertice of the triangle
			//cout << Triangles[l][1]<<" ";	//second one
			//cout << Triangles[l][2]<<" ";	//third one
			//cout << Triangles[l][3]<<endl;	//region label
		//}

		//for(unsigned int l=0; l<n_s; l++){
			//cout << BoundEdges[l][0]<<" ";	//index of the first vertice of the edge on a bound
			//cout << BoundEdges[l][1]<<" ";	//second one
			//cout << BoundEdges[l][2]<<endl;	//boundary label
		//}
	//}

void Mesh::build_edges(){
	vector<Edge> Edges;
	for(unsigned int k=0; k<n_t; k++){
		vector<Vertice> V(3);
		V = Triangles[k].get_vertices();
		vector<Vertice> V0({V[0], V[1]});
		vector<Vertice> V1({V[0], V[2]});
		vector<Vertice> V2({V[1], V[2]});
		vector<Edge> Ek(3);
		Ek[0] = Edge(0, V0);
		Ek[1] = Edge(0, V1);
		Ek[2] = Edge(0, V2);

		if (!count(Edges.begin(), Edges.end(), Ek[0])){
			Edges.push_back(Ek[0]);
		}
		if (!count(Edges.begin(), Edges.end(), Ek[1])){
			Edges.push_back(Ek[1]);
		}
		if (!count(Edges.begin(), Edges.end(), Ek[2])){
			Edges.push_back(Ek[2]);
		}
	}
	n_e = Edges.size();
	AllEdges = Edges;
}

void Mesh::build_half_bridges(){
	for(unsigned int k=0; k<n_t; k++){		//boucle sur les triangles
		Triangle T_k = Triangles[k];
		vector<Vertice> V_half_Tk;
		for(int i=0; i<3; i++){				//boucle sur les sommet
			Edge E_opp_i = T_k.edge_opp(i);
			Vertice V_half_i = E_opp_i.Compute_middle();
			int label_half_i = E_opp_i.get_label();

			V_half_i.set_label(label_half_i);
			if (!count(Vertices.begin(), Vertices.end(), V_half_i)){
				V_half_i.set_index(n_v);
				int new_label = int( E_opp_i.get_vertices()[0].get_label() == 1 && E_opp_i.get_vertices()[1].get_label()==1 );
				V_half_i.set_label(new_label);
				Vertices.push_back(V_half_i);
				V_half_Tk.push_back(V_half_i);
				T_k.Add_Half_Vertices(V_half_i);
				n_v = n_v + 1;

			}
			else{
				auto it = find(Vertices.begin(), Vertices.end(), V_half_i);
				int index_V = std::distance(Vertices.begin(), it);
				int new_index = Vertices[index_V].get_index();
				int new_label = Vertices[index_V].get_label();
				V_half_i.set_index(new_index);
				V_half_i.set_label(new_label);
				V_half_Tk.push_back(V_half_i);
				T_k.Add_Half_Vertices(V_half_i);
			}
		set_triangle(k, T_k);
		}
	}
	cout<<"n__________v  "<<n_v<<endl;

}
