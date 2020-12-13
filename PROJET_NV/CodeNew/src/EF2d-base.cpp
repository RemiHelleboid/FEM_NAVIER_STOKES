#include "EF2d-base.hpp"
#include <map>
#include <fstream>
#include <cstdlib>
#include <cassert>
using namespace std;
extern  int debug;

void Mesh2d::BuildAdjacent()
{
    int fait = std::numeric_limits<int>::min() ; //-2^31
    map< pair<int,int>, int> arete;
    int err=0; // nb of error
    for(int k=0; k< nbe; ++k)
    {
        int i0 = be(k,0);
        int i1 = be(k,1);
        if(i1<i0) std::swap(i0,i1);
        pair<int,int>  a(i0,i1);
        if(debug)
            cout << " b " << k << ":  (" << i0 << "," << i1 << ") " << b[k].lab <<endl;
        if(arete.find(a) != arete.end())
        {
            err++;
            std::cerr << " Erreur l'arete de bord existe 2 fois "<< i0 << " " << i1 << std::endl;
        }
        else
            arete[a]= -1-k;
    }
    
    for(int k=0; k< nt; ++k)
        for(int i=0; i< 3; ++i)
        {
            int i0 = operator()(k,(i+1)%3);
            int i1 = operator()(k,(i+2)%3);
            if(i1<i0) std::swap(i0,i1);
            pair<int,int>  a(i0,i1);
            if(arete.find(a) == arete.end())
            {
                arete[a]= 3*k + i;
            }
            else
            {
                int cadj = arete[a];
                if( cadj >= 0) { // arete interne
                    int kk = cadj/3, ii=cadj%3;
                    t[k ].dataadj[i ]=kk*3+ii;
                    t[kk].dataadj[ii]=k *3+i ;
                    arete[a]=fait; // fait !!
                }
                else if( cadj == fait)
                {
                   std::cerr << " Erreur l'arete entre plus de 2 objet  "<< i0 << " " << i1 << std::endl;
                   err++;
                }
                else
                {
                   if(debug)
                       cout << " AdjB " << k << " " << i << " :  " << -cadj-1 << " (" << i0 << " , "<< i1 << ")" <<  endl;
                  t[k].dataadj[i]=cadj;
                }
                
            }
        }
    // verif to le adj on ete defini ...
    int erra =0;
    for(int k=0; k< nt; ++k)
        for(int i=0; i< 3; ++i)
            if(t[k].dataadj[i]== fait) erra++;
    err += erra;
    if( err)
    {
        std::cout << " nb arete element sans adjacent " << erra << std::endl;
        std::cout << " erreur dans la construction des adjacent "<< std::endl;
        abort();
    }
        
}
Mesh2d::Mesh2d(const char * filename)
{
  std::ifstream  f(filename); 
  assert( f); 
  int I[4] ; 
  f >> nv >> nt >> nbe ;
  assert( f.good());
  t = new Simplex[nt];
  v = new Vertex[nv];
  b= new BorderSimplex[nv];
  assert( t && v);
  double mes =0, mesb=0;
  for(int i=0;i<nv;++i)
    { 
      f >> v[i] ; 
      assert( f.good());
    }

  for(int k=0;k<nt;++k)
    { 
      for(int i=0;i< 4; ++i)
	f >> I[i] ; 
      assert( f.good());
      t[k].build(v,I,-1);
      mes += t[k].mes; 
    }
    for(int k=0;k<nbe;++k)
    {
        for(int i=0;i< 3; ++i)
            f >> I[i] ;
        assert( f.good());
        b[k].build(v,I,-1);
        mesb += b[k].mes;
    }
  BuildAdjacent();
  std::cout<< " End read " << nv << " " << nt << " mes =" << mes << " mesb = " << mesb <<  std::endl;
}
