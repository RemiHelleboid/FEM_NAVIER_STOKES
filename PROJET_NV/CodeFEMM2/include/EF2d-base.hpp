#pragma once
#include "R1.hpp"
#include "R2.hpp"
#include "clock.hpp"
#include <cassert>
#include <limits>
class Label {
public:
  int lab;
  int  OnGamma() const { return lab;}
  Label(int l=0) : lab(l) {}
};

class Vertex :public R2,public  Label
{
    typedef R2 Rd;
 public:
  Vertex() {}

private: // not copy operator
    Vertex(const Vertex &);
    Vertex& operator=(const Vertex &);

};

inline std::ostream& operator <<(std::ostream& f, const  Vertex & P )
{ return  f << P.x << ',' << P.y  << ',' << P.lab  ;}
inline  std::istream& operator >>(std::istream& f,  Vertex & P)
{ return f >>  P.x >>  P.y >> P.lab ;  }


class Simplex {
public:

  static const int nbv =3, d=2;
  Vertex * v[nbv];
  int dataadj[nbv]; // Attention les valeurs sont code ..
  int reg;
  double mes;
  Simplex(){ (v[0]=(v[1]=(v[2]=0)));}
  void build(Vertex *v0,int * I,int offset=0)
  {// I array of vertex number
      int bad = std::numeric_limits<int>::min() ; //-2^31
    for(int i=0; i < nbv; ++i)
      v[i] =  v0 + I[i]+offset;
    reg=I[nbv];
    mes = det(*v[0], *v[1], *v[2]) * 0.5;
    dataadj[0]=dataadj[1]=dataadj[2]=bad ; // 2^31-1..
    assert(mes>0) ;
  }
  void GradLambdaK(R2 *G) const
  {
    double K2 = mes*2;
    G[0] = R2(*v[1],*v[2]).perp()/K2;
    G[1] = R2(*v[2],*v[0]).perp()/K2;
    G[2] = R2(*v[0],*v[1]).perp()/K2;
  }

  Vertex & operator[](int i) { assert(i>=0 && i < nbv); return *(v[i]); }
  const Vertex & operator[](int i) const { assert(i>=0 && i < nbv); return *(v[i]);}

  const R2 &P(int i) const { return *v[i];} //
  int SimplexAdj(int i,int & ii) const // Element adjacent et par quel face
    { if (dataadj[i]<0) { return ii=-1; }
        ii=dataadj[i]%nbv;
        return dataadj[i]/nbv;}
  int BAdj(int i) //  Bord adjacent
    { if (dataadj[i]>=0) { return -1; }
        else return -dataadj[i]-1;}

  R2 operator()(R2 Phat) const { return (1-Phat.x-Phat.y)* P(0) + Phat.x* P(1) +  Phat.y* P(2) ; }
private: // not copy operator
    Simplex(const Simplex &);
    Simplex& operator=(const Simplex &);

};

class BorderSimplex {
public:

    static const int nbv =2, d=2;
    Vertex * v[nbv];
    int lab;
    double mes;
    BorderSimplex(){ (v[0]=(v[1]));}
    void build(Vertex *v0,int * I,int offset=0)
    {// I array of vertex number
        for(int i=0; i < nbv; ++i)
            v[i] =  v0 + I[i]+offset;
        lab=I[nbv];
        R2 AB(*v[0],*v[1]);
        mes = Norme2(AB);
        assert(mes>0) ;
    }

    Vertex & operator[](int i) { assert(i>=0 && i < nbv); return *(v[i]); }
    const Vertex & operator[](int i) const { assert(i>=0 && i < nbv); return *(v[i]);}

    const R2 &P(int i) const { return *v[i];}
    void GradLambdaK(R2 *G) const {  R2 AB(*v[0],*v[1]); G[0] = - (G[1] = AB);}
    R2 operator()(R1 Ph) const { return (1-Ph.x)* P(0) + Ph.x* P(1)  ; }
    int OnGamma() const {return lab;}
private: // not copy operator
    BorderSimplex(const BorderSimplex &);
    BorderSimplex& operator=(const BorderSimplex &);

};


class Mesh2d
{
public:
  static const  int NotonGamma=std::numeric_limits<int>::min() ;
  int nv,nt,nbe;
  Vertex * v;
  Simplex *t;
  BorderSimplex *b;
  Mesh2d(const char *  filename);
  ~Mesh2d() { delete [] v; delete [] t; delete [] b;}
  // destuctor => careful with copie operator
  // no copy operator
  // chech index number
  int CheckV(int i) const { assert( i>=0 && i < nv); return i; }
  int CheckT(int i) const { assert( i>=0 && i < nt); return i; }
  int CheckBE(int i) const { assert( i>=0 && i < nbe); return i; }
  int operator()(const Vertex & vv) const { return CheckV(&vv-v);}
  int operator()(const  Simplex & tt) const  { return CheckT(&tt-t);}
  int operator()(const Vertex * vv)const  { return CheckV(vv-v);}  // (1)
  int operator()(const  Simplex * tt) const { return CheckT(tt-t);}
  int operator()(const BorderSimplex * vv)const  { return CheckBE(vv-b);}  // (1)
  int operator()(const  BorderSimplex & tt) const { return CheckBE(&tt-b);}
  void getnum(int k, int *Ik) const { Ik[0]=t[k].v[0]-v; Ik[1]=t[k].v[1]-v; Ik[2]=t[k].v[2]-v;}
  Simplex & operator[](int k) { return t[CheckT(k)]; }
  const Simplex & operator[](int k) const { return t[CheckT(k)]; }
  const Vertex & operator()(int i) const { return v[CheckV(i)]; }
  const BorderSimplex & be(int i) const { return b[CheckBE(i)]; }
  int    be(int i, int j ) const { return operator()(b[CheckBE(i)][j]); }
  int  operator()(int k, int i) const { return  operator()(t[k].v[i]); }// call (1)
  void  BuildAdjacent();
  int OnGamma(int k, int a) const { int ie = t[k].BAdj(a); return ie >=0 ? b[ie].OnGamma(): NotonGamma;}
  static const int d=2;
private:// no copy too hard to build...
  Mesh2d(const Mesh2d &);
  Mesh2d & operator=(const Mesh2d &);

};
