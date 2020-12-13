#pragma once
#include "EF2d-base.hpp"
#include "QuadratureFormular.hpp"
#include <map>
#include <vector>
#include <cassert>
#include "SparseLinearSolver.hpp"
typedef  HashMatrix<int,double> hashMatrix;
//#include "SparseMat.hpp"
using  std::vector;
extern int debug ;
typedef GQuadratureFormular<R2> QF2;
const int opV=0,opDX=1,opDY=2;
int  SetWP0(const R2 &P,double *w);// Value and derivative of B. Func
int  SetWP1(const R2 &P,double *w);
int  SetWP2(const R2 &P,double *w);
int  SetWPP2(const R2 &P,int k,double *WP2,const Mesh2d & Th);
int  SetWP2(const QF2 & qf,vector<double>& WP2,vector<double>& D);
int  SetWP1(const QF2 & qf,vector<double>& WP1,vector<double>& D);
int  SetWP0(const QF2 & qf,vector<double>& WP0,vector<double>& D);



void SetWK(const Simplex &K,const vector<double>& WP1,vector<double>& WP1K);

void SetDK(const Simplex & K,const vector<double>& D,vector<double>& DK);

struct ItemOP {
    double c;
    int opu,opv; // Operation
    int cu,cv;   // composante
    double *v;
};

 ostream & operator<<(ostream & f ,const ItemOP & o);
ostream & operator<<(ostream & f ,const vector<ItemOP> & vop);
void AddMass(vector<ItemOP> & vop,double alpha,int ccu=0,int ccv=-1,double *vv=0);
void AddLap(vector<ItemOP> & vop,double beta,int ccu=0,int ccv=-1,double *vv=0);
void AddPdivV(vector<ItemOP> & vop,double beta,int cp,int cv1,int cv2,double *v0=0,double *v1=0);
void AdddivUQ(vector<ItemOP> & vop,double beta,int cp,int cv1,int cv2,double *v0=0,double *v1=0);

double normsum(double *u,double *v, int n);

double min(double *v, int n);
double max(double *v, int n);
int BuildNum(const Mesh2d & Th,int ndofs,int ndofa,int ndoft, vector<int> & num);

void AddMatSymConst(hashMatrix &A,int npq,int ndofKi,const int *iK,
                    const vector<ItemOP> & vop, const double *D,
                    const double *Wi,int ndofKj=0,const int *jK=0,const double *Wj=0,int oi=0,int oj=0);
/*
void AddMatSymConst(MatrixMap &A,int npq,int ndofK,const int *iK,
                    const vector<ItemOP> & vop, const double *D,const double *W)
;
*/
void Add2MatLap2QF(hashMatrix &A, const vector<double> & Wh,const vector<double> & Dh, vector<ItemOP>  vop,const vector<int> & num, const Mesh2d & Th);
void Add2MatLap2QF( const vector<ItemOP> & vop,hashMatrix &A, const vector<double> ** pWh,const vector<double> Dh,const vector<int> ** pnum,int *offset,const Mesh2d & Th);
double Integral(const Mesh2d & Th, const vector<double> &D, const vector<double> & F);

void SetF(const Mesh2d &Th,const QF2 &qf,const vector<double> &WP2,const vector<int>& num,const vector<double> & u, vector<double> &F, double (*FF)(double x,double y,double v));

template <typename T>
ostream & operator <<(ostream & f,const  vector<T> & v)
{
    for(int i=0; i< v.size();++i)
        f << i<< " " <<v[i] <<endl;
    return f;
}

double CalcUp(R2 Pp,int Kp,vector<double> u,const  vector<int> & num ,int ndofu,  int c,const Mesh2d & Th);
void Bartocart(const Mesh2d & Th,R2 & Phat,int Kp);
//double CalcUp2(R2 Pp,int Kp,vector<double> u,const  vector<int> & num ,int ndofu,  int c);
