#pragma once
inline void MATERROR(int i,const char *cmm)
{
    std::cerr << " MATERROR " << i << " : " << cmm << std::endl;
    std::abort();
}

template<class TypeIndex=int,class TypeScalar=double>
class VirtualMatrix {
public:
    typedef TypeIndex I;
    typedef TypeScalar R;

    I n,m; // size of matrix
    long state,codeini,codesym,codenum;
    VirtualMatrix(I NN,I MM=-1) : n(NN),m(MM<0 ? n : MM),
      state(0),codeini(0),codesym(0),codenum(0) {}
   static R *Set2Const(I n,R *x,R c=R()) { std::fill(x,x+n,c); return x;}

  virtual void dosolver(R *x,R*b,int N=1)
    {  MATERROR(1,"VirtualMatrix:: no solver ?????"); }
  virtual void fac_init(){}  // n, nzz fixe
  virtual void fac_symbolic(){} //  i,j fixe
  virtual void fac_numeric(){}   // a fixe
  virtual void SetState(){}
  virtual R* addMatMul(R *x,R*Ax) const {  MATERROR(1,"VirtualMatrix:: no AddMatMul ?????"); return 0;}
  virtual R* addMatTransMul(R *x,R*Atx) const {  MATERROR(1,"VirtualMatrix:: no addMatTransMul ?????");return 0;}
  R* MatMul(R *x,R*Ax) const { return addMatMul(x, Set2Const(n,Ax)); }
  R* MatTransMul(R *x,R*Atx) const { return addMatMul(x, Set2Const(m,Atx));}
  virtual bool WithSolver() const {return false;} // by default no solver

  void CheckState(long ci=0,long cs=0, long cn=0)
    {
      if(ci &&  ci != codeini) { codeini=ci; state=0;}// n, nzz fixe
      else if(cs &&  cs != codesym) { codesym=cs; state=1;}//  i,j fixe
      else if(cn &&  cn != codenum) { codenum=cn; state=2;}// a fixe
    };

  R* solve(R *x,R *b,int N=1)
  {
    SetState();
    if( state==0) {fac_init(); state=1;}
    if( state==1) {fac_symbolic(); state=2;}
    if( state==2) {fac_numeric();state=3;}
    dosolver(x,b,N);
    return x;
  }

  virtual ~VirtualMatrix(){}
};


template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> *A,TypeScalar *x, TypeScalar *Ax) { return A->MatMul(x,Ax);}
template<class TypeIndex=int,class TypeScalar=double>
inline double * ProduitMatVec(const VirtualMatrix<TypeIndex,TypeScalar> &A,TypeScalar *x, TypeScalar *Ax) { return A.MatMul(x,Ax);}
