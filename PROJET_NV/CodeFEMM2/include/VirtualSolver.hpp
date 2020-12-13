#pragma once
#ifndef __VirtualSolver_HPP__
#define __VirtualSolver_HPP__
#include <map>
#include <vector>
#include <string>
#include "cholmod.h"

struct Data_Sparse_Solver {
    typedef void * pcommworld;
    static std::map<std::string,int> *  mds;
    typedef std::map<std::string,int>::const_iterator IMDS;
    bool initmat;  // 0
    //  TypeSolveMat* typemat;
    double epsilon; //1
    const void * precon; //3
    int NbSpace;//4
    int strategy;// 5
    double tgv;//6
    bool factorize;
    double tol_pivot;//8
    double tol_pivot_sym;//9
    int itmax ; //10
    int verb;
    std::string data_filename;
    std::vector<long> lparams;  //  copy arry more secure ...
    std::vector<double> dparams;

    std::map<std::string,std::string> * smap;

    std::vector<long> perm_r;
    std::vector<long> perm_c;
    std::vector<double> scale_r;
    std::vector<double> scale_c;
    std::string sparams;
    pcommworld commworld;  // pointeur sur le commworld
    int master; //  master rank in comm add FH 02/2013 for MUMPS ... => VDATASPARSESOLVER exist
    std::vector<double> * rinfo;
    std::vector<long> * info;
    Data_Sparse_Solver()
    :
    initmat(1),
    epsilon(1e-6),
    // typemat(0),
     precon(0),
    NbSpace(50),
    strategy(0),
    tgv(1e30),
    factorize(0),
    tol_pivot(-1),
    tol_pivot_sym(-1),
    itmax(0),
    verb(verbosity),
    smap(0) ,
    commworld(0),
    master(0),
    rinfo(0),
    info(0)
    {}

    static   std::map<std::string,int> * Set_mds()
    {
        using namespace std;
        static int init=0;
        static map<string,int> ms;

        if( init==0)
        {
            ms["eps"]=1;
            ms["precon"]=3;
            ms["nbkrilov"]=4;
            ms["nbspace"]=4;
            ms["strategy"]=5;
            ms["tgv"]=6;
            ms["tol_pivot"]=8;
            ms["tol_pivot_sym"]=9;
            ms["itmax"]=10;
            ms["data_filename"]=11;
            ms["verbosity"]=12;
        }
        return & ms;
    }
    void  set(va_list ap)
    {
        int k=0;
        const char *what;

        while (( what=va_arg(ap,const  char *) ))
        {
            k++;
            assert(k<20);
            if( mds==0) mds = Set_mds();
            if( what==0 ) break;
            IMDS iw = mds->find(what);
            if(iw == mds->end()) {break;}
            int cas = iw->second;
            switch (cas)
            {
                case 1: epsilon=va_arg(ap,double);
                    cout << " ds : eps = " << epsilon << endl; break;
                case 3: precon=va_arg(ap,void *);
                    cout << " ds : precon = " << precon << endl; break;
                case 4: NbSpace=va_arg(ap,int);
                    cout << " ds : nbspace (krilov)  = " << NbSpace << endl; break;
                case 5: strategy=va_arg(ap,int);
                    cout << " ds : strategy  = " << strategy << endl; break;
                case 6: tgv=va_arg(ap,double);
                     cout << " ds : tgv  = " << tgv  << endl; break;
                case 8: tol_pivot=va_arg(ap,double);  cout << " ds : tol_pivot  = " << tol_pivot << endl; break;
                case 9: tol_pivot_sym=va_arg(ap,double); cout << " ds : tol_pivot_sym  = " << tol_pivot_sym << endl;  break;
                case 10: itmax=va_arg(ap,int);
                    cout << " ds : itmax  = " << itmax << endl; break;
                case 12: verbosity=va_arg(ap,int);
                    cout << " ds : itmax  = " << itmax << endl; break;
                    // case 11: data_filename=va_arg(ap,string); break;
            }
        }
    }
private:
    Data_Sparse_Solver(const Data_Sparse_Solver& ); // pas de copie

};


template<class I, class R>
class VirtualSolver {
public:
    int state;
    VirtualSolver() : state(0),codeini(0),codesym(0),codenum(0) {}
    long codeini,codesym,codenum;
    virtual void dosolver(R *x,R*b,int N=1) =0;

    virtual void fac_init(){}  // n, nzz fixe
    virtual void fac_symbolic(){} //  i,j fixe
    virtual void fac_numeric(){}   // a fixe
    virtual void SetState(){}

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

    virtual ~VirtualSolver(){}
};

#endif
