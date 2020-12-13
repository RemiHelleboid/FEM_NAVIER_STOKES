#include "umfpack.h"
#include "EF2d.hpp"
#include <cmath>
#include <cstdlib>

#include <cassert>
#include <fstream>

using namespace std;
int debug =0;

extern double *pdebugF;
double g(R2 P) { double x=P.x*M_PI, y=P.y*M_PI*2; return  cos(x)*cos(y);}

double FF1(double x,double y,double u,double ux, double uy) { return 1./sqrt(1+ux*ux+uy*uy);}
double FF(double x,double y,double u,double ux, double uy) { return sqrt(1+ux*ux+uy*uy);}



int  main(int argc, const char** argv)
{
    // int2d(Th) ( (grad(u),grad(v))/sqrt( 1+(grad(u),grad(u)) = O, u= g sur Gamma
    const char sumfpack[]="UMFPACK", *ts=sumfpack;
    if(argc>3) ts=argv[3];
    debug =0;
    cout << " lecture de " << argv[1] << endl;
    assert(argc>=2);
    if(argc>3) debug =atoi(argv[4]);;
    int niter=100;
    if( argc>2) niter = atoi(argv[2]);
    cout << " lecture de " << argv[1] << " niter = "<< niter << " debug =" << debug << " ts = "<<*ts <<   endl;
    Mesh2d Th(argv[1]);

    vector<int> num(Th.nt*6);
    int n=BuildNum(Th,1,1,0,num);
    cout << " ndof = " << n << endl;
    const QF2 &qf=QuadratureFormular_T_5;
    vector<double> F(qf.n*Th.nt);
    if(debug>2)pdebugF= & F[0]; // for debug of function quadarature

    vector<double> WP2(1),D(1);
    SetWP2(qf,WP2,D);
    vector<double> u(n),b(n),uo(n);
    for(int i=0; i< n; ++i)
    {
        u[i]=0;
        b[i]=0;
    }


    vector<ItemOP> OPA;

    AddLap(OPA,1,0,0,&F[0]);
    double err=0;
    for(int iter =0; iter< niter; ++iter)
    {
        uo = u;
        SetF(Th,qf,WP2,num,u,F,FF1);
        hashMatrix  A(n);
        Add2MatLap2QF(A,WP2,D, OPA, num, Th);


        if( debug) cout << " A: min " << min(&A.aij[0],A.nnz)  << " max " << max(&A.aij[0],A.nnz)   << endl;



        double tgv = 1e30;
        if(debug) cout << " b: min " << min(&b[0],n)  << " max " << max(&b[0],n)   << endl;
        int label;
        for(int k=0; k< Th.nt; ++k)
            for(int e=0; e<3; ++e)
                if ( (label=Th.OnGamma(k,e))!= Mesh2d::NotonGamma)
                {
                    int ip0=(e+1)%3,ip1=(e+2)%3,ipa=e+3;//  les 3 DL de l'aretes
                    int i0=Th(k,ip0), i1 = Th(k,ip1), ia=Th[k].BAdj(e);
                    int iAp[]={ip0,ip1,ipa};
                    if(debug)
                        cout << " on  : " << label << " / " << k << " " << e << " " << "ia =" << ia  << " (" << i0 << " " << i1 <<")"<< endl;
                    R2 P[3];
                    P[0]=Th(i0);
                    P[1]=Th(i1);
                    P[2]=(P[0]+P[1])*.5;

                    for(int ii=0; ii<3; ++ii)
                    {
                        double bc=g(P[ii]);
                        int i=  num[k*6+iAp[ii]];
                        A(i,i)=tgv;
                        b[i]= tgv*bc;
                        u[i]=bc;
                    }
                }
        double eps=1e-8;
        SparseLinearSolver<int,double> ALS(A,ts,"eps",eps,nullptr);

        if(debug)
        {
            cout << " A  " << A << endl;
            cout << " b  " << b << endl;
        }
        ALS.solve(&u[0],&b[0]);

        if(debug>3)
        {
            cout << " A  " << A << endl;
            cout << " b  " << b << endl;
        }

        err = normsum(&u[0],&uo[0],n);
        SetF(Th,qf,WP2,num,u,F,FF);
        cout << iter << " min u " << min(&u[0],n)  << " max " << max(&u[0],n)   << " err = "<< err << " J = " << Integral(Th,D,F) << endl;
        if(err < 1e-8) break;
    }
    {
        ofstream f("gpminsurf.txt");
        for(int k=0; k<Th.nt; ++k)
        {  const Simplex & K=Th[k];
            for(int ip=0; ip<=3; ip++)
                f <<  (const R2 &) K[ip%3]<< " "<< u[num[6*k+(ip%3)]]<< endl;
            f << endl << endl;
        }

    }
    {
        ofstream f("minsuf.txt");
        for(int k=0; k<Th.nt; ++k)
        {
            for(int ip=0; ip<6; ip++)
                f << u[num[6*k+ip]] << " ";
            f << endl ;
        }

    }
    //  pour lance la verification freefem++  unix
    char ff[256];

    sprintf(ff,"FreeFem++ freefem/plotSurfMin.edp -ns -Th %s ",argv[1]);
    printf("%s\n",ff);
    system(ff);


}
