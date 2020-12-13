#include <suitesparse/umfpack.h>
#include "EF2d.hpp"

#include <cassert>
#include <fstream>
#include <vector>
using namespace std;
int debug =0;





double g(const R2 & P){ return 2.;}
int  main(int argc, const char** argv)
{

    assert(argc>=2);
    const char sumfpack[]="UMFPACK", *ts=sumfpack;
    if(argc>2) ts=argv[2];
     if(argc>3) debug =1;
     cout << " lecture de " << argv[1] << endl;
    Mesh2d Th(argv[1]);

    vector<int> num(Th.nt*6);
    int n=BuildNum(Th,1,1,0,num);
    cout << " ndof = " << n << endl;
    const QF2 &qf=QuadratureFormular_T_5;
    vector<double> WP2(1),D(1);
    SetWP2(qf,WP2,D);
    vector<ItemOP> OPA,OPM;
    AddLap(OPA,1);
    AddMass(OPM,1);

    hashMatrix  A(n),M(n);
    Add2MatLap2QF(A,WP2,D, OPA, num, Th);
    Add2MatLap2QF(M,WP2,D, OPM, num, Th);
    vector<double> u(n),b(n),fh(n);
    for(int i=0; i< n; ++i)
    {
        u[i]=0;
        fh[i]=1;

    }


    double tgv = 1e30;
    cout << " fh: min " << min(&fh[0],n)  << " max " << max(&fh[0],n)   << endl;
    ProduitMatVec(M,&fh[0],&b[0]);
    cout << " b: min " << min(&b[0],n)  << " max " << max(&b[0],n)   << endl;

    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            if ( Th.OnGamma(k,e)==2)
            {
                int ip0=(e+1)%3,ip1=(e+2)%3,ipa=e+3;//  les 3 DL de l'aretes
                int i0=Th(k,ip0), i1 = Th(k,ip1), ia=Th[k].BAdj(e);
                int iAp[]={ip0,ip1,ipa};
                if(debug)
                    cout << " on 2 : " << k << " " << e << " " << "ia =" << ia  << " (" << i0 << " " << i1 <<")"<< endl;
                R2 P[3];
                P[0]=Th(i0);
                P[1]=Th(i1);
                P[2]=(P[0]+P[1])*.5;

                for(int ii=0; ii<3; ++ii)
                    {
                        int i=  num[k*6+iAp[ii]];
                        A(i,i)=tgv;
                        double gii =g(P[ii]);
                        b[i]= tgv*gii;
                        u[i]=gii;
                    }
            }
    double eps=1e-8;
    SparseLinearSolver<int,double> ALS(A,ts,"eps",eps,nullptr);

    if(debug)
    {
    cout << " M  " << M << endl;
    cout << " A  " << A << endl;
    cout << " b  " << b << endl;
    }
    ALS.solve(&u[0],&b[0]);

    ProduitMatVec(A,&u[0],&fh[0]);
    for(int i=0; i< n; ++i)
    {
        fh[i] -= b[i];
    }
    double dmin=min(&fh[0],n) , dmax=max(&fh[0],n);
    cout << " diff min " << dmin  << " max " << dmax<< endl;
    assert( (abs(dmin) < 1e-10) && (abs(dmax) < 1e-10)) ;

    cout << " min u " << min(&u[0],n)  << " max " << max(&u[0],n)   << endl;

    {
        ofstream f("gp.txt");
        for(int k=0; k<Th.nt; ++k)
        {  const Simplex & K=Th[k];
            for(int ip=0; ip<=3; ip++)
                f <<  (const R2 &) K[ip%3]<< " "<< u[num[6*k+(ip%3)]]<< endl;
            f << endl << endl;
        }

    }
    {
        ofstream f("uk2.txt");
        for(int k=0; k<Th.nt; ++k)
        {
            for(int ip=0; ip<6; ip++)
                f << u[num[6*k+ip]] << " ";
            f << endl ;
        }

    }
    //  pour lance la verification freefem++  unix
    char ff[256];
    sprintf(ff,"FreeFem++ freefem/plotuk2.edp -ns -Th %s ",argv[1]);
    system(ff);


}
