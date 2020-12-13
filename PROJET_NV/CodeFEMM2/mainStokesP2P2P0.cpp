
#include "EF2d.hpp"

#include <cassert>
#include <fstream>

using namespace std;
int debug =0;
// fonction CL
double g(const R2 & P, int label,int c)
{
    if( (label == 3)  && (c ==0) ) return P.x*(1.-P.x)*4;
    return 0.;
}
int  main(int argc, const char** argv)
{

    assert(argc>=2);
    const char sumfpack[]="UMFPACK", *ts=sumfpack;
    if(argc>2) ts=argv[2];
    if(argc>3) debug =1;
    cout << " lecture de " << argv[1] << endl;
    Mesh2d Th(argv[1]);

    vector<int> nump(Th.nt),numu(Th.nt*6);
    int n2=BuildNum(Th,1,1,0,numu);// p2
    int n1=BuildNum(Th,0,0,1,nump);// p0
    int n=n2+n2+n1;
    cout << " n1 = " << n1 << endl;
    cout << " n2 = " << n2 << endl;
    cout << " ndof = " << n2+n2+n1 << endl;
    const QF2 &qf=QuadratureFormular_T_5;
    vector<double> WPp(1),WPu(1),D(1);
    double nu = 0.01, epsp = - 1e-7;
    int ndofp,ndofu;
    ndofu=SetWP2(qf,WPu,D);
    ndofp=SetWP0(qf,WPp,D);
    vector<ItemOP> OPA;
    AddLap(OPA,nu,0,0);
    AddLap(OPA,nu,1,1);
    AddMass(OPA,epsp,2,2);
    AddPdivV(OPA,-1,2,0,1);
    AdddivUQ(OPA,-1,0,1,2);
    cout << " OP A ="<< OPA << endl;

    //  Data block
    const vector<double> * BWP[] = { &WPu,&WPu,&WPp};
    const vector<int> * Bnum[]  = { & numu,&numu,&nump};
    int offset[]={ 0, n2, n2+n2, n2+n2+n1 };
    int BDofK[]={ndofu,ndofu,ndofp};
    hashMatrix  DataA(n);
    Add2MatLap2QF(OPA,DataA,BWP,D, Bnum,offset, Th);
    vector<double> u(n),b(n),fh(n);
    for(int i=0; i< n; ++i)
    {
        u[i]=0;
        fh[i]=1;

    }





    double tgv = 1e30;
    cout << " fh: min " << min(&fh[0],n)  << " max " << max(&fh[0],n)   << endl;
    int label;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            if ( (label=Th.OnGamma(k,e))  != Mesh2d::NotonGamma) // sur le bord
            {
                int ip0=(e+1)%3,ip1=(e+2)%3,ipa=e+3;//  le 3 DL de l'aretes
                int i0=Th(k,ip0), i1 = Th(k,ip1), ia=Th[k].BAdj(e);
                if(debug)
                    cout << " on 2 : " << k << " " << e << " " << "ia =" << ia  << " (" << i0 << " " << i1 <<")"<< endl;
                R2 P[3];
                P[0]=Th(i0);
                P[1]=Th(i1);
                P[2] = (P[0]+P[1])*.5;
                int NumOnEdge[6]={-1,-1,-1,-1,-1,-1};
                NumOnEdge[ip0]=0;
                NumOnEdge[ip1]=1;
                NumOnEdge[ipa]=2;

                for(int c=0; c<2; ++c) // CL sur  2 composante ...
                {
                    const  vector<int> & num= *Bnum[c];
                    int dofk=BDofK[c];
                    int oc=offset[c], j;
                    for(int ip=0; ip<dofk; ++ip)
                        if( (j=NumOnEdge[ip])>=0)
                        {

                          int i =  num[k*dofk+ip]+oc;
                            DataA(i,i)=tgv;
                            double bc =g(P[j],label,c);
                            b[i]=bc*tgv;
                            u[i]=bc;
                        }
                }
            }
    if(debug)
    {
    cout << " A  " << DataA << endl;
    cout << " b  " << b << endl;
    }
    double eps=1e-8;
    SparseLinearSolver<int,double> A(DataA,ts,"eps",eps,nullptr);

    if(debug)
    {
        cout << " A  " << DataA << endl;
        cout << " b  " << b << endl;
    }
    A.solve(&u[0],&b[0]);

    cout << " A: min " << min(&DataA.aij[0],DataA.nnz)  << " max " << max(&DataA.aij[0],DataA.nnz)   << endl;



    cout << " min u1:  " << min(&u[0],n2)  << " max " << max(&u[0],n2)   << endl;
    cout << " min u2:  " << min(&u[n2],n2)  << " max " << max(&u[n2],n2)   << endl;
    cout << " min p:  " << min(&u[2*n2],n1)  << " max " << max(&u[n2*2],n1)   << endl;


    {
        ofstream f("u1u2p0.txt");
        for(int k=0; k<Th.nt; ++k)
        {
            for(int c=0; c<3; ++c)
            {
                const  vector<int> & num=*Bnum[c];
                int dofk=BDofK[c];
                int oc=offset[c];

                for(int ip=0; ip<dofk; ip++)
                  f << u[num[dofk*k+ip]+oc] << " ";
            }
            f << endl ;
        }

    }
    //  pour lance la verification freefem++  unix
    char ff[256];
    sprintf(ff,"FreeFem++ freefem/plotu1u2p0.edp -ns -Th %s ",argv[1]);
    system(ff);

}
