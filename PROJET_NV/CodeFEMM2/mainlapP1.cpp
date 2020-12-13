#include <suitesparse/umfpack.h>
#include "EF2d.hpp"
#include "QuadratureFormular.hpp"
#include <map>
#include <vector>
#include <cassert>
#include <fstream>

using namespace std;
int debug =0;


struct ItemOper {
    double c;
    int opu,opv;
};
void AddMass(vector<ItemOper> & vop,double alpha)
{
    ItemOper op={alpha,0,0};
    vop.push_back(op);
}
void AddLap(vector<ItemOper> & vop,double beta)
{
    ItemOper op1={beta,1,1};
    ItemOper op2={beta,2,2};
    vop.push_back(op1);
    vop.push_back(op2);
}
void AddMatSymConst(hashMatrix &A,int npq,int ndofK,const int *iK,
                    const vector<ItemOper> & vop, const double *D,const double *W)
{
    int ndofK3=ndofK*3;
    for(int ip=0; ip<ndofK; ++ip)
        for(int jp=0; jp<ndofK; ++jp)
        {
            double aKij =0;
            for(int op=0; op<vop.size(); ++op)
                for(int p=0,k=0; p<npq; ++p, k+=ndofK3)
                {
                    aKij += D[p]*vop[op].c*W[k+ip*3+vop[op].opv]*W[k+jp*3+vop[op].opu];
                }
            int i = iK[ip];
            int j = iK[jp];
            //cout << i << " " << j << " " << aKij << endl;
            A(i,j) +=aKij;

        }
}


void Add2MatLap2QF(hashMatrix &A, const vector<double> & Wh,const vector<double> & Dh, const vector<ItemOper> & vop,const vector<int> & num, const Mesh2d & Th)
{
    int npq= Dh.size();
    int ndofK=num.size()/Th.nt;

    vector<double> WK(Wh.size()),DK(Dh.size());

    for(int k=0; k< Th.nt; ++k)
    {
        const Simplex &K=Th[k];
        SetWK(K,Wh,WK); SetDK(K,Dh,DK);
        AddMatSymConst(A,npq,ndofK,&num[k*ndofK],vop,&DK[0],&WK[0]);
    }

}



int  main(int argc, const char** argv)
{
    const char sumfpack[]="UMFPACK", *ts=sumfpack;
    if(argc>2) ts=argv[2];
    assert(argc>=2);
    if(argc>3) debug =1;
    cout << " lecture de " << argv[1] << endl;
    Mesh2d Th(argv[1]);
    int n = Th.nv;
    vector<int> num(Th.nt*3);
    BuildNum(Th,1,0,0,num);
    if(debug) //  pour avoir la meme numerotaion que freefem++
    for(int k=0,j=0; k<Th.nt; ++k)
        for(int i=0; i<3; ++i)
            num[j++]=Th(k,i);

    const QF2 &qf=QuadratureFormular_T_2;
    vector<double> WP1(1),D(1);
    SetWP1(qf,WP1,D);
    vector<ItemOper> OPA,OPM;
    AddLap(OPA,1);
    AddMass(OPM,1);

    hashMatrix  A(n),M(n);
    Add2MatLap2QF(A,WP1,D, OPA, num, Th);
    Add2MatLap2QF(M,WP1,D, OPM, num, Th);
    vector<double> u(n),b(n),fh(n);
    for(int i=0; i< n; ++i)
    {
        u[i]=0;
        fh[i]=1;

    }


    cout << " A: min " << min(&A.aij[0],A.nnz)  << " max " << max(&A.aij[0],A.nnz)   << endl;
    cout << " M: min " << min(&M.aij[0],M.nnz)  << " max " << max(&M.aij[0],M.nnz)   << endl;

    ProduitMatVec(&M,&fh[0],&b[0]);


    double tgv = 1e30;
    cout << " fh: min " << min(&fh[0],n)  << " max " << max(&fh[0],n)   << endl;

    cout << " b: min " << min(&b[0],n)  << " max " << max(&b[0],n)   << endl;
    //  Ici c'est beaucoup plus complique car nos avons pas l'information
    // il faut passez par les triangle ...
    // car on a renumerote le dof sommet..
  /*
    for(int i=0; i< n; ++i)

        if ( Th.v[i].OnGamma()==2)
        {
            A(i,i)=tgv;
            b[i]= tgv*2;
        }
   */
    // ok pour du P1

    for(int k=0; k< Th.nt; ++k)
        for(int ip=0; ip<3; ++ip)
          if ( Th[k][ip].OnGamma()==2)
          {
              int i=  num[k*3+ip];
              cout << i << " " << Th(k,ip) <<  " " << k << " " << ip << endl;

              A(i,i)=tgv;
              b[i]= tgv*2;
              u[i]=2;
          }
    if(debug)
    {
        A.CSR();
        M.CSR();
    cout << " M  " << A << endl;
    cout << " A  " << A << endl;
    cout << " b  " << b << endl;
    }
    double eps=1e-8;
    SparseLinearSolver<int,double> ALS(A,ts,"eps",eps,nullptr);
    ALS.solve(&u[0],&b[0]);

    ProduitMatVec(&A,&u[0],&fh[0]);
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
                f <<  (const R2 &) K[ip%3]<< " "<< u[num[3*k+(ip%3)]]<< endl;
            f << endl << endl;
        }

    }
    {
        ofstream f("uk.txt");
        for(int k=0; k<Th.nt; ++k)
        {
            for(int ip=0; ip<3; ip++)
                f << u[num[3*k+ip]] << " ";
            f << endl ;
        }

    }

    //  pour lance la verification freefem++  unix
    char ff[256];
    sprintf(ff,"FreeFem++ freefem/plotuk.edp -ns -Th %s ",argv[1]);
    system(ff);

}
