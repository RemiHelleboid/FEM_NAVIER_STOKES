/* ----------- README ----------------


  The modified parts are at lines:
  32-33 --> I define the function that calculates P' i.e. x - dt * u
  35 --> I define a function that puts into a vector the x and y coordinates of P'
  157 --> I coded the part about calculating the function defined in lines 26-27 starting from P
  186 --> I coded the part that finds the triangle in which P' goes

*/

#include "EF2d.hpp"
#include "Walk.hpp"
#include <cassert>
#include <fstream>

using namespace std;
int debug =0;

// fonction CL
double g(const R2 & P, int label,int c)
{
  //  if( (label == 4)  && (c ==0) ) {return -16*(1-P.y)*(0.5-P.y);}        // Condition of Th.msh
    if( (label == 4)  && (c ==0) ) return P.y*(1.-P.y)*4;                   // Condition of 1.msh


    return 0.;
}

const double dt=0.01;
// definition of x- dt * u
double FFx(double x,double y,double ux, double uy) { return x - dt * ux;}
double FFy(double x,double y,double ux, double uy) { return y - dt * uy;}

void CreatePp(vector<R2> &Pp, vector<double> Fx,vector<double> Fy){         // puts all the P' calculated in a vector
  for(int k=0; k<Fx.size(); k++){
    R2 A(Fx[k],Fy[k]);
    Pp[k]=A;}
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
    int n2=BuildNum(Th,1,1,0,numu);         // p2
    int n1=BuildNum(Th,1,0,0,nump);         // p1
    int n=n2+n2+n1;
    cout << " n1 = " << n1 << endl;
    cout << " n2 = " << n2 << endl;
    cout << " ndof = " << n2+n2+n1 << endl;
    const QF2 &qf=QuadratureFormular_T_5;
    vector<double> WPp(1),WPu(1),D(1);
    double nu = 0.01*dt, epsp = - 1e-7;
    int ndofp,ndofu;
    ndofu=SetWP2(qf,WPu,D);
    ndofp=SetWP1(qf,WPp,D);
    vector<ItemOP> OPA;
    AddLap(OPA,nu,0,0);
    AddLap(OPA,nu,1,1);
    AddMass(OPA,1,0,0);
    AddMass(OPA,1,1,1);
    AddMass(OPA,epsp,2,2);
    AddPdivV(OPA,-dt,2,0,1);
    AdddivUQ(OPA,-1,0,1,2);
    cout << " OP A ="<< OPA << endl;

    //  Data block
    const vector<double> * BWP[] = { &WPu,&WPu,&WPp};
    const vector<int> * Bnum[]  = { & numu,&numu,&nump};
    int offset[]={ 0, n2, n2+n2, n2+n2+n1 };
    int BDofK[]={ndofu,ndofu,ndofp};
    hashMatrix  DataA(n);
    Add2MatLap2QF(OPA,DataA,BWP,D, Bnum,offset, Th);
    vector<double> u(n),b(n),fh(n),ux(n2), uy(n2);






    for(int i=0; i< n; ++i)
    {

        u[i]=0;
        fh[i]=1;

    }



    double tgv = 1e30;
    cout << " fh: min " << min(&fh[0],n) << " max " << max(&fh[0],n) << endl;
    int label;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            if ( ((label=Th.OnGamma(k,e))  != Mesh2d::NotonGamma)  && (label != 2) )  // sur le bord
            {
                int ip0=(e+1)%3,ip1=(e+2)%3,ipa=e+3;                                  //  le 3 DL de l'aretes
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

    // ------------------------------- CALCULATION OF P' -----------------------------


    vector<double> F1(qf.n*Th.nt);
    vector<double> F2(qf.n*Th.nt);
    vector<R2> Pp(qf.n*Th.nt);
    for (int k=0; k<n2; k++){           // I divide u in the x and y parts defining ux and uy
      ux[k]=u[k];
      uy[k]=u[n2+k];

    }

    SetF(Th,qf,WPu,numu,ux,F1,FFx);     // calcualates the x coordinates of P'
    SetF(Th,qf,WPu,numu,uy,F2,FFy);     // calculates the y coordinates of P'

    // SETF calculates the functions F1 and F2 just in the quadrature points, so our P'
    // will be the application of F1 and F2 of the quadrature points

    CreatePp(Pp, F1, F2);

    // Prints all for verification

     std::cout << "F1 =  " << F1  << '\n' << "F2 =  " <<  F2 << endl;
     std::cout << "P' =  " << '\n' <<Pp  << endl;
  // std::cout << "u =  " << '\n'<< u<< endl;
     std::cout << "ux =  " << '\n'<< ux<< endl;
     std::cout << "uy =  " << '\n'<< uy<< endl;



// ----------------------------------- FINDING OF THE TRIANGLE ----------------------------------

    R2 Phat;                              // ??? I don't know what it is, but it is overwritten in the function FIND, so I initialize it at 0;
    const Simplex * Khat;
                                        // is the pointer at the found triangle
    bool outside;                         // tells you if P' goes to the border
    const Simplex * tstart = & Th[0];     // is the triangle from which the function start searching
    int it;

    for(int k=0; k<qf.n*Th.nt; k+=qf.n){

      Khat=Find(Th,Pp[k],Phat,outside,tstart,it);
      std::cout << "number of qp = " << qf.n <<'\n';
      std::cout << "Pp = " << Pp[k] <<'\n';
      std::cout << "K = " << it <<'\n';
      //std::cout << "Vertex  = " << (*Khatprev)[0]<< " " << (*Khatprev)[1]<< " "<< (*Khatprev)[2]<<'\n';
      std::cout << "Vertex prime = " << (*Khat)[0]<< " " << (*Khat)[1]<< " "<< (*Khat)[2]<<'\n';
      std::cout << "Phat = " << Phat <<'\n'<< "outside = " << outside <<'\n' <<endl;

// -------------------USING WP2 ----------------------
      vector<double> WPprime(qf.n*6*3);
      for(int q=0; q<qf.n; q++){

        int ndofk;
        ndofk=SetWP2(Pp[k*qf.n+q],WPprime);
          for(int pp=0; pp<ndofk; pp++){
            for(int c=0; c<3; ++c)
            {
                const  vector<int> & num=*Bnum[c];
                int dofk=BDofK[c];
                int oc=offset[c];

                for(int ip=0; ip<dofk; ip++)
                   u[num[dofk*k+ip]+oc] << " ";


          }


      }

/*      for(int c=0; c<2; ++c)
        {
            const  vector<int> & num=*Bnum[c];
            int dofk=BDofK[c];
            int oc=offset[c];

            for(int ip=0; ip<dofk; ip++)
              u[num[dofk*k+ip]+oc]=u[num[dofk*it+ip]+oc];
                    std::cout << "umentre=  " << '\n'<< u<< endl;

        }

*/
  }
  //std::cout << "udopo =  " << '\n'<< u<< endl;
  //std::cout << "num " << numu <<'\n';

    {
        ofstream f("u1u2p1.txt");
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
    sprintf(ff,"FreeFem++ freefem/plotu1u2p1.edp -ns -Th %s ",argv[1]);
    system(ff);

}
