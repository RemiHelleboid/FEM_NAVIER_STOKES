/* ----------- README ----------------


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
    if( (label == 4)  && (c ==0) ) {return -16*(1-P.y)*(0.5-P.y);}        // Condition of Th.msh
  //  if( (label == 4)  && (c ==0) ) return P.y*(1.-P.y)*4;                   // Condition of 1.msh
  //  if( (label == 2)  && (c ==0) ) {return -2*P.y*P.y+2*P.y+10;}

    return 0.;
}

// ---------------------------------------------------------
// for k...

// from p to Pp
//               // Pp[k*qf.n+p]                // we put *Bnum[0]      //BDofK[0]  // offset[0]
// double CalcUp(R2 Pp,int Kp,vector<double> u,const  vector<int> num ,int ndofu,  int c){
//
//   double  WPprime;
//   double uprime=0.0;
//   int sizeWP;
//   sizeWP=SetWP22(Pp,WPprime);
//
//   for(int i=0; i<ndofu; i++){
//
//     uprime+=u[num[ndofu*Kp+i]+c]*WPprime[3*i];
//
//
//   }
//
//   return uprime;
// }



// ---------------------------------------------------------
 int numb;
const double dt=0.1;
const double T=0.1;
// definition of x- dt * u
double FFx(double x,double y,double ux) { return ux ;}
double FFy(double x,double y,double uy) { return uy ;}

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
    double nu = 0.05*dt, epsp = - 1e-7;
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
    vector<double> u(n),b(n),fh(n),ux(n2), uy(n2), unew(n);

    //------------------------_TIME CICLE_----------------------------------

    for(double t=0; t<T; t+=dt){


      b=unew;

      double tgv = 1e30;
      cout << " fh: min " << min(&fh[0],n) << " max " << max(&fh[0],n) << endl;
      int label;
      for(int k=0; k< Th.nt; ++k)
          for(int e=0; e<3; ++e)
              if ( ((label=Th.OnGamma(k,e))  != Mesh2d::NotonGamma)  )  // sur le bord
              {
                  int ip0=(e+1)%3,ip1=(e+2)%3,ipa=e+3;                                  //  le 3 DL de l'aretes
                  int i0=Th(k,ip0), i1 = Th(k,ip1), ia=Th[k].BAdj(e); // (i0 ,i1) bord arret , ia number of the arret
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
                      for(int ip=0; ip<dofk; ++ip){
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

 


      int k1=1799;
      const Simplex &K=Th[k1];
      std::cout << "Triangle K =" <<k1 << '\n'<<endl;
      for(int p=0; p<qf.n; p++){
        //int p=0;
        R2 P=K(qf[p]);
        double ttt, tttt;

        ttt=CalcUp(P,k1, u, *Bnum[0] ,BDofK[0] , offset[0],Th);
        tttt=CalcUp(P,k1, u, *Bnum[1] ,BDofK[1] , offset[1],Th);
        std::cout << "Quadrature point number " << p<<'\n';
        std::cout << " Calcupx = " << ttt << " Calcupy = " << tttt <<'\n';
        std::cout << " F1 = " << F1[p+qf.n*k1] << " F2 = " << F2[p+qf.n*k1] <<'\n'<<endl;

      }

      CreatePp(Pp, F1, F2);

      // Prints all for verification

      //std::cout << "F1 =  " << F1  << '\n' << "F2 =  " <<  F2 << endl;
      // std::cout << "P' =  " << '\n' <<Pp[0]  << endl;
      // std::cout << "u =  " << '\n'<< u<< endl;
      // std::cout << "ux =  " << '\n'<< ux<< endl;
      // std::cout << "uy =  " << '\n'<< uy<< endl;



  // ----------------------------------- FINDING OF THE TRIANGLE ----------------------------------

      // is the triangle from which the function start searching

        vector<double> WPfin(qf.n*6*3), WPKfin(qf.n*6*3);
        int aaa;
        aaa=SetWP2(qf,WPfin,D);
        // std::cout << "Bnum[0] = " << Bnum[0] <<'\n';
        // std::cout << "BDofK[0]*k = " << BDofK[0]*k<<'\n';
        // std::cout << "BDofK[0]*k = " << BDofK[0]*k<<'\n';


       for(int k=0; k<Th.nt; k++){
          const  Simplex & K= Th[k];

          int *jK= &numu[BDofK[0]*k];

          SetWK(K,WPfin,WPKfin);
          double mes = K.mes;

          for(int j=0; j<6; j++){
            double uxj=0;

            for(int p=0; p<qf.n; p++){
              double upx;
              int it;
              R2 Phat;                              // ??? I don't know what it is, but it is overwritten in the function FIND, so I initialize it at 0;
              const Simplex * Khat;
                                            // is the pointer at the found triangle
              bool outside;                 // tells you if P' goes to the border
              const Simplex * tstart = & Th[0];
              Khat=Find(Th,Pp[k*qf.n+p],Phat,outside,tstart,it);

              if(outside==1){
                Pp[k*qf.n+p]=Th[it](Phat);
                numb++;

              }

              upx=CalcUp(Pp[k*qf.n+p],it, u, *Bnum[0] ,BDofK[0] , offset[0],Th);
              // std::cout << "number of qp = " << qf.n <<'\n';
              // std::cout << "Pp = " << Pp[k*qf.n+p] <<'\n';
              // std::cout << "K = " << it <<'\n';
              // std::cout << "Vertex prime = " << (*Khat)[0]<< " " << (*Khat)[1]<< " "<< (*Khat)[2]<<'\n';
              // std::cout << "Phat = " << Phat <<'\n'<< "outside = " << outside <<'\n' <<endl;

              uxj+= D[p] * upx * WPKfin[6*3*p+3*j];
            }
            uxj*=mes;
            int jj=jK[j];
            //take care about boundary conditions
            if(1)
            unew[jj+offset[0]]+=uxj;

          }

        }

  //--------------------Y coordinates of U ------------------------

  for(int k=0; k<Th.nt; k++){
     const  Simplex & K= Th[k];

     int *jK= &numu[BDofK[1]*k];


     SetWK(K,WPfin,WPKfin);
     double mes = K.mes;

     for(int j=0; j<6; j++){
       double uyj=0;

       for(int p=0; p<qf.n; p++){
         double upy;
         int it;
         R2 Phat;                              // ??? I don't know what it is, but it is overwritten in the function FIND, so I initialize it at 0;
         const Simplex * Khat;
                                       // is the pointer at the found triangle
         bool outside;                 // tells you if P' goes to the border
         const Simplex * tstart = & Th[0];
         Khat=Find(Th,Pp[k*qf.n+p],Phat,outside,tstart,it);

         if(outside==1){
             Pp[k*qf.n+p]=Th[it](Phat);
             numb++;
           }
         upy=CalcUp(Pp[k*qf.n+p],it, u, *Bnum[1] ,BDofK[1] , offset[1],Th);

         uyj+= D[p] * upy * WPKfin[6*3*p+3*j];
         //std::cout << "line = "<< 6*3*p+3*j<< " "<< "WPKfin = " << WPKfin[6*3*p+3*j] <<'\n';

       }
       uyj*=mes;
       int jj=jK[j];
       //take care about boundary conditions
       if(1)
       unew[jj+offset[1]]+=uyj;

     }
   }



 //       std::cout << "count = " << numb <<'\n';

/*

////--------------------- PLOT OF DATA ------------------
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
    sprintf(ff,"FreeFem++ plotu1u2p1.edp -ns -Th %s ",argv[1]);
    system(ff);
*/
}
}
