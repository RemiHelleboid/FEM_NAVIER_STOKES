
#include "EF2d.hpp"
#include "Walk.hpp"
#include <cassert>
#include <fstream>
#include <omp.h>

using namespace std;
int debug =0;


vector<double> DIFF_VECTOR(vector<double> &U, vector<double> &V){
    int Nu = U.size();
    int Nv = V.size();
    if (Nu != Nv){
        assert("ERROR DIFF OF TWO VECTOR : SIZES DON'T MATCH");
    }
    else{
        vector<double> D(Nu);
        #pragma omp parallel for
            for(int i = 0; i < Nu; ++i){
                D[i] = U[i] - V[i];
            }
        return(D);
    }
}

double Entry_condition(const R2 & P, int label, int c)
{
    if( (label == 4)  && (c ==0) )
        return -16*(1-P.y)*(0.5-P.y); //Th.msh
    return 0.;
}

int Construction_Matrix_Operation(vector<ItemOP> &Matrix_M_Construction_Operation,vector<ItemOP> &Matrix_A_Construction_Operation, double epsp, double nu, double dt){
    AddLap(Matrix_A_Construction_Operation,nu,0,0);
    AddLap(Matrix_A_Construction_Operation,nu,1,1);
    AddMass(Matrix_A_Construction_Operation,epsp,2,2);
    AddPdivV(Matrix_A_Construction_Operation,-dt,2,0,1);
    AdddivUQ(Matrix_A_Construction_Operation,-1,0,1,2);
    AddMass(Matrix_M_Construction_Operation,1,0,0);
    AddMass(Matrix_M_Construction_Operation,1,1,1);
    cout << " Opération to construct matrix A : "<< Matrix_A_Construction_Operation << endl;
    cout << " Opération to construct matrix M : "<< Matrix_M_Construction_Operation << endl;
    return(1);
}

int ExportSolutionCSV(Mesh2d &Th, int* offset, int* BDofK, const vector<int>* Bnum[], vector<double> &Pressure, vector<double> &Velocity, std::string filename){
    ofstream f(filename);
    f<<"Index,X,Y,label,Vx,Vy,P"<<std::endl;
    for(int index_triangle=0; index_triangle<Th.nt; ++index_triangle)
    {
        for(int ip=0; ip<3; ip++){
        {
                f<<Th(index_triangle, ip)<<","<<Th(Th(index_triangle, ip))<< ",";
                for(int index_field=0; index_field < 3; ++index_field){
                    int dofk = BDofK[index_field];
                    int oc = offset[index_field];
                    const  vector<int> & num=  *Bnum[index_field];
                    // f << Velocity[num[dofk * index_triangle + ip] + oc] << " ";
                if (ip<3)
                     f<<Velocity[num[dofk * index_triangle + ip] + oc]<<",";
                }
            }
            f << endl;
        }
        f << endl ;
    }

    return(1);
}


int Set_boundary_condition(Mesh2d &Th, vector<double> &Pressure, vector<double> &Velocity, hashMatrix &DataMatrixA, hashMatrix &DataMatrixM, const vector<int>* Bnum[] ,int * offset, int * BDofK){
    double TGV = 1e30;
    int label = 0;
    for(int index_triangle=0; index_triangle< Th.nt; ++index_triangle)
        for(int index_vertex=0; index_vertex<3; ++index_vertex){
            label = Th.OnGamma(index_triangle, index_vertex);
            if ( label != Mesh2d::NotonGamma && label!=2 )  // sur le bord
            {
                int ip0=(index_vertex+1)%3;
                int ip1=(index_vertex+2)%3;
                int ipa=index_vertex+3;
                int i0=Th(index_triangle,ip0);
                int i1 = Th(index_triangle,ip1);
                int ia=Th[index_triangle].BAdj(index_vertex);

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

                          int i =  num[index_triangle*dofk+ip]+oc;
                          DataMatrixA(i,i) = TGV;
                          double bc = Entry_condition(P[j],label,c);
                          Pressure[i] = bc * TGV;
                          Velocity[i] = bc;
                      }
                    }
              }
          }
      }


    return(1);
}

int mcount=0;
int ccount=0;
const double dt=1.0;
const double Final_Time=10.0;
// definition of x- dt * Velocity
double FunctionForSetF_x(double x,double y,double ux) { return x;}
double FunctionForSetF_y(double x,double y,double uy) { return y;}

double mysdot(int dimension_total,const double *x,const double *y)
{
    double s=0;
    for(int i=0;i<dimension_total;++i)
        s+= x[i]*y[i];
    return  s;
}

int  main(int argc, const char** argv)
{

    assert(argc>=2);
    const char SolverSparse[]="UMFPACK", *ts=SolverSparse;
    if(argc>2) ts=argv[2];
    if(argc>3) debug =1;
    //cout << " lecture de " << argv[1] << endl;
    Mesh2d Th(argv[1]);

    vector<int> ListIndexP(Th.nt);
    vector<int> ListIndexU(Th.nt*6);
    int dimension_Uh=BuildNum(Th,1,1,0,ListIndexU);         // p2
    int dimension_Ph=BuildNum(Th,1,0,0,ListIndexP);         // p1
    int dimension_total=dimension_Uh+dimension_Uh+dimension_Ph;
    cout << " dimension_Ph = " << dimension_Ph << endl;
    cout << " dimension_Uh = " << dimension_Uh << endl;
    cout << " ndof = " << dimension_Uh+dimension_Uh+dimension_Ph << endl;
    const QF2 &Quadrature_formula=QuadratureFormular_T_5;
    vector<double> P1_basis_func_elem_ref(1),P2_basis_func_elem_ref(1),D(1);
    double nu = 0.0025*dt, epsp = - 1e-7;
    int Nb_P2_func_per_elem=SetWP2(Quadrature_formula,P2_basis_func_elem_ref,D);
    int Nb_P1_func_per_elem=SetWP1(Quadrature_formula,P1_basis_func_elem_ref,D);
    vector<ItemOP> Matrix_A_Construction_Operation;
    vector<ItemOP> Matrix_M_Construction_Operation;

    Construction_Matrix_Operation(Matrix_M_Construction_Operation, Matrix_A_Construction_Operation, epsp, nu, dt);


    //  Data block
    const vector<double> * BWP[] = { &P2_basis_func_elem_ref,&P2_basis_func_elem_ref,&P1_basis_func_elem_ref};
    const vector<int> * Bnum[]  = { & ListIndexU,&ListIndexU,&ListIndexP};
    int offset[]={ 0, dimension_Uh, dimension_Uh+dimension_Uh, dimension_Uh+dimension_Uh+dimension_Ph };
    int BDofK[]={Nb_P2_func_per_elem,Nb_P2_func_per_elem,Nb_P1_func_per_elem};
    hashMatrix  DataMatrixA(dimension_total);
    hashMatrix DataMatrixM(dimension_total);
    Add2MatLap2QF(Matrix_A_Construction_Operation,DataMatrixA,BWP,D, Bnum,offset, Th);
    Add2MatLap2QF(Matrix_M_Construction_Operation,DataMatrixM,BWP,D, Bnum,offset, Th);
    vector<double> Velocity(dimension_total);
    vector<double> Pressure(dimension_total);
    vector<double> Second_member(dimension_total);
    vector<double> Velocity_x(dimension_Uh);
    vector<double> Velocity_y(dimension_Uh);
    vector<double> Prec_Velocity(dimension_total);
    vector<double> Error(dimension_total);
    vector<double> Me(dimension_total);
    vector<double> Stokes_Velocity(dimension_total);

    //------------------------_TIME CICLE_----------------------------------
    int epoch=1;
    for(double time=0; time<Final_Time; time+=dt){
        cout << " Second_member: min " << min(&Second_member[0],dimension_total) << " max " << max(&Second_member[0],dimension_total) << endl;
        Set_boundary_condition(Th, Pressure, Velocity, DataMatrixA, DataMatrixM, Bnum, offset, BDofK);
        vector<double> New_Velocity(dimension_total);

        double eps=1e-6;
        SparseLinearSolver<int,double> Linear_system_solver(DataMatrixA,ts,"eps",eps,nullptr);

        Linear_system_solver.solve(&Velocity[0],&Pressure[0]);
        cout << " min u1:  " << min(&Velocity[0],dimension_Uh)  << " max " << max(&Velocity[0],dimension_Uh)   << endl;
        cout << " min u2:  " << min(&Velocity[dimension_Uh],dimension_Uh)  << " max " << max(&Velocity[dimension_Uh],dimension_Uh)   << endl;

      // ------------------------------- CALCULATION OF P' -----------------------------


        vector<double> F1(Quadrature_formula.n*Th.nt);
        vector<double> F2(Quadrature_formula.n*Th.nt);
        vector<R2> Pp(Quadrature_formula.n*Th.nt);
        #pragma omp parallel for
        for (int k=0; k<dimension_Uh; k++){
            Velocity_x[k]=Velocity[k];
            Velocity_y[k]=Velocity[dimension_Uh+k];
      }

      SetF(Th,Quadrature_formula,P2_basis_func_elem_ref,ListIndexU,Velocity_x,F1,FunctionForSetF_x);
      SetF(Th,Quadrature_formula,P2_basis_func_elem_ref,ListIndexU,Velocity_y,F2,FunctionForSetF_y);

      for(int k=0; k<F2.size(); k++){
          R2 point(F1[k],F2[k]);
          Pp[k]=point;}


  // ----------------------------------- FINDING OF THE TRIANGLE ----------------------------------

      // is the triangle from which the function start searching

        vector<double> WPfin(Quadrature_formula.n*6*3), WPKfin(Quadrature_formula.n*6*3);
        int aaa;
        aaa=SetWP2(Quadrature_formula,WPfin,D);


       for(int k=0; k<Th.nt; k++){

          const  Simplex & K= Th[k];

          int *jK= &ListIndexU[BDofK[0]*k];
          // int *jKy= &ListIndexU[BDofK[1]*k];

          SetWK(K,WPfin,WPKfin);
          double mes = K.mes;


          for(int j=0; j<6; j++){
            double uxj=0;
            double uyj=0;

            for(int p=0; p<Quadrature_formula.n; p++){
              double upx;
              double upy;
              int it;
              R2 Phat;                              // ??? I don't know what it is, but it is overwritten in the function FIND, so I initialize it at 0;
              const Simplex * Khat;
                                            // is the pointer at the found triangle
              bool outside;                 // tells you if P' goes to the border
              const Simplex * tstart = & Th[k];
              Khat=Find(Th,Pp[k*Quadrature_formula.n+p],Phat,outside,tstart,it);

              if(outside==1){
                Pp[k*Quadrature_formula.n+p]=Th[it](Phat);

              }

              upx = CalcUp(Pp[k*Quadrature_formula.n+p],it, Velocity, *Bnum[0] ,BDofK[0] , offset[0],Th);
              upy = CalcUp(Pp[k*Quadrature_formula.n+p],it, Velocity, *Bnum[1] ,BDofK[1] , offset[1],Th);

              uxj+= D[p] * upx * WPKfin[6*3*p+3*j];
              uyj+= D[p] * upy * WPKfin[6*3*p+3*j];

            }
            uxj*=mes;
            uyj*=mes;

            int jj=jK[j];// jjy=jKy[j];
              New_Velocity[jj+offset[0]]+=uxj;
              New_Velocity[jj+offset[1]]+=uyj;

          }

        }

    Pressure=New_Velocity;

    Error = DIFF_VECTOR(Velocity, Stokes_Velocity);



//--------  ------ L2 ERROR --------------------
   ProduitMatVec(DataMatrixM,&Error[0],&Me[0]);
  double errL2 = sqrt(mysdot(dimension_total,&Error[0],&Me[0]));

//-------------- Linf  ERROR -------------------

double errLinf=normsum(&Velocity[0],&Stokes_Velocity[0],dimension_Uh+dimension_Uh);


  cout << " err L2= "<<  " " << errL2  << endl;
  cout << " err Linf= "<<  " " << errLinf << endl;


  if(mcount==0){
    Stokes_Velocity=Velocity;
     AddMass(Matrix_A_Construction_Operation,1,0,0);
     AddMass(Matrix_A_Construction_Operation,1,1,1);
     Add2MatLap2QF(Matrix_A_Construction_Operation,DataMatrixA,BWP,D, Bnum,offset, Th);
     mcount++;
  }


////--------------------- PLOT OF DATA ------------------
if(epoch==1){

        ofstream f("u1u2p1.txt");
        for(int k=0; k<Th.nt; ++k)
        {
            for(int c=0; c<3; ++c)
            {
                const  vector<int> & num=*Bnum[c];
                int dofk=BDofK[c];
                int oc=offset[c];

                for(int ip=0; ip<dofk; ip++)
                  f << Velocity[num[dofk*k+ip]+oc] << " ";
            }
            f << endl ;
        }

    std::string FileCSV = "SteadyStokes" + to_string(time) + ".csv";
    ExportSolutionCSV(Th, offset, BDofK, Bnum, Pressure, Velocity, FileCSV);

    //  pour lance la verification freefem++  unix
    char ff[256];
    sprintf(ff,"FreeFem++ freefem/plotu1u2p1.edp -ns -Th %s ",argv[1]);
    system(ff);

}

std::cout << "Iteration = " << epoch << '\dimension_total';
std::cout << "TIME = " << second() << '\dimension_total';
epoch++;
}
}
