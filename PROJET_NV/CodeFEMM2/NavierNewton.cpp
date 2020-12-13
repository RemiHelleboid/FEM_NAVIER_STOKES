#include "EF2d.hpp"
#include <string>
#include <cassert>
#include <fstream>

using namespace std;



#include "EF2d.hpp"
#include <string>
#include <cassert>
#include <fstream>

using namespace std;
int debug =0;

// fonction CL
double Initial_State(const R2 & P, int label, int c)
{
    if( (label == 4)  && (c ==0) )
        return -16*(1-P.y)*(0.5-P.y); //Th.msh
    return 0.;
}

int Construction_Matrix_Operation(vector<ItemOP> &Matrix_Construction_Operation, double epsp, double nu){
    AddLap(Matrix_Construction_Operation, nu, 1, 1);
    AddPdivV(Matrix_Construction_Operation, -1, 2, 0, 1);
    AddMass(Matrix_Construction_Operation, epsp, 2, 2);
    AdddivUQ(Matrix_Construction_Operation, -1, 0, 1, 2);
    AddLap(Matrix_Construction_Operation, nu, 0, 0);
    cout << " Opération to construct matrix A : "<< Matrix_Construction_Operation << endl;
    return(1);
}

int Set_boundary_condition(Mesh2d &Th, vector<double> &Pressure, vector<double> &Velocity, hashMatrix &DataMatrix, const vector<int>* Bnum[] ,int * offset, int * BDofK){
    const double TGV = 1e20;
    int label = 0;
    for(int index_triangle=0; index_triangle < Th.nt; ++index_triangle){
        for(int index_vertex=0; index_vertex<3; ++index_vertex){
            label = Th.OnGamma(index_triangle, index_vertex);
            if ((label != Mesh2d::NotonGamma) && (label != 2)) // sur le bord
            {
                int ip0 = (index_vertex + 1) % 3;
                int ip1 = (index_vertex + 2) % 3;
                int ipa = (index_vertex + 3);

                //On  vérifie que les index sont corrects
                int i0 = Th(index_triangle, ip0);
                int i1 = Th(index_triangle, ip1);
                int ia = Th[index_triangle].BAdj(index_vertex);

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
                        if( (j=NumOnEdge[ip])>=0){
                            int i =  num[index_triangle*dofk+ip]+oc;
                            DataMatrix(i,i) = TGV;
                            double bc = Initial_State(P[j], label,c);
                            Pressure[i] = bc * TGV;
                            Velocity[i] = bc;
                        }
                    }
                }
            }
        }
    }
    return(1);
}

int ExportSolution(Mesh2d &Th, int* offset, int* BDofK, const vector<int>* Bnum[], vector<double> &Pressure, vector<double> &Velocity, std::string filename){
    ofstream f(filename);
    for(int index_triangle=0; index_triangle<Th.nt; ++index_triangle)
    {
        for(int ip=0; ip<3; ip++){
        {
                for(int index_field=0; index_field < 3; ++index_field){
                    int dofk = BDofK[index_field];
                    int oc = offset[index_field];
                    const  vector<int> & num=  *Bnum[index_field];
                    f << Velocity[num[dofk * index_triangle + ip] + oc] << " ";
                if (ip<3)
                    std::cout<<"Vertex : "<<Th(index_triangle, ip)<<"   "<<Th(Th(index_triangle, ip))<< "   "<< Velocity[num[dofk * index_triangle + ip] + oc]<<std::endl;
                }
            }
            f << endl;
        }
        f << endl ;
    }
    std::cout<<Velocity.size()<<std::endl;
    std::cout<<Pressure.size()<<std::endl;
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
            std::cout<<endl ;
        }
        f << endl ;
        std::cout<<endl;
    }
    std::cout<<Velocity.size()<<std::endl;
    std::cout<<Pressure.size()<<std::endl;
    return(1);
}






int  main(int argc, const char** argv)
{
    assert(argc>=2);
    const char SolverSparse[]="UMFPACK", *TypeOfSolver=SolverSparse;
    if(argc>2) TypeOfSolver=argv[2];
    if(argc>3) debug =1;
    cout << " lecture de " << argv[1] << endl;
    Mesh2d Th(argv[1]);

    vector<int> ListIndexP(Th.nt);
    vector<int> ListIndexU(Th.nt*6);
    int dimension_Uh = BuildNum(Th, 1, 1, 0, ListIndexU);// p2
    int dimension_Ph = BuildNum(Th, 1, 0, 0, ListIndexP);// p1
    int dimension_total = dimension_Uh + dimension_Uh + dimension_Ph;
    cout << " dimension_Ph = " << dimension_Ph << endl;
    cout << " dimension_Uh = " << dimension_Uh << endl;
    cout << " dimension_total = " << dimension_total << endl;
    const QF2 &Quadrature_formula = QuadratureFormular_T_5;
    vector<double> P1_basis_func_elem_ref(1), P2_basis_func_elem_ref(1), D(1);
    double nu = 1.0e-3;
    double epsp = - 1e-7;
    int Nb_P2_func_per_elem = SetWP2(Quadrature_formula, P2_basis_func_elem_ref, D);
    int Nb_P1_func_per_elem = SetWP1(Quadrature_formula, P1_basis_func_elem_ref, D);
    vector<ItemOP> Matrix_Construction_Operation;
    Construction_Matrix_Operation(Matrix_Construction_Operation, epsp, nu);

    const vector<double> * BWP[] = { &P2_basis_func_elem_ref, &P2_basis_func_elem_ref, &P1_basis_func_elem_ref};
    const vector<int> *Bnum[]  = {&ListIndexU, &ListIndexU, &ListIndexP};
    int offset[] = {0, dimension_Uh, dimension_Uh + dimension_Uh, dimension_Uh + dimension_Uh + dimension_Ph};
    int BDofK[] = {Nb_P2_func_per_elem, Nb_P2_func_per_elem, Nb_P1_func_per_elem};
    hashMatrix  DataMatrix(dimension_total);
    hashMatrix DataM(dimension_total);
    Add2MatLap2QF(Matrix_Construction_Operation,DataMatrix, BWP, D, Bnum, offset, Th);
    vector<double> Velocity(dimension_total);
    vector<double> Pressure(dimension_total);
    vector<double> Second_member(dimension_total);

    Set_boundary_condition(Th, Pressure, Velocity, DataMatrix, Bnum, offset, BDofK);

    double tolerance=1e-6;
    SparseLinearSolver<int,double> Linear_system_solver(DataMatrix, TypeOfSolver, "eps", tolerance, nullptr);
    Linear_system_solver.solve(&Velocity[0], &Pressure[0]);

    cout << " A: min " << min(&DataMatrix.aij[0],DataMatrix.nnz)  << " max " << max(&DataMatrix.aij[0],DataMatrix.nnz)   << endl;
    cout << " min u1:  " << min(&Velocity[0], dimension_Uh)  << " max " << max(&Velocity[0],dimension_Uh) << endl;
    cout << " min u2:  " << min(&Velocity[dimension_Uh], dimension_Uh)  << " max " << max(&Velocity[dimension_Uh],dimension_Uh)   << endl;
    cout << " min p:  " << min(&Pressure[dimension_Ph], dimension_Ph)  << " max " << max(&Pressure[dimension_Ph],dimension_Ph)   << endl;

    ExportSolution(Th, offset, BDofK, Bnum, Pressure, Velocity, "SteadyStokes.data");
    ExportSolutionCSV(Th, offset, BDofK, Bnum, Pressure, Velocity, "SteadyStokes.csv");

    char ff[256];
    sprintf(ff,"FreeFem++ freefem/check41.edp -ns -Th %s ",argv[1]);
    system(ff);

}
