
#include "EF2d.hpp"
#include <string>
#include <cassert>
#include <fstream>

int debug =0;

std::vector<double> DIFF_VECTOR(std::vector<double> &U, std::vector<double> &V){
    int Nu = U.size();
    int Nv = V.size();
    if (Nu != Nv){
        assert("ERROR DIFF OF TWO std::vector : SIZES DON'T MATCH");
    }
    else{
        std::vector<double> D(Nu);
        #pragma omp parallel for
            for(int i = 0; i < Nu; ++i){
                D[i] = U[i] - V[i];
            }
        return(D);
    }
}

// fonction CL
double Initial_State(const R2 & P, int label, int c)
{
    if( (label == 4)  && (c ==0) )
        return -16*(1-P.y)*(0.5-P.y); //Mesh.msh
    return 0.;
}

int Construction_Matrix_Operation(std::vector<ItemOP> &Matrix_Construction_Operation, double epsp, double nu){
    AddLap(Matrix_Construction_Operation, nu, 0, 0);
    AddLap(Matrix_Construction_Operation, nu, 1, 1);
    AddPdivV(Matrix_Construction_Operation, -1, 2, 0, 1);
    AddMass(Matrix_Construction_Operation, epsp, 2, 2);
    AdddivUQ(Matrix_Construction_Operation, -1, 0, 1, 2);
    std::cout << " Opération to construct matrix A : "<< Matrix_Construction_Operation << std::endl;
    return(1);
}

int Set_boundary_condition(Mesh2d &Mesh, std::vector<double> &Pressure, std::vector<double> &Velocity, hashMatrix &DataMatrix, const std::vector<int>* Bnum[] ,int * offset, int * BDofK){
    const double TGV = 1e20;
    int label = 0;
    for(int index_triangle=0; index_triangle < Mesh.nt; ++index_triangle){
        for(int index_vertex=0; index_vertex<3; ++index_vertex){
            label = Mesh.OnGamma(index_triangle, index_vertex);
            if ((label != Mesh2d::NotonGamma) && (label != 2)) // sur le bord
            {
                int ip0 = (index_vertex + 1) % 3;
                int ip1 = (index_vertex + 2) % 3;
                int ipa = (index_vertex + 3);

                //On  vérifie que les index sont corrects
                int i0 = Mesh(index_triangle, ip0);
                int i1 = Mesh(index_triangle, ip1);
                int ia = Mesh[index_triangle].BAdj(index_vertex);

                R2 P[3];
                P[0]=Mesh(i0);
                P[1]=Mesh(i1);
                P[2] = (P[0]+P[1])*.5;
                int NumOnEdge[6]={-1,-1,-1,-1,-1,-1};
                NumOnEdge[ip0]=0;
                NumOnEdge[ip1]=1;
                NumOnEdge[ipa]=2;
                for(int c=0; c<2; ++c) // CL sur  2 composante ...
                {
                    const  std::vector<int> & num= *Bnum[c];
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

int ExportSolution(Mesh2d &Mesh, int* offset, int* BDofK, const std::vector<int>* Bnum[], std::vector<double> &Pressure, std::vector<double> &Velocity, std::string filename){
    ofstream f(filename);
    for(int index_triangle=0; index_triangle<Mesh.nt; ++index_triangle)
    {
        for(int ip=0; ip<3; ip++){
        {
            if (ip<3)
                for(int index_field=0; index_field < 3; ++index_field){
                    int dofk = BDofK[index_field];
                    int oc = offset[index_field];
                    const  std::vector<int> & num=  *Bnum[index_field];
                    f << Velocity[num[dofk * index_triangle + ip] + oc] << " ";
                }
            }
            f << std::endl;
        }
        f << std::endl ;
    }
    std::cout<<Velocity.size()<<std::endl;
    std::cout<<Pressure.size()<<std::endl;
    return(1);
}

int ExportSolutionCSV(Mesh2d &Mesh, int* offset, int* BDofK, const std::vector<int>* Bnum[], std::vector<double> &Pressure, std::vector<double> &Velocity, std::string filename){
    ofstream f(filename);
    f<<"Index,X,Y,label,Vx,Vy,P"<<std::endl;
    for(int index_triangle=0; index_triangle<Mesh.nt; ++index_triangle)
    {
        for(int ip=0; ip<3; ip++){
        {
                f<<Mesh(index_triangle, ip)<<","<<Mesh(Mesh(index_triangle, ip))<< ",";
                for(int index_field=0; index_field < 3; ++index_field){
                    int dofk = BDofK[index_field];
                    int oc = offset[index_field];
                    const  std::vector<int> & num=  *Bnum[index_field];
                    // f << Velocity[num[dofk * index_triangle + ip] + oc] << " ";
                if (ip<3)
                     f<<Velocity[num[dofk * index_triangle + ip] + oc]<<",";
                }
            }
            f << std::endl;
            std::cout<<std::endl ;
        }
        f << std::endl ;
        std::cout<<std::endl;
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
    Mesh2d Mesh(argv[1]);
    int NbTriangle = Mesh.nt;

    std::vector<int> ListIndexP(Mesh.nt);
    std::vector<int> ListIndexU(Mesh.nt*6);
    int dimension_Uh = BuildNum(Mesh, 1, 1, 0, ListIndexU);// p2
    int dimension_Ph = BuildNum(Mesh, 1, 0, 0, ListIndexP);// p1
    int dimension_total = dimension_Uh + dimension_Uh + dimension_Ph;
    std::cout << " dimension_Ph = " << dimension_Ph << std::endl;
    std::cout << " dimension_Uh = " << dimension_Uh << std::endl;
    std::cout << " dimension_total = " << dimension_total << std::endl;
    const QF2 &Quadrature_formula = QuadratureFormular_T_5;
    std::vector<double> P1_basis_func_elem_ref(1), P2_basis_func_elem_ref(1), D(1);
    double nu = 0.0025;
    double epsp = - 1e-7;
    int Nb_P2_func_per_elem = SetWP2(Quadrature_formula, P2_basis_func_elem_ref, D);
    int Nb_P1_func_per_elem = SetWP1(Quadrature_formula, P1_basis_func_elem_ref, D);
    std::vector<ItemOP> Matrix_Construction_Operation;
    Construction_Matrix_Operation(Matrix_Construction_Operation, epsp, nu);

    const std::vector<double> * BWP[] = { &P2_basis_func_elem_ref, &P2_basis_func_elem_ref, &P1_basis_func_elem_ref};
    const std::vector<int> *Bnum[]  = {&ListIndexU, &ListIndexU, &ListIndexP};
    int offset[] = {0, dimension_Uh, dimension_Uh + dimension_Uh, dimension_Uh + dimension_Uh + dimension_Ph};
    int BDofK[] = {Nb_P2_func_per_elem, Nb_P2_func_per_elem, Nb_P1_func_per_elem};
    hashMatrix  DataMatrix(dimension_total);
    hashMatrix DataM(dimension_total);
    Add2MatLap2QF(Matrix_Construction_Operation,DataMatrix, BWP, D, Bnum, offset, Mesh);
    std::vector<double> Velocity(dimension_total);
    std::vector<double> Pressure(dimension_total);
    std::vector<double> Second_member(dimension_total);

    Set_boundary_condition(Mesh, Pressure, Velocity, DataMatrix, Bnum, offset, BDofK);

    double tolerance=1e-8;
    SparseLinearSolver<int,double> Linear_system_solver(DataMatrix, TypeOfSolver, "eps", tolerance, nullptr);
    Linear_system_solver.solve(&Velocity[0], &Pressure[0]);

    std::cout << " A: min " << min(&DataMatrix.aij[0],DataMatrix.nnz)  << " max " << max(&DataMatrix.aij[0],DataMatrix.nnz)   << std::endl;
    std::cout << " min u1:  " << min(&Velocity[0], dimension_Uh)  << " max " << max(&Velocity[0],dimension_Uh) << std::endl;
    std::cout << " min u2:  " << min(&Velocity[dimension_Uh], dimension_Uh)  << " max " << max(&Velocity[dimension_Uh],dimension_Uh)   << std::endl;
    std::cout << " min p:  " << min(&Pressure[dimension_Ph], dimension_Ph)  << " max " << max(&Pressure[dimension_Ph],dimension_Ph)   << std::endl;

    ExportSolution(Mesh, offset, BDofK, Bnum, Pressure, Velocity, "SteadyStokes.data");
    ExportSolutionCSV(Mesh, offset, BDofK, Bnum, Pressure, Velocity, "SteadyStokes.csv");

    char ff[256];
    sprintf(ff,"FreeFem++ freefem/check41.edp -ns -Mesh %s ",argv[1]);
    system(ff);

}
