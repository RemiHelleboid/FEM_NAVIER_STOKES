#include "EF2d.hpp"
#include "Walk.hpp"
#include <cassert>
#include <fstream>

int debug=0;
double time_step = 0.5;
int MaxNbEpoch = 20.;

std::vector<double> DIFF_VECTOR(std::vector<double> &U, std::vector<double> &V){
    int Nu = U.size();
    int Nv = V.size();
    if (Nu != Nv){
        assert("ERROR DIFF OF TWO VECTOR : SIZES DON'T MATCH");
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

double Entry_condition(const R2 & P, int label, int c)
{
    double value = 0.0;
    if( (label == 4)  && (c ==0) ){
        value = -16*(1-P.y)*(0.5-P.y);
    }
    return(value);
}

int Construction_Matrix_Operation(std::vector<ItemOP> &Matrix_M_Construction_Operation,std::vector<ItemOP> &Matrix_A_Construction_Operation, double RELAXATION, double nu, double time_step){
    AddLap(Matrix_A_Construction_Operation, nu , 0, 0);
    AddLap(Matrix_A_Construction_Operation, nu, 1, 1);
    AddMass(Matrix_A_Construction_Operation, 1, 0, 0);
    AddMass(Matrix_A_Construction_Operation, 1, 1, 1);
    AddMass(Matrix_A_Construction_Operation, RELAXATION, 2, 2);
    AddPdivV(Matrix_A_Construction_Operation, -time_step, 2, 0, 1);
    AdddivUQ(Matrix_A_Construction_Operation, -1, 0, 1, 2);
    AddMass(Matrix_M_Construction_Operation, 1, 0, 0);
    AddMass(Matrix_M_Construction_Operation, 1, 1, 1);
    std::cout << " Opération to construct matrix A : "<< Matrix_A_Construction_Operation << std::endl;
    std::cout << " Opération to construct matrix M : "<< Matrix_M_Construction_Operation << std::endl;
    return(1);
}

int Set_boundary_condition(Mesh2d &Mesh, std::vector<double> &Pressure, std::vector<double> &Velocity, hashMatrix &DataMatrixA, hashMatrix &DataMatrixM, const std::vector<int>* Block_to_global_index[] ,int * BlockSizes, int * Nb_basis_func_per_Triangle){
    double TGV = 1e20;
    int label = -100;
    for(int k=0; k < Mesh.nt; ++k){
        for(int e=0; e<3; ++e){
            label = Mesh.OnGamma(k,e);
            if (label != Mesh2d::NotonGamma && label!=2 ){
                int index_0 = (e+1)%3;
                int index_1 = (e+2)%3;
                int index_adj = e+3;
                int i0 = Mesh(k, index_0);
                int i1 = Mesh(k, index_1);
                int ia = Mesh[k].BAdj(e);
                R2 P[3];
                P[0] = Mesh(i0);
                P[1] = Mesh(i1);
                P[2] = (P[0]+P[1])*.5;
                int NumOnEdge[6] = {-1,-1,-1,-1,-1,-1};
                NumOnEdge[index_0] = 0;
                NumOnEdge[index_1] = 1;
                NumOnEdge[index_adj] = 2;
                for(int c=0; c<2; ++c){
                    const  std::vector<int> & Numerotation= *Block_to_global_index[c];
                    int Basis_size_triangle = Nb_basis_func_per_Triangle[c];
                    int Current_block_size = BlockSizes[c];
                    int j = 0;
                    for(int ip=0; ip<Basis_size_triangle; ++ip){
                        int j = NumOnEdge[ip];
                        if (j>=0){
                            int i =  Numerotation[k * Basis_size_triangle + ip] + Current_block_size;
                            DataMatrixA(i,i) = TGV;
                            double bound_cond = Entry_condition(P[j], label, c);
                            Velocity[i] = bound_cond;
                            Pressure[i] = bound_cond * TGV;
                        }
                    }
                }
            }
        }
    }
    return(1);
}


int ExportSolutionCSV(Mesh2d &Mesh, int* BlockSizes, int* Nb_basis_func_per_Triangle, const std::vector<int>* Block_to_global_index[], std::vector<double> &Pressure, std::vector<double> &Velocity, std::string filename){
    ofstream f(filename);
    f<<"Index,X,Y,label,Vx,Vy,P"<<std::endl;
    for(int index_triangle=0; index_triangle<Mesh.nt; ++index_triangle)
    {
        for(int ip=0; ip<3; ip++){
        {
                f<<Mesh(index_triangle, ip)<<","<<Mesh(Mesh(index_triangle, ip))<< ",";
                for(int index_field=0; index_field < 3; ++index_field){
                    int dofk = Nb_basis_func_per_Triangle[index_field];
                    int oc = BlockSizes[index_field];
                    const  std::vector<int> & num=  *Block_to_global_index[index_field];
                if (ip<3)
                     f<<Velocity[num[dofk * index_triangle + ip] + oc]<<",";
                }
            }
            f << std::endl;
        }
        f << std::endl ;
    }
    return(1);
}

double FunctionForSetF_x(double x,double y,double Velocity_x) { return x - time_step * Velocity_x;}
double FunctionForSetF_y(double x,double y,double Velocity_y) { return y - time_step * Velocity_y;}

int  main(int argc, const char** argv)
{

    assert(argc>=2);
    const char SolverSparse[]="UMFPACK", *KindSolver=SolverSparse;
    if(argc>2) KindSolver=argv[2];
    Mesh2d Mesh(argv[1]);
    int NbTriangle = Mesh.nt;

    std::vector<int> Local_to_Global_P(NbTriangle);
    std::vector<int> Local_to_Global_U(NbTriangle * 6);
    const int dimension_Uh = BuildNum(Mesh, 1, 1, 0, Local_to_Global_U);
    const int dimension_Ph = BuildNum(Mesh, 1, 0, 0, Local_to_Global_P);
    const int dimension_total = dimension_Uh + dimension_Uh + dimension_Ph;

    const QF2 &Quadrature_formula = QuadratureFormular_T_5;
    std::vector<double> P1_basis_func_elem_ref(1);
    std::vector<double> P2_basis_func_elem_ref(1);
    std::vector<double> Weigth_quadra(1);
    double NU = 1.05e-3 * time_step;    //WATER
    const double RELAXATION = - 1e-7;
    const int Nb_P2_func_per_elem = SetWP2(Quadrature_formula, P2_basis_func_elem_ref, Weigth_quadra);
    const int Nb_P1_func_per_elem = SetWP1(Quadrature_formula, P1_basis_func_elem_ref, Weigth_quadra);
    std::vector<ItemOP> Matrix_A_Construction_Operation;
    std::vector<ItemOP> Matrix_M_Construction_Operation;

    Construction_Matrix_Operation(Matrix_M_Construction_Operation, Matrix_A_Construction_Operation, RELAXATION, NU, time_step);

    const std::vector<double> * Block_to_basis_func[] = {&P2_basis_func_elem_ref, &P2_basis_func_elem_ref, &P1_basis_func_elem_ref};
    const std::vector<int> * Block_to_global_index[]  = {& Local_to_Global_U, &Local_to_Global_U, &Local_to_Global_P};
    int BlockSizes[]={0, dimension_Uh, dimension_Uh + dimension_Uh, dimension_Uh + dimension_Uh + dimension_Ph};
    int Nb_basis_func_per_Triangle[] = {Nb_P2_func_per_elem, Nb_P2_func_per_elem, Nb_P1_func_per_elem};
    hashMatrix DataMatrixA(dimension_total);
    hashMatrix DataMatrixM(dimension_total);
    Add2MatLap2QF(Matrix_A_Construction_Operation, DataMatrixA, Block_to_basis_func, Weigth_quadra, Block_to_global_index, BlockSizes, Mesh);
    Add2MatLap2QF(Matrix_M_Construction_Operation, DataMatrixM, Block_to_basis_func, Weigth_quadra, Block_to_global_index, BlockSizes, Mesh);
    std::vector<double> Velocity(dimension_total);
    std::vector<double> Pressure(dimension_total);
    std::vector<double> Second_member(dimension_total);
    std::vector<double> Velocity_x(dimension_Uh);
    std::vector<double> Velocity_y(dimension_Uh);
    std::vector<double> Prec_Velocity(dimension_total);
    std::vector<double> Error(dimension_total);
    std::vector<double> MatrixError(dimension_total);
    std::vector<double> Stokes_Velocity(dimension_total);


    for(int epoch = 1; epoch <= MaxNbEpoch; epoch++){
        std::cout << "EPOCH" << epoch << '\n';
        std::vector<double> New_Velocity(dimension_total);
        Set_boundary_condition(Mesh, Pressure, Velocity, DataMatrixA, DataMatrixM, Block_to_global_index , BlockSizes, Nb_basis_func_per_Triangle);
        double tolerance = 1e-6;
        SparseLinearSolver<int,double> Linear_system_solver(DataMatrixA, KindSolver, "eps", tolerance, nullptr);
        Linear_system_solver.solve(&Velocity[0], &Pressure[0]);
        std::vector<double> F1(Quadrature_formula.n * NbTriangle);
        std::vector<double> F2(Quadrature_formula.n * NbTriangle);
        std::vector<R2> Pp(Quadrature_formula.n * NbTriangle);
        #pragma omp parallel for
        for (int k=0; k<dimension_Uh; k++){
            Velocity_x[k] = Velocity[k];
            Velocity_y[k] = Velocity[dimension_Uh + k];
        }
        SetF(Mesh, Quadrature_formula, P2_basis_func_elem_ref, Local_to_Global_U, Velocity_x, F1, FunctionForSetF_x);
        SetF(Mesh, Quadrature_formula, P2_basis_func_elem_ref, Local_to_Global_U, Velocity_y, F2, FunctionForSetF_y);
        for(int k=0; k<F2.size(); k++){
            R2 point(F1[k],F2[k]);
            Pp[k] = point;
        }

        std::vector<double> WPfin(Quadrature_formula.n*6*3), WPKfin(Quadrature_formula.n*6*3);
        int aaa;
        aaa=SetWP2(Quadrature_formula,WPfin,Weigth_quadra);
        for(int k=0; k<NbTriangle; k++){
            const  Simplex & K= Mesh[k];
            int *jK= &Local_to_Global_U[Nb_basis_func_per_Triangle[0]*k];
            SetWK(K,WPfin,WPKfin);
            double Area_K = K.mes;
            for(int j=0; j<6; j++){
                double uxj=0;
                double uyj=0;
                for(int p=0; p<Quadrature_formula.n; p++){
                    double upx;
                    double upy;
                    int it;
                    R2 Phat;
                    const Simplex * Khat;
                    bool outside;
                    const Simplex * tstart = & Mesh[k];
                    Khat=Find(Mesh, Pp[k*Quadrature_formula.n+p],Phat,outside,tstart,it);
                    if(outside==1){
                        Pp[k*Quadrature_formula.n+p]=Mesh[it](Phat);
                    }
                    upx = CalcUp(Pp[k*Quadrature_formula.n+p],it, Velocity, *Block_to_global_index[0] ,Nb_basis_func_per_Triangle[0] , BlockSizes[0], Mesh);
                    upy = CalcUp(Pp[k*Quadrature_formula.n+p],it, Velocity, *Block_to_global_index[1] ,Nb_basis_func_per_Triangle[1] , BlockSizes[1], Mesh);
                    uxj+= Weigth_quadra[p] * upx * WPKfin[6*3*p+3*j];
                    uyj+= Weigth_quadra[p] * upy * WPKfin[6*3*p+3*j];
                }
                uxj *= Area_K;
                uyj *= Area_K;
                int jj = jK[j];
                New_Velocity[jj + BlockSizes[0]] += uxj;
                New_Velocity[jj + BlockSizes[1]] += uyj;
            }
        }

        Pressure = New_Velocity;
        Error = DIFF_VECTOR(Velocity, Prec_Velocity);
        std::cout<<Weigth_quadra.size()<<std::endl;
        std::cout<<Weigth_quadra[0]<<std::endl;
        std::cout << "TIME = " << epoch*time_step << '\n';
        double DiffIter=normsum(&Velocity[0], &Prec_Velocity[0], 2*dimension_Uh);
        std::cout << " Diff btw two itteration = "<<  " " << DiffIter  << std::endl;
        Prec_Velocity=Velocity;
        std::string FileCSV = "NavierStokesDATA/NavierStokes" + to_string(int(epoch)) + ".csv";
        ExportSolutionCSV(Mesh, BlockSizes, Nb_basis_func_per_Triangle, Block_to_global_index, Pressure, Velocity, FileCSV);

        if(epoch==MaxNbEpoch){
            {
                ofstream f("NS_1.txt");
                for(int k=0; k<NbTriangle; ++k)
                {
                    for(int c=0; c<3; ++c)
                    {
                        const  std::vector<int> & num=*Block_to_global_index[c];
                        int dofk=Nb_basis_func_per_Triangle[c];
                        int oc=BlockSizes[c];

                        for(int ip=0; ip<dofk; ip++)
                          f << Velocity[num[dofk*k+ip]+oc] << " ";
                    }
                    f << std::endl ;
                }
            }
        char ff[256];
        sprintf(ff,"FreeFem++ freefem/FreeFemNavierStokes.edp -Th %s ",argv[1]);
        system(ff);
        std::cout << "Epoch " << epoch << '\n';
        std::cout << "Diff btw two epoch= "<<  " " << DiffIter  << std::endl;
        }
    }
    return(1);
}
