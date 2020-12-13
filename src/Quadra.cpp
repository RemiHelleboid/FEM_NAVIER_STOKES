#include "Quadra.hpp"


using namespace std;



double Quad_Triangle_Ref(Triangle &T, double (*f)(double, double)){
    double x1, y1, x2, y2, x3, y3;
    double area_K = T.area();
    T.compute_middle();
    vector V_half = T.get_vertices_middle();
    x1 = V_half[0].x();
    y1 = V_half[0].y();
    x2 = V_half[1].x();
    y2 = V_half[1].y();
    x3 = V_half[2].x();
    y3 = V_half[2].y();
    double I = (area_K / 3) * (f(x1, y1) + f(x2, y2) + f(x3, y3));
    return(I);
}
