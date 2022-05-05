
#include <CGAL/Epick_d.h>
//#include <CGAL/point_generators_d.h>
//#include <CGAL/predicates_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/property_map.h>
#include <CGAL/algorithm.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/assertions.h>

#include <Eigen/Dense>

#include "CLI11.hpp"

#include <random>
#include <iostream>
#include <iterator>
#include <vector>

#include <math.h>       /* acos */

#include "stetvol.h"

typedef CGAL::Epick_d< CGAL::Dimension_tag<4> >         K;
typedef CGAL::Delaunay_triangulation<K>                  T;
typedef K::Point_d                                      Point;
typedef CGAL::Spatial_sort_traits_adapter_d<K,
  CGAL::First_of_pair_property_map<std::pair<Point,int>>
> Search_traits;

Eigen::Matrix<double,Eigen::Dynamic,4> V,CC;
Eigen::Matrix<int,Eigen::Dynamic,4> F,N;
std::vector<bool> degtet;
std::vector<int> po2vo,vo2po;


const double s3a = 2.0*M_PI*M_PI;
const double offset = 0.5;

constexpr std::array<std::array<int, 4>, 6> eo{
    0, 1, 2, 3,
    0, 2, 1, 3,
    0, 3, 1, 2,
    1, 2, 0, 3,
    1, 3, 0, 2,
    2, 3, 0, 1
};

// edge connectivity
struct edge
{
    std::array<int, 2> v;
    bool operator==(edge const &o) const
    {
        return ((v[0] == o.v[0]) || (v[0] == o.v[1])) && ((v[1] == o.v[0]) || (v[1] == o.v[1]));
    }
    int &operator[](int idx)
    {
        return v[idx];
    }
    const int &operator[](int idx) const
    {
        return v[idx];
    }
};

template <>
struct std::hash<edge>
{
    size_t operator()(edge const &key) const
    {
        return (size_t)key[0] * (size_t)key[1];
    }
};

double root(int k, int l, int c)
{
//    std::cerr << "poly: x^"<< k << " - " << l << " x - " << c << " = 0" << std::endl;;
    double p = 1.0;
    double s = (std::pow(p,k)-(double)l*p-(double)c) /
        ((double)k*std::pow(p,k-1)-(double)l);
    while (std::fabs(s) > 1e-15)
    {
        p = p - s;
        s = (std::pow(p,k)-(double)l*p-(double)c) / ((double)k*std::pow(p,k-1)-(double)l);
        //std::cerr << p << ", " << r << ", " << d << std::endl;

    }
//    std::cerr << "x = " << p << std::endl;
    return p;
}

double rootp(std::vector<int> coef)
{
    std::cerr << "poly: " << coef[0];
    for (int i = 1; i < coef.size(); i++)
    std::cerr << " + " << coef[i] << " x^" << i;
    std::cerr << std::endl;;
    
    double x = 1.0;
    double f = (double)coef[0], df = 0.0;
    for (int i = 1; i < coef.size(); i++)
    {
        f += (double)coef[i] * std::pow(x,i);
        df += (double)i * (double)coef[i] * std::pow(x,i-1);
    }
    double s = f/df;
    int r = 0;
    while (std::fabs(s) > 1e-15 && std::fabs(x) < 1e3 && r < 20)
    {
        x = x - s;
        f = (double)coef[0]; df = 0.0;
        for (int i = 1; i < coef.size(); i++)
        {
            f += (double)coef[i] * std::pow(x,i);
            df += (double)i * (double)coef[i] * std::pow(x,i-1);
        }
        s = f/df;
        r++;
    }
    std::cerr << "x = " << x;
    if (r >= 20) std::cerr << " (not converged)";
    std::cerr << std::endl;
    return x;
}



Eigen::Matrix<double,Eigen::Dynamic,4> double_fib(int n, double mc0, double mc1)
{
//    std::cerr << "Creating " << n << " points by the fibonacci method." << std::endl;
//    auto tp0 = std::chrono::steady_clock::now();
    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    for (int i = 0; i < n; i++)
    {
//        double s = ((double)i+offset)/((double)n-1.0+2.0*offset);
        double s = ((double)i+offset)/((double)n);
        // create x and y based on fib spiral in the plane
        double ab = 2.0 * M_PI * s * (double)n;
        double r = sqrt(s);
        double theta = ab * mc0;
        double R = sqrt(1.0-s);
        double phi = ab * mc1;
        V.row(i) << r*sin(theta),r*cos(theta),R*sin(phi),R*cos(phi);
    }
//    auto tp1 = std::chrono::steady_clock::now();
//    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    return V;
}

void create_delaunay()
{
    int np = V.rows();
    std::vector<Point> points(np);
    for (int i = 0; i < np; i++)
        points[i] = Point(V(i,0),V(i,1),V(i,2),V(i,3));
    
//    t.insert(Point(CGAL::Origin()));
//    t.insert(points.begin(), points.end());
    
    T t(4);
    T::Full_cell_handle hint;
    hint = t.insert(Point(CGAL::Origin()))->full_cell();
    for (int i = 0; i < np; i++) hint = t.insert(points[po2vo[i]],hint)->full_cell();
    
    CGAL_assertion( t.is_valid() );

//    std::cerr << "Triangulation successfully computed: "
//    << t.number_of_vertices() << " vertices, "
//    << t.number_of_finite_full_cells() << " finite cells.\n Establishing tetrahedralization... "
//    << std::endl;

//    V = Eigen::Matrix<double,Eigen::Dynamic,4>(np,4);
    CGAL::Unique_hash_map<T::Vertex_handle, int> VH;
    T::Vertex_iterator vit = t.vertices_begin();
    int i = 0;
    for (; vit != t.vertices_end(); ++vit)
        if (!t.is_infinite(vit) &&
            vit->point() != Point(CGAL::Origin()) )
        {
            VH[vit] = po2vo[i];
//            VH[vit] = i;
//            for (int j = 0; j < 4; j++)
//                V(i,j) = vit->point()[j];
            i++;
        }
    assert(np==i);

    std::vector<T::Full_cell_handle> surface_cells;
    std::back_insert_iterator<std::vector<T::Full_cell_handle>> out(surface_cells);
    t.incident_full_cells(t.infinite_vertex(), out);
    
    int nf = surface_cells.size();
    CGAL::Unique_hash_map<T::Full_cell_handle, int> FH;
    i = 0;
    for (auto sit = surface_cells.begin();
         sit != surface_cells.end(); ++sit)
    {
        FH[*sit] = i;
        i++;
    }
    
    F = Eigen::Matrix<int,Eigen::Dynamic,4>(nf,4);
    N = Eigen::Matrix<int,Eigen::Dynamic,4>(nf,4);
    i = 0;
    for (auto sit = surface_cells.begin();
         sit != surface_cells.end(); ++sit)
    {
        int k = 0, ii = -1;
        for (int j = 0; j < 5; j++)
        {
            if ((*sit)->vertex(j) != t.infinite_vertex())
            {
                F(i,k) = VH[(*sit)->vertex(j)];
                N(i,k) = FH[(*sit)->neighbor(j)];
                k++;
            }
            else ii = j;
        }
        assert(k==4 && ii >= 0);
        if (ii % 2 == 0)
        {
            std::swap(F(i,0),F(i,1));
            std::swap(N(i,0),N(i,1));
        }
        i++;
    }
//    std::cerr << "done.\nThere are " << nf << " facets on the convex hull."<< std::endl;
}

Eigen::Vector4d areavector(Eigen::Matrix4d P)
{
    Eigen::Matrix4d Q;
    Q.row(0) = P.row(0);
    for (int j = 1; j < 4; j++)
        Q.row(j) = P.row(j)-P.row(0);
    Eigen::FullPivLU<Eigen::Matrix4d> lu(Q);
    double d = fabs(lu.determinant());
    if (d < 1e-15) return Eigen::RowVector4d::Zero();
    Eigen::Vector4d a = lu.solve(Eigen::Vector4d(1.0,0.0,0.0,0.0));
    // ||a|| * area = d (up to prop constants)
    // ||a|| should be area
    double as = a.squaredNorm();
    if (as < 1e-15) return Eigen::RowVector4d::Zero();
    return (d/as) * a;
}

void analyze_tets()
{
    int nd = 0;
    const Eigen::Vector4d b(1.0,0.0,0.0,0.0);
    CC = Eigen::Matrix<double,Eigen::Dynamic,4>(F.rows(),4);
    degtet.resize(F.rows());
    for (int f = 0; f < F.rows(); f++)
    {
        Eigen::Matrix4d P;
        for (int j = 0; j < 4; j++)
            P.row(j) = V.row(F(f,j));
        for (int j = 1; j < 4; j++)
            P.row(j) -= P.row(0);
        Eigen::FullPivLU<Eigen::Matrix4d> lu(P);
        degtet[f] = (lu.rank() < 4);
        if (degtet[f])
            nd++;
        else
            CC.row(f) = (lu.solve(b)).normalized();
    }
}

double vb(double a)
{
    double v = vol(Eigen::VectorXd::Constant(6,a));
    double c = 2.0*M_PI/a;
    double t = s3a/v;
    //    std::cerr << "vb( " << a << " ): " << t*(6.0/c-1.0) << std::endl;
    return t*(6.0/c-1.0);
}


double ab(int v)
{
    double al = 2.0*M_PI/4.0;
    double am = std::acos(1.0/3.0);
    double au = 0.5*(al+am);
    double va = vb(au);
    while ( va > 0.0 && va < (double)v )
    {
        al = au;
        au = 0.5*(al+am);
        va = vb(au);
    }
    am = 0.5*(al+au);
    va = vb(am);
    while (al-au > 1e-15 && std::fabs((double)v-va) >= 1e-12)
    {
        //        std::cerr << std::setprecision (15) << (al-au) << ", " << fabs((double)v-va) << std::endl;
        if (va > (double)v)
            au = am;
        else
            al = am;
        am = 0.5*(al+au);
        va = vb(am);
    }
    return am;
}

double radius(double a)
{
    Eigen::Matrix<double,3,4> R;
    R.col(0) <<  1, 0, -1/sqrt(2.0);
    R.col(1) << -1, 0, -1/sqrt(2.0);
    R.col(2) << 0,  1,  1/sqrt(2.0);
    R.col(3) << 0, -1,  1/sqrt(2.0);
    Eigen::Matrix4d P;
    
    double s = 0.5;
    double h = sqrt(1.0-1.5*s*s);
    P.topRows<3>() = s*R;
    P.row(3) = Eigen::RowVector4d::Constant(h);
    
    while (dihedrals(P)[0] > a)
    {
        s *= 0.5;
        h = sqrt(1.0-1.5*s*s);
        P.topRows<3>() = s*R;
        P.row(3) = Eigen::RowVector4d::Constant(h);
    }
    double su = 2.0*s, sl = s;
    s = 0.5*(sl+su);
    h = sqrt(1.0-1.5*s*s);
    P.topRows<3>() = s*R;
    P.row(3) = Eigen::RowVector4d::Constant(h);
    double di = dihedrals(P)[0];
    while (su-sl > 1e-15 && fabs(di-a) > 1e-15)
    {
        if (di > a)
            su = s;
        else
            sl = s;
        s = 0.5*(sl+su);
        h = sqrt(1.0-1.5*s*s);
        P.topRows<3>() = s*R;
        P.row(3) = Eigen::RowVector4d::Constant(h);
        di = dihedrals(P)[0];
    }
    return std::acos(h);
}





int main(int argc, char *argv[])
{
    CLI::App app{"Exploring irrational numbers for double fib sampling"};

    int np = 64;
    app.add_option("-n", np, "number of points to start with");

    int npb = 1000000;
    app.add_option("-b", npb, "bound for number of points");

    double magic_ratio = 1.43;
    app.add_option("-m", magic_ratio, "ratio bound for dispersion");

    
    bool rot = false;
    app.add_flag("-r,--rotation,--so3", rot, "treat points as rotations");
    
    CLI11_PARSE(app, argc, argv);


    std::vector<std::vector<int>> param;
    std::vector<double> psi;
    
    int nc = 0;
    std::vector<int> c(8,0);
    for (c[0] = -4; c[0] <= -1; c[0]++)
    for (c[1] = -1; c[1] <= 1; c[1]++)
    for (c[2] = -1; c[2] <= 1; c[2]++)
    for (c[3] = -1; c[3] <= 1; c[3]++)
    for (c[4] = -1; c[4] <= 1; c[4]++)
    for (c[5] = -1; c[5] <= 1; c[5]++)
    for (c[6] = -1; c[6] <= 1; c[6]++)
    for (c[7] = 0; c[7] <= 1; c[7]++)
    {
        int cs = 0;
        for (int i = 0; i < 6; i++) if (c[i] != 0) cs++;
        if (cs > 5) continue;
        double phi = rootp(c);
        if (phi > 1.0 && phi < 2.0)
        {
                param.push_back(c);
                psi.push_back(phi);
                nc++;
        }
    }
    std::cout << "Number of candidates: " << nc << std::endl;
    
//    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(nc,nc);

    struct cp
    {
        int idx0, idx1, exp0, exp1;
        std::vector<double> value;
    };
//    std::vector<std::pair<std::pair<std::pair<int,int>,int>,std::vector<double>>> pv(nc*(nc-1)/2);
    std::vector<cp> pv;
    for (int ri = 0; ri < nc-1; ri++)
        for (int rj = ri+1; rj < nc; rj++)
        {
            cp p0 = {ri,rj,-1,-1,std::vector<double>()};
            cp p1 = {ri,rj,1,1,std::vector<double>()};
            pv.push_back(p0);
            pv.push_back(p1);
        }
    for (int ri = 0; ri < nc-1; ri++)
        {
            cp p0 = {ri,ri,-1,-2,std::vector<double>()};
            cp p1 = {ri,ri,1,2,std::vector<double>()};
            pv.push_back(p0);
            pv.push_back(p1);
        }

    int pc = pv.size();

            
    std::cerr << "Starting with " << pc << " pairs." << std::endl;
    const Eigen::Vector4d b(1.0,0.0,0.0,0.0);

    for (; np < npb; np *= 2 )
    {
        std::cout << "n: " << np << std::endl;
        std::vector<Point> points(2*np);
        double disp_bound = radius(ab(2*np));
        
        for (int p = 0; p < pc; p++)
        {
            V = double_fib(np,
                           std::pow(psi[pv[p].idx0],pv[p].exp0),
                           std::pow(psi[pv[p].idx1],pv[p].exp1));
        
        int pidx = 0;
        for (int i = 0; i < np; i++)
        {
            points[pidx] = Point(V(i,0),V(i,1),V(i,2),V(i,3));
            pidx++;
            points[pidx] = Point(-V(i,0),-V(i,1),-V(i,2),-V(i,3));
            pidx++;
        }
        
        T t(4);
        t.insert(Point(CGAL::Origin()))->full_cell();
        t.insert(points.begin(),points.end());
        
            std::vector<T::Full_cell_handle> surface_cells;
            std::back_insert_iterator<std::vector<T::Full_cell_handle>> out(surface_cells);
            t.incident_full_cells(t.infinite_vertex(), out);

            double mr = 0.0;
            for (auto sit = surface_cells.begin();
                 sit != surface_cells.end(); ++sit)
            {
                Eigen::Matrix4d P;
                int k = 0;
                for (int j = 0; j < 5; j++)
                {
                    if ((*sit)->vertex(j) != t.infinite_vertex())
                    {
                        // P.row(k) = (*sit)->point();
                        for (int e = 0; e < 4; e++)
                            P(k,e) = (*sit)->vertex(j)->point()[e];
                        k++;
                    }
                }
                for (int j = 1; j < 4; j++)
                    P.row(j) -= P.row(0);
                Eigen::FullPivLU<Eigen::Matrix4d> lu(P);
                if (lu.rank() >= 4)
                {
                    Eigen::Vector4d cc = (lu.solve(b)).normalized();
                    double r = acos(fabs(cc.dot(P.row(0))));
                    if (r > mr) mr = r;
                }
            }
            
            pv[p].value.push_back(mr/disp_bound);
       
        }
        
//        double br = pv[0].second.back();
//        for (int i = 1; i < pc; i++)
//            if (pv[i].second.back() < br) br = pv[i].second.back();
//        br *= 1.2;
        
        int pcr = 0;
        for (int i = 0; i < pc; i++)
            if (pv[i].value.back() < magic_ratio)
            {
                pv[pcr] = pv[i];
                pcr++;
            }
        pc = pcr;
        std::cout << "Retaining " << pc << " pairs."  << std::endl;

    }
    
    std::sort(pv.begin(),pv.begin()+pc,[](auto &left, auto &right) {
            return std::accumulate(left.value.begin(),
                                   left.value.end(),0.0) < std::accumulate(right.value.begin(),
                                                right.value.end(),0.0);
    });

    std::cerr << "Top 100 sums:" << std::endl;
    for (int i = 0; i < std::min(100,pc); i++)
    {
        std::cout << std::accumulate(pv[i].value.begin(),pv[i].value.end(),0.0) << " " << psi[pv[i].idx0] << " " << pv[i].exp0 << " " << std::pow(psi[pv[i].idx0],pv[i].exp0) << " " << psi[pv[i].idx1] << " " << pv[i].exp1 << " " << std::pow(psi[pv[i].idx1],pv[i].exp1)
        << std::endl;
        for (auto ce : param[pv[i].idx0]) std::cerr << ce << " ";
        std::cerr << " exp: " << pv[i].exp0 << std::endl;
        for (auto ce : param[pv[i].idx1]) std::cerr << ce << " ";
        std::cerr << " exp: " << pv[i].exp1 << std::endl;
    }

    return 0;
}
