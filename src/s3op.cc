
#include <CGAL/Epick_d.h>
//#include <CGAL/point_generators_d.h>
//#include <CGAL/predicates_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_d.h>
#include <CGAL/property_map.h>
#include <CGAL/algorithm.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Timer.h>


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


Eigen::Matrix<double,Eigen::Dynamic,4> V,AV;
Eigen::Matrix<int,Eigen::Dynamic,4> F,N;
std::vector<int> po2vo,vo2po;

const double s3a = 2.0*M_PI*M_PI;


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



bool read_points(std::basic_istream<char>& vf)
{
    int n;
    vf >> n;
    std::cerr << "Reading " << n << " points... " << std::flush;
    V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    for (int i = 0; i < n; i++)
    {
        vf >> V(i,0) >> V(i,1) >> V(i,2) >> V(i,3);
        if ( fabs(1.0-V.row(i).squaredNorm()) > 1e-12 )
        {
            V.row(i).normalize();
            //std::cerr << "Point in row " << i << " not unit distance: " << V.row(i) << "\nDistance: " << std::setprecision(16) << V.row(i).squaredNorm() << std::endl;
        }
    }
    std::cerr << "done!" << std::endl;
    return true;
}

void create_delaunay()
{
    int np = V.rows();
    CGAL::Timer cost;
    std::vector<Point> points(np);

    
    for (int i = 0; i < np; i++)
        points[i] = Point(V(i,0),V(i,1),V(i,2),V(i,3));
    
    T t(4);
    std::cerr << "Inserting... " << std::flush;
    T::Full_cell_handle hint;
    cost.reset(); cost.start();
    hint = t.insert(Point(CGAL::Origin()))->full_cell();
    for (int i = 0; i < np; i++) hint = t.insert(points[po2vo[i]],hint)->full_cell();
    cost.stop(); std::cerr << "done in "<<cost.time()<<" seconds." << std::endl;
    
//    T t(4);
//    std::cerr << "Delaunay with standard inserter... " << std::flush;
//    cost.reset(); cost.start();
//    t.insert(Point(CGAL::Origin()));
//    t.insert(points.begin(), points.end());
//    cost.stop(); std::cerr << "done in "<<cost.time()<<" seconds." << std::endl;
//    CGAL_assertion( t.is_valid() );
    

    std::cerr << "Triangulation successfully computed: "
    << t.number_of_vertices() << " vertices, "
    << t.number_of_finite_full_cells() << " finite cells.\n Establishing tetrahedralization... "
    << std::endl;

    CGAL::Unique_hash_map<T::Vertex_handle, int> VH;
    T::Vertex_iterator vit = t.vertices_begin();
    int i = 0;
    for (; vit != t.vertices_end(); ++vit)
        if (!t.is_infinite(vit) &&
            vit->point() != Point(CGAL::Origin()) )
        {
            VH[vit] = po2vo[i];
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
    std::cerr << "done.\nThere are " << nf << " facets on the convex hull."<< std::endl;
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

double compute_AV()
{
    int nd = 0;
    const Eigen::Vector4d b(1.0,0.0,0.0,0.0);
    AV = Eigen::MatrixXd::Zero(F.rows(),4);
    Eigen::VectorXd cr = Eigen::VectorXd::Ones(F.rows());
    
    for (int f = 0; f < F.rows(); f++)
    {
        Eigen::Matrix4d P;
        for (int j = 0; j < 4; j++)
            P.row(j) = V.row(F(f,j));
        for (int j = 1; j < 4; j++)
            P.row(j) -= P.row(0);
        
        Eigen::FullPivLU<Eigen::Matrix4d> lu(P);
        double d = fabs(lu.determinant());
        if (d < 1e-15) {AV.row(f) = Eigen::RowVector4d::Zero(); continue;}
        
        Eigen::Vector4d a = lu.solve(Eigen::Vector4d(1.0,0.0,0.0,0.0));
        // a, for now, is on the hyperplane through the tet.
        // use this to compute smallest circle.

        double an = a.norm();
        if (an < 1e-12) continue;

        a /= an;
               
        cr(f) = a.dot(P.row(0));
        
        // ||a|| * area = d (up to prop constants)
        // ||a|| should be area
 
        // scale by 'area'
        //a *= d/an;
        
        // scale by circumradius
        // a *= std::pow(std::acos(cr(f)),4);
        
//        a *= d*d/an;

        AV.row(f) = a;
    }
    double disp = cr.minCoeff();
    std::cerr << "Dispersion: " << acos(disp)*180.0/M_PI << std::endl;
    return disp;
}

double shortest_edge()
{
    // per edge measures: edge lengths
    std::unordered_map<edge,int> eidx;
    int ne = 0;
    for (int f = 0; f < F.rows(); f++)
        for (int j = 0; j < 6; j++)
        {
            edge e = edge{F(f,eo[j][0]),F(f,eo[j][1])};
            auto ei = eidx.find(e);
            if (ei == eidx.end())
                eidx[e] = ne++;
        }
    Eigen::VectorXd el(ne);
    for (auto e: eidx)
        el[e.second] = acos(V.row(e.first[0]).dot(V.row(e.first[1])));
    std::cerr << "Shortest edge: " << el.minCoeff() << std::endl;
    return el.minCoeff();
}

double vol()
{
    double vol = 0.0;
    for (int f = 0; f < F.rows(); f++)
    {
        Eigen::Matrix4d P;
        for (int j = 0; j < 4; j++)
            P.row(j) = V.row(F(f,j));
        double d = P.determinant();
        if (d < 0.0) return 1e12;
        vol += P.determinant();
    }
    std::cerr << "Volume: " << (vol/(12.0*M_PI*M_PI)) << std::endl;
    return vol;
}


double disp()
{
    int nd = 0;
    const Eigen::Vector4d b(1.0,0.0,0.0,0.0);
    AV = Eigen::MatrixXd::Zero(F.rows(),4);
    Eigen::VectorXd cr = Eigen::VectorXd::Ones(F.rows());
    
    for (int f = 0; f < F.rows(); f++)
    {
        Eigen::Matrix4d P;
        for (int j = 0; j < 4; j++)
            P.row(j) = V.row(F(f,j));
        for (int j = 1; j < 4; j++)
            P.row(j) -= P.row(0);
        
        Eigen::FullPivLU<Eigen::Matrix4d> lu(P);
        double d = fabs(lu.determinant());
        if (d < 1e-15) {AV.row(f) = Eigen::RowVector4d::Zero(); continue;}
        
        Eigen::Vector4d a = lu.solve(Eigen::Vector4d(1.0,0.0,0.0,0.0));
        
        double an = a.norm();
        if (an < 1e-12) continue;
        
        a /= an;
        
        cr(f) = a.dot(P.row(0));
    }
    double disp = cr.minCoeff();
    std::cerr << "Dispersion: " << acos(disp)*180.0/M_PI << std::endl;
    return disp;
}

int main(int argc, char *argv[])
{
    CLI::App app{"Optimizing points on the 3-sphere"};

    
//    int seed = 0;
//    app.add_option("-s,--seed", seed, "random number seed");

    int n = 1;
    app.add_option("-n", n, "number of iterations");
    
    bool rot = false;
    app.add_flag("-r,--rotation,--so3", rot, "treat points as rotations");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");

    CLI11_PARSE(app, argc, argv);

    if (*fileopt)
    {
        std::ifstream fin(filename);
        read_points(fin);
        fin.close();
    }
    else
        read_points(std::cin);
    
    int np = V.rows();
    std::cerr << "Points: ";
    std::cerr << np;
    std::cerr << std::endl;

    if (rot)
    {
        V.conservativeResize(2*np,4);
        for (int i = 0; i < np; i++)
            V.row(i+np) = -V.row(i);
       // np *= 2;
    }
    
    std::vector<std::pair<Point,int>> pi(V.rows());
    for (int i = 0; i < V.rows(); i++)
        pi[i] = std::make_pair(Point(V(i,0),V(i,1),V(i,2),V(i,3)),i);
    Search_traits traits;
    CGAL::spatial_sort(pi.begin(),pi.end(),traits);
    po2vo.resize(V.rows()); vo2po.resize(V.rows());
    auto it = pi.begin();
    for (int i = 0; i < V.rows(); i++, ++it)
    {
        po2vo[i] = it->second;
        vo2po[it->second] = i;
    }
    
    for (int i = 0; i < n; i++)
    {
        create_delaunay();
        double odisp = compute_AV();
        double ovol = vol();
        Eigen::MatrixXd DV = Eigen::MatrixXd::Zero(V.rows(),4);
        for (int f = 0; f < F.rows(); f++)
        {
            for (int j = 0; j < 4; j++)
            {
                Eigen::RowVector4d va = AV.row(f)-V.row(F(f,j));
                DV.row(F(f,j)) += /* va.squaredNorm() * */ va;
            }
        }
//        for (int v = 0; v < V.rows(); v++)
//        {
//            DV.row(v) -= (DV.row(v).dot(V.row(v))) * V.row(v);
//        }
        double m = DV.rowwise().squaredNorm().sum();
        double mdv = sqrt(DV.rowwise().squaredNorm().maxCoeff());
        double ms = /*0.25*/shortest_edge()/mdv;
        // backtrack gradient descent with ms as starting step
        //double ovol = vol();
        Eigen::MatrixXd OV = V;
        double cdisp, ldisp = odisp;
        double cvol, lvol = ovol;
        while (ms > 1e-12)
        {
            std::cerr << "ms: " << ms << std::endl;
            V = OV+ms*DV;
            V.rowwise().normalize();
            cvol = vol();
            cdisp = disp();
//            if (cdisp >= odisp && cdisp <= ldisp) //+ 0.25*ms*m)
//                break;
            if (cvol <= ovol && cvol >= lvol) //+ 0.25*ms*m)
                break;
            ldisp = cdisp;
            lvol = cvol;
            ms *= 0.5;
        }
        ms *= 2.0;
        V = OV+ms*DV;
        V.rowwise().normalize();
        if (rot)
        {
            for (int v = 0; v < np; v++)
            {
                Eigen::RowVector4d p = V.row(v);
                p -= V.row(v+np);
                p *= 0.5;
                V.row(v) = p;
                V.row(v+np) = -p;
            }
        }
        if (ms <= 1e-12)
        {
            std::cerr << "Breaking after " << i << " iterations, no more improvement" << std::endl;
            break;
        }
    }
    
    std::cerr << "Final: " << std::endl;
    create_delaunay();
    compute_AV();
    vol();
    
    if (rot)
    {
        Eigen::VectorXd dd(np);
        for (int i = 0; i < np; i++)
        {
            dd(i) = (V.row(i)+V.row(i+np)).norm();
            //std::cerr << V.row(i) << " / " << V.row(i+np) << std::endl;
        }
        std::cerr << "Asymmetry:  " << dd.maxCoeff() << std::endl;
    }
        
    return 0;
}
