#include "s3io.h"

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

#include <boost/math/special_functions/beta.hpp>

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

    std::cerr << "Triangulation successfully computed: "
    << t.number_of_vertices() << " vertices, "
    << t.number_of_finite_full_cells() << " finite cells.\n Establishing tetrahedralization... "
    << std::endl;

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
    
    F.resize(0,4); F.resize(nf,4);
    N.resize(0,4); N.resize(nf,4);
    
    //F = Eigen::Matrix<int,Eigen::Dynamic,4>(nf,4);
    //N = Eigen::Matrix<int,Eigen::Dynamic,4>(nf,4);
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

void analyze_tets()
{
    int nd = 0;
    const Eigen::Vector4d b(1.0,0.0,0.0,0.0);
    CC.resize(0,4); CC.resize(F.rows(),4);
    //CC = Eigen::Matrix<double,Eigen::Dynamic,4>(F.rows(),4);
    degtet.resize(0);
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


int main(int argc, char *argv[])
{
    CLI::App app{"Statistics based on the convex hull (Delaunay / Voronoi) of points on the 3-sphere"};

    
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
        V = read_points(fin);
        fin.close();
    }
    else
        V = read_points(std::cin);
    
    int np = V.rows();
    std::cerr << "Points: ";
    std::cout << np;
    std::cerr << std::endl;

    if (rot)
    {
        V.conservativeResize(2*np,4);
        for (int i = 0; i < np; i++)
            V.row(i+np) = -V.row(i);
        //np += np;
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
        analyze_tets();

        Eigen::VectorXd cr(F.rows());
        for (int f = 0; f < F.rows(); f++)
        {
            //circumradius of this tet
            if (degtet[f])
                cr[f] = 0.0;
            else
                cr[f] = acos(fabs(CC.row(f).dot(V.row(F(f,0)))));
        }
        std::cerr << "Dispersion:";
        int fid;
        std::cout << " " << (cr.maxCoeff(&fid)*180.0/M_PI) << std::endl;
 
        // per edge indices
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

        // compute per edge dual face centroids
        Eigen::Matrix<double,Eigen::Dynamic,4> EC = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(ne,4);
        for (int f = 0; f < F.rows(); f++)
            for (int j = 0; j < 6; j++)
                EC.row(eidx[edge{F(f,eo[j][0]),F(f,eo[j][1])}]) += CC.row(f);
        EC.rowwise().normalize();
        
        
    
        // volume of voronoi cells, distributed from tets
        Eigen::VectorXd v = Eigen::VectorXd::Zero(V.rows());
        Eigen::MatrixXd CM = Eigen::MatrixXd::Zero(V.rows(),4);
        Eigen::Matrix4d P;
        for (int f = 0; f < F.rows(); f++)
        {
            P.col(0) = CC.row(f);
            for (int j = 0; j < 6; j++)
            {
                int v0 = F(f,eo[j][0]);
                int v1 = F(f,eo[j][1]);
                P.col(1) = EC.row(eidx[edge{v0,v1}]);
                int fn = N(f,eo[j][2]);
                for (int k = 2; k < 4; k++)
                {
                    P.col(2) = CC.row(N(f,eo[j][k]));
                    P.col(3) = V.row(v0);
                    if (fabs(P.determinant()) > 1e-12)
                    {
                        double vv = 0.5*vol(dihedrals(P));
                        v[v0] += vv;
                        CM.row(v0) += vv * (P.rowwise().sum()).normalized();
                    }
                    P.col(3) = V.row(v1);
                    if (fabs(P.determinant()) > 1e-12)
                    {
                        double vv = 0.5*vol(dihedrals(P));
                        v[v1] += vv;
                        CM.row(v1) += vv * (P.rowwise().sum()).normalized();
                    }
                }
            }
        }
     //   std::cerr << v.transpose() << std::endl;
     
        std::cerr << "Total volume (sanity check): " << v.sum() << "  ( expected: " << (2.0*M_PI*M_PI) << " ) " <<  std::endl;

        double vm = v.mean();
        Eigen::VectorXd vv = v - vm*Eigen::VectorXd::Ones(V.rows());
        // std::cerr << "zero mean? : " << vv.mean() << std::endl;
        std::cerr << "Volume variance:";
        std::cout << " " << std::sqrt(vv.squaredNorm()/(double)V.rows()) << std::endl;

        CM.rowwise().normalize();
       // std::cerr << (V-CM) << std::endl;
        
        
        V = CM;
        
//        if (rot)
//        {
//            for (int v = 0; v < np; v++)
//            {
//                Eigen::RowVector4d p = V.row(v);
//                p -= V.row(v+np);
//                p *= 0.5;
//                V.row(v) = p;
//                V.row(v+np) = -p;
//            }
//        }
//        V.rowwise().normalize();
    }
    return 0;
}
