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

    int nvit = t.number_of_vertices();
    
    std::cerr << "Triangulation successfully computed: "
    << nvit << " vertices, "
    << t.number_of_finite_full_cells() << " finite cells.\n Establishing tetrahedralization... "
    << std::endl;

    nvit--;
    if (nvit < np)
    {
        V = Eigen::Matrix<double,Eigen::Dynamic,4>(nvit,4);
    }
    CGAL::Unique_hash_map<T::Vertex_handle, int> VH;
    T::Vertex_iterator vit = t.vertices_begin();
    int i = 0;
    for (; vit != t.vertices_end(); ++vit)
        if (!t.is_infinite(vit) &&
            vit->point() != Point(CGAL::Origin()) )
        {
            if (nvit < np)
            {
                VH[vit] = i;
                if (nvit < np)
                    for (int j = 0; j < 4; j++)
                        V(i,j) = vit->point()[j];
            }
            else
                VH[vit] = po2vo[i];
            i++;
        }
    //assert(np==i);
    np = i;
    
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
    std::cerr << "done.\nThere are " << nf << " facets on the convex hull of the Delaunay triangulation."<< std::endl;
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
    CLI::App app{"Statistics based on the convex hull (Delaunay / Voronoi) of points on the 3-sphere"};

    
    bool rot = false;
    app.add_flag("-r,--rotation,--so3", rot, "treat points as rotations");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");
    
    bool dispersion = false;
    app.add_flag("-d,--dispersion", dispersion, "report dispersion");

    bool volume_variance = false;
    app.add_flag("-v,--volume_variance", volume_variance, "report volume variance");

    int nvb = 0;
    app.add_option("--nvb", nvb, "number of bins for volume histogram");

    bool shortest_edge = false;
    app.add_flag("-e,--shortest_edge", shortest_edge, "report shortest edge");

    bool sphericity = false;
    app.add_flag("-s,--sphericity", sphericity, "report shortest edge");
    
    bool tetcount = false;
    app.add_flag("-t,--tetcount", tetcount, "report number of tets (relative to lower bound)");


    
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
    if (nvb == 0) std::cout << np;
    else std::cerr << np;
    std::cerr << std::endl;

    if (rot)
    {
        V.conservativeResize(2*np,4);
        for (int i = 0; i < np; i++)
            V.row(i+np) = -V.row(i);
        np += np;
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


    create_delaunay();
    analyze_tets();

    // per tet measures: circumradii and vertex degrees
    
    Eigen::VectorXi degree;
    
    if (sphericity) // || degree stuff
    {
        degree = Eigen::VectorXi::Zero(V.rows());
        for (int f = 0; f < F.rows(); f++)
             for (int j = 0; j < 4; j++)
                degree[F(f,j)]++;
    }
    
    if (dispersion )
    {
        Eigen::VectorXd cr(F.rows());
        for (int f = 0; f < F.rows(); f++)
        {
            //circumradius of this tet
            if (degtet[f])
                cr[f] = 0.0;
            else
                cr[f] = acos(fabs(CC.row(f).dot(V.row(F(f,0)))));
        }
        int fid;
        double disp = cr.maxCoeff(&fid);
        double disp_bound = radius(ab(V.rows()));
        std::cerr << "Dispersion:";
         std::cout << " " << disp;
        std::cerr << std::endl << "Bound:";
        std::cout << " " << disp_bound;
        std::cerr << std::endl << "Ratio:";
        std::cout << " " << (disp/disp_bound) << std::endl;
    }
    
    if (tetcount)
    {
        double v = vol(Eigen::VectorXd::Constant(6,ab(V.rows())));
        double ltb = s3a/v;
        std::cerr << "Number of facets relative to lower bound:";
        std::cout << " " << ((double)F.rows() / ltb) << std::endl;
    }
    
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

    std::cerr << "V: " << V.rows() << ", E: " << ne << ", F: " << (2*F.rows()) << ", T: " << F.rows() << std::endl;
    std::cerr << "V-E+F-T = " << (V.rows() - ne + F.rows()) << std::endl;
    
    
    if (shortest_edge)
    {
        // per edge measures: edge lengths
        Eigen::VectorXd el(ne);
        for (auto e: eidx)
            el[e.second] = acos(V.row(e.first[0]).dot(V.row(e.first[1])));
        std::cerr << "Shortest point-point distance:";
        std::cout << " " << el.minCoeff() << std::endl;
    }
    
    if (volume_variance || nvb > 0)
    {
        // compute per edge dual face centroids
        Eigen::Matrix<double,Eigen::Dynamic,4> EC(ne,4);
        for (int f = 0; f < F.rows(); f++)
            for (int j = 0; j < 6; j++)
                EC.row(eidx[edge{F(f,eo[j][0]),F(f,eo[j][1])}]) += CC.row(f);
        EC.rowwise().normalize();
    
        // volume of voronoi cells, distributed from tets
        Eigen::VectorXd v = Eigen::VectorXd::Zero(np);
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
                    if (fabs(P.determinant()) > 1e-8)
                        v[v0] += 0.5*vol(dihedrals(P));
                    P.col(3) = V.row(v1);
                    if (fabs(P.determinant()) > 1e-8)
                        v[v1] += 0.5*vol(dihedrals(P));
                }
            }
        }
     
        std::cerr << "Total volume (sanity check): " << v.sum() << "  ( expected: " << (2.0*M_PI*M_PI) << " ) " <<  std::endl;

        if (nvb > 0)
        {
            double ap = s3a/(double)np; // hyperarea per sample
            double bw = 2.0*ap / (double)nvb;
            std::vector<int> hv(nvb);
            for (int i = 0; i < np; i++)
            {
                int b = (int)(v[i]/bw);
                if (b > nvb-1) continue;
                hv[std::max(b,0)]++;
            }
            for (int i = 0; i < nvb; i++)
            {
                std::cout << ((double)i+0.5)*bw/ap << " " << (hv[i]/2) << std::endl;
            }
        }
        
        if (volume_variance)
        {
            double vm = v.mean();
            Eigen::VectorXd vv = v - vm*Eigen::VectorXd::Ones(np);
           // std::cerr << "zero mean? : " << vv.mean() << std::endl;
            std::cerr << "Volume variance:";
            std::cout << " " << std::sqrt(vv.squaredNorm()/(double)np) << std::endl;
        }
    }

    if (sphericity)
    {
        // roundness of voronoi cells, distributed from tets
        Eigen::VectorXd ah = Eigen::VectorXd::Zero(np);
        for (int f = 0; f < F.rows(); f++)
        {
            for (int j = 0; j < 4; j++)
            {
                int vi = F(f,j);
                ah[vi] += CC.row(f).dot(V.row(vi));
            }
        }
        for (int i = 0; i < np; i++)
            ah[i] /= (double)degree[i];

        Eigen::VectorXd vh = Eigen::VectorXd::Zero(np);
        for (int f = 0; f < F.rows(); f++)
        {
            for (int j = 0; j < 4; j++)
            {
                int vi = F(f,j);
                double d = CC.row(f).dot(V.row(vi)) - ah[vi];
                vh[vi] += d*d;
            }
        }
        for (int i = 0; i < np; i++)
            vh[i] /= (double)degree[i];
        
        double f = boost::math::ibeta_inv(1.5,0.5,s3a/(double)np);
        vh /= f;
        
        std::cerr << "Sphericity (min,max): " << vh.minCoeff() << " ,";
        std::cout << " " << vh.maxCoeff() << std::endl;
    }
    
//
//  degree stuff
//    int maxd = degree.maxCoeff();
//    Eigen::VectorXi hd = Eigen::VectorXi::Zero(maxd+1);
//    for (int i = 0; i < np; i++)
//        hd[degree[i]]++;
//    std::cerr << "Vertex degrees: " << std::endl;
//    for (int i = 0; i <= maxd; i++)
//        if (hd[i] > 0)
//            std::cerr << i << ": " << hd[i] << ", ";//std::endl;

    return 0;
}
