#include "s3io.h"

#include <boost/math/special_functions/beta.hpp>
#include <boost/sort/spreadsort/spreadsort.hpp>

#include <Eigen/Dense>

#include "CLI11.hpp"

#include <random>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
//#include <execution>

#include <math.h>       /* acos */


const double s3a = 2.0*M_PI*M_PI;
const double k = 40.0;
Eigen::VectorXd ibs;
int ns;


double zpf(Eigen::RowVector4d q)
{
    double s = 0.0;
    if (fabs(q[1]) > 1e-16) s += q[1]/sqrt(2.0*(q[0]*q[0]+q[1]*q[1]));
    if (fabs(q[3]) > 1e-16) s += q[3]/sqrt(2.0*(q[2]*q[2]+q[3]*q[3]));
    double r = std::acos(s/sqrt(2.0));
//    double r = std::acos((q[1]+q[3])/sqrt(2.0));
    return 0.5*(1.0+cos(k*r*r));
}

int main(int argc, char *argv[])
{
    CLI::App app{"Clifford torus visualization for points on the 3-sphere"};

    int nx = 256;
    int ny = 256;
    
    int seed = 0;
    app.add_option("-s,--seed", seed, "random number seed (for rotation)");

    double t = 1.0;
    app.add_option("-t,-w,-k", t, "thickness, kernel width, in radians");

    bool rot = false;
    app.add_flag("-r,--rotation,--so3", rot, "treat points as rotations");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");

    bool zps = false;
    app.add_flag("-z,--sample", zps, "use points to sample zone plate");

    
    CLI11_PARSE(app, argc, argv);

    Eigen::Matrix<double,Eigen::Dynamic,4> V;
    if (*fileopt)
    {
        std::ifstream fin(filename);
        V = read_points(fin);
        fin.close();
    }
    else
        V = read_points(std::cin);
    
    double d = 0.0;
    int np = V.rows();
    std::cerr << "Points: ";
//    std::cout << np;
    std::cerr << np;
    std::cerr << std::endl;

    if (rot)
    {
        V.conservativeResize(2*np,4);
        for (int i = 0; i < np; i++)
            V.row(i+np) = -V.row(i);
        np += np;
    }
    
    if (seed > 0)
    {
        std::default_random_engine rng(seed);
        std::normal_distribution<double> dn;
        Eigen::Matrix4d X;
        for (int r = 0; r < 4; r++)
            X.row(r) << dn(rng),dn(rng),dn(rng),dn(rng);
        std::cerr << X << std::endl;
        Eigen::Matrix4d Q = X.fullPivHouseholderQr().matrixQ();
        Eigen::Matrix4d L = Eigen::Matrix4d::Zero();
        for (int r = 0; r < 4; r++)
            L(r,r) = Q(r,r)/std::abs(Q(r,r));
    
        // std::cerr << (Q*L) << std::endl;
    
        V *= (Q*L);
    }
    
    double ap = s3a/(double)np;
    
    double ct = std::sqrt(2.0)*std::cos(t);

    std::cerr << "ct: " << ct << std::endl;
    int nd = 0;

    if (zps)
    {
        Eigen::MatrixXd F = Eigen::MatrixXd::Zero(nx,ny);
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nx,ny);
        int r = 10;

        double sigma = 2.0*M_PI/sqrt(nx*nx+ny*ny);
        double wm = 1.0;
        for (int i = 0; i < np; i++)
        {
            
            Eigen::RowVector4d q = V.row(i);
            double s0 = std::sqrt(q[0]*q[0]+q[1]*q[1]);
            double s1 = std::sqrt(q[2]*q[2]+q[3]*q[3]);
            if (s0+s1<ct) continue;

            nd++;
            //std::cerr << i << std::endl;
            //std::cerr << q << std::endl;

            double x = std::atan2(q[0],q[1]);
            double y = std::atan2(q[2],q[3]);
            int px = (x+M_PI)*(double)nx/(2.0*M_PI);
            int py = (y+M_PI)*(double)ny/(2.0*M_PI);
            double z = zpf(q);
            
            for (int cx = std::max(0,px-r); cx <= std::min(px+r,nx-1); cx++)
                for (int cy = std::max(0,py-r); cy <= std::min(py+r,ny-1); cy++)
                {
                    x = (double)cx*(2.0*M_PI)/(double)nx - M_PI;
                    y = (double)cy*(2.0*M_PI)/(double)ny - M_PI;
                    Eigen::RowVector4d c;
                    c << sin(x), cos(x), sin(y), cos(y);
                    c *= 1/std::sqrt(2.0);
                    double d = std::acos(fabs(c.dot(q)));
                    double w = exp(-d*d / (sigma*sigma));
                    if (cx == px && cy == py)
                        if (w < wm) wm = w;
                    K(cx,cy) += w;
                    F(cx,cy) += w*z;
                    //if (w < wm) wm = w;
                }
         }
       std::cerr << wm << std::endl;
       std::cout << "P5\n" << nx << " " << ny << "\n255\n";
        for (int cy = 0; cy < ny; cy++)
            for (int cx = 0; cx < nx; cx++)
            {
                int g = 0;
                if (K(cx,cy) > 0.0)
                    g = (int)(255.0*F(cx,cy)/K(cx,cy));
                g = std::max(0,std::min(255,g));
                std::cout << (char)g;
            }
    }
    else
    {
        std::cout << R"(<svg viewBox="0 0 100 100" xmlns="http://www.w3.org/2000/svg">)" << std::endl;
    
        std::cout.setf( std::ios::fixed );
        std::cout.precision( 2 );
        for (int i = 0; i < np; i++)
        {
            Eigen::RowVector4d q = V.row(i);
            double s0 = std::sqrt(q[0]*q[0]+q[1]*q[1]);
            double s1 = std::sqrt(q[2]*q[2]+q[3]*q[3]);
            if (s0+s1<ct) continue;

            double x = std::atan2(q[0],q[1]);
            double y = std::atan2(q[2],q[3]);
            x += M_PI; y += M_PI;
            x *= 100.0/(2.0*M_PI);
            y *= 100.0/(2.0*M_PI);

            std::cout << R"(<circle cx=")" << x << R"(" cy=")" << y << R"(" r = ")"<< ((s0+s1-ct)/(std::sqrt(2.0)-ct)) << R"("/>)" << std::endl;
            nd++;
        }
        std::cout << R"(</svg>)" << std::endl;
    }
    std::cerr << "Points used: " << nd << " / " << np << std::endl;

    return 0;
}
