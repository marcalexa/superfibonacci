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
#include <execution>

#include <cmath>

const double s3a = 2.0*M_PI*M_PI;
const double fm = 100.0;
Eigen::VectorXd ibs;
int ns;



int main(int argc, char *argv[])
{
    CLI::App app{"Pair correlations on the 3-sphere"};

    
    int seed = 0;
    app.add_option("-s,--seed", seed, "random number seed");

    int n = 1;
    app.add_option("-n", n, "number of samples");

    int nh = 100;
    app.add_option("-b", nh, "number of histogram bins");

    double variance = 0.0;
    app.add_option("-v", variance, "variance of kernel");

    bool auto_variance = false;
    app.add_flag("-a,--auto", auto_variance, "estimate variance based sample count");
    
    bool rot = false;
    app.add_flag("-r,--rotation,--so3", rot, "treat points as rotations");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");

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
    
    double ap = s3a/(double)np; // hyperarea per sample
    double bw = (M_PI/2.0)/(double)nh; // width per bin in [0.pi/2]

    
    if (auto_variance)
    {
        // ap contains volume per sample, if sphere, radius is
        double r = std::pow(ap*(3.0/4.0)/M_PI,1.0/3.0);
        double fwxm = 2.0*r;
        double c = fwxm / (2.0*std::sqrt(2.0*std::log(fm)));
        variance = c*c;
        std::cerr << "Auto-variance: " << variance << std::endl;
    }
    int offset = 0;
    if (variance > 0.0)
    {
        offset = 1;
        while (exp( -0.5 * (double)(offset*offset) * bw*bw / variance ) > 1e-8 )
            offset++;
    }
    std::cerr << "offset: " << offset << std::endl;
    
    std::vector<double> h(nh,0);
    
    std::vector<int> a(np);
    std::iota(std::begin(a), std::end(a), 0);
    std::shuffle(a.begin(),a.end(), std::default_random_engine(seed));
    a.resize(n);
    
    std::for_each(
                  std::execution::par,
                  std::begin(a), std::end(a), [&](int i)
                  {
                      Eigen::VectorXd c = V*V.row(i).transpose();
                      
                      for (int j = 0; j < np; j++)
                      {
                          if (j==i || c[j] <= 1e-16 ) continue;
                          double d = std::acos(c[j]);
                          double dw = std::sin(d+0.5*bw); // added half a bin width to avoid division by zero below
                          dw *= dw; dw = 1.0/dw;
                          int b = (int)(d/bw);
                          if (offset > 0)
                          {
                              int k = std::max(0,b-offset);
                              int u = std::min(b+offset,nh-1);
                              for (; k <= u; k++)
                              {
                                  double dd = (double)k*bw - d;
                                  h[k] += dw*std::exp( -0.5*dd*dd/variance );
                              }
                          }
                          else
                          {
                              h[std::min(std::max(b,0),nh-1)] += 1.0;
                          }
                      }
                  });

    std::cerr << "Point correlation, " << nh << " bins"<< std::endl;
    for (int i = 0; i < nh; i++)
    {
        double leb = (double)i*bw;
        double ueb = (double)(i+1)*bw;
        double lh = std::sin(leb);
        double uh = std::sin(ueb);
        double bm = boost::math::ibeta(1.5,0.5,std::max(0.0,uh*uh)) - boost::math::ibeta(1.5,0.5,std::max(0.0,lh*lh));
        std::cout << ((double)i/((double)nh)*(M_PI/2.0)) << " " <<
        ( h[i]*ap*bw ) << std::endl;
    }
    


    return 0;
}
