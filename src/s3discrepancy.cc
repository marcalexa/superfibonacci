#include "s3io.h"

//#include <boost/math/special_functions/beta.hpp>
#include <boost/sort/spreadsort/spreadsort.hpp>

#include <Eigen/Dense>

#include "CLI11.hpp"

#include <random>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>    // std::for_each
#include <execution>

#include <math.h>       /* acos */


const double s3a = 2.0*M_PI*M_PI;
Eigen::VectorXd ibs;
int ns;


Eigen::Matrix<double,Eigen::Dynamic,4> uniform_random(int n, int seed = 0)
{
    std::default_random_engine rng(seed);
    std::normal_distribution<double> dn;
    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < 4; j++) V(i,j) = dn(rng);
        V.row(i).normalize();
    }
    
    return V;
}



int main(int argc, char *argv[])
{
    CLI::App app{"Sampling discrepancy on the 3-sphere"};

    
    int seed = 0;
    app.add_option("-s,--seed", seed, "random number seed");

    int n = 1;
    app.add_option("-n", n, "number of sample directions");

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
    std::cout << np;
    std::cerr << std::endl;

    if (rot)
    {
        V.conservativeResize(2*np,4);
        for (int i = 0; i < np; i++)
            V.row(i+np) = -V.row(i);
        np += np;
    }
    Eigen::Matrix<double,Eigen::Dynamic,4> N = uniform_random(n,seed);
    
    double ap = s3a/(double)np;


    Eigen::VectorXd dd = Eigen::VectorXd::Zero(n);
    
    std::vector<int> a(n);
    std::iota(std::begin(a), std::end(a), 0);
    std::for_each(
                  std::execution::par,
                  std::begin(a), std::end(a), [&](int i)
                  {
                      Eigen::VectorXd h = V*N.row(i).transpose();
                      std::vector<double> lh(h.data(),h.data()+np);
                      boost::sort::spreadsort::spreadsort(lh.begin(), lh.begin()+np);

                      for (int j = 0; j < np; j++)
                      {
                          if (lh[j] > 0.0) break;
                          double h = -lh[j];
                          double a = 2.0*M_PI * (std::acos(h) - std::sqrt(1-h*h)*h);
                          double cd = std::max(fabs(a-(double)j * ap),
                                               fabs(a-(double)(j+1) * ap));
                          if (cd > dd[i])
                          {
                              dd[i] = cd;
                          }
                      }
                  });
    
    d = dd.maxCoeff();

    std::cerr << "Discrepancy (" << n << " random centers, parallel, all critical circles): ";
    std::cout << " " << d << "\n";
    
    return 0;
}
