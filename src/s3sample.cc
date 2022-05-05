
#include <Eigen/Dense>

#include "CLI11.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iterator>
#include <vector>

const double offset = 0.5;
const double sn = 1.533751168755204288118041;

double rootp(std::vector<int> coef)
{
    std::cerr << "poly: " << coef[0];
    for (int i = 1; i < coef.size(); i++)
    std::cerr << " + " << coef[i] << " x^" << i;
    std::cerr << std::endl;;
    
    double x = 2.0;
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


void write_points(std::basic_ostream<char>& vf, Eigen::Matrix<double,Eigen::Dynamic,4> V)
{
    int n = V.rows();
    vf.precision(16);
    vf << n << "\n";
    for (int i = 0; i < n; i++)
    {
        vf << V(i,0) << " " << V(i,1) << " " << V(i,2) << " " << V(i,3) << "\n";
    }
}

Eigen::Matrix<double,Eigen::Dynamic,4> double_fib(int n, const double mc0, const double mc1, bool time2cout = false)
{
    std::cerr << "Creating ";
    if (time2cout)
        std::cout << n << " ";
    else
        std::cerr << n << " ";
    std::cerr << "points by the fibonacci method." << std::endl;
    std::cerr << std::setprecision(16) << "Using phi = " << mc0 << " , psi = " << mc1 << std::endl;

    auto tp0 = std::chrono::steady_clock::now();
    Eigen::Matrix<double,Eigen::Dynamic,4> V(n,4);
    double s,ab,r,R,theta,phi;
    double dn = 1.0/(double)n;
    for (int i = 0; i < n; i++)
    {
        s = (double)i+offset;
        ab = 2.0 * M_PI * s;
        theta = ab * mc0;
        phi = ab * mc1;
        s *= dn;
        r = sqrt(s);
        R = sqrt(1.0-s);
        V.row(i) << r*sin(theta),r*cos(theta),R*sin(phi),R*cos(phi);
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in ";
    if (time2cout)
    {
        std::cout << (std::chrono::duration <double, std::milli> (tp1-tp0).count());
        std::cerr << " ms";
        std::cout << std::endl;
    }
    else
        std::cerr << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    return V;
}

Eigen::Matrix<double,Eigen::Dynamic,4> uniform_random(int n, bool time2cout = false, int seed = 0)
{
    std::cerr << "Creating ";
    if (time2cout)
        std::cout << n << " ";
    else
        std::cerr << n << " ";
    std::cerr << "points by uniform random sampling." << std::endl;
    std::default_random_engine rng(seed);
    std::normal_distribution<double> dn;
    //std::uniform_real_distribution<double> du;

    auto tp0 = std::chrono::steady_clock::now();
    Eigen::Matrix<double,Eigen::Dynamic,4> V(n,4);
    for (int i = 0; i < n; i++)
    {
       V.row(i) << dn(rng),dn(rng),dn(rng),dn(rng);
       V.row(i).normalize();
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in ";
    if (time2cout)
    {
        std::cout << (std::chrono::duration <double, std::milli> (tp1-tp0).count());
        std::cerr << " ms";
        std::cout << std::endl;
    }
    else
        std::cerr << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;

    return V;
}


int main(int argc, char *argv[])
{
    CLI::App app{"Super-Fibonacci-Sampling the 3-sphere"};

    
    int seed = 0;
    app.add_option("-s,--seed", seed, "random number seed");

    int n = 1;
    app.add_option("-n", n, "number of points");
    
    bool uniform = false;
    app.add_flag("-u,--uniform0", uniform, "Use uniform sampling instead of Super-Fibonacci");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");

    bool sout = false;
    app.add_flag("-o",sout,"Output points to cout");

    int centering = 0;
    app.add_option ("-c", centering, "Centering the point set (multiple times)");
    
    bool reflect = false;
    app.add_flag("-r,--reflect", reflect, "Generate half of n points, then reflect them");

    std::vector<int> phi_coef, psi_coef;
    app.add_option("--phicoef", phi_coef, "Coefficients for polynomial, whose root around 1 is phi");
    app.add_option("--psicoef", psi_coef, "Coefficients for polynomial, whose root around 1 is psi");

    int phi_exp = -1, psi_exp = -1;
    app.add_option("--phiexp", phi_exp, "Exponent for phi (-1)");
    app.add_option("--psiexp", psi_exp, "Exponent for psi (-1)");

    
    CLI11_PARSE(app, argc, argv);

    Eigen::Matrix<double,Eigen::Dynamic,4> V;
    std::cerr << "Sampling " << n << " points on S3" << std::endl;
    if (reflect) n /= 2;

    if (uniform)
    {
        std::cerr << "Using uniform sampling!" << std::endl;
        V = uniform_random(n,!sout,seed); 
    }
    else
    {
        double iphi = 1.0/sqrt(2.0);
        if (phi_exp != -1) iphi = std::pow(sqrt(2.0),phi_exp);
        if (phi_coef.size() > 0)
            iphi = std::pow(rootp(phi_coef),phi_exp);
        double ipsi = 1.0/sn;
        if (psi_exp != -1) ipsi = std::pow(sn,psi_exp);
        if (psi_coef.size() > 0)
            ipsi = std::pow(rootp(psi_coef),psi_exp);
        V = double_fib(n,iphi,ipsi,!sout);
    }

    for (int i = 0; i < centering; i++)
    {
        Eigen::RowVector4d m = V.colwise().mean();
        V.rowwise() -= m;
        V.rowwise().normalize();
    }
    
    if (reflect)
    {
        V.conservativeResize(2*n,Eigen::NoChange);
        for (int i = 0; i < n; i++)
        {
            V(i+n,0) = V(i,2);
            V(i+n,1) = V(i,3);
            V(i+n,2) = V(i,0);
            V(i+n,3) = V(i,1);
        }
    }
    
    if (*fileopt)
    {
        std::ofstream fout(filename);
        write_points(fout,V);
        fout.close();
    }
    if (sout)
        write_points(std::cout,V);

    return 0;
}
