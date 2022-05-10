
#include <iostream>
#include <vector>
#include <array>
#include <cmath>

const double sn = 1.533751168755204288118041;

int n;

std::vector<std::array<double,4> > Q;


void write_points()
{
    std::cout.precision(16);
    for (int i = 0; i < n; i++)
    {
        std::cout << Q[i][0] << " " << Q[i][1] << " " << Q[i][2] << " " << Q[i][3] << "\n";
    }
}

void double_fib(double mc0, double mc1)
{
    std::cerr << "Creating " << n << " points by the fibonacci method." << std::endl;

    double s,ab,r,R,theta,phi;
    double dn = 1.0/(double)n;
    for (int i = 0; i < n; i++)
    {
        s = (double)i+0.5;
        ab = 2.0 * M_PI * s;
        theta = ab * mc0;
        phi = ab * mc1;
        s *= dn;
        r = sqrt(s);
        R = sqrt(1.0-s);
        Q[i][0] = r*sin(theta);
        Q[i][1] = r*cos(theta);
        Q[i][2] = R*sin(phi);
        Q[i][3] = R*cos(phi);
    }
    std::cerr << "Done.\n";
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Expects a single integer argument, specifying the number of points!\n";
        exit(0);
    }
    n = atoi(argv[1]);
    if (n <= 0)
    {
        std::cerr << "Couldn't parse integer argument!\n";
        exit(0);
    }

    Q.resize(n);

    double iphi = 1.0/sqrt(2.0);
    double ipsi = 1.0/sn;
    double_fib(iphi,ipsi);

    write_points();

    return 0;
}
