#include "s3io.h"

#include <Eigen/Dense>
#include "CLI11.hpp"

#include <math.h>

const double s3a = 2.0*M_PI*M_PI;

int main(int argc, char *argv[])
{
    CLI::App app{"Discrepancy on the 3-sphere"};

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
    
    int np = V.rows();
    std::cerr << "Points: ";
    std::cout << np;
    std::cerr << std::endl;

    // assume points are rotations and always add opposite point
    V.conservativeResize(2*np,4);
    for (int i = 0; i < np; i++)
        V.row(i+np) = -V.row(i);
    
    double ap = s3a/(double)(2*np);
    
    long n = np;
    n = n*(n-1)/2 + n*(n-1)*(n-2)/6 + n*(n-1)*(n-2)*(n-3)/24;
    std::cerr << "Number of circles: " << n << std::endl;
    
    double md2 = 0, md3 = 0, md4 = 0;

    n = 0;
    for (int i0 = 0; i0 < np-1; i0++)
        for (int i1 = i0+1; i1 < np; i1++)
        {
            Eigen::Vector4d v0 = V.row(i0);
            Eigen::Vector4d v1 = V.row(i1);
            Eigen::Vector4d c;
            double h = v0.dot(v1);
            if (h > 0)
                c = (v0+v1).normalized();
            else
            {
                c = (v0-v1).normalized();
                h *= -1;
            }
            int ni = 0;
            for (int j = 0; j < 2*np; j++)
                if (j != i0 && j != i1 && j != (i0+np) && j != (i1+np) && V.row(j).dot(c) > h)
                    ni++;
            double a = 2.0*M_PI * (std::acos(h) - std::sqrt(1-h*h)*h);
            double cd = std::max(fabs(a-(double)ni * ap),
                                 fabs(a-(double)(ni+2) * ap));
            if (cd > md2) md2 = cd;
            n++;
        }
    std::cerr << "Checked " << n << " pairs. Max discrepancy: " << md2 << std::endl;
    n = 0;
    for (int i0 = 0; i0 < np-2; i0++)
        for (int i1 = i0+1; i1 < np-1; i1++)
            for (int i2 = i1+1; i2 < np; i2++)
            {
                Eigen::Vector4d v0 = V.row(i0);
                Eigen::Vector4d v1 = V.row(i1);
                Eigen::Vector4d v2 = V.row(i2);
                if (v0.dot(v1) < 0) v1 *= -1;
                if (v0.dot(v2) < 0) v2 *= -1;
                Eigen::Matrix3d M(3,3);
                M.row(0) << (v0-v1).dot(v0), (v0-v1).dot(v1), (v0-v1).dot(v2);
                M.row(1) << (v0-v2).dot(v0), (v0-v2).dot(v1), (v0-v2).dot(v2);
                M.row(2) << 1,1,1;
                Eigen::Vector3d y = M.partialPivLu().solve(Eigen::Vector3d(0,0,1));
                Eigen::Vector4d c = (y[0]*v0+y[1]*v1+y[2]*v2).transpose();
                c.normalize();
                double h = c.dot(v0);
                int ni = 0;
                for (int j = 0; j < 2*np; j++)
                    if (j != i0 && j != i1 && j != i2 && j != (i0+np) && j != (i1+np) && j != (i2+np) && V.row(j).dot(c) > h)
                        ni++;
                double a = 2.0*M_PI * (std::acos(h) - std::sqrt(1-h*h)*h);
                double cd = std::max(fabs(a-(double)ni * ap),
                                     fabs(a-(double)(ni+3) * ap));
                if (cd > md3) md3 = cd;
                n++;
            }
    std::cerr << "Checked " << n << " triples. Max discrepancy: " << md3 << std::endl;

    n = 0;
    for (int i0 = 0; i0 < np-3; i0++)
        for (int i1 = i0+1; i1 < np-2; i1++)
            for (int i2 = i1+1; i2 < np-1; i2++)
                for (int i3 = i2+1; i3 < np; i3++)
                {
                    Eigen::Vector4d v0 = V.row(i0);
                    Eigen::Vector4d v1 = V.row(i1);
                    Eigen::Vector4d v2 = V.row(i2);
                    Eigen::Vector4d v3 = V.row(i3);
                    if (v0.dot(v1) < 0) v1 *= -1;
                    if (v0.dot(v2) < 0) v2 *= -1;
                    if (v0.dot(v3) < 0) v3 *= -1;
                    Eigen::Matrix4d M(4,4);
                    M.row(0) << (v0-v1).dot(v0), (v0-v1).dot(v1), (v0-v1).dot(v2), (v0-v1).dot(v3);
                    M.row(1) << (v0-v2).dot(v0), (v0-v2).dot(v1), (v0-v2).dot(v2), (v0-v2).dot(v3);
                    M.row(2) << (v0-v3).dot(v0), (v0-v3).dot(v1), (v0-v3).dot(v2), (v0-v3).dot(v3);
                    M.row(3) << 1,1,1,1;
                    Eigen::Vector4d y = M.partialPivLu().solve(Eigen::Vector4d(0,0,0,1));
                    Eigen::Vector4d c = y[0]*v0+y[1]*v1+y[2]*v2+y[3]*v3;
                    c.normalize();
                    double h = c.dot(v0);
                    int ni = 0;
                    for (int j = 0; j < 2*np; j++)
                        if (j != i0 && j != i1 && j != i2 && j != i3 && j != (i0+np) && j != (i1+np) && j != (i2+np) && j != (i3+np) && V.row(j).dot(c) > h)
                            ni++;
                    double a = 2.0*M_PI * (std::acos(h) - std::sqrt(1-h*h)*h);
                    double cd = std::max(fabs(a-(double)ni * ap),
                                         fabs(a-(double)(ni+4) * ap));
                    if (cd > md4) md4 = cd;
                   n++;
                }
    
    std::cerr << "Checked " << n << " 4-tuples. Max discrepancy: " << md4 << std::endl;

    std::cerr << "Overall: ";
    std::cout << std::max(std::max(md2,md3),md4) << "\n";
        
    return 0;
}
