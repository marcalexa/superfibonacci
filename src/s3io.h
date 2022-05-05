
#ifndef INCLUDED_S3IO_HEADER
#define INCLUDED_S3IO_HEADER

#include <Eigen/Dense>

#include <iostream>


Eigen::Matrix<double,Eigen::Dynamic,4> read_points(std::basic_istream<char>& vf)
{
    std::string line;
    // first number found on a line that is not a comment is assumed to be number of points
    int n = -1;
    while (n < 0 && std::getline(vf >> std::ws, line))
    {
        if (std::isdigit(line[0]))
        {
            std::istringstream iss(line);
            iss >> n;
        }
    }
    std::cerr << "Trying to read " << n << " points... " << std::flush;
    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    int i = 0;
    while (i < n && std::getline(vf >> std::ws, line))
    {
        if (std::isdigit(line[0]) || line[0] == '-' || line[0] == '.' )
        {
            std::istringstream iss(line);
            iss >> V(i,0) >> V(i,1) >> V(i,2) >> V(i,3);
            if ( fabs(1.0-V.row(i).squaredNorm()) > 1e-12 )
            {
                V.row(i).normalize();
            }
            i++;
        }
    }
    std::cerr << "done!" << std::endl;
    if (i < n)
    {
        n = i;
        std::cerr << "Warning: Could only find " << n << "points." << std::endl;
        V.conservativeResize(n,Eigen::NoChange);
    }
    if (std::getline(vf >> std::ws, line))
    {
        std::cerr << "Warning: there are additional lines. The first ones is:\n " << line << std::endl;
    }
    
    return V;
}

#endif
