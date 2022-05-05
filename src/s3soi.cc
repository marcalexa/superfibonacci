
#include <Eigen/Dense>

#include "CLI11.hpp"

#include <chrono>
#include <random>
#include <iostream>
#include <iterator>
#include <vector>

#include "SOI/rotutils.h"


const int Nstep=1000;
const double minimal_step=1.e-10;
Eigen::Matrix<double,Eigen::Dynamic,3> P;
double tt;

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

Eigen::Matrix<double,Eigen::Dynamic,4> double_fib(int n, double mc0, double mc1)
{
    std::cerr << "Creating " << n << " points by the fibonacci method." << std::endl;
    std::cerr << std::setprecision(16) << "Using phi = " << mc0 << " , psi = " << mc1 << std::endl;

    auto tp0 = std::chrono::steady_clock::now();
    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    for (int i = 0; i < n; i++)
    {
        double s = ((double)i+0.5)/((double)n);
        double ab = 2.0 * M_PI * s * (double)n;
        double r = sqrt(s);
        double theta = ab * mc0;
        double R = sqrt(1.0-s);
        double phi = ab * mc1;
        V.row(i) << r*sin(theta),r*cos(theta),R*sin(phi),R*cos(phi);
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    return V;
}

Eigen::Matrix<double,Eigen::Dynamic,4> uniform_random(int n, int seed = 0)
{
    std::cerr << "Creating " << n << " points by uniform random sampling." << std::endl;
    auto tp0 = std::chrono::steady_clock::now();
    std::default_random_engine rng(seed);
    //std::normal_distribution<double> dn;
    std::uniform_real_distribution<double> du;

    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(n,4);
    for (int i = 0; i < n; i++)
    {
//        double s0 = sqrt(-2.0*log(du(rng)));
        double s0 = sqrt(2.0*du(rng));
        double a0 = 2.0*M_PI*du(rng);
//        double s1 = sqrt(-2.0*log(du(rng)));
        double s1 = sqrt(2.0*du(rng));
        double a1 = 2.0*M_PI*du(rng);
       //V.row(i) << dn(rng),dn(rng),dn(rng),dn(rng);
       // V.row(i) << 1.0,0.0,0.0,0.0;
        V.row(i) << s0*cos(a0),s0*sin(a0),s1*cos(a1),s1*sin(a1);
        V.row(i).normalize();
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    return V;
}


int num_icosamp(int L)
{
    return (12 + L*30 + 20*(L-1)*L/2);
}

int num_unifsamp(int n)
{
    return (floor(((double) n*n) / M_PI));
}

double dot(vec v1,vec v2)
{
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

double length(vec v)
{
    return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

double vecdist(vec v1,vec v2)
{
    vec v;
    v[0] = v2[0] - v1[0]; v[1] = v2[1] - v1[1]; v[2] = v2[2] - v1[2];
    return length(v);
}

double get_coulomb_energy(int N,vec p[])
{
    double e = 0;
    int i,j;
    for(i = 0;i<N;i++)
        for(j = i+1; j<N; j++ )
        {
            e += 1/ vecdist(p[i],p[j]);
        }
    return e;
}

void get_forces(int N,vec f[], vec p[])
{
    int i,j;
    vec r;
    double rr, ff, l;
    for(i = 0;i<N;i++)
    {
        f[i][0] = 0;f[i][1] = 0;f[i][2] = 0;
    }
    for(i = 0;i<N;i++)
    {
        for(j = i+1; j<N; j++ )
        {
            r[0] = p[i][0]-p[j][0];
            r[1] = p[i][1]-p[j][1];
            r[2] = p[i][2]-p[j][2];
            l = length(r); l = 1/(l*l*l);
            ff = l*r[0]; f[i][0] += ff; f[j][0] -= ff;
            ff = l*r[1]; f[i][1] += ff; f[j][1] -= ff;
            ff = l*r[2]; f[i][2] += ff; f[j][2] -= ff;
        }
    }
    return;
}


void set_icosa( vec *iverts, edge *iedges, sim *ifaces )
{
    iverts[0][0] = -0.72361; iverts[0][1] = 0.52573;   iverts[0][2] = 0.44721;
    iverts[1][0] = -0.72361; iverts[1][1] = -0.52573;  iverts[1][2] = 0.44721;
    iverts[2][0] = -0.27639; iverts[2][1] = 0.85065;   iverts[2][2] = -0.44721;
    iverts[3][0] = -0.27639; iverts[3][1] = -0.85065;  iverts[3][2] = -0.44721;
    iverts[4][0] = 0.00000;  iverts[4][1] = 0.00000;   iverts[4][2] = 1.00000;
    iverts[5][0] = 0.00000;  iverts[5][1] = 0.00000;   iverts[5][2] = -1.00000;
    iverts[6][0] = 0.27639;  iverts[6][1] = 0.85065;   iverts[6][2] = 0.44721;
    iverts[7][0] = 0.27639;  iverts[7][1] = -0.85065;  iverts[7][2] = 0.44721;
    iverts[8][0] = 0.72361;  iverts[8][1] = 0.52573;   iverts[8][2] = -0.44721;
    iverts[9][0] = 0.72361;  iverts[9][1] = -0.52573;  iverts[9][2] = -0.44721;
    iverts[10][0] = 0.89443; iverts[10][1] = 0.00000;  iverts[10][2] = 0.44721;
    iverts[11][0] = -0.89443;iverts[11][1] = 0.00000;  iverts[11][2] = -0.44721;
    
    iedges[0][0] = 0;  iedges[0][1] = 1;
    iedges[1][0] = 0;  iedges[1][1] = 2;
    iedges[2][0] = 0;  iedges[2][1] = 4;
    iedges[3][0] = 0;  iedges[3][1] = 6;
    iedges[4][0] = 0;  iedges[4][1] = 11;
    iedges[5][0] = 1;  iedges[5][1] = 3;
    iedges[6][0] = 1;  iedges[6][1] = 4;
    iedges[7][0] = 1;  iedges[7][1] = 7;
    iedges[8][0] = 1;  iedges[8][1] = 11;
    iedges[9][0] = 2;  iedges[9][1] = 5;
    iedges[10][0] = 2; iedges[10][1] = 6;
    iedges[11][0] = 2; iedges[11][1] = 8;
    iedges[12][0] = 2; iedges[12][1] = 11;
    iedges[13][0] = 3; iedges[13][1] = 5;
    iedges[14][0] = 3; iedges[14][1] = 7;
    iedges[15][0] = 3; iedges[15][1] = 9;
    iedges[16][0] = 3; iedges[16][1] = 11;
    iedges[17][0] = 4; iedges[17][1] = 6;
    iedges[18][0] = 4; iedges[18][1] = 7;
    iedges[19][0] = 4; iedges[19][1] = 10;
    iedges[20][0] = 5; iedges[20][1] = 8;
    iedges[21][0] = 5; iedges[21][1] = 9;
    iedges[22][0] = 5; iedges[22][1] = 11;
    iedges[23][0] = 6; iedges[23][1] = 8;
    iedges[24][0] = 6; iedges[24][1] = 10;
    iedges[25][0] = 7; iedges[25][1] = 9;
    iedges[26][0] = 7; iedges[26][1] = 10;
    iedges[27][0] = 8; iedges[27][1] = 9;
    iedges[28][0] = 8; iedges[28][1] = 10;
    iedges[29][0] = 9; iedges[29][1] = 10;
    
    ifaces[0][0] = 0;  ifaces[0][1] = 1;  ifaces[0][2] = 4;
    ifaces[1][0] = 0;  ifaces[1][1] = 1;  ifaces[1][2] = 11;
    ifaces[2][0] = 0;  ifaces[2][1] = 2;  ifaces[2][2] = 6;
    ifaces[3][0] = 0;  ifaces[3][1] = 2;  ifaces[3][2] = 11;
    ifaces[4][0] = 0;  ifaces[4][1] = 4;  ifaces[4][2] = 6;
    ifaces[5][0] = 1;  ifaces[5][1] = 3;  ifaces[5][2] = 7;
    ifaces[6][0] = 1;  ifaces[6][1] = 3;  ifaces[6][2] = 11;
    ifaces[7][0] = 1;  ifaces[7][1] = 4;  ifaces[7][2] = 7;
    ifaces[8][0] = 2;  ifaces[8][1] = 5;  ifaces[8][2] = 8;
    ifaces[9][0] = 2;  ifaces[9][1] = 5;  ifaces[9][2] = 11;
    ifaces[10][0] = 2; ifaces[10][1] = 6; ifaces[10][2] = 8;
    ifaces[11][0] = 3; ifaces[11][1] = 5; ifaces[11][2] = 9;
    ifaces[12][0] = 3; ifaces[12][1] = 5; ifaces[12][2] = 11;
    ifaces[13][0] = 3; ifaces[13][1] = 7; ifaces[13][2] = 9;
    ifaces[14][0] = 4; ifaces[14][1] = 6; ifaces[14][2] = 10;
    ifaces[15][0] = 4; ifaces[15][1] = 7; ifaces[15][2] = 10;
    ifaces[16][0] = 5; ifaces[16][1] = 8; ifaces[16][2] = 9;
    ifaces[17][0] = 6; ifaces[17][1] = 8; ifaces[17][2] = 10;
    ifaces[18][0] = 7; ifaces[18][1] = 9; ifaces[18][2] = 10;
    ifaces[19][0] = 8; ifaces[19][1] = 9; ifaces[19][2] = 10;
    
    return;
}


void uvec(vec *iverts, int *thisface, vec bcen, vec v)
{
    int i, j;
    
    for (i=0;i<3;i++)
    {
        v[i] = 0.0;
        for (j=0;j<3;j++)
        {
            v[i] += bcen[j] * iverts[thisface[j]][i];
        }
    }
    
    return;
}


void anginterp(vec vec1, vec vec2, double a, vec p)

{
    int i, j;
    vec ax, R[3];
    double angle;
    
    CROSS(vec1, vec2, ax);
    NORMAL(ax);
    
    angle = 0.0;
    for (i=0; i<3; i++)
    {
        angle += vec1[i] * vec2[i];
        p[i] = vec1[i];
    }
    angle = a * acos(angle);
    
    if (fabs(angle-M_PI) < ZERO)
    {
        std::cerr << "Warning: your matrix is not uniquely defined when start and end points are antipodal." << std::endl;
    }
    
    AA_MAT(ax,angle,R);
    
    ROTATE(R,p);
    NORMAL(p);
    
    return;
}

void sample_icosa(int L)
{
    int np = num_icosamp(L);
    P = Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(np,3);
    int pi = 0;
    
    std::cerr << "Creating " << np << " points on S2 by subdividing Icosahedron " << L << " times..." << std::endl;
    auto tp0 = std::chrono::steady_clock::now();

    
    int i, j, k, ii, jj, kk;
    int T;
    
    vec  iverts[12], vert1, vert2, bcen;
    sim  ifaces[20];
    edge iedges[30];
    
    double w, a, step, twist;
    vec p, p1, p2, R[3];
    int flat = 0;
    
    /* define icosahedron */
    set_icosa(iverts, iedges, ifaces);
    
    /* loop over all vertices */
    for (i=0; i<12; i++)
    {
        P.row(pi) << iverts[i][0], iverts[i][1], iverts[i][2];
        pi++;
    }
    
    /* loop over all edges */
    for (i=0; i<30; i++)
    {
        /* calculate subdivisions */
        for (j=1; j<L+1; j++)
        {
            /* take a linear combination of the two endpoints */
            
            a =  (double) j/L;
            
            bcen[0] = a;
            bcen[1] = 1.0-a;
            bcen[2] = 0.0;
            
            if (!flat)
            {
                anginterp(iverts[iedges[i][0]],iverts[iedges[i][1]],a,p);   /* note: this could be made more efficient */
            }
            else
            {
                for (k=0; k<3; k++)
                {
                    p[k] = a * iverts[iedges[i][0]][k] + (1.0-a) * iverts[iedges[i][1]][k];
                }
            }
            
            P.row(pi) << p[0],p[1],p[2];
            pi++;
        }
    }
    
    /* loop over all simplices */
    for (i=0;i<20;i++)
    {
        
       
        /* step along simplex edge using edge barystep j/L */
        for (j=1;j<L;j++)
        {
            a = (double) j/L;
            
            if (!flat)
            {
                anginterp(iverts[ifaces[i][1]],iverts[ifaces[i][0]],a,p1);
                anginterp(iverts[ifaces[i][2]],iverts[ifaces[i][0]],a,p2);
                for (k=1;k<L-j+1;k++)
                {
                    w = (double) k / (L-j);
                    anginterp(p1,p2,w,p);  /* note: this could be made more efficient */
                    P.row(pi) << p[0],p[1],p[2];
                    pi++;
                }
            }
            else
            {
                vert1[0] = 0.0;
                vert1[1] = a;
                vert1[2] = 1.0 - a;
                
                vert2[0] = 1.0 - a;
                vert2[1] = a;
                vert2[2] = 0.0;
                
                /* step across simplex interior using edge barystep k/(L-j) */
                for (k=1;k<L-j;k++)
                {
                    w = (double) k / (L-j);
                    
                    /* define Barycentric coords relative to simplex */
                    for (ii=0; ii<3; ii++)
                    {
                        bcen[ii] = w * vert1[ii] + (1.0 - w) * vert2[ii];
                    }
                    
                    for (jj=0; jj<3; jj++)
                    {
                        p[jj] = 0.0;
                        for (kk=0;kk<3;kk++)
                        {
                            p[jj] += bcen[kk] * iverts[ifaces[i][kk]][jj];
                        }
                    }
                    /* print out results */
                    P.row(pi) << p[0],p[1],p[2];
                    pi++;
                }
            }
        }
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    tt += (std::chrono::duration <double, std::milli> (tp1-tp0).count());
    std::cerr << pi << std::endl;

}

void sample_sphere(int N)
{
    
    std::cerr << "Creating " << N << " points on S2 by repulsion..." << std::endl;
    auto tp0 = std::chrono::steady_clock::now();

    int i,k;
    double l, e, e0, d;
    vec *p0, *p1, *f, *pp0, *pp1;
    double step=0.01;

   
    p0 = (vec *) calloc(N,sizeof(vec));
    p1 = (vec *) calloc(N,sizeof(vec));
    f = (vec *) calloc(N,sizeof(vec));
    pp0 = p0;
    pp1 = p1;

    std::default_random_engine rng;
//    std::normal_distribution<double> dn;
    std::uniform_real_distribution<double> du(-1.0,1.0);

    for( i = 0; i<N; i++ )
    {
        p0[i][0] = du(rng);
        p0[i][1] = du(rng);
        p0[i][2] = du(rng);
        l = length(p0[i]);
        if(l!=0.0)
        {
            p0[i][0] /= l;
            p0[i][1] /= l;
            p0[i][2] /= l;
        }
        else
            i--;
    }

    e0 = get_coulomb_energy(N,p0);
    for(k = 0;k<Nstep;k++)
    {
        get_forces(N,f,p0);
        for(i=0; i < N;i++)
        {
            d = dot(f[i],pp0[i]);
            f[i][0]  -= pp0[i][0]*d;
            f[i][1]  -= pp0[i][1]*d;
            f[i][2]  -= pp0[i][2]*d;
            pp1[i][0] = pp0[i][0]+f[i][0]*step;
            pp1[i][1] = pp0[i][1]+f[i][1]*step;
            pp1[i][2] = pp0[i][2]+f[i][2]*step;
            l = length(pp1[i]);
            pp1[i][0] /= l;
            pp1[i][1] /= l;
            pp1[i][2] /= l;
        }
        e = get_coulomb_energy(N,pp1);
        if(e >= e0)
        {
            /* not successful step */
            step /= 2;
            if(step < minimal_step)
                break;
            continue;
        }
        else
        {
            /* successful step */
            vec *t = pp0;      pp0 = pp1; pp1 = t;
            e0 = e;
            step*=2;
        }
    }

    P = Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(N,3);

    for (i = 0;i < N;i++)
    {
        P.row(i) << p0[i][0],p0[i][1],p0[i][2];
    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    tt += (std::chrono::duration <double, std::milli> (tp1-tp0).count());

}
    

Eigen::Matrix<double,Eigen::Dynamic,4>
sample_SOI(int n)
{
    int i = 0;
    double  s, w[3];        /* current point from vector list */
    double  alpha;          /* increment angle for rotation sample */
    double  angles[3];      /* ordered angles */

    alpha = 2.0 * M_PI/((double) n);

    int np = P.rows();
    int nq = np*n;
    Eigen::Matrix<double,Eigen::Dynamic,4> V = Eigen::Matrix<double,Eigen::Dynamic,4>::Zero(nq,4);
    int qi = 0;

    std::cerr << "Creating " << nq << " points on S3 by SOI from " << np << " points on S2..." << std::endl;
    auto tp0 = std::chrono::steady_clock::now();

    
    for (int pi = 0; pi < np; pi++)
    {
        w[0] = P(pi,0); w[1] = P(pi,1); w[2] = P(pi,2);
        
        if (w[2] > 1.0)
        {
            w[2] = 1.0;
        }
        
        if (w[2] < -1.0)
        {
            w[2] = -1.0;
        }
        
        angles[1] = acos(w[2]);
        
        // printf("%lf %lf %lf %lf\n",w[0],w[1],w[2],angles[1]);
        if (angles[1] == 0.0)
        {
            angles[2] = 0.0;
        }
        else
        {
            s = -w[0]/sin(angles[1]);
            if (s > 1.0)
            {
                s = 1.0;
            }
            
            if (s < -1.0)
            {
                s = -1.0;
            }
            
            angles[2] = acos(s);
            
            if(w[1] < 0)
            {
                angles[2] = 2.0 * M_PI - angles[2];
            }
        }
        if (angles[2]  >  M_PI)
        {
            angles[2] += - 2.0 * M_PI;
        }
        
        for (i=0; i<n; i++ )
        {
            angles[0] = - M_PI + i * alpha;
            
            double M[9];
            EA_MAT(M,angles);
            Eigen::Map<Eigen::Matrix<double,3,3> > R(M);
            Eigen::Quaterniond q(R);
            V.row(qi) << q.x(),q.y(),q.z(),q.w();
            qi++;
        }

    }
    auto tp1 = std::chrono::steady_clock::now();
    std::cerr << "Done in " << (std::chrono::duration <double, std::milli> (tp1-tp0).count()) << " ms" << std::endl;
    tt += (std::chrono::duration <double, std::milli> (tp1-tp0).count());

    return V;
}
    
int main(int argc, char *argv[])
{
    CLI::App app{"SOI-Sampling the 3-sphere"};

    
    int seed = 0;
    app.add_option("-s,--seed", seed, "random number seed");

    int n = 1;
    app.add_option("-n", n, "number of points");
    
    bool repulsion = false;
    app.add_flag("-r,--repulsion", repulsion, "Use repulsion on S2");
    
    std::string filename = "default";
    CLI::Option *fileopt = app.add_option("-f,--filename", filename, "Point set filename");

    bool sout = false;
    app.add_flag("-o",sout,"Output to cout");

   
    CLI11_PARSE(app, argc, argv);

    Eigen::Matrix<double,Eigen::Dynamic,4> V;
    std::cerr << "Trying to aample " << n << " points on S3" << std::endl;

    
    //int i = 0;

    tt = 0.0;
    if (repulsion)
    {
        int l = 0;
        while (l*num_unifsamp(l) < n)
            l++;
        
        sample_sphere(num_unifsamp(l));
        V = sample_SOI(l);
    }
    else
    {
        int l = 0;
        while ( ((int) floor(5.675 * (double) (l+1))) * num_icosamp(l) < n )
            l++;
  
        sample_icosa(l);
        V = sample_SOI((int) floor(5.675 * (double) (l+1)));
    }
    if (!sout)
    {
        std::cerr << "Created ";
        std::cout << V.rows() << " ";
        std::cerr << "points on S3 in ";
        std::cout << tt;
        std::cerr << " ms";
        std::cout << std::endl;
    }
    else
        std::cerr << "Created " << V.rows() << " points on S3 in " << tt << " ms" << std::endl;

    
    
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
