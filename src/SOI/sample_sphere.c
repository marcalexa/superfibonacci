/* Author:              V.Bulatov@ic.ac.uk */
/* Revised:             Julie C. Mitchell  */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "rotutils.h"

double    frand(void);
double    dot(vec v1,vec v2);
double    length(vec v);
double    vecdist(vec v1,vec v2);
double    get_coulomb_energy(int N,vec p[]);
void      get_forces(int N,vec f[], vec p[]);

int main(int argc, char *argv[])
{
    static int N=100,Nstep=1000;
    static double step=0.01;
    static double minimal_step=1.e-10;
    int i,k;
    double l, e, e0, d;
    vec *p0, *p1, *f, *pp0, *pp1;
    char inBuf[80], *str;

    if(argc < 2)
    {
        fprintf(stderr,"Usage:  sample_sphere   number_of_points   [maximal_number_of_steps]   >  outfile \n");
        exit(1);
    }

    N = atoi(argv[1]);

    if(argc > 2)
	{
    	Nstep = atoi(argv[2]);
	}
	
    p0 = (vec *) calloc(N,sizeof(vec));
    p1 = (vec *) calloc(N,sizeof(vec));
    f = (vec *) calloc(N,sizeof(vec));
    pp0 = p0;
    pp1 = p1;
	
    srand(time(NULL));
	
    for(i = 0; i<N; i++ )
    {
        p0[i][0] = 2*frand();
        p0[i][1] = 2*frand();
        p0[i][2] = 2*frand();
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
        /* fprintf(stderr,"\rn: %5d, e = %18.8lf step = %12.10lf",k,e,step); */
        fflush(stderr);
    }
    for (i = 0;i < N;i++)
    {
        printf("%15.10lf %15.10lf %15.10lf\n",p0[i][0],p0[i][1],p0[i][2]);
    }
    return 0;
}
	
double frand(void)
{
    return ((rand()-(RAND_MAX/2))/(RAND_MAX/2.));
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
