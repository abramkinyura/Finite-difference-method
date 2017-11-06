#include <stdio.h>
#include <math.h>
#include <stdlib.h>



int main(void)
{
int N;

double dist;
double alphacoeff;
double *a;
double *b;
double *c;
double *f;
double *u1;
double *u2;
double *alpha;
double *beta;
int i;
printf("Input N\n");
scanf("%d",&N);
printf("Input alpha\n");
scanf("%lf",&alphacoeff);
//N=N-1;
double h=1.0/(double)(N);
 a = (double*)malloc((N+1) * sizeof(double));
 b = (double*)malloc((N+1) * sizeof(double));
 c = (double*)malloc((N+1) * sizeof(double));
 f = (double*)malloc((N+1) * sizeof(double));
 u1 = (double*)malloc((N+1) * sizeof(double)); 
 u2 = (double*)malloc((N+1) * sizeof(double));
 alpha=(double*)malloc((N+1)*sizeof(double));
 beta=(double*)malloc((N+1)*sizeof(double));

FILE *out=fopen("data.txt","w");

for (i=0;i<=N;i++)
u1[i]=0.0;

for (i=1;i<=(N-1);i++)
b[i]=1.0/(h*h);

for (i=1;i<=(N-1);i++)
c[i]=2.0/(h*h);

for (i=1;i<=(N-1);i++)
a[i]=1.0/(h*h);





dist=1.0;

while (dist>=1e-6)

{
dist=0.0;
for (i=1;i<=(N-1);i++)
f[i]=alphacoeff*exp(i*h)-u1[i]/(1+fabs(u1[i])); 

alpha[1]=0;
beta[1]=0;

for (i=1;i<=(N-1);i++)
alpha[i+1]=a[i]/(c[i]-b[i]*alpha[i]);

for (i=1;i<=(N-1);i++)
beta[i+1]=(f[i]+b[i]*beta[i])/(c[i]-b[i]*alpha[i]);

u2[N]=0;

for(i=(N-1);i>=0;i--)
u2[i]=alpha[i+1]*u2[i+1]+beta[i+1];

for(i=0;i<=N;i++)
printf("%lf\n",u2[i]);

for(i=0;i<=N;i++)
dist+=h*fabs(u1[i]-u2[i])*fabs(u1[i]-u2[i]);

printf("dist=%lf\n",dist);

printf("\n");

for(i=0;i<=N;i++)
u1[i]=u2[i];
}

for (i=0;i<=N;i++)
fprintf(out,"%lf\t%lf\n",i*h,u2[i]);


//printf("dist=%lf\n",dist);


return 0;
    }	








