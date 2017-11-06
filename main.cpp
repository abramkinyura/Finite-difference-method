include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>

using namespace std;

double nextu(double u1, double u2, double fn, double h, double a)
{
	return 2*u2 - u1 + a*h*h*u2*u2*u2 - fn*h*h;
}

double uN(double u1, double* u, double* f, int N, double a)
{
	u[1]=u1;
	double h=1.0/(double)N;
	for(int i=2; i<=N; ++i){ u[i]=nextu(u[i-2],u[i-1],f[i-1],h,a); if(fabs(u[i])>1e50) return u[i];}
	return u[N];
}

double bisec(double left, double right, double* u, double* f, int N, double a, int iter)
{
	double u1=(left+right)/2;
	double u_N=uN(u1,u,f,N,a);
	if(fabs(u_N-1)<1e-10)
	{
		cout<<"iter="<<iter<<endl;
		 return u1;
	}
	double u_Nl=uN(left,u,f,N,a);
	//double u_Nr=uN(right,u,f,N,a);
	//cout<<"left="<<left<<" "<<u_Nl<<" right="<<right<<" "<<u_Nr<<" u1="<<u1<<" "<<u_N<<endl;
	if((u_N-1)*(u_Nl-1)<0) return bisec(left,u1,u,f,N,a,iter+1);
	else return bisec(u1,right,u,f,N,a,iter+1);
}
int main()
{
	double a;
	double left, right;
	ofstream fout("data.txt");
	if(!fout)
	{
	    cout<<"Can't open file\n";
	     return 1;
	}
	int varf;
	int N,N0;
//	double h;
	cout<<"N=";
	cin>>N;
	N0=N;
	double h0=(double)1/(double)N0;
	double* u0=(double*)malloc((N+1)*sizeof(double));
//	h=1.0/(double)N;
	cout<<"f(x)=sin(x) (1) or sin(x^2) (2) ?\n";
	cin>>varf;
	cout<<"a=";
	cin>>a;
	int coef=1;
	for(int j=0; j<5; ++j)
	{
	  cout<<"N="<<N<<endl;
	double* f=(double*)malloc((N+1)*sizeof(double));
	double* u=(double*)malloc((N+1)*sizeof(double));
	for(int i=0; i<=N; ++i)
	{
		double x=(double)i/double(N);
		switch(varf)
		{
			case 1:
			f[i]=sin(x);
			break;
			case 2:
			f[i]=sin(x*x);
			break;
			default:
			cout<<"Error\n";
			return 1;
		}
	}
	u[0]=0;
	bool a1,a2;
	a1=a2=false;
	for(;;)
	{
		double u1;
		cin>>u1;
		double u_N=uN(u1,u,f,N,a);
		if(fabs(u_N)>1e50) {cout<<"Too big value"<<endl; continue;}
		cout<<u_N<<endl;
		
		if(!a1 && uN(u1,u,f,N,a)<1) 
		{
			left=u1;
			a1=true;
		}
		if(!a2 && uN(u1,u,f,N,a)>1)
		{
			right=u1;
			a2=true;
		}
		if(a1&&a2)
		{
			cout<<"u1^(0)="<<left<<" uN^(0)="<<uN(left,u,f,N,a)<<endl<<"u1^(1)="<<right<<" uN^(1)="<<uN(right,u,f,N,a)<<endl;
			break;
		}
	}
	double u1=bisec(left, right, u, f, N, a, 1);
	u[N]=uN(u1,u,f,N,a);
	cout<<fabs(u[N]-1)<<endl;
	for(int i=0; i<=N; ++i)
	{
	  double x=(double)i/(double)N;
	   fout<<x<<"\t"<<u[i]<<endl;
	}
	if(j==0)for(int k=0; k<N0+1; ++k) u0[k]=u[coef*k];
	if(j>0){for(int k=0; k<N0+1; ++k)
	{
	    cout<<"N="<<N/2<<" u("<<(double)k*h0<<")~"<<u0[k]<<endl;
	    cout<<"N="<<N<<" u("<<(double)k*h0<<")~"<<u[coef*k]<<endl;
	    cout<<"res="<<fabs(u0[k]-u[coef*k])<<endl;
	}}
	for(int k=0; k<N0+1; ++k) u0[k]=u[coef*k];
	N=2*N;
	coef*=2;
	free(u);
	free(f);
	}
	free(u0);
	
}
	
			

