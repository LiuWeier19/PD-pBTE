#include<iostream>
#include <iomanip>
#include<fstream>
#include<math.h>
#include<cstdio>
#include<cmath>
#include<string>
#include<cstring>
#include<sstream>
#include <stdio.h>
#include "omp.h"
using namespace std;

double max(double a,double b,double c,double d){
	a=max(a,b);
	c=max(c,d);
	return max(a,c);
} 
double min_d(double a,double b,double c,double d){
	a=min(a,b);
	c=min(c,d);
	return min(a,c);
} 
int min_i(int a,int b,int c,int d){
  int e=a;
	//e=min(a,b);
	//f=min(c,d);
	//return min(e,f);
  if(b<e){
    e=b;
  }
  if(c<e){
    e=c;
  }
  if(d<e){
    e=d;
  }
  return e;
}
double Trape(double *a,int n,double delta,double *theta){
	double sum=a[0]*theta[0]+a[n]*theta[n];
	for(int i=1;i<n;i++){
		sum+=2*a[i]*theta[i];
	}
	sum=sum*delta/2;
	return sum;	
}
double Trape(double *a,int n,double delta){
	double sum=a[0]+a[n];
	for(int i=1;i<n;i++){
		sum+=2*a[i];
	}
	sum=sum*delta/2;
	return sum;	
}

double simpson(double *a,int n,double delta,double *theta){ //simpson rule to get the intergration of  a*sin(theta) 
	double sum=a[0]*sin(theta[0])+a[n]*sin(theta[n]);
	for(int i=1;i<=n/2.;i++){
		if(i<n/2.)
            sum+=2*a[2*i]*sin(theta[2*i]);
        sum+=4*a[2*i-1]*sin(theta[2*i-1]);	
	}
	sum=delta*sum/3.;		
return sum;			
} 
double simpson_2(double *a,int n,double delta,double *sintheta){ //simpson rule to get the intergration of  a*sin(theta) 
	double sum=a[0]*sintheta[0]+a[n]*sintheta[n];
	for(int i=1;i<=n/2.;i++){
		if(i<n/2.)
            sum+=2*a[2*i]*sintheta[2*i];
        sum+=4*a[2*i-1]*sintheta[2*i-1];	
	}
	sum=delta*sum/3.;		
return sum;			
}
double simpson(double *a,int n,double delta){ //simpson rule to get the intergration  
	double sum=a[0]+a[n];
	for(int i=1;i<=n/2.;i++){
		if(i<n/2.)
            sum+=2*a[2*i];
        sum+=4*a[2*i-1];	
	}
	sum=delta*sum/3.;		
	return sum;			
} 
double simpson_3(double *a,int n,double delta){ //simpson rule to get the intergration  
	double sum=a[0]+a[n];
    double nn2=n/2.;
    int i;
	for(int i=1;i<(n>>1);i++){
		//if(i<n/2.)
            sum+=2*a[i<<1];
        sum+=4*a[(i<<1)-1];	
	}
    if(i==nn2){
        sum+=4*a[(i<<1)-1];
    }
	sum=delta*sum/3.;		
	return sum;			
}

double simpson_I0(double *a[],double *theta,double delta1,double delta2,int n1, int n2){  //I0
	double sum1[n2+1],sum2;
    for(int i=0;i<n2+1;i++){
    	sum1[i]=simpson(a[i],n1,delta1,theta);
    } 
   sum2=simpson(sum1,n2,delta2);   		
   return sum2;	
}
double simpson_I0_2(double *a[],double *theta,double delta1,double delta2,int n1, int n2){  //I0
	double sum1[n2+1],sum2;
    for(int i=0;i<n2+1;i++){
    	sum1[i]=simpson_2(a[i],n1,delta1,theta);
    } 
   sum2=simpson(sum1,n2,delta2);   		
   return sum2;	
}
double qx(double *I[],double *theta,double delta1,double delta2,int n1, int n2){
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*cos(theta[0])*sin(theta[0])+I[i][n1]*cos(theta[n1])*sin(theta[n1]);
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*cos(theta[2*j])*sin(theta[2*j]);
      sum1[i]+=4*I[i][2*j-1]*cos(theta[2*j-1])*sin(theta[2*j-1]);
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
}
double qx_2(double *I[],double *sintheta,double *costheta,double delta1,double delta2,int n1, int n2){
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*costheta[0]*sintheta[0]+I[i][n1]*costheta[n1]*sintheta[n1];
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*costheta[2*j]*sintheta[2*j];
      sum1[i]+=4*I[i][2*j-1]*costheta[2*j-1]*sintheta[2*j-1];
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
}
double qx_3(double *I[],double *sincostheta,double delta1,double delta2,int n1, int n2){
  //cout<<1<<endl;
  
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*sincostheta[0]+I[i][n1]*sincostheta[n1];
  //cout<<sum1[i]<<" , i:"<<i<<endl;
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*sincostheta[2*j];
      sum1[i]+=4*I[i][2*j-1]*sincostheta[2*j-1];
      //cout<<sum1[i]<<" , i:"<<i<<" , j:"<<j<<endl;
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
  //return 1.;
}

double qz(double *I[],double *theta,double *phi,double delta1,double delta2,int n1, int n2){
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*sin(theta[0])*sin(theta[0])*sin(phi[i])+I[i][n1]*sin(theta[n1])*sin(theta[n1])*sin(phi[i]);
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*sin(theta[2*j])*sin(theta[2*j])*sin(phi[i]);
      sum1[i]+=4*I[i][2*j-1]*sin(theta[2*j-1])*sin(theta[2*j-1])*sin(phi[i]);
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
}
double qz_2(double *I[],double *sintheta,double *sinphi,double delta1,double delta2,int n1, int n2){
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*sintheta[0]*sintheta[0]*sinphi[i]+I[i][n1]*sintheta[n1]*sintheta[n1]*sinphi[i];
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*sintheta[2*j]*sintheta[2*j]*sinphi[i];
      sum1[i]+=4*I[i][2*j-1]*sintheta[2*j-1]*sintheta[2*j-1]*sinphi[i];
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
}

double qy_2(double *I[],double *sintheta,double *sinphi,double delta1,double delta2,int n1, int n2){
  double sum1[n2+1],sum2;
  for(int i=0;i<n2+1;i++){
  sum1[i]=I[i][0]*sintheta[0]*sintheta[0]*sinphi[i]+I[i][n1]*sintheta[n1]*sintheta[n1]*sinphi[i];
    for(int j=1;j<=n1/2.;j++){
      if(j<n1/2.)
        sum1[i]+=2*I[i][2*j]*sintheta[2*j]*sintheta[2*j]*sinphi[i];
      sum1[i]+=4*I[i][2*j-1]*sintheta[2*j-1]*sintheta[2*j-1]*sinphi[i];
    }
    sum1[i]=sum1[i]*delta1/3.;
  }
  sum2=simpson(sum1,n2,delta2);
  return sum2;
}

double GL2(double *a,double minx,double maxx,int n){
    double sum=0,aa=0,bb=0;
    double deltaGL=(maxx-minx)/(double)n;
    double deltaGL2=deltaGL/2.;
    for(int i=0;i<n;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
        sum+=(a[i+i]+a[i+i+1])*deltaGL2;
    }

    return sum;
}

double GL2(double *a,double minx,double maxx,int n,double *sintheta){
    double sum=0,aa=0,bb=0;
    double deltaGL=(maxx-minx)/(double)n;
    double deltaGL2=deltaGL/2.;
    for(int i=0;i<n;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
        sum+=(a[i+i]*sintheta[i+i]+a[i+i+1]*sintheta[i+i+1])*deltaGL2;
    }

    return sum;
}

double GL3(double *a,double minx,double maxx,int n){
    double sum=0,aa=0,bb=0;
    double deltaGL=(maxx-minx)/(double)n;
    double deltaGL2=deltaGL/2.;
    for(int i=0;i<n;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
        sum+=(0.55555556*a[i*3]+0.88888889*a[i*3+1]+0.55555556*a[i*3+2])*deltaGL2;
    }

    return sum;
}

double GL4(double *a,double minx,double maxx,int n){
    double sum=0,aa=0,bb=0;
    double deltaGL=(maxx-minx)/(double)n;
    double deltaGL2=deltaGL/2.;
    for(int i=0;i<n;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
        sum+=(0.3478548451*a[i*4]+0.6521451549*a[i*4+1]+0.6521451549*a[i*4+2]+0.3478548451*a[i*4+3])*deltaGL2;
    }

    return sum;
}

double GL4(double *a,double minx,double maxx,int n,double *sintheta){
    double sum=0,aa=0,bb=0;
    double deltaGL=(maxx-minx)/(double)n;
    double deltaGL2=deltaGL/2.;
    for(int i=0;i<n;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
        sum+=(0.3478548451*a[i*4]*sintheta[4*i]+0.6521451549*a[i*4+1]*sintheta[4*i+1]+0.6521451549*a[i*4+2]*sintheta[4*i+2]+0.3478548451*a[i*4+3]*sintheta[4*i+3])*deltaGL2;
    }

    return sum;
}

double GL2_I0(double *a[],double *theta,double minx1,double maxx1,int n1,double minx2,double maxx2,int n2){  //I0
	double sum1[2*n2],sum2;
    for(int i=0;i<2*n2;i++){
    	sum1[i]=GL2(a[i],minx1,maxx1,n1,theta);
    } 
   sum2=GL2(sum1,minx2,maxx2,n2);   		
   return sum2;	
} 

double GL2_q_2(double *a[],double *theta,double *phi,double minx1,double maxx1,int n1,double minx2,double maxx2,int n2){  //I0
	double sum1[2*n2],sum2;
    for(int i=0;i<2*n2;i++){
    	sum1[i]=GL2(a[i],minx1,maxx1,n1,theta)*phi[i];
    } 
   sum2=GL2(sum1,minx2,maxx2,n2);   		
   return sum2;	
} 

double GL4_I0(double *a[],double *theta,double minx1,double maxx1,int n1,double minx2,double maxx2,int n2){  //I0
	double sum1[4*n2],sum2;
    for(int i=0;i<4*n2;i++){
    	sum1[i]=GL4(a[i],minx1,maxx1,n1,theta);
    } 
   sum2=GL4(sum1,minx2,maxx2,n2);   		
   return sum2;	
} 

double GL4_q_2(double *a[],double *theta,double *phi,double minx1,double maxx1,int n1,double minx2,double maxx2,int n2){  //I0
	double sum1[4*n2],sum2;
    for(int i=0;i<4*n2;i++){
    	sum1[i]=GL4(a[i],minx1,maxx1,n1,theta)*phi[i];
    } 
   sum2=GL4(sum1,minx2,maxx2,n2);   		
   return sum2;	
} 

double deter33(double D[3][3]){
  double s=D[0][0]*D[1][1]*D[2][2]+D[0][1]*D[1][2]*D[2][0]+D[0][2]*D[1][0]*D[2][1]-D[0][0]*D[1][2]*D[2][1]-D[0][1]*D[1][0]*D[2][2]-D[0][2]*D[1][1]*D[2][0];
  return s;
}

double deter22(double D[2][2]){
  double s=D[0][0]*D[1][1]-D[0][1]*D[1][0];
  return s;
}