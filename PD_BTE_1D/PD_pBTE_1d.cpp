//1d
#include"../calFunction.h"
//the global variable
const double kb = 1.3806488e-23;// Planck constant, J s.
const double h = 6.62606957e-34;// Dirac constant, J s.
const double hbar = 1.054571726e-34;
const double pi =3.1415926535897932384626433832795028841953;
const double Tl=301.,Tr=300.,C=0.930e6,v=1804.;// Tl-left temperature,Tr-right temperature,v-velocity,C-heat capacity
const int nx=100.,nu=25.; //nx-x parts,nu-Number of angle integration intervals
const double ll=260.4e-9;//mean free path of Si
const double ll_1=1./ll;
const double deltau=1./nu;
const double kn=1000.;
const int gln=4;//The number of points of Gauss-Legendre integral
const int glnu=gln*nu;


double eprt_pddo(int x,int n,double u2,double *I0,double deltax,double L,double I[nx+1][nu+1]){
    double If1;
    if(x==0){
        If1=C*v*Tl/4./pi;
    }else if(x==nx){
        If1=C*v*Tr/4./pi;
    }else{
        If1=I0[x]-ll*u2*0.5*(1/deltax)*(I[x+1][n]-I[x-1][n]);
    }
    return If1;
}

void out_file(string file,double *T,double *x,int nx,double *q,double *k){   //output file 
    ofstream outfile;
    outfile.open(file);
    outfile<<"x,"<<"T*,"<<"q,"<<"q*"<<endl;
	for(int i=0;i<nx+1;i++){
		outfile<<x[i]<<","<<T[i]<<","<<q[i]<<","<<k[i]<<endl;
	}
	outfile.close();
}

void makeparaf(double u[], double xl[nx+1], double xlong[nx+1],double deltaxh, double** para1, double** para_d1, int** HN){
	for(int x=0;x<nx+1;x++){
		for(int j=0;j<glnu;j++){
			if(x==0){
				//para1[x][j]=0;
				//HN[x][j]=0;
			}else if(xl[x]-xl[0]<deltaxh){
				double sum1=0,sum2=0;
				int i=1;
				while(i<=x){
					sum1+=xlong[x-i];
					sum2+=u[j]*xlong[x-i]/(xl[x-i]-xl[x]);
					i++;
				}

				sum2=(sum2/sum1)-1./ll;
				para_d1[x][j]=1./sum1;
				para1[x][j]=1./sum2;
				HN[x][j]=i-1;

				//cout<<HN[x][j]<<endl;

			}else{
				double sum1=0,sum2=0;
				int i=1;
				while(((xl[x]-xl[x-i])<=deltaxh)&&((x-i)>=0)){
					sum1+=xlong[x-i];
					sum2+=u[j]*xlong[x-i]/(xl[x-i]-xl[x]);
					i++;
				}

				sum2=(sum2/sum1)-1./ll;
				para_d1[x][j]=1./sum1;
				para1[x][j]=1./sum2;
				HN[x][j]=i-1;
			}
		}
	}
}

void makeparab(double u[], double xl[nx+1], double xlong[nx+1],double deltaxh, double** para1, double** para_d1, int** HN){
	for(int x=0;x<nx+1;x++){
		for(int j=0;j<glnu;j++){
			if(x==nx){
				//para1[x][j]=0;
				//HN[x][j]=0;
			}else if(xl[nx]-xl[x]<deltaxh){
				double sum1=0,sum2=0;
				int i=1;
				while(i<=(nx-x)){
					sum1+=xlong[x+i];
					sum2+=u[j]*xlong[x+i]/(xl[x+i]-xl[x]);
					i++;
				}
				sum2=(sum2/sum1)-1./ll;
				para_d1[x][j]=1./sum1;
				para1[x][j]=1./sum2;
				HN[x][j]=i-1;
			}else{
				double sum1=0,sum2=0;
				int i=1;
				while(((xl[x+i]-xl[x])<=deltaxh)&&(x+i<=nx)){
					sum1+=xlong[x+i];
					sum2+=u[j]*xlong[x+i]/(xl[x+i]-xl[x]);
					i++;
				}
				sum2=(sum2/sum1)-1./ll;
				para_d1[x][j]=1./sum1;
				para1[x][j]=1./sum2;
				HN[x][j]=i-1;
			}
		}
	}
}

double eprt_pddo_f(int x,int n,double u2,double *I0_ll,double deltax,double L,double** I,double xl[nx+1],double xlong[nx+1],double deltaxh, double para1,double para_d1,int HN){
    double If1;
	if(x==0){
		If1=C*v*Tl/4./pi;
	}else{
		double sum1=0;//sum1-If
		int i=1;
		while(i<=HN){
			sum1+=I[x-i][n]*u2*xlong[x-i]/(xl[x-i]-xl[x]);
			i++;
		}
		sum1=sum1*para_d1-I0_ll[x];

		If1=sum1*para1;
	}
    return If1;
}

double eprt_pddo_b(int x,int n,double u2,double *I0_ll,double deltax,double L,double** I,double xl[nx+1],double xlong[nx+1],double deltaxh, double para1,double para_d1,int HN){
    double If1;
	if(x==nx){
		If1=C*v*Tr/4./pi;
	}else{
		double sum1=0;//sum1-If
		int i=1;
		while(i<=HN){
			sum1+=I[x+i][n]*u2*xlong[x+i]/(xl[x+i]-xl[x]);
			i++;
		}
		sum1=sum1*para_d1-I0_ll[x];

		If1=sum1*para1;
	}

    return If1;
}



int main(){
    clock_t start,finish;
    double totaltime;
    start=clock();  //timing
    double beta;
	double u1[glnu]={-1},u2[glnu]={0},u[glnu],x[nx+1],xlong[nx+1];
    double q[nx+1]={0},k[nx+1]={0},qf[nx+1],qb[nx+1],k2[nx+1]={0};
    double sum1[nx+1]={0},sum2[nx+1]={0},I0[nx+1],T[nx+1],T2[nx+1],I0_ll[nx+1],I00[nx+1];
    double L=ll/kn,deltax=L/nx;   
    double deltaxh=3.*deltax;
	deltaxh=deltaxh*(1+1e-9);

	//If,Ib
	double** If=new double*[nx+1];
	double** Ib=new double*[nx+1];
	double** Iff=new double*[nx+1];
	double** Ibb=new double*[nx+1];

	double** para1f=new double*[nx+1];double** para1b=new double*[nx+1];
	double** para_d1f=new double*[nx+1];double** para_d1b=new double*[nx+1];
	int** HN1f=new int*[nx+1];int** HN1b=new int*[nx+1];

	for(int l=0;l<nx+1;l++){
		If[l]=new double[glnu];
		Ib[l]=new double[glnu];
		Iff[l]=new double[glnu];
		Ibb[l]=new double[glnu];

		para1f[l]=new double[glnu];
		para1b[l]=new double[glnu];
		para_d1f[l]=new double[glnu];
		para_d1b[l]=new double[glnu];
		HN1f[l]=new int[glnu];
		HN1b[l]=new int[glnu];
	}
	
	for(int l=0;l<nx+1;l++)
       { x[l]=l*deltax;
	   	 xlong[l]=deltax;
		 beta=double(l)/double(nx);
         T[l]=(1-beta)*Tl+beta*Tr;  //guess the temperature
         }     
		


	double aa=0,bb=0,maxx=1.,minx=0.;
	double deltaGL=(maxx-minx)/(double)nu;
	for(int i=0;i<nu;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
		u2[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u2[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u2[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u2[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
    }

	aa=0;bb=0;maxx=-1.;minx=0.;
	deltaGL=(maxx-minx)/(double)nu;
	for(int i=0;i<nu;i++){
        aa=(double)i*deltaGL;bb=((double)(i+1))*deltaGL;
		u1[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u1[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u1[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u1[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
    }

	for (int l=0;l<nx+1;l++){	
		for(int n=0;n<glnu;n++){
				If[l][n]=C*v*T[l]/4./pi;
				Ib[l][n]=C*v*T[l]/4./pi;
		}        
	}

	for(int l=0;l<nx+1;l++){	
			sum1[l]=GL4(If[l],0.,1.,nu);
		    sum2[l]=GL4(Ib[l],-1.,0.,nu);
		    I0[l]=(sum1[l]+sum2[l])/2.; 
			I0_ll[l]=I0[l]/ll; 
	} 

	double qa1=C*v*Tl/4.;
	double qb2=C*v*Tr/4.;

	makeparaf(u2, x,xlong,deltaxh, para1f,para_d1f, HN1f);
	makeparab(u1, x,xlong,deltaxh, para1b,para_d1b, HN1b);

	for(int l=0;l<nx+1;l++){
		sum1[l]=GL4(If[l],0.,1.,nu);
		sum2[l]=GL4(Ib[l],-1.,0.,nu);

		I0[l]=(sum1[l]+sum2[l])/2.;   
		I0_ll[l]=I0[l]/ll;

	}

    int j=0;double residual=10.;

    while(residual>1e-6&&j<200000){
		double c=0.;
		for(int l=0;l<nx+1;l++){
			I00[l]=I0[l];
		}       

		for (int l=0;l<nx+1;l++){
			for(int n=0;n<glnu;n++){
				If[l][n]=eprt_pddo_f(l,n,u2[n],I0_ll,deltax,L,If,x,xlong,deltaxh,para1f[l][n],para_d1f[l][n],HN1f[l][n]);
			}  
		}

        for (int l=nx;l>=0;l--){
			for(int n=0;n<glnu;n++){
				Ib[l][n]=eprt_pddo_b(l,n,u1[n],I0_ll,deltax,L,Ib,x,xlong,deltaxh,para1b[l][n],para_d1b[l][n],HN1b[l][n]);
			}  
		}

		for(int l=0;l<nx+1;l++){
			sum1[l]=GL4(If[l],0.,1.,nu);
		    sum2[l]=GL4(Ib[l],-1.,0.,nu);

		    I0[l]=(sum1[l]+sum2[l])/2.;   
			I0_ll[l]=I0[l]/ll;
			c=max(fabs(I00[l]-I0[l]),c);
	    }

		residual=c;
		j++;                             
	} 

	cout<<j<<endl;  
    cout<<residual<<endl;  


    for(int l=0;l<nx+1;l++){
		sum1[l]=GL4(If[l],0.,1.,nu);
		sum2[l]=GL4(Ib[l],-1.,0.,nu);

		I0[l]=(sum1[l]+sum2[l])/2.;

		T[l]=4.*pi*I0[l]/C/v;    
		T2[l]=T[l]-Tr;      

		sum1[l]=GL4(If[l],0.,1.,nu,u2);
		sum2[l]=GL4(Ib[l],-1.,0.,nu,u1); 

		q[l]=2*pi*(sum2[l]+sum1[l]);

	    k[l]=q[l]*(L);
		k2[l]=q[l]/(qa1-qb2);
	}
	cout<<j<<endl;

	string file0="PD_1D_kn_"+to_string(kn)+".xml";
	out_file(file0,T2,x,nx,q,k2);

    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;

    cout<<"\n time is "<<showpoint<<setprecision(8)<<totaltime<<" s"<<endl;
    return 0;
}