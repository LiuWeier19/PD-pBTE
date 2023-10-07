//  The following code is the simulation code for "Two-dimensional phonon heat conduction along thin films" as presented in the Manuscript. 
//  To meet the approximation conditions for the analytical solution, deltaT needs to be significantly smaller than T0. Additionally, for small Knudsen numbers kn, 
//  using a deltaT that is too small can lead to significant computational errors due to the substantial difference between the deviation function g and e^0. 
//  Therefore, deltaT should be adjusted for different Knudsen numbers. When the computed results approach periodicity, we use the results at x/nx=0.5 as the outcomes presented in the Manuscript.

#include"../calFunction.h"
using namespace std;
const double kb = 1.3806488e-23;// Planck constant, J s.
const double h = 6.62606957e-34;// Dirac constant, J s.
const double hbar = 1.054571726e-34;
const double pi =3.1415926535897932384626433832795028841953;
const double C=0.930e6,v=1804.;//v-velocity,C-heat capacity
const double T0=300.;// T0- Temperature at low temperature boundary
const double deltaT=0.00005;// Temperature difference between two boundaries
const int nx=10,n1=8,n2=16,nz=1000;//nx-x ,n1-theta ,n2-phi,nz-z
const double llsi=260.4e-9;//mean free path of Si
const double para3=1/llsi;
const double kn=10.;
const double D=llsi/kn,deltaz=D/nz;
const double delta1=pi/2./n1,delta2=pi/n2,deltax=deltaz;
const double L=nx*deltax;
const int nr=3;
const double areaSpace=deltax*deltax;
const double areaSpace2=deltax;

const double kbu=C*v*llsi/3.;
const int gln=4;
const int glnn1=gln*n1;
const int glnn2=gln*n2;

struct realHlist{
  int nodeNumX;
  int nodeNumZ;

  double XiX;
  double XiZ;

  realHlist* next;
};

void out_file(string file,double T[][nz+1],double *x,int nx,int nz,double *q,double *k){
    ofstream outfile;
    outfile.open(file);
    outfile<<"x,"<<"qbar,"<<"k*,"<<"q"<<endl;
	for(int i=0;i<nx+2*nr;i++){
		outfile<<showpoint<<setprecision(8)<<x[i]<<","<<showpoint<<setprecision(8)<<q[i]<<","<<showpoint<<setprecision(8)<<k[i]<<",";
		for(int j=0;j<nz;j++){
			outfile<<showpoint<<setprecision(10)<<T[i][j]<<",";
		}
     outfile<<showpoint<<setprecision(10)<<T[i][nz]<<endl;
	}
	outfile.close();
}
double simpson_qx_z(double *qx,int nz,double deltaz){
  return simpson(qx,nz,deltaz);
}
double simpson_deltaT(double *T1,double *T0,int nz,double deltaz){
  return simpson(T1,nz,deltaz)-simpson(T0,nz,deltaz);
}


double PD_eprtff(double para1,double para2,int i,int j,int m,int n,double I0[nx+2*nr][nz+1],double deltax,double deltaz,double L,double D,double p,double qf[nx+2*nr][nz+1],double qb[nx+2*nr][nz+1],double qu[nx+2*nr][nz+1],double qd[nx+2*nr][nz+1],double **tempff[nx+2*nr][nz+1],double deltaxh,double paraS2ff,double paraC2ff,double paraSCff, double para1ff, double para21ff, double para22ff,int HNnode,int *HrnX,int *HrnZ,double *HrnXiX,double *HrnXiZ){
    double temp;
    int nrr;
    
    if(i<nr){
      temp=tempff[i][j][m][n];
    }
    else if(i>nx+nr-1){
      temp=tempff[i][j][m][n];
    }
    else if(j==0){
      temp=(-qd[i][j]/pi);
    }

    else{
        double sumx=0;
        double sumz=0;
        double dI1=0;
        
        for(int l=0;l<HNnode;l++){
          if((i-HrnX[l])*para1+(j-HrnZ[l])*para2>0.){
            sumx+=tempff[HrnX[l]][HrnZ[l]][m][n]*HrnXiX[l];
            sumz+=tempff[HrnX[l]][HrnZ[l]][m][n]*HrnXiZ[l];
          }
        }
        dI1=(sumx*para21ff)+(sumz*para22ff)-(I0[i][j]*para3);
        temp=dI1*para1ff;

      }

	return temp;
}


double PD_eprtbf(double para1,double para2,int i,int j,int m,int n,double I0[nx+2*nr][nz+1],double deltax,double deltaz,double L,double D,double p,double qf[nx+2*nr][nz+1],double qb[nx+2*nr][nz+1],double qu[nx+2*nr][nz+1],double qd[nx+2*nr][nz+1],double **tempff[nx+2*nr][nz+1],double deltaxh,double paraS2bf,double paraC2bf,double paraSCbf, double para1bf, double para21bf, double para22bf,int HNnode,int *HrnX,int *HrnZ,double *HrnXiX,double *HrnXiZ){
    double temp;
    int nrr;
    if(i>nx+nr-1){
      temp=tempff[i][j][m][n];
    }
    else if(i<nr){
      temp=tempff[i][j][m][n];
    }
    else if(j==0){
        temp=(-qd[i][j]/pi);
    }
    
    else{
        double sumx=0;
        double sumz=0;
        double dI1=0;
        
        for(int l=0;l<HNnode;l++){
          if((i-HrnX[l])*para1+(j-HrnZ[l])*para2>0.){
            sumx+=tempff[HrnX[l]][HrnZ[l]][m][n]*HrnXiX[l];
            sumz+=tempff[HrnX[l]][HrnZ[l]][m][n]*HrnXiZ[l];
          }
        }
        dI1=(sumx*para21bf)+(sumz*para22bf)-(I0[i][j]*para3);
        temp=dI1*para1bf;
      }
  
	return temp;
}

void makeparaff(double **paraS2ff[nx+2*nr][nz+1],double **paraC2ff[nx+2*nr][nz+1],double **paraSCff[nx+2*nr][nz+1],double *sintheta1,double *sintheta2,double *costheta1,double *costheta2,double *sinu1,double *sinu2, double **para1ff[nx+2*nr][nz+1],double **para21ff[nx+2*nr][nz+1],double **para22ff[nx+2*nr][nz+1]){
  double area=areaSpace;
  int nrr;

  double sum1=0;
  double sum2=0;
  double xi=0;
  for(int l=0;l<nx+2*nr;l++){
  	for(int m=0;m<nz+1;m++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          paraS2ff[l][m][n][i]=0;
          paraC2ff[l][m][n][i]=0;
          paraSCff[l][m][n][i]=0;
          para1ff[l][m][n][i]=0;
          sum1=0;
          sum2=0;
            for(int j=0;j<2*nr+1;j++){
              for(int k=0;k<2*nr+1;k++){
                if(j==nr&&k==nr){
                }else if(abs(j-nr)*abs(j-nr)+abs(k-nr)*abs(k-nr)<=nr*nr){
                  if((nr-j)*costheta1[i]+(nr-k)*sintheta1[i]*sinu1[n]>0.){
                    if(l-nr+j>=0&&l-nr+j<=nx+2*nr-1&&m-nr+k>=0&&m-nr+k<=nz){
                      paraS2ff[l][m][n][i]+=area*(k-nr)*deltaz*(k-nr)*deltaz/(deltax*deltax*(j-nr)*(j-nr)+deltaz*deltaz*(k-nr)*(k-nr));
                      paraC2ff[l][m][n][i]+=area*(j-nr)*deltax*(j-nr)*deltax/(deltax*deltax*(j-nr)*(j-nr)+deltaz*deltaz*(k-nr)*(k-nr));
                      paraSCff[l][m][n][i]+=area*(j-nr)*deltax*(k-nr)*deltaz/(deltax*deltax*(j-nr)*(j-nr)+deltaz*deltaz*(k-nr)*(k-nr));

                      xi=areaSpace2/((j-nr)*(j-nr)+(k-nr)*(k-nr));
                      sum1+=(xi*(j-nr));
                      sum2+=(xi*(k-nr));
                    }
                    
                  }     
                }
              }
            }

          if(paraSCff[l][m][n][i]==0){
            if(paraC2ff[l][m][n][i]==0){
              if(paraS2ff[l][m][n][i]!=0){
                para21ff[l][m][n][i]=0.;
                para22ff[l][m][n][i]=sintheta1[i]*sinu1[n]/paraS2ff[l][m][n][i];
                para1ff[l][m][n][i]=sum2*para22ff[l][m][n][i];
              }
              
            }else if(paraS2ff[l][m][n][i]==0){
              para21ff[l][m][n][i]=costheta1[i]/paraC2ff[l][m][n][i];
              para22ff[l][m][n][i]=0.;
              para1ff[l][m][n][i]=sum1*para21ff[l][m][n][i];
            }else{
              para21ff[l][m][n][i]=costheta1[i]/paraC2ff[l][m][n][i];
              para22ff[l][m][n][i]=sintheta1[i]*sinu1[n]/paraS2ff[l][m][n][i];
              para1ff[l][m][n][i]=sum1*para21ff[l][m][n][i]+sum2*para22ff[l][m][n][i];
            }
          }else{
            para21ff[l][m][n][i]=((costheta1[i]*(paraS2ff[l][m][n][i]/paraSCff[l][m][n][i])/(paraC2ff[l][m][n][i]*(paraS2ff[l][m][n][i]/paraSCff[l][m][n][i])-paraSCff[l][m][n][i]))-(sintheta1[i]*sinu1[n]/(paraS2ff[l][m][n][i]*(paraC2ff[l][m][n][i]/paraSCff[l][m][n][i])-paraSCff[l][m][n][i])));
            para22ff[l][m][n][i]=((sintheta1[i]*sinu1[n]*(paraC2ff[l][m][n][i]/paraSCff[l][m][n][i])/(paraS2ff[l][m][n][i]*(paraC2ff[l][m][n][i]/paraSCff[l][m][n][i])-paraSCff[l][m][n][i]))-(costheta1[i]/(paraC2ff[l][m][n][i]*(paraS2ff[l][m][n][i]/paraSCff[l][m][n][i])-paraSCff[l][m][n][i])));
            para1ff[l][m][n][i]=sum1*para21ff[l][m][n][i]+sum2*para22ff[l][m][n][i];
          }
          para1ff[l][m][n][i]+=-1./llsi;
          para1ff[l][m][n][i]=1./para1ff[l][m][n][i];
        }
      }
  	}
  } 
}

void makeparabf(double **paraS2bf[nx+2*nr][nz+1],double **paraC2bf[nx+2*nr][nz+1],double **paraSCbf[nx+2*nr][nz+1],double *sintheta1,double *sintheta2,double *costheta1,double *costheta2,double *sinu1,double *sinu2, double **para1bf[nx+2*nr][nz+1], double **para21bf[nx+2*nr][nz+1], double **para22bf[nx+2*nr][nz+1]){
  double area=areaSpace;
  int nrr;

  double sum1=0;
  double sum2=0;
  double xi=0;
  for(int l=0;l<nx+2*nr;l++){
  	for(int m=0;m<nz+1;m++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          paraS2bf[l][m][n][i]=0;
          paraC2bf[l][m][n][i]=0;
          paraSCbf[l][m][n][i]=0;
          para1bf[l][m][n][i]=0;
          sum1=0;
          sum2=0;

            for(int j=0;j<2*nr+1;j++){
              for(int k=0;k<2*nr+1;k++){
                if(j==nr&&k==nr){
                }else if(abs(j-nr)*abs(j-nr)+abs(k-nr)*abs(k-nr)<=nr*nr){
                  if((nr-j)*costheta2[i]+(nr-k)*sintheta2[i]*sinu1[n]>0.){
                    if(l-nr+j>=0&&l-nr+j<=nx+2*nr-1&&m-nr+k>=0&&m-nr+k<=nz){
                      paraS2bf[l][m][n][i]+=area*(k-nr)*deltax*(k-nr)*deltax/(deltax*deltax*(j-nr)*(j-nr)+deltax*deltax*(k-nr)*(k-nr));
                      paraC2bf[l][m][n][i]+=area*(j-nr)*deltax*(j-nr)*deltax/(deltax*deltax*(j-nr)*(j-nr)+deltax*deltax*(k-nr)*(k-nr));
                      paraSCbf[l][m][n][i]+=area*(j-nr)*deltax*(k-nr)*deltax/(deltax*deltax*(j-nr)*(j-nr)+deltax*deltax*(k-nr)*(k-nr));

                      xi=areaSpace2/((j-nr)*(j-nr)+(k-nr)*(k-nr));
                      sum1+=(xi*(j-nr));
                      sum2+=(xi*(k-nr));
                    }
                    
                  }     
                }
              }
            }

          if(paraSCbf[l][m][n][i]==0){
            if(paraC2bf[l][m][n][i]==0){
              if(paraS2bf[l][m][n][i]!=0){
                para21bf[l][m][n][i]=0.;
                para22bf[l][m][n][i]=sintheta2[i]*sinu1[n]/paraS2bf[l][m][n][i];
                para1bf[l][m][n][i]=para22bf[l][m][n][i]*sum2;
              } 
            }else if(paraS2bf[l][m][n][i]==0){
              para21bf[l][m][n][i]=costheta2[i]/paraC2bf[l][m][n][i];
              para22bf[l][m][n][i]=0.;
              para1bf[l][m][n][i]=para21bf[l][m][n][i]*sum1;
            }else{
              para21bf[l][m][n][i]=costheta2[i]/paraC2bf[l][m][n][i];
              para22bf[l][m][n][i]=sintheta2[i]*sinu1[n]/paraS2bf[l][m][n][i];
              para1bf[l][m][n][i]=sum1*para21bf[l][m][n][i]+sum2*para22bf[l][m][n][i];
            }
          }else{
            para21bf[l][m][n][i]=((costheta2[i]*(paraS2bf[l][m][n][i]/paraSCbf[l][m][n][i])/(paraC2bf[l][m][n][i]*(paraS2bf[l][m][n][i]/paraSCbf[l][m][n][i])-paraSCbf[l][m][n][i]))-(sintheta2[i]*sinu1[n]/(paraS2bf[l][m][n][i]*(paraC2bf[l][m][n][i]/paraSCbf[l][m][n][i])-paraSCbf[l][m][n][i])));
            para22bf[l][m][n][i]=((sintheta2[i]*sinu1[n]*(paraC2bf[l][m][n][i]/paraSCbf[l][m][n][i])/(paraS2bf[l][m][n][i]*(paraC2bf[l][m][n][i]/paraSCbf[l][m][n][i])-paraSCbf[l][m][n][i]))-(costheta2[i]/(paraC2bf[l][m][n][i]*(paraS2bf[l][m][n][i]/paraSCbf[l][m][n][i])-paraSCbf[l][m][n][i])));
            para1bf[l][m][n][i]=sum1*para21bf[l][m][n][i]+sum2*para22bf[l][m][n][i];
          }
          para1bf[l][m][n][i]+=-1./llsi;

          para1bf[l][m][n][i]=1./para1bf[l][m][n][i];
        }
      }
  	}
  } 
}

void initializeRealNodes(int HNnode[nx+2*nr][nz+1],realHlist *Hrn[nx+2*nr][nz+1],int *HrnX[nx+2*nr][nz+1],int *HrnZ[nx+2*nr][nz+1],double *HrnXiX[nx+2*nr][nz+1],double *HrnXiZ[nx+2*nr][nz+1]){
  realHlist *lastnode;
  realHlist *thisnode;
  for(int l=0;l<nx+2*nr;l++){
  	for(int m=0;m<nz+1;m++){
      HNnode[l][m]=0;

      for(int j=0;j<2*nr+1;j++){
        for(int k=0;k<2*nr+1;k++){
          if(j==nr&&k==nr){
          }else if(abs(j-nr)*abs(j-nr)+abs(k-nr)*abs(k-nr)<=nr*nr){
            if(l-nr+j>=0&&l-nr+j<=nx+2*nr-1&&m-nr+k>=0&&m-nr+k<=nz){
              if(HNnode[l][m]==0){
                Hrn[l][m]=new realHlist;
                Hrn[l][m]->nodeNumX=l-nr+j;
                Hrn[l][m]->nodeNumZ=m-nr+k;

                Hrn[l][m]->XiX=(j-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                Hrn[l][m]->XiZ=(k-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                lastnode=Hrn[l][m];
              }else{
                thisnode=new realHlist;
                thisnode->nodeNumX=l-nr+j;
                thisnode->nodeNumZ=m-nr+k;

                thisnode->XiX=(j-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                thisnode->XiZ=(k-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                lastnode->next=thisnode;
                lastnode=thisnode;
              }
               HNnode[l][m]++;
               
            }     
          }
        }
      }

      HrnX[l][m]=new int[HNnode[l][m]];
      HrnZ[l][m]=new int[HNnode[l][m]];
      HrnXiX[l][m]=new double[HNnode[l][m]];
      HrnXiZ[l][m]=new double[HNnode[l][m]];

      thisnode=Hrn[l][m];
      for(int ii=0;ii<HNnode[l][m];ii++){
        HrnX[l][m][ii]=thisnode->nodeNumX;HrnZ[l][m][ii]=thisnode->nodeNumZ;HrnXiX[l][m][ii]=thisnode->XiX;HrnXiZ[l][m][ii]=thisnode->XiZ;

        lastnode=thisnode;
        thisnode=lastnode->next;

        delete lastnode;
      }
  	}
  }
}


int main(){
    clock_t start,finish;
    cout<<kn<<endl;
    double beta,residual=10.;
    double p=0.0; 
    double theta1[glnn1],theta2[glnn1],u1[glnn2],u2[glnn2],Iff[nx+2*nr][nz+1],Ifb[nx+2*nr][nz+1],Ibf[nx+2*nr][nz+1],Ibb[nx+2*nr][nz+1];

    double qf[nx+2*nr][nz+1],qb[nx+2*nr][nz+1],qu[nx+2*nr][nz+1],qd[nx+2*nr][nz+1];
    double T[nx+2*nr][nz+1],k[nx+2*nr],x[nx+2*nr],y[nz+1],qxx[nx+2*nr];
    double q[nx+2*nr][nz+1],TT,I0[nx+2*nr][nz+1],I00[nx+2*nr][nz+1],hf[nx+2*nr][nz+1];
    double totaltime;
    double Ifx[nx+2*nr][nz+1],Ibx[nx+2*nr][nz+1];

    int HNnode[nx+2*nr][nz+1];
    realHlist *Hrn[nx+2*nr][nz+1];

    int *HrnX[nx+2*nr][nz+1], *HrnZ[nx+2*nr][nz+1];
    double *HrnXiX[nx+2*nr][nz+1],*HrnXiZ[nx+2*nr][nz+1];


    double sintheta1[glnn1],sintheta2[glnn1],sinu1[glnn2],sinu2[glnn2];
    double costheta1[glnn1],costheta2[glnn1],cosu1[glnn2],cosu2[glnn2];
    double para11[glnn1],para21[glnn1][glnn2];
    double para12[glnn1],para22[glnn1][glnn2];
    double sincostheta1[glnn1],sincostheta2[glnn1];
    double sinsintheta1[glnn1],sinsintheta2[glnn1];

    double **tempff[nx+2*nr][nz+1],**tempbf[nx+2*nr][nz+1];
    double **ltempff[nx+2*nr][nz+1],**ltempbf[nx+2*nr][nz+1];
    double **paraS2ff[nx+2*nr][nz+1],**paraS2bf[nx+2*nr][nz+1];
    double **paraC2ff[nx+2*nr][nz+1],**paraC2bf[nx+2*nr][nz+1];
    double **paraSCff[nx+2*nr][nz+1],**paraSCbf[nx+2*nr][nz+1];

    double **para1ff[nx+2*nr][nz+1],**para1bf[nx+2*nr][nz+1];

    double **para21ff[nx+2*nr][nz+1],**para21bf[nx+2*nr][nz+1];
    double **para22ff[nx+2*nr][nz+1],**para22bf[nx+2*nr][nz+1];

    start=clock();  //timing

    for(int i=0;i<nx+2*nr;i++){
      for(int j=0;j<nz+1;j++){
        tempff[i][j]=new double*[glnn2];
        tempbf[i][j]=new double*[glnn2];
        ltempff[i][j]=new double*[glnn2];
        ltempbf[i][j]=new double*[glnn2];
        paraS2ff[i][j]=new double*[glnn2];
        paraS2bf[i][j]=new double*[glnn2];
        paraC2ff[i][j]=new double*[glnn2];
        paraC2bf[i][j]=new double*[glnn2];
        paraSCff[i][j]=new double*[glnn2];
        paraSCbf[i][j]=new double*[glnn2];

        para1ff[i][j]=new double*[glnn2];
        para1bf[i][j]=new double*[glnn2];

        para21ff[i][j]=new double*[glnn2];
        para21bf[i][j]=new double*[glnn2];
        para22ff[i][j]=new double*[glnn2];
        para22bf[i][j]=new double*[glnn2];
        for(int k=0;k<glnn2;k++){
          tempff[i][j][k]=new double[glnn1];
          tempbf[i][j][k]=new double[glnn1];
          ltempff[i][j][k]=new double[glnn1];
          ltempbf[i][j][k]=new double[glnn1];
          paraS2ff[i][j][k]=new double[glnn1];
          paraS2bf[i][j][k]=new double[glnn1];
          paraC2ff[i][j][k]=new double[glnn1];
          paraC2bf[i][j][k]=new double[glnn1];
          paraSCff[i][j][k]=new double[glnn1];
          paraSCbf[i][j][k]=new double[glnn1];

          para1ff[i][j][k]=new double[glnn1];
          para1bf[i][j][k]=new double[glnn1];

          para21ff[i][j][k]=new double[glnn1];
          para21bf[i][j][k]=new double[glnn1];
          para22ff[i][j][k]=new double[glnn1];
          para22bf[i][j][k]=new double[glnn1];
        }
      }
    }


    double *Isffd[glnn2],*Isfbd[glnn2],*Isbfd[glnn2],*Isbbd[glnn2];

    for(int i=0;i<glnn2;i++){
      Isffd[i]=new double[glnn1];
      Isfbd[i]=new double[glnn1];
      Isbfd[i]=new double[glnn1];
      Isbbd[i]=new double[glnn1];
    }

    double deltaxh=deltax*nr;

    double ddeltaT=deltaT/nx;

    for(int l=0;l<nx+2*nr;l++){ 
		  x[l]=(l-nr+1)*deltax;
      for(int i=0;i<nz+1;i++){
        T[l][i]=T0+ddeltaT*(l-nr+1)-ddeltaT*0.5;
      }
           //guess the temputure
    }  


    double aa=0,bb=0,maxx=pi/2.,minx=0.;
	  double deltaGL=(maxx-minx)/(double)n1;
	  for(int i=0;i<n1;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;
		  theta1[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;theta1[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;theta1[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;theta1[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sintheta1[4*i]=sin(theta1[4*i]);sintheta1[4*i+1]=sin(theta1[4*i+1]);sintheta1[4*i+2]=sin(theta1[4*i+2]);sintheta1[4*i+3]=sin(theta1[4*i+3]);
      costheta1[4*i]=cos(theta1[4*i]);costheta1[4*i+1]=cos(theta1[4*i+1]);costheta1[4*i+2]=cos(theta1[4*i+2]);costheta1[4*i+3]=cos(theta1[4*i+3]);
      sincostheta1[4*i]=sintheta1[4*i]*costheta1[4*i];sincostheta1[4*i+1]=sintheta1[4*i+1]*costheta1[4*i+1];sincostheta1[4*i+2]=sintheta1[4*i+2]*costheta1[4*i+2];sincostheta1[4*i+3]=sintheta1[4*i+3]*costheta1[4*i+3];
      sinsintheta1[4*i]=sintheta1[4*i]*sintheta1[4*i];sinsintheta1[4*i+1]=sintheta1[4*i+1]*sintheta1[4*i+1];sinsintheta1[4*i+2]=sintheta1[4*i+2]*sintheta1[4*i+2];sinsintheta1[4*i+3]=sintheta1[4*i+3]*sintheta1[4*i+3];
      para11[4*i]=costheta1[4*i]/deltax;para11[4*i+1]=costheta1[4*i+1]/deltax;para11[4*i+2]=costheta1[4*i+2]/deltax;para11[4*i+3]=costheta1[4*i+3]/deltax;
    }

	  aa=0;bb=0;maxx=pi/2.;minx=pi;
	  deltaGL=(maxx-minx)/(double)n1;
	  for(int i=0;i<n1;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;
      theta2[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;theta2[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;theta2[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;theta2[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sintheta2[4*i]=sin(theta2[4*i]);sintheta2[4*i+1]=sin(theta2[4*i+1]);sintheta2[4*i+2]=sin(theta2[4*i+2]);sintheta2[4*i+3]=sin(theta2[4*i+3]);
      costheta2[4*i]=cos(theta2[4*i]);costheta2[4*i+1]=cos(theta2[4*i+1]);costheta2[4*i+2]=cos(theta2[4*i+2]);costheta2[4*i+3]=cos(theta2[4*i+3]);
      sincostheta2[4*i]=sintheta2[4*i]*costheta2[4*i];sincostheta2[4*i+1]=sintheta2[4*i+1]*costheta2[4*i+1];sincostheta2[4*i+2]=sintheta2[4*i+2]*costheta2[4*i+2];sincostheta2[4*i+3]=sintheta2[4*i+3]*costheta2[4*i+3];
      sinsintheta2[4*i]=sintheta2[4*i]*sintheta2[4*i];sinsintheta2[4*i+1]=sintheta2[4*i+1]*sintheta2[4*i+1];sinsintheta2[4*i+2]=sintheta2[4*i+2]*sintheta2[4*i+2];sinsintheta2[4*i+3]=sintheta2[4*i+3]*sintheta2[4*i+3];
      para12[4*i]=costheta2[4*i]/deltax;para12[4*i+1]=costheta2[4*i+1]/deltax;para12[4*i+2]=costheta2[4*i+2]/deltax;para12[4*i+3]=costheta2[4*i+3]/deltax;
    }

    aa=0;bb=0;maxx=pi;minx=0.;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;
		  u1[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u1[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u1[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u1[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu1[4*i]=sin(u1[4*i]);sinu1[4*i+1]=sin(u1[4*i+1]);sinu1[4*i+2]=sin(u1[4*i+2]);sinu1[4*i+3]=sin(u1[4*i+3]);
      cosu1[4*i]=cos(u1[4*i]);cosu1[4*i+1]=cos(u1[4*i+1]);cosu1[4*i+2]=cos(u1[4*i+2]);cosu1[4*i+3]=cos(u1[4*i+3]);
    }

    aa=0;bb=0;maxx=pi;minx=2.*pi;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;
		  u2[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u2[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u2[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u2[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu2[4*i]=sin(u2[4*i]);sinu2[4*i+1]=sin(u2[4*i+1]);sinu2[4*i+2]=sin(u2[4*i+2]);sinu2[4*i+3]=sin(u2[4*i+3]);
      cosu2[4*i]=cos(u2[4*i]);cosu2[4*i+1]=cos(u2[4*i+1]);cosu2[4*i+2]=cos(u2[4*i+2]);cosu2[4*i+3]=cos(u2[4*i+3]);
    }

    for (int n=0;n<glnn1;n++){
        for (int m=0;m<glnn2;m++){
            para21[n][m]=sintheta1[n]*sinu1[m]/deltaz;
            para22[n][m]=sintheta2[n]*sinu1[m]/deltaz;
        }
    }

    initializeRealNodes(HNnode,Hrn,HrnX,HrnZ,HrnXiX,HrnXiZ);

    makeparaff(paraS2ff,paraC2ff,paraSCff,sintheta1,sintheta2,costheta1,costheta2,sinu1,sinu2,para1ff,para21ff,para22ff);
    makeparabf(paraS2bf,paraC2bf,paraSCbf,sintheta1,sintheta2,costheta1,costheta2,sinu1,sinu2,para1bf,para21bf,para22bf);

        for(int l=0;l<nx+2*nr;l++){
            for(int m=0;m<nz+1;m++){
                I0[l][m]=C*v*T[l][m]/4./pi;
                qd[l][m]=-I0[l][m]*pi;
            }
        }

        finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\n time is "<<totaltime<<" s"<<endl;

        int j=0;
        residual=10.;
        while(residual>1e-4&&j<50000){
            double a=0;

            memcpy(I00,I0,sizeof(I0));//I00

            if(j==0){
              for(int l=0;l<nx+2*nr;l++){
  		  	      for(int m=0;m<nz+1;m++){

                  for(int n=0;n<glnn2;n++){
          for(int i=0;i<glnn1;i++){
            if(m==0){
              tempff[l][m][n][i]=I0[l][m];
              tempbf[l][m][n][i]=I0[l][m];
            }
            else if(i==0||n==0||n==n2){
              tempff[l][m][n][i]=I0[l][m]-llsi*costheta1[i]*C*v*deltaT*(1.)/4./pi/L;
              tempbf[l][m][n][i]=I0[l][m]-llsi*costheta2[i]*C*v*deltaT*(1.)/4./pi/L;
            }else{
              tempff[l][m][n][i]=I0[l][m]-llsi*costheta1[i]*C*v*deltaT*(1.-exp(-m*deltaz/llsi/(sintheta1[i]*sinu1[n])))/4./pi/L;
              tempbf[l][m][n][i]=I0[l][m]-llsi*costheta2[i]*C*v*deltaT*(1.-exp(-m*deltaz/llsi/(sintheta2[i]*sinu1[n])))/4./pi/L;
            }
            
          }
        }
  			        }
  		        } 
            }else{
              for(int l=0;l<nx+2*nr;l++){
  		  	      for(int m=0;m<nz+1;m++){
                  if(l<nr){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                            tempff[l][m][n][i]=tempff[l+nx][m][n][i]-I0[l+nx][m]+I0[l][m];
                            tempbf[l][m][n][i]=tempbf[l+nx][m][n][i]-I0[l+nx][m]+I0[l][m];
                        }
                    }
                  }

                  if(l>nx+nr-1){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                          tempff[l][m][n][i]=tempff[l-nx][m][n][i]-I0[l-nx][m]+I0[l][m];
                          tempbf[l][m][n][i]=tempbf[l-nx][m][n][i]-I0[l-nx][m]+I0[l][m];
                        }
                    }
                  
                  }
  			        }
  		        } 
            }

            for(int l=0;l<nx+2*nr;l++){
  		  	    for(int m=0;m<nz+1;m++){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                            tempff[l][m][n][i]=PD_eprtff(costheta1[i],sintheta1[i]*sinu1[n],l,m,n,i,I0,deltax,deltaz,L,D,p,qf,qb,qu,qd,tempff,deltaxh,paraS2ff[l][m][n][i],paraC2ff[l][m][n][i],paraSCff[l][m][n][i],para1ff[l][m][n][i],para21ff[l][m][n][i],para22ff[l][m][n][i],HNnode[l][m],HrnX[l][m],HrnZ[l][m],HrnXiX[l][m],HrnXiZ[l][m]);
                        }
                    }
  			    }
  		    }  


            for(int l=nx+2*nr-1;l>=0;l--){
  		  	    for(int m=0;m<nz+1;m++){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                            tempbf[l][m][n][i]=PD_eprtbf(costheta2[i],sintheta2[i]*sinu1[n],l,m,n,i,I0,deltax,deltaz,L,D,p,qf,qb,qu,qd,tempbf,deltaxh,paraS2bf[l][m][n][i],paraC2bf[l][m][n][i],paraSCbf[l][m][n][i],para1bf[l][m][n][i],para21bf[l][m][n][i],para22bf[l][m][n][i],HNnode[l][m],HrnX[l][m],HrnZ[l][m],HrnXiX[l][m],HrnXiZ[l][m]);
                        }
                    }
  			    }
  		    }

            for(int l=0;l<nx+2*nr;l++){
  		  	    for(int m=0;m<nz+1;m++){
                    if(l>=nr&&l<nx+nr){
                        I0[l][m]=(GL4_I0(tempff[l][m],sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(tempff[l][nz-m],sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(tempbf[l][m],sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(tempbf[l][nz-m],sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;

                      if(m==0){
                        qd[l][m]=(GL4_q_2(tempff[l][nz-m],sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(tempbf[l][nz-m],sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2));
                      }
                    }

                    a=max(fabs(I0[l][m]-I00[l][m]),a);  
                }
            }

            residual=a;
            j++;
            cout<<j<<endl; 
            cout<<residual<<endl;
        }

    for(int l=0;l<nx+2*nr;l++){
  		  	    for(int m=0;m<nz+1;m++){
                      qf[l][m]=(GL4_I0(tempff[l][m],sincostheta1,0.,pi/2.,n1,0.,pi,n2)+GL4_I0(tempff[l][nz-m],sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2));
                      qb[l][m]=(GL4_I0(tempbf[l][m],sincostheta2,pi/2.,pi,n1,0.,pi,n2)+GL4_I0(tempbf[l][nz-m],sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2));
                      q[l][m]= qf[l][m]+qb[l][m]; 

                    if(l>=nr&&l<nx+nr){
                      I0[l][m]=(GL4_I0(tempff[l][m],sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(tempff[l][nz-m],sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(tempbf[l][m],sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(tempbf[l][nz-m],sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;
                    }

                    if(m==0){
                      qd[l][m]=(GL4_q_2(tempff[l][nz-m],sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(tempbf[l][nz-m],sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2));
                    }  

                    T[l][m]=4.*pi*I0[l][m]/C/v;
                }
            }

            for(int l=0;l<nx+2*nr;l++){
                TT=simpson_deltaT(T[0],T[nx],nz,deltaz); 
                qxx[l]=simpson_qx_z(q[l],nz,deltaz);  

                k[l]=qxx[l]*L/(deltaT*D)/kbu;
            }

        string file0="PD_kn_"+to_string(kn)+"_Q_r_"+to_string(nr)+"_deltaT_"+to_string(deltaT)+".xml";
        out_file(file0,q,x,nx, nz,qxx,k);

    //}
    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\n time is "<<totaltime<<" s"<<endl; 
}
