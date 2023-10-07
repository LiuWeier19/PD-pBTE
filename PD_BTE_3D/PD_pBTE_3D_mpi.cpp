//  The following code is the simulation code for "Three-dimensional steady-state phonon transport" as presented in the Manuscript.
//  Please note that this code does not provide a widely applicable angle allocation function. If you need to reproduce the examples 
//  in the manuscript, the angle allocation part in this code should suffice. If you encounter any errors, please write a separate angle allocation function for each process.
#include "mpi.h"
#include"../calFunction.h"

const double kb = 1.3806488e-23;// Planck constant, J s.
const double h = 6.62606957e-34;// Dirac constant, J s.
const double hbar = 1.054571726e-34;
const double pi =3.1415926535897932384626433832795028841953;
const double C=0.930e6,v=1804.;
const double T0=300.;
const double deltaT=1.;
const double Tx0=301.,Txn=T0, Ty0=T0,Tyn=T0,Tz0=T0,Tzn=T0;
const int nx=50,ny=50,nz=50,mn1=8,mn2=8;//nx-x ,n1-theta ,n2-phi,nz-z,ny-y,nz-z
const double llsi=260.4e-9;
const double para3=1/llsi;
const double retime=llsi/v;
const double kbulk=C*v*llsi/3.;
const double kn=10.;
const double Lx=llsi/kn;
const double deltax=Lx/nx,deltay=deltax,deltaz=deltax;
const double Ly=ny*deltay,Lz=nz*deltaz;

const double delta1=pi/2./mn1,delta2=pi/2./mn2;
const double ITx0=C*v*Tx0/4./pi,ITxn=C*v*T0/4./pi, ITy0=C*v*T0/4./pi,ITyn=C*v*T0/4./pi,ITz0=C*v*T0/4./pi,ITzn=C*v*T0/4./pi;

const int nr=3;
const double deltaxh=deltax*nr;
const double volumeSpace=deltax*deltay*deltaz;
const int gln=4;

struct NodePara{
  double nodeX;
  double nodeY;
  double nodeZ;
  double temp;
  double ltemp;
  double para11;
  double para12;
  double para13;
  double para22;
  double para23;
  double para33;
  double cpara;
  double vgpara;
  int HN;

  double vgradientI;

  double AD;
};

struct Hlist{
  NodePara* Hnode;
  Hlist* next;
};

struct realNode{
  double nodeX;
  double nodeY;
  double nodeZ;
  double I0;
  double I00;
  double T;
  double Tl;
  int HN;

  double I0ll;
  double I0im;
  double I0b;

  double qx;
  double qy;
  double qz;

  double qxim;
  double qyim;
  double qzim;

  double qfx;
  double qfy;
  double qfz;

  double gradientTx;
  double gradientTy;
  double gradientTz;
};
struct realHlist{
  int nodeNumX;
  int nodeNumY;
  int nodeNumZ;

  double XiX;
  double XiY;
  double XiZ;

  realHlist* next;
};


double simpson_qxx(double *q[],int ny,int nz,double deltay,double deltaz){
  double sum1[ny+1],sum2;
  for(int i=0;i<ny+1;i++){
    sum1[i]=simpson(q[i],nz,deltaz);
  }
  sum2=simpson(sum1,ny,deltay);
  return sum2;
}

double simpson_meanValue_qx(realNode **rn[],int ny,int nz,double deltay,double deltaz){
  double sum1[ny+1],sum2,q1[nz+1];
  for(int i=0;i<ny+1;i++){
    for(int j=0;j<nz+1;j++){
      q1[j]=rn[i][j]->qx;
    }
    sum1[i]=simpson(q1,nz,deltaz);
  }
  sum2=simpson(sum1,ny,deltay);
  return sum2;
}

double simpson_meanValue_T(realNode **rn[],int ny,int nz,double deltay,double deltaz){
  double sum1[ny+1],sum2,q1[nz+1];
  for(int i=0;i<ny+1;i++){
    for(int j=0;j<nz+1;j++){
      q1[j]=rn[i][j]->T;
    }
    sum1[i]=simpson(q1,nz,deltaz);
  }
  sum2=simpson(sum1,ny,deltay);
  return sum2;
}

void out_file3(string file,double **T[],double qxx[],double k[],double x[]){
    ofstream outfile;
    outfile.open(file);
    outfile<<"x,"<<"qbar,"<<"k,"<<"q"<<endl;
    for(int i=0;i<nx+1;i++){
      outfile<<showpoint<<setprecision(8)<<x[i]<<","<<showpoint<<setprecision(8)<<qxx[i]<<","<<showpoint<<setprecision(8)<<k[i]<<endl;
      for(int j=0;j<ny+1;j++){
        for(int k=0;k<nz;k++){
          outfile<<showpoint<<setprecision(10)<<T[i][j][k]<<",";
        }
        outfile<<showpoint<<setprecision(10)<<T[i][j][nz]<<endl;
      }
    }
    outfile.close();
}
void out_file_T(string file,realNode ***rn[],double qxx[],double k[],double x[]){
    ofstream outfile;
    outfile.open(file);
    outfile<<"x,"<<"q,"<<"k,"<<"T"<<endl;
    for(int i=0;i<nx+1;i++){
      outfile<<showpoint<<setprecision(8)<<x[i]<<","<<showpoint<<setprecision(8)<<qxx[i]<<","<<showpoint<<setprecision(8)<<k[i]<<endl;
      for(int j=0;j<ny+1;j++){
        for(int k=0;k<nz;k++){
          outfile<<showpoint<<setprecision(10)<<rn[i][j][k]->T<<",";
        }
        outfile<<showpoint<<setprecision(10)<<rn[i][j][nz]->T<<endl;
      }
    }
    outfile.close();
}

void makepara_mpi(double ****para11[],double ****para12[],double ****para13[],double ****para22[],double ****para23[],double ****para33[],double sintheta[],double costheta[],double sinu[],double cosu[],double ****cpara[],double ****ADpara[],double ****para11_1[],double ****para12_1[],double ****para13_1[],double ****para22_1[],double ****para23_1[],double ****para33_1[],int ****dj[],int n1,int n2){
  double volume=volumeSpace;
  int nrr;
  double xi=0;
  double pdIxD=0,pdIyD=0,pdIzD=0;
  double pdIx=0,pdIy=0,pdIz=0;
  double AD=0;
  for(int i=0;i<nx+1;i++){
      for(int j=0;j<ny+1;j++){
        for(int k=0;k<nz+1;k++){
          for(int l=0;l<gln*n2;l++){
            for(int m=0;m<gln*n1;m++){
              para11[i][j][k][l][m]=0;para12[i][j][k][l][m]=0;para13[i][j][k][l][m]=0;para22[i][j][k][l][m]=0;para23[i][j][k][l][m]=0;para33[i][j][k][l][m]=0;
              pdIxD=0,pdIyD=0,pdIzD=0;pdIx=0,pdIy=0,pdIz=0;

                for(int pr1=0;pr1<2*nr+1;pr1++){
                  for(int pr2=0;pr2<2*nr+1;pr2++){
                    for(int pr3=0;pr3<2*nr+1;pr3++){
                      if(pr1==nr&&pr2==nr&&pr3==nr){
                      }else if((nr-pr1)*(nr-pr1)+(nr-pr2)*(nr-pr2)+(nr-pr3)*(nr-pr3)<=nr*nr){
                        if((nr-pr1)*costheta[m]+(nr-pr2)*sintheta[m]*cosu[l]+(nr-pr3)*sintheta[m]*sinu[l]>0.){
                          if(i+pr1-nr>=0&&i+pr1-nr<=nx&&j+pr2-nr>=0&&j+pr2-nr<=ny&&k+pr3-nr>=0&&k+pr3-nr<=nz){
                            xi=volume/((pr1-nr)*(pr1-nr)+(pr2-nr)*(pr2-nr)+(pr3-nr)*(pr3-nr));

                            para11[i][j][k][l][m]+=(pr1-nr)*(pr1-nr)*xi;
                            para12[i][j][k][l][m]+=(pr1-nr)*(pr2-nr)*xi;
                            para13[i][j][k][l][m]+=(pr1-nr)*(pr3-nr)*xi;
                            para22[i][j][k][l][m]+=(pr2-nr)*(pr2-nr)*xi;
                            para23[i][j][k][l][m]+=(pr2-nr)*(pr3-nr)*xi;
                            para33[i][j][k][l][m]+=(pr3-nr)*(pr3-nr)*xi;

                            pdIx+=(pr1-nr)*xi/deltax;
                            pdIy+=(pr2-nr)*xi/deltax;
                            pdIz+=(pr3-nr)*xi/deltax;
                          }
                        }
                      }
                    }
                  }
                }
              if(para22[i][j][k][l][m]==0.&&para33[i][j][k][l][m]==0.){
                cpara[i][j][k][l][m]=1./(costheta[m]*pdIx/para11[i][j][k][l][m]-1./llsi);
                ADpara[i][j][k][l][m]=0.;

                para11_1[i][j][k][l][m]=1./para11[i][j][k][l][m];
              }else if(para11[i][j][k][l][m]==0.&&para33[i][j][k][l][m]==0.){
                cpara[i][j][k][l][m]=1./(sintheta[m]*cosu[l]*(pdIy/para22[i][j][k][l][m])-1./llsi);
                ADpara[i][j][k][l][m]=0.;

                para22_1[i][j][k][l][m]=1./para22[i][j][k][l][m];
              }else if(para11[i][j][k][l][m]==0.&&para22[i][j][k][l][m]==0.){
                cpara[i][j][k][l][m]=1./(sintheta[m]*sinu[l]*(pdIz/para33[i][j][k][l][m])-1./llsi);
                ADpara[i][j][k][l][m]=0.;

                para33_1[i][j][k][l][m]=1./para33[i][j][k][l][m];
              }else if(para33[i][j][k][l][m]==0.){
                double paraA[2][2];
                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=para12[i][j][k][l][m];
                paraA[1][0]=para12[i][j][k][l][m];paraA[1][1]=para22[i][j][k][l][m];
                ADpara[i][j][k][l][m]=deter22(paraA);

                paraA[0][0]=pdIx;paraA[0][1]=para12[i][j][k][l][m];
                paraA[1][0]=pdIy;paraA[1][1]=para22[i][j][k][l][m];
                pdIxD=deter22(paraA);

                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=pdIx;
                paraA[1][0]=para12[i][j][k][l][m];paraA[1][1]=pdIy;
                pdIyD=deter22(paraA);

                cpara[i][j][k][l][m]=1./(((costheta[m]*pdIxD+sintheta[m]*cosu[l]*pdIyD)/ADpara[i][j][k][l][m])-1./llsi);

                if(para12[i][j][k][l][m]==0.){
                  dj[i][j][k][l][m]=1;
                  para11_1[i][j][k][l][m]=1./para11[i][j][k][l][m];
                  para22_1[i][j][k][l][m]=1./para22[i][j][k][l][m];
                  para12_1[i][j][k][l][m]=0.;
                }else{
                  dj[i][j][k][l][m]=0;
                  para11_1[i][j][k][l][m]=para22[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para22_1[i][j][k][l][m]=para11[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para12_1[i][j][k][l][m]=para12[i][j][k][l][m]/ADpara[i][j][k][l][m];
                }
              }else if(para22[i][j][k][l][m]==0.){
                double paraA[2][2];
                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=para13[i][j][k][l][m];
                paraA[1][0]=para13[i][j][k][l][m];paraA[1][1]=para33[i][j][k][l][m];
                ADpara[i][j][k][l][m]=deter22(paraA);

                paraA[0][0]=pdIx;paraA[0][1]=para13[i][j][k][l][m];
                paraA[1][0]=pdIz;paraA[1][1]=para33[i][j][k][l][m];
                pdIxD=deter22(paraA);

                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=pdIx;
                paraA[1][0]=para13[i][j][k][l][m];paraA[1][1]=pdIz;
                pdIzD=deter22(paraA);

                cpara[i][j][k][l][m]=1./(((costheta[m]*pdIxD+sintheta[m]*sinu[l]*pdIzD)/ADpara[i][j][k][l][m])-1./llsi);

                if(para13[i][j][k][l][m]==0.){
                  dj[i][j][k][l][m]=1;
                  para11_1[i][j][k][l][m]=1./para11[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=1./para33[i][j][k][l][m];
                  para13_1[i][j][k][l][m]=0.;
                }else{
                  dj[i][j][k][l][m]=0;
                  para11_1[i][j][k][l][m]=para33[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=para11[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para13_1[i][j][k][l][m]=para13[i][j][k][l][m]/ADpara[i][j][k][l][m];
                }
              }else if(para11[i][j][k][l][m]==0.){
                double paraA[2][2];
                paraA[0][0]=para22[i][j][k][l][m];paraA[0][1]=para23[i][j][k][l][m];
                paraA[1][0]=para23[i][j][k][l][m];paraA[1][1]=para33[i][j][k][l][m];
                ADpara[i][j][k][l][m]=deter22(paraA);

                paraA[0][0]=pdIy;paraA[0][1]=para23[i][j][k][l][m];
                paraA[1][0]=pdIz;paraA[1][1]=para33[i][j][k][l][m];
                pdIyD=deter22(paraA);

                paraA[0][0]=para22[i][j][k][l][m];paraA[0][1]=pdIy;
                paraA[1][0]=para23[i][j][k][l][m];paraA[1][1]=pdIz;
                pdIzD=deter22(paraA);

                cpara[i][j][k][l][m]=1./(((sintheta[m]*cosu[l]*pdIyD+sintheta[m]*sinu[l]*pdIzD)/ADpara[i][j][k][l][m])-1./llsi);

                if(para23[i][j][k][l][m]==0.){
                  dj[i][j][k][l][m]=1;
                  para22_1[i][j][k][l][m]=1./para22[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=1./para33[i][j][k][l][m];
                  para23_1[i][j][k][l][m]=0.;
                }else{
                  dj[i][j][k][l][m]=0;
                  para22_1[i][j][k][l][m]=para33[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=para22[i][j][k][l][m]/ADpara[i][j][k][l][m];
                  para23_1[i][j][k][l][m]=para23[i][j][k][l][m]/ADpara[i][j][k][l][m];
                }
              }else{
                double paraA[3][3];
                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=para12[i][j][k][l][m];paraA[0][2]=para13[i][j][k][l][m];
                paraA[1][0]=para12[i][j][k][l][m];paraA[1][1]=para22[i][j][k][l][m];paraA[2][2]=para23[i][j][k][l][m];
                paraA[2][0]=para13[i][j][k][l][m];paraA[2][1]=para23[i][j][k][l][m];paraA[2][2]=para33[i][j][k][l][m];
                ADpara[i][j][k][l][m]=deter33(paraA);

                paraA[0][0]=pdIx;paraA[0][1]=para12[i][j][k][l][m];paraA[0][2]=para13[i][j][k][l][m];
                paraA[1][0]=pdIy;paraA[1][1]=para22[i][j][k][l][m];paraA[2][2]=para23[i][j][k][l][m];
                paraA[2][0]=pdIz;paraA[2][1]=para23[i][j][k][l][m];paraA[2][2]=para33[i][j][k][l][m];
                pdIxD=deter33(paraA);

                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=pdIx;paraA[0][2]=para13[i][j][k][l][m];
                paraA[1][0]=para12[i][j][k][l][m];paraA[1][1]=pdIy;paraA[2][2]=para23[i][j][k][l][m];
                paraA[2][0]=para13[i][j][k][l][m];paraA[2][1]=pdIz;paraA[2][2]=para33[i][j][k][l][m];
                pdIyD=deter33(paraA);

                paraA[0][0]=para11[i][j][k][l][m];paraA[0][1]=para12[i][j][k][l][m];paraA[0][2]=pdIx;
                paraA[1][0]=para12[i][j][k][l][m];paraA[1][1]=para22[i][j][k][l][m];paraA[2][2]=pdIy;
                paraA[2][0]=para13[i][j][k][l][m];paraA[2][1]=para23[i][j][k][l][m];paraA[2][2]=pdIz;
                pdIzD=deter33(paraA);

                cpara[i][j][k][l][m]=1./(((costheta[m]*pdIxD+sintheta[m]*cosu[l]*pdIyD+sintheta[m]*sinu[l]*pdIzD)/ADpara[i][j][k][l][m])-1./llsi);

                if(para12[i][j][k][l][m]==0.&&para13[i][j][k][l][m]==0.&&para23[i][j][k][l][m]==0.){
                  dj[i][j][k][l][m]=1;
                  para11_1[i][j][k][l][m]=1./para11[i][j][k][l][m];
                  para22_1[i][j][k][l][m]=1./para22[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=1./para33[i][j][k][l][m];
                  para12_1[i][j][k][l][m]=0.;
                  para13_1[i][j][k][l][m]=0.;
                  para23_1[i][j][k][l][m]=0.;
                }else{
                  dj[i][j][k][l][m]=0;

                  para11_1[i][j][k][l][m]=(para22[i][j][k][l][m]*para33[i][j][k][l][m]-para23[i][j][k][l][m]*para23[i][j][k][l][m])/ADpara[i][j][k][l][m];
                  para22_1[i][j][k][l][m]=(para11[i][j][k][l][m]*para33[i][j][k][l][m]-para13[i][j][k][l][m]*para13[i][j][k][l][m])/ADpara[i][j][k][l][m];
                  para33_1[i][j][k][l][m]=(para22[i][j][k][l][m]*para11[i][j][k][l][m]-para12[i][j][k][l][m]*para12[i][j][k][l][m])/ADpara[i][j][k][l][m];
                  para12_1[i][j][k][l][m]=(para13[i][j][k][l][m]*para23[i][j][k][l][m]-para12[i][j][k][l][m]*para33[i][j][k][l][m])/ADpara[i][j][k][l][m];
                  para13_1[i][j][k][l][m]=(para12[i][j][k][l][m]*para23[i][j][k][l][m]-para22[i][j][k][l][m]*para13[i][j][k][l][m])/ADpara[i][j][k][l][m];
                  para23_1[i][j][k][l][m]=(para12[i][j][k][l][m]*para13[i][j][k][l][m]-para11[i][j][k][l][m]*para23[i][j][k][l][m])/ADpara[i][j][k][l][m];

                  ADpara[i][j][k][l][m]=1./ADpara[i][j][k][l][m];
                }
              }

            }
          }
        }
      }
    }

}

double PD_BTE_tran_GS_3(string type,int xx,int yy,int zz,int l,int m,double ****temp1[],realNode ***rn[],double costheta,double sinthetacosu,double sinthetasinu,realHlist *Hnode,double para11,double para12,double para13,double para22,double para23,double para33,double cpara,double ADpara,double para11_1,double para12_1,double para13_1,double para22_1,double para23_1,double para33_1, int dj,double I0ll){
  double temp=0;
  int nrr;
  double volume=volumeSpace;
  
  if((xx==nx&&type=="bff")||(xx==nx&&type=="bfb")||(xx==nx&&type=="bbf")||(xx==nx&&type=="bbb")){
    temp=ITxn;
  }else if((yy==0&&type=="fff")||(yy==0&&type=="ffb")||(yy==0&&type=="bff")||(yy==0&&type=="bfb")){
    temp=ITy0;
  }else if((yy==ny&&type=="fbf")||(yy==ny&&type=="fbb")||(yy==ny&&type=="bbf")||(yy==ny&&type=="bbb")){
    temp=ITyn;
  }else if((zz==0&&type=="fff")||(zz==0&&type=="fbf")||(zz==0&&type=="bff")||(zz==0&&type=="bbf")){
    temp=ITz0;
  }else if((zz==nz&&type=="ffb")||(zz==nz&&type=="fbb")||(zz==nz&&type=="bfb")||(zz==nz&&type=="bbb")){
    temp=ITzn;
  }
  else if((xx==0&&type=="fff")||(xx==0&&type=="ffb")||(xx==0&&type=="fbf")||(xx==0&&type=="fbb")){
    temp=ITx0;
  }
  else{
    double vgt=0;
    double dIxD=0,dIyD=0,dIzD=0,pdIxD=0,pdIyD=0,pdIzD=0;
    double dIx=0,dIy=0,dIz=0,pdIx=0,pdIy=0,pdIz=0;
    double xi=0,AD=0;
    realHlist* thisnode=Hnode;
    for(int pr=0;pr<rn[xx][yy][zz]->HN;pr++){
      if((xx-thisnode->nodeNumX)*costheta+(yy-thisnode->nodeNumY)*sinthetacosu+(zz-thisnode->nodeNumZ)*sinthetasinu>0.){
        dIx+=temp1[thisnode->nodeNumX][thisnode->nodeNumY][thisnode->nodeNumZ][l][m]*thisnode->XiX;
        dIy+=temp1[thisnode->nodeNumX][thisnode->nodeNumY][thisnode->nodeNumZ][l][m]*thisnode->XiY;
        dIz+=temp1[thisnode->nodeNumX][thisnode->nodeNumY][thisnode->nodeNumZ][l][m]*thisnode->XiZ;
      }
      thisnode=thisnode->next;
    }

    if(para22==0.&&para33==0.){
      temp=((costheta*dIx*para11_1)-I0ll)*cpara;
    }else if(para11==0.&&para33==0.){
      temp=(sinthetacosu*(dIy*para22_1)-I0ll)*cpara;
    }else if(para11==0.&&para22==0.){
      temp=(sinthetasinu*(dIz*para33_1)-I0ll)*cpara;
    }else if(para33==0.){
      if(dj==1){
        dIxD=dIx*para11_1;
        dIyD=dIy*para22_1;
        temp=(((costheta*dIxD+sinthetacosu*dIyD))-I0ll)*cpara;
      }else{
        dIxD=dIx*para11_1+dIy*para12_1;
        dIyD=dIx*para12_1+dIy*para22_1;
        temp=(((costheta*dIxD+sinthetacosu*dIyD))-I0ll)*cpara;
      }
    }else if(para22==0.){
      if(dj==1){
        dIxD=dIx*para11_1;
        dIzD=dIz*para33_1;
        temp=((costheta*dIxD+sinthetasinu*dIzD)-I0ll)*cpara;
      }else{
        dIxD=dIx*para11_1+dIz*para13_1;
        dIzD=dIx*para13_1+dIz*para33_1;
        temp=((costheta*dIxD+sinthetasinu*dIzD)-I0ll)*cpara;
      }
    }else if(para11==0.){
      if(dj==1){
        dIyD=dIy*para22_1;
        dIzD=dIz*para33_1;
        temp=((sinthetacosu*dIyD+sinthetasinu*dIzD)-I0ll)*cpara;
      }else{
        dIyD=dIy*para22_1+dIz*para23_1;
        dIzD=dIy*para23_1+dIz*para33_1;
        temp=((sinthetacosu*dIyD+sinthetasinu*dIzD)-I0ll)*cpara;
      }
    }else{
      if(dj==1){
        dIxD=dIx*para11_1;
        dIyD=dIy*para22_1;
        dIzD=dIz*para33_1;
        temp=((costheta*dIxD+sinthetacosu*dIyD+sinthetasinu*dIzD)-I0ll)*cpara;
      }else{
        double paraA[3][3];
        paraA[0][0]=dIx;paraA[0][1]=para12;paraA[0][2]=para13;
        paraA[1][0]=dIy;paraA[1][1]=para22;paraA[2][2]=para23;
        paraA[2][0]=dIz;paraA[2][1]=para23;paraA[2][2]=para33;
        dIxD=deter33(paraA);

        paraA[0][0]=para11;paraA[0][1]=dIx;paraA[0][2]=para13;
        paraA[1][0]=para12;paraA[1][1]=dIy;paraA[2][2]=para23;
        paraA[2][0]=para13;paraA[2][1]=dIz;paraA[2][2]=para33;
        dIyD=deter33(paraA);

        paraA[0][0]=para11;paraA[0][1]=para12;paraA[0][2]=dIx;
        paraA[1][0]=para12;paraA[1][1]=para22;paraA[2][2]=dIy;
        paraA[2][0]=para13;paraA[2][1]=para23;paraA[2][2]=dIz;
        dIzD=deter33(paraA);

        temp=(((costheta*dIxD+sinthetacosu*dIyD+sinthetasinu*dIzD)*ADpara)-I0ll)*cpara;
      }

    }
  }
  return temp;
}

double ND(int x,int y,double L,int c1,int c2){
  double s=(exp(-1.*0.5*((deltax*(x-c1)*deltax*(x-c1)/L/L)+(deltax*(y-c2)*deltax*(y-c2)/L/L))));
  return s;
}

void initializeRealNodes(realNode ***rn[],realHlist ***Hrn[]){
  int nrr;
  realHlist *lastnode;
  realHlist *thisnode;
    for(int i=0;i<nx+1;i++){
      for(int j=0;j<ny+1;j++){
        for(int k=0;k<nz+1;k++){
          rn[i][j][k]->HN=0;
   
                for(int pr1=0;pr1<2*nr+1;pr1++){
                  for(int pr2=0;pr2<2*nr+1;pr2++){
                    for(int pr3=0;pr3<2*nr+1;pr3++){
                      if(pr1==nr&&pr2==nr&&pr3==nr){
                      }else if((nr-pr1)*(nr-pr1)+(nr-pr2)*(nr-pr2)+(nr-pr3)*(nr-pr3)<=nr*nr){
                          if(i+pr1-nr>=0&&i+pr1-nr<=nx&&j+pr2-nr>=0&&j+pr2-nr<=ny&&k+pr3-nr>=0&&k+pr3-nr<=nz){

                            if(rn[i][j][k]->HN==0){
                              Hrn[i][j][k]=new realHlist;
                              Hrn[i][j][k]->nodeNumX=i+pr1-nr;
                              Hrn[i][j][k]->nodeNumY=j+pr2-nr;
                              Hrn[i][j][k]->nodeNumZ=k+pr3-nr;

                              Hrn[i][j][k]->XiX=volumeSpace*(Hrn[i][j][k]->nodeNumX-i)/(((Hrn[i][j][k]->nodeNumX-i)*(Hrn[i][j][k]->nodeNumX-i)+(Hrn[i][j][k]->nodeNumY-j)*(Hrn[i][j][k]->nodeNumY-j)+(Hrn[i][j][k]->nodeNumZ-k)*(Hrn[i][j][k]->nodeNumZ-k))*deltax);
                              Hrn[i][j][k]->XiY=volumeSpace*(Hrn[i][j][k]->nodeNumY-j)/(((Hrn[i][j][k]->nodeNumX-i)*(Hrn[i][j][k]->nodeNumX-i)+(Hrn[i][j][k]->nodeNumY-j)*(Hrn[i][j][k]->nodeNumY-j)+(Hrn[i][j][k]->nodeNumZ-k)*(Hrn[i][j][k]->nodeNumZ-k))*deltax);
                              Hrn[i][j][k]->XiZ=volumeSpace*(Hrn[i][j][k]->nodeNumZ-k)/(((Hrn[i][j][k]->nodeNumX-i)*(Hrn[i][j][k]->nodeNumX-i)+(Hrn[i][j][k]->nodeNumY-j)*(Hrn[i][j][k]->nodeNumY-j)+(Hrn[i][j][k]->nodeNumZ-k)*(Hrn[i][j][k]->nodeNumZ-k))*deltax);

                              lastnode=Hrn[i][j][k];
                            }else{
                              thisnode=new realHlist;
                              thisnode->nodeNumX=i+pr1-nr;
                              thisnode->nodeNumY=j+pr2-nr;
                              thisnode->nodeNumZ=k+pr3-nr;

                              thisnode->XiX=volumeSpace*(thisnode->nodeNumX-i)/(((thisnode->nodeNumX-i)*(thisnode->nodeNumX-i)+(thisnode->nodeNumY-j)*(thisnode->nodeNumY-j)+(thisnode->nodeNumZ-k)*(thisnode->nodeNumZ-k))*deltax);
                              thisnode->XiY=volumeSpace*(thisnode->nodeNumY-j)/(((thisnode->nodeNumX-i)*(thisnode->nodeNumX-i)+(thisnode->nodeNumY-j)*(thisnode->nodeNumY-j)+(thisnode->nodeNumZ-k)*(thisnode->nodeNumZ-k))*deltax);
                              thisnode->XiZ=volumeSpace*(thisnode->nodeNumZ-k)/(((thisnode->nodeNumX-i)*(thisnode->nodeNumX-i)+(thisnode->nodeNumY-j)*(thisnode->nodeNumY-j)+(thisnode->nodeNumZ-k)*(thisnode->nodeNumZ-k))*deltax);

                              lastnode->next=thisnode;
                              lastnode=thisnode;
                            }
                            rn[i][j][k]->HN++;
                          }
                      }
                    }
                  }
                }

        }
      }
    }
}

int main(int argc,char** argv){
    int iam,iam2,np,np2;
    MPI_Comm comm;

    MPI_Init(&argc,&argv);
    MPI_Comm_dup(MPI_COMM_WORLD,&comm);
    MPI_Comm_rank(comm,&iam2);
    MPI_Comm_size(comm,&np2);
    MPI_Status st;

    int tag=1;

    int n1,n2;
    int sn1,sn2,msn1,msn2;

    if(np2==(mn1*mn2)*2){
      np=np2/2;
      iam=iam2/2;
    }else{
      np=np2;
      iam=iam2;
    }
    
    if(np<mn2){
        n1=mn1;sn1=0;msn1=mn1-1;
        if((mn2/np)*np==mn2){
            n2=(mn2/np);
            sn2=iam*(mn2/np);
        }else if(iam<mn2%np){
            n2=(mn2/np)+1;
            sn2=iam*((mn2/np)+1);
        }else{
            n2=(mn2/np);
            sn2=iam*((mn2/np)+1)-(iam-mn2%np);
        }
        msn2=sn2+n2-1;
    }else{
        n2=1;
        if(np<=mn1*mn2){
            if(np%mn2==0){
                sn2=iam/(np/mn2);msn2=sn2;
                if((mn1*mn2/np)*np==mn1*mn2){
                    n1=mn1*mn2/np;
                    sn1=(iam%(np/mn2))*(mn1*mn2/np);
                }else if(iam%(np/mn2)<(mn1%(np/mn2))){
                    n1=(mn1*mn2/np)+1;
                    sn1=(iam%(np/mn2))*((mn1*mn2/np)+1);
                }else{
                    n1=mn1*mn2/np;
                    sn1=(iam%(np/mn2))*((mn1*mn2/np)+1)-(iam%(np/mn2)-mn1%(np/mn2));
                }
                msn1=sn1+n1-1;
            }else{
              cout<<"error:Please manually assign angles to each process."<<endl;
              return 0;
            }
        }else{
            n1=1;
            n2=1;

            cout<<"Warning: np > mn1*mn2, The output may not be correct."<<endl;
        }
    }
    np=np2;
    iam=iam2;
	if(iam==0)
		cout<<np<<endl;

    clock_t start,finish;
    double totaltime;

    double beta,residual=10.;
    double TT=0;

    realNode ***rn[nx+1];
    realHlist ***Hrn[nx+1];

    double **qx2[nx+1],**qy2[nx+1],**qz2[nx+1];

    double *Isfffd[gln*n2],*Isffbd[gln*n2],*Isfbfd[gln*n2],*Isfbbd[gln*n2],*Isbffd[gln*n2],*Isbfbd[gln*n2],*Isbbfd[gln*n2],*Isbbbd[gln*n2];
    double theta1[gln*n1],theta2[gln*n1],u1[gln*n2],u2[gln*n2],u3[gln*n2],u4[gln*n2];
    double sintheta1[gln*n1],sintheta2[gln*n1],sinu1[gln*n2],sinu2[gln*n2],sinu3[gln*n2],sinu4[gln*n2];
    double costheta1[gln*n1],costheta2[gln*n1],cosu1[gln*n2],cosu2[gln*n2],cosu3[gln*n2],cosu4[gln*n2];
    double sincostheta1[gln*n1],sincostheta2[gln*n1],sinsintheta1[gln*n1],sinsintheta2[gln*n1];

    double sintheta1cosu1[gln*n2][gln*n1],sintheta1cosu2[gln*n2][gln*n1],sintheta1cosu3[gln*n2][gln*n1],sintheta1cosu4[gln*n2][gln*n1];
    double sintheta2cosu1[gln*n2][gln*n1],sintheta2cosu2[gln*n2][gln*n1],sintheta2cosu3[gln*n2][gln*n1],sintheta2cosu4[gln*n2][gln*n1];
    double sintheta1sinu1[gln*n2][gln*n1],sintheta1sinu2[gln*n2][gln*n1],sintheta1sinu3[gln*n2][gln*n1],sintheta1sinu4[gln*n2][gln*n1];
    double sintheta2sinu1[gln*n2][gln*n1],sintheta2sinu2[gln*n2][gln*n1],sintheta2sinu3[gln*n2][gln*n1],sintheta2sinu4[gln*n2][gln*n1];

    double x[nx+1],y[ny+1],z[nz+1],kh[nx+1],qxx[nx+1];

    double **I0[nx+1],**I00[nx+1],**I0im[nx+1],**I0ll[nx+1];
    double *I0ll_1[nx+1],*I0im_1[nx+1];

    double **qxim[nx+1],**qyim[nx+1],**qzim[nx+1];
    double *qxim_1[nx+1],*qyim_1[nx+1],*qzim_1[nx+1];

    string file0;

    if(iam==0){
      start=clock();  //timing
    }
    

    for(int i=0;i<nx+1;i++){
      rn[i]=new realNode**[ny+1];
      Hrn[i]=new realHlist**[ny+1];
      qx2[i]=new double*[ny+1];qy2[i]=new double*[ny+1];qz2[i]=new double*[ny+1];

      I0[i]=new double*[ny+1];I00[i]=new double*[ny+1];
      
      I0im[i]=new double*[ny+1];I0ll[i]=new double*[ny+1];qxim[i]=new double*[ny+1];qyim[i]=new double*[ny+1];qzim[i]=new double*[ny+1];

      for(int j=0;j<ny+1;j++){
        I0[i][j]=new double[nz+1];I00[i][j]=new double[nz+1];
        
        I0im[i][j]=new double[nz+1];I0ll[i][j]=new double[nz+1];qxim[i][j]=new double[nz+1];qyim[i][j]=new double[nz+1];qzim[i][j]=new double[nz+1];        
        
        rn[i][j]=new realNode*[nz+1];
        Hrn[i][j]=new realHlist*[nz+1];
        qx2[i][j]=new double[nz+1];qy2[i][j]=new double[nz+1];qz2[i][j]=new double[nz+1];

      }
    }

    double *I0b[nx+1];double *qxb[nx+1];double *qyb[nx+1];double *qzb[nx+1];
      for(int i=0;i<nx+1;i++){
        I0b[i]=new double[np*(nz+1)*(ny+1)];qxb[i]=new double[np*(nz+1)*(ny+1)];qyb[i]=new double[np*(nz+1)*(ny+1)];qzb[i]=new double[np*(nz+1)*(ny+1)];

        I0ll_1[i]=new double[(nz+1)*(ny+1)];I0im_1[i]=new double[(nz+1)*(ny+1)];qxim_1[i]=new double[(nz+1)*(ny+1)];qyim_1[i]=new double[(nz+1)*(ny+1)];qzim_1[i]=new double[(nz+1)*(ny+1)];
      }

    double ****tempfff[nx+1];
      double ****ltempfff[nx+1];
      double ****para11fff[nx+1],****para22fff[nx+1],****para33fff[nx+1],****para12fff[nx+1],****para13fff[nx+1],****para23fff[nx+1];
      double ****para11fff1[nx+1],****para22fff1[nx+1],****para33fff1[nx+1],****para12fff1[nx+1],****para13fff1[nx+1],****para23fff1[nx+1];
      double ****cparafff[nx+1],****ADparafff[nx+1];
      int ****djfff[nx+1];

      for(int i=0;i<nx+1;i++){
        tempfff[i]=new double***[ny+1];
        ltempfff[i]=new double***[ny+1];

        para11fff[i]=new double***[ny+1];para22fff[i]=new double***[ny+1];para33fff[i]=new double***[ny+1];para12fff[i]=new double***[ny+1];para13fff[i]=new double***[ny+1];para23fff[i]=new double***[ny+1];
        cparafff[i]=new double***[ny+1];
        ADparafff[i]=new double***[ny+1];

        para11fff1[i]=new double***[ny+1];para22fff1[i]=new double***[ny+1];para33fff1[i]=new double***[ny+1];para12fff1[i]=new double***[ny+1];para13fff1[i]=new double***[ny+1];para23fff1[i]=new double***[ny+1];
        djfff[i]=new int***[ny+1];

        for(int j=0;j<ny+1;j++){
          tempfff[i][j]=new double**[nz+1];
          ltempfff[i][j]=new double**[nz+1];

          para11fff[i][j]=new double**[nz+1];para22fff[i][j]=new double**[nz+1];para33fff[i][j]=new double**[nz+1];para12fff[i][j]=new double**[nz+1];para13fff[i][j]=new double**[nz+1];para23fff[i][j]=new double**[nz+1];
          cparafff[i][j]=new double**[nz+1];
          ADparafff[i][j]=new double**[nz+1];

          para11fff1[i][j]=new double**[nz+1];para22fff1[i][j]=new double**[nz+1];para33fff1[i][j]=new double**[nz+1];para12fff1[i][j]=new double**[nz+1];para13fff1[i][j]=new double**[nz+1];para23fff1[i][j]=new double**[nz+1];
          djfff[i][j]=new int**[nz+1];
          for(int k=0;k<nz+1;k++){
            tempfff[i][j][k]=new double*[gln*n2];
            ltempfff[i][j][k]=new double*[gln*n2];

            para11fff[i][j][k]=new double*[gln*n2];para22fff[i][j][k]=new double*[gln*n2];para33fff[i][j][k]=new double*[gln*n2];para12fff[i][j][k]=new double*[gln*n2];para13fff[i][j][k]=new double*[gln*n2];para23fff[i][j][k]=new double*[gln*n2];
            cparafff[i][j][k]=new double*[gln*n2];
            ADparafff[i][j][k]=new double*[gln*n2];

            para11fff1[i][j][k]=new double*[gln*n2];para22fff1[i][j][k]=new double*[gln*n2];para33fff1[i][j][k]=new double*[gln*n2];para12fff1[i][j][k]=new double*[gln*n2];para13fff1[i][j][k]=new double*[gln*n2];para23fff1[i][j][k]=new double*[gln*n2];
            djfff[i][j][k]=new int*[gln*n2];
            for(int l=0;l<gln*n2;l++){
              tempfff[i][j][k][l]=new double[gln*n1];
              ltempfff[i][j][k][l]=new double[gln*n1];

              para11fff[i][j][k][l]=new double[gln*n1];para22fff[i][j][k][l]=new double[gln*n1];para33fff[i][j][k][l]=new double[gln*n1];para12fff[i][j][k][l]=new double[gln*n1];para13fff[i][j][k][l]=new double[gln*n1];para23fff[i][j][k][l]=new double[gln*n1];
              cparafff[i][j][k][l]=new double[gln*n1];
              ADparafff[i][j][k][l]=new double[gln*n1];

              para11fff1[i][j][k][l]=new double[gln*n1];para22fff1[i][j][k][l]=new double[gln*n1];para33fff1[i][j][k][l]=new double[gln*n1];para12fff1[i][j][k][l]=new double[gln*n1];para13fff1[i][j][k][l]=new double[gln*n1];para23fff1[i][j][k][l]=new double[gln*n1];
              djfff[i][j][k][l]=new int[gln*n1];
            }
          }
        }
      }

      double ****tempbff[nx+1];
      double ****ltempbff[nx+1];
      double ****para11bff[nx+1],****para22bff[nx+1],****para33bff[nx+1],****para12bff[nx+1],****para13bff[nx+1],****para23bff[nx+1];

      double ****para11bff1[nx+1],****para22bff1[nx+1],****para33bff1[nx+1],****para12bff1[nx+1],****para13bff1[nx+1],****para23bff1[nx+1];

      double ****cparabff[nx+1],****ADparabff[nx+1];
      int ****djbff[nx+1];
      for(int i=0;i<nx+1;i++){
        tempbff[i]=new double***[ny+1];
        ltempbff[i]=new double***[ny+1];

        para11bff[i]=new double***[ny+1];para22bff[i]=new double***[ny+1];para33bff[i]=new double***[ny+1];para12bff[i]=new double***[ny+1];para13bff[i]=new double***[ny+1];para23bff[i]=new double***[ny+1];
        cparabff[i]=new double***[ny+1];
        ADparabff[i]=new double***[ny+1];

        para11bff1[i]=new double***[ny+1];para22bff1[i]=new double***[ny+1];para33bff1[i]=new double***[ny+1];para12bff1[i]=new double***[ny+1];para13bff1[i]=new double***[ny+1];para23bff1[i]=new double***[ny+1];
        djbff[i]=new int***[ny+1];

        for(int j=0;j<ny+1;j++){
          tempbff[i][j]=new double**[nz+1];
          ltempbff[i][j]=new double**[nz+1];

          para11bff[i][j]=new double**[nz+1];para22bff[i][j]=new double**[nz+1];para33bff[i][j]=new double**[nz+1];para12bff[i][j]=new double**[nz+1];para13bff[i][j]=new double**[nz+1];para23bff[i][j]=new double**[nz+1];
          cparabff[i][j]=new double**[nz+1];
          ADparabff[i][j]=new double**[nz+1];

          para11bff1[i][j]=new double**[nz+1];para22bff1[i][j]=new double**[nz+1];para33bff1[i][j]=new double**[nz+1];para12bff1[i][j]=new double**[nz+1];para13bff1[i][j]=new double**[nz+1];para23bff1[i][j]=new double**[nz+1];
          djbff[i][j]=new int**[nz+1];

          for(int k=0;k<nz+1;k++){
            tempbff[i][j][k]=new double*[gln*n2];
            ltempbff[i][j][k]=new double*[gln*n2];

            para11bff[i][j][k]=new double*[gln*n2];para22bff[i][j][k]=new double*[gln*n2];para33bff[i][j][k]=new double*[gln*n2];para12bff[i][j][k]=new double*[gln*n2];para13bff[i][j][k]=new double*[gln*n2];para23bff[i][j][k]=new double*[gln*n2];
            cparabff[i][j][k]=new double*[gln*n2];
            ADparabff[i][j][k]=new double*[gln*n2];

            para11bff1[i][j][k]=new double*[gln*n2];para22bff1[i][j][k]=new double*[gln*n2];para33bff1[i][j][k]=new double*[gln*n2];para12bff1[i][j][k]=new double*[gln*n2];para13bff1[i][j][k]=new double*[gln*n2];para23bff1[i][j][k]=new double*[gln*n2];
            djbff[i][j][k]=new int*[gln*n2];
            for(int l=0;l<gln*n2;l++){
              tempbff[i][j][k][l]=new double[gln*n1];
              ltempbff[i][j][k][l]=new double[gln*n1];

              para11bff[i][j][k][l]=new double[gln*n1];para22bff[i][j][k][l]=new double[gln*n1];para33bff[i][j][k][l]=new double[gln*n1];para12bff[i][j][k][l]=new double[gln*n1];para13bff[i][j][k][l]=new double[gln*n1];para23bff[i][j][k][l]=new double[gln*n1];
              cparabff[i][j][k][l]=new double[gln*n1];
              ADparabff[i][j][k][l]=new double[gln*n1];

              para11bff1[i][j][k][l]=new double[gln*n1];para22bff1[i][j][k][l]=new double[gln*n1];para33bff1[i][j][k][l]=new double[gln*n1];para12bff1[i][j][k][l]=new double[gln*n1];para13bff1[i][j][k][l]=new double[gln*n1];para23bff1[i][j][k][l]=new double[gln*n1];
              djbff[i][j][k][l]=new int[gln*n1];
            }
          }
        }
      }

    for(int i=0;i<gln*n2;i++){
      Isfffd[i]=new double[gln*n1];Isffbd[i]=new double[gln*n1];Isfbfd[i]=new double[gln*n1];Isfbbd[i]=new double[gln*n1];
      Isbffd[i]=new double[gln*n1];Isbfbd[i]=new double[gln*n1];Isbbfd[i]=new double[gln*n1];Isbbbd[i]=new double[gln*n1];
    }

    double aa=0,bb=0,maxx=(double)(msn1+1)*(pi/2.)/mn1,minx=sn1*(pi/2.)/mn1;
	  double deltaGL=(maxx-minx)/(double)n1;
	  for(int i=0;i<n1;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;

		  theta1[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;theta1[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;theta1[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;theta1[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sintheta1[4*i]=sin(theta1[4*i]);sintheta1[4*i+1]=sin(theta1[4*i+1]);sintheta1[4*i+2]=sin(theta1[4*i+2]);sintheta1[4*i+3]=sin(theta1[4*i+3]);
      costheta1[4*i]=cos(theta1[4*i]);costheta1[4*i+1]=cos(theta1[4*i+1]);costheta1[4*i+2]=cos(theta1[4*i+2]);costheta1[4*i+3]=cos(theta1[4*i+3]);
      sincostheta1[4*i]=sintheta1[4*i]*costheta1[4*i];sincostheta1[4*i+1]=sintheta1[4*i+1]*costheta1[4*i+1];sincostheta1[4*i+2]=sintheta1[4*i+2]*costheta1[4*i+2];sincostheta1[4*i+3]=sintheta1[4*i+3]*costheta1[4*i+3];
      sinsintheta1[4*i]=sintheta1[4*i]*sintheta1[4*i];sinsintheta1[4*i+1]=sintheta1[4*i+1]*sintheta1[4*i+1];sinsintheta1[4*i+2]=sintheta1[4*i+2]*sintheta1[4*i+2];sinsintheta1[4*i+3]=sintheta1[4*i+3]*sintheta1[4*i+3];

    }

    aa=0;bb=0;maxx=pi-(double)(msn1+1)*(pi/2.)/mn1;minx=pi-sn1*(pi/2.)/mn1;
	  deltaGL=(maxx-minx)/(double)n1;
	  for(int i=0;i<n1;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;

      theta2[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;theta2[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;theta2[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;theta2[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sintheta2[4*i]=sin(theta2[4*i]);sintheta2[4*i+1]=sin(theta2[4*i+1]);sintheta2[4*i+2]=sin(theta2[4*i+2]);sintheta2[4*i+3]=sin(theta2[4*i+3]);
      costheta2[4*i]=cos(theta2[4*i]);costheta2[4*i+1]=cos(theta2[4*i+1]);costheta2[4*i+2]=cos(theta2[4*i+2]);costheta2[4*i+3]=cos(theta2[4*i+3]);
      sincostheta2[4*i]=sintheta2[4*i]*costheta2[4*i];sincostheta2[4*i+1]=sintheta2[4*i+1]*costheta2[4*i+1];sincostheta2[4*i+2]=sintheta2[4*i+2]*costheta2[4*i+2];sincostheta2[4*i+3]=sintheta2[4*i+3]*costheta2[4*i+3];
      sinsintheta2[4*i]=sintheta2[4*i]*sintheta2[4*i];sinsintheta2[4*i+1]=sintheta2[4*i+1]*sintheta2[4*i+1];sinsintheta2[4*i+2]=sintheta2[4*i+2]*sintheta2[4*i+2];sinsintheta2[4*i+3]=sintheta2[4*i+3]*sintheta2[4*i+3];
    }

    aa=0;bb=0;maxx=(double)(msn2+1)*(pi/2.)/mn2;minx=sn2*(pi/2.)/mn2;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;

		  u1[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u1[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u1[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u1[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu1[4*i]=sin(u1[4*i]);sinu1[4*i+1]=sin(u1[4*i+1]);sinu1[4*i+2]=sin(u1[4*i+2]);sinu1[4*i+3]=sin(u1[4*i+3]);
      cosu1[4*i]=cos(u1[4*i]);cosu1[4*i+1]=cos(u1[4*i+1]);cosu1[4*i+2]=cos(u1[4*i+2]);cosu1[4*i+3]=cos(u1[4*i+3]);
    }

    aa=0;bb=0;maxx=pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2;minx=pi/2.+(mn2-msn2-1)*(pi/2.)/mn2;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;

		  u2[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u2[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u2[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u2[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu2[4*i]=sin(u2[4*i]);sinu2[4*i+1]=sin(u2[4*i+1]);sinu2[4*i+2]=sin(u2[4*i+2]);sinu2[4*i+3]=sin(u2[4*i+3]);
      cosu2[4*i]=cos(u2[4*i]);cosu2[4*i+1]=cos(u2[4*i+1]);cosu2[4*i+2]=cos(u2[4*i+2]);cosu2[4*i+3]=cos(u2[4*i+3]);
    }

    aa=0;bb=0;maxx=1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2;minx=1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;
      
		  u3[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u3[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u3[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u3[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu3[4*i]=sin(u3[4*i]);sinu3[4*i+1]=sin(u3[4*i+1]);sinu3[4*i+2]=sin(u3[4*i+2]);sinu3[4*i+3]=sin(u3[4*i+3]);
      cosu3[4*i]=cos(u3[4*i]);cosu3[4*i+1]=cos(u3[4*i+1]);cosu3[4*i+2]=cos(u3[4*i+2]);cosu3[4*i+3]=cos(u3[4*i+3]);
    }

    aa=0;bb=0;maxx=2.*pi-(double)(msn2+1)*(pi/2.)/mn2;minx=2.*pi-sn2*(pi/2.)/mn2;
	  deltaGL=(maxx-minx)/(double)n2;
	  for(int i=0;i<n2;i++){
      aa=(double)i*deltaGL+minx;bb=((double)(i+1))*deltaGL+minx;

		  u4[4*i]=(aa+bb)/2.-0.8611363116*(bb-aa)/2.;u4[4*i+1]=(aa+bb)/2.-0.3399810436*(bb-aa)/2.;u4[4*i+2]=(aa+bb)/2.+0.3399810436*(bb-aa)/2.;u4[4*i+3]=(aa+bb)/2.+0.8611363116*(bb-aa)/2.;
      sinu4[4*i]=sin(u4[4*i]);sinu4[4*i+1]=sin(u4[4*i+1]);sinu4[4*i+2]=sin(u4[4*i+2]);sinu4[4*i+3]=sin(u4[4*i+3]);
      cosu4[4*i]=cos(u4[4*i]);cosu4[4*i+1]=cos(u4[4*i+1]);cosu4[4*i+2]=cos(u4[4*i+2]);cosu4[4*i+3]=cos(u4[4*i+3]);
    }

    for(int i=0;i<gln*n2;i++){
      for(int j=0;j<gln*n1;j++){
        sintheta1cosu1[i][j]=cosu1[i]*sintheta1[j];sintheta1cosu2[i][j]=cosu2[i]*sintheta1[j];sintheta1cosu3[i][j]=cosu3[i]*sintheta1[j];sintheta1cosu4[i][j]=cosu4[i]*sintheta1[j];
        sintheta2cosu1[i][j]=cosu1[i]*sintheta2[j];sintheta2cosu2[i][j]=cosu2[i]*sintheta2[j];sintheta2cosu3[i][j]=cosu3[i]*sintheta2[j];sintheta2cosu4[i][j]=cosu4[i]*sintheta2[j];
        sintheta1sinu1[i][j]=sinu1[i]*sintheta1[j];sintheta1sinu2[i][j]=sinu2[i]*sintheta1[j];sintheta1sinu3[i][j]=sinu3[i]*sintheta1[j];sintheta1sinu4[i][j]=sinu4[i]*sintheta1[j];
        sintheta2sinu1[i][j]=sinu1[i]*sintheta2[j];sintheta2sinu2[i][j]=sinu2[i]*sintheta2[j];sintheta2sinu3[i][j]=sinu3[i]*sintheta2[j];sintheta2sinu4[i][j]=sinu4[i]*sintheta2[j];
      }
    }

    for(int i=0;i<nx+1;i++){
      for(int j=0;j<ny+1;j++){
        for(int k=0;k<nz+1;k++){
          rn[i][j][k]=new realNode;
          rn[i][j][k]->nodeX=i*deltax;rn[i][j][k]->nodeY=j*deltay;rn[i][j][k]->nodeZ=k*deltaz;
          rn[i][j][k]->T=T0;

          I0[i][j][k]=C*v*T0/4./pi;
          I00[i][j][k]=I0[i][j][k];
          I0ll[i][j][k]=I0[i][j][k]/llsi;
          I0ll_1[i][j*(nz+1)+k]=I0[i][j][k]/llsi;

          for(int l=0;l<gln*n2;l++){
            for(int m=0;m<gln*n1;m++){

                tempfff[i][j][k][l][m]=I0[i][j][k];
                ltempfff[i][j][k][l][m]=I0[i][j][k];

                tempbff[i][j][k][l][m]=I0[i][j][k];
                ltempbff[i][j][k][l][m]=I0[i][j][k];
            }
          }
        }
      }
    }

    for(int i=0;i<nx+1;i++){
      x[i]=i*deltax;
      qxx[i]=0;
      kh[i]=C*v*llsi/3.;
    }

    for(int i=0;i<ny+1;i++){
      y[i]=i*deltay;
    }
    for(int i=0;i<ny+1;i++){
      z[i]=i*deltaz;
    }

    if(np<(mn1*mn2)*2){
      makepara_mpi(para11fff,para12fff,para13fff,para22fff,para23fff,para33fff,sintheta1,costheta1,sinu1,cosu1,cparafff,ADparafff,para11fff1,para12fff1,para13fff1,para22fff1,para23fff1,para33fff1,djfff,n1,n2);
      makepara_mpi(para11bff,para12bff,para13bff,para22bff,para23bff,para33bff,sintheta2,costheta2,sinu1,cosu1,cparabff,ADparabff,para11bff1,para12bff1,para13bff1,para22bff1,para23bff1,para33bff1,djbff,n1,n2);
    }else{
      if(iam%2==0){
        makepara_mpi(para11fff,para12fff,para13fff,para22fff,para23fff,para33fff,sintheta1,costheta1,sinu1,cosu1,cparafff,ADparafff,para11fff1,para12fff1,para13fff1,para22fff1,para23fff1,para33fff1,djfff,n1,n2);
      }else if(iam%2==1){
        makepara_mpi(para11bff,para12bff,para13bff,para22bff,para23bff,para33bff,sintheta2,costheta2,sinu1,cosu1,cparabff,ADparabff,para11bff1,para12bff1,para13bff1,para22bff1,para23bff1,para33bff1,djbff,n1,n2);
      }
    }

    initializeRealNodes(rn,Hrn);

    if(iam==0){
      finish=clock();
      totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
      cout<<"\n time is "<<totaltime<<" s"<<endl;
    }

     
    int nd=0;
    double rrr=0.;
    double tempfffff=0;
    while(residual>1e-4&&nd<10000){
      rrr=0.;
      for(int i=0;i<nx+1;i++){
          for(int j=0;j<ny+1;j++){
            for(int k=0;k<nz+1;k++){  
              I00[i][j][k]=I0[i][j][k];
            }
          }
        }

        if(np<(mn1*mn2)*2){
          for(int i=0;i<nx+1;i++){
            for(int j=0;j<ny+1;j++){
              for(int k=0;k<nz+1;k++){ 
                for(int l=0;l<gln*n2;l++){
                  for(int m=0;m<gln*n1;m++){
                    tempfff[i][j][k][l][m]=PD_BTE_tran_GS_3("fff",i,j,k,l,m,tempfff,rn,costheta1[m],sintheta1cosu1[l][m],sintheta1sinu1[l][m],Hrn[i][j][k],para11fff[i][j][k][l][m],para12fff[i][j][k][l][m],para13fff[i][j][k][l][m],para22fff[i][j][k][l][m],para23fff[i][j][k][l][m],para33fff[i][j][k][l][m],cparafff[i][j][k][l][m],ADparafff[i][j][k][l][m],para11fff1[i][j][k][l][m],para12fff1[i][j][k][l][m],para13fff1[i][j][k][l][m],para22fff1[i][j][k][l][m],para23fff1[i][j][k][l][m],para33fff1[i][j][k][l][m],djfff[i][j][k][l][m],I0ll_1[i][j*(nz+1)+k]);
                  }
                }
              }
            }
          }

          for(int i=nx;i>=0;i--){
            for(int j=0;j<ny+1;j++){
              for(int k=0;k<nz+1;k++){  
                for(int l=0;l<gln*n2;l++){
                  for(int m=0;m<gln*n1;m++){
                    tempbff[i][j][k][l][m]=PD_BTE_tran_GS_3("bff",i,j,k,l,m,tempbff,rn,costheta2[m],sintheta2cosu1[l][m],sintheta2sinu1[l][m],Hrn[i][j][k],para11bff[i][j][k][l][m],para12bff[i][j][k][l][m],para13bff[i][j][k][l][m],para22bff[i][j][k][l][m],para23bff[i][j][k][l][m],para33bff[i][j][k][l][m],cparabff[i][j][k][l][m],ADparabff[i][j][k][l][m],para11bff1[i][j][k][l][m],para12bff1[i][j][k][l][m],para13bff1[i][j][k][l][m],para22bff1[i][j][k][l][m],para23bff1[i][j][k][l][m],para33bff1[i][j][k][l][m],djbff[i][j][k][l][m],I0ll_1[i][j*(nz+1)+k]);
                  }
                }
              }
            }
          }
        }else{
          if(iam%2==0){
            for(int i=0;i<nx+1;i++){
              for(int j=0;j<ny+1;j++){
                for(int k=0;k<nz+1;k++){ 
                  for(int l=0;l<gln*n2;l++){
                    for(int m=0;m<gln*n1;m++){
                      tempfff[i][j][k][l][m]=PD_BTE_tran_GS_3("fff",i,j,k,l,m,tempfff,rn,costheta1[m],sintheta1cosu1[l][m],sintheta1sinu1[l][m],Hrn[i][j][k],para11fff[i][j][k][l][m],para12fff[i][j][k][l][m],para13fff[i][j][k][l][m],para22fff[i][j][k][l][m],para23fff[i][j][k][l][m],para33fff[i][j][k][l][m],cparafff[i][j][k][l][m],ADparafff[i][j][k][l][m],para11fff1[i][j][k][l][m],para12fff1[i][j][k][l][m],para13fff1[i][j][k][l][m],para22fff1[i][j][k][l][m],para23fff1[i][j][k][l][m],para33fff1[i][j][k][l][m],djfff[i][j][k][l][m],I0ll_1[i][j*(nz+1)+k]);
                    }
                  }
                }
              }
            }
          }else if(iam%2==1){
            for(int i=nx;i>=0;i--){
              for(int j=0;j<ny+1;j++){
                for(int k=0;k<nz+1;k++){  
                  for(int l=0;l<gln*n2;l++){
                    for(int m=0;m<gln*n1;m++){
                      tempbff[i][j][k][l][m]=PD_BTE_tran_GS_3("bff",i,j,k,l,m,tempbff,rn,costheta2[m],sintheta2cosu1[l][m],sintheta2sinu1[l][m],Hrn[i][j][k],para11bff[i][j][k][l][m],para12bff[i][j][k][l][m],para13bff[i][j][k][l][m],para22bff[i][j][k][l][m],para23bff[i][j][k][l][m],para33bff[i][j][k][l][m],cparabff[i][j][k][l][m],ADparabff[i][j][k][l][m],para11bff1[i][j][k][l][m],para12bff1[i][j][k][l][m],para13bff1[i][j][k][l][m],para22bff1[i][j][k][l][m],para23bff1[i][j][k][l][m],para33bff1[i][j][k][l][m],djbff[i][j][k][l][m],I0ll_1[i][j*(nz+1)+k]);
                    }
                  }
                }
              }
            }
          }
        }

        for(int i=0;i<nx+1;i++){
          for(int j=0;j<ny+1;j++){
            for(int k=0;k<nz+1;k++){
              for(int l=0;l<gln*n2;l++){
                for(int m=0;m<gln*n1;m++){
                  Isfffd[l][m]=tempfff[i][j][k][l][m];Isffbd[l][m]=tempfff[i][j][nz-k][l][m];Isfbfd[l][m]=tempfff[i][ny-j][k][gln*n2-1-l][m];Isfbbd[l][m]=tempfff[i][ny-j][nz-k][gln*n2-1-l][m];
                  Isbffd[l][m]=tempbff[i][j][k][l][m];Isbfbd[l][m]=tempbff[i][j][nz-k][l][m];Isbbfd[l][m]=tempbff[i][ny-j][k][gln*n2-1-l][m];Isbbbd[l][m]=tempbff[i][ny-j][nz-k][gln*n2-1-l][m];
                }
              }

              if(np<(mn1*mn2)*2){
                I0im_1[i][j*(nz+1)+k]=(
                  GL4_I0(Isfffd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isffbd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_I0(Isfbfd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isfbbd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbffd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbfbd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbbfd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbbbd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                )/4./pi;
              }else{
                if(iam%2==0){
                  I0im_1[i][j*(nz+1)+k]=(
                    GL4_I0(Isfffd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                    +GL4_I0(Isffbd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                    +GL4_I0(Isfbfd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                    +GL4_I0(Isfbbd,sintheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                  )/4./pi;
                }else if(iam%2==1){
                  I0im_1[i][j*(nz+1)+k]=(
                    GL4_I0(Isbffd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                    +GL4_I0(Isbfbd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                    +GL4_I0(Isbbfd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                    +GL4_I0(Isbbbd,sintheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                  )/4./pi;
                }
              }

            }
          }
        }

      for(int i=0;i<nx+1;i++){
        MPI_Gather(I0im_1[i],(ny+1)*(nz+1),MPI_DOUBLE,I0b[i],(ny+1)*(nz+1),MPI_DOUBLE,0,comm);
      }

      if(iam==0){
        for(int i=0;i<nx+1;i++){
          for(int j=0;j<ny+1;j++){
            for(int k=0;k<nz+1;k++){
              I0[i][j][k]=0;
              for(int llq=0;llq<np;llq++){
                I0[i][j][k]+=I0b[i][llq*(nz+1)*(ny+1)+j*(nz+1)+k];
              }

              I0ll_1[i][j*(nz+1)+k]=I0[i][j][k]/llsi;
              rrr=max(abs(I0[i][j][k]-I00[i][j][k]),rrr);
            }
          }
        }
      }
        
        for(int i=0;i<nx+1;i++){
          MPI_Bcast(I0ll_1[i], (ny+1)*(nz+1), MPI_DOUBLE, 0, comm);
        }

        MPI_Bcast(&rrr, 1, MPI_DOUBLE, 0, comm);

      residual=rrr;
      nd++;
      if(iam==0){
        cout<<nd<<endl;
        cout<<residual<<endl;
      }
    }

      for(int i=0;i<nx+1;i++){
        for(int j=0;j<ny+1;j++){
          for(int k=0;k<nz+1;k++){
            for(int l=0;l<gln*n2;l++){
              for(int m=0;m<gln*n1;m++){
                Isfffd[l][m]=tempfff[i][j][k][l][m];Isffbd[l][m]=tempfff[i][j][nz-k][l][m];Isfbfd[l][m]=tempfff[i][ny-j][k][gln*n2-1-l][m];Isfbbd[l][m]=tempfff[i][ny-j][nz-k][gln*n2-1-l][m];
                Isbffd[l][m]=tempbff[i][j][k][l][m];Isbfbd[l][m]=tempbff[i][j][nz-k][l][m];Isbbfd[l][m]=tempbff[i][ny-j][k][gln*n2-1-l][m];Isbbbd[l][m]=tempbff[i][ny-j][nz-k][gln*n2-1-l][m];
              }
            }

            if(np<(mn1*mn2)*2){
              qxim_1[i][j*(nz+1)+k]=GL4_I0(Isfffd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_I0(Isffbd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_I0(Isfbfd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_I0(Isfbbd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                +GL4_I0(Isbffd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_I0(Isbfbd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_I0(Isbbfd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_I0(Isbbbd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);

              qyim_1[i][j*(nz+1)+k]=GL4_q_2(Isfffd,sinsintheta1,cosu1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isffbd,sinsintheta1,cosu4,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_q_2(Isfbfd,sinsintheta1,cosu2,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isfbbd,sinsintheta1,cosu3,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbffd,sinsintheta2,cosu1,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbfbd,sinsintheta2,cosu4,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbbfd,sinsintheta2,cosu2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbbbd,sinsintheta2,cosu3,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);
              
              qzim_1[i][j*(nz+1)+k]=GL4_q_2(Isfffd,sinsintheta1,sinu1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isffbd,sinsintheta1,sinu4,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_q_2(Isfbfd,sinsintheta1,sinu2,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isfbbd,sinsintheta1,sinu3,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbffd,sinsintheta2,sinu1,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbfbd,sinsintheta2,sinu4,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbbfd,sinsintheta2,sinu2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                +GL4_q_2(Isbbbd,sinsintheta2,sinu3,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);
            }else{
              if(iam%2==0){
                 qxim_1[i][j*(nz+1)+k]=GL4_I0(Isfffd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isffbd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_I0(Isfbfd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isfbbd,sincostheta1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                 ;

                qyim_1[i][j*(nz+1)+k]=GL4_q_2(Isfffd,sinsintheta1,cosu1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isffbd,sinsintheta1,cosu4,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isfbfd,sinsintheta1,cosu2,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isfbbd,sinsintheta1,cosu3,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                 ;
                
                qzim_1[i][j*(nz+1)+k]=GL4_q_2(Isfffd,sinsintheta1,sinu1,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isffbd,sinsintheta1,sinu4,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isfbfd,sinsintheta1,sinu2,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isfbbd,sinsintheta1,sinu3,sn1*(pi/2.)/mn1,(double)(msn1+1)*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2)
                  ;
              }else{
                qxim_1[i][j*(nz+1)+k]=
                  GL4_I0(Isbffd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbfbd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbbfd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_I0(Isbbbd,sincostheta2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);

                qyim_1[i][j*(nz+1)+k]=
                  GL4_q_2(Isbffd,sinsintheta2,cosu1,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbfbd,sinsintheta2,cosu4,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbbfd,sinsintheta2,cosu2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbbbd,sinsintheta2,cosu3,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);
                
                qzim_1[i][j*(nz+1)+k]=
                  GL4_q_2(Isbffd,sinsintheta2,sinu1,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,sn2*(pi/2.)/mn2,(double)(msn2+1)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbfbd,sinsintheta2,sinu4,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,2.*pi-(double)(msn2+1)*(pi/2.)/mn2,2.*pi-sn2*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbbfd,sinsintheta2,sinu2,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,pi/2.+(mn2-msn2-1)*(pi/2.)/mn2,pi/2.+(double)(mn2-sn2)*(pi/2.)/mn2,n2)
                  +GL4_q_2(Isbbbd,sinsintheta2,sinu3,pi-(double)(msn1+1)*(pi/2.)/mn1,pi-sn1*(pi/2.)/mn1,n1,1.5*pi-(double)(mn2-sn2)*(pi/2.)/mn2,1.5*pi-(mn2-msn2-1)*(pi/2.)/mn2,n2);
              }
            }
          }
        }
      }

      for(int i=0;i<nx+1;i++){
        MPI_Gather(qxim_1[i],(ny+1)*(nz+1),MPI_DOUBLE,qxb[i],(ny+1)*(nz+1),MPI_DOUBLE,0,comm);
        MPI_Gather(qyim_1[i],(ny+1)*(nz+1),MPI_DOUBLE,qyb[i],(ny+1)*(nz+1),MPI_DOUBLE,0,comm);
        MPI_Gather(qzim_1[i],(ny+1)*(nz+1),MPI_DOUBLE,qzb[i],(ny+1)*(nz+1),MPI_DOUBLE,0,comm);
      }

      if(iam==0){
        for(int i=0;i<nx+1;i++){
          for(int j=0;j<ny+1;j++){
            for(int k=0;k<nz+1;k++){
              rn[i][j][k]->qx=0;
              rn[i][j][k]->qy=0;
              rn[i][j][k]->qz=0;
              for(int llq=0;llq<np;llq++){

                rn[i][j][k]->qx+=qxb[i][llq*(nz+1)*(ny+1)+j*(nz+1)+k];
                rn[i][j][k]->qy+=qyb[i][llq*(nz+1)*(ny+1)+j*(nz+1)+k];
                rn[i][j][k]->qz+=qzb[i][llq*(nz+1)*(ny+1)+j*(nz+1)+k];
              }

              qx2[i][j][k]=3.*rn[i][j][k]->qx*Lx/C/v/llsi/(Tx0-T0);
              qy2[i][j][k]=3.*rn[i][j][k]->qy*Lx/C/v/llsi/(Tx0-T0);
              qz2[i][j][k]=3.*rn[i][j][k]->qz*Lx/C/v/llsi/(Tx0-T0);

              rn[i][j][k]->T=I0[i][j][k]*4.*pi/C/v;
            }
          }
        }

        for(int i=0;i<nx+1;i++){
          qxx[i]=simpson_meanValue_qx(rn[i],ny,nz,deltay,deltaz);
          TT=simpson_meanValue_T(rn[0],ny,nz,deltay,deltaz)-simpson_meanValue_T(rn[nx],ny,nz,deltay,deltaz);
          kh[i]=qxx[i]*Lx/TT;
        }

      file0="3D_T_kn_"+to_string(kn)+"_np_"+to_string(np)+".xml";
      out_file_T(file0,rn,qxx,kh,x);

      file0="3D_qx_kn_"+to_string(kn)+"_np_"+to_string(np)+".xml";
      out_file3(file0,qx2,qxx,kh,x);

      file0="3D_qy_kn_"+to_string(kn)+"_np_"+to_string(np)+".xml";
      out_file3(file0,qy2,qxx,kh,x);

      file0="3D_qz_kn_"+to_string(kn)+"_np_"+to_string(np)+".xml";
      out_file3(file0,qz2,qxx,kh,x);
    

      finish=clock();
      totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
      cout<<"\n time is "<<totaltime<<" s"<<endl;

      }

  MPI_Finalize();
     
}