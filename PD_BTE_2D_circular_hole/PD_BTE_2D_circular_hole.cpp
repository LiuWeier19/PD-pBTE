#include"calFunction.h"

using namespace std;
const double kb = 1.3806488e-23;// Planck constant, J s.
const double h = 6.62606957e-34;// Dirac constant, J s.
const double hbar = 1.054571726e-34;
const double pi =3.1415926535897932384626433832795028841953;
const double C=0.930e6,v=1804.;
const double T0=0.;

const double Tx0=1.,Txn=0.;
const double deltaT=1.;

const int nx=50,n1=8,n2=16,nz=50,n_boundry=200;//nx-x ,n1-theta ,n2-phi,nz-z
const double llsi=260.4e-9;
const double para3=1/llsi;
const double kn=0.1;
const double D=llsi/kn,deltaz=D/nz;
const double delta1=pi/2./n1,delta2=pi/n2,deltax=deltaz;
const double L=nx*deltax;
const int nr=3;
const double areaSpace=deltax*deltax;
const double areaSpace2=deltax;
const double ri=0.5;

const double kbu=C*v*llsi/3.;
const int gln=4;
const int glnn1=gln*n1;
const int glnn2=gln*n2;

const double deltaxh=deltax*nr;

const int ncx=(nx/2)+nr-1,ncz=nz/2,ncr=nx*0.3;
const int xl1=ncx-ncr-nr,xl2=ncx+ncr+nr,zl1=ncz-ncr-nr,zl2=ncz+ncr+nr;

const double cx=ncx*deltax,cz=ncz*deltaz,cr=ncr*deltax;
const double deltaI0=C*v*deltaT/4./pi;

struct realHlist{
  int nodeNumX;
  int nodeNumZ;

  double XiX;
  double XiZ;

  double cXiX;
  double cXiZ;

  realHlist* next;

  double **Xitempff;
  double **Xitempbf;
};
struct realNode{
  double X;
  double Z;

  double area;

  double I0;
  double I00;
  double T;

  double qx;
  double qz;

  double qf;
  double qb;
  double qu;
  double qd;

  int HNnode;int *HrnX;int *HrnZ;int *HrnN;int *HrnN_b;
  double *HrnXiX;double *HrnXiZ;double *cXiX;double *cXiZ;

  double ***Xitempff;
  double ***Xitempbf;

  double **tempff;
  double **tempbf;

  double **paraS2ff;
  double **paraC2ff;
  double **paraSCff;

  double **paraS2bf;
  double **paraC2bf;
  double **paraSCbf;

  double **para1ff;
  double **para1bf;

  double **para21ff;
  double **para22ff;

  double **para21bf;
  double **para22bf;
};

struct realNode_boundry{
  realNode *thisRnode;

  double qdmm;
  int clot[glnn2];
};

void out_file(string file,double T[][nz+1],double *x,int nx,int nz,double *q,double *k){
    ofstream outfile;
    outfile.open(file);
    outfile<<"x,"<<"qxbar,"<<"k,"<<"qx"<<endl;
	for(int i=0;i<nx+2*nr-1;i++){
		outfile<<showpoint<<setprecision(8)<<x[i]<<","<<showpoint<<setprecision(8)<<q[i]<<","<<showpoint<<setprecision(8)<<k[i]<<",";
		for(int j=0;j<nz;j++){
			outfile<<showpoint<<setprecision(10)<<T[i][j]<<",";
		}
     outfile<<showpoint<<setprecision(10)<<T[i][nz]<<endl;
	}
	outfile.close();
}

void out_file_T_ast_2(string file,realNode *rn[][nz+1],realNode_boundry *rn_b_z0[],realNode_boundry *rn_b_zn[],int nx,int nz){
    ofstream outfile;
    outfile.open(file);
    outfile<<"n, "<<"x, "<<"z, "<<"T*"<<endl;
    for(int i=0;i<nx+2*nr-1;i++){
      for(int j=0;j<nz+1;j++){
        outfile<<i*(nz+1)+j<<", "<<rn[i][j]->X/L<<", "<<rn[i][j]->Z/D<<", "<<showpoint<<setprecision(8)<<rn[i][j]->T<<endl;
      }
    }

    for(int i=0;i<n_boundry+1;i++){
      outfile<<(nx+2*nr-1)*(nz+1)+i<<", "<<rn_b_z0[i]->thisRnode->X/L<<", "<<rn_b_z0[i]->thisRnode->Z/D<<", "<<showpoint<<setprecision(8)<<rn_b_z0[i]->thisRnode->T<<endl;
    }
    for(int i=0;i<n_boundry-1;i++){
      outfile<<(nx+2*nr-1)*(nz+1)+i<<", "<<rn_b_zn[i]->thisRnode->X/L<<", "<<rn_b_zn[i]->thisRnode->Z/D<<", "<<showpoint<<setprecision(8)<<rn_b_zn[i]->thisRnode->T<<endl;
    }
	outfile.close();
}


double simpson_qx_z(double *qx,int nz,double deltaz){
  return simpson(qx,nz,deltaz);
}
double simpson_deltaT(double *T1,double *T0,int nz,double deltaz){
  return simpson(T1,nz,deltaz)-simpson(T0,nz,deltaz);
}


void initializeRealNodes_c(realNode *rn[][nz+1],realHlist *Hrn[nx+2*nr-1][nz+1],realNode_boundry *rn_b_z0[],realNode_boundry *rn_b_zn[],int col[][nz+1]){
  realHlist *firstnode;
  realHlist *lastnode;
  realHlist *thisnode;

  for(int l=0;l<nx+2*nr-1;l++){
  	for(int m=0;m<nz+1;m++){
      rn[l][m]->HNnode=0;

      if((col[l][m]==1)||(col[l][m]==2)){
        for(int j=0;j<2*nr+1;j++){
          for(int k=0;k<2*nr+1;k++){
            if(j==nr&&k==nr){
            }else if((abs(j-nr)*abs(j-nr)+abs(k-nr)*abs(k-nr)<=nr*nr)){
              if(l-nr+j>=0&&l-nr+j<=nx+2*nr-2){
                int za=0;
                if(m-nr+k>=0&&m-nr+k<=nz){
                  za=(m-nr+k);
                }else if(m-nr+k<0){
                  za=(m-nr+k)+nz;
                }else{
                  za=(m-nr+k)-(nz);
                }
                if((col[l-nr+j][za]==1)||(col[l-nr+j][za]==2)){
                  if(rn[l][m]->HNnode==0){
                    firstnode=new realHlist;
                    firstnode->Xitempff=rn[l-nr+j][za]->tempff;
                    firstnode->Xitempbf=rn[l-nr+j][za]->tempbf;
                
                    firstnode->cXiX=(j-nr)*deltax;
                    firstnode->cXiZ=(k-nr)*deltaz;

                    firstnode->XiX=(j-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                    firstnode->XiZ=(k-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                    lastnode=firstnode;
                  }else{
                    thisnode=new realHlist;

                    thisnode->Xitempff=rn[l-nr+j][za]->tempff;
                    thisnode->Xitempbf=rn[l-nr+j][za]->tempbf;

                    thisnode->cXiX=(j-nr)*deltax;
                    thisnode->cXiZ=(k-nr)*deltaz;

                    thisnode->XiX=(j-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                    thisnode->XiZ=(k-nr)*areaSpace2/(((j-nr)*(j-nr)+(k-nr)*(k-nr)));
                    lastnode->next=thisnode;
                    lastnode=thisnode;
                  }
                  rn[l][m]->HNnode++;
                }
              }
            }
          }
        }

        if(col[l][m]==2){
          for(int j=0;j<n_boundry+1;j++){
            if((rn_b_z0[j]->thisRnode->X-rn[l][m]->X)*(rn_b_z0[j]->thisRnode->X-rn[l][m]->X)+(rn_b_z0[j]->thisRnode->Z-rn[l][m]->Z)*(rn_b_z0[j]->thisRnode->Z-rn[l][m]->Z)<=deltaxh*deltaxh*(1+1e-12)){
              if(rn[l][m]->HNnode==0){
                firstnode=new realHlist;
                firstnode->Xitempff=rn_b_z0[j]->thisRnode->tempff;
                firstnode->Xitempbf=rn_b_z0[j]->thisRnode->tempbf;
                
                firstnode->cXiX=rn_b_z0[j]->thisRnode->X-rn[l][m]->X;
                firstnode->cXiZ=rn_b_z0[j]->thisRnode->Z-rn[l][m]->Z;

                firstnode->XiX=rn_b_z0[j]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                firstnode->XiZ=rn_b_z0[j]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                lastnode=firstnode;
              }else{
                thisnode=new realHlist;

                thisnode->Xitempff=rn_b_z0[j]->thisRnode->tempff;
                thisnode->Xitempbf=rn_b_z0[j]->thisRnode->tempbf;

                thisnode->cXiX=rn_b_z0[j]->thisRnode->X-rn[l][m]->X;
                thisnode->cXiZ=rn_b_z0[j]->thisRnode->Z-rn[l][m]->Z;

                thisnode->XiX=rn_b_z0[j]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                thisnode->XiZ=rn_b_z0[j]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                lastnode->next=thisnode;
                lastnode=thisnode;
              }
              rn[l][m]->HNnode++;
            }
          }

          for(int j=0;j<n_boundry-1;j++){
            if((rn_b_zn[j]->thisRnode->X-rn[l][m]->X)*(rn_b_zn[j]->thisRnode->X-rn[l][m]->X)+(rn_b_zn[j]->thisRnode->Z-rn[l][m]->Z)*(rn_b_zn[j]->thisRnode->Z-rn[l][m]->Z)<=deltaxh*deltaxh*(1+1e-12)){
              if(rn[l][m]->HNnode==0){
                firstnode=new realHlist;
                firstnode->Xitempff=rn_b_zn[j]->thisRnode->tempff;
                firstnode->Xitempbf=rn_b_zn[j]->thisRnode->tempbf;
                
                firstnode->cXiX=rn_b_zn[j]->thisRnode->X-rn[l][m]->X;
                firstnode->cXiZ=rn_b_zn[j]->thisRnode->Z-rn[l][m]->Z;

                firstnode->XiX=rn_b_zn[j]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                firstnode->XiZ=rn_b_zn[j]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                lastnode=firstnode;
              }else{
                thisnode=new realHlist;

                thisnode->Xitempff=rn_b_zn[j]->thisRnode->tempff;
                thisnode->Xitempbf=rn_b_zn[j]->thisRnode->tempbf;

                thisnode->cXiX=rn_b_zn[j]->thisRnode->X-rn[l][m]->X;
                thisnode->cXiZ=rn_b_zn[j]->thisRnode->Z-rn[l][m]->Z;

                thisnode->XiX=rn_b_zn[j]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                thisnode->XiZ=rn_b_zn[j]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                lastnode->next=thisnode;
                lastnode=thisnode;
              }
              rn[l][m]->HNnode++;
            }
          }
        }

        rn[l][m]->HrnXiX=new double[rn[l][m]->HNnode];
        rn[l][m]->HrnXiZ=new double[rn[l][m]->HNnode];
        rn[l][m]->cXiX=new double[rn[l][m]->HNnode];
        rn[l][m]->cXiZ=new double[rn[l][m]->HNnode];

        rn[l][m]->Xitempff=new double**[rn[l][m]->HNnode];
        rn[l][m]->Xitempbf=new double**[rn[l][m]->HNnode];

        thisnode=firstnode;
        for(int ii=0;ii<rn[l][m]->HNnode;ii++){
          rn[l][m]->HrnXiX[ii]=thisnode->XiX;rn[l][m]->HrnXiZ[ii]=thisnode->XiZ;
          rn[l][m]->cXiX[ii]=thisnode->cXiX;rn[l][m]->cXiZ[ii]=thisnode->cXiZ;

          rn[l][m]->Xitempff[ii]=thisnode->Xitempff;
          rn[l][m]->Xitempbf[ii]=thisnode->Xitempbf;

          lastnode=thisnode;
          thisnode=lastnode->next;

          delete lastnode;
        }
      }
  	}
  }

  for(int l=0;l<n_boundry+1;l++){
    rn_b_z0[l]->thisRnode->HNnode=0;

    for(int j=xl1;j<=xl2;j++){
      for(int k=zl1;k<=zl2;k++){
        if(col[j][k]==2){
          if((rn[j][k]->X-rn_b_z0[l]->thisRnode->X)*(rn[j][k]->X-rn_b_z0[l]->thisRnode->X)+(rn[j][k]->Z-rn_b_z0[l]->thisRnode->Z)*(rn[j][k]->Z-rn_b_z0[l]->thisRnode->Z)<=deltaxh*deltaxh*(1+1e-12)){
            if(rn_b_z0[l]->thisRnode->HNnode==0){
                firstnode=new realHlist;
                firstnode->Xitempff=rn[j][k]->tempff;
                firstnode->Xitempbf=rn[j][k]->tempbf;
                
                firstnode->cXiX=rn[j][k]->X-rn_b_z0[l]->thisRnode->X;
                firstnode->cXiZ=rn[j][k]->Z-rn_b_z0[l]->thisRnode->Z;

                firstnode->XiX=rn[j][k]->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                firstnode->XiZ=rn[j][k]->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                lastnode=firstnode;
              }else{
                thisnode=new realHlist;

                thisnode->Xitempff=rn[j][k]->tempff;
                thisnode->Xitempbf=rn[j][k]->tempbf;

                thisnode->cXiX=rn[j][k]->X-rn_b_z0[l]->thisRnode->X;
                thisnode->cXiZ=rn[j][k]->Z-rn_b_z0[l]->thisRnode->Z;

                thisnode->XiX=rn[j][k]->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                thisnode->XiZ=rn[j][k]->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                lastnode->next=thisnode;
                lastnode=thisnode;
              }
              rn_b_z0[l]->thisRnode->HNnode++;


          }
        }
      }
    }

    if(l==0){
      if(rn_b_z0[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_z0[l+1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_z0[l+1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_z0[l+1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_z0[l+1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        firstnode->XiX=rn_b_z0[l+1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_z0[l+1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_z0[l+1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_z0[l+1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_z0[l+1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_z0[l+1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        thisnode->XiX=rn_b_z0[l+1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_z0[l+1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_z0[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_zn[0]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_zn[0]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_zn[0]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_zn[0]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

      thisnode->XiX=rn_b_zn[0]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_zn[0]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_z0[l]->thisRnode->HNnode++;

    }else if(l==n_boundry){
      if(rn_b_z0[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_z0[l-1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_z0[l-1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_z0[l-1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_z0[l-1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        firstnode->XiX=rn_b_z0[l-1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_z0[l-1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_z0[l-1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_z0[l-1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_z0[l-1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_z0[l-1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        thisnode->XiX=rn_b_z0[l-1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_z0[l-1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_z0[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_zn[n_boundry-2]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_zn[n_boundry-2]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_zn[n_boundry-2]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_zn[n_boundry-2]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

      thisnode->XiX=rn_b_zn[n_boundry-2]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_zn[n_boundry-2]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_z0[l]->thisRnode->HNnode++;

    }else{
      if(rn_b_z0[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_z0[l+1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_z0[l+1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_z0[l+1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_z0[l+1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        firstnode->XiX=rn_b_z0[l+1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_z0[l+1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_z0[l+1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_z0[l+1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_z0[l+1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_z0[l+1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

        thisnode->XiX=rn_b_z0[l+1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_z0[l+1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_z0[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_z0[l-1]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_z0[l-1]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_z0[l-1]->thisRnode->X-rn_b_z0[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_z0[l-1]->thisRnode->Z-rn_b_z0[l]->thisRnode->Z;

      thisnode->XiX=rn_b_z0[l-1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_z0[l-1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_z0[l]->thisRnode->HNnode++;
    }

    rn_b_z0[l]->thisRnode->HrnXiX=new double[rn_b_z0[l]->thisRnode->HNnode];
    rn_b_z0[l]->thisRnode->HrnXiZ=new double[rn_b_z0[l]->thisRnode->HNnode];
    rn_b_z0[l]->thisRnode->cXiX=new double[rn_b_z0[l]->thisRnode->HNnode];
    rn_b_z0[l]->thisRnode->cXiZ=new double[rn_b_z0[l]->thisRnode->HNnode];

    rn_b_z0[l]->thisRnode->Xitempff=new double**[rn_b_z0[l]->thisRnode->HNnode];
    rn_b_z0[l]->thisRnode->Xitempbf=new double**[rn_b_z0[l]->thisRnode->HNnode];

    thisnode=firstnode;
    for(int ii=0;ii<rn_b_z0[l]->thisRnode->HNnode;ii++){
      rn_b_z0[l]->thisRnode->HrnXiX[ii]=thisnode->XiX;rn_b_z0[l]->thisRnode->HrnXiZ[ii]=thisnode->XiZ;
      rn_b_z0[l]->thisRnode->cXiX[ii]=thisnode->cXiX;rn_b_z0[l]->thisRnode->cXiZ[ii]=thisnode->cXiZ;

      rn_b_z0[l]->thisRnode->Xitempff[ii]=thisnode->Xitempff;
      rn_b_z0[l]->thisRnode->Xitempbf[ii]=thisnode->Xitempbf;

      lastnode=thisnode;
      thisnode=lastnode->next;

      delete lastnode;
    }

  }

  for(int l=0;l<n_boundry-1;l++){
    rn_b_zn[l]->thisRnode->HNnode=0;

    for(int j=xl1;j<=xl2;j++){
      for(int k=zl1;k<=zl2;k++){
        if(col[j][k]==2){
          if((rn[j][k]->X-rn_b_zn[l]->thisRnode->X)*(rn[j][k]->X-rn_b_zn[l]->thisRnode->X)+(rn[j][k]->Z-rn_b_zn[l]->thisRnode->Z)*(rn[j][k]->Z-rn_b_zn[l]->thisRnode->Z)<=deltaxh*deltaxh*(1+1e-12)){
            if(rn_b_zn[l]->thisRnode->HNnode==0){
                firstnode=new realHlist;
                firstnode->Xitempff=rn[j][k]->tempff;
                firstnode->Xitempbf=rn[j][k]->tempbf;
                
                firstnode->cXiX=rn[j][k]->X-rn_b_zn[l]->thisRnode->X;
                firstnode->cXiZ=rn[j][k]->Z-rn_b_zn[l]->thisRnode->Z;

                firstnode->XiX=rn[j][k]->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                firstnode->XiZ=rn[j][k]->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
                lastnode=firstnode;
              }else{
                thisnode=new realHlist;

                thisnode->Xitempff=rn[j][k]->tempff;
                thisnode->Xitempbf=rn[j][k]->tempbf;

                thisnode->cXiX=rn[j][k]->X-rn_b_zn[l]->thisRnode->X;
                thisnode->cXiZ=rn[j][k]->Z-rn_b_zn[l]->thisRnode->Z;

                thisnode->XiX=rn[j][k]->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                thisnode->XiZ=rn[j][k]->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
                lastnode->next=thisnode;
                lastnode=thisnode;
              }
              rn_b_zn[l]->thisRnode->HNnode++;

          }
        }
      }
    }

    if(l==0){
      if(rn_b_zn[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_zn[l+1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_zn[l+1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_zn[l+1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_zn[l+1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        firstnode->XiX=rn_b_zn[l+1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_zn[l+1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_zn[l+1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_zn[l+1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_zn[l+1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_zn[l+1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        thisnode->XiX=rn_b_zn[l+1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_zn[l+1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_zn[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_z0[0]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_z0[0]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_z0[0]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_z0[0]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

      thisnode->XiX=rn_b_z0[0]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_z0[0]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_zn[l]->thisRnode->HNnode++;

    }else if(l==n_boundry-2){
      if(rn_b_zn[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_zn[l-1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_zn[l-1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_zn[l-1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_zn[l-1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        firstnode->XiX=rn_b_zn[l-1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_zn[l-1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_zn[l-1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_zn[l-1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_zn[l-1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_zn[l-1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        thisnode->XiX=rn_b_zn[l-1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_zn[l-1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_zn[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_z0[n_boundry]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_z0[n_boundry]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_z0[n_boundry]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_z0[n_boundry]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

      thisnode->XiX=rn_b_z0[n_boundry]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_z0[n_boundry]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_zn[l]->thisRnode->HNnode++;

    }else{
      if(rn_b_zn[l]->thisRnode->HNnode==0){
        firstnode=new realHlist;
        firstnode->Xitempff=rn_b_zn[l+1]->thisRnode->tempff;
        firstnode->Xitempbf=rn_b_zn[l+1]->thisRnode->tempbf;
                
        firstnode->cXiX=rn_b_zn[l+1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        firstnode->cXiZ=rn_b_zn[l+1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        firstnode->XiX=rn_b_zn[l+1]->thisRnode->area*firstnode->cXiX/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        firstnode->XiZ=rn_b_zn[l+1]->thisRnode->area*firstnode->cXiZ/((firstnode->cXiX)*(firstnode->cXiX)+(firstnode->cXiZ)*(firstnode->cXiZ));
        lastnode=firstnode;
      }else{
        thisnode=new realHlist;

        thisnode->Xitempff=rn_b_zn[l+1]->thisRnode->tempff;
        thisnode->Xitempbf=rn_b_zn[l+1]->thisRnode->tempbf;

        thisnode->cXiX=rn_b_zn[l+1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
        thisnode->cXiZ=rn_b_zn[l+1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

        thisnode->XiX=rn_b_zn[l+1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        thisnode->XiZ=rn_b_zn[l+1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
        lastnode->next=thisnode;
        lastnode=thisnode;
      }
      rn_b_zn[l]->thisRnode->HNnode++;

      thisnode=new realHlist;

      thisnode->Xitempff=rn_b_zn[l-1]->thisRnode->tempff;
      thisnode->Xitempbf=rn_b_zn[l-1]->thisRnode->tempbf;

      thisnode->cXiX=rn_b_zn[l-1]->thisRnode->X-rn_b_zn[l]->thisRnode->X;
      thisnode->cXiZ=rn_b_zn[l-1]->thisRnode->Z-rn_b_zn[l]->thisRnode->Z;

      thisnode->XiX=rn_b_zn[l-1]->thisRnode->area*thisnode->cXiX/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      thisnode->XiZ=rn_b_zn[l-1]->thisRnode->area*thisnode->cXiZ/((thisnode->cXiX)*(thisnode->cXiX)+(thisnode->cXiZ)*(thisnode->cXiZ));
      lastnode->next=thisnode;
      lastnode=thisnode;
      rn_b_zn[l]->thisRnode->HNnode++;
    }

    rn_b_zn[l]->thisRnode->HrnXiX=new double[rn_b_zn[l]->thisRnode->HNnode];
    rn_b_zn[l]->thisRnode->HrnXiZ=new double[rn_b_zn[l]->thisRnode->HNnode];
    rn_b_zn[l]->thisRnode->cXiX=new double[rn_b_zn[l]->thisRnode->HNnode];
    rn_b_zn[l]->thisRnode->cXiZ=new double[rn_b_zn[l]->thisRnode->HNnode];

    rn_b_zn[l]->thisRnode->Xitempff=new double**[rn_b_zn[l]->thisRnode->HNnode];
    rn_b_zn[l]->thisRnode->Xitempbf=new double**[rn_b_zn[l]->thisRnode->HNnode];

    thisnode=firstnode;
    for(int ii=0;ii<rn_b_zn[l]->thisRnode->HNnode;ii++){
      rn_b_zn[l]->thisRnode->HrnXiX[ii]=thisnode->XiX;rn_b_zn[l]->thisRnode->HrnXiZ[ii]=thisnode->XiZ;
      rn_b_zn[l]->thisRnode->cXiX[ii]=thisnode->cXiX;rn_b_zn[l]->thisRnode->cXiZ[ii]=thisnode->cXiZ;

      rn_b_zn[l]->thisRnode->Xitempff[ii]=thisnode->Xitempff;
      rn_b_zn[l]->thisRnode->Xitempbf[ii]=thisnode->Xitempbf;

      lastnode=thisnode;
      thisnode=lastnode->next;

      delete lastnode;
    }
  }
}

void makeparaff_c(realNode *rn[][nz+1],double *sintheta1,double *sintheta2,double *costheta1,double *costheta2,double *sinu1,double *sinu2,realNode_boundry *rn_b_z0[],realNode_boundry *rn_b_zn[]){
  int nrr;
  double sum1=0;
  double sum2=0;
  double xi=0;
  for(int l=0;l<nx+2*nr-1;l++){
  	for(int m=0;m<nz+1;m++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn[l][m]->paraS2ff[n][i]=0;
          rn[l][m]->paraC2ff[n][i]=0;
          rn[l][m]->paraSCff[n][i]=0;

          rn[l][m]->para1ff[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn[l][m]->HNnode;ii++){
            if((-1.*rn[l][m]->cXiX[ii])*costheta1[i]+(-1.*rn[l][m]->cXiZ[ii])*sintheta1[i]*sinu1[n]>0.){
              sum1+=rn[l][m]->HrnXiX[ii];
              sum2+=rn[l][m]->HrnXiZ[ii];

              rn[l][m]->paraS2ff[n][i]+=rn[l][m]->cXiZ[ii]*rn[l][m]->HrnXiZ[ii];
              rn[l][m]->paraC2ff[n][i]+=rn[l][m]->cXiX[ii]*rn[l][m]->HrnXiX[ii];
              rn[l][m]->paraSCff[n][i]+=rn[l][m]->cXiZ[ii]*rn[l][m]->HrnXiX[ii];
            }
          }

          if(rn[l][m]->paraSCff[n][i]==0){
            if(rn[l][m]->paraC2ff[n][i]==0){
              if(rn[l][m]->paraS2ff[n][i]!=0){
                rn[l][m]->para21ff[n][i]=0.;
                rn[l][m]->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn[l][m]->paraS2ff[n][i];
                rn[l][m]->para1ff[n][i]=sum2*rn[l][m]->para22ff[n][i];
              }
              
            }else if(rn[l][m]->paraS2ff[n][i]==0){
              rn[l][m]->para21ff[n][i]=costheta1[i]/rn[l][m]->paraC2ff[n][i];
              rn[l][m]->para22ff[n][i]=0.;
              rn[l][m]->para1ff[n][i]=sum1*rn[l][m]->para21ff[n][i];
            }else{
              rn[l][m]->para21ff[n][i]=costheta1[i]/rn[l][m]->paraC2ff[n][i];
              rn[l][m]->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn[l][m]->paraS2ff[n][i];
              rn[l][m]->para1ff[n][i]=sum1*rn[l][m]->para21ff[n][i]+sum2*rn[l][m]->para22ff[n][i];
            }
          }else{
            rn[l][m]->para21ff[n][i]=((costheta1[i]*(rn[l][m]->paraS2ff[n][i]/rn[l][m]->paraSCff[n][i])/(rn[l][m]->paraC2ff[n][i]*(rn[l][m]->paraS2ff[n][i]/rn[l][m]->paraSCff[n][i])-rn[l][m]->paraSCff[n][i]))-(sintheta1[i]*sinu1[n]/(rn[l][m]->paraS2ff[n][i]*(rn[l][m]->paraC2ff[n][i]/rn[l][m]->paraSCff[n][i])-rn[l][m]->paraSCff[n][i])));
            rn[l][m]->para22ff[n][i]=((sintheta1[i]*sinu1[n]*(rn[l][m]->paraC2ff[n][i]/rn[l][m]->paraSCff[n][i])/(rn[l][m]->paraS2ff[n][i]*(rn[l][m]->paraC2ff[n][i]/rn[l][m]->paraSCff[n][i])-rn[l][m]->paraSCff[n][i]))-(costheta1[i]/(rn[l][m]->paraC2ff[n][i]*(rn[l][m]->paraS2ff[n][i]/rn[l][m]->paraSCff[n][i])-rn[l][m]->paraSCff[n][i])));
            rn[l][m]->para1ff[n][i]=sum1*rn[l][m]->para21ff[n][i]+sum2*rn[l][m]->para22ff[n][i];
          }

          rn[l][m]->para1ff[n][i]+=-1./llsi;

          rn[l][m]->para1ff[n][i]=1./rn[l][m]->para1ff[n][i];
          
        }
      }
  	}
  } 

  for(int l=0;l<n_boundry+1;l++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn_b_z0[l]->thisRnode->paraS2ff[n][i]=0;
          rn_b_z0[l]->thisRnode->paraC2ff[n][i]=0;
          rn_b_z0[l]->thisRnode->paraSCff[n][i]=0;

          rn_b_z0[l]->thisRnode->para1ff[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn_b_z0[l]->thisRnode->HNnode;ii++){
            if((-1.*rn_b_z0[l]->thisRnode->cXiX[ii])*costheta1[i]+(-1.*rn_b_z0[l]->thisRnode->cXiZ[ii])*sintheta1[i]*sinu1[n]>0.){
              sum1+=rn_b_z0[l]->thisRnode->HrnXiX[ii];
              sum2+=rn_b_z0[l]->thisRnode->HrnXiZ[ii];

              rn_b_z0[l]->thisRnode->paraS2ff[n][i]+=rn_b_z0[l]->thisRnode->cXiZ[ii]*rn_b_z0[l]->thisRnode->HrnXiZ[ii];
              rn_b_z0[l]->thisRnode->paraC2ff[n][i]+=rn_b_z0[l]->thisRnode->cXiX[ii]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
              rn_b_z0[l]->thisRnode->paraSCff[n][i]+=rn_b_z0[l]->thisRnode->cXiZ[ii]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
            }
          }

          if(rn_b_z0[l]->thisRnode->paraSCff[n][i]==0){
            if(rn_b_z0[l]->thisRnode->paraC2ff[n][i]==0){
              if(rn_b_z0[l]->thisRnode->paraS2ff[n][i]!=0){
                rn_b_z0[l]->thisRnode->para21ff[n][i]=0.;
                rn_b_z0[l]->thisRnode->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn_b_z0[l]->thisRnode->paraS2ff[n][i];
                rn_b_z0[l]->thisRnode->para1ff[n][i]=sum2*rn_b_z0[l]->thisRnode->para22ff[n][i];
              }
              
            }else if(rn_b_z0[l]->thisRnode->paraS2ff[n][i]==0){
              rn_b_z0[l]->thisRnode->para21ff[n][i]=costheta1[i]/rn_b_z0[l]->thisRnode->paraC2ff[n][i];
              rn_b_z0[l]->thisRnode->para22ff[n][i]=0.;
              rn_b_z0[l]->thisRnode->para1ff[n][i]=sum1*rn_b_z0[l]->thisRnode->para21ff[n][i];
            }else{
              rn_b_z0[l]->thisRnode->para21ff[n][i]=costheta1[i]/rn_b_z0[l]->thisRnode->paraC2ff[n][i];
              rn_b_z0[l]->thisRnode->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn_b_z0[l]->thisRnode->paraS2ff[n][i];
              rn_b_z0[l]->thisRnode->para1ff[n][i]=sum1*rn_b_z0[l]->thisRnode->para21ff[n][i]+sum2*rn_b_z0[l]->thisRnode->para22ff[n][i];
            }
          }else{
            rn_b_z0[l]->thisRnode->para21ff[n][i]=((costheta1[i]*(rn_b_z0[l]->thisRnode->paraS2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])/(rn_b_z0[l]->thisRnode->paraC2ff[n][i]*(rn_b_z0[l]->thisRnode->paraS2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])-rn_b_z0[l]->thisRnode->paraSCff[n][i]))-(sintheta1[i]*sinu1[n]/(rn_b_z0[l]->thisRnode->paraS2ff[n][i]*(rn_b_z0[l]->thisRnode->paraC2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])-rn_b_z0[l]->thisRnode->paraSCff[n][i])));
            rn_b_z0[l]->thisRnode->para22ff[n][i]=((sintheta1[i]*sinu1[n]*(rn_b_z0[l]->thisRnode->paraC2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])/(rn_b_z0[l]->thisRnode->paraS2ff[n][i]*(rn_b_z0[l]->thisRnode->paraC2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])-rn_b_z0[l]->thisRnode->paraSCff[n][i]))-(costheta1[i]/(rn_b_z0[l]->thisRnode->paraC2ff[n][i]*(rn_b_z0[l]->thisRnode->paraS2ff[n][i]/rn_b_z0[l]->thisRnode->paraSCff[n][i])-rn_b_z0[l]->thisRnode->paraSCff[n][i])));
            rn_b_z0[l]->thisRnode->para1ff[n][i]=sum1*rn_b_z0[l]->thisRnode->para21ff[n][i]+sum2*rn_b_z0[l]->thisRnode->para22ff[n][i];
          }

          rn_b_z0[l]->thisRnode->para1ff[n][i]+=-1./llsi;

          rn_b_z0[l]->thisRnode->para1ff[n][i]=1./rn_b_z0[l]->thisRnode->para1ff[n][i];
          
        }
      }
  }

  for(int l=0;l<n_boundry-1;l++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn_b_zn[l]->thisRnode->paraS2ff[n][i]=0;
          rn_b_zn[l]->thisRnode->paraC2ff[n][i]=0;
          rn_b_zn[l]->thisRnode->paraSCff[n][i]=0;

          rn_b_zn[l]->thisRnode->para1ff[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn_b_zn[l]->thisRnode->HNnode;ii++){
            if((-1.*rn_b_zn[l]->thisRnode->cXiX[ii])*costheta1[i]+(-1.*rn_b_zn[l]->thisRnode->cXiZ[ii])*sintheta1[i]*sinu1[n]>0.){
              sum1+=rn_b_zn[l]->thisRnode->HrnXiX[ii];
              sum2+=rn_b_zn[l]->thisRnode->HrnXiZ[ii];

              rn_b_zn[l]->thisRnode->paraS2ff[n][i]+=rn_b_zn[l]->thisRnode->cXiZ[ii]*rn_b_zn[l]->thisRnode->HrnXiZ[ii];
              rn_b_zn[l]->thisRnode->paraC2ff[n][i]+=rn_b_zn[l]->thisRnode->cXiX[ii]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
              rn_b_zn[l]->thisRnode->paraSCff[n][i]+=rn_b_zn[l]->thisRnode->cXiZ[ii]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
            }
          }

          if(rn_b_zn[l]->thisRnode->paraSCff[n][i]==0){
            if(rn_b_zn[l]->thisRnode->paraC2ff[n][i]==0){
              if(rn_b_zn[l]->thisRnode->paraS2ff[n][i]!=0){
                rn_b_zn[l]->thisRnode->para21ff[n][i]=0.;
                rn_b_zn[l]->thisRnode->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn_b_zn[l]->thisRnode->paraS2ff[n][i];
                rn_b_zn[l]->thisRnode->para1ff[n][i]=sum2*rn_b_zn[l]->thisRnode->para22ff[n][i];
              }
              
            }else if(rn_b_zn[l]->thisRnode->paraS2ff[n][i]==0){
              rn_b_zn[l]->thisRnode->para21ff[n][i]=costheta1[i]/rn_b_zn[l]->thisRnode->paraC2ff[n][i];
              rn_b_zn[l]->thisRnode->para22ff[n][i]=0.;
              rn_b_zn[l]->thisRnode->para1ff[n][i]=sum1*rn_b_zn[l]->thisRnode->para21ff[n][i];
            }else{
              rn_b_zn[l]->thisRnode->para21ff[n][i]=costheta1[i]/rn_b_zn[l]->thisRnode->paraC2ff[n][i];
              rn_b_zn[l]->thisRnode->para22ff[n][i]=sintheta1[i]*sinu1[n]/rn_b_zn[l]->thisRnode->paraS2ff[n][i];
              rn_b_zn[l]->thisRnode->para1ff[n][i]=sum1*rn_b_zn[l]->thisRnode->para21ff[n][i]+sum2*rn_b_zn[l]->thisRnode->para22ff[n][i];
            }
          }else{
            rn_b_zn[l]->thisRnode->para21ff[n][i]=((costheta1[i]*(rn_b_zn[l]->thisRnode->paraS2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])/(rn_b_zn[l]->thisRnode->paraC2ff[n][i]*(rn_b_zn[l]->thisRnode->paraS2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])-rn_b_zn[l]->thisRnode->paraSCff[n][i]))-(sintheta1[i]*sinu1[n]/(rn_b_zn[l]->thisRnode->paraS2ff[n][i]*(rn_b_zn[l]->thisRnode->paraC2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])-rn_b_zn[l]->thisRnode->paraSCff[n][i])));
            rn_b_zn[l]->thisRnode->para22ff[n][i]=((sintheta1[i]*sinu1[n]*(rn_b_zn[l]->thisRnode->paraC2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])/(rn_b_zn[l]->thisRnode->paraS2ff[n][i]*(rn_b_zn[l]->thisRnode->paraC2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])-rn_b_zn[l]->thisRnode->paraSCff[n][i]))-(costheta1[i]/(rn_b_zn[l]->thisRnode->paraC2ff[n][i]*(rn_b_zn[l]->thisRnode->paraS2ff[n][i]/rn_b_zn[l]->thisRnode->paraSCff[n][i])-rn_b_zn[l]->thisRnode->paraSCff[n][i])));
            rn_b_zn[l]->thisRnode->para1ff[n][i]=sum1*rn_b_zn[l]->thisRnode->para21ff[n][i]+sum2*rn_b_zn[l]->thisRnode->para22ff[n][i];
          }

          rn_b_zn[l]->thisRnode->para1ff[n][i]+=-1./llsi;

          rn_b_zn[l]->thisRnode->para1ff[n][i]=1./rn_b_zn[l]->thisRnode->para1ff[n][i];
          
        }
      }
  }

}

void makeparabf_c(realNode *rn[][nz+1],double *sintheta1,double *sintheta2,double *costheta1,double *costheta2,double *sinu1,double *sinu2,realNode_boundry *rn_b_z0[],realNode_boundry *rn_b_zn[]){
  int nrr;

  double sum1=0;
  double sum2=0;
  double xi=0;
  for(int l=0;l<nx+2*nr-1;l++){
  	for(int m=0;m<nz+1;m++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn[l][m]->paraS2bf[n][i]=0;
          rn[l][m]->paraC2bf[n][i]=0;
          rn[l][m]->paraSCbf[n][i]=0;

          rn[l][m]->para1bf[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn[l][m]->HNnode;ii++){
            if((-1.*rn[l][m]->cXiX[ii])*costheta2[i]+(-1.*rn[l][m]->cXiZ[ii])*sintheta2[i]*sinu1[n]>0.){
              sum1+=rn[l][m]->HrnXiX[ii];
              sum2+=rn[l][m]->HrnXiZ[ii];

              rn[l][m]->paraS2bf[n][i]+=rn[l][m]->cXiZ[ii]*rn[l][m]->HrnXiZ[ii];
              rn[l][m]->paraC2bf[n][i]+=rn[l][m]->cXiX[ii]*rn[l][m]->HrnXiX[ii];
              rn[l][m]->paraSCbf[n][i]+=rn[l][m]->cXiZ[ii]*rn[l][m]->HrnXiX[ii];
            }

          }

          if(rn[l][m]->paraSCbf[n][i]==0){
            if(rn[l][m]->paraC2bf[n][i]==0){
              if(rn[l][m]->paraS2bf[n][i]!=0){
                rn[l][m]->para21bf[n][i]=0.;
                rn[l][m]->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn[l][m]->paraS2bf[n][i];
                rn[l][m]->para1bf[n][i]=rn[l][m]->para22bf[n][i]*sum2;
              } 
            }else if(rn[l][m]->paraS2bf[n][i]==0){
              rn[l][m]->para21bf[n][i]=costheta2[i]/rn[l][m]->paraC2bf[n][i];
              rn[l][m]->para22bf[n][i]=0.;
              rn[l][m]->para1bf[n][i]=rn[l][m]->para21bf[n][i]*sum1;
            }else{
              rn[l][m]->para21bf[n][i]=costheta2[i]/rn[l][m]->paraC2bf[n][i];
              rn[l][m]->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn[l][m]->paraS2bf[n][i];
              rn[l][m]->para1bf[n][i]=sum1*rn[l][m]->para21bf[n][i]+sum2*rn[l][m]->para22bf[n][i];
            }
          }else{
            rn[l][m]->para21bf[n][i]=((costheta2[i]*(rn[l][m]->paraS2bf[n][i]/rn[l][m]->paraSCbf[n][i])/(rn[l][m]->paraC2bf[n][i]*(rn[l][m]->paraS2bf[n][i]/rn[l][m]->paraSCbf[n][i])-rn[l][m]->paraSCbf[n][i]))-(sintheta2[i]*sinu1[n]/(rn[l][m]->paraS2bf[n][i]*(rn[l][m]->paraC2bf[n][i]/rn[l][m]->paraSCbf[n][i])-rn[l][m]->paraSCbf[n][i])));
            rn[l][m]->para22bf[n][i]=((sintheta2[i]*sinu1[n]*(rn[l][m]->paraC2bf[n][i]/rn[l][m]->paraSCbf[n][i])/(rn[l][m]->paraS2bf[n][i]*(rn[l][m]->paraC2bf[n][i]/rn[l][m]->paraSCbf[n][i])-rn[l][m]->paraSCbf[n][i]))-(costheta2[i]/(rn[l][m]->paraC2bf[n][i]*(rn[l][m]->paraS2bf[n][i]/rn[l][m]->paraSCbf[n][i])-rn[l][m]->paraSCbf[n][i])));
            rn[l][m]->para1bf[n][i]=sum1*rn[l][m]->para21bf[n][i]+sum2*rn[l][m]->para22bf[n][i];
          }

          rn[l][m]->para1bf[n][i]+=-1./llsi;

          rn[l][m]->para1bf[n][i]=1./rn[l][m]->para1bf[n][i];
        }
      }
  	}
  } 

  for(int l=0;l<n_boundry+1;l++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn_b_z0[l]->thisRnode->paraS2bf[n][i]=0;
          rn_b_z0[l]->thisRnode->paraC2bf[n][i]=0;
          rn_b_z0[l]->thisRnode->paraSCbf[n][i]=0;

          rn_b_z0[l]->thisRnode->para1bf[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn_b_z0[l]->thisRnode->HNnode;ii++){
            if((-1.*rn_b_z0[l]->thisRnode->cXiX[ii])*costheta2[i]+(-1.*rn_b_z0[l]->thisRnode->cXiZ[ii])*sintheta2[i]*sinu1[n]>0.){
              sum1+=rn_b_z0[l]->thisRnode->HrnXiX[ii];
              sum2+=rn_b_z0[l]->thisRnode->HrnXiZ[ii];

              rn_b_z0[l]->thisRnode->paraS2bf[n][i]+=rn_b_z0[l]->thisRnode->cXiZ[ii]*rn_b_z0[l]->thisRnode->HrnXiZ[ii];
              rn_b_z0[l]->thisRnode->paraC2bf[n][i]+=rn_b_z0[l]->thisRnode->cXiX[ii]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
              rn_b_z0[l]->thisRnode->paraSCbf[n][i]+=rn_b_z0[l]->thisRnode->cXiZ[ii]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
            }

          }

          if(rn_b_z0[l]->thisRnode->paraSCbf[n][i]==0){
            if(rn_b_z0[l]->thisRnode->paraC2bf[n][i]==0){
              if(rn_b_z0[l]->thisRnode->paraS2bf[n][i]!=0){
                rn_b_z0[l]->thisRnode->para21bf[n][i]=0.;
                rn_b_z0[l]->thisRnode->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn_b_z0[l]->thisRnode->paraS2bf[n][i];
                rn_b_z0[l]->thisRnode->para1bf[n][i]=rn_b_z0[l]->thisRnode->para22bf[n][i]*sum2;
              } 
            }else if(rn_b_z0[l]->thisRnode->paraS2bf[n][i]==0){
              rn_b_z0[l]->thisRnode->para21bf[n][i]=costheta2[i]/rn_b_z0[l]->thisRnode->paraC2bf[n][i];
              rn_b_z0[l]->thisRnode->para22bf[n][i]=0.;
              rn_b_z0[l]->thisRnode->para1bf[n][i]=rn_b_z0[l]->thisRnode->para21bf[n][i]*sum1;
            }else{
              rn_b_z0[l]->thisRnode->para21bf[n][i]=costheta2[i]/rn_b_z0[l]->thisRnode->paraC2bf[n][i];
              rn_b_z0[l]->thisRnode->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn_b_z0[l]->thisRnode->paraS2bf[n][i];
              rn_b_z0[l]->thisRnode->para1bf[n][i]=sum1*rn_b_z0[l]->thisRnode->para21bf[n][i]+sum2*rn_b_z0[l]->thisRnode->para22bf[n][i];
            }
          }else{
            rn_b_z0[l]->thisRnode->para21bf[n][i]=((costheta2[i]*(rn_b_z0[l]->thisRnode->paraS2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])/(rn_b_z0[l]->thisRnode->paraC2bf[n][i]*(rn_b_z0[l]->thisRnode->paraS2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])-rn_b_z0[l]->thisRnode->paraSCbf[n][i]))-(sintheta2[i]*sinu1[n]/(rn_b_z0[l]->thisRnode->paraS2bf[n][i]*(rn_b_z0[l]->thisRnode->paraC2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])-rn_b_z0[l]->thisRnode->paraSCbf[n][i])));
            rn_b_z0[l]->thisRnode->para22bf[n][i]=((sintheta2[i]*sinu1[n]*(rn_b_z0[l]->thisRnode->paraC2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])/(rn_b_z0[l]->thisRnode->paraS2bf[n][i]*(rn_b_z0[l]->thisRnode->paraC2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])-rn_b_z0[l]->thisRnode->paraSCbf[n][i]))-(costheta2[i]/(rn_b_z0[l]->thisRnode->paraC2bf[n][i]*(rn_b_z0[l]->thisRnode->paraS2bf[n][i]/rn_b_z0[l]->thisRnode->paraSCbf[n][i])-rn_b_z0[l]->thisRnode->paraSCbf[n][i])));
            rn_b_z0[l]->thisRnode->para1bf[n][i]=sum1*rn_b_z0[l]->thisRnode->para21bf[n][i]+sum2*rn_b_z0[l]->thisRnode->para22bf[n][i];
          }

          rn_b_z0[l]->thisRnode->para1bf[n][i]+=-1./llsi;

          rn_b_z0[l]->thisRnode->para1bf[n][i]=1./rn_b_z0[l]->thisRnode->para1bf[n][i];
        }
      }
  }

  for(int l=0;l<n_boundry-1;l++){
      for(int n=0;n<glnn2;n++){
        for(int i=0;i<glnn1;i++){
          rn_b_zn[l]->thisRnode->paraS2bf[n][i]=0;
          rn_b_zn[l]->thisRnode->paraC2bf[n][i]=0;
          rn_b_zn[l]->thisRnode->paraSCbf[n][i]=0;

          rn_b_zn[l]->thisRnode->para1bf[n][i]=0;

          sum1=0;
          sum2=0;

          for(int ii=0;ii<rn_b_zn[l]->thisRnode->HNnode;ii++){
            if((-1.*rn_b_zn[l]->thisRnode->cXiX[ii])*costheta2[i]+(-1.*rn_b_zn[l]->thisRnode->cXiZ[ii])*sintheta2[i]*sinu1[n]>0.){
              sum1+=rn_b_zn[l]->thisRnode->HrnXiX[ii];
              sum2+=rn_b_zn[l]->thisRnode->HrnXiZ[ii];
              rn_b_zn[l]->thisRnode->paraS2bf[n][i]+=rn_b_zn[l]->thisRnode->cXiZ[ii]*rn_b_zn[l]->thisRnode->HrnXiZ[ii];
              rn_b_zn[l]->thisRnode->paraC2bf[n][i]+=rn_b_zn[l]->thisRnode->cXiX[ii]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
              rn_b_zn[l]->thisRnode->paraSCbf[n][i]+=rn_b_zn[l]->thisRnode->cXiZ[ii]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
            }

          }

          if(rn_b_zn[l]->thisRnode->paraSCbf[n][i]==0){
            if(rn_b_zn[l]->thisRnode->paraC2bf[n][i]==0){
              if(rn_b_zn[l]->thisRnode->paraS2bf[n][i]!=0){
                rn_b_zn[l]->thisRnode->para21bf[n][i]=0.;
                rn_b_zn[l]->thisRnode->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn_b_zn[l]->thisRnode->paraS2bf[n][i];
                rn_b_zn[l]->thisRnode->para1bf[n][i]=rn_b_zn[l]->thisRnode->para22bf[n][i]*sum2;
              } 
            }else if(rn_b_zn[l]->thisRnode->paraS2bf[n][i]==0){
              rn_b_zn[l]->thisRnode->para21bf[n][i]=costheta2[i]/rn_b_zn[l]->thisRnode->paraC2bf[n][i];
              rn_b_zn[l]->thisRnode->para22bf[n][i]=0.;
              rn_b_zn[l]->thisRnode->para1bf[n][i]=rn_b_zn[l]->thisRnode->para21bf[n][i]*sum1;
            }else{
              rn_b_zn[l]->thisRnode->para21bf[n][i]=costheta2[i]/rn_b_zn[l]->thisRnode->paraC2bf[n][i];
              rn_b_zn[l]->thisRnode->para22bf[n][i]=sintheta2[i]*sinu1[n]/rn_b_zn[l]->thisRnode->paraS2bf[n][i];
              rn_b_zn[l]->thisRnode->para1bf[n][i]=sum1*rn_b_zn[l]->thisRnode->para21bf[n][i]+sum2*rn_b_zn[l]->thisRnode->para22bf[n][i];
            }
          }else{
            rn_b_zn[l]->thisRnode->para21bf[n][i]=((costheta2[i]*(rn_b_zn[l]->thisRnode->paraS2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])/(rn_b_zn[l]->thisRnode->paraC2bf[n][i]*(rn_b_zn[l]->thisRnode->paraS2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])-rn_b_zn[l]->thisRnode->paraSCbf[n][i]))-(sintheta2[i]*sinu1[n]/(rn_b_zn[l]->thisRnode->paraS2bf[n][i]*(rn_b_zn[l]->thisRnode->paraC2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])-rn_b_zn[l]->thisRnode->paraSCbf[n][i])));
            rn_b_zn[l]->thisRnode->para22bf[n][i]=((sintheta2[i]*sinu1[n]*(rn_b_zn[l]->thisRnode->paraC2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])/(rn_b_zn[l]->thisRnode->paraS2bf[n][i]*(rn_b_zn[l]->thisRnode->paraC2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])-rn_b_zn[l]->thisRnode->paraSCbf[n][i]))-(costheta2[i]/(rn_b_zn[l]->thisRnode->paraC2bf[n][i]*(rn_b_zn[l]->thisRnode->paraS2bf[n][i]/rn_b_zn[l]->thisRnode->paraSCbf[n][i])-rn_b_zn[l]->thisRnode->paraSCbf[n][i])));
            rn_b_zn[l]->thisRnode->para1bf[n][i]=sum1*rn_b_zn[l]->thisRnode->para21bf[n][i]+sum2*rn_b_zn[l]->thisRnode->para22bf[n][i];
          }

          rn_b_zn[l]->thisRnode->para1bf[n][i]+=-1./llsi;

          rn_b_zn[l]->thisRnode->para1bf[n][i]=1./rn_b_zn[l]->thisRnode->para1bf[n][i];
        }
      }
  }

}

int dyeing(int xl,int yl,int r,int col[nx+2*nr-1][nz+1],int nr){
  int nnc=0;
  for(int i=0;i<nx+2*nr-1;i++){
		for(int j=0;j<nz+1;j++){
      if(((double)i-(double)xl)*((double)i-(double)xl)+((double)j-(double)yl)*((double)j-(double)yl)<=((double)r)*((double)r)){
        col[i][j]=0;
      }else if((i-xl)*(i-xl)+(j-yl)*(j-yl)<=(r+nr)*(r+nr)){
        col[i][j]=2;
        nnc++;
      }else{
        col[i][j]=1;
        nnc++;
      }
		}
	}
  return nnc;
}

void makeClot(double X,double Z,double xl,double zl, int clot[],double u1[glnn2],double u2[glnn2],double theta1[glnn1],double theta2[glnn1]){
    double Xi=fabs(X-xl);
    double Zi=fabs(Z-zl);
    double tanxz0=Xi/Zi;
    if((X>xl&&Z>zl)||(X<xl&&Z<zl)){
        for(int i=0;i<glnn2;i++){
                int mn=-1;
                for(int jj=0;jj<glnn1;jj++){
                    double tan1=-sin(u1[i])*tan(theta2[jj]);
                    if(tan1>tanxz0&&mn==-1){
                        mn=jj;
                        break;
                    }
                }
                if(mn!=-1){
                    clot[i]=mn;
                }else{
                    clot[i]=glnn1;
                }
        }
    }else{
        for(int i=0;i<glnn2;i++){
                int mn=-1;
                for(int jj=0;jj<glnn1;jj++){
                    double tan1=sin(u1[i])*tan(theta1[jj]);
                    if(tan1>tanxz0&&mn==-1){
                        mn=jj;
                        break;
                    }
                }
                if(mn!=-1){
                    clot[i]=mn;
                }else{
                    clot[i]=glnn1;
                }
        }
    }
}

double PD_eprtff_c(double para1,double para2,int i,int j,int m,int n,double I0,realNode *rn[][nz+1],double deltaxh,double para1ff, double para21ff, double para22ff){
    double temp;
    int nrr;
    if((i<nr)||(i>nx+nr-1)){
      temp=rn[i][j]->tempff[m][n];
    }else if(j==0){
      temp=rn[i][j]->tempff[m][n];
    }    
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;
        
        for(int l=0;l<rn[i][j]->HNnode;l++){
          if((-1.*rn[i][j]->cXiX[l])*para1+(-1.*rn[i][j]->cXiZ[l])*para2>0.){
            sumx+=rn[i][j]->Xitempff[l][m][n]*rn[i][j]->HrnXiX[l];
            sumz+=rn[i][j]->Xitempff[l][m][n]*rn[i][j]->HrnXiZ[l];
          }
        }
        dI1=(sumx*para21ff)+(sumz*para22ff)-(I0*para3);
        temp=dI1*para1ff;
      }

	return temp;
}

double PD_eprtff_bc_z0(double para1,double para2,int l,int m,int n,double I0,realNode_boundry *rn_b_z0[],double deltaxh,double para1ff, double para21ff, double para22ff){
    double temp;
    int nrr;

    if(l==n_boundry){
      temp=(-1.*rn_b_z0[l]->thisRnode->qb)/pi;
    }
    else if((l>(n_boundry/2))&&(n<rn_b_z0[l]->clot[m])&&(l!=n_boundry)){
      temp=rn_b_z0[l]->qdmm;
    }
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;

        for(int ii=0;ii<rn_b_z0[l]->thisRnode->HNnode;ii++){
          if((-1.*rn_b_z0[l]->thisRnode->cXiX[ii])*para1+(-1.*rn_b_z0[l]->thisRnode->cXiZ[ii])*para2>0.){
            sumx+=rn_b_z0[l]->thisRnode->Xitempff[ii][m][n]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
            sumz+=rn_b_z0[l]->thisRnode->Xitempff[ii][m][n]*rn_b_z0[l]->thisRnode->HrnXiZ[ii];
          }
        }
        dI1=(sumx*para21ff)+(sumz*para22ff)-(I0*para3);
        temp=dI1*para1ff;
      }
	return temp;
}

double PD_eprtff_bc_zn(double para1,double para2,int l,int m,int n,double I0,realNode_boundry *rn_b_zn[],double deltaxh,double para1ff, double para21ff, double para22ff){
    double temp;
    int nrr;

    if(l==(n_boundry-2)/2){
      temp=(-1.*rn_b_zn[l]->thisRnode->qd)/pi;
    }
    else if(l>((n_boundry-2)/2)){
      temp=rn_b_zn[l]->qdmm;
    }

    else if((l<((n_boundry-2)/2))&&(n>rn_b_zn[l]->clot[m])){
      temp=rn_b_zn[l]->qdmm;
    }
    
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;

        for(int ii=0;ii<rn_b_zn[l]->thisRnode->HNnode;ii++){
          if((-1.*rn_b_zn[l]->thisRnode->cXiX[ii])*para1+(-1.*rn_b_zn[l]->thisRnode->cXiZ[ii])*para2>0.){
            sumx+=rn_b_zn[l]->thisRnode->Xitempff[ii][m][n]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
            sumz+=rn_b_zn[l]->thisRnode->Xitempff[ii][m][n]*rn_b_zn[l]->thisRnode->HrnXiZ[ii];
          }
        }
        dI1=(sumx*para21ff)+(sumz*para22ff)-(I0*para3);
        temp=dI1*para1ff;
      }

	return temp;
}

double PD_eprtbf_c(double para1,double para2,int i,int j,int m,int n,double I0,realNode *rn[][nz+1],double deltaxh,double para1bf, double para21bf, double para22bf){
    double temp;
    int nrr;
    if((i<nr-1)||(i>nx+nr-2)){
      temp=rn[i][j]->tempbf[m][n];
    }else if(j==0){
      temp=rn[i][j]->tempbf[m][n];
    }
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;

        for(int l=0;l<rn[i][j]->HNnode;l++){
          if((-1.*rn[i][j]->cXiX[l])*para1+(-1.*rn[i][j]->cXiZ[l])*para2>0.){
            sumx+=rn[i][j]->Xitempbf[l][m][n]*rn[i][j]->HrnXiX[l];
            sumz+=rn[i][j]->Xitempbf[l][m][n]*rn[i][j]->HrnXiZ[l];
          }
        }

        dI1=(sumx*para21bf)+(sumz*para22bf)-(I0*para3);
        temp=dI1*para1bf;
      }
	return temp;
}

double PD_eprtbf_bc_z0(double para1,double para2,int l,int m,int n,double I0,realNode_boundry *rn_b_z0[],double deltaxh,double para1bf, double para21bf, double para22bf){
    double temp;
    int nrr;
    if(l==0){
      temp=(rn_b_z0[l]->thisRnode->qf)/pi;
    }else if((l<(n_boundry/2))&&(n<rn_b_z0[l]->clot[m])){
      temp=rn_b_z0[l]->qdmm;
    }
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;

        for(int ii=0;ii<rn_b_z0[l]->thisRnode->HNnode;ii++){
          if((-1.*rn_b_z0[l]->thisRnode->cXiX[ii])*para1+(-1.*rn_b_z0[l]->thisRnode->cXiZ[ii])*para2>0.){
            sumx+=rn_b_z0[l]->thisRnode->Xitempbf[ii][m][n]*rn_b_z0[l]->thisRnode->HrnXiX[ii];
            sumz+=rn_b_z0[l]->thisRnode->Xitempbf[ii][m][n]*rn_b_z0[l]->thisRnode->HrnXiZ[ii];
          }
        }
        dI1=(sumx*para21bf)+(sumz*para22bf)-(I0*para3);
        temp=dI1*para1bf;
      }
	return temp;
}

double PD_eprtbf_bc_zn(double para1,double para2,int l,int m,int n,double I0,realNode_boundry *rn_b_zn[],double deltaxh,double para1bf, double para21bf, double para22bf){
    double temp;
    int nrr;
    if(l==(n_boundry-2)/2){
      temp=(-1.*rn_b_zn[l]->thisRnode->qd)/pi;
    }else if(l<((n_boundry-2)/2)){
      temp=rn_b_zn[l]->qdmm;
    }else if((l>(n_boundry/2))&&(n>rn_b_zn[l]->clot[m])){
      temp=rn_b_zn[l]->qdmm;
    }
    else{
        double sum1=0;
        double sum2=0;
        double area=areaSpace;
        double sumx=0;
        double sumz=0;
        double dI1=0;
        double dI2=0;
        double xi=0;

        for(int ii=0;ii<rn_b_zn[l]->thisRnode->HNnode;ii++){
          if((-1.*rn_b_zn[l]->thisRnode->cXiX[ii])*para1+(-1.*rn_b_zn[l]->thisRnode->cXiZ[ii])*para2>0.){
            sumx+=rn_b_zn[l]->thisRnode->Xitempbf[ii][m][n]*rn_b_zn[l]->thisRnode->HrnXiX[ii];
            sumz+=rn_b_zn[l]->thisRnode->Xitempbf[ii][m][n]*rn_b_zn[l]->thisRnode->HrnXiZ[ii];
          }
        }

        dI1=(sumx*para21bf)+(sumz*para22bf)-(I0*para3);
        temp=dI1*para1bf;
      }
	return temp;
}

double makeQDMM(double sintheta1[],double costheta1[],double sintheta2[],double costheta2[],double sinu1[],double sinu2[],double xl,double yl,int clot[],double l,double m,double *Isffd[],double *Isfbd[],double *Isbfd[],double *Isbbd[],double sincostheta1[],double sincostheta2[],double sinsintheta1[],double sinsintheta2[]){
    double qout=0.,qoutx=0.,qoutz=0.;
    double xi=fabs(l-xl);
    double yi=fabs(m-yl);
    double ci=sqrt((xi*xi)+(yi*yi));
    double sinz0=yi/ci;
    double cosz0=xi/ci;

    double **Idmmff,**Idmmfb,**Idmmbf,**Idmmbb;
    Idmmff=new double*[glnn2];Idmmfb=new double*[glnn2];Idmmbf=new double*[glnn2];Idmmbb=new double*[glnn2];

    for(int ii=0;ii<glnn2;ii++){
      Idmmff[ii]=new double[glnn1];Idmmfb[ii]=new double[glnn1];Idmmbf[ii]=new double[glnn1];Idmmbb[ii]=new double[glnn1];
    }
    if(l>xl&&m>yl){
        qoutx+=GL4_I0(Isbbd,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Isbbd,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);
        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmbf[ii][jj]=Isbfd[ii][jj];
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmbf[ii][jj]=0;
            }
        }
        qoutx+=GL4_I0(Idmmbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Idmmbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmfb[ii][jj]=0;
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmfb[ii][jj]=Isfbd[ii][jj];
            }
        }

        qoutx+=GL4_I0(Idmmfb,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Idmmfb,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2);

        qout=-cosz0*qoutx-sinz0*qoutz;
    }else if(l<xl&&m<yl){

        qoutx+=GL4_I0(Isffd,sincostheta1,0.,pi/2.,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Isffd,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmbf[ii][jj]=0;
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmbf[ii][jj]=Isbfd[ii][jj];
            }
        }

        qoutx+=GL4_I0(Idmmbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Idmmbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmfb[ii][jj]=Isfbd[ii][jj];
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmfb[ii][jj]=0;
            }
        }

        qoutx+=GL4_I0(Idmmfb,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Idmmfb,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2);

        qout=cosz0*qoutx+sinz0*qoutz;
    }else if(l<xl&&m>yl){

        qoutx+=GL4_I0(Isfbd,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Isfbd,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmff[ii][jj]= Isffd[ii][jj];
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmff[ii][jj]=0;
            }
        }

        qoutx+=GL4_I0(Idmmff,sincostheta1,0.,pi/2.,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Idmmff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmbb[ii][jj]=0;
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmbb[ii][jj]=Isbbd[ii][jj];
            }
        }

        qoutx+=GL4_I0(Idmmbb,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Idmmbb,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

        qout=cosz0*qoutx-sinz0*qoutz;
    }else if(l>xl&&m<yl){
        qoutx+=GL4_I0(Isbfd,sincostheta2,pi/2.,pi,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Isbfd,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);

        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmff[ii][jj]= 0;
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmff[ii][jj]=Isffd[ii][jj];
            }
        }
        qoutx+=GL4_I0(Idmmff,sincostheta1,0.,pi/2.,n1,0.,pi,n2);
        qoutz+=GL4_q_2(Idmmff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2);
        
        for(int ii=0;ii<glnn2;ii++){
            for(int jj=0;jj<clot[ii];jj++){
                Idmmbb[ii][jj]=Isbbd[ii][jj];
            }
            for(int jj=clot[ii];jj<glnn1;jj++){
                Idmmbb[ii][jj]=0;
            }
        }
        qoutx+=GL4_I0(Idmmbb,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
        qoutz+=GL4_q_2(Idmmbb,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

        qout=-cosz0*qoutx+sinz0*qoutz;
    }

    for(int ii=0;ii<glnn2;ii++){
      delete Idmmff[ii];delete Idmmfb[ii];delete Idmmbf[ii];delete Idmmbb[ii];
    }
    delete Idmmff;delete Idmmfb;delete Idmmbf;delete Idmmbb;
    
    return qout;
}


int main(){
    clock_t start,finish;
    cout<<"kn: "<<kn<<endl;
    double beta,residual=10.;
    double theta1[glnn1],theta2[glnn1],u1[glnn2],u2[glnn2],Iff[nx+2*nr-1][nz+1],Ifb[nx+2*nr-1][nz+1],Ibf[nx+2*nr-1][nz+1],Ibb[nx+2*nr-1][nz+1];

    double qf[nx+2*nr-1][nz+1],qb[nx+2*nr-1][nz+1],qu[nx+2*nr-1][nz+1],qd[nx+2*nr-1][nz+1];
    double T[nx+2*nr-1][nz+1],k[nx+2*nr-1],x[nx+2*nr-1],y[nz+1],qxx[nx+2*nr-1];
    double q[nx+2*nr-1][nz+1],TT,I0[nx+2*nr-1][nz+1],I00[nx+2*nr-1][nz+1],hf[nx+2*nr-1][nz+1];
    double totaltime;
    double Ifx[nx+2*nr-1][nz+1],Ibx[nx+2*nr-1][nz+1];

    int HNnode[nx+2*nr-1][nz+1];
    realHlist *Hrn[nx+2*nr-1][nz+1];

    double sintheta1[glnn1],sintheta2[glnn1],sinu1[glnn2],sinu2[glnn2];
    double costheta1[glnn1],costheta2[glnn1],cosu1[glnn2],cosu2[glnn2];
    double para11[glnn1],para21[glnn1][glnn2];
    double para12[glnn1],para22[glnn1][glnn2];
    double sincostheta1[glnn1],sincostheta2[glnn1];
    double sinsintheta1[glnn1],sinsintheta2[glnn1];

    realNode *rn[nx+2*nr-1][nz+1];

    int col[nx+2*nr-1][nz+1];

    double qx0,kll,kl;

    start=clock();  //timing

    int nc=0;

    nc=dyeing(ncx,ncz,ncr,col,nr);

  double areaB=((double)(nx+1)*(double)(nz+1)*areaSpace)-(double)nc*areaSpace-pi*cr*cr;
  double areaBn=areaB/(2.*(double)n_boundry);

    for(int l=0;l<nx+2*nr-1;l++){
      for(int m=0;m<nz+1;m++){
        rn[l][m]=new realNode;
        rn[l][m]->X=(double)l*deltax;
        rn[l][m]->Z=(double)m*deltaz;

        rn[l][m]->area=areaSpace;

        rn[l][m]->tempff=new double*[glnn2];
        rn[l][m]->tempbf=new double*[glnn2];
        rn[l][m]->paraS2ff=new double*[glnn2];
        rn[l][m]->paraS2bf=new double*[glnn2];
        rn[l][m]->paraC2ff=new double*[glnn2];
        rn[l][m]->paraC2bf=new double*[glnn2];
        rn[l][m]->paraSCff=new double*[glnn2];
        rn[l][m]->paraSCbf=new double*[glnn2];

        rn[l][m]->para1ff=new double*[glnn2];
        rn[l][m]->para1bf=new double*[glnn2];

        rn[l][m]->para21ff=new double*[glnn2];
        rn[l][m]->para21bf=new double*[glnn2];
        rn[l][m]->para22ff=new double*[glnn2];
        rn[l][m]->para22bf=new double*[glnn2];
        for(int k=0;k<glnn2;k++){
          rn[l][m]->tempff[k]=new double[glnn1];
          rn[l][m]->tempbf[k]=new double[glnn1];
          rn[l][m]->paraS2ff[k]=new double[glnn1];
          rn[l][m]->paraS2bf[k]=new double[glnn1];
          rn[l][m]->paraC2ff[k]=new double[glnn1];
          rn[l][m]->paraC2bf[k]=new double[glnn1];
          rn[l][m]->paraSCff[k]=new double[glnn1];
          rn[l][m]->paraSCbf[k]=new double[glnn1];

          rn[l][m]->para1ff[k]=new double[glnn1];
          rn[l][m]->para1bf[k]=new double[glnn1];

          rn[l][m]->para21ff[k]=new double[glnn1];
          rn[l][m]->para21bf[k]=new double[glnn1];
          rn[l][m]->para22ff[k]=new double[glnn1];
          rn[l][m]->para22bf[k]=new double[glnn1];
        }
      }
    }

    realNode_boundry *rn_b_z0[n_boundry+1];
    realNode_boundry *rn_b_zn[n_boundry-1];
    for(int l=0;l<n_boundry+1;l++){
      rn_b_z0[l]=new realNode_boundry;
      rn_b_z0[l]->thisRnode=new realNode;
      rn_b_z0[l]->thisRnode->X=((double)ncx*deltax)-(deltax*(double)ncr*cos(pi*l/n_boundry));
      rn_b_z0[l]->thisRnode->Z=((double)ncz*deltaz)-(deltaz*(double)ncr*sin(pi*l/n_boundry));

      rn_b_z0[l]->thisRnode->area=areaSpace;
      rn_b_z0[l]->thisRnode->T=((rn_b_z0[l]->thisRnode->X-(double)(nr-1)*deltax)/L)*(Txn-Tx0)+Tx0;

      rn_b_z0[l]->thisRnode->tempff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->tempbf=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraS2ff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraS2bf=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraC2ff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraC2bf=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraSCff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->paraSCbf=new double*[glnn2];

      rn_b_z0[l]->thisRnode->para1ff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->para1bf=new double*[glnn2];

      rn_b_z0[l]->thisRnode->para21ff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->para21bf=new double*[glnn2];
      rn_b_z0[l]->thisRnode->para22ff=new double*[glnn2];
      rn_b_z0[l]->thisRnode->para22bf=new double*[glnn2];
      for(int k=0;k<glnn2;k++){
        rn_b_z0[l]->thisRnode->tempff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->tempbf[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraS2ff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraS2bf[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraC2ff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraC2bf[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraSCff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->paraSCbf[k]=new double[glnn1];

        rn_b_z0[l]->thisRnode->para1ff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->para1bf[k]=new double[glnn1];

        rn_b_z0[l]->thisRnode->para21ff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->para21bf[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->para22ff[k]=new double[glnn1];
        rn_b_z0[l]->thisRnode->para22bf[k]=new double[glnn1];
      }
    }
    rn_b_z0[0]->thisRnode->X=((double)(ncx-ncr)*deltax);
    rn_b_z0[0]->thisRnode->Z=((double)ncz*deltaz);

    rn_b_z0[n_boundry/2]->thisRnode->X=((double)ncx*deltax);
    rn_b_z0[n_boundry/2]->thisRnode->Z=((double)(ncz-ncr)*deltaz);

    rn_b_z0[n_boundry]->thisRnode->X=((double)(ncx+ncr)*deltax);
    rn_b_z0[n_boundry]->thisRnode->Z=((double)ncz*deltaz);

    for(int l=0;l<n_boundry-1;l++){
      rn_b_zn[l]=new realNode_boundry;
      rn_b_zn[l]->thisRnode=new realNode;
      rn_b_zn[l]->thisRnode->X=rn_b_z0[l+1]->thisRnode->X;
      rn_b_zn[l]->thisRnode->Z=2.*((double)ncz*deltaz)-(rn_b_z0[l+1]->thisRnode->Z);

      rn_b_zn[l]->thisRnode->area=areaSpace;
      rn_b_z0[l]->thisRnode->T=((rn_b_z0[l]->thisRnode->X-(double)(nr-1)*deltax)/L)*(Txn-Tx0)+Tx0;

      rn_b_zn[l]->thisRnode->tempff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->tempbf=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraS2ff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraS2bf=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraC2ff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraC2bf=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraSCff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->paraSCbf=new double*[glnn2];

      rn_b_zn[l]->thisRnode->para1ff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->para1bf=new double*[glnn2];

      rn_b_zn[l]->thisRnode->para21ff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->para21bf=new double*[glnn2];
      rn_b_zn[l]->thisRnode->para22ff=new double*[glnn2];
      rn_b_zn[l]->thisRnode->para22bf=new double*[glnn2];
      for(int k=0;k<glnn2;k++){
        rn_b_zn[l]->thisRnode->tempff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->tempbf[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraS2ff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraS2bf[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraC2ff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraC2bf[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraSCff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->paraSCbf[k]=new double[glnn1];

        rn_b_zn[l]->thisRnode->para1ff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->para1bf[k]=new double[glnn1];

        rn_b_zn[l]->thisRnode->para21ff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->para21bf[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->para22ff[k]=new double[glnn1];
        rn_b_zn[l]->thisRnode->para22bf[k]=new double[glnn1];
      }
    }

    for(int l=0;l<nx+2*nr-1;l++){ 
		  x[l]=(l)*deltax;
      for(int i=0;i<nz+1;i++){
        T[l][i]=(((double)(l-nr+1)*deltax)/L)*(Txn-Tx0)+Tx0;
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

    initializeRealNodes_c(rn,Hrn,rn_b_z0,rn_b_zn,col);

    makeparaff_c(rn,sintheta1,sintheta2,costheta1,costheta2,sinu1,sinu2,rn_b_z0,rn_b_zn);
    makeparabf_c(rn,sintheta1,sintheta2,costheta1,costheta2,sinu1,sinu2,rn_b_z0,rn_b_zn);
        for(int l=0;l<nx+2*nr-1;l++){
            for(int m=0;m<nz+1;m++){

                I0[l][m]=C*v*T[l][m]/4./pi;
                qd[l][m]=-I0[l][m]*pi;

                for(int n=0;n<glnn2;n++){
                  for(int i=0;i<glnn1;i++){
                    rn[l][m]->tempff[n][i]=I0[l][m];
                    rn[l][m]->tempbf[n][i]=I0[l][m];
                  }
                }
            }
        }

        for(int l=0;l<n_boundry+1;l++){
          rn_b_z0[l]->thisRnode->I0=C*v*rn_b_z0[l]->thisRnode->T/4./pi;
          rn_b_z0[l]->thisRnode->qf=rn_b_z0[l]->thisRnode->I0*pi;
          rn_b_z0[l]->thisRnode->qb=rn_b_z0[l]->thisRnode->I0*pi;
          rn_b_z0[l]->thisRnode->qu=rn_b_z0[l]->thisRnode->I0*pi;
          rn_b_z0[l]->thisRnode->qd=rn_b_z0[l]->thisRnode->I0*pi;

          makeClot(rn_b_z0[l]->thisRnode->X,rn_b_z0[l]->thisRnode->Z,cx,cz,rn_b_z0[l]->clot,u1,u2,theta1,theta2);

          rn_b_z0[l]->qdmm=rn_b_z0[l]->thisRnode->I0;

          for(int n=0;n<glnn2;n++){
            for(int i=0;i<glnn1;i++){
              rn_b_z0[l]->thisRnode->tempff[n][i]=rn_b_z0[l]->thisRnode->I0;
              rn_b_z0[l]->thisRnode->tempbf[n][i]=rn_b_z0[l]->thisRnode->I0;
            }
          }
        }

        for(int l=0;l<n_boundry-1;l++){
          rn_b_zn[l]->thisRnode->I0=C*v*rn_b_zn[l]->thisRnode->T/4./pi;
          rn_b_zn[l]->thisRnode->qf=rn_b_zn[l]->thisRnode->I0*pi;
          rn_b_zn[l]->thisRnode->qb=rn_b_zn[l]->thisRnode->I0*pi;
          rn_b_zn[l]->thisRnode->qu=rn_b_zn[l]->thisRnode->I0*pi;
          rn_b_zn[l]->thisRnode->qd=rn_b_zn[l]->thisRnode->I0*pi;

          makeClot(rn_b_zn[l]->thisRnode->X,rn_b_zn[l]->thisRnode->Z,cx,cz,rn_b_zn[l]->clot,u1,u2,theta1,theta2);

          rn_b_zn[l]->qdmm=rn_b_zn[l]->thisRnode->I0;

          for(int n=0;n<glnn2;n++){
            for(int i=0;i<glnn1;i++){
              rn_b_zn[l]->thisRnode->tempff[n][i]=rn_b_zn[l]->thisRnode->I0;
              rn_b_zn[l]->thisRnode->tempbf[n][i]=rn_b_zn[l]->thisRnode->I0;
            }
          }
        }
        kl=C*v*llsi/3.;
        kll=C*v*llsi/3.;

        int j=0;
        residual=100000.;
        while(residual>1e-3&&j<100){
            double a=0;

            memcpy(I00,I0,sizeof(I0));

            for(int l=0;l<nx+2*nr-1;l++){
              for(int n=0;n<glnn2;n++){
                for(int i=0;i<glnn1;i++){
                  rn[l][0]->tempff[n][i]=rn[l][nz]->tempff[n][i];
                  rn[l][0]->tempbf[n][i]=rn[l][nz]->tempbf[n][i];
                }
              }
            }

            for(int l=0;l<nr;l++){
              for(int m=0;m<nz+1;m++){
                for(int n=0;n<glnn2;n++){
                  for(int i=0;i<glnn1;i++){
                    rn[l][m]->tempff[n][i]=rn[l+nx][m]->tempff[n][i]+deltaI0;
                  }
                }
              }
            }

            for(int l=nx+nr;l<nx+2*nr-1;l++){
              for(int m=0;m<nz+1;m++){
                for(int n=0;n<glnn2;n++){
                  for(int i=0;i<glnn1;i++){
                    rn[l][m]->tempff[n][i]=rn[l-nx][m]->tempff[n][i]-deltaI0;
                  }
                }
              }
            }

            for(int l=0;l<nr-1;l++){
              for(int m=0;m<nz+1;m++){
                for(int n=0;n<glnn2;n++){
                  for(int i=0;i<glnn1;i++){
                    rn[l][m]->tempbf[n][i]=rn[l+nx][m]->tempbf[n][i]+deltaI0;
                  }
                }
              }
            }

            for(int l=nx+nr-1;l<nx+2*nr-1;l++){
              for(int m=0;m<nz+1;m++){
                for(int n=0;n<glnn2;n++){
                  for(int i=0;i<glnn1;i++){
                    rn[l][m]->tempbf[n][i]=rn[l-nx][m]->tempbf[n][i]-deltaI0;
                  }
                }
              }
            }

            for(int l=0;l<n_boundry+1;l++){
              rn_b_z0[l]->thisRnode->I00=rn_b_z0[l]->thisRnode->I0;
            }

            for(int l=0;l<n_boundry-1;l++){
              rn_b_zn[l]->thisRnode->I00=rn_b_zn[l]->thisRnode->I0;
            }
            for(int l=nr;l<nx+nr;l++){
  		  	    for(int m=0;m<nz+1;m++){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                            rn[l][m]->tempff[n][i]=PD_eprtff_c(costheta1[i],sintheta1[i]*sinu1[n],l,m,n,i,I0[l][m],rn,deltaxh,rn[l][m]->para1ff[n][i],rn[l][m]->para21ff[n][i],rn[l][m]->para22ff[n][i]);
                        }
                    }
  			    }
  		    }  

            for(int l=nx+nr-2;l>=nr-1;l--){
  		  	    for(int m=0;m<nz+1;m++){
                    for(int n=0;n<glnn2;n++){
                        for(int i=0;i<glnn1;i++){
                            rn[l][m]->tempbf[n][i]=PD_eprtbf_c(costheta2[i],sintheta2[i]*sinu1[n],l,m,n,i,I0[l][m],rn,deltaxh,rn[l][m]->para1bf[n][i], rn[l][m]->para21bf[n][i], rn[l][m]->para22bf[n][i]);
                        }
                    }
  			    }
  		    }

          for(int l=0;l<n_boundry+1;l++){
            for(int n=0;n<glnn2;n++){
              for(int i=0;i<glnn1;i++){
                rn_b_z0[l]->thisRnode->tempff[n][i]=PD_eprtff_bc_z0(costheta1[i],sintheta1[i]*sinu1[n],l,n,i,rn_b_z0[l]->thisRnode->I0,rn_b_z0,deltaxh,rn_b_z0[l]->thisRnode->para1ff[n][i],rn_b_z0[l]->thisRnode->para21ff[n][i],rn_b_z0[l]->thisRnode->para22ff[n][i]);
                rn_b_z0[l]->thisRnode->tempbf[n][i]=PD_eprtbf_bc_z0(costheta2[i],sintheta2[i]*sinu1[n],l,n,i,rn_b_z0[l]->thisRnode->I0,rn_b_z0,deltaxh,rn_b_z0[l]->thisRnode->para1bf[n][i],rn_b_z0[l]->thisRnode->para21bf[n][i],rn_b_z0[l]->thisRnode->para22bf[n][i]);
              }
            }
          }

          for(int l=0;l<n_boundry-1;l++){
            for(int n=0;n<glnn2;n++){
              for(int i=0;i<glnn1;i++){
                rn_b_zn[l]->thisRnode->tempff[n][i]=PD_eprtff_bc_zn(costheta1[i],sintheta1[i]*sinu1[n],l,n,i,rn_b_zn[l]->thisRnode->I0,rn_b_zn,deltaxh,rn_b_zn[l]->thisRnode->para1ff[n][i],rn_b_zn[l]->thisRnode->para21ff[n][i],rn_b_zn[l]->thisRnode->para22ff[n][i]);
                rn_b_zn[l]->thisRnode->tempbf[n][i]=PD_eprtbf_bc_zn(costheta2[i],sintheta2[i]*sinu1[n],l,n,i,rn_b_zn[l]->thisRnode->I0,rn_b_zn,deltaxh,rn_b_zn[l]->thisRnode->para1bf[n][i],rn_b_zn[l]->thisRnode->para21bf[n][i],rn_b_zn[l]->thisRnode->para22bf[n][i]);
              }
            }
          }


            for(int l=0;l<nx+2*nr-1;l++){
  		  	    for(int m=0;m<nz+1;m++){
                      if(m<=(nz/2)){
                        I0[l][m]=(GL4_I0(rn[l][m]->tempff,sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(rn[l][nz-m]->tempff,sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(rn[l][m]->tempbf,sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(rn[l][nz-m]->tempbf,sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;
                      }else{
                        I0[l][m]=I0[l][nz-m];
                      }

                      if(l==(nx/2)+nr-1){
                        for(int n=0;n<glnn2;n++){
                          for(int i=0;i<glnn1;i++){
                            rn[l][m]->tempff[n][i]=rn[l][m]->tempff[n][i]-I0[l][m]+C*v*(T0+deltaT/2.)/4./pi;
                            rn[l][m]->tempbf[n][i]=rn[l][m]->tempbf[n][i]-I0[l][m]+C*v*(T0+deltaT/2.)/4./pi;
                          }
                        }
                        I0[l][m]=C*v*(T0+deltaT/2.)/4./pi;
                      }

                      if(col[l][m]==0){
                        T[l][m]=0;
                      }else{
                        T[l][m]=4.*pi*I0[l][m]/C/v;
                      }

                    a=max(fabs(I0[l][m]-I00[l][m]),a);  
                }
            }

            rn_b_z0[0]->thisRnode->I0=(GL4_I0(rn_b_z0[0]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(rn_b_z0[0]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(rn_b_z0[0]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(rn_b_z0[0]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;
            
            rn_b_z0[0]->thisRnode->qf=GL4_I0(rn_b_z0[0]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,0.,pi,n2)+GL4_I0(rn_b_z0[0]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
            rn_b_z0[0]->thisRnode->qb=GL4_I0(rn_b_z0[0]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2)+GL4_I0(rn_b_z0[0]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
            rn_b_z0[0]->thisRnode->qu=GL4_q_2(rn_b_z0[0]->thisRnode->tempff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2)+GL4_q_2(rn_b_z0[0]->thisRnode->tempbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);
            rn_b_z0[0]->thisRnode->qd=GL4_q_2(rn_b_z0[0]->thisRnode->tempff,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(rn_b_z0[0]->thisRnode->tempbf,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

            rn_b_z0[0]->thisRnode->qx=rn_b_z0[0]->thisRnode->qf+rn_b_z0[0]->thisRnode->qb;

            rn_b_z0[n_boundry]->thisRnode->I0=(GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;
            
            rn_b_z0[n_boundry]->thisRnode->qf=GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,0.,pi,n2)+GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
            rn_b_z0[n_boundry]->thisRnode->qb=GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2)+GL4_I0(rn_b_z0[n_boundry]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
            rn_b_z0[n_boundry]->thisRnode->qu=GL4_q_2(rn_b_z0[n_boundry]->thisRnode->tempff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2)+GL4_q_2(rn_b_z0[n_boundry]->thisRnode->tempbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);
            rn_b_z0[n_boundry]->thisRnode->qd=GL4_q_2(rn_b_z0[n_boundry]->thisRnode->tempff,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(rn_b_z0[n_boundry]->thisRnode->tempbf,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

            rn_b_z0[n_boundry]->thisRnode->qx=rn_b_z0[n_boundry]->thisRnode->qf+rn_b_z0[n_boundry]->thisRnode->qb;

            for(int l=1;l<n_boundry;l++){
              rn_b_z0[l]->thisRnode->I0=(GL4_I0(rn_b_z0[l]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(rn_b_zn[l-1]->thisRnode->tempff,sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(rn_b_z0[l]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(rn_b_zn[l-1]->thisRnode->tempbf,sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;

              rn_b_z0[l]->thisRnode->qf=GL4_I0(rn_b_z0[l]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,0.,pi,n2)+GL4_I0(rn_b_zn[l-1]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
              rn_b_z0[l]->thisRnode->qb=GL4_I0(rn_b_z0[l]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2)+GL4_I0(rn_b_zn[l-1]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
              rn_b_z0[l]->thisRnode->qu=GL4_q_2(rn_b_z0[l]->thisRnode->tempff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2)+GL4_q_2(rn_b_z0[l]->thisRnode->tempbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);
              rn_b_z0[l]->thisRnode->qd=GL4_q_2(rn_b_zn[l-1]->thisRnode->tempff,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(rn_b_zn[l-1]->thisRnode->tempbf,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

              rn_b_z0[l]->thisRnode->qx=rn_b_z0[l]->thisRnode->qf+rn_b_z0[l]->thisRnode->qb;

              rn_b_z0[l]->qdmm=makeQDMM(sintheta1,costheta1,sintheta2,costheta2,sinu1,sinu2,cx,cz,rn_b_z0[l]->clot,rn_b_z0[l]->thisRnode->X,rn_b_z0[l]->thisRnode->Z,rn_b_z0[l]->thisRnode->tempff,rn_b_zn[l-1]->thisRnode->tempff,rn_b_z0[l]->thisRnode->tempbf,rn_b_zn[l-1]->thisRnode->tempbf,sincostheta1,sincostheta2,sinsintheta1,sinsintheta2)/pi;
            }

            for(int l=0;l<n_boundry-1;l++){
              rn_b_zn[l]->thisRnode->I0=rn_b_z0[l+1]->thisRnode->I0;

              rn_b_zn[l]->thisRnode->qf=GL4_I0(rn_b_zn[l]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,0.,pi,n2)+GL4_I0(rn_b_z0[l+1]->thisRnode->tempff,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2);
              rn_b_zn[l]->thisRnode->qb=GL4_I0(rn_b_zn[l]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2)+GL4_I0(rn_b_z0[l+1]->thisRnode->tempbf,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);
              rn_b_zn[l]->thisRnode->qu=GL4_q_2(rn_b_zn[l]->thisRnode->tempff,sinsintheta1,sinu1,0.,pi/2.,n1,0.,pi,n2)+GL4_q_2(rn_b_zn[l]->thisRnode->tempbf,sinsintheta2,sinu1,pi/2.,pi,n1,0.,pi,n2);
              rn_b_zn[l]->thisRnode->qd=GL4_q_2(rn_b_z0[l+1]->thisRnode->tempff,sinsintheta1,sinu2,0.,pi/2.,n1,pi,2.*pi,n2)+GL4_q_2(rn_b_z0[l+1]->thisRnode->tempbf,sinsintheta2,sinu2,pi/2.,pi,n1,pi,2.*pi,n2);

              rn_b_zn[l]->thisRnode->qx=rn_b_zn[l]->thisRnode->qf+rn_b_zn[l]->thisRnode->qb;

              rn_b_zn[l]->qdmm=makeQDMM(sintheta1,costheta1,sintheta2,costheta2,sinu1,sinu2,cx,cz,rn_b_zn[l]->clot,rn_b_zn[l]->thisRnode->X,rn_b_zn[l]->thisRnode->Z,rn_b_zn[l]->thisRnode->tempff,rn_b_z0[l+1]->thisRnode->tempff,rn_b_zn[l]->thisRnode->tempbf,rn_b_z0[l+1]->thisRnode->tempbf,sincostheta1,sincostheta2,sinsintheta1,sinsintheta2)/pi;
            }

            residual=a;
            j++;  
            cout<<j<<endl; 
            cout<<residual<<endl;
        }

    for(int l=0;l<nx+2*nr-1;l++){
  		  	    for(int m=0;m<nz+1;m++){

                    if(m<=(nz/2)){

                        I0[l][m]=(GL4_I0(rn[l][m]->tempff,sintheta1,0.,pi/2.,n1,0.,pi,n2)
                          +GL4_I0(rn[l][nz-m]->tempff,sintheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                          +GL4_I0(rn[l][m]->tempbf,sintheta2,pi/2.,pi,n1,0.,pi,n2)
                          +GL4_I0(rn[l][nz-m]->tempbf,sintheta2,pi/2.,pi,n1,pi,2.*pi,n2))/4./pi;


                      q[l][m]=GL4_I0(rn[l][m]->tempff,sincostheta1,0.,pi/2.,n1,0.,pi,n2)
                        +GL4_I0(rn[l][nz-m]->tempff,sincostheta1,0.,pi/2.,n1,pi,2.*pi,n2)
                        +GL4_I0(rn[l][m]->tempbf,sincostheta2,pi/2.,pi,n1,0.,pi,n2)
                        +GL4_I0(rn[l][nz-m]->tempbf,sincostheta2,pi/2.,pi,n1,pi,2.*pi,n2);

                    }else{
                      I0[l][m]=I0[l][nz-m];
                      q[l][m]=q[l][nz-m];
                    }
                      rn[l][m]->I0=I0[l][m];
                    if(col[l][m]==0){
                      T[l][m]=0;
                    }else{
                      T[l][m]=4.*pi*I0[l][m]/C/v;
                    }
                    rn[l][m]->T=T[l][m];
                }
            }

            for(int l=0;l<nx+2*nr-1;l++){
                qxx[l]=simpson_qx_z(q[l],nz,deltaz);  
            }
            double qxzx[nx+1];

            for(int l=0;l<nx+1;l++){
                qxzx[l]=q[l+nr-1][0];
            }

            double qx1=simpson(qxzx,nx,deltax);

            qx0=(qxx[nr-1]+qxx[nx+nr-1])/2.;
            TT=simpson_deltaT(T[nr-1],T[nx+nr-1],nz,deltaz); 

            for(int l=0;l<nx+2*nr-1;l++){
                k[l]=qx0*L/TT;
            }
            cout<<TT<<endl;

        string file0="PD_kn_"+to_string(kn)+"_Qx_p0.3.xml";
        out_file(file0,q,x,nx, nz,qxx,k);

        file0="PD_kn_"+to_string(kn)+"_T_p0.3.xml";
        out_file_T_ast_2(file0,rn,rn_b_z0,rn_b_zn,nx, nz);

    finish=clock();
    totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
    cout<<"\n time is "<<totaltime<<" s"<<endl; 

    return 0;
}
