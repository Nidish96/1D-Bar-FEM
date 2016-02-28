#include"System.h"
#include<cstdio>
#include<cstdlib>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_integration.h>

extern int interval;
extern Integration Type;

/* SYSTEM FUNCTIONS */
double _A_(double x)
{
  return 1.0;
}

double _E_(double x)
{
  return 1.0;
}

void MatrixPrint(FILE* fid,gsl_matrix *M)
{
  int i,j;
  for( i=0;i<M->size1;i++ )
  {
    for( j=0;j<M->size2;j++ )
      fprintf(fid,"%lf\t",gsl_matrix_get(M,i,j));
    fprintf(fid,"\n");
  }
}

/* CLASS BOUNDARYCONDS */
void BOUNDARYCONDS::ChCondIndex(int old_id,int new_id)
{
  int i;
  for( i=0;i<F;i++ )
    if( FF[i].index==old_id )
    {
      FF[i].index=new_id;
      break;
    }
  for( i=0;i<D;i++ )
    if( DD[i].index==old_id )
    {
      DD[i].index=new_id;
      break;
    }
}

/* CLASS ELEMENT */
void ELEMENT::PushNode(NODE *N)
{
  O++;
  NL = (NODE**)realloc(NL,O*sizeof(NODE*));
  NL[O-1] = N;
}

void ELEMENT::setstartnode(int i)
{
  startnode = i;
}

void ELEMENT::setendnode(int i)
{
  endnode = i;
}

void ELEMENT::LPrint(FILE* fid)
{
  fprintf(fid,"ELEMENT %d\n",id);
  for( int i=0;i<O;i++ )
  {
    fprintf(fid,"\t");
    NL[i]->NPrint(fid);
  }
}

void ELEMENT::SetupK_BF_tr(double (*forcing)(double))
{
  int i,j;
  double x,h,xs,xf;
  gsl_vector *NN = gsl_vector_alloc(O);
  gsl_vector *NNP = gsl_vector_alloc(O);
  gsl_vector *NNPP = gsl_vector_alloc(O);
  gsl_vector *BB = gsl_vector_alloc(O);
  gsl_vector *BBP = gsl_vector_alloc(O);

  KL = gsl_matrix_calloc(O,O);
  BForce = gsl_vector_calloc(O);

  xs = NL[0]->retX();
  xf = NL[O-1]->retX();
  h = (double) (xf-xs)/interval;

  for( x=xs;x<xf;x+=h )
  {
    for( i=0;i<O;i++ )
    {
      gsl_vector_set(NN,i,1.0);
      gsl_vector_set(NNP,i,1.0);
      gsl_vector_set(NNPP,i,1.0);
      for( j=0;j<O;j++ )
      {
        if( i!=j )
        {
          gsl_vector_set(NN,i,gsl_vector_get(NN,i)*(x-NL[j]->retX())/(NL[i]->retX()-NL[j]->retX()));
          gsl_vector_set(NNP,i,gsl_vector_get(NNP,i)*(x+h-NL[j]->retX())/(NL[i]->retX()-NL[j]->retX()));
          gsl_vector_set(NNPP,i,gsl_vector_get(NNPP,i)*(x+h+h-NL[j]->retX())/(NL[i]->retX()-NL[j]->retX()));
        }
      }
      gsl_vector_set(BB,i,(gsl_vector_get(NNP,i)-gsl_vector_get(NN,i))/h);
      gsl_vector_set(BBP,i,(gsl_vector_get(NNPP,i)-gsl_vector_get(NNP,i))/h);
    }

    for( i=0;i<O;i++ )
    {
      gsl_vector_set(BForce,i,gsl_vector_get(BForce,i)+
        (gsl_vector_get(NN,i)*forcing(x)+gsl_vector_get(NNP,i)*forcing(x+h))*h/2.0);
      for( j=0;j<O;j++ )
        gsl_matrix_set(KL,i,j,gsl_matrix_get(KL,i,j)+
          (gsl_vector_get(BB,i)*gsl_vector_get(BB,j)*_A_(x)*_E_(x)+
            gsl_vector_get(BBP,i)*gsl_vector_get(BBP,j)*_A_(x+h)*_E_(x+h))*h/2.0);
    }
  }

  gsl_vector_free(NN);
  gsl_vector_free(NNP);
  gsl_vector_free(BB);
}

/* ISOPARAMETRIC MAPPING FUNCTIONS */
/* N Vector */
gsl_vector* NN(double iso_s,int O)
{
  int i,j;
  gsl_vector *ret = gsl_vector_alloc(O);
  for(i=0;i<O;i++)
  {
    gsl_vector_set(ret,i,1.0);
    for(j=0;j<O;j++)
      if(j!=i)
        gsl_vector_set(ret,i,gsl_vector_get(ret,i)*(double)((iso_s+1.0)*(O-1.0)/2.0-j)/(i-j));
  }
  return ret;
}
/* B Vector */
gsl_vector* BB(gsl_vector* NV,double iso_s,int O)
{
  int i,j;
  gsl_vector *ret = gsl_vector_calloc(O);
  for( i=0;i<O;i++ )
  {
    for( j=0;j<O;j++ )
      if( j!=i )
        gsl_vector_set(ret,i,gsl_vector_get(ret,i)+1.0/(iso_s+1.0-2.0*j/(O-1)));
    gsl_vector_set(ret,i,gsl_vector_get(ret,i)*gsl_vector_get(NV,i));
  }
  return ret;
}
/* x Map */
double Map_X(gsl_vector* XX,double iso_s,int O)
{
  double ret;
  gsl_vector* NV = NN(iso_s,O);
  gsl_blas_ddot(NV,XX,&ret);
  gsl_vector_free(NV);
  return ret;
}
/* Jacobian Map */
double Map_Jac(gsl_vector* XX,double iso_s,int O)
{
  double ret;
  gsl_vector* NV = NN(iso_s,O);
  gsl_vector* BV = BB(NV,iso_s,O);

  gsl_blas_ddot(BV,XX,&ret);
  gsl_vector_free(NV);
  gsl_vector_free(BV);
  return ret;
}

void ELEMENT::SetupK_BF_gl_iso(double (*forcing)(double))
{
  KL = gsl_matrix_calloc(O,O);
  BForce = gsl_vector_calloc(O);
  gsl_vector *XX = gsl_vector_alloc(O);
  gsl_vector *Btmp,*Ntmp;

  int intg_i,i,j,n_GL;
  n_GL = (O+1.0)/2.0;
  gsl_integration_glfixed_table *Tab = gsl_integration_glfixed_table_alloc(n_GL);
  double si,wi;
  double x,jac;

  for( i=0;i<O;i++ )
    gsl_vector_set(XX,i,NL[i]->retX());

  for( intg_i=0;intg_i<n_GL;intg_i++ )
  {
    gsl_integration_glfixed_point(-1.0,1.0,intg_i,&si,&wi,Tab);
    x = Map_X(XX,si,O);
    jac = Map_Jac(XX,si,O);
    Ntmp = NN(si,O);
    Btmp = BB(Ntmp,si,O);
    for( i=0;i<O;i++ )
    {
      for( j=0;j<O;j++ )
      {
        gsl_matrix_set(KL,i,j,gsl_matrix_get(KL,i,j)+
        gsl_vector_get(Btmp,i)*gsl_vector_get(Btmp,j)*_A_(x)*_E_(x)/jac*wi);
      }
      gsl_vector_set(BForce,i,gsl_vector_get(BForce,i)+
      gsl_vector_get(Ntmp,i)*forcing(x)*jac*wi);
    }
    gsl_vector_free(Btmp);
    gsl_vector_free(Ntmp);
  }

  gsl_integration_glfixed_table_free(Tab);
  gsl_vector_free(XX);
}

void ELEMENT::SetupK_BF_gl(double (*forcing)(double))
{
  KL = gsl_matrix_calloc(O,O);
  BForce = gsl_vector_calloc(O);
  gsl_vector *XX = gsl_vector_alloc(O);
  gsl_vector *Btmp = gsl_vector_alloc(O),*Ntmp = gsl_vector_alloc(O);

  int intg_i,i,j,n_GL;
  n_GL = (O+1.0)/2.0;
  gsl_integration_glfixed_table *Tab = gsl_integration_glfixed_table_alloc(n_GL);
  double xi,wi;

  for( i=0;i<O;i++ )
    gsl_vector_set(XX,i,NL[i]->retX());

  for( intg_i=0;intg_i<n_GL;intg_i++ )
  {
    gsl_integration_glfixed_point(NL[0]->retX(),NL[O-1]->retX(),intg_i,&xi,&wi,Tab);
    gsl_vector_set_all(Ntmp,1.0);
    gsl_vector_set_zero(Btmp);
    for( i=0;i<O;i++ )
    {
      for( j=0;j<O;j++ )
      {
        if( j!=i )
        {
          gsl_vector_set(Ntmp,i,gsl_vector_get(Ntmp,i)*
          (xi-gsl_vector_get(XX,j))/(gsl_vector_get(XX,i)-gsl_vector_get(XX,j)));
          gsl_vector_set(Btmp,i,gsl_vector_get(Btmp,i)+
          1.0/(xi-gsl_vector_get(XX,j)));
        }
      }
      gsl_vector_set(Btmp,i,gsl_vector_get(Btmp,i)*gsl_vector_get(Ntmp,i));
    }

    for( i=0;i<O;i++ )
    {
      for( j=0;j<O;j++ )
      {
        gsl_matrix_set(KL,i,j,gsl_matrix_get(KL,i,j)+
        gsl_vector_get(Btmp,i)*gsl_vector_get(Btmp,j)*_A_(xi)*_E_(xi)*wi);
      }
      gsl_vector_set(BForce,i,gsl_vector_get(BForce,i)+
      gsl_vector_get(Ntmp,i)*forcing(xi)*wi);
    }
  }

  gsl_integration_glfixed_table_free(Tab);
  gsl_vector_free(XX);
}

/* CLASS SYSTEM */
SYSTEM::SYSTEM(int ii,int ND,double LL)
{
  int i;
  id = ii;
  Length = LL;
  NDNUM = ND;

  N = (NODE*)malloc(NDNUM*sizeof(NODE));
  for( i=0;i<NDNUM;i++ )
  {
    N[i].setid(i);
    N[i].setF(0.0);
    N[i].setU(0.0);
    N[i].setX((double)i*Length/(NDNUM-1));
  }
}

void SYSTEM::SETBC(char a,double xbc,double vbc)
{
  int i;
  for(i=0;i<NDNUM;i++)
    if( xbc==N[i].retX() )
      break;
  if( i!=NDNUM )  /* NODE FOUND AT REQUIRED POINT - INDEX i*/
  {
    switch(a){
      case 'f': N[i].setF(vbc); BC.PushForce(i,vbc);
                break;  /* Setting Force BC */
      case 'd': N[i].setU(vbc); BC.PushDisp(i,vbc);
                break;  /* Setting Disp BC */
      default: printf("INVALID BC CHAR CODE\nQUITTING\n");
                    exit(1);
    }
  }
  else  /* NO NODE FOUND AT REQUIRED POINT */
  {
    int setxid;
    for( i=0;i<NDNUM-1;i++ )
    {
      if( N[i].retX()<xbc && N[i+1].retX()>xbc )
      {
        setxid = i+1; /* x must be set at a node with id i+1 */
        break;
      }
    }

    NDNUM++;
    N = (NODE*)realloc(N,NDNUM*sizeof(NODE));
    for(i=NDNUM-1;i>setxid;i--)
    {
      N[i] = N[i-1];
      N[i].setid(i);
      BC.ChCondIndex(i-1,i);
    }
    N[setxid].setid(setxid);
    N[setxid].setX(xbc);

    switch(a){
      case 'f': N[setxid].setF(vbc);  BC.PushForce(setxid,vbc);
                break;  /* Force BC */
      case 'd': N[setxid].setU(vbc);  BC.PushDisp(setxid,vbc);
                break;  /* Displacement BC */
      default:  printf("INVALID BC CHAR CODE\nQUITTING\n");
                exit(1);
    }
  }
}

void SYSTEM::InitELs(int ER,double (*forcing)(double))
{
  int i,j,k,nsperel;
  ELNUM = (ER!=-1)?ER:NDNUM-1;
  if( ELNUM>=NDNUM )
  {
    printf("INVALID ELEMENT NUMBER\nQUITTING\n");
    exit(1);
  }

  L = (ELEMENT*)malloc(ELNUM*sizeof(ELEMENT));
  NPE = (int*)malloc(ELNUM*sizeof(int));

  double tmp = (double)(NDNUM+(ELNUM-1))/ELNUM+0.5;
  nsperel = (int)tmp;

  k = 0;
  for( i=0;i<ELNUM-1;i++ )
    NPE[i] = nsperel;
  NPE[ELNUM-1] = (NDNUM-(nsperel*(ELNUM-1)-(ELNUM-1)));
  for( i=0;i<ELNUM;i++ )
  {
    L[i].setid(i);
    L[i].setstartnode(k);
    for(j=0;j<NPE[i];j++)
      L[i].PushNode(&N[k++]);
    k--;
    L[i].setendnode(k);

    /* SET UP ELEMENT STIFFNESS MATRIX */
    if( Type==Trapezoidal )
    {
      L[i].SetupK_BF_tr(forcing);
      // fprintf(stderr,"\nTRAPEZOIDAL\n");
    }
    else
    {
      L[i].SetupK_BF_gl(forcing);
      // fprintf(stderr,"\nGauss-Legendre\n");
    }
  }
}

void SYSTEM::StitchK_BF()
{
  int i,j,l,ls,le;
  K = gsl_matrix_calloc(NDNUM,NDNUM);
  BODYFORCE = gsl_vector_calloc(NDNUM);

  for( l=0;l<ELNUM;l++ )
  {
    ls = L[l].retstartnode();
    le = L[l].retendnode();

    for( i=ls;i<=le;i++ )
    {
      gsl_vector_set(BODYFORCE,i,gsl_vector_get(BODYFORCE,i) +
              L[l].retBFi(i-ls));
      for( j=ls;j<=le;j++ )
        gsl_matrix_set(K,i,j,gsl_matrix_get(K,i,j) +
          L[l].retKLij(i-ls,j-ls));
    }
  }
  // printf("STITCHED MATRIX\n");
  // MatrixPrint(stdout,K);
  //
  // printf("STITCHED BODY FORCE VECTOR\n");
  // gsl_vector_fprintf(stdout,BODYFORCE,"%lf");
}

void SYSTEM::Solve()
{
  int dnum = BC.retD(),fnum = BC.retF();
  int i,j,k,flag;
  gsl_matrix *AMX = gsl_matrix_alloc(NDNUM,NDNUM);
  gsl_vector *BMX = gsl_vector_calloc(NDNUM);
  gsl_vector *XMX = gsl_vector_alloc(NDNUM);
  gsl_permutation *PP = gsl_permutation_calloc(NDNUM);

  for( i=0;i<NDNUM;i++ )
    for( j=0;j<NDNUM;j++ )
    {
      flag = 1;
      for( k=0;k<dnum;k++ )
        if( BC.retDid(k)==j )
        {
          flag = 0;
          break;
        }
      if( flag==0 )
      {
        if( i==j )
          gsl_matrix_set(AMX,i,j,-1);
        else
          gsl_matrix_set(AMX,i,j,0);
      }
      else
        gsl_matrix_set(AMX,i,j,gsl_matrix_get(K,i,j));

      flag = 1;
      for( k=0;k<dnum;k++ )
        if( BC.retDid(k)==i )
        {
          flag = 0;
          break;
        }
      if( flag==1 )
        gsl_vector_set(BMX,i,N[i].retF());
    }

  for( k=0;k<dnum;k++ )
    for( i=0;i<NDNUM;i++ )
      gsl_vector_set(BMX,i,gsl_vector_get(BMX,i)-gsl_matrix_get(K,i,BC.retDid(k))*N[i].retU()
                      +gsl_vector_get(BODYFORCE,i));

  i = 1;
  gsl_linalg_LU_decomp(AMX,PP,&i);
  gsl_linalg_LU_solve(AMX,PP,BMX,XMX);

  gsl_permutation_free(PP);
  gsl_matrix_free(AMX);
  gsl_vector_free(BMX);

  for( i=0;i<NDNUM;i++ )
  {
    flag = 1;
    for( k=0;k<dnum;k++ )
    {
      if( BC.retDid(k)==i )
      {
        flag = 0;
        break;
      }
    }
    if( flag==0 )
      N[i].setF(gsl_vector_get(XMX,i));
    else
      N[i].setU(gsl_vector_get(XMX,i));
  }

  gsl_vector_free(XMX);
}

void SYSTEM::SYSPRINT(FILE* fid)
{
  int i,dnum=BC.retD(),fnum=BC.retF();
  fprintf(fid,"1D AXIAL ELEMENT FINITE ELEMENT MODELLING\n\n");
  fprintf(fid,"SYSTEM SPECIFICATIONS\nLength : %lf\nNodes : %d\tElements : %d\n\n",Length,NDNUM,ELNUM);
  fprintf(fid,"BOUNDARY CONDITIONS\n");
  fprintf(fid,"DISPLACEMENT BCs : %d\n",dnum);
  for( i=0;i<dnum;i++ )
    fprintf(fid,"\tDISP_BC %d : %lf @ %d\n",i,BC.retDval(i),BC.retDid(i));
  fprintf(fid,"\nFORCE BCs : %d\n",fnum);
  for( i=0;i<fnum;i++ )
    fprintf(fid,"\tFORCE_BC %d : %lf @ %d\n",i,BC.retFval(i),BC.retFid(i));
  fprintf(fid,"\n");

  fprintf(fid,"SOLUTION\n");
  fprintf(fid,"BY ELEMENT\n");
  for( i=0;i<ELNUM;i++ )
    L[i].LPrint(fid);

  fprintf(fid,"\nBY NODE\n\tid\tX\t\tF\t\tU\n");
  for( i=0;i<NDNUM;i++ )
    fprintf(fid,"\t%d\t%lf\t%lf\t%lf\n",N[i].retid(),N[i].retX(),
                                      N[i].retF(),N[i].retU());

  printf("\nid\tX\tF\tU\n");
  for( i=0;i<NDNUM;i++ )
    printf("%d\t%lf\t%lf\t%lf\n",N[i].retid(),N[i].retX(),
                                  N[i].retF(),N[i].retU());
}
