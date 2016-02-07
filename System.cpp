#include"System.h"
#include<cstdio>
#include<cstdlib>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_blas.h>

/* SYSTEM FUNCTIONS */
double _A_(double x)
{
  return 1;
}

double _E_(double x)
{
  return 1;
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
    NL[i]->NPrint(fid);
}

void ELEMENT::SetupK()
{
  int i,j;
  double s1=-1,s2=+1,*S,s,x,h,xp,xd;
  gsl_matrix *THTA = gsl_matrix_alloc(O,O);
  gsl_vector *BMX = gsl_vector_alloc(O);
  gsl_vector *XMX = gsl_vector_alloc(O);
  gsl_permutation *PP = gsl_permutation_calloc(O);

  /* TO GET COEFFICIENTS OF ISOPARAMETRIC MAPPING x = sum(ai*s^i),i,0,O-1 */
  S = (double*)malloc(O*sizeof(double));
  for( i=0;i<O;i++ )
  {
    S[i] = s1+i*(s2-s1)/(O-1);
    for( j=0;j<O;j++ )
      gsl_matrix_set(THTA,i,j,pow(S[i],j));
    gsl_vector_set(BMX,i,NL[i]->retX());
  }

  i = 1;
  gsl_linalg_LU_decomp(THTA,PP,&i);
  gsl_linalg_LU_solve(THTA,PP,BMX,XMX);

  gsl_matrix_free(THTA);
  gsl_permutation_free(PP);

  /* XMX Stores the coefficients in, x = a0+a1*s1+a2*s2^2+a3*s3^2+... ,
    isoparametric mapping of x onto s */
  gsl_vector *BMMX = gsl_vector_alloc(O);
  gsl_vector *BPMX = gsl_vector_alloc(O);
  gsl_vector *BPPMX = gsl_vector_alloc(O);

  KL = gsl_matrix_calloc(O,O);

  h = 0.001;
  for( s=s1;s<s2;s+=h )
  {
    x = gsl_vector_get(XMX,0);
    xp = x;
    xd = 0;
    for( j=1;j<O;j++ )
    {
      x += gsl_vector_get(XMX,j)*pow(s,j);
      xp += gsl_vector_get(XMX,j)*pow(s+h,j);
      xd += gsl_vector_get(XMX,j)*j*pow(s,j-1);
    }

    /* Casting BMX as column of first derivative of weights */
    for( i=0;i<O;i++ )
    {
      gsl_vector_set(BMMX,i,1);
      gsl_vector_set(BMX,i,1);
      gsl_vector_set(BPMX,i,1);
      gsl_vector_set(BPPMX,i,1);
      for( j=0;j<O;j++ )
        if( j!=i )
        {
          gsl_vector_set(BMMX,i,gsl_vector_get(BMMX,i)*(s-h-S[i])/(S[i]-S[j]));
          gsl_vector_set(BMX,i,gsl_vector_get(BMX,i)*(s-S[i])/(S[i]-S[j]));
          gsl_vector_set(BPMX,i,gsl_vector_get(BPMX,i)*(s+h-S[i])/(S[i]-S[j]));
          gsl_vector_set(BPPMX,i,gsl_vector_get(BPPMX,i)*(s+2*h-S[i])/(S[i]-S[j]));
        }
      gsl_vector_set(BMMX,i,(gsl_vector_get(BPMX,i)-gsl_vector_get(BMMX,i))/(2*h));
      gsl_vector_set(BPMX,i,(gsl_vector_get(BPPMX,i)-gsl_vector_get(BMX,i))/(2*h));
      gsl_vector_set(BMX,i,gsl_vector_get(BMMX,i));
    }

    for( i=0;i<O;i++ )
      for( j=0;j<O;j++ )
        gsl_matrix_set( KL,i,j,gsl_matrix_get(KL,i,j)+
          ( gsl_vector_get(BMX,i)*gsl_vector_get(BMX,j)*_A_(x)*_E_(x) +
        gsl_vector_get(BPMX,i)*gsl_vector_get(BPMX,j)*_A_(xp)*_E_(xp) )*h/(2*xd));
  }
  printf("ELEMENT %d STIFFNESS MATRIX\n",id);
  MatrixPrint(stdout,KL);

  gsl_vector_free(BMMX);
  gsl_vector_free(BMX);
  gsl_vector_free(BPMX);
  gsl_vector_free(BPPMX);
  gsl_vector_free(XMX);
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
    N[setxid].setid(setxid);
    N[setxid].setX(xbc);
    for(i=NDNUM-1;i>setxid;i--)
    {
      N[i] = N[i-1];
      N[i].setid(i);
      BC.ChCondIndex(i-1,i);
    }

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

void SYSTEM::InitELs(int ER)
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
    L[i].SetupK();
  }
}

void SYSTEM::StitchK()
{
  int i,j,l,ls,le;
  K = gsl_matrix_calloc(NDNUM,NDNUM);

  for( l=0;l<ELNUM;l++ )
  {
    ls = L[l].retstartnode();
    le = L[l].retendnode();

    for( i=ls;i<=le;i++ )
      for( j=ls;j<=le;j++ )
        gsl_matrix_set(K,i,j,gsl_matrix_get(K,i,j) +
          L[l].retKLij(i-ls,j-ls));
  }
  printf("STITCHED MATRIX\n");
  MatrixPrint(stdout,K);
}
