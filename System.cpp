#include"System.h"
#include<cstdio>
#include<cstdlib>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_blas.h>

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

void ELEMENT::setSENodeIDs(int i)
{
  startnode = i;
  endnode = startnode + (O-1);
}

void ELEMENT::LPrint(FILE* fid)
{
  fprintf(fid,"ELEMENT %d\n",id);
  for( int i=0;i<O;i++ )
    NL[i]->NPrint(fid);
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

  nsperel = (NDNUM+(ELNUM-1))/ELNUM;

  k = 0;
  for( i=0;i<ELNUM-1;i++ )
    NPE[i] = nsperel;
  NPE[ELNUM-1] = (NDNUM-k);
  for( i=0;i<ELNUM;i++ )
  {
    L[i].setid(i);
    L[i].setSENodeIDs(k);
    for(j=0;j<NPE[i];j++)
      L[i].PushNode(&N[k++]);
    k--;

    /* SET UP ELEMENT STIFFNESS MATRIX */
  }
}
