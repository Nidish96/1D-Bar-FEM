#include"System.h"
#include<cstdio>
#include<cstdlib>

/* CLASS ELEMENT */
void ELEMENT::setNodes(int o,NODE N[])
{
  O=o;
  NL = (NODE**)malloc(O*sizeof(NODE*));
  for( int i=0;i<O;i++ )
    NL[i] = &N[i];
  startnode = NL[0]->retid();
  endnode = NL[O-1]->retid();
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

  }
}
