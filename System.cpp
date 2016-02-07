#include"System.h"
#include<cstdio>
#include<cstdlib>

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
  for( int i=0;i<O;i++ )
    NL[i]->NPrint(fid);
}
