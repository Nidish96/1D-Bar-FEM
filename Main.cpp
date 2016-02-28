#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cctype>
#include"System.h"

int interval;
Integration Type;

double bodyforce(double x)
{
  return x;
}

int main(int args,char* argc[])
{
  double P = 1.0;
  double L = 1.0;
  int No_Nodes = (args>1)?atoi(argc[1]):20;
  int No_Elements = (args>2)?atoi(argc[2]):-1;
  interval = (args>3)?atoi(argc[3]):100;
  if( args>4 )
  {
    // fprintf(stderr,"\n%s\n\n",argc[4]);
    for( int i=0;i<sizeof(argc[4]);i++ )
      argc[4][i] = tolower(argc[4][i]);
    if( strcmp(argc[4],"gl") )
      Type = Trapezoidal;
    else
      Type = GaussLegendre;
  }

  SYSTEM SYS(1,No_Nodes,L);

  SYS.SETBC('d',0.0,0.0);
  SYS.SETBC('f',1.0,P);

  SYS.InitELs(No_Elements,bodyforce);
  SYS.StitchK_BF();
  SYS.Solve();
  SYS.SYSPRINT(fopen("SYSTEM.FEM","w+"));

  return 0;
}
