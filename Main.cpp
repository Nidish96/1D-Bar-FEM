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
  int i;
  double P = 1.0;
  double L = 1.0;
  int No_Nodes = 20;
  int No_Elements = -1;
  interval = 100;
  int intgn = 1;

  for( i=1;i<args;i++ )
  {
    argc[i][1] = tolower(argc[i][1]);
    switch(argc[i][1]){
      case 'n': sscanf(argc[i],"-n%d",&No_Nodes); break;
      case 'e': sscanf(argc[i],"-e%d",&No_Elements);  break;
      case 'p': sscanf(argc[i],"-p%lf",&P); break;
      case 'l': sscanf(argc[i],"-l%lf",&L); break;
      case 'i': sscanf(argc[i],"-i%d",&intgn);
                switch(intgn){
                  case 0: Type = Trapezoidal; break;
                  case 1: Type = GaussLegendre; break;
                  case 2: Type = GaussLegendreIso;  break;} break;
      case 'v': sscanf(argc[i],"-v%d",&interval);  break;
      case 'r': break;
      case 'h': i=100;  break;
    }
  }
  if( i==101||args==1 )
  {
    fprintf(stderr,"\n\tFinite Element Solver for 1D Bar Element\n");
    fprintf(stderr,"\tArgs  \tType    \tDescription\n");
    fprintf(stderr,"\t-n    \tint     \tNumber of Nodes [%d]\n",No_Nodes);
    fprintf(stderr,"\t-e    \tint     \tNumber of Elements (-1 implies Nodes-1) [%d]\n",No_Elements);
    fprintf(stderr,"\t-l    \tfl      \tLength of Element in metres [%lf]\n",L);
    fprintf(stderr,"\t-i    \t<0,1,2> \tType of Integration [%d]\n",intgn);
    fprintf(stderr,"\t      \t        \t0 - Trapezoidal\n");
    fprintf(stderr,"\t      \t        \t1 - Gauss-Legendre \n");
    fprintf(stderr,"\t      \t        \t2 - Gauss-Legendre with Isoparametric Mapping\n");
    fprintf(stderr,"\t-v    \tint     \tNumber of Intervals (for Trapezoidal integration) [%d]\n",interval);
    fprintf(stderr,"\t-r    \t-       \tRun without default Parameters\n");
    fprintf(stderr,"\t-h    \t-       \tPrint Usage Message\n");
    exit(0);
  }

  SYSTEM SYS(1,No_Nodes,L);

  SYS.SETBC('d',0.0,0.0);
  SYS.SETBC('f',L,P);

  SYS.InitELs(No_Elements,bodyforce);
  SYS.StitchK_BF();
  SYS.Solve();
  SYS.SYSPRINT(fopen("SYSTEM.FEM","w+"));

  return 0;
}
