#include<cstdio>
#include<cstdlib>
#include"System.h"

double bodyforce(double x)
{
  return x;
}

int main(int args,char* argc[])
{
  double P = 1.0;
  double L = 1.0;

  SYSTEM SYS(1,20,L);

  SYS.SETBC('d',0.0,0.0);
  SYS.SETBC('f',1.0,P);

  SYS.InitELs(-1,bodyforce);
  SYS.StitchK_BF();
  SYS.Solve();
  SYS.SYSPRINT(fopen("SYSTEM.FEM","w+"));

  return 0;
}
