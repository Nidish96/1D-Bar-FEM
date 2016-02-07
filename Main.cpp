#include<cstdio>
#include<cstdlib>
#include"System.h"

double bodyforce(double x)
{
  return 1;
}

int main(int args,char* argc[])
{
  SYSTEM SYS(1,2,1);
  SYS.SETBC('d',0.0,0);

  SYS.InitELs(-1,bodyforce);
  SYS.StitchK_BF();
  SYS.Solve();
  SYS.SOLNPRINT(stdout);

  return 0;
}
