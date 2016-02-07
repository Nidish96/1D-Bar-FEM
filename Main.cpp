#include<cstdio>
#include<cstdlib>
#include"System.h"

int main(int args,char* argc[])
{
  SYSTEM SYS(1,8,1);
  SYS.InitELs(-1);
  SYS.StitchK();

  return 0;
}
