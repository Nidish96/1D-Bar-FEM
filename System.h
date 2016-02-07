#ifndef SYSTEM_DEFD
#define SYSTEM_DEFD

#include<cstdio>

class NODE{
  int id;
  double X, F, U; /* Position, Force, Displacement */

public:
  void setid(int i){id=i;}
  void setX(double x){X=x;}
  void setF(double f){F=f;}
  void setU(double u){U=u;}

  int retid(){return id;}
  double retX(){return X;}
  double retF(){return F;}
  double retU(){return U;}

  void NPrint(FILE* fid){fprintf(fid,"NODE %d: %lf\t%lf\t%lf\n",id,X,F,U);}
};

class ELEMENT{
  int id;
  int O;  /* No. of Nodes in this element */
  int startnode,endnode;  /* id of start and end nodes */
  NODE **NL;  /* Pointers to Nodes in the element */

public:
  void setid(int i){id=i;}
  void setNodes(int,NODE[]);

  int retid(){return id;}
  int retO(){return O;}
  int retstartnode(){return startnode;}
  int retendnode(){return endnode;}

  void LPrint(FILE*);
};

#endif
