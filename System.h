#ifndef SYSTEM_DEFD
#define SYSTEM_DEFD

#include<cstdio>
#include<cstdlib>

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

struct COND{
  int index;
  double Val;
};

class BOUNDARYCONDS{
  int F,D;  /* No. of BCs in Force & Displ */
  COND *FF,*DD;

public:
  BOUNDARYCONDS(){F=0;D=0;
    FF = (COND*)malloc(sizeof(COND));
    DD = (COND*)malloc(sizeof(COND));}

  int retD(){return D;}
  int retF(){return F;}
  void PrintFF(FILE* fid,int i){
    if(i<F) fprintf(fid,"BC %d: %d\t%lf\n",i,FF[i].index,FF[i].Val);
    else fprintf(fid,"ERROR - BC.F NOT FOUND");}
  void PrintDD(FILE* fid,int i){
    if(i<D) fprintf(fid,"BC %d: %d\t%lf\n",i,DD[i].index,DD[i].Val);
    else fprintf(fid,"ERROR - BC.D NOT FOUND");}

  void PushForce(int i,double v){F++;
    FF = (COND*)realloc(FF,F*sizeof(COND));
    FF[F-1].index = i;  FF[F-1].Val = v;}
  void PushDisp(int i,double v){D++;
    DD = (COND*)realloc(DD,D*sizeof(COND));
    DD[D-1].index = i;  DD[D-1].Val = v;}
  void ChCondIndex(int,int);
};

class SYSTEM{
  int id;
  double Length;    /* Total Length of bar */
  int ELNUM,NDNUM;  /* No. of Elements & Nodes */
  ELEMENT *L;
  int *NPE;         /* No. of nodes in each element */
  NODE *N;
  BOUNDARYCONDS BC;

public:
  SYSTEM(int,int,double);

  void SETBC(char,double,double);
};

#endif
