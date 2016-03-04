#ifndef SYSTEM_DEFD
#define SYSTEM_DEFD

#include<cstdio>
#include<cstdlib>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>

enum Integration{ Trapezoidal,GaussLegendre,GaussLegendreIso };

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

  ~NODE(){}
};

class ELEMENT{
  int id;
  int O;  /* No. of Nodes in this element */
  int startnode,endnode;  /* id of start and end nodes */
  NODE **NL;  /* Pointers to Nodes in the element */
  gsl_matrix *KL; /* Element Stiffness Matrix */
  gsl_vector *BForce;

public:
  ELEMENT(){O = 0;  NL = (NODE**)malloc(sizeof(NODE*));}
  void setid(int i){id=i;}
  void PushNode(NODE*);
  void setstartnode(int);
  void setendnode(int);

  void SetupK_BF_gl(double (*forcing)(double));  /* Gauss-Legendre Without Isoparametric Mapping */
  void SetupK_BF_gl_iso(double (*forcing)(double));  /* Gauss-Legendre With Isoparametric Mapping */
  void SetupK_BF_tr(double (*forcing)(double)); /* Trapezoidal Without Isoparametric Mapping */

  int retid(){return id;}
  int retO(){return O;}
  int retstartnode(){return startnode;}
  int retendnode(){return endnode;}
  double retKLij(int i,int j){return gsl_matrix_get(KL,i,j);}
  double retBFi(int i){return gsl_vector_get(BForce,i);}

  void LPrint(FILE*);

  void freeall(){free(NL);  gsl_matrix_free(KL);  gsl_vector_free(BForce);}

  ~ELEMENT(){}
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
  int retDid(int i){return DD[i].index;}
  int retFid(int i){return FF[i].index;}
  double retDval(int i){return DD[i].Val;}
  double retFval(int i){return FF[i].Val;}
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

  ~BOUNDARYCONDS(){free(FF); free(DD);}
};

class SYSTEM{
  int id;
  double Length;    /* Total Length of bar */
  int ELNUM,NDNUM;  /* No. of Elements & Nodes */
  ELEMENT *L;
  NODE *N;
  int *NPE;         /* No. of nodes in each element */
  BOUNDARYCONDS BC;
  gsl_matrix *K;     /* Stitched System Stiffness matrix */
  gsl_vector *BODYFORCE;

public:
  SYSTEM(int,int,double); /* id,nodes,length */

  void SETBC(char,double,double); /* Function to take BCs */

  void InitELs(int,double (*forcing)(double)); /* Initialize Elements */
  void StitchK_BF();
  void Solve();
  void SYSPRINT(FILE*);

  void ConditionNumber();

  ~SYSTEM(){gsl_matrix_free(K);
            gsl_vector_free(BODYFORCE);
            free(N);  free(NPE);
            for( int i=0;i<ELNUM;i++ )
              L[i].freeall();
            free(L);}
};

void MatrixPrint(FILE*,gsl_matrix*);
double _A_(double);
double _E_(double);

#endif
