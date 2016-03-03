#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cctype>
#include<math.h>
#include"System.h"

int interval;
Integration Type;

double *cf;
int cfn;

double *af;
int afn;

double *ef;
int efn;

double bodyforce(double x)
{
  double ret=0;
  for( int i=0;i<cfn;i++ )
    ret += cf[i]*pow(x,i);
  return ret;
}

int main(int args,char* argc[])
{
  int i,j,k;
  double P = 1.0;
  double L = 1.0;
  int No_Nodes = 20;
  int No_Elements = -1;
  interval = 100;
  int intgn = 1;

  /* Body Force Polynomial Coefficients */
  cfn = 2;
  cf = (double*)malloc(cfn*sizeof(double));
  cf[0] = 0;  cf[1] = 1.;

  /* Area Polynomial Coefficients */
  afn = 1;
  af = (double*)malloc(afn*sizeof(double));
  af[0] = 1;

  /* Elasticity Modulus Polynomial Coefficients */
  efn = 1;
  ef = (double*)malloc(efn*sizeof(double));
  ef[0] = 1;

  double **bf;
  int bfn=1;
  bf = (double**)malloc(bfn*sizeof(double*));
  for( i=0;i<bfn;i++ )
    bf[i] = (double*)malloc(2*sizeof(double));  /* first column - position;
                                                   second column - force */
  bf[0][0] = L;  bf[0][1] = P;
  double **bd;
  int bdn=1;
  bd = (double**)malloc(bdn*sizeof(double*));
  for( i=0;i<bfn;i++ )
    bd[i] = (double*)malloc(2*sizeof(double));  /* first column - position;
                                                   second column - disp */
  bd[0][0] = 0;  bd[0][1] = 0;

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
      case 'a': sscanf(argc[i],"-a%d,",&afn); af = (double*)realloc(af,afn*sizeof(double));
                for( j=0;j<afn;j++ )
                {
                  for( k=0;k<strlen(argc[i]);k++ )
                    if( argc[i][k]==',' ){ argc[i] = &argc[i][k];  break;}
                  sscanf(argc[i],",%lf",&af[j]);
                  argc[i] = &argc[i][1];
                }                break;
      case 'y': sscanf(argc[i],"-y%d,",&efn); ef = (double*)realloc(ef,efn*sizeof(double));
                for( j=0;j<efn;j++ )
                {
                  for( k=0;k<strlen(argc[i]);k++ )
                    if( argc[i][k]==',' ){ argc[i] = &argc[i][k];  break;}
                  sscanf(argc[i],",%lf",&ef[j]);
                  argc[i] = &argc[i][1];
                }                break;
      case 'c': sscanf(argc[i],"-c%d,",&cfn); cf = (double*)realloc(cf,cfn*sizeof(double));
                for( j=0;j<cfn;j++ )
                {
                  for( k=0;k<strlen(argc[i]);k++ )
                    if( argc[i][k]==',' ){ argc[i] = &argc[i][k];  break;}
                  sscanf(argc[i],",%lf",&cf[j]);
                  argc[i] = &argc[i][1];
                }                break;
      case 'b': argc[i][2] = tolower(argc[i][2]);
                switch(argc[i][2]){
                  case 'f': bfn++;  bf = (double**)realloc(bf,bfn*sizeof(double*));
                                    bf[bfn-1] = (double*)malloc(2*sizeof(double));
                            sscanf(argc[i],"-bf%lf,%lf",&bf[bfn-1][0],&bf[bfn-1][1]); break;
                  case 'd': bdn++;  bd = (double**)realloc(bd,bdn*sizeof(double*));
                                    bd[bdn-1] = (double*)malloc(2*sizeof(double));
                            sscanf(argc[i],"-bd%lf,%lf",&bd[bdn-1][0],&bd[bdn-1][1]); break;
                  case 'i': printf("\n\tInteractive Boundary Condition Setup \n");

                  break;
                }
      case 'r': break;
      case 'h': i=100;  break;
    }
  }
  if( i==101||args==1 )
  {
    fprintf(stderr,"\n\tFinite Element Solver for 1D Bar Element\n");
    fprintf(stderr,"\tArgs  \tType      \tDescription\n");
    fprintf(stderr,"\t-n    \tint       \tNumber of Nodes [%d]\n",No_Nodes);
    fprintf(stderr,"\t-e    \tint       \tNumber of Elements (-1 implies Nodes-1) [%d]\n",No_Elements);
    fprintf(stderr,"\t-l    \tfl        \tTotal Length in metres [%.2f]\n",L);
    fprintf(stderr,"\t-a    \tint,fl,.. \tCoefficients for Area polynomial [%d",afn);
    for( k=0;k<afn;k++ )  fprintf(stderr,",%.2f",af[k]);  fprintf(stderr,"]\n");
    fprintf(stderr,"\t-y    \tint,fl,.. \tCoefficients for Elasticity Modulus polynomial [%d",efn);
    for( k=0;k<efn;k++ )  fprintf(stderr,",%.2f",ef[k]);  fprintf(stderr,"]\n");
    fprintf(stderr,"\t-c    \tint,fl,.. \tCoefficients for body force polynomial [%d",cfn);
    for( k=0;k<cfn;k++ )  fprintf(stderr,",%.2f",cf[k]);  fprintf(stderr,"]\n");
    fprintf(stderr,"\t-bf   \tfl,fl     \tForce Boundary Condition <pos,force> [%.2f,%.2f]\n",bf[0][0],bf[0][1]);
    fprintf(stderr,"\t-bd   \tfl,fl     \tDispl Boundary Condition <pos,disp> [%.2f,%.2f]\n",bd[0][0],bd[0][1]);
    fprintf(stderr,"\t      \t          \tThe above two may be called multiple times to save multiple boundary conditions\n");
    fprintf(stderr,"\t      \t          \t<num,c[0],c[1],c[2]..,c[num-1]> applied as Sum(c[i]*x^i)\n");
    fprintf(stderr,"\t-i    \t<0,1,2>   \tType of Integration [%d]\n",intgn);
    fprintf(stderr,"\t      \t          \t0 - Trapezoidal\n");
    fprintf(stderr,"\t      \t          \t1 - Gauss-Legendre \n");
    fprintf(stderr,"\t      \t          \t2 - Gauss-Legendre with Isoparametric Mapping\n");
    fprintf(stderr,"\t-v    \tint       \tNumber of Intervals (for Trapezoidal integration) [%d]\n",interval);
    fprintf(stderr,"\t-r    \t-         \tRun without default Parameters\n");
    fprintf(stderr,"\t-h    \t-         \tPrint Usage Message\n");
    exit(0);
  }

  /* System Initialization */
  SYSTEM SYS(1,No_Nodes,L);

  /* Setting Boundary Conditions */
  if( bdn>1 )
  {
    bdn--;  bd = &bd[1];
  }
  if( bfn>1 )
  {
    bfn--;  bf = &bf[1];
  }
  for( i=0;i<bfn;i++ )
    for( j=0;j<bdn;j++ )
      if( bf[i][0]==bd[j][0] )
      {
        fprintf(stderr,"\n\tBOUNDARY DATA ERROR  - Force & Displacement BCs Set at one Position\n");
        exit(1);
      }
  for( i=0;i<bdn;i++ )
  {
    SYS.SETBC('d',bd[i][0],bd[i][1]);
    free(bd[i]);
  }
  for( i=0;i<bfn;i++ )
  {
    SYS.SETBC('f',bf[i][0],bf[i][1]);
    free(bf[i]);
  }

  SYS.InitELs(No_Elements,bodyforce);
  SYS.StitchK_BF();
  SYS.Solve();
  SYS.SYSPRINT(fopen("SYSTEM.FEM","w+"));

  free(af);
  free(cf);
  free(ef);
  return 0;
}
