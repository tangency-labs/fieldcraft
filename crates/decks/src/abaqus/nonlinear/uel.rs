//! Original generated C user element; the public UEL calling contract supplies the ABI.
//! No vendor example implementation is used. Nodal forces are assembled externally.
pub const SOURCE: &str = r#"/* FieldCraft StVK rod2: constant reference mass, no physical damping.
   Public Abaqus UEL and C interoperability contracts; original implementation.
   FIELDCRAFT_UEL_STANDALONE only selects the Linux ABI for independent C tests. */
#include <math.h>
#include <stdio.h>
#ifdef FIELDCRAFT_UEL_STANDALONE
#define FOR_NAME(lower,upper) lower##_
#else
#include <omi_for_c.h>
#endif
#ifdef __cplusplus
extern "C" {
#endif
extern void FOR_NAME(xit,XIT)(void);
void FOR_NAME(uel,UEL)(double *rhs,double *amat,double *state,double *energy,
 int *ndof,int *nrhs,int *nsvars,double *props,int *nprops,double *coords,
 int *mcrd,int *nnode,double *u,double *du,double *v,double *a,int *jtype,
 double *time,double *dt,int *kstep,int *kinc,int *jelem,double *params,
 int *ndload,int *jdltyp,double *adlmag,double *predef,int *npredf,int *flags,
 int *mlvarx,double *ddlmag,int *mdload,double *pnewdt,int *jprops,
 int *njprop,double *period) {
 (void)du;(void)jtype;(void)time;(void)kstep;(void)kinc;(void)jelem;
 (void)jdltyp;(void)adlmag;(void)predef;(void)npredf;(void)mlvarx;
 (void)ddlmag;(void)mdload;(void)pnewdt;(void)jprops;(void)njprop;(void)period;
 if(*ndof!=6||*nnode!=2||*nrhs!=1||*nsvars!=6||*nprops!=4||*mcrd<3||*ndload!=0||flags[3]!=0) {
  fprintf(stderr,"FieldCraft UEL: unsupported dimensions, distributed loads or perturbation\n");
  FOR_NAME(xit,XIT)();return;
 }
 double chord[3],l2=0,d2=0,f[6],k[36],m[36]={0};
 for(int i=0;i<3;i++) {
  double x=coords[i+*mcrd]-coords[i];l2+=x*x;
  chord[i]=x+u[i+3]-u[i];d2+=chord[i]*chord[i];
 }
 if(!isfinite(l2)||!isfinite(d2)||l2<=0||d2<=1e-20*l2||
    props[0]<=0||props[1]<=0||props[2]<=0||(props[3]!=0&&props[3]!=1)) {
  fprintf(stderr,"FieldCraft UEL: invalid rod geometry or properties\n");FOR_NAME(xit,XIT)();return;
 }
 double l=sqrt(l2),g=(d2-l2)/(2*l2),ea=props[0]*props[1],mt=props[2]*props[1]*l;
 for(int i=0;i<6;i++) {
  double si=i<3?-1:1;rhs[i]=0;
  f[i]=si*ea/l*g*chord[i%3];
  for(int j=0;j<6;j++) {
   int ij=i+6*j;double sj=j<3?-1:1;
   k[ij]=si*sj*ea/l*(chord[i%3]*chord[j%3]/l2+(i%3==j%3?g:0));
   if(i%3==j%3) m[ij]=props[3]==1?(i==j?mt/2:0):(i==j?mt/3:mt/6);
   amat[ij]=0;
  }
 }
 int request=flags[2];
 if(request==2||request==4||request==6) {
  for(int i=0;i<36;i++) amat[i]=request==2?k[i]:m[i];
  if(request==6) for(int i=0;i<6;i++) rhs[i]=-f[i];
  return;
 }
 if(request==3) return; /* no viscous or structural damping */
 if(request!=1) {
  fprintf(stderr,"FieldCraft UEL: unsupported request; DIRECT,NOHAF required\n");FOR_NAME(xit,XIT)();return;
 }
 for(int i=0;i<8;i++) energy[i]=0;
 energy[1]=0.5*ea*l*g*g;
 if(*kinc==0) {
  for(int i=0;i<6;i++) {
   rhs[i]=-f[i];state[i]=f[i];
   for(int j=0;j<6;j++) {
    amat[i+6*j]=k[i+6*j];
    if(flags[0]==11||flags[0]==12) energy[0]+=0.5*v[i]*m[i+6*j]*v[j];
   }
  }
  return; /* output-only call: Abaqus does not commit zero-increment state */
 }
 if(flags[0]==1||flags[0]==2) {
  for(int i=0;i<6;i++) rhs[i]=-f[i];
  for(int i=0;i<36;i++) amat[i]=k[i];
 } else if(flags[0]==11||flags[0]==12) {
  double alpha=params[0],beta=params[1];
  if(*dt<=0||beta<=0||alpha>0||alpha<(-1.0/3.0)) {
   fprintf(stderr,"FieldCraft UEL: unsupported time operator\n");FOR_NAME(xit,XIT)();return;
  }
  for(int i=0;i<6;i++) {
   rhs[i]=-(1+alpha)*f[i]+alpha*state[i];
   for(int j=0;j<6;j++) {
    int ij=i+6*j;rhs[i]-=m[ij]*a[j];
    amat[ij]=m[ij]/(beta*(*dt)*(*dt))+(1+alpha)*k[ij];
    energy[0]+=0.5*v[i]*m[ij]*v[j];
   }
  }
 } else {
  fprintf(stderr,"FieldCraft UEL: unsupported procedure\n");FOR_NAME(xit,XIT)();return;
 }
 for(int i=0;i<6;i++) state[i]=f[i];
 /* Abaqus supplies start-of-increment state and commits only accepted increments. */
}
#ifdef __cplusplus
}
#endif
"#;
