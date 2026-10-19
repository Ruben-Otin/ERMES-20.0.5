// Unit test for linear_solvers/complex_matrix_operations ( all sparse products vs dense references,
// several thread counts, sizes smaller than the thread count ). Build and run with: make test ( in unix/ )
#include "../linear_solvers/complex_matrix_operations.h"
#include <cstdio>
#include <random>
using namespace Kratos;
typedef std::vector<std::vector<ComplexType> > D;
static std::mt19937 rng(1); static std::uniform_real_distribution<double> U(-1,1);
ComplexType rc(){return ComplexType(U(rng),U(rng));}
double maxdiff(ComplexVectorType&a,ComplexVectorType&b){double m=0;for(size_t i=0;i<a.size();i++)m=std::max(m,std::abs(a[i]-b[i]));return m;}
void dense(const D&M,ComplexVectorType&x,ComplexVectorType&y){for(size_t i=0;i<M.size();i++){y[i]=0;for(size_t j=0;j<M.size();j++)y[i]+=M[i][j]*x[j];}}
int main(){
  int fails=0;
  for(int n: {1,2,7,100,1003}) for(int th: {1,2,3,8}){
    // upper-stored matrix
    ComplexMatrixType A(n), F(n);
    D S(n,std::vector<ComplexType>(n)),H=S,Fd=S;
    for(int i=0;i<n;i++) for(int j=i;j<n;j++) if(j==i || U(rng)>0.7){ ComplexType a=rc(); A(i,j)=a; S[i][j]=a;S[j][i]=a; H[i][j]=a; if(j!=i) H[j][i]=std::conj(a); else H[i][i]=a;}
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) if(U(rng)>0.7){ComplexType a=rc(); F(i,j)=a; Fd[i][j]=a;}
    ComplexVectorType x(n),y(n),r(n); for(auto&v:x)v=rc();
    auto conjD=[&](D M){for(auto&r:M)for(auto&v:r)v=std::conj(v);return M;};
    auto trD=[&](D M){D T=M;for(int i=0;i<n;i++)for(int j=0;j<n;j++)T[i][j]=M[j][i];return T;};
    // Herm with possibly complex diagonal: old semantics for transpose uses diag as stored
    struct C{const char*nm; void(*f)(ComplexMatrixType&,ComplexVectorType&,ComplexVectorType&,int); D M; ComplexMatrixType* A;};
    D Ht=trD(H); // cTransMultHerm: diag as stored -> trD keeps diag same
    D Hc=conjD(H); // cConjgMultHerm: conj incl diag
    std::vector<C> cs={{"MultSymm",cMultSymm,S,&A},{"TransMultSymm",cTransMultSymm,S,&A},{"MultHerm",cMultHerm,H,&A},
      {"TransMultHerm",cTransMultHerm,Ht,&A},{"ConjgMultSymm",cConjgMultSymm,conjD(S),&A},{"ConjgMultHerm",cConjgMultHerm,Hc,&A},
      {"MultFull",cMultFull,Fd,&F},{"ConjgMultFull",cConjgMultFull,conjD(Fd),&F},{"TransMultFull",cTransMultFull,trD(Fd),&F}};
    for(auto&c:cs){ for(auto&v:y)v=ComplexType(99,99); c.f(*c.A,x,y,th); dense(c.M,x,r); double d=maxdiff(y,r); if(d>1e-12){printf("FAIL %s n=%d th=%d d=%g\n",c.nm,n,th,d);fails++;} }
    ComplexVectorType z(n); cMult(F,x,z); dense(Fd,x,r); if(maxdiff(z,r)>1e-12){printf("FAIL cMult\n");fails++;}
    cTransposeMult(F,x,z); dense(trD(Fd),x,r); if(maxdiff(z,r)>1e-12){printf("FAIL cTransposeMult\n");fails++;}
  }
  // vector ops large
  int n=50000; ComplexVectorType a(n),b(n),c(n); for(int i=0;i<n;i++){a[i]=rc();b[i]=rc();}
  ComplexType d(0,0),dc(0,0); double nn=0; for(int i=0;i<n;i++){d+=a[i]*b[i];dc+=a[i]*std::conj(b[i]);nn+=std::norm(a[i]);}
  if(std::abs(cDot(a,b)-d)>1e-9||std::abs(cConjgDot(a,b)-dc)>1e-9||std::abs(cTwoNorm(a)-nn)>1e-9){printf("FAIL dots\n");fails++;}
  c=b; cScaleFirstAndAdd(ComplexType(2,1),a,c); double m=0; for(int i=0;i<n;i++)m=std::max(m,std::abs(c[i]-(ComplexType(2,1)*a[i]+b[i]))); if(m>1e-12){printf("FAIL axpy\n");fails++;}
  printf("fails=%d\n",fails); return fails;
}
