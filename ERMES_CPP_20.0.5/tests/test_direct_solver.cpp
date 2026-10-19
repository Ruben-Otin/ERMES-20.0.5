// Unit test for linear_solvers/direct_complex_solver ( sparse LU ): all ERMES storage modes,
// matrices needing off-diagonal pivoting, tiny sizes and singular detection.
// Build and run with: make test ( in unix/ )
#include "../linear_solvers/direct_complex_solver.h"
#include <cstdio>
#include <random>
using namespace Kratos;
typedef std::complex<double> C;
static std::mt19937 rng(7); static std::uniform_real_distribution<double> U(-1,1);
static C rc(){ return C(U(rng),U(rng)); }

// 3D grid graph ( 7-point stencil ) with complex weights: FEM-like sparsity
static void Grid(int m, std::vector<std::pair<int,int> >& edges){
  auto id=[m](int i,int j,int k){return (i*m+j)*m+k;};
  for(int i=0;i<m;i++)for(int j=0;j<m;j++)for(int k=0;k<m;k++){
    if(i+1<m)edges.push_back({id(i,j,k),id(i+1,j,k)}); if(j+1<m)edges.push_back({id(i,j,k),id(i,j+1,k)}); if(k+1<m)edges.push_back({id(i,j,k),id(i,j,k+1)});}
}
static double RelRes(int n, std::vector<std::vector<C> >& D, Vector<C>& x, Vector<C>& b){
  double r=0,bn=0; for(int i=0;i<n;i++){C s=-b[i]; for(int j=0;j<n;j++) s+=D[i][j]*x[j]; r+=std::norm(s); bn+=std::norm(b[i]);} return std::sqrt(r/bn);
}
int main(){
  int fails=0; String pre("None");
  for(int m: {1,2,5,9}){
    int n=m*m*m; std::vector<std::pair<int,int> > E; Grid(m,E);
    std::vector<std::vector<C> > D(n,std::vector<C>(n,0.0));
    // --- complex symmetric, upper storage ( full-wave mode ); zero diagonal on some rows -> pivoting
    CSRMatrix<C> S(n);
    for(int i=0;i<n;i++){ C d=(i%7==3)?C(0,0):C(6,0.3)+0.1*rc(); if(d!=C(0,0)) S(i,i)=d; D[i][i]=d; }
    for(auto e:E){ C a=C(-1,0)+0.2*rc(); S(e.first,e.second)+=a; D[e.first][e.second]+=a; D[e.second][e.first]+=a; }
    Vector<C> b(n),x(n); for(int i=0;i<n;i++) b[i]=rc();
    { DirectComplexSolver s(1,0,0,1e-10,pre,false); s.Solve(S,b,x); double r=RelRes(n,D,x,b); printf("symm n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
    // same with 4 threads ( parallel subtrees + parallel dense kernels )
    { DirectComplexSolver s(4,0,0,1e-10,pre,false); s.Solve(S,b,x); double r=RelRes(n,D,x,b); printf("symm 4 threads n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
    // --- full storage ( plasma Full_matrix ): non-symmetric
    CSRMatrix<C> F(n); std::vector<std::vector<C> > DF(n,std::vector<C>(n,0.0));
    for(int i=0;i<n;i++){ C d=C(6,1)+rc(); F(i,i)=d; DF[i][i]=d; }
    for(auto e:E){ C a=rc(), c=rc(); F(e.first,e.second)=a; F(e.second,e.first)=c; DF[e.first][e.second]=a; DF[e.second][e.first]=c; }
    { DirectComplexSolver s(1,0,0,1e-10,pre,true); s.Solve(F,b,x); double r=RelRes(n,DF,x,b); printf("full n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
    { DirectComplexSolver s(4,0,0,1e-10,pre,true); s.Solve(F,b,x); double r=RelRes(n,DF,x,b); printf("full 4 threads n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
    // --- Hermitian ( upper ) + symmetric aux ( upper ) and Hermitian + full aux
    CSRMatrix<C> H(n), As(n), Af(n); std::vector<std::vector<C> > DH(n,std::vector<C>(n,0.0)), DHF;
    for(int i=0;i<n;i++){ double d=6+U(rng); H(i,i)=C(d,0); DH[i][i]+=d; }
    for(auto e:E){ C a=rc(); H(e.first,e.second)=a; DH[e.first][e.second]+=a; DH[e.second][e.first]+=std::conj(a); }
    DHF=DH;
    for(int i=0;i<n;i+=3){ C a=C(0,0.5)+0.1*rc(); As(i,i)=a; DH[i][i]+=a; Af(i,i)=a; DHF[i][i]+=a; }
    for(size_t k=0;k<E.size();k+=4){ auto e=E[k]; C a=0.3*rc(), c=0.3*rc(); As(e.first,e.second)=a; DH[e.first][e.second]+=a; DH[e.second][e.first]+=a;
      Af(e.first,e.second)=a; Af(e.second,e.first)=c; DHF[e.first][e.second]+=a; DHF[e.second][e.first]+=c; }
    { DirectComplexSolver s(1,0,0,1e-10,pre,false); s.Solve_HermSymm(H,As,b,x); double r=RelRes(n,DH,x,b); printf("herm-symm n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
    { DirectComplexSolver s(1,0,0,1e-10,pre,false); s.Solve_HermFull(H,Af,b,x); double r=RelRes(n,DHF,x,b); printf("herm-full n=%d relres=%.2e\n",n,r); if(!(r<1e-10)) fails++; }
  }
  // --- singular matrix detection ( two identical rows )
  { CSRMatrix<C> S(3); S(0,0)=1; S(0,1)=1; S(1,1)=1; S(2,2)=0.0; std::vector<int> Ap,Ai; std::vector<C> Ax;
    BuildFullCSC(S,DirectComplexSolver::SYMMETRIC_UPPER,NULL,DirectComplexSolver::FULL,3,Ap,Ai,Ax);
    ComplexSparseLU LU; LU.SetMatrix(3,Ap,Ai,Ax); std::string msg; bool ok=LU.Factorize(0.01,msg); printf("singular detected: %s\n", ok?"NO":"yes"); if(ok) fails++; }
  printf("fails=%d\n",fails); return fails;
}
