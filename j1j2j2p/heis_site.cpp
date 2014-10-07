#include <dmtk/dmtk.h>

using namespace dmtk;

static dmtk::Vector<Block<double> > site_block;

template<class T>
void build_heis_site(Block<T> &r, double b = 0) 
{
  int dim = 2;
  Basis basis(2);
 
  QN qn;
  qn["N"] = 1;
  qn["Sz"] = -1;
  basis(0) = State(0,qn);
  qn["Sz"] = 1;
  basis(1) = State(1,qn);
  basis.reorder();

  r.clear();
  r.resize(basis);
  r.set_lattice(Lattice(1,OBC));

  Matrix<T> h(dim,dim);
  Matrix<T> sz(dim,dim), sp(dim,dim), sm(dim,dim);

// tables for Hamiltonian

  sz(0,0) = -0.5;
  sz(1,1) = 0.5;

  sp(0,1) = 1;
  sm(1,0) = 1;

  h = T(-b)*sz;

  r.push_back(Sz<T>(sz,r.basis(),0));
  r.push_back(Splus<T>(sp,r.basis(),0));
  r.push_back(Sminus<T>(sm,r.basis(),0));
  r.push_back(H<T>(h,r.basis(),0));
}

template<class T>
Hami<T>
heis_chain(int ls, size_t bc, double j1z, double j1x, double j2z, double j2x, double j2pz, double j2px)
{
  Hami<T> chain(Lattice(ls,bc));
  Lattice l = chain.lattice();
  chain.set_use_hc(true);

  for(int is = 0; is < ls; is++){
    
    if(l(is)[NN1_RIGHT].pos() != l(is).pos()){
      chain += Sz<T>(l(is).pos()) * Sz<T>(l[NN1_RIGHT].pos()) * T(j1z);
      chain += Splus<T>(l(is).pos()) * Sminus<T>(l[NN1_RIGHT].pos()) * T(0.5 * j1x);
    }

    if(is%2) { // Odd case. 1,3,5,...
      if(j2pz != 0.0 && l(is)[NN3_RIGHT].pos() != l(is).pos())
        chain += Sz<T>(l(is).pos()) * Sz<T>(l[NN3_RIGHT].pos()) * T(j2pz);
      if(j2px != 0.0 && l(is)[NN3_RIGHT].pos() != l(is).pos())
        chain += Splus<T>(l(is).pos()) * Sminus<T>(l[NN3_RIGHT].pos()) * T(0.5 * j2px);
    } else { // Even case. 0,2,4,...
      if(j2z != 0.0 && l(is)[NN3_RIGHT].pos() != l(is).pos())
        chain += Sz<T>(l(is).pos()) * Sz<T>(l[NN3_RIGHT].pos()) * T(j2z);
      if(j2x != 0.0 && l(is)[NN3_RIGHT].pos() != l(is).pos())
        chain += Splus<T>(l(is).pos()) * Sminus<T>(l[NN3_RIGHT].pos()) * T(0.5 * j2x);
    }

    
  }

  build_heis_site<T>(chain.site, 0);

  cout << chain.description() << endl;

  return chain;
}

