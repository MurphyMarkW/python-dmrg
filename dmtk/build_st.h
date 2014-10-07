#ifndef __DMTK_BUILD_ST_H__
#define __DMTK_BUILD_ST_H__

#include "site.h"

namespace dmtk
{

BasicOp<double>
build_h12 (const Hami<double> &hij, double tstep, bool imaginary_time)
{
// Build composite Hamiltonian
  Block<double> hij_block;
  hij_block = build_complex_site(hij,Hami<double>(),false,false);  

  BasicOp<double> *_real_hij = hij_block(H<double>().set_name(hij.name().c_str()));  
  if(!_real_hij) { cout << "ERROR\n"; exit(-1); }
  BasicOp<double> real_hij = *_real_hij;
  Vector<double> wk(hij_block.basis().dim());
  real_hij = real_hij.diagonalize(wk);
//  for(int i = 0; i < wk.size(); i++) cout << i << " " << wk[i] << endl;

  int i = 0;
  BMatrix<double>::iterator biter;
  for(biter = real_hij.begin(); biter != real_hij.end(); biter++){
    Matrix<double> &m = (*biter);

    if(fabs(tstep) < 1.e-10){
      m = I<double>();
    } else {
      Matrix<double> u = (*biter);
      Matrix<double> ut(u.rows(),u.cols());
      ut = ctranspose(u);
  
  
      m = double(0.);
      for(int j = 0; j < u.rows(); j++){
        m(j,j) = std::exp(-wk[i++]*tstep);
      }
    
      Matrix<double> aux(u.cols(),m.rows());
      aux = product(m,ut);
      m = product(u,aux);
    }

    for(int l = 0; l < m.rows(); l++)
     for(int j = 0; j < m.rows(); j++) cout << l << " " << j << " " << m(l,j) << endl;
  }
  
  return real_hij;
}


BasicOp<complex<double> >
build_h12 (const Hami<complex<double> > &hij, double tstep, bool imaginary_time)
{
// Build composite Hamiltonian
  Block<complex<double> > hij_block;
  hij_block = build_complex_site(hij, Hami<complex<double> >(), false,false);  

  BasicOp<complex<double> > *_real_hij = hij_block(H<complex<double> >().set_name(hij.name().c_str()));  
  if(!_real_hij) { cout << "ERROR\n"; exit(-1); }
  BasicOp<complex<double> > real_hij = *_real_hij;

  BMatrix<complex<double> >::iterator biter;

  for(biter = real_hij.begin(); biter != real_hij.end(); biter++){
    SubMatrix<complex<double> > &block = (*biter);
    for(int l = 0; l < block.rows(); l++)
     for(int j = 0; j < block.rows(); j++) cout << l << " " << j << " " << block(l,j) << endl;
  }


  Vector<double> wk(hij_block.basis().dim());
  real_hij = real_hij.diagonalize(wk);
//  for(int i = 0; i < wk.size(); i++) cout << i << " " << wk[i] << endl;

  int i = 0;
  for(biter = real_hij.begin(); biter != real_hij.end(); biter++){
    Matrix<complex<double> > &m = (*biter);

    if(fabs(tstep) < 1.e-10){
      m = I<complex<double> >();
    } else {
      Matrix<complex<double> > u = (*biter);
      Matrix<complex<double> > ut(u.rows(),u.cols());
      ut = ctranspose(u);
  
  
      m = complex<double>(0.,0.);
      for(int j = 0; j < u.rows(); j++){
        if(imaginary_time)
          m(j,j) = exp(complex<double>(-wk[i++]*tstep));
        else 
          m(j,j) = exp(complex<double>(0.,-wk[i++]*tstep));
      }
    
      Matrix<complex<double> > aux(u.cols(),m.rows());
      aux = product(m,ut);
      m = product(u,aux);
    }

/*
    for(int l = 0; l < m.rows(); l++)
     for(int j = 0; j < m.rows(); j++) cout << l << " " << j << " " << m(l,j) << endl;
*/
  }
  
  return real_hij;
}


template <class T>
std::vector<BasicOp<T> >
build_st_ops(const Hami<T> &hami, double tstep, bool imaginary_time)
{
  const Lattice &l = hami.lattice();
  std::vector<BasicOp<T> > st_ops;

  for(int i = 0; i < l.size()-1; i++){
    Hami<T> h12(Lattice(2,OBC));
    h12.clear();
    h12.set_use_hc(hami.use_hc());
//    h12.set_use_hc(true);
    h12.set_name(hami.name().c_str());
    Block<T> rl = hami.get_site(i);
    Block<T> rr = hami.get_site(i+1);
    h12.sites[0] = &rl;
    h12.sites[1] = &rr;
    typename Hami<T>::const_iterator hiter;
    for(hiter = hami.begin(); hiter != hami.end(); hiter++){
      const Term<T> &t = *hiter;
      if(t.size() == 1){
        BasicOp<T> ref_op = t[0];
        if(ref_op.site() == i || ref_op.site() == i+1){
          T coef = t.coef();
          coef *= ((ref_op.site() == 0 || ref_op.site() == l.size()-1) ? 1. : 0.5);
          ref_op.set_site((ref_op.site() == i ? 0 : 1));
          h12 += ref_op * coef;
        }
      } else {
        BasicOp<T> ref_op1 = t[0];
        BasicOp<T> ref_op2 = t[1];
        int site1 = ref_op1.site();
        int site2 = ref_op2.site();
        if((site1 == i && site2 == i+1) || (site1 == i+1 && site2 == i)){
          ref_op1.set_site((ref_op1.site() == i ? 0 : 1));
          ref_op2.set_site((ref_op2.site() == i ? 0 : 1));
          h12 += ref_op2*ref_op1 * t.coef();
        }
      }
    }
    h12 += H<T>(0).set_name(h12.name().c_str())*T((i == 0) ? 1. : 0.5);
    h12 += H<T>(1).set_name(h12.name().c_str())*T((i == l.size()-2) ? 1. : 0.5);
    cout << "BUILD ST OPS " << i << " " << h12.description() << endl;
    BasicOp<T> st_op;
    st_op = build_h12(h12, tstep, imaginary_time); 
//    st_ops[i] = st_op;
    st_ops.push_back(st_op);
  }
  return st_ops;
}

template<class T>
std::vector<BasicOp<T> >
default_st_ops(const Hami<T> &st_hami, double dt = 0.)
{
  int lx = st_hami.lattice().size();
  std::vector<BasicOp<T> > st_ops(lx-1);
  for(int i = 0; i < lx-1; i++){
    const Block<T> &rl = st_hami.get_site(i);
    const Block<T> &rr = st_hami.get_site(i+1);
    Basis aux_basis(rl.basis(),rr.basis());
    aux_basis.reorder();

    BasicOp<T> hij;

    PackedBasis::const_iterator pbiter;
    for(pbiter = aux_basis.subspaces().begin(); pbiter != aux_basis.subspaces().end(); pbiter++){
      SubMatrix<T> block(pbiter->qn(),*pbiter,*pbiter);
      block = I<T>();
      hij.push_back(block);
    }
    st_ops[i] = hij;
  }
  return st_ops;
}

} // namespace dmtk

#endif // __DMTK_BUILD_ST_H__
