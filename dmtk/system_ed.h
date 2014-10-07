#ifndef __DMTK_SYSTEM_ED__
#define __DMTK_SYSTEM_ED__

#include <complex>
#include <iostream>
#include <fstream>
#include <dmtk/dmtk.h>
#define USE_RAMSEY
//#define NO_SAVE
using namespace std;

namespace dmtk
{

template<class T>
class 
SystemED:public System<T>
{
  private:
    typedef System<T> _S;

    Vector<int> nbits0;
    Vector<long int> mask;
    Vector<int> dim0;
    Vector<int> bit0;
    bool _imaginary_time;
    int _n_exp_lanczos;
    long int maxdim;
    int dimbasis;
    bool _use_file;
    bool _use_hash;
    bool _expand_basis;
    int _calc_gap;

  public:
    Vector<long int> basis;
    Vector<int> hash_index;
    Vector<int> collision;

    Vector<T> gs;
    Vector<T> gs1;
    double phase;
    double dt;
    double my_time;

    bool _use_targets;

    SystemED():_imaginary_time(false),_n_exp_lanczos(10),dt(0),my_time(0),dimbasis(0),_use_file(false), _use_hash(true),_expand_basis(false),_calc_gap(0) {}
    SystemED(const Hami<T> &hh, const Lattice &l, const char *name): System<T>(hh,l,name),_imaginary_time(false),_n_exp_lanczos(10),dt(0),my_time(0),dimbasis(0),_use_file(false), _use_hash(true),_expand_basis(false),_calc_gap(0) {}

    void build_masks(void);
    void build_basis(void);
    BasicOp<T> operator()(const BasicOp<T> &ref_op) const
    {
      const Block<T>& b = _S::h.get_site(ref_op.site());
      BasicOp<T>op = ref_op;
      op = *b(op.set_site(0));
      op.set_site(ref_op.site());
      return op;
    }

    void product(const BasicOp<T> &op, long int i, Vector<long int> &res, Vector<T> &coef) const;
    void product(const BasicOp<T> &op1, const BasicOp<T> &op2, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const;
    void product(const BasicOp<T> &op1, const BasicOp<T> &op2, const BasicOp<T> &op3, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const;
    void product(const BasicOp<T> &op1, const BasicOp<T> &op2, const BasicOp<T> &op3, const BasicOp<T> &op4, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const;

    Vector<T> product(const BasicOp<T> &op, const Vector<T> &state, T tcoef = T(1)) const;
    Vector<T> product(const BasicOp<T> &op1, const BasicOp<T> &op2, const Vector<T> &state, T tcoef = T(1)) const;

    Vector<T> product(const Vector<T> &state);
    Vector<T> product(long int state);

    T measure(const BasicOp<T> &op) const;
    T measure(const BasicOp<T> &op1, const BasicOp<T> &op2) const;

    void save_hami();

    int get_state_index(const long int &i) const;
    bool imaginary_time() { return _imaginary_time; } const
    SystemED & set_imaginary_time(bool b) { _imaginary_time = b; return *this; }
    int n_exp_lanczos() { return _n_exp_lanczos; } const
    SystemED & set_n_exp_lanczos(int b) { _n_exp_lanczos = b; return *this; }
    Matrix<T> density_matrix(int pos) const;

    int state(long int i, int site) const
    {
      long int kbits = (i & mask(site));
      kbits = kbits >> bit0(site);
      return kbits;
    }

    void add_state(long int state, int &idx);
    int find_state(long int state) const;
    void set_maxdim(int dim) 
    {
      maxdim = dim;
      hash_index.resize(dim+1);
      collision.resize(dim+1);
      basis.resize(dim+1);
      hash_index = 0;
      collision = 0;
      basis = (long int)(0);
    }
    int get_maxdim() const { return maxdim; }
    int get_dimbasis() const { return dimbasis; }
    void clear_basis() 
    {
      basis.resize(0);
      hash_index.resize(0);
      collision.resize(0);
      dimbasis = 0;
    }

    SystemED & set_use_file(bool b) { _use_file = b; return *this; }
    bool use_file() const { return _use_file; }
    SystemED & set_use_hash(bool b) { _use_hash = b; return *this; }
    bool use_hash() const { return _use_hash; }
    SystemED & set_expand_basis(bool b) { _expand_basis = b; return *this; }
    bool expand_basis() const { return _expand_basis; }
    SystemED & set_calc_gap(int calc_gap) { _calc_gap = calc_gap; return *this; }
    int calc_gap() const { return _calc_gap; }
};


template<class T>
int
SystemED<T>::find_state(long int state) const
{
  long int idum = 1234567*(state+1);
  if(idum < 0) idum = -idum;
  long int index = (long int)(quickran(idum)*maxdim)+1;
  if(index < 0) index = -index;

  while (true){
    if(collision(index) == 0) { index=-index; break; }
    if(basis(hash_index(index)) == state) break;
    if(collision(index) == -1) { index = -index; break; }
    index=collision(index);
  }

  return index;
}

template<class T>
void
SystemED<T>::add_state(long int state, int &idx)
{
  idx = find_state(state);
  if(idx >= 0) return;

  idx = -idx;

  int i = idx;
  if(collision(idx) == -1) { 
    i = maxdim;
    while(i>=0 && collision(i) != 0) i--;
    if(i == -1) { cout << "ERROR: two many states\n"; exit(-1); }
  } 

  basis(dimbasis) = state;
  collision(idx)=i;
  collision(i)=-1;
  hash_index(i)=dimbasis;
  dimbasis++;
  idx=i;
}


template<class T>
void
SystemED<T>::build_masks(void)
{
  int nbits = 0;
  nbits0.resize(_S::h.lattice().size());
  mask.resize(_S::h.lattice().size());
  dim0.resize(_S::h.lattice().size());
  bit0.resize(_S::h.lattice().size());

  for(int i = 0; i < _S::h.lattice().size(); i++){
    const Block<T> &_b = _S::h.get_site(i);
    int d = _b.basis().dim();
    for(int j = 0; j < d; j++)
      cout << "SITE: " << i << " .BASIS STATE " << j << " " << _b.basis()(j).qn().n() << " " << _b.basis()(j).qn().sz() << endl;
    dim0(i) = d;
    for(int j = 6; j >= 0; j--){
      if(IBITS(dim0(i)-1,j) == 1) {
        bit0(i) = nbits;
        nbits += (j+1); 
        nbits0(i) = j+1;
        break;
      }
    }
    cout << "SITE " << i << " " << dim0(i) << " " << bit0(i) << " " << nbits0(i) << endl;
  }

  maxdim = (1 << nbits);
  cout << "MAXDIM = " << nbits << " " << maxdim << endl;

  for(int site = 0; site < _S::h.lattice().size(); site++){
    for(int bit = 0; bit < nbits0(site); bit++) mask(site) |= (1 << bit);
    mask(site) = mask(site) << bit0(site);
  }
}

template<class T>
void
SystemED<T>::build_basis(void)
{
  int dim = 0;
  build_masks();
  basis.resize(maxdim); 
  for(long int i = 0; i < maxdim; i++){
    bool found = true;
    for(int site = 0; site < _S::h.lattice().size(); site++){
      const Block<T> &b = _S::h.get_site(site);
      if(dim0(site) <= state(i,site)) { found = false; break; }
    }
    if(found) {
      QN dqn(0,0);
      for(int site = 0; site < _S::h.lattice().size(); site++){
        const Block<T> &b = _S::h.get_site(site);
// cout << i << " " << site << " " << state(i,site) << " " << b.basis()(state(i,site)).qn().n() << " " << b.basis()(state(i,site)).qn().sz() << endl;
        dqn += b.basis()(state(i,site)).qn();
      }
//cout << dqn.n() << " " << dqn.sz() << " " << this->qnt.n() << " " << this->qnt.sz() << endl;
//cout << "--------------------------------------------------------\n";
      if(dqn.n() == this->qnt.n() && dqn.sz() == this->qnt.sz()) basis(dim++) = i;
//      if(dqn.equal(this->qnt, this->_grand_canonical)) basis(dim++) = i;
    }
  }
  cout << "DIM = " << dim << endl;
  dimbasis = dim;
  basis.resize(dim);
  gs.resize(dim);
}

template<class T>
void
SystemED<T>::product(const BasicOp<T> &op, long int i, Vector<long int> &res, Vector<T> &coef) const
{
  res.resize(op.dim());
  coef.resize(op.dim());
  int index = 0;
  int site = op.site();
  int col = state(i, site);
  Vector<int> dqn(_S::h.lattice().size());
  for(int bit = 1; bit < _S::h.lattice().size(); bit++){
    const Block<T>& b = _S::h.get_site(bit-1);
    dqn(bit) = dqn(bit-1)+b.basis()(state(i,bit-1)).qn().n();
  }
  int sign = 1;
  if(op.fermion()) sign *= SGN(dqn(site));
  for(int row = 0; row < op.dim(); row++){
    if(op(col,row) != T(0)){
      coef(index) = op(col,row)*T(sign);
      res(index) = i;
      for(int bit = 0; bit < nbits0(site); bit++){
        if(IBITS(row, bit) == 0)
          res(index) = IBCLR(res(index), bit+bit0(site));
        else
          res(index) = IBSET(res(index), bit+bit0(site));
      }
      index++;
    } 
  }
  res.resize(index);
  coef.resize(index);
}

template<class T>
void
SystemED<T>::product(const BasicOp<T> &op1, const BasicOp<T> &op2, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const
{
  const BasicOp<T> &m1 = op1;
  const BasicOp<T> &m2 = op2;
  res.resize(m1.dim()*m2.dim());
  coef.resize(m1.dim()*m2.dim());
  int index = 0;
  int site1 = op1.site();
  int site2 = op2.site();
  Vector<int> dqn(_S::h.lattice().size());
  for(int bit = 1; bit < _S::h.lattice().size(); bit++){
    const Block<T>& b = _S::h.get_site(bit-1);
    dqn(bit) = dqn(bit-1)+b.basis()(state(i,bit-1)).qn().n();
  }
  int sign = 1;
  if(op1.fermion()) sign *= SGN(dqn(site1));
  if(op2.fermion()) sign *= SGN(dqn(site2));
// cout << op1.description() << " " << op1.site() << " " << op2.description() << " " << op2.site() << " " << sign << endl;
  if(!use_hc){
    int col1 = state(i, site1);
    int col2 = state(i, site2);
    for(int row1 = 0; row1 < m1.dim(); row1++){
      if(m1(col1,row1) != T(0))
      {
        for(int row2 = 0; row2 < m2.dim(); row2++){
          if(m2(col2,row2) != T(0))
          {
//cout << "M1 " << op1.description() << " " << col1 << " " << row1 << " " << m1(col1,row1) << endl;
//cout << "M2 " << op2.description() << " " << col2 << " " << row2 << " " << m2(col2,row2) << endl;
            coef(index) = T(sign)*m1(col1,row1)*m2(col2,row2);
            res(index) = i;
            for(int bit = 0; bit < nbits0(site1); bit++){
              if(IBITS(row1, bit) == 0)
                res(index) = IBCLR(res(index), bit+bit0(site1));
              else
                res(index) = IBSET(res(index), bit+bit0(site1));
            }
            for(int bit = 0; bit < nbits0(site2); bit++){
              if(IBITS(row2, bit) == 0)
                res(index) = IBCLR(res(index), bit+bit0(site2));
              else
                res(index) = IBSET(res(index), bit+bit0(site2));
            }
            index++;
          }
        }
      } 
    }
  } else {
    int row1 = state(i, site1);
    int row2 = state(i, site2);
    if(op1.fermion() && op2.fermion()) sign = -sign;
    for(int col1 = 0; col1 < m1.dim(); col1++){
      if(m1(col1,row1) != T(0)){
        for(int col2 = 0; col2 < m2.dim(); col2++){
          if(m2(col2,row2) != T(0)){
            coef(index) = std::conj(T(sign)*m1(col1,row1)*m2(col2,row2));
            res(index) = i;
            for(int bit = 0; bit < nbits0(site1); bit++){
              if(IBITS(col1, bit) == 0)
                res(index) = IBCLR(res(index), bit+bit0(site1));
              else
                res(index) = IBSET(res(index), bit+bit0(site1));
            }
            for(int bit = 0; bit < nbits0(site2); bit++){
              if(IBITS(col2, bit) == 0)
                res(index) = IBCLR(res(index), bit+bit0(site2));
              else
                res(index) = IBSET(res(index), bit+bit0(site2));
            }
            index++;
          }
        }
      } 
    }
  }

/*
 if(!use_hc)
 cout << op1.description() << " " << op1.site() << " " << op2.description() << " " << op2.site() << endl;
 else
 cout << op1.description() << " " << op1.site() << " " << op2.description() << " " << op2.site() << " (H.C.) "<< endl;
 for(int j = 0; j < index; j++) cout << i << " " << res(j) << " " << coef(j) << endl;
*/

  res.resize(index);
  coef.resize(index);
}

template<class T>
void
SystemED<T>::product(const BasicOp<T> &op1, const BasicOp<T> &op2, const BasicOp<T> &op3, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const
{
  const BasicOp<T> &m1 = op1;
  const BasicOp<T> &m2 = op2;
  const BasicOp<T> &m3 = op3;
  res.resize(m1.dim()*m2.dim()*m3.dim());
  coef.resize(m1.dim()*m2.dim()*m3.dim());
  int index = 0;
  int site1 = op1.site();
  int site2 = op2.site();
  int site3 = op3.site();
  Vector<int> dqn(_S::h.lattice().size());
  for(int bit = 1; bit < _S::h.lattice().size(); bit++){
    const Block<T>& b = _S::h.get_site(bit-1);
    dqn(bit) = dqn(bit-1)+b.basis()(state(i,bit-1)).qn().n();
  }
  int sign = 1;
  if(op1.fermion()) sign *= SGN(dqn(site1));
  if(op2.fermion()) sign *= SGN(dqn(site2));
  if(op3.fermion()) sign *= SGN(dqn(site3));
// cout << op1.description() << " " << op1.site() << " " << op2.description() << " " << op2.site() << endl;
  if(!use_hc){
    int col1 = state(i, site1);
    int col2 = state(i, site2);
    int col3 = state(i, site3);
    for(int row1 = 0; row1 < m1.dim(); row1++){
      if(m1(col1,row1) != T(0))
      {
        for(int row2 = 0; row2 < m2.dim(); row2++){
          if(m2(col2,row2) != T(0))
          {
            for(int row3 = 0; row3 < m3.dim(); row3++){
              if(m3(col3,row3) != T(0))
              {
                coef(index) = T(sign)*m1(col1,row1)*m2(col2,row2)*m3(col3,row3);
                res(index) = i;
                for(int bit = 0; bit < nbits0(site1); bit++){
                  if(IBITS(row1, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site1));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site1));
                }
                for(int bit = 0; bit < nbits0(site2); bit++){
                  if(IBITS(row2, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site2));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site2));
                }
                for(int bit = 0; bit < nbits0(site3); bit++){
                  if(IBITS(row3, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site3));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site3));
                }
                index++;
              }
            }
          }
        }
      } 
    }
  } else {
    int row1 = state(i, site1);
    int row2 = state(i, site2);
    int row3 = state(i, site3);
    int nfermions = op1.fermion() + op2.fermion() + op3.fermion();
    if(nfermions > 1) sign = -sign;
    for(int col1 = 0; col1 < m1.dim(); col1++){
      if(m1(col1,row1) != T(0)){
        for(int col2 = 0; col2 < m2.dim(); col2++){
          if(m2(col2,row2) != T(0)){
            for(int col3 = 0; col3 < m3.dim(); col3++){
              if(m3(col3,row3) != T(0)){

                coef(index) = std::conj(T(sign)*m1(col1,row1)*m2(col2,row2)*m3(col3,row3));
                res(index) = i;
                for(int bit = 0; bit < nbits0(site1); bit++){
                  if(IBITS(col1, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site1));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site1));
                }
                for(int bit = 0; bit < nbits0(site2); bit++){
                  if(IBITS(col2, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site2));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site2));
                }
                for(int bit = 0; bit < nbits0(site3); bit++){
                  if(IBITS(col3, bit) == 0)
                    res(index) = IBCLR(res(index), bit+bit0(site3));
                  else
                    res(index) = IBSET(res(index), bit+bit0(site3));
                }
                index++;
              }
            }
          }
        }
      } 
    }
  }

  res.resize(index);
  coef.resize(index);
}

template<class T>
void
SystemED<T>::product(const BasicOp<T> &op1, const BasicOp<T> &op2, const BasicOp<T> &op3, const BasicOp<T> &op4, long int i, Vector<long int> &res, Vector<T> &coef, bool use_hc) const
{
  const BasicOp<T> &m1 = op1;
  const BasicOp<T> &m2 = op2;
  const BasicOp<T> &m3 = op3;
  const BasicOp<T> &m4 = op4;
  res.resize(m1.dim()*m2.dim()*m3.dim()*m4.dim());
  coef.resize(m1.dim()*m2.dim()*m3.dim()*m4.dim());
  int index = 0;
  int site1 = op1.site();
  int site2 = op2.site();
  int site3 = op3.site();
  int site4 = op4.site();
  Vector<int> dqn(_S::h.lattice().size());
  for(int bit = 1; bit < _S::h.lattice().size(); bit++){
    const Block<T>& b = _S::h.get_site(bit-1);
    dqn(bit) = dqn(bit-1)+b.basis()(state(i,bit-1)).qn().n();
//cout << i << " " << bit << " " << state(i,bit-1) << " " << b.basis()(state(i,bit-1)).qn().n() << endl;
  }
  int sign = 1;
  if(op1.fermion()) sign *= SGN(dqn(site1));
  if(op2.fermion()) sign *= SGN(dqn(site2));
  if(op3.fermion()) sign *= SGN(dqn(site3));
  if(op4.fermion()) sign *= SGN(dqn(site4));

// cout << op1.description() << " " << op1.site() << " " << op2.description() << " " << op2.site() << " " << op3.description() << " " << op3.site() << " " << op4.description() << " " << op4.site() << " " << sign << endl;
  if(!use_hc){
    int col1 = state(i, site1);
    int col2 = state(i, site2);
    int col3 = state(i, site3);
    int col4 = state(i, site4);
    for(int row1 = 0; row1 < m1.dim(); row1++){
      if(m1(col1,row1) != T(0))
      {
        for(int row2 = 0; row2 < m2.dim(); row2++){
          if(m2(col2,row2) != T(0))
          {
            for(int row3 = 0; row3 < m3.dim(); row3++){
              if(m3(col3,row3) != T(0))
              {
                for(int row4 = 0; row4 < m4.dim(); row4++){
                  if(m4(col4,row4) != T(0))
                  {
                    coef(index) = T(sign)*m1(col1,row1)*m2(col2,row2)*m3(col3,row3)*m4(col4,row4);
                    res(index) = i;
                    for(int bit = 0; bit < nbits0(site1); bit++){
                      if(IBITS(row1, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site1));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site1));
                    }
                    for(int bit = 0; bit < nbits0(site2); bit++){
                      if(IBITS(row2, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site2));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site2));
                    }
                    for(int bit = 0; bit < nbits0(site3); bit++){
                      if(IBITS(row3, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site3));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site3));
                    }
                    for(int bit = 0; bit < nbits0(site4); bit++){
                      if(IBITS(row4, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site4));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site4));
                    }
                    index++;
                  }
                }
              }
            }
          }
        }
      } 
    }
  } else {
    int row1 = state(i, site1);
    int row2 = state(i, site2);
    int row3 = state(i, site3);
    int row4 = state(i, site4);
    int nfermions = op1.fermion() + op2.fermion() + op3.fermion() + op4.fermion();
    if(nfermions == 2 || nfermions == 3) sign = -sign;
    for(int col1 = 0; col1 < m1.dim(); col1++){
      if(m1(col1,row1) != T(0)){
        for(int col2 = 0; col2 < m2.dim(); col2++){
          if(m2(col2,row2) != T(0)){
            for(int col3 = 0; col3 < m3.dim(); col3++){
              if(m3(col3,row3) != T(0)){
                for(int col4 = 0; col4 < m4.dim(); col4++){
                  if(m4(col4,row4) != T(0)){

                    coef(index) = std::conj(T(sign)*m1(col1,row1)*m2(col2,row2)*m3(col3,row3)*m4(col4,row4));
                    res(index) = i;
                    for(int bit = 0; bit < nbits0(site1); bit++){
                      if(IBITS(col1, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site1));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site1));
                    }
                    for(int bit = 0; bit < nbits0(site2); bit++){
                      if(IBITS(col2, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site2));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site2));
                    }
                    for(int bit = 0; bit < nbits0(site3); bit++){
                      if(IBITS(col3, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site3));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site3));
                    }
                    for(int bit = 0; bit < nbits0(site4); bit++){
                      if(IBITS(col4, bit) == 0)
                        res(index) = IBCLR(res(index), bit+bit0(site4));
                      else
                        res(index) = IBSET(res(index), bit+bit0(site4));
                    }
                    index++;
                  }
                }
              }
            }
          }
        }
      } 
    }
  }

  res.resize(index);
  coef.resize(index);
}


template<class T>
int
SystemED<T>::get_state_index(const long int &i) const
{
  if(use_hash()) {
    int index = find_state(i);
    if(index >= 0) index = hash_index(index);
    return index;
  }

  int origin = 0;
  int end = basis.size() - 1;
  int index = -1;

/*
  for(int index = 0; index < basis.size(); index++){
    if(i == basis(index)) return index;
  }
  return -1;
*/

  while(origin <= end){
    int index_old = index;
    index = (origin + end) / 2;

    long int s = basis(index);

    if(s == i) {
//cout << "SEARCH " << i << " " << index << endl;
      return index;
    }
    if(s > i)
      end = index;
    else
      origin = index;

    if(index == index_old){
      if(index == end)
        end = end - 1;
      else
        origin = origin + 1;
    }
  }

  return -1;
}

template<class T>
Vector<T>
SystemED<T>::product(const BasicOp<T> &op, const Vector<T> &state, T tcoef) const
{
  Vector<long int>::const_iterator iter;
  Vector<T> coef;
  Vector<T> res = state;
  Vector<long int> this_res;
  res = T(0);
  for(long int i = 0; i < state.size(); i++){
    product(op, basis(i), this_res, coef);
    for(int j = 0; j < this_res.size(); j++){
      int index = get_state_index(this_res(j));
//      if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
      if(index >= 0) res(index) += coef(j)*state(i)*tcoef;
    }
  }

  return res;
}


template<class T>
Vector<T>
SystemED<T>::product(const BasicOp<T> &op1, const BasicOp<T> &op2, const Vector<T> &state, T tcoef) const
{
  Vector<long int>::const_iterator iter;
  Vector<T> coef;
  Vector<T> res = state;
  Vector<long int> this_res;
  res = T(0);

  bool use_hc = _S::use_hc();
  if(op1.is_diagonal() && op2.is_diagonal()) use_hc = false;
//if(!use_hc) cout << op1.description() << " " << op2.description() << endl;

  for(long int i = 0; i < state.size(); i++){
    product(op1, op2, basis(i), this_res, coef, false);
    for(int j = 0; j < this_res.size(); j++){
      int index = get_state_index(this_res(j));
//      if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
      if(index >= 0) res(index) += coef(j)*state(i)*tcoef;
    }
    if(use_hc){
      product(op1, op2, basis(i), this_res, coef, true);
      for(int j = 0; j < this_res.size(); j++){
        int index = get_state_index(this_res(j));
//        if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
        if(index >= 0) res(index) += coef(j)*state(i)*std::conj(tcoef);
      }
    }
  }

  return res;
}

template<class T>
Vector<T>
SystemED<T>::product(long int state)
{
  Vector<T> res(basis.size());
  Vector<T> coef;
  Vector<long int> this_res;
  res = T(0);

  typename Hami<T>::const_iterator iter;
  for(iter = _S::h.begin(); iter != _S::h.end(); iter++){
    const Term<T>& t = (*iter);
    if(t.size() == 1){
// cout << t.description() << endl;
      BasicOp<T> ref_op = t[0];
      const Block<T>& b = _S::h.get_site(ref_op.site());
      BasicOp<T>op = *b(ref_op.set_site(0));
      op.set_site(t[0].site());
      product(op, state, this_res, coef);
      for(int j = 0; j < this_res.size(); j++){
        int index = get_state_index(this_res(j));
if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
        if(_expand_basis && dimbasis < maxdim) {
           add_state(this_res(j), index); 
           index = hash_index(index);
        }
        if(index >= 0) res(index) += coef(j)*t.coef();
      }
    } else if(t.size() == 2) {
      BasicOp<T> ref_op1 = t[0];
      BasicOp<T> ref_op2 = t[1];
      const Block<T>& b1 = _S::h.get_site(ref_op1.site());
      const Block<T>& b2 = _S::h.get_site(ref_op2.site());
//cout << ref_op1.description() << " " << ref_op2.description() << endl;
      BasicOp<T>op1 = *b1(ref_op1.set_site(0));
      BasicOp<T>op2 = *b2(ref_op2.set_site(0));
      bool use_hc = _S::use_hc();
      if(op1.is_diagonal() && op2.is_diagonal()) use_hc = false;
//if(!use_hc) cout << t.description() << endl; else cout << t.description() << " + h.c\n";
      op1.set_site(t[0].site());
      op2.set_site(t[1].site());

      product(op2, op1, state, this_res, coef, false);
      for(int j = 0; j < this_res.size(); j++){
        int index = get_state_index(this_res(j));
//        if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
        if(_expand_basis && dimbasis < maxdim) {
           add_state(this_res(j), index); 
           index = hash_index(index);
        }
        if(index >= 0) res(index) += coef(j)*t.coef();
      }
      if(use_hc){
        product(op2, op1, state, this_res, coef, true);
        for(int j = 0; j < this_res.size(); j++){
          int index = get_state_index(this_res(j));
//          if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
          if(_expand_basis && dimbasis < maxdim) {
             add_state(this_res(j), index); 
             index = hash_index(index);
          }
          if(index >= 0) res(index) += coef(j)*std::conj(t.coef());
        }
      }
    } else if(t.size() == 3) {
      BasicOp<T> ref_op1 = t[0];
      BasicOp<T> ref_op2 = t[1];
      BasicOp<T> ref_op3 = t[2];
      const Block<T>& b1 = _S::h.get_site(ref_op1.site());
      const Block<T>& b2 = _S::h.get_site(ref_op2.site());
      const Block<T>& b3 = _S::h.get_site(ref_op3.site());
      BasicOp<T>op1 = *b1(ref_op1.set_site(0));
      BasicOp<T>op2 = *b2(ref_op2.set_site(0));
      BasicOp<T>op3 = *b3(ref_op3.set_site(0));
      bool use_hc = _S::use_hc();
      if(op1.is_diagonal() && op2.is_diagonal() && op3.is_diagonal()) use_hc = false;
      op1.set_site(t[0].site());
      op2.set_site(t[1].site());
      op3.set_site(t[2].site());

      product(op3, op2, op1, state, this_res, coef, false);
      for(int j = 0; j < this_res.size(); j++){
        int index = get_state_index(this_res(j));
//        if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
        if(_expand_basis && dimbasis < maxdim) {
           add_state(this_res(j), index); 
           index = hash_index(index);
        }
        if(index >= 0) res(index) += coef(j)*t.coef();
      }
      if(use_hc){
        product(op3, op2, op1, state, this_res, coef, true);
        for(int j = 0; j < this_res.size(); j++){
          int index = get_state_index(this_res(j));
//          if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
          if(_expand_basis && dimbasis < maxdim) {
             add_state(this_res(j), index); 
             index = hash_index(index);
          }
          if(index >= 0) res(index) += coef(j)*std::conj(t.coef());
        }
      }
    } else if(t.size() == 4) {
      BasicOp<T> ref_op1 = t[0];
      BasicOp<T> ref_op2 = t[1];
      BasicOp<T> ref_op3 = t[2];
      BasicOp<T> ref_op4 = t[3];
      const Block<T>& b1 = _S::h.get_site(ref_op1.site());
      const Block<T>& b2 = _S::h.get_site(ref_op2.site());
      const Block<T>& b3 = _S::h.get_site(ref_op3.site());
      const Block<T>& b4 = _S::h.get_site(ref_op4.site());
      BasicOp<T>op1 = *b1(ref_op1.set_site(0));
      BasicOp<T>op2 = *b2(ref_op2.set_site(0));
      BasicOp<T>op3 = *b3(ref_op3.set_site(0));
      BasicOp<T>op4 = *b4(ref_op4.set_site(0));
      bool use_hc = _S::use_hc();
      if(op1.is_diagonal() && op2.is_diagonal() && op3.is_diagonal() && op4.is_diagonal()) use_hc = false;
      op1.set_site(t[0].site());
      op2.set_site(t[1].site());
      op3.set_site(t[2].site());
      op4.set_site(t[3].site());

      product(op4, op3, op2, op1, state, this_res, coef, false);
      for(int j = 0; j < this_res.size(); j++){
        int index = get_state_index(this_res(j));
//        if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
        if(_expand_basis && dimbasis < maxdim) {
           add_state(this_res(j), index); 
           index = hash_index(index);
        }
        if(index >= 0) res(index) += coef(j)*t.coef();
      }
      if(use_hc){
        product(op4, op3, op2, op1, state, this_res, coef, true);
        for(int j = 0; j < this_res.size(); j++){
          int index = get_state_index(this_res(j));
//          if(index == -1) cout << "ERROR: state not found " << this_res(j) << endl;
          if(_expand_basis && dimbasis < maxdim) {
             add_state(this_res(j), index); 
             index = hash_index(index);
          }
          if(index >= 0) res(index) += coef(j)*std::conj(t.coef());
        }
      }
    }
  }
  return res;
}


template<class T>
Vector<T>
SystemED<T>::product(const Vector<T> &state) 
{
  Vector<T> res = state;
  Vector<T> this_res = state;
  res = T(0);
  this_res = T(0);


  typename Hami<T>::const_iterator iter;

/*
  for(iter = _S::h.begin(); iter != _S::h.end(); iter++){
    const Term<T>& t = (*iter);
//cout << t.description() << endl;
    if(t.size() == 1){
      BasicOp<T> ref_op = t[0];
      const Block<T>& b = _S::h.get_site(ref_op.site());
      BasicOp<T>op = *b(ref_op.set_site(0));
      op.set_site(t[0].site());
      res += product(op, state, t.coef()); 
    } else {
      BasicOp<T> ref_op1 = t[0];
      BasicOp<T> ref_op2 = t[1];
      const Block<T>& b1 = _S::h.get_site(ref_op1.site());
      const Block<T>& b2 = _S::h.get_site(ref_op2.site());
//cout << ref_op1.description() << " " << ref_op2.description() << endl;
      BasicOp<T>op1 = *b1(ref_op1.set_site(0));
      BasicOp<T>op2 = *b2(ref_op2.set_site(0));
      op1.set_site(t[0].site());
      op2.set_site(t[1].site());
      res += product(op1, op2, state, t.coef()); 
    }

  }
*/

  for(int i = 0; i < basis.size(); i++){
    this_res = product(basis(i));
    this_res *= state(i);
    res += this_res;
  }
  if(_calc_gap > 0){
    T x = dmtk::product(gs1,res);
    res += x*gs1;
  }

  return res;
}


template<class T>
void
SystemED<T>::save_hami() 
{
  Vector<T> res(basis.size());
  Vector<T> data(basis.size());
  Vector<int> col(basis.size());

  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", _S::name());
  ofstream real_file(file,std::ios::out|std::ios::binary);
  if(!real_file) {
    cout << "*** ERROR: Could not open file " << file << endl;
    return;
  }
#ifdef WITH_BZIP2
  bzip2_stream::bzip2_ostream outputfile(real_file);
#else
  ofstream &outputfile = real_file;
#endif
//////////////////////////////////////////////////////////////////////
  cout << "GENERATING HAMILTONIAN " << file << endl;

  for(int i = 0; i < basis.size(); i++){
    res = product(basis(i));
    int index = 0;
    for(int j = 0; j < basis.size(); j++){
      if(res(j) != T(0) && j <= i){
        data(index) = res(j);
        col(index++) = j;
      }
    }
    col.resize(index);
    col.write(outputfile);
    data.resize(index);
    data.write(outputfile);
  }


  outputfile.close();
}


template<class T>
Vector<T>
product_hami(SystemED<T> &s, Vector<T> &state)
{
  Vector<T> res(s.basis.size());
  res = T(0);
  Vector<T> data(s.basis.size());
  Vector<int> col(s.basis.size());

// Reading file
  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", s.name());
  ifstream real_file(file,std::ios::in|std::ios::binary);
  if(!real_file) {
    cout << "*** ERROR: Could not open file " << file << endl;
    return res;
  }
#ifdef WITH_BZIP2
  bzip2_stream::bzip2_istream inputfile(real_file);
#else
  ifstream &inputfile = real_file;
#endif

  for(int i = 0; i < s.basis.size(); i++){
    col.read(inputfile);    
    data.read(inputfile);    
    for(int j = 0; j < data.size(); j++){
      res(i) += state(col(j))*data(j);
      if(col(j) != i) res(col(j)) += state(i)*std::conj(data(j));
    }
  }

  inputfile.close();

  if(s.calc_gap() > 0){
    T x = dmtk::product(s.gs1,state);
    res += T(x*100000.)*s.gs1;
  }

  return res;
}

template<class T>
Vector<T>
product(SystemED<T> &s, Vector<T> &state)
{
  if(!s.use_file()) return s.product(state); 
  return product_hami(s, state);
}


template<class T>
void
hmatrix_verif(SystemED<T>& s)
{
  Vector<T> aux(s.basis.size());
  Vector<T> aux2(s.basis.size());
  Vector<T> &auxv = aux;
  for(int i = 0; i < aux.size(); i++){
    auxv = T(0);
    auxv(i) = T(1);
    aux2 = T(0);
    aux2 = product(s,aux);
    for(int j = 0; j < auxv.size(); j++){
      cout << "VERIF " << i << " " << j << " " << aux2(j) << endl;
    } 
  }
}

template<class T>
void
full_diag(SystemED<T>& s)
{
  Vector<T> aux(s.basis.size());
  Vector<T> aux2(s.basis.size());
  Vector<T> &auxv = aux;
  Matrix<T> hami(s.basis.size(),s.basis.size());
  for(int i = 0; i < aux.size(); i++){
    auxv = T(0);
    auxv(i) = T(1);
    aux2 = T(0);
    aux2 = product(s,aux);
    for(int j = 0; j < auxv.size(); j++){
      hami(i,j) = aux2(j);
    } 
  }

  double w[s.basis.size()];
  double work[3*s.basis.size()];
  DMTK_int info = 0;
  DMTK_int basis_size = s.basis().size();
  dsyev_('V', 'L', basis_size, hami.array(), basis_size, w, work, 3*basis_size, info);

  std::cout << info << std::endl;

  for(int i = 0; i < s.basis.size(); i++)
    std::cout << "Eigenvalue " << i << " " << w[i] << std::endl;

}


/////////////////////////////////////////////////////////////////////
// TIME EVOLUTION
/////////////////////////////////////////////////////////////////////
Vector<double>
exp_h(SystemED<double> &S)
{
  Vector<Vector<double> >aux(1);
  Vector<double> _e(1);
  Vector<double> a(1000), b(1000);
  Vector<double> d(1000), e(1000);
  int tridiag_size;
  Vector<double>aux_vector = S.gs;
  aux[0] = aux_vector;

  cout << ">>>>>>> LANCZOS <<<<<<<\n";
  tridiag_size = S.n_exp_lanczos() == -1 ? 10 : S.n_exp_lanczos();
  lanczos<double, SystemED<double>, Vector<double> >(S, aux, aux_vector, _e, 1, a, b, tridiag_size, 1.e-16, true, true, "vectors_test.dat");
  Matrix<double> z(1000,1000);
  z = I<double>();
  d(Range(1,tridiag_size)) = a(Range(1,tridiag_size));
  e(Range(2,tridiag_size+1)) = b(Range(1,tridiag_size));
  tqli(d, e, tridiag_size, z, true);

  size_t gs_index;
  for(int i = 1; i <= tridiag_size; i++){
    if(d[i] == _e[0]) {gs_index = i-1; break;}
  }

  Matrix<double> m_exp_h(tridiag_size, tridiag_size);
  m_exp_h = 0.;
  for(int i = 0; i < tridiag_size; i++)
    m_exp_h(i,i) = std::exp(-(d[i+1]-S.phase)*0.5*S.dt);

  Matrix<double> u(tridiag_size, tridiag_size);
  for(int col = 0; col < tridiag_size; col++)
   for(int row = 0; row < tridiag_size; row++)
     u(col,row) = z(col+1,row+1);

  Matrix<double> ut(u.rows(),u.cols());
  ut = ctranspose(u);
  Matrix<double> aux_m(u.cols(),m_exp_h.rows());
  aux_m = product(m_exp_h,ut);
  m_exp_h = product(u,aux_m);

  Vector<double> aux_gs(tridiag_size);
  aux_gs = m_exp_h.column(0);

  ifstream inputfile("vectors_test.dat",std::ios::in|std::ios::binary);
  aux_vector = 0.;
  for(int i = 0; i < tridiag_size; i++){
    aux[0].read(inputfile);
    aux_vector += aux_gs[i] * aux[0];
  }
  inputfile.close();
                               
  return aux_vector;
}

Vector<complex<double> >
exp_h(SystemED<complex<double> > &S)
{
  Vector<Vector<complex<double> > >aux(1);
  Vector<double> _e(1);
  Vector<double> a(1000), b(1000);
  Vector<double> d(1000), e(1000);
  int tridiag_size;
  Vector<complex<double> >aux_vector = S.gs;
  aux[0] = aux_vector;

  cout << ">>>>>>> LANCZOS <<<<<<<\n";
  tridiag_size = S.n_exp_lanczos() == -1 ? 10 : S.n_exp_lanczos();
  lanczos<complex<double>, SystemED<complex<double> >, Vector<complex<double> > >(S, aux, aux_vector, _e, 1, a, b, tridiag_size, 1.e-16, true, true, "vectors_test.dat");
  Matrix<double> z(1000,1000);
  z = I<double>();
  d(Range(1,tridiag_size)) = a(Range(1,tridiag_size));
  e(Range(2,tridiag_size+1)) = b(Range(1,tridiag_size));
  tqli(d, e, tridiag_size, z, true);

  size_t gs_index;
  for(int i = 1; i <= tridiag_size; i++){
    if(d[i] == _e[0]) {gs_index = i-1; break;}
  }

  Matrix<complex<double> > m_exp_h(tridiag_size, tridiag_size);
  m_exp_h = complex<double>(0);
  for(int i = 0; i < tridiag_size; i++){
    complex<double> coef = S.imaginary_time() ? complex<double>(-(d[i+1]-S.phase)*0.5*S.dt) : complex<double>(0.,-(d[i+1]-S.phase)*S.dt);
    m_exp_h(i,i) = std::exp(coef);
  }

  Matrix<complex<double> > u(tridiag_size, tridiag_size);
  for(int col = 0; col < tridiag_size; col++)
   for(int row = 0; row < tridiag_size; row++)
     u(col,row) = z(col+1,row+1);

  Matrix<complex<double> > ut(u.rows(),u.cols());
  ut = ctranspose(u);
  Matrix<complex<double> > aux_m(u.cols(),m_exp_h.rows());
  aux_m = product(m_exp_h,ut);
  m_exp_h = product(u,aux_m);

  Vector<complex<double> > aux_gs(tridiag_size);
  aux_gs = m_exp_h.column(0);

  ifstream inputfile("vectors_test.dat",std::ios::in|std::ios::binary);
  aux_vector = complex<double>(0);
  for(int i = 0; i < tridiag_size; i++){
    aux[0].read(inputfile);
    aux_vector += aux_gs[i] * aux[0];
  }
  inputfile.close();
                               
  return aux_vector;
}

Vector<double>
exp_h_rk(SystemED<double> &s)
{
  Vector<double> new_state, k_state;
  Vector<double> hphi, aux_hphi;
  Vector<double> &gs = s.gs;
  double dt = s.dt * 0.05;
  double my_time = s.my_time;
  double one_sixth(1./6.);
  double one_third(1./3.);
  double one_half(0.5);
  double c_energy(s.phase);

  new_state = gs;
  k_state = gs;

  for(int i = 1; i <= 10; i++){
//    hphi = product_default(s, gs);
    hphi = product(s, gs);
    hphi -= c_energy*gs;

    k_state = hphi;
    k_state *= -dt;

    new_state += one_sixth*k_state;

//    aux_hphi = product_default(s, k_state);
    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt;

    new_state += one_third*k_state;

//    aux_hphi = product_default(s, k_state);
    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt;

    new_state += one_third*k_state;

//    aux_hphi = product_default(s, k_state);
    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += aux_hphi;
    k_state *= -dt;

    new_state += one_sixth*k_state;
    gs = new_state;
    my_time += dt;
  }

  return new_state;
}

Vector<complex<double> >
exp_h_rk(SystemED<complex<double> > &s)
{
  Vector<complex<double> > new_state, k_state;
  Vector<complex<double> > hphi, aux_hphi;
  Vector<complex<double> > &gs = s.gs;
  double dt = s.imaginary_time() ? s.dt*0.05 : s.dt * 0.1;
  double my_time = s.my_time;
  complex<double> one_sixth(1./6.,0.);
  complex<double> one_third(1./3.,0.);
  complex<double> one_half(0.5,0.);
  complex<double> c_energy(s.phase,0.);

  new_state = gs;
  k_state = gs;
  complex<double> coef = s.imaginary_time() ? -dt : complex<double>(0,-dt);

  for(int i = 1; i <= 10; i++){
    hphi = product(s, gs);
    hphi -= c_energy*gs;

    k_state = hphi;
    k_state *= coef;
    new_state += one_sixth*k_state;

    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*(aux_hphi);
    k_state *= coef;

    new_state += one_third*k_state;

    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*(aux_hphi);
    k_state *= coef;

    new_state += one_third*k_state;

    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += aux_hphi;
    k_state *= coef;

    new_state += one_sixth*k_state;
    gs = new_state;
    my_time += dt;
  }

  return new_state;
}


template<class T>
Matrix<T>
SystemED<T>::density_matrix(int pos) const
{
  Matrix<T> dm(dim0(pos),dim0(pos));
  Matrix<T> coef(dim0(pos),basis.size());
  for(int i = 0; i < basis.size(); i++){
    long int state1 = state(basis(i),pos);
    long int state2 = basis(i);
    for(int bit = 0; bit < nbits0(pos); bit++){
      state2 = IBCLR(state2, bit+bit0(pos));
    }
    long int index = get_state_index(state2);
    coef(state1,index) += gs(i);
  }
  for(int i = 0; i < dm.rows(); i++){
    slice_iter<T> coli = coef.column(i);
    for(int j = 0; j < dm.rows(); j++){
      slice_iter<T> colj = coef.column(j);
      dm(i,j) = T(0);
      for(int k = 0; k < basis.size(); k++) dm(i,j) += coli(k)*conj(colj(k)); 
      cout << "DM(" << i << "," << j << ") = " << dm(i,j) << endl;
    }
  }


  return dm;
}

template<class T>
T
SystemED<T>::measure(const BasicOp<T> &ref_op) const
{
  BasicOp<T> op = operator()(ref_op);
  Vector<T> res = gs;
  res = product(op, gs);
  T x = 0;
  for(int i = 0; i < gs.size(); i++) x += std::conj(gs(i))*res(i);
  cout << ref_op.description() << " = " << x << endl;
  return x;
}

template<class T>
T
SystemED<T>::measure(const BasicOp<T> &_op1, const BasicOp<T> &_op2) const
{
  BasicOp<T> op1 = operator()(_op1);
  BasicOp<T> op2 = operator()(_op2);
  Vector<T> res = gs;
  res = product(op1, op2, gs);
  T x = 0;
  for(int i = 0; i < gs.size(); i++) x += std::conj(gs(i))*res(i);
  cout << _op1.description() << _op2.description() << " = " << x << endl;
  return x;
}

/////////////////////////////////////////////////////////////////////
} // namespace dmtk

#endif // __DMTK_SYSTEM_ED__
