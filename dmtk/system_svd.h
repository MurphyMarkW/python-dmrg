#ifndef __DMTK_SYSTEM_SVD_H__
#define __DMTK_SYSTEM_SVD_H__
//*******************************************************************
using namespace std;

namespace dmtk
{

template<class T>
class SystemSVD
{
  protected:
    BMatrix<T> rho_l;
    BMatrix<T> rho_r;
    Basis rho_l_basis;
    Basis rho_r_basis;
  public:

    SystemSVD() {};

    SystemSVD(const Hami<T> &_h, const Lattice& lattice, const char *the_name): System(_h,lattice,the_name) {}

    void write(std::ostream &s) const
    {
      s.write((const char *)&_numsweeps, sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        s.write((const char *)&_sweeps(0,i), sizeof(double));
        s.write((const char *)&_sweeps(1,i), sizeof(double));
      }
      s.write((const char *)&_in_warmup, sizeof(bool));
      s.write((const char *)&_sweep, sizeof(int));
      s.write((const char *)&dir, sizeof(size_t));
      s.write((const char *)&iter, sizeof(int));
//      s.write((const char *)&_ntargets, sizeof(int));
      s.write((const char *)&_nstates, sizeof(int));
      int n = _propagate_state.size();
      s.write((const char *)&n, sizeof(int));
      qnt.write(s);
      qn.write(s);
      gs.write(s);
      n = _target.size();
      s.write((const char *)&n, sizeof(int));
      for(int i = 0; i < n; i++){
        s.write((const char *)&_target_weight[i], sizeof(double));
        _target[i].write(s);
      }
      for(int i = 0; i < _propagate_state.size(); i++){
        _propagate_state[i].write(s);
      }
      rho_l.write(s);
      rho_l_basis.write(s);
      rho_r.write(s);
      rho_r_basis.write(s);
      rightblock.write(s);
      leftblock.write(s);
    }

    void read(std::istream &s)
    {
      s.read((char *)&_numsweeps, sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        s.read((char *)&_sweeps(0,i), sizeof(double));
        s.read((char *)&_sweeps(1,i), sizeof(double));
      }
      s.read((char *)&_in_warmup, sizeof(bool));
      s.read((char *)&_sweep, sizeof(int));
      s.read((char *)&dir, sizeof(size_t));
      s.read((char *)&iter, sizeof(int));
      s.read((char *)&_ntargets, sizeof(int));
      int n;
      s.read((char *)&n, sizeof(int));
      _propagate_state.resize(n);
      qnt.read(s);
      qn.read(s);
      gs.read(s);
      s.read((char *)&n, sizeof(int));
      _target.resize(n);
      _target_weight.resize(n);
      for(int i = 0; i < n; i++){
        s.read((char *)&_target_weight[i], sizeof(double));
        _target[i].read(s);
      }
      for(int i = 0; i < _propagate_state.size(); i++){
        _propagate_state[i].read(s);
      }
      rho_l.read(s);
      rho_r.read(s);
      rho_l_basis.read(s);
      rho_r_basis.read(s);
      rightblock.read(s);
      leftblock.read(s);
    }

};

/////////////////////////////////////////////////////////////////////////
// truncate:
// build and diagonalize density matrix. 
/////////////////////////////////////////////////////////////////////////
template<class T>
static Vector<double>
diagonal(BMatrix<T> &m)
{
  Vector<double> w(m.subspaces().dim());

  typename BMatrix<T>::iterator iter;
  for(iter = m.begin(); iter != m.end(); iter++){
    SubMatrix<T> &b = (*iter);

    Vector<double> d(b.rows());
    for(int i = 0; i < d.size(); i++) d(i) = real(b(i,i));
    w(b.row_range()) = d(Range(0,b.rows()-1));
    b = I<T>();
  }
  return w;
}

template<class T>
BMatrix<T>
SystemSVD<T>::density_matrix(int position)
{
  CTimer clock;
  clock.Start();
  signal_emit(SYSTEM_SIGNAL_BUILD_DM);
//******* Calculate density matrix (block by block)
  BMatrix<T> dm;

  Basis basis(_b1->basis(),_b2->basis());
  if(position == RIGHT) basis = Basis(_b3->basis(),_b4->basis());
  basis.reorder();

  if(verbose() == 3) 
  {
    for(int i = 0; i < basis.size(); i++)
      cout << "NEW BASIS " << i << " " << basis[i].qn().n() << " " << basis[i].qn().sz() << " " << basis[i].qn().kx1() << "/" << basis[i].qn().kx2() << " " << basis[i].qn().ky1() << "/" << basis[i].qn().ky2() << endl;
  }

  dm.resize(basis);

  cout << "Building Density Matrix in blocks \n";

  PackedBasis::iterator basis_iter;
  for(basis_iter = basis.subspace_begin(); basis_iter != basis.subspace_end(); basis_iter++){
    SubSpace s = *basis_iter;
    SubMatrix<T> block(s.qn(),s,s);
    dm.push_back(block);
  }

  for(int nt = 0; nt < _ntargets; nt++){
    if(_target_weight[nt] != 0.0) {
      VectorState<T> &s = _target[nt]; 
#ifdef USE_K
      VectorState<T> s2 = s;
      s2.randomize();
//      s2 *= T(0.000001);
      s2 *= T(_random_coef);
      s += s2; 
      gs = s;
//    if(_in_warmup && lattice().size() != size()) s.randomize();
#endif // USE_K
      T norm = product(s,s);
      cout << "Target vector " << nt << " " << _target_weight[nt] << " " << norm << endl;
      s /= sqrt(norm);
//      if(s.qn_mask() != QN::get_qn_mask()){
//        VectorState<T> s2 = s;
//        s2.set_qn_mask(s2.qn(),QN::get_qn_mask());
//        s2.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis()); 
//        dm += s2.density_matrix(position)*T(_target_weight[nt]);
//      } else 
      {
        dm += s.density_matrix(position)*T(_target_weight[nt]);
      }
    }
  }

  cout << "Lap: " << clock.LapTime().c_str() << endl;
  clock.Lap();

  T trace = T(0);
  typename BMatrix<T>::iterator biter;
  for(biter = dm.begin(); biter != dm.end(); biter++){
    SubMatrix<T> &block = (*biter);
/*
    for(int i = 0; i < block.size1(); i++)
      for(int j = 0; j < block.size2(); j++)
        cout << i << " " << j << " " << block(i,j) << endl;
    cout << "---------------------------\n"; 
*/
    trace += block.trace();
  }
  cout << "TRACE = " << trace << endl; 
  signal_emit(SYSTEM_SIGNAL_DM_READY);
  cout << "Density Matrix time: " << clock.TotalTime().c_str() << endl;

  return dm;
}


template<class T>
void
SystemSVD<T>::truncate(int position, int new_size, bool build_targets, bool diagonalize)
{
  int ls = lattice().ls();

  CTimer clock;
  clock.Start();

//******* Create target states
  if(build_targets) build_target_states(position);
  
//******* Measure control operators
  typename Hami<T>::iterator op_iter;
  for(op_iter = control_ops.begin(); op_iter != control_ops.end(); op_iter++){
    Term<T> &term = *op_iter;
    BasicOp<T> this_op = term[0].internals();
    bool info;
    if(dir == LEFT2RIGHT){
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(iter),NULL,info);
      } else if(iter == 1) {
        measure_operator(this_op.set_site(0),NULL,info);
        measure_operator(this_op.set_site(1),NULL,info);
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(ls-3),NULL,info);
        measure_operator(this_op.set_site(ls-2),NULL,info);
        measure_operator(this_op.set_site(ls-1),NULL,info);
      }
    } else {
      if(iter > 1 && iter < ls-3){
        measure_operator(this_op.set_site(ls-2-iter),NULL,info);
      } else if(iter == 1) {
        measure_operator(this_op.set_site(ls-1),NULL,info);
        measure_operator(this_op.set_site(ls-2),NULL,info);
        measure_operator(this_op.set_site(ls-3),NULL,info);
      } else if(iter == ls-3){
        measure_operator(this_op.set_site(1),NULL,info);
        measure_operator(this_op.set_site(0),NULL,info);
      }
    }
  }

#ifndef USE_PRODUCT_DEFAULT
  if(_rotate_terms){
    H<T> haux(h23);
    hint(MASK_BLOCK2|MASK_BLOCK3);
    VectorState<T> v23 = gs.condense(MASK_BLOCK2|MASK_BLOCK3);
    VectorState<T> res23(v23);
    res23 = T(0);
    product(h23, v23, res23, BLOCK2, T(1), false);
    cout << "LOCAL ENERGY " << (dir == LEFT2RIGHT ? iter : ls-2-iter) << " " << product(res23,v23)/product(v23,v23) << endl;
    h23 = haux;
  }
#endif // USE_PRODUCT_DEFAULT
//********************************************************************
//                   DENSITY MATRIX
//********************************************************************
  signal_emit(SYSTEM_SIGNAL_TRUNCATE);

  Basis basis(_b1->basis(),_b2->basis());
  const Block<T> *the_block = _b1;
  if(position == RIGHT) {
    basis = Basis(_b3->basis(),_b4->basis());
    the_block = _b4;
  }
  basis.reorder();

/////////////////////////////////////////////////////////////////////
  if(basis.dim() <= new_size){
// We use indentity matrix for the new operators since
// we don't want to rotate or truncate them

    BMatrix<T> aux_rho(basis);
    aux_rho.clear();
    typename PackedBasis::const_iterator biter;
    for(biter = basis.subspaces().begin(); biter != basis.subspaces().end(); biter++){
      SubMatrix<T> block(biter->qn(),*biter,*biter);
      block=(I<T>());
      aux_rho.push_back(block);
    }

    rho_basis = basis;
    rho = aux_rho;

    if(dm_ops.size() > 0){
      _target_weight[0] = 1.0;
      _target.resize(1);
      _target_weight.resize(1);
    }

    return;
  }
/////////////////////////////////////////////////////////////////////

  Vector<double> w(basis.dim()); // eigenvalues
  double entropy = 0.0;

  if(QN::get_qn_mask() != _grand_canonical){

//******* Calculate entropy
    size_t aux_mask = QN::get_qn_mask();
    QN::set_qn_mask(_grand_canonical);
    rho = density_matrix(position); 
    if(diagonalize){
      rho.diagonalize(w);
      for(int i = 0; i < w.size(); i++) cout << i << " " << w[i] << endl;
      BMatrix<T> u, v;
      rho.svd(w, u, v);
    }
    else 
      w = diagonal(rho);
    
    for(int i = 0; i < w.size(); i++){
      if(w[i] >= 1.e-10) {
//        cout << "WARNING: w[i] <= 0 " << w[i] << endl;
        entropy -= w[i]*log(w[i]);
      }
      if(_verbose > 0) 
        cout << "RHO EIGENVALUE " << i << " " << w[i] << endl;
    }
    cout << "ITER = " << iter << " ENTROPY = " << entropy << endl;
    QN::set_qn_mask(aux_mask);
  }

//******* Diagonalize Density Matrix in blocks
  cout << "Diagonalizing Density Matrix in blocks\n";
// Test with energy eigenvectors
  Vector<size_t> idx(basis.dim());
  rho_w.resize(basis.dim());

  if(use_energy_rg()){
    if(position == LEFT)
      rho = h12;
    else
      rho = h34;
  } else {
    rho = density_matrix(position);

    BMatrix<T> aux_rho;
    aux_rho = rho;
    if(diagonalize) {
      aux_rho.diagonalize(w);
      for(int i = 0; i < w.size(); i++) cout << i << " " << w[i] << endl;
      BMatrix<T> u, v;
      aux_rho.svd(w, u, v);
    } else
      w = diagonal(aux_rho);
    indexx<double, Vector<double> >(basis.dim(), w, idx, use_energy_rg());
    double wtotal = 1.;
    for(int i = 0; i < std::min(int(new_size),int(basis.dim())); i++){
      wtotal -= w(idx(i));
      if(verbose())
        cout << "REAL RHO EIGENVALUE " << i << " " << w[idx(i)] << endl;
    }
    for(int i = 0; i < basis.dim(); i++) rho_w(i) = w(idx(i));

    cout << "REAL TRUNCATION ERROR = " << iter << " " << basis.dim() << " " << wtotal << endl;

    if(dm_boundary_ops.size() != 0 && QN::get_qn_mask() == _grand_canonical){
      BMatrix<T> rho_copy; rho_copy = rho;
      if(diagonalize) 
        rho_copy.diagonalize(w);
      else
        w = diagonal(rho_copy);
      double entropy = 0.0;
      for(int i = 0; i < w.size(); i++){
        if(w[i] >= 1.e-10) {
//          cout << "WARNING: w[i] <= 0 " << w[i] << endl;
          entropy -= w[i]*log(w[i]);
        }
        if(_verbose > 0) 
          cout << "RHO EIGENVALUE " << i << " " << w[i] << endl;
      }
      cout << "ITER = " << iter << " ENTROPY = " << entropy << endl;
    }

//  Apply boundary operators to calculate delta_rho
    delta_rho();
  }


// Now we diagonalize the final d.m.
  if(diagonalize) {
    rho.diagonalize(w);
    for(int i = 0; i < w.size(); i++) cout << i << " " << w[i] << endl;
    BMatrix<T> u, v;
    rho.svd(w, u, v);
  } else
    w = diagonal(rho);

  cout << "Lap: " << clock.LapTime().c_str() << endl;

//******* Calculate entropy

  if(dm_boundary_ops.size() == 0 && QN::get_qn_mask() == _grand_canonical){
    double entropy = 0.0;
    for(int i = 0; i < w.size(); i++){
      if(w[i] >= 1.e-10) {
//        cout << "WARNING: w[i] <= 0 " << w[i] << endl;
        entropy -= w[i]*log(w[i]);
      }
      if(_verbose > 0) 
        cout << "RHO EIGENVALUE " << i << " " << w[i] << endl;
    }
    cout << "ITER = " << iter << " ENTROPY = " << entropy << endl;
  }
//*******************************************************************
// Entanglement spectrum
// Make sure that you are not using boundary operators 
//*******************************************************************
  if(_print_dm_spectrum){
    typename BMatrix<T>::iterator rho_iter;
    int i = 0;
    for(rho_iter = rho.begin(); rho_iter != rho.end(); rho_iter++){
      SubMatrix<T> &sm = (*rho_iter);
      cout << "RHO QN " << sm.qn() << endl;
//      cout << "===============================================" << endl;
      for(int j = 0; j < sm.size1(); j++){
//        cout << "RHO QN " << sm.qn().n() << " " << sm.qn().sz() << " " << sm.qn().kx() << " " << sm.qn().ky();
        cout << " " << w[i+j] << endl;
      }
      i += sm.size1();
    }
  }

//******* Reorder eigenvalues and eigenvectors, and truncate
  cout << "Reordering eigenvalues and eigenvectors, and truncating\n";

/*
  if(_target_subspaces && _in_warmup && lattice().size() != size()){
    for(int i = 0; i < basis.dim(); i++){
      if(basis[i].qn().n() == qnt.n()/2) w[i] = 1.;
      if(basis[i].qn().n() == ls/2 && ls < qnt.n()) w[i] = 1.;
    } 
  }
*/

  indexx<double, Vector<double> >(basis.dim(), w, idx, use_energy_rg());

  Vector<size_t> new_idx(basis.dim());
  for(int i = 0; i < basis.dim(); i++) new_idx(idx(i)) = i;

  if(_use_error){
    double error = 1.;
    for(int i = 0; i < basis.dim(); i++){
      error -= w(idx(i));
      if(error <= _error && i >= new_size){
//        new_size = std::min(new_size, (int)((i+1)*1.1));
        new_size = (int)((i+1)*1.1);
        if(_error_max_size != -1) 
          new_size = std::min(new_size, _error_max_size);
        cout << "ERROR = " << _error << " " << error << endl;
        cout << "NEW SIZE = " << i+1 << " " << new_size << endl;
        break;
      }
    }
  }

  Vector<double> vaux = w;
  rho_basis = basis;
  int new_dim = 0;

  typename BMatrix<T>::iterator biter = rho.begin();
//  if(_target_subspaces && _in_warmup && lattice().size() != size())
  if(_target_subspaces && _in_warmup)
//  if(_target_subspaces)
  {
    Vector<double> bweight(rho.size());
    Vector<size_t> bidx(rho.size());
    Vector<SubMatrix<T>*> bvector(rho.size());
    size_t ib = 0;
    while(biter != rho.end()){
      SubMatrix<T> &block = (*biter);
//      if(block.qn().n() > qnt.n()) { biter++; continue; }
//      if(block.qn().n() > qnt.n()/2*3) { biter++; continue; }
//      if(qnt.n()-block.qn().n() > h.lattice().size()-the_block->lattice().size()) { biter++; continue; }
//      if(block.qn().n() == qnt.n() && block.qn().kx() != qnt.kx()) { biter++; continue; }
      bvector[ib] = &block;
      Range col_range = block.col_range();
      bweight[ib] = 0;
      for(size_t col = col_range.begin(); col <= col_range.end(); col++)
        bweight[ib] += w(col);
      ib++;
      biter++;
    }
    int bsize = ib;
    indexx<double, Vector<double> >(bsize, bweight, bidx, use_energy_rg());
    int ngroup = std::max(_nsub,int(rho.size()/bsize));
    int nmax = 0;
    for(ib = 0; ib < bsize; ib++){
      SubMatrix<T> &block(*bvector(bidx(ib)));
      Range col_range = block.col_range();
      nmax += std::min(int(col_range.size()),ngroup);
// FIXME      if(nmax > new_size) break;   
    }
//    bsize = std::min(ib,rho.size());

    BMatrix<T> new_rho;
    Vector<double> waux(rho.size());
    for(biter = rho.begin(); biter != rho.end(); biter++){
      SubMatrix<T> &block = (*biter);
      bool found = false;
      for(ib = 0; ib < bsize; ib++){
        SubMatrix<T> *_block = bvector(bidx(ib));
        if(_block == &block) { found = true; break; }
      }
      if(!found) continue; 

// cout << "NEW BLOCK " << " " << block.qn().n() << " " << block.qn().sz() << " " << block.qn().kx() << endl;
      SubMatrix<T> aux_block(block); // copy of the original
      Range col_range = block.col_range();
      int first_col = new_dim;
      Vector<double> waux(col_range.size());
      waux(Range(0,col_range.size()-1)) = vaux(col_range);
      indexx<double, Vector<double> >(col_range.size(), waux, idx, use_energy_rg());

      int maxcol = std::min(size_t(ngroup), col_range.size());
 
      for(size_t col = 0; col < maxcol; col++){
        block.column(col) = aux_block.column(idx(col));
        w(new_dim) = waux(idx(col));
        rho_basis(new_dim) = basis(idx(col)+col_range.begin());
        new_dim++;
      }
      if(first_col != new_dim){
        block.resize(Range(first_col,new_dim-1),block.row_range());
        new_rho.push_back(*biter);
      }
      if(new_dim >= new_size) break;
    }
  }
  else
  {
    while(biter != rho.end()){
      SubMatrix<T> &block = (*biter);
      SubMatrix<T> aux_block(block); // copy of the original
      Range col_range = block.col_range();
      int new_col = 0;
      int first_col = new_dim;
      for(size_t col = col_range.begin(); col <= col_range.end(); col++){
        if(new_idx(col) < new_size){
// cout << "NEW STATE " << col << " " << basis(col).qn().n() << " " << basis(col).qn().sz() << endl;
           block.column(new_col) = aux_block.column(col-col_range.begin());
           w(new_dim) = vaux(col);
           rho_basis(new_dim) = basis(col);
           new_col++;
           new_dim++;
        }
      }
      if(first_col == new_dim){
#ifndef USE_K
        biter = rho.erase(biter);
#else // USE_K
//        if(!_in_warmup && iter <= lattice().size()/2){
        if((_in_warmup && iter <= lattice().size()/2) || _target_subspaces){
          int add_states = std::min(block.size1(),size_t(2));
          block.resize(add_states,block.size2());
          rho_basis(new_dim) = basis(col_range.begin());
          if(add_states == 2) rho_basis(new_dim+1) = basis(col_range.begin()+1);
          new_dim += add_states;
          biter++;
        } else {
          biter = rho.erase(biter);
        }
#endif
      }else{
        block.resize(Range(first_col,new_dim-1),block.row_range());
        biter++;
      }
//      if(new_dim >= new_size) break;
    }
  }

  for(biter = rho.begin(); biter != rho.end(); biter++) cout << "ADDING BLOCK " << (*biter).qn() << endl;
  cout << "NEW DIM " << new_dim << endl;
  rho_basis.resize(new_dim); // truncate basis
  rho.repack(rho_basis);

//******* Check precision after truncation

  double xweight = 0;
  for(int i = 0; i < new_dim; i++){
//    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
    xweight += w(i);
  }

  cout << "TOTAL WEIGHT = " << xweight << endl;

//  xweight = 0;
//  for(int i = new_size; i < w.size(); i++){
////    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
//    xweight += w(i);
//  }
//
//  cout << "DISCARDED WEIGHT = " << xweight << endl;

  precision = 1.;
  for(int i = 0; i < new_dim; i++){
//    cout << i << " " << w(i) << " " << vaux(idx(i)) << endl;
    precision -= w(i);
  }

//////////////////////////////////////////////////////////////////
// dm_ops
//////////////////////////////////////////////////////////////////
  if(dm_ops.size() > 0){
    _target_weight[0] = 1.0;
    _target.resize(1);
    _target_weight.resize(1);
  }
//////////////////////////////////////////////////////////////////

  cout << "-------------------------------------------\n";
  cout << "Truncation error = " << precision << endl;
  cout << "Truncation time: " << clock.TotalTime().c_str() << endl;
  cout << "-------------------------------------------\n";

  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;

  outputfile.precision(10);

  outputfile << "Truncation error = " << precision << endl;
  outputfile << "Truncation time: " << clock.TotalTime().c_str() << endl;
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// rotate:
// build new block and operators, truncating with the density matrix.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
SystemSVD<T>::rotate(int position, Block<T>& b)
{
  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
  }

  CTimer clock;
  clock.Start();

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

// We build new lattice for the new block 

  if(verbose() > 0){
    for(int i = 0; i < rho_basis.size(); i++)
      cout << "RHO BASIS " << i << " " << rho_basis[i].qn().n() << " " << rho_basis[i].qn().sz() << " " << rho_basis[i].qn().kx1() << endl;
  }

  b = Block<T>(rho_basis, b1.n_orbitals()+b2.n_orbitals());
  b.set_single_site(false);
  Lattice l, l2;

  if(position == LEFT){
    l = _b1->lattice();
    l2 = _b2->lattice();
  } else {
    l = _b4->lattice();
    l2 = _b3->lattice();
  }
  for(int is = 0; is < l2.size(); is++){
    l.push_back(Site());
    l.nx() += 1;
  }
  b.set_lattice(l);

  b = T(0);
  Basis basis(b1.basis(),b2.basis());
  basis.reorder();

/////////////////////////////////////////////////////////////
  if(verbose() == 3)
  {
    Vector<const Block<T>*> _b(4);
    _b[0] = _b1; 
    _b[1] = _b2; 
    _b[2] = _b3; 
    _b[3] = _b4; 
    for(int i = 0; i < 4; i++){
      cout << "BLOCK " << i+1 << endl;
      typename Block<T>::const_iterator biter;
      for(biter = _b[i]->begin(); biter != _b[i]->end(); biter++){
        cout << "OPERATOR " << (*biter).name() << "(" << (*biter).site() << "," << (*biter).internal_site() << ")" << endl;
      }
    }
    cout << "-------------\n";
  }
/////////////////////////////////////////////////////////////
//********************************************************************
// Test using identity
//********************************************************************
// TODO: comment this lines 
/*
  typename BMatrix<T>::iterator biter = rho.begin();

  for(biter = rho.begin(); biter != rho.end(); biter++){
    SubMatrix<T> &block(*biter);
    block = I<T>(); 
  }
*/
// ***************************************************************************
// Hamiltonian terms for interactions between blocks
// ***************************************************************************
  signal_emit(SYSTEM_SIGNAL_ROTATE);
  if(_rotate_terms) rotate_terms(position, b, basis, rho_basis, NULL);
// ***************************************************************************
// Local Hamiltonian
// ***************************************************************************
  cout << "NEW BLOCK: HAMILTONIAN TERMS: " << this->h.name() << endl;
  if(_rotate_terms) rotate_hami(position, b, basis, rho_basis, NULL); 
// ***************************************************************************
// Other Hamiltonians 
// ***************************************************************************
  if(_rotate_terms){
    typename std::vector<Hami<T> >::iterator hiter;
    for(hiter = hs.begin(); hiter != hs.end(); hiter++){
      rotate_terms(position, b, basis, rho_basis, &(*hiter));
      rotate_hami(position, b, basis, rho_basis, &(*hiter)); 
    }
  }
// ***************************************************************************
// Correlations
// ***************************************************************************
  if(_rotate_corr){
    rotate_corr(position, b, basis, rho_basis, NULL); 
    rotate_corr(position, b, basis, rho_basis, &ops); 
    rotate_corr(position, b, basis, rho_basis, &dm_ops); 
  }
// ***************************************************************************
  rho_basis = basis;
  cout << "-------------------------------------------\n";
  cout << "Rotation time: " << clock.TotalTime().c_str() << endl;
  cout << "-------------------------------------------\n";

  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;

  outputfile << "Rotation time: " << clock.TotalTime().c_str() << endl;
  outputfile.close();

}

/////////////////////////////////////////////////////////////////////////
// rotate_tems:
// rotate hamiltonian terms 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_terms(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami)
{
  if(!this_hami) this_hami = &this->h;

  int _mask = MASK_BLOCK1 | MASK_BLOCK2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;
  size_t the_block = BLOCK1;
  size_t site_block = BLOCK2;

  if(position == RIGHT){
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
    pos1 = BLOCK3;
    pos2 = BLOCK4;
    site_block = BLOCK3;
    the_block = BLOCK4;
  }

  CTimer clock;
  clock.Start();

//********************************************************************
//                   NEW BLOCK
//********************************************************************
//******* Change of basis

  cout << "NEW BLOCK: SINGLE OPERATORS\n";

// ***************************************************************************
// We add to the outer block the operators necessary to build the interaction
// between blocks. We take them from the single-site blocks and we add them 
// to the outer block. 
// We also add local terms of the hamiltonian on site block
// ***************************************************************************

  typename Hami<T>::const_iterator titer;
  for(titer = this_hami->begin(); titer != this_hami->end(); titer++){
    const Term<T>& t = (*titer);
//    if(T(t.coef()) == T(0)) continue;
    if(t.size() == 1 && t[0].name() == "I") continue;
  
    bool found = false; 
    const BasicOp<T> *_op = NULL;

    typename Term<T>::const_iterator oiter;
    int bmask = 0;

    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(oiter->is_hami()) {
        rotate_hami(position, b, basis, rho_basis, &(oiter->hami()));
        rotate_terms(position, b, basis, rho_basis, &(oiter->hami()));
      } else {
        bmask |= mask(block(oiter->site()));
      }
    }
    if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1) && bmask != MASK_BLOCK2)
       found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4) && bmask != MASK_BLOCK3)
       found = true;

    if(found){
      BasicOp<T> real_op;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(block(oiter->site()) == site_block && !oiter->is_hami()){
          _op = operator()(*oiter);
          if(!_op) cout << "WARNING : Operator " << (*oiter).description() << " not found\n";
          real_op = oiter->internals();
          real_op.dqn = _op->dqn;
          break;
        }
      }
      bool add = b.contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        cout << "ADDING NEW OPERATOR " << real_op.description() << endl;
        clock.Lap();

        new_operator(real_op, *_op, rho, basis, 1-position);
        b.push_back(real_op);
        cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// ***************************************************************************
// Composite operators for interactions 
// ***************************************************************************
// 1- We look for terms that involve pieces in one block and pieces
//    in another block (not a site block) and we rotate them 
//    We also add TERM_LOCAL and TERM_EXTERN that are in the block

  for(titer = this_hami->begin(); titer != this_hami->end(); titer++){
    const Term<T>& t = (*titer);
    bool found = false;
    if((t.type() == TERM_PRODUCT && t.size() >= 2) || t.type() == TERM_LOCAL || t.type() == TERM_EXTERN){
      Term<T> aux_term;
      BasicOp<T> top2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        if(block(top.site()) == the_block){
          aux_term *= top;
          found = true;
        }
        bmask |= mask(block(top.site()));
      }
      bool doit = false;
      if(found && position == LEFT && (bmask & MASK_BLOCK1) && ((!(bmask & MASK_BLOCK2) && bmask != MASK_BLOCK1) || (bmask == MASK_BLOCK1 && (t.type() == TERM_LOCAL || t.type() == TERM_EXTERN))))
         doit = true;
      if(found && position == RIGHT && (bmask & MASK_BLOCK4) && ((!(bmask & MASK_BLOCK3) && bmask != MASK_BLOCK4) || (bmask == MASK_BLOCK4 && (t.type() == TERM_LOCAL || t.type() == TERM_EXTERN))))
         doit = true;

      if(doit){ // we found a piece of a composite operator
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0].internals();
        else
          top1 = aux_term; // the piece contains more than one operator

        const BasicOp<T>* op1 = operator()(top1);

        if(!op1) {
          cout << "ERROR 1: Operator " << top1.description() << " not found\n";
          continue;
        }

        BasicOp<T> new_op(top1.internals()); 
        new_op.dqn = op1->dqn;
        bool add = b.contains(new_op);

        if(!add && op1){
          new_op.resize(rho_basis);
          clock.Lap();
  
          cout << "NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, rho, basis, position);
          b.push_back(new_op);
          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }

// 2- We look for terms that involve pieces in one block and pieces
//    in a site block and we put the pieces together

  for(titer = this_hami->begin(); titer != this_hami->end(); titer++){
    const Term<T>& t = (*titer);
    bool found1 = false;
    bool found2 = false;
    if(t.size() >= 2){
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        bmask |= mask(block(top.site()));
        
        if(block(top.site()) == pos1){
          aux_term1 *= top;
          found1 = true;
        }
        if(block(top.site()) == pos2){
          aux_term2 *= top;
          found2 = true;
        }
      }
//    if bmask == _mask, this term will go into the Hamiltonian
//    unless it is TERM_EXTERN
      if(((bmask != _mask) || (bmask == _mask && t.type() == TERM_EXTERN)) && found1 && found2){ // we have a new composite operator
        BasicOp<T> top1(aux_term1);
        if(aux_term1.size() == 1)
          top1 = aux_term1[0].internals();
        BasicOp<T> top2(aux_term2);
        if(aux_term2.size() == 1)
          top2 = aux_term2[0].internals();
        const BasicOp<T>* op1 = operator()(top1);
        const BasicOp<T>* op2 = operator()(top2);

        if(!op1) {
          cout << "ERROR 2: Operator " << top1.description() << " in term " << t.description() << " not found\n";
          continue;
        }
        if(!op2) {
          cout << "ERROR 2: Operator " << top2.description() << " in term " << t.description() << " not found\n";
          continue;
        }

        Term<T> new_term; 
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i].internals();
          if(block(top.site()) == pos1 || block(top.site()) == pos2){
            new_term *= top;
          }
        }
        BasicOp<T> new_op(new_term); 
        new_op.dqn = op1->dqn + op2->dqn;
        bool add = b.contains(new_op);

        if(!add && op1 && op2){
          new_op.resize(rho_basis);
          clock.Lap();
  
          cout << "ADDING NEW OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
          b.push_back(new_op);
          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }
}

// ***************************************************************************
// Measurement operators 
// ***************************************************************************
/////////////////////////////////////////////////////////////////////////
// rotate_corr:
// rotate measurement operators 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_corr(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami)
{
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;
  size_t the_block = BLOCK1;
  size_t site_block = BLOCK2;

  if(position == RIGHT){
    pos1 = BLOCK3;
    pos2 = BLOCK4;
    site_block = BLOCK3;
    the_block = BLOCK4;
  }

  CTimer clock;
  clock.Start();


  const Hami<T> *ph = &corr;
  if(this_hami) ph = this_hami;
  const Hami<T> &_h = *ph;

// 1a-Look for operators in S.ops involving the other blocks and site block 
// We also rotate individual operators

  typename Hami<T>::const_iterator titer;
  if(this_hami == &ops){
    for(titer = ops.begin(); titer != ops.end(); titer++){
      const Term<T>& t = (*titer);

      if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
 
      bool found = false; 
      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(!oiter->is_hami()) {
          bmask |= mask(block(oiter->site()));
        } else {
          continue;
        }
      }
      if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1)) found = true;
      if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4)) found = true;
  
      if(found){
        Term<T> aux_term;
        BasicOp<T> top2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i].internals();
          if(block(top.site()) == site_block){
            aux_term *= top;
          }
        }
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0].internals();
        else
          top1 = aux_term; // the piece contains more than one operator
  
        const BasicOp<T>* _op = operator()(top1);
  
        if(!_op) {
          cout << "ERROR 1a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> real_op(top1.internals()); 
        real_op.dqn = _op->dqn;
        bool add = b.contains(real_op);

        if(!add && _op){ // If we haven't added it yet, we do it now
          real_op.resize(rho_basis);
  
          cout << "ADDING NEW MEAS. OPERATOR " << real_op.description() << endl;
          clock.Lap();
  
          new_operator(real_op, *_op, rho, basis, 1-position);
          b.push_back(real_op);
          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  } 

// 1b-Look for operators in S.corr involving the other blocks and site block 
// We also rotate individual operators

  for(titer = _h.begin(); titer != _h.end(); titer++){
    const Term<T>& t = (*titer);
//    if(T(t.coef()) == T(0)) continue;
 
/* 
    if(t.type() == TERM_PRODUCT && t.size() == 2 && t[1].is_hami()){
      rotate_additive(position, b, basis, rho_basis, t); 
    }
*/
 
    bool found = false; 

    typename Term<T>::const_iterator oiter;
    int bmask = 0;
    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(!oiter->is_hami()) {
        bmask |= mask(block(oiter->site()));
      } else {
        rotate_terms(position, b, basis, rho_basis, &(oiter->hami()));
        rotate_hami(position, b, basis, rho_basis, &(oiter->hami()));
        continue;
      }
    }
    if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
    if(position == LEFT && (bmask & MASK_BLOCK2) && !(bmask & MASK_BLOCK1)) found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK3) && !(bmask & MASK_BLOCK4)) found = true;

    if(found){
      Term<T> aux_term;
      BasicOp<T> top2;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        if(block(top.site()) == site_block){
          aux_term *= top;
        }
      }
      if(aux_term.size() == t.size()) aux_term.coef() = t.coef();

      BasicOp<T> top1;
      if(aux_term.size() == 1) // the piece contains a single operator
        top1 = aux_term[0].internals();
      else
        top1 = aux_term; // the piece contains more than one operator

      const BasicOp<T>* _op = operator()(top1);

      if(!_op) {
        cout << "ERROR 1b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> real_op(top1.internals()); 
      real_op.dqn = _op->dqn;
      bool add = b.contains(real_op);

      if(!add && _op){ // If we haven't added it yet, we do it now
        real_op.resize(rho_basis);

        cout << "ADDING NEW MEAS. OPERATOR " << real_op.description() << endl;
        clock.Lap();

        new_operator(real_op, *_op, rho, basis, 1-position);
        b.push_back(real_op);
        cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// 2a- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  if(!this_hami){
    for(titer = ops.begin(); titer != ops.end(); titer++){
      const Term<T>& t = (*titer);
      bool found = false;
  
      if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;

      typename Term<T>::const_iterator oiter;
      int bmask = 0;
      for(oiter = t.begin(); oiter != t.end(); oiter++){
        if(!oiter->is_hami()) {
          bmask |= mask(block(oiter->site()));
        } else {
          continue;
        }
      }
      if(position == LEFT && (bmask & MASK_BLOCK1) && !(bmask & MASK_BLOCK2)) found = true;
      if(position == RIGHT && (bmask & MASK_BLOCK4) && !(bmask & MASK_BLOCK3)) found = true;
  
      if(found){ // we found a piece of a composite operator
        Term<T> aux_term;
        BasicOp<T> top2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i].internals();
          if(block(top.site()) == the_block){
            aux_term *= top;
          }
        }
        if(aux_term.size() == t.size()) aux_term.coef() = t.coef();
  
        BasicOp<T> top1;
        if(aux_term.size() == 1) // the piece contains a single operator
          top1 = aux_term[0].internals();
        else
          top1 = aux_term; // the piece contains more than one operator
  
        const BasicOp<T>* _op = operator()(top1);
  
        if(!_op) {
          cout << "ERROR 2a: Operator " << top1.description() << " not found\n";
          continue;
        }
  
        BasicOp<T> new_op(top1.internals()); 
        new_op.dqn = _op->dqn;
        bool add = b.contains(new_op);
  
        if(!add && _op){
          new_op.resize(rho_basis);
          clock.Lap();
  
          cout << "NEW MEAS. OPERATOR " << new_op.description() << endl;
          new_op = T(0);
          new_operator(new_op, *_op, rho, basis, position);
          b.push_back(new_op);
          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }

// 2b- We look for terms that involve pieces in old block and pieces
// in another block (not site blocks) and we rotate them 
// We also rotate individual operators

  for(titer = _h.begin(); titer != _h.end(); titer++){
    const Term<T>& t = (*titer);
    if(t.size() == 1 && (!_store_products && !t[0].is_hami())) continue;
    bool found = false;

    typename Term<T>::const_iterator oiter;
    int bmask = 0;
    for(oiter = t.begin(); oiter != t.end(); oiter++){
      if(!oiter->is_hami()) {
        bmask |= mask(block(oiter->site()));
      } else {
        continue;
      }
    }
    if(position == LEFT && (bmask & MASK_BLOCK1) && !(bmask & MASK_BLOCK2)) found = true;
    if(position == RIGHT && (bmask & MASK_BLOCK4) && !(bmask & MASK_BLOCK3)) found = true;

    if(found){ // we found a piece of a composite operator
      Term<T> aux_term;
      BasicOp<T> top2;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        if(block(top.site()) == the_block){
          aux_term *= top;
        }
      }

      BasicOp<T> top1;
      if(aux_term.size() == 1) // the piece contains a single operator
        top1 = aux_term[0].internals();
      else
        top1 = aux_term; // the piece contains more than one operator

      const BasicOp<T>* _op = operator()(top1);

      if(!_op) {
        cout << "ERROR 2b: Operator " << top1.description() << " not found\n";
        continue;
      }

      BasicOp<T> new_op(top1.internals()); 
      new_op.dqn = _op->dqn;
      bool add = b.contains(new_op);

      if(!add && _op){
        new_op.resize(rho_basis);
        clock.Lap();

        cout << "NEW MEAS. OPERATOR " << new_op.description() << endl;
        new_op = T(0);
        new_operator(new_op, *_op, rho, basis, position);
        b.push_back(new_op);
        cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
  }

// 3- Composite operators for correlations
  if(_store_products){
    for(titer = _h.begin(); titer != _h.end(); titer++){
      const Term<T>& t = (*titer);
      bool found1 = false;
      bool found2 = false;
      if(t.type() == TERM_PRODUCT && t.size() >= 2){
        Term<T> aux_term1, aux_term2;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> top = t[i].internals();
          if(block(top.site()) == pos1){
            aux_term1 *= top;
            found1 = true;
          }
          if(block(top.site()) == pos2){
            aux_term2 *= top;
            found2 = true;
          }
        } 
        if(found1 && found2){ // we have a new composite operator
          BasicOp<T> top1, top2;
          if(aux_term1.size() == 1) // the piece contains a single operator
            top1 = aux_term1[0].internals();
          else
            top1 = aux_term1; // the piece contains more than one operator
          if(aux_term2.size() == 1) // the piece contains a single operator
            top2 = aux_term2[0].internals();
          else
            top2 = aux_term2; // the piece contains more than one operator
          const BasicOp<T>* op1 = operator()(top1);
          const BasicOp<T>* op2 = operator()(top2);
  
          if(!op1) {
            cout << "ERROR 3: Operator " << top1.description() << " not found\n";
            continue;
          }
          if(!op2) {
            cout << "ERROR 3: Operator " << top2.description() << " not found\n";
            continue;
          }
 
          Term<T> new_term; 
          for(int i = 0; i < t.size(); i++){
            BasicOp<T> top = t[i].internals();
            if(block(top.site()) == pos1 || block(top.site()) == pos2){
              new_term *= top;
            }
          }

          BasicOp<T> new_op(new_term); 
          new_op.dqn = op1->dqn + op2->dqn;
          bool add = b.contains(new_op);
  
          if(!add && op1 && op2){
            new_op.resize(rho_basis);
            clock.Lap();
    
            cout << "ADDING NEW CORR. TERM " << new_op.description() << endl;
            new_op = T(0);
            new_operator(new_op, *op1, *op2, pos1, pos2, rho, basis, T(1));
  
            b.push_back(new_op);
            cout << "Lap: " << clock.LapTime().c_str() << endl;
          }
        }
      }
    }

  }

}

/////////////////////////////////////////////////////////////////////////
// rotate_hami:
// rotate hamiltonian terms 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_hami(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami)
{
  if(!this_hami) this_hami = &this->h;

#ifndef USE_PRODUCT_DEFAULT

  if(this_hami == &this->h && !_use_hamis){
    BasicOp<T> _h = BasicOp<T>(this->h);
//    _h.dqn.kx2() = qnt.kx2();
//    _h.dqn.ky2() = qnt.ky2();
    _h.resize(rho_basis);
    if(position == LEFT){
      new_operator(_h, h12, rho, basis);
      b.push_back(_h);
    } else { // position == RIGHT
      new_operator(_h, h34, rho, basis);
      b.push_back(_h);
    }

    return;
  }

#endif // USE_PRODUCT_DEFAULT

  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;
  int _mask = MASK_BLOCK1 | MASK_BLOCK2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
    pos1 = BLOCK3;
    pos2 = BLOCK4;
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

  CTimer clock;
  clock.Start();

  cout << "HAMILTONIAN TERMS: " << this_hami->name() << endl;

  BasicOp<T> _h = BasicOp<T>(*this_hami);
//  _h.dqn.kx2() = qnt.kx2();
//  _h.dqn.ky2() = qnt.ky2();
  _h.resize(rho_basis);

  bool found2 = false;
  bool found1 = false;
  const BasicOp<T> *local_op2 = b2(BasicOp<T>(*this_hami));
  if(local_op2 && local_op2->name() == this_hami->name()){ 
    cout << "OPERATOR " << this_hami->name() << "2" << endl;
    clock.Lap();
    found2 = true;
    new_operator(_h, *local_op2, rho, basis, RIGHT);
    cout << "Lap: " << clock.LapTime().c_str() << endl;
  } else {
#ifdef WITH_WARNINGS
    if(this_hami->name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
    else
      cout << "You can ignore the warning\n";
#endif
  }

  const BasicOp<T> *local_op1 = b1(BasicOp<T>(*this_hami));
  if(local_op1 && local_op1->name() == this_hami->name()){ 
    cout << "OPERATOR " << this_hami->name() << "1" << endl;
    clock.Lap();
    found1 = true;
    new_operator(_h, *local_op1, rho, basis, LEFT);
    cout << "Lap: " << clock.LapTime().c_str() << endl;
  } else {
#ifdef WITH_WARNINGS
    if(this_hami->name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
    else
      cout << "You can ignore the warning\n";
#endif
  }

  typename Hami<T>::const_iterator hiter;
  for(hiter = this_hami->begin(); hiter != this_hami->end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TERM_EXTERN) continue;
    if(t.type() == TERM_LOCAL) continue;
    if(t.size() == 1 && t[0].name() == "I") continue;
    if(t.size() == 1 && t[0].name() != this_hami->name()) {
      const BasicOp<T>&top = t[0];

/*
      // Rotate Hamiltonians terms in the original Hamiltonian
      // This is now done inside rotate_terms
      if(top.is_hami()){
        rotate_hami(position, b, basis, rho_basis, &(top.hami()));
        continue;
      }
*/
      if(top.is_hami()) continue;

      size_t ib = block(top.site());
      if((ib == BLOCK2 && position == LEFT) || (ib == BLOCK3 && position == RIGHT)) {
        const BasicOp<T>* op= operator()(top);
        bool calc_hc = this_hami->use_hc();
        if (op->is_diagonal()) calc_hc = false;
                                                                                
        if(!calc_hc)
          cout << "TERM " << t.name(true) << endl;
        else
          cout << "TERM " << t.name(true) << " + h.c." << endl;
                                                                                
        clock.Lap();
        new_operator(_h, *op, rho, basis, 1-position, T(t.coef()));
        cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
//if we don't have a site Hamiltonian we try to find the missing term
      if((ib == BLOCK1 && position == LEFT && !found1) || (ib == BLOCK4 && position == RIGHT && !found2) || (position == LEFT && ib == BLOCK1 && b1.lattice().size() == 1) || (position == RIGHT && ib == BLOCK4 && b2.lattice().size() == 1)){
        const BasicOp<T>* op= operator()(top);
        bool calc_hc = this_hami->use_hc();
        if (op->is_diagonal()) calc_hc = false;
                      
        if(!calc_hc)
          cout << "MISSING TERM " << t.name(true) << endl;
        else
          cout << "MISSING TERM " << t.name(true) << " + h.c." << endl;
                                                                                
        clock.Lap();
        
        new_operator(_h, *op, rho, basis, position, T(t.coef()));
        cout << "Lap: " << clock.LapTime().c_str() << endl;
      }
    }
// Multiple products
    if(t.size() >= 2 && t.type() == TERM_PRODUCT){
      bool found1 = false;
      bool found2 = false;
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        bmask |= mask(block(top.site()));
        if(block(top.site()) == pos1){
          aux_term1 *= top;
          found1 = true;
        }
        if(block(top.site()) == pos2){
          aux_term2 *= top;
          found2 = true;
        }
      } 
      if((bmask & _mask) == bmask && found1 && found2){ // we have a new composite operator

        BasicOp<T> top1, top2;
        if(aux_term1.size() == 1) // the piece contains a single operator
          top1 = aux_term1[0].internals();
        else
          top1 = aux_term1; // the piece contains more than one operator
        if(aux_term2.size() == 1) // the piece contains a single operator
          top2 = aux_term2[0].internals();
        else
          top2 = aux_term2; // the piece contains more than one operator
        const BasicOp<T>* op1 = operator()(top1);
        const BasicOp<T>* op2 = operator()(top2);

        if(!op1) {
          cout << "ERROR : Operator " << top1.description() << " not found\n";
          continue;
        }
        if(!op2) {
          cout << "ERROR : Operator " << top2.description() << " not found\n";
          continue;
        }

        bool calc_hc = this_hami->use_hc();
        if (op1->is_diagonal() && op2->is_diagonal()) calc_hc = false;
//        if (t.size() > 2) {
////          calc_hc = false; // you have to add the h.c. terms explicitly
//          bool is_diag = true;
//          for(int i = 0; i < t.size(); i++) 
//            if(t[i].is_diagonal()) { is_diag = false; break; }
//          if(is_diag) calc_hc = false;
//        }

        if(op1 && op2){
          clock.Lap();
  
          if(!calc_hc)
            cout << "TERM " << t.description() << endl;
          else
            cout << "TERM " << t.description() << " + h.c." << endl;
          new_operator(_h, *op1, *op2, pos1, pos2, rho, basis, T(t.coef()), calc_hc);

          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
  }

  b.push_back(_h);

}

/////////////////////////////////////////////////////////////////////////
// rotate_additive:
// rotate additive operators 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate_additive(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Term<T> &term)
{
  const BasicOp<T>& top1 = term[0];
  const BasicOp<T>& top2 = term[1];

  if(!top2.is_hami()) {
    cout << "WARNING: (rotate_additive) Operator is not additive\n";
  }

  cout << "ROTATE ADDITIVE " << term.name() << endl;

  const Block<T> *pb1 = _b1;
  const Block<T> *pb2 = _b2;
  size_t pos1 = BLOCK1;
  size_t pos2 = BLOCK2;
  int _mask = MASK_BLOCK1 | MASK_BLOCK2;

  if(position == RIGHT){
    pb1 = _b3;
    pb2 = _b4;
    pos1 = BLOCK4;
    pos2 = BLOCK3;
    _mask = MASK_BLOCK3 | MASK_BLOCK4;
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);

  CTimer clock;
  clock.Start();

  BasicOp<T> this_op(top2);
  this_op.resize(rho_basis);
//////////////////////////////////////////////////////////////////
// NOW WE ADD THE SECOND TERM WHICH HAS THE SUM
/////////////////////////////////////////////////////////////////////////

  if(b.contains(top2.name(), 0, 0, OP_ADDITIVE)) return;

  bool found = operator[](pos1).contains(top2.name(), 0, 0, OP_ADDITIVE);
  bool add = false;
  
  if(found){ // we transform old operator 
    clock.Lap();
    const BasicOp<T>* op = operator[](pos1)(top2); 
    if(op){
      cout << "NEW ADDITIVE OPERATOR 1 " << this_op.description().c_str() << endl;
      new_operator(this_op, *op, rho, basis, position);
      cout << "Lap: " << clock.LapTime().c_str() << endl;
      add = true;
    } else {
      cout << "WARNING: " << op->name() << " not found\n";
    }
  }

  typename Hami<T>::iterator hiter;
  Hami<T> this_hami = top2.hami();
  for(hiter = this_hami.begin(); hiter != this_hami.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;
    BasicOp<T>top = t[0].internals();

    size_t ib = block(top.site());
    if((ib == BLOCK2 && position == LEFT) || (ib == BLOCK3 && position == RIGHT)){
      const BasicOp<T>* op = operator[](ib)(top.set_site(0));

      cout << "TERM " << t.name(true) << endl;

      clock.Lap();
      if(op){
        cout << "ADDING NEW ADDITIVE OPERATOR 2 " << this_op.name() << endl;
        new_operator(this_op, *op, rho, basis, 1-position, T(t.coef()));
        add = true;
      } else {
        cout << "WARNING: " << top.name() << " not found\n";
      } 
      cout << "Lap: " << clock.LapTime().c_str() << endl;
    }
    if(!found && ((ib == BLOCK1 && position == LEFT) || (ib == BLOCK4 && position == RIGHT))){
      const BasicOp<T>* op = operator()(top);

      cout << "TERM " << t.name(true) << endl;

      clock.Lap();
      if(op){
        cout << "ADDING NEW ADDITIVE OPERATOR 1 " << this_op.name() << endl;
        new_operator(this_op, *op, rho, basis, position, T(t.coef()));
        add = true;
      } else {
        cout << "WARNING: " << top.name() << " not found\n";
      } 
      cout << "Lap: " << clock.LapTime().c_str() << endl;
    }
  }

  if(add) b.push_back(this_op);
}

} // namespace dmtk

#endif // __DMTK_SYSTEM_SVD_H__
