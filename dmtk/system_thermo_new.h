#ifndef __DMTK_SYSTEM_TIME_H__
#define __DMTK_SYSTEM_TIME_H__

#include <dmtk/dmtk.h>
#include <iostream>
#include <iosfwd>
#include <list>
#include <string>
#include <iomanip>

using namespace dmtk;

#define NTARGETS 4

namespace dmtk
{

template<class T>
class SystemThermo:public System<T>
{
  private:
    SystemThermo() {} // no default constructor
    bool _apply_product_time;
    bool time_evolve;
    bool _use_rk;
    int _n_lanczos;
    bool _imaginary_time;
    bool _in_empty_sweep;

    typedef dmtk::System<T> _S;
  protected:

  public:
    BasicOp<T> hij; 
    Vector<BasicOp<T>* > hsites;
    bool _propagate_gs;
    bool _evolve_gs;
    double dt;
    double my_time;
    double phase;

    SystemThermo(const Hami<T> &_h, const Lattice& lattice, double _dt, const char *name):System<T>(_h, lattice, name), _evolve_gs(false), _propagate_gs(false), _apply_product_time(false), time_evolve(false), _imaginary_time(true), my_time(0), phase(0), _n_lanczos(-1), _use_rk(true), _in_empty_sweep(false)
      { hsites.resize(lattice.size()); dt = _dt; }

    void init_time_iteration(const Block<T>&b1, const Block<T>&b2, const Block<T>&b3, const Block<T>&b4, bool rotate_states = true, bool build_operators = true);
    void build_target_states(int pos = 0);
    void empty_sweep(size_t _dir, size_t _start, size_t _end, bool build_operators = true);
    void start_sweep(const BasicOp<T> *op, size_t t, size_t _dir = RIGHT2LEFT, size_t _start = 0, size_t _block = BLOCK2);
    void start_sweep(const Hami<T> *h, size_t t, size_t _dir = RIGHT2LEFT, size_t _start = 0);
    void time_sweep(size_t t, size_t _dir, size_t _start, bool build_operators);
    void time_sweep1(size_t t, size_t _dir, size_t _iter, bool build_operators);
    void evolve(bool update_gs = true);
    void set_evolve_gs(bool _evolve) { _evolve_gs = _evolve; }

    SystemThermo &set_use_rk(bool b) { _use_rk = b; return *this; }
    bool use_rk() { return _use_rk; }

    void thermo_field(size_t state = 0, Vector<size_t> *states = NULL, size_t first = 0);

    BasicOp<T>& get_hij(int pos)
      { if(hsites[pos]) return *hsites[pos]; else return hij; }
    const BasicOp<T>& get_hij(int pos) const
      { if(hsites[pos]) return *hsites[pos]; else return hij; }

    void propagate_gs(const VectorState<T>& _gs, bool evolve = false) { _propagate_gs = true; _evolve_gs = evolve, _S::_target.resize(2); _S::_target[0] = _gs; _S::_target[1] = _gs; _S::_target_weight.resize(2); _S::_target_weight[0] = 0.5; _S::_target_weight[1] = 0.5; _S::write_gs(_S::gs, _S::iter, _S::dir); }

    SystemThermo &set_apply_product_time(bool b) { _apply_product_time = b; return *this; }
    bool apply_product_time() const { return _apply_product_time; }

    int n_exp_lanczos() const { return _n_lanczos; }
    SystemThermo &set_n_exp_lanczos(int n_lanczos)
      { _n_lanczos = n_lanczos; return *this; }

    bool imaginary_time() const { return _imaginary_time; }
    SystemThermo &set_imaginary_time(bool use)
      { _imaginary_time = use; return *this; }

    void write_status() const
    {
      _S::write_status();

      char file[255];
      snprintf(file, 255, "system_thermo_%s.dat", _S::_name);
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
      outputfile.write((const char *)&my_time, sizeof(double));
      outputfile.close();
    }

    virtual void read_status()
    {
      char file[255];
      snprintf(file, 255, "system_thermo_%s.dat", _S::_name);
      ifstream real_file(file,std::ios::in|std::ios::binary);
      if(!real_file) {
        cout << "*** ERROR: Could not open file " << file << endl;
        return;
      }
#ifdef WITH_BZIP2
      bzip2_stream::bzip2_istream inputfile(real_file);
#else
      ifstream &inputfile = real_file;
#endif
      inputfile.read((char *)&my_time, sizeof(double));
      inputfile.close();
    }
};

#ifdef USE_RK
#warning USING RK - if you do not want this undef USE_RK
// Declaration of function "product_time"
// You have to provide the implementation
VectorState<complex<double> >
product_time(SystemThermo<complex<double> > &s, const VectorState<complex<double> > &vv, double xtime);

#ifndef USE_CUSTOM_PRODUCT
VectorState<complex<double> >
product(SystemThermo<complex<double> > &s, const VectorState<complex<double> > &vv)
{
  VectorState<complex<double> >aux(vv);
  System<complex<double> > &_s = s;
  aux = product(_s, vv);
  if(s.apply_product_time()) aux += product_time(s, vv, s.my_time);
  return aux;
}

VectorState<double>
product(SystemThermo<double> &s, const VectorState<double> &vv)
{
  System<double> &_s = s;
  return product(_s, vv);
}
#endif // USE_CUSTOM_PRODUCT

VectorState<double>
exp_h(SystemThermo<double> &S, const VectorState<double> & evolve_vector)
{
  Vector<VectorState<double> >aux(1);
  Vector<double> _e(1);
  Vector<double> a(100), b(100);
  Vector<double> d(100), e(100);
  int tridiag_size;
  VectorState<double>aux_vector = evolve_vector;
  aux[0] = aux_vector;

  S.set_apply_product_time(true);

  cout << ">>>>>>> LANCZOS <<<<<<<\n";
  tridiag_size = S.n_exp_lanczos() == -1 ? 15 : S.n_exp_lanczos();
  lanczos<double, SystemThermo<double>, VectorState<double> >(S, aux, aux_vector, _e, 1, a, b, tridiag_size, 1.e-16, true, true, "vectors_test.dat", true);
  Matrix<double> z(100,100);
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
    m_exp_h(i,i) = std::exp(-d[i+1]*S.dt);

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

VectorState<complex<double> >
exp_h(SystemThermo<complex<double> > &S, const VectorState<complex<double> > &evolve_vector)
{
  Vector<VectorState<complex<double> > >aux(1);
  Vector<double> _e(1);
  Vector<double> a(100), b(100);
  Vector<double> d(100), e(100);
  int tridiag_size;
  VectorState<complex<double> >aux_vector = evolve_vector;
  aux[0] = aux_vector;

  S.set_apply_product_time(true);

  cout << ">>>>>>> LANCZOS <<<<<<<\n";
  tridiag_size = S.n_exp_lanczos() == -1 ? 15 : S.n_exp_lanczos();
  lanczos<complex<double>, SystemThermo<complex<double> >, VectorState<complex<double> > >(S, aux, aux_vector, _e, 1, a, b, tridiag_size, 1.e-16, true, true, "vectors_test.dat", true);
  Matrix<double> z(100,100);
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
    complex<double> coef = S.imaginary_time() ? complex<double>(-d[i+1]*S.dt) : complex<double>(0.,-(d[i+1]-S.phase)*S.dt);
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

VectorState<double>
exp_h_rk(SystemThermo<double> &s, const VectorState<double> &evolve_vector)
{
  VectorState<double> new_state, k_state;
  VectorState<double> hphi, aux_hphi;
  VectorState<double> vphi, aux_vphi;
  VectorState<double> gs = evolve_vector;
  double dt = s.dt * 0.1;
  double my_time = s.my_time;
  double one_sixth(1./6.);
  double one_third(1./3.);
  double one_half(0.5);
  double c_energy(s.phase);

  new_state = gs;
  k_state = gs;

  s.set_apply_product_time(false);

  for(int i = 1; i <= 10; i++){
    hphi = product_default(s, gs);
//    hphi = product(s, gs);
    hphi -= c_energy*gs;

    k_state = hphi;
    k_state *= -dt;

    new_state += one_sixth*k_state;

    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt;

    new_state += one_third*k_state;

    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt;

    new_state += one_third*k_state;

    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += aux_hphi;
    k_state *= -dt;

    new_state += one_sixth*k_state;
    gs = new_state;
    my_time += dt;
  }

  s.set_apply_product_time(true);

  return new_state;
}

VectorState<complex<double> >
exp_h_rk(SystemThermo<complex<double> > &s, const VectorState<complex<double> > &evolve_vector)
{
  VectorState<complex<double> > new_state, k_state;
  VectorState<complex<double> > hphi, aux_hphi;
  VectorState<complex<double> > vphi, aux_vphi;
  VectorState<complex<double> > gs = evolve_vector;
  double dt = s.imaginary_time() ? s.dt*0.1 : s.dt * 0.1;
  double my_time = s.my_time;
  complex<double> one_sixth(1./6.,0.);
  complex<double> one_third(1./3.,0.);
  complex<double> one_half(0.5,0.);
  complex<double> c_energy(s.phase,0.);

  new_state = gs;
  k_state = gs;
  complex<double> coef = s.imaginary_time() ? -dt : complex<double>(0,-dt);

  s.set_apply_product_time(false);

  for(int i = 1; i <= 10; i++){
    hphi = product_default(s, gs);
//    hphi = product(s, gs);
    hphi -= c_energy*gs;

    vphi = product_time(s, gs, my_time);

    k_state = hphi + vphi;
    k_state *= coef;

    new_state += one_sixth*k_state;

    vphi = product_time(s, gs, my_time+0.5*dt);
    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(s, k_state, my_time+0.5*dt);
    k_state = hphi + vphi;
    k_state += one_half*(aux_hphi + aux_vphi);
    k_state *= coef;

    new_state += one_third*k_state;

    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(s, k_state, my_time+0.5*dt);
    k_state = hphi + vphi;
    k_state += one_half*(aux_hphi + aux_vphi);
    k_state *= coef;

    new_state += one_third*k_state;

    vphi = product_time(s, gs, my_time+dt);
    aux_hphi = product_default(s, k_state);
//    aux_hphi = product(s, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(s, k_state, my_time+dt);
    k_state = hphi + vphi;
    k_state += aux_hphi + aux_vphi;
    k_state *= coef;

    new_state += one_sixth*k_state;
    gs = new_state;
    my_time += dt;
  }

  s.set_apply_product_time(true);

  return new_state;
}

template<>
void
SystemThermo<double>::evolve(bool update_gs)
{
  VectorState<double> new_state, k_state;
  VectorState<double> hphi, aux_hphi;
  VectorState<double> vphi, aux_vphi;
  double one_sixth(1./6.);
  double one_third(1./3.);
  double one_half(0.5);
  double thirty_one(31./162.);
  double fourteen(14./162.);
  double five(5./162.);
  double sixteen(16./81.);
  double twenty(20./81.);
  double two(2./81.);
  double c_energy(phase);
  CTimer clock;

  clock.Start();

  if(_S::verbose()) cout << ">>> EVOLVE <<<\n";

  set_apply_product_time(false);


  int maxstate = _evolve_gs ? 1 : 0;
  int _first_vector = _propagate_gs ? 1 : 0;

  for(int nstate = 0; nstate <= maxstate; nstate++){

    int first_vector = _first_vector + nstate*(NTARGETS-1);
 
    int evolve_vector = 1-nstate;
    if(!_propagate_gs) evolve_vector = 0;
 
    _S::gs = _S::_target[evolve_vector];
  
    new_state = _S::gs;
    k_state = _S::gs;
    _S::_target[first_vector+1] = _S::gs;
    _S::_target[first_vector+2] = _S::gs;
  
    hphi = product_default(*this, gs);
//    hphi = product(*this, _S::gs);
    hphi -= c_energy*_S::gs;
  
    k_state = hphi;
    k_state *= -dt*0.5;
    cout << "NEW KSTATE 1 " << sqrt(product(k_state,k_state)) << endl;
  
  //  new_state += k_state; // Euler - TEST
  
    new_state += one_sixth*k_state;
    _S::_target[first_vector+1] += thirty_one*k_state;
    _S::_target[first_vector+2] += sixteen*k_state;
  
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt*0.5;
    cout << "NEW KSTATE 2 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_third*k_state;
    _S::_target[first_vector+1] += fourteen*k_state;
    _S::_target[first_vector+2] += twenty*k_state;
  
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += one_half*aux_hphi;
    k_state *= -dt*0.5;
    cout << "NEW KSTATE 3 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_third*k_state;
    _S::_target[first_vector+1] += fourteen*k_state;
    _S::_target[first_vector+2] += twenty*k_state;
  
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    k_state = hphi;
    k_state += aux_hphi;
    k_state *= -dt*0.5;
    cout << "NEW KSTATE 4 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_sixth*k_state;
    _S::_target[first_vector+1] -= five*k_state;
    _S::_target[first_vector+2] -= two*k_state;
    _S::_target[first_vector+3] = new_state;
  
    if(update_gs){
      clock.Lap();
  #ifdef USE_EXP_RK
      new_state = exp_h_rk(*this, _S::_target[evolve_vector]);
  #else
      new_state = exp_h(*this, _S::_target[evolve_vector]);
  #endif
      my_time += dt;
      _S::_target[evolve_vector] = new_state;
    }
  } 

  _S::gs = _S::_target[_first_vector];
  _S::write_gs(_S::gs, _S::iter, _S::dir);
  cout << "NORM GS:" << sqrt(product(_S::gs,_S::gs)) << endl;
  cout << "CPU time(new_gs):" << clock.LapTime() << endl;

/*
  double wsum = 0;
  if(_propagate_gs)
    _S::_target_weight[first_vector] = 0.5;
  else
    _S::_target_weight[first_vector] = 0.4;

  int last_vector = _S::_target.size()-1;
  _S::_target_weight[last_vector] = 0.4;
  if(!_propagate_gs){
    for(int i = 0; i < _S::_target.size(); i++)
      _S::_target_weight[i] = 1./_S::_target.size();
  } else {
    for(int i = 1; i < _S::_target.size(); i++)
      _S::_target_weight[i] = (1.-_S::_target_weight[0])/(_S::_target.size()-1);
  }
*/

  if(_propagate_gs){
    if(!_evolve_gs){
      _S::_target_weight[0] = 0.5;
      _S::_target_weight[1] = 0.25;
      _S::_target_weight[2] = _S::_target_weight[3] = _S::_target_weight[4] = 0.25/3;
    } else {
      _S::_target_weight[0] = 0.25;
      _S::_target_weight[1] = 0.25;
      for(int i = 2; i < 8; i++) _S::_target_weight[i] = 0.5/6;
    }
  } else {
    _S::_target_weight[0] = 0.5;
    _S::_target_weight[1] = _S::_target_weight[2] = _S::_target_weight[3] = 0.5/3;
  }

  for(int i = 0; i < _S::_target.size(); i++)
    cout << "TARGET WEIGHT " << i << " " << _S::_target_weight[i] << endl;
  

//  for(int i = 1; i < 4; i++){
//    T overlap = product(_S::_target[0],_S::_target[i]);
//    double x = product(_S::_target[first_vector],_S::_target[i]).real();
//    cout << "TARGET NORM = " << sqrt(x) << " " << overlap << endl;
//    _S::_target[i] -= (overlap/x)*_S::_target[0];
//  }

  set_apply_product_time(true);
  cout << "CPU time(evolve):" << clock.TotalTime() << endl;
}

template<>
void
SystemThermo<complex<double> >::evolve(bool update_gs)
{
  VectorState<complex<double> > new_state, k_state;
  VectorState<complex<double> > hphi, aux_hphi;
  VectorState<complex<double> > vphi, aux_vphi;
  complex<double> one_sixth(1./6.,0.);
  complex<double> one_third(1./3.,0.);
  complex<double> one_half(0.5,0.);
  complex<double> thirty_one(31./162.,0.);
  complex<double> fourteen(14./162.,0.);
  complex<double> five(5./162.,0.);
  complex<double> sixteen(16./81.,0.);
  complex<double> twenty(20./81.,0.);
  complex<double> two(2./81.,0.);
  complex<double> c_energy(phase,0.);
  complex<double> coef = _imaginary_time ? -dt : complex<double>(0.,-dt);
  CTimer clock;

  clock.Start();

  if(_S::verbose()) cout << ">>> EVOLVE <<<\n";

  set_apply_product_time(false);

  int maxstate = _evolve_gs ? 1 : 0;
  int _first_vector = _propagate_gs ? 1 : 0;

  for(int nstate = 0; nstate <= maxstate; nstate++){

    int first_vector = _first_vector + nstate*(NTARGETS-1);
 
    int evolve_vector = 1-nstate;
    if(!_propagate_gs) evolve_vector = 0;
 
    _S::gs = _S::_target[evolve_vector];
  
    new_state = _S::gs;
    k_state = _S::gs;
    _S::_target[first_vector+1] = _S::gs;
    _S::_target[first_vector+2] = _S::gs;
  
    hphi = product_default(*this, gs);
//    hphi = product(*this, _S::gs);
    hphi -= c_energy*_S::gs;
  
    vphi = product_time(*this, _S::gs, my_time);
  
    k_state = hphi + vphi;
    k_state *= coef;
    cout << "NEW KSTATE 1 " << sqrt(product(k_state,k_state)) << endl;

  //  new_state += k_state; // Euler - TEST
  
    new_state += one_sixth*k_state;
    _S::_target[first_vector+1] += thirty_one*k_state;
    _S::_target[first_vector+2] += sixteen*k_state;
  
    vphi = product_time(*this, _S::gs, my_time+0.5*dt);
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(*this, k_state, my_time+0.5*dt);
    k_state = hphi + vphi;
    k_state += one_half*(aux_hphi + aux_vphi);
    k_state *= coef;
    cout << "NEW KSTATE 2 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_third*k_state;
    _S::_target[first_vector+1] += fourteen*k_state;
    _S::_target[first_vector+2] += twenty*k_state;
  
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(*this, k_state, my_time+0.5*dt);
    k_state = hphi + vphi;
    k_state += one_half*(aux_hphi + aux_vphi);
    k_state *= coef;
    cout << "NEW KSTATE 3 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_third*k_state;
    _S::_target[first_vector+1] += fourteen*k_state;
    _S::_target[first_vector+2] += twenty*k_state;
  
    vphi = product_time(*this, _S::gs, my_time+dt);
    aux_hphi = product_default(*this, k_state);
//    aux_hphi = product(*this, k_state);
    aux_hphi -= c_energy*k_state;
    aux_vphi = product_time(*this, k_state, my_time+dt);
    k_state = hphi + vphi;
    k_state += aux_hphi + aux_vphi;
    k_state *= coef;
    cout << "NEW KSTATE 4 " << sqrt(product(k_state,k_state)) << endl;
  
    new_state += one_sixth*k_state;
    _S::_target[first_vector+1] -= five*k_state;
    _S::_target[first_vector+2] -= two*k_state;
    _S::_target[first_vector+3] = new_state;
  
    if(update_gs){
      clock.Lap();
  #ifdef USE_EXP_RK
      new_state = exp_h_rk(*this, _S::_target[evolve_vector]);
  #else
      new_state = exp_h(*this, _S::_target[evolve_vector]);
  #endif
      my_time += dt;
      _S::_target[evolve_vector] = new_state;
    }
  }
  _S::gs = _S::_target[_first_vector];
   _S::write_gs(_S::gs, _S::iter, _S::dir);
   cout << "NORM GS:" << sqrt(product(_S::gs,_S::gs)) << endl;
   cout << "CPU time(new_gs):" << clock.LapTime() << endl;
 
/* 
  double wsum = 0;
  if(_propagate_gs)
    _S::_target_weight[first_vector] = 0.5;
  else
    _S::_target_weight[first_vector] = 0.4;

  int last_vector = _S::_target.size()-1;
  _S::_target_weight[last_vector] = 0.4;
  if(!_propagate_gs){
    for(int i = 0; i < _S::_target.size(); i++)
      _S::_target_weight[i] = 1./_S::_target.size();
  } else {
    for(int i = 1; i < _S::_target.size(); i++)
      _S::_target_weight[i] = (1.-_S::_target_weight[0])/(_S::_target.size()-1);
  }
*/

  if(_propagate_gs){
    if(!_evolve_gs){
      _S::_target_weight[0] = 0.5;
      _S::_target_weight[1] = 0.25;
      _S::_target_weight[2] = _S::_target_weight[3] = _S::_target_weight[4] = 0.25/3;
    } else {
      _S::_target_weight[0] = 0.25;
      _S::_target_weight[1] = 0.25;
      for(int i = 2; i < 8; i++) _S::_target_weight[i] = 0.5/6;
    }
  } else {
    _S::_target_weight[0] = 0.5;
    _S::_target_weight[1] = _S::_target_weight[2] = _S::_target_weight[3] = 0.5/3;
  }

  for(int i = 0; i < _S::_target.size(); i++)
    cout << "TARGET WEIGHT " << i << " " << _S::_target_weight[i] << endl;

/*
  for(int i = 1; i < 4; i++){
    T overlap = product(_S::_target[0],_S::_target[i]);
    double x = product(_S::_target[first_vector],_S::_target[i]).real();
    cout << "TARGET NORM = " << sqrt(x) << " " << overlap << endl;
    _S::_target[i] -= (overlap/x)*_S::_target[0];
  }
*/
  set_apply_product_time(true);
  cout << "CPU time(evolve):" << clock.TotalTime() << endl;
}
#else // !USE_RK
#warning USING S-T - if you do not want this define USE_RK
template <class T>
void
SystemThermo<T>::evolve(bool update_gs)
{
  my_time += dt;
}
#endif //USE_RK

template<class T>
void
SystemThermo<T>::thermo_field(size_t state, Vector<size_t> *states, size_t first)
{
  int ls = _S::lattice().ls();
  _S::_timer.Start();
  bool save_use_rk = _use_rk; // We shouldn't use rk in this part
  _use_rk = false;
  size_t save_qn_mask = _S::_grand_canonical;
//  _S::_grand_canonical = 0; 
  _S::_in_warmup = true;
  bool save_use_single_site = _S::_use_single_site;
  bool save_save_hmatrix = _S::_save_hmatrix;
  _S::_save_hmatrix = false;
  bool save_use_composite = _S::_use_composite;
  _S::_use_composite = false;


  _S::start(false);
  _S::init_correlations();

  Vector<size_t> _states(ls);
  _states = state;
  if(states) _states = *states;

  QN qn_iter;

//****************************************************************
// We build the thermo field
//****************************************************************
  for(_S::iter = 1; _S::iter <= ls-3; _S::iter++)
  {
    _S::read_block(_S::leftblock, _S::iter, LEFT);
    const Block<T>& site1 = _S::h.get_site(_S::iter);
    const Block<T>& site2 = _S::h.get_site(_S::iter+1);
    _S::rightblock = _S::h.get_site(_S::iter+2);

    cout << "THERMO-FIELD BUILDING ITERATION (LEFT-TO-RIGHT) " << _S::iter << endl;

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, false);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

    cout << _S::_b1->basis().size() << " " << _S::_b2->basis().size() << " " << _S::_b3->basis().size() << " " << _S::_b4->basis().size() << endl;
    qn_iter = QN();
    if(_S::iter == 1){
      QN qn1 = _S::_b1->basis()[_states[0]].qn();
      QN qn2 = _S::_b2->basis()[_states[1]].qn();
      QN qn3 = _S::_b3->basis()[_states[2]].qn();
      QN qn4 = _S::_b4->basis()[_states[3]].qn();
      qn_iter = qn1 + qn2 + qn3 + qn4;
      cout << _states[0] << " " << qn1.n() << " " << qn1.sz() << endl;
      cout << _states[1] << " " << qn2.n() << " " << qn2.sz() << endl;
      cout << _states[2] << " " << qn3.n() << " " << qn3.sz() << endl;
      cout << _states[3] << " " << qn4.n() << " " << qn4.sz() << endl;
      _S::gs.set_qn_mask(qn_iter, _S::_grand_canonical);
      _S::gs.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      cout << "GS SIZE " << _S::gs.size() << endl;
      _S::gs = T(0);
      _S::gs(_states[0],_states[1],_states[2],_states[3]) = T(1.);
      _S::measure();
    }else{
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      QN qn1 = _S::_b1->basis()[first].qn();
      QN qn2 = _S::_b2->basis()[_states[_S::iter]].qn();
      QN qn3 = _S::_b3->basis()[_states[_S::iter+1]].qn();
      QN qn4 = _S::_b4->basis()[_states[_S::iter+2]].qn();
/*
      qn_iter = qn1 + qn2 + qn3 + qn4;
      cout << first << " " << qn1.n() << " " << qn1.sz() << endl;
      cout << _states[iter] << " " << qn2.n() << " " << qn2.sz() << endl;
      cout << _states[iter+1] << " " << qn3.n() << " " << qn3.sz() << endl;
      cout << _states[iter+2] << " " << qn4.n() << " " << qn4.sz() << endl;
      _S::gs.set_qn_mask(qn_iter, _S::_grand_canonical);
      _S::gs.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      cout << "GS SIZE " << _S::gs.size() << endl;
      _S::gs = T(0);
      _S::gs(first,_states[iter],_states[iter+1],_states[iter+2]) = T(1.);
*/
      read_rho(rho1, basis1, _S::iter, LEFT);
      read_gs(_S::gs, _S::iter-1, LEFT);
      _S::seed = _S::gs;
      _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      _S::seed = T(0);

      new_seed(_S::gs, _S::seed, rho1, basis1, LEFT);
      qn_iter = _S::gs.qn() + qn4;
      cout << "QN " << qn_iter.n() << " " << qn_iter.sz() << endl;
      _S::gs.set_qn_mask(qn_iter, _S::_grand_canonical);
      _S::gs.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      _S::gs = T(0);

      typename VectorState<T>::iterator siter;
      for(siter = _S::seed.subspace_begin(); siter != _S::seed.subspace_end(); siter++){
        StateSpace s = *siter;
        StateSpace sgs = _S::gs.get_qn_space(s[1].qn(),s[3].qn(),s[4].qn(),qn4);
        state_slice<T> ss = _S::seed(s);
        state_slice<T> ssgs = _S::gs(sgs);
        for(int i1 = 0; i1 < s[1].dim(); i1++){
          for(int i2 = 0; i2 < s[3].dim(); i2++){
            for(int i3 = 0; i3 < s[4].dim(); i3++){
              ssgs(i1,i2,i3,_states[_S::iter+2]-sgs[4].begin()) = ss(i1,0,i2,i3);  
            }
          }
        }

      }
    }
    _S::measure();

    _S::truncate(LEFT, site1.basis().size()*site1.basis().size());
    _S::rotate(LEFT, _S::newblock);
    _S::write_iter(LEFT);
  }

  for(_S::iter = 1; _S::iter <= ls-3; _S::iter++)
  {
    _S::read_block(_S::leftblock, ls-_S::iter-2, LEFT);
    _S::read_block(_S::rightblock, _S::iter, RIGHT);
    const Block<T>& site1 = _S::h.get_site(ls-2-_S::iter);
    const Block<T>& site2 = _S::h.get_site(ls-2-_S::iter+1);

    cout << "THERMO-FIELD RIGHT-TO-LEFT ITERATION " << _S::iter << endl;
    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, false);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

    cout << _S::_b1->basis().size() << " " << _S::_b2->basis().size() << " " << _S::_b3->basis().size() << " " << _S::_b4->basis().size() << endl;
    cout << "GS SIZE " << _S::gs.size() << endl;

    if(_S::iter != 1){
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      cout << "QN GS " << _S::gs.qn().n() << " " << _S::gs.qn().sz() << endl;
      cout << "QN " << _S::qnt.n() << " " << _S::qnt.sz() << endl;
      read_rho(rho1, basis1, _S::iter, RIGHT);
      read_rho(rho2, basis2, _S::lattice().ls()-_S::iter-2+1, LEFT);
      read_gs(_S::gs, _S::iter-1, RIGHT);
      cout << "NEW SEED" << endl;
      _S::seed.set_qn_mask(_S::gs.qn(), _S::_grand_canonical);
      _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      new_seed(_S::gs, _S::seed, rho1, rho2, basis1, basis2, RIGHT);
      _S::gs = _S::seed;
    }
    _S::measure();

    _S::truncate(RIGHT, site2.basis().size()*site2.basis().size());
    _S::rotate(RIGHT, _S::newblock);
    _S::write_iter(RIGHT);
  }
  for(_S::iter = 1; _S::iter <= ls-3; _S::iter++)
  {
    _S::read_block(_S::leftblock, _S::iter, LEFT);
    _S::read_block(_S::rightblock, ls-_S::iter-2, RIGHT);
    const Block<T>& site1 = _S::h.get_site(_S::iter);
    const Block<T>& site2 = _S::h.get_site(_S::iter+1);

    cout << "THERMO-FIELD LEFT-TO-RIGHT ITERATION " << _S::iter << endl;

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, false);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

    cout << _S::_b1->basis().size() << " " << _S::_b2->basis().size() << " " << _S::_b3->basis().size() << " " << _S::_b4->basis().size() << endl;
    cout <<"GS SIZE " << _S::gs.size() << endl;

    if(_S::iter != 1){
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      read_rho(rho1, basis1, _S::iter, LEFT);
      read_rho(rho2, basis2, _S::lattice().ls()-_S::iter-2+1, RIGHT);
      read_gs(_S::gs, _S::iter-1, LEFT);
      _S::gs.resize(_S::_grand_canonical);
      cout << "NEW SEED" << endl;
      _S::seed.set_qn_mask(_S::gs.qn(), _S::_grand_canonical);
      _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
      new_seed(_S::gs, _S::seed, rho1, rho2, basis1, basis2, LEFT);
      _S::gs = _S::seed;
    }
    _S::measure();

    if(_S::iter < ls-3){
      _S::truncate(LEFT,site1.basis().size()*site1.basis().size());
      _S::rotate(LEFT, _S::newblock);
    }
    _S::write_iter(LEFT);
  }

  _use_rk = save_use_rk;

  _S::write_gs(_S::gs, ls-2, LEFT);
  _S::_target[0] = product_default(*this, _S::gs);
  cout << "THERMO FIELD VARIATIONAL ENERGY = "  << product(_S::gs,_S::_target[0])/product(_S::gs,_S::gs) << endl;
  _S::_use_single_site = save_use_single_site;
  _S::_save_hmatrix = save_save_hmatrix;
  _S::_use_composite = save_use_composite;
}

template<class T>
void
SystemThermo<T>::time_sweep1(size_t t, size_t _dir = RIGHT2LEFT, size_t _iter = 1, bool build_operators = true)
{
  int ls = _S::lattice().ls();
  CTimer clock;

  clock.Start();

  char file[255];
  snprintf(file,255,"iter_%s.dat",_S::_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "Time sweep\n";

  _S::dir = _dir;
  _S::iter = _iter;
  _S::signal_emit(SYSTEM_SIGNAL_START_ITER);
//****************************************************************
  cout << "NEW TIME SWEEP\n";

  size_t first_state = _propagate_gs ? 1 : 0;
  if(_evolve_gs) first_state = 0;

  if(_dir == LEFT2RIGHT){
    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
    int pos = _iter;
    if(_iter > ls -3) pos--;
    const Block<T>& site1 = _S::h.get_site(pos);
    const Block<T>& site2 = _S::h.get_site(pos+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    _S::_timer.Lap();
    cout << "TIME SWEEP ITERATION " << endl;
    outputfile << "TIME SWEEP ITERATION " << endl;
    cout << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, true, build_operators);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;


    for(int n = first_state; n < _S::_target.size(); n++){
      if(_S::iter != 0 && _S::iter != ls-2){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK2" << endl;
        VectorState<T> v23 = _S::_target[n].condense(MASK_BLOCK2|MASK_BLOCK3);
        VectorState<T> res23(v23);
        res23 = T(0);
        product(get_hij(_S::iter), v23, res23, BLOCK2);
        _S::_target[n] = res23.decondense(MASK_BLOCK2|MASK_BLOCK3, _S::_target[n]);
      } 
      if(_S::iter == 0){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK1" << endl;
        VectorState<T> v12 = _S::_target[n].condense(MASK_BLOCK1|MASK_BLOCK2);
        VectorState<T> res12(v12);
        res12 = T(0);
        product(get_hij(_S::iter), v12, res12, BLOCK1);
        _S::_target[n] = res12.decondense(MASK_BLOCK1|MASK_BLOCK2, _S::_target[n]);
      } 
      if(_S::iter == ls-2){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK4" << endl;
        VectorState<T> v34 = _S::_target[n].condense(MASK_BLOCK3|MASK_BLOCK4);
        VectorState<T> res34(v34);
        res34 = T(0);
        product(get_hij(_S::iter), v34, res34, BLOCK4);
        _S::_target[n] = res34.decondense(MASK_BLOCK3|MASK_BLOCK4, _S::_target[n]);
      } 
    }
    
    _S::gs = _S::_target[first_state];
   
    if(_S::iter > 0 && _S::iter < ls-3){
      _S::m = t;
      _S::m1 = std::min(_S::m1*_S::m2,_S::m);

      _S::truncate(LEFT, _S::m1);
      _S::rotate(LEFT, _S::newblock);
      _S::write_iter(LEFT);
    } else {
//      build_S::_target_states(0);
      _S::write_gs(_S::gs,_S::iter,LEFT);
    }

    cout << "===========================================\n";
    cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
    cout << "===========================================\n";
    outputfile << "===========================================\n";
    outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
    outputfile << "===========================================\n";
  } else {
  
    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
    int pos = _iter;
    if(_iter > ls -3) pos--;
    const Block<T>& site1 = _S::h.get_site(ls-2-pos);
    const Block<T>& site2 = _S::h.get_site(ls-2-pos+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    _S::_timer.Lap();
    cout << "TIME SWEEP ITERATION " << endl;
    outputfile << "TIME SWEEP ITERATION " << endl;
    cout << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;
    outputfile << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, true, build_operators);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

    for(int n = first_state; n < _S::_target.size(); n++){
      if(_S::iter != 0 && _S::iter != ls-2){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK2" << endl;
        VectorState<T> v23 = _S::_target[n].condense(MASK_BLOCK2|MASK_BLOCK3);
        VectorState<T> res23(v23);
        res23 = T(0);
        product(get_hij(ls-2-_S::iter), v23, res23, BLOCK2);
        _S::_target[n] = res23.decondense(MASK_BLOCK2|MASK_BLOCK3, _S::_target[n]);
      }
      if(_S::iter == ls-2){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK1" << endl;
        VectorState<T> v12 = _S::_target[n].condense(MASK_BLOCK1|MASK_BLOCK2);
        VectorState<T> res12(v12);
        res12 = T(0);
        product(get_hij(ls-2-_S::iter), v12, res12, BLOCK1);
        _S::_target[n] = res12.decondense(MASK_BLOCK1|MASK_BLOCK2, _S::_target[n]);
      } 
      if(_S::iter == 0){
        cout << "PRODUCT H" << ls-2-_S::iter << " BLOCK4" << endl;
        VectorState<T> v34 = _S::_target[n].condense(MASK_BLOCK3|MASK_BLOCK4);
        VectorState<T> res34(v34);
        res34 = T(0);
        product(get_hij(ls-2-_S::iter), v34, res34, BLOCK4);
        _S::_target[n] = res34.decondense(MASK_BLOCK3|MASK_BLOCK4, _S::_target[n]);
      } 
    }

    _S::gs = _S::_target[first_state];
 
    if(_S::iter > 0 && _S::iter < ls-3){
      _S::m = t;
      _S::m4 = std::min(_S::m4*_S::m3,_S::m);

      _S::truncate(RIGHT, _S::m4);
      _S::rotate(RIGHT, _S::newblock);
      _S::write_iter(RIGHT);
    } else {
//      build_S::_target_states(0);
      _S::write_gs(_S::gs,_S::iter,RIGHT);
    }

    cout << "===========================================\n";
    cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
    cout << "===========================================\n";
    outputfile << "===========================================\n";
    outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
    outputfile << "===========================================\n";
  }

  _S::signal_emit(SYSTEM_SIGNAL_END_ITER);
  outputfile.close();
  cout << "CPU time (time_sweep1): " << clock.TotalTime() << endl;
}


template<class T>
void
SystemThermo<T>::time_sweep(size_t t, size_t _dir = RIGHT2LEFT, size_t _start = 0, bool build_operators = true)
{
  int ls = _S::lattice().ls();
  CTimer clock;

  clock.Start();

  char file[255];
  snprintf(file,255,"iter_%s.dat",_S::_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "Time sweep\n";
  outputfile.close();

//****************************************************************
  cout << "NEW TIME SWEEP\n";
  for(int _iter = _start; _iter < ls-1; _iter++){
    time_sweep1(t, _dir, _iter, build_operators);
    if (_S::signal_emit(SYSTEM_SIGNAL_MEASURE)) _S::measure();
  }
//****************************************************************
  _S::signal_emit(SYSTEM_SIGNAL_END_SWEEP);
  cout << "CPU time (time_sweep): " << clock.TotalTime() << endl;
}

template<class T>
void
SystemThermo<T>::start_sweep(const BasicOp<T> *op, size_t t, size_t _dir, size_t _start, size_t _block)
{
  int ls = _S::lattice().ls();
  _in_empty_sweep = true;

  time_evolve = true;

  _S::iter = _start;
  _S::dir = _dir;
  _S::m = t;
  if(_S::iter == 1)
#ifdef USE_RK
    read_gs(_S::gs, ls-2-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#else
    read_gs(_S::gs, ls-1-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#endif
  else
    read_gs(_S::gs, _S::iter-1, _dir == LEFT2RIGHT ? LEFT : RIGHT);

  size_t first_state = _propagate_gs ? 1 : 0;
#ifdef USE_RK
  if(_use_rk){
    _S::_target.resize(NTARGETS+first_state);
    _S::_target_weight.resize(NTARGETS+first_state);
    _S::set_ntargets(NTARGETS+first_state);
    if(_propagate_gs && _evolve_gs){
      _S::_target.resize(2*NTARGETS+first_state-1);
      _S::_target_weight.resize(2*NTARGETS+first_state-1);
      _S::set_ntargets(2*NTARGETS+first_state-1);
    }
  }
#else
  _S::_target.resize(1+first_state);
  _S::_target_weight.resize(1+first_state);
  _S::set_ntargets(1+first_state);
#endif

  if(op){
    for(int i = 0; i < _S::_ntargets; i++) _S::_target[i] = _S::gs;
    cout << "APPLYING OPERATOR " << op->description() << " " << _block << endl;
    cout << "GS NORM = " << product(_S::gs,_S::gs) << endl;
//    cout << _S::qn.n() << " " << _S::qn.sz() << endl;
//    cout << op->dqn.n() << " " << op->dqn.sz() << endl;

    _S::qn = _S::qn+op->dqn;
    VectorState<T> _seed(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis(),_S::gs.qn()+op->dqn);
    cout << "SEED NORM = " << product(_S::seed,_S::seed) << endl;
                                                                                
    product(*op, _S::gs, _seed, _block);
    cout << _seed.qn().n() << " " << _seed.qn().sz() << endl;
    _S::qnt = _S::qnt+op->dqn;
    cout << "NORM = " << product(_seed,_seed) << endl;
//  _seed = _seed / sqrt(product(_seed,_seed));
    
    _S::_target[first_state] = _seed;
    _S::gs = _S::_target[first_state];
  } else {
    _S::_target[first_state] = _S::gs;
  }
  if(_start == 1)
#ifdef USE_RK
    write_gs(_S::gs, ls-2-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#else
    write_gs(_S::gs, ls-1-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#endif
  else
   _S::write_gs(_S::gs, _S::iter-1, _dir == LEFT2RIGHT ? LEFT : RIGHT);

/*
  if(_dir == LEFT2RIGHT){

    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
    const Block<T>& site1 = _S::h.get_site(_S::iter);
    const Block<T>& site2 = _S::h.get_site(_S::iter+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);

  } else {
   
    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
    const Block<T>& site1 = _S::h.get_site(ls-2-_S::iter);
    const Block<T>& site2 = _S::h.get_site(ls-2-_S::iter+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
  }
*/
/*
  _S::m = t;
  if(_S::dir == LEFT2RIGHT){
    _S::m1 = std::min(_S::m1*_S::m2,_S::m);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;
    _S::truncate(LEFT, _S::m1);
    _S::rotate(LEFT, _S::newblock);
    _S::write_iter(LEFT);
  } else {
    _S::m4 = std::min(_S::m4*_S::m3,_S::m);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;
    _S::truncate(RIGHT, _S::m4);
    _S::rotate(RIGHT, _S::newblock);
    _S::write_iter(RIGHT);
  }
  _S::gs = _S::_target[0]; 
//  _S::measure_n(1, &_S::_target[first_state]);
//  _S::measure(&_S::_target[first_state]);
  _S::gs = _S::_target[first_state]; 
  _S::write_gs(_S::gs, _S::iter, _S::dir == LEFT2RIGHT ? LEFT : RIGHT);
*/

#ifdef USE_RK
/*
  _S::read_block(_S::rightblock, ls-_start-2, RIGHT);
  _S::read_block(_S::leftblock, _start, LEFT);
  const Block<T>& site1 = _S::h.get_site(_start);
  const Block<T>& site2 = _S::h.get_site(_start+1);

  _S::_b1 = &_S::leftblock;
  _S::_b2 = &site1;
  _S::_b3 = &site2;
  _S::_b4 = &_S::rightblock;

  cout << "|" << _S::_b1->lattice().size() << "|-|" << _S::_b2->lattice().size() << "|-|" << _S::_b3->lattice().size() << "|-|" << _S::_b4->lattice().size() << "|" << endl;

  _S::hint(MASK_BLOCK1|MASK_BLOCK2);
  _S::hint(MASK_BLOCK3|MASK_BLOCK4);
  _S::hint(MASK_BLOCK2|MASK_BLOCK3, false);
*/
  return;
#endif

  _S::dir = _dir;
  _S::iter = _start; 
  cout << "MOVING TO THE END TO START TIME SLICES\n";

  char file[255];
  snprintf(file,255,"iter_%s.dat",_S::_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "Time sweep\n";

  _S::dir = _dir;
//****************************************************************
  cout << "START TIME SWEEP\n";
  if(_dir == LEFT2RIGHT){
    for(int _iter = _start; _iter < ls-1; _iter++)
    {
      _S::iter = _iter;
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
      int pos = _iter;
      if(_iter > ls -3) pos--;
      const Block<T>& site1 = _S::h.get_site(pos);
      const Block<T>& site2 = _S::h.get_site(pos+1);

      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();

      _S::_timer.Lap();
      cout << "START TIME SWEEP ITERATION " << endl;
      outputfile << "START TIME SWEEP ITERATION " << endl;
      cout << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;
      outputfile << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;
 
      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

      _S::gs = _S::_target[first_state];
      if(_S::iter < ls-3){
        _S::m = t;
        _S::m1 = std::min(_S::m1*_S::m2,_S::m);
  
        _S::truncate(LEFT, _S::m1);
        _S::rotate(LEFT, _S::newblock);
        _S::write_iter(LEFT);
      } else {
//        build_S::_target_states(0);
        _S::write_gs(_S::gs,_S::iter,LEFT);
        _S::write_gs(_S::gs,_S::iter+1,LEFT);
        _S::write_gs(_S::gs,1,RIGHT);
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  } else {
    for(int _iter = _start; _iter < ls-1; _iter++) 
    {
      _S::iter = _iter;
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
      int pos = _iter;
      if(_iter > ls-3) pos--;
      const Block<T>& site1 = _S::h.get_site(ls-2-pos);
      const Block<T>& site2 = _S::h.get_site(ls-2-pos+1);

      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();

      _S::_timer.Lap();
      cout << "START TIME SWEEP ITERATION " << endl;
      outputfile << "START TIME SWEEP ITERATION " << endl;
      cout << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;

      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

      _S::gs = _S::_target[first_state];
      if(_S::iter < ls-3){
        _S::m = t;
        _S::m4 = std::min(_S::m4*_S::m3,_S::m);

        _S::truncate(RIGHT, _S::m4);
        _S::rotate(RIGHT, _S::newblock);
        _S::write_iter(RIGHT);
      } else {
//        build_S::_target_states(0);
        _S::write_gs(_S::gs,_S::iter,RIGHT);
        _S::write_gs(_S::gs,_S::iter+1,RIGHT);
        _S::write_gs(_S::gs,1,LEFT);
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  }

  _in_empty_sweep = false;
  outputfile.close();
}

template<class T>
void
SystemThermo<T>::start_sweep(const Hami<T> *h, size_t t, size_t _dir, size_t _start)
{
  int ls = _S::lattice().ls();
  _in_empty_sweep = true;

  time_evolve = true;

  _S::iter = _start;
  _S::dir = _dir;
  _S::m = t;
  if(_S::iter == 1)
#ifdef USE_RK
    read_gs(_S::gs, ls-2-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#else
    read_gs(_S::gs, ls-1-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#endif
  else
    read_gs(_S::gs, _S::iter-1, _dir == LEFT2RIGHT ? LEFT : RIGHT);

  size_t first_state = _propagate_gs ? 1 : 0;
#ifdef USE_RK
  if(_use_rk){
    _S::_target.resize(NTARGETS+first_state);
    _S::_target_weight.resize(NTARGETS+first_state);
    _S::set_ntargets(NTARGETS+first_state);
    if(_propagate_gs && _evolve_gs){
      _S::_target.resize(2*NTARGETS+first_state-1);
      _S::_target_weight.resize(2*NTARGETS+first_state-1);
      _S::set_ntargets(2*NTARGETS+first_state-1);
    }
  }
#else
  _S::_target.resize(1+first_state);
  _S::_target_weight.resize(1+first_state);
  _S::set_ntargets(1+first_state);
#endif

  if(h){ // ONLY WORKS IN THE GRAND CANONICAL
    for(int i = 0; i < _S::_ntargets; i++) _S::_target[i] = _S::gs;
    typename Hami<T>::const_iterator iter;
    VectorState<T> _seed(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis(),_S::gs.qn());
    for(iter = h->begin(); iter != h->end(); iter++){
      Term<T> t = *iter;
      BasicOp<T> op = t[0];
      const BasicOp<T> *real_op = _S::operator()(op);
      if(!real_op) cout << "ERROR: Operator " << op.description() << endl;
      size_t _block = _S::block(real_op->site());
      cout << "APPLYING OPERATOR " << real_op->description() << " " << _block << endl;
      cout << "GS NORM = " << product(_S::gs,_S::gs) << endl;

      product(*real_op, _S::gs, _seed, _block);
      cout << _seed.qn().n() << " " << _seed.qn().sz() << endl;
      cout << "NORM = " << product(_seed,_seed) << endl;
//      _seed = _seed / sqrt(product(_seed,_seed));
    }
    _S::_target[first_state] = _seed;
    _S::gs = _S::_target[first_state];
    _S::measure();
  } else {
    _S::_target[first_state] = _S::gs;
  }
  if(_start == 1)
#ifdef USE_RK
    write_gs(_S::gs, ls-2-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#else
    write_gs(_S::gs, ls-1-_S::iter, _dir == LEFT2RIGHT ? RIGHT : LEFT);
#endif
  else
   _S::write_gs(_S::gs, _S::iter-1, _dir == LEFT2RIGHT ? LEFT : RIGHT);

/*
  if(_dir == LEFT2RIGHT){

    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
    const Block<T>& site1 = _S::h.get_site(_S::iter);
    const Block<T>& site2 = _S::h.get_site(_S::iter+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);

  } else {
   
    _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
    _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
    const Block<T>& site1 = _S::h.get_site(ls-2-_S::iter);
    const Block<T>& site2 = _S::h.get_site(ls-2-_S::iter+1);

    _S::m1 = _S::leftblock.dim();
    _S::m4 = _S::rightblock.dim();
    _S::m2 = site1.dim();
    _S::m3 = site2.dim();

    init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
  }
*/
/*
  _S::m = t;
  if(_S::dir == LEFT2RIGHT){
    _S::m1 = std::min(_S::m1*_S::m2,_S::m);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;
    _S::truncate(LEFT, _S::m1);
    _S::rotate(LEFT, _S::newblock);
    _S::write_iter(LEFT);
  } else {
    _S::m4 = std::min(_S::m4*_S::m3,_S::m);
    cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;
    _S::truncate(RIGHT, _S::m4);
    _S::rotate(RIGHT, _S::newblock);
    _S::write_iter(RIGHT);
  }
  _S::gs = _S::_target[0]; 
//  _S::measure_n(1, &_S::_target[first_state]);
//  _S::measure(&_S::_target[first_state]);
  _S::gs = _S::_target[first_state]; 
  _S::write_gs(_S::gs, _S::iter, _S::dir == LEFT2RIGHT ? LEFT : RIGHT);
*/

#ifdef USE_RK
/*
  _S::read_block(_S::rightblock, ls-_start-2, RIGHT);
  _S::read_block(_S::leftblock, _start, LEFT);
  const Block<T>& site1 = _S::h.get_site(_start);
  const Block<T>& site2 = _S::h.get_site(_start+1);

  _S::_b1 = &_S::leftblock;
  _S::_b2 = &site1;
  _S::_b3 = &site2;
  _S::_b4 = &_S::rightblock;

  cout << "|" << _S::_b1->lattice().size() << "|-|" << _S::_b2->lattice().size() << "|-|" << _S::_b3->lattice().size() << "|-|" << _S::_b4->lattice().size() << "|" << endl;

  _S::hint(MASK_BLOCK1|MASK_BLOCK2);
  _S::hint(MASK_BLOCK3|MASK_BLOCK4);
  _S::hint(MASK_BLOCK2|MASK_BLOCK3, false);
*/
  return;
#endif

  _S::dir = _dir;
  _S::iter = _start; 
  cout << "MOVING TO THE END TO START TIME SLICES\n";

  char file[255];
  snprintf(file,255,"iter_%s.dat",_S::_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cerr << "*** ERROR: could not open " << file << endl;
  outputfile << "Time sweep\n";

  _S::dir = _dir;
//****************************************************************
  cout << "START TIME SWEEP\n";
  if(_dir == LEFT2RIGHT){
    for(int _iter = _start; _iter < ls-1; _iter++)
    {
      _S::iter = _iter;
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
      int pos = _iter;
      if(_iter > ls -3) pos--;
      const Block<T>& site1 = _S::h.get_site(pos);
      const Block<T>& site2 = _S::h.get_site(pos+1);

      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();

      _S::_timer.Lap();
      cout << "START TIME SWEEP ITERATION " << endl;
      outputfile << "START TIME SWEEP ITERATION " << endl;
      cout << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;
      outputfile << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;
 
      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

      _S::gs = _S::_target[first_state];
      if(_S::iter < ls-3){
        _S::m = t;
        _S::m1 = std::min(_S::m1*_S::m2,_S::m);
  
        _S::truncate(LEFT, _S::m1);
        _S::rotate(LEFT, _S::newblock);
        _S::write_iter(LEFT);
      } else {
//        build_S::_target_states(0);
        _S::write_gs(_S::gs,_S::iter,LEFT);
        _S::write_gs(_S::gs,_S::iter+1,LEFT);
        _S::write_gs(_S::gs,1,RIGHT);
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  } else {
    for(int _iter = _start; _iter < ls-1; _iter++) 
    {
      _S::iter = _iter;
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
      int pos = _iter;
      if(_iter > ls-3) pos--;
      const Block<T>& site1 = _S::h.get_site(ls-2-pos);
      const Block<T>& site2 = _S::h.get_site(ls-2-pos+1);

      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();

      _S::_timer.Lap();
      cout << "START TIME SWEEP ITERATION " << endl;
      outputfile << "START TIME SWEEP ITERATION " << endl;
      cout << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;

      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

      _S::gs = _S::_target[first_state];
      if(_S::iter < ls-3){
        _S::m = t;
        _S::m4 = std::min(_S::m4*_S::m3,_S::m);

        _S::truncate(RIGHT, _S::m4);
        _S::rotate(RIGHT, _S::newblock);
        _S::write_iter(RIGHT);
      } else {
//        build_S::_target_states(0);
        _S::write_gs(_S::gs,_S::iter,RIGHT);
        _S::write_gs(_S::gs,_S::iter+1,RIGHT);
        _S::write_gs(_S::gs,1,LEFT);
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _S::_timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  }

  _in_empty_sweep = false;
  outputfile.close();
}

template<class T>
void
SystemThermo<T>::empty_sweep(size_t _dir, size_t _start, size_t _end, bool build_operators)
{
  int ls = _S::lattice().ls();
  size_t first_state = _propagate_gs ? 1 : 0;
  CTimer clock;
  _in_empty_sweep = true;

  clock.Start();

  _S::dir = _dir;
//****************************************************************
  if(_dir == LEFT2RIGHT){
    for(int _iter = _start; _iter <= _end; _iter++)
    {
      _S::iter = _iter;
      _S::signal_emit(SYSTEM_SIGNAL_START_ITER);
      clock.Lap();
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,_S::iter)), LEFT);
      int pos = _iter;
      if(_iter > ls -3) pos--;
      const Block<T>& site1 = _S::h.get_site(pos);
      const Block<T>& site2 = _S::h.get_site(pos+1);
      cout << "EMPTY SWEEP " << _S::iter << " " << ls-_S::iter-2 << endl;
      cout << "LEFT-TO-RIGHT ITERATION " << _S::iter << endl;

      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();

      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, true, build_operators);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;

      _S::gs = _S::_target[first_state];

      if(_S::iter != ls-3){
        _S::m1 = std::min(_S::m1*_S::m2,_S::m);
  
        if(build_operators) {
          _S::truncate(LEFT, _S::m1);
          _S::rotate(LEFT, _S::newblock);
          _S::write_iter(LEFT);
        } else {
          _S::write_gs(_S::gs,_S::iter,LEFT);
        }
      } else {
//        build_S::_target_states(0);
        _S::write_gs(_S::gs,_S::iter,LEFT);
        _S::write_gs(_S::gs,_S::iter+1,LEFT);
      }
      if (_S::signal_emit(SYSTEM_SIGNAL_MEASURE)) _S::measure();

      cout << "===========================================\n";
      cout << "CPU time (empty_sweep1): " << clock.LapTime() << endl;
      _S::signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
  } else {
    for(int _iter = _start; _iter <= _end; _iter++) 
    {   
      _S::iter = _iter;
      _S::signal_emit(SYSTEM_SIGNAL_START_ITER);
      clock.Lap();
      _S::read_block(_S::rightblock, std::min(ls-3,std::max(1,_S::iter)), RIGHT);
      _S::read_block(_S::leftblock, std::min(ls-3,std::max(1,ls-_S::iter-2)), LEFT);
      int pos = _iter;
      if(_iter > ls -3) pos--;
      const Block<T>& site1 = _S::h.get_site(ls-2-pos);
      const Block<T>& site2 = _S::h.get_site(ls-2-pos+1);
      cout << "EMPTY SWEEP " << _S::iter << " " << ls-_S::iter-2 << endl;
      cout << "RIGHT-TO-LEFT ITERATION " << _S::iter << endl;
  
      _S::m1 = _S::leftblock.dim();
      _S::m4 = _S::rightblock.dim();
      _S::m2 = site1.dim();
      _S::m3 = site2.dim();
  
      init_time_iteration(_S::leftblock, site1, site2, _S::rightblock, true, build_operators);
      cout << _S::m1 << " " << _S::m2 << " " << _S::m3 << " " << _S::m4 << endl;
 
      _S::gs = _S::_target[first_state];

      if(_S::iter != ls-3){
        _S::m4 = std::min(_S::m4*_S::m3,_S::m);
 
        if(build_operators){
          _S::truncate(RIGHT, _S::m4);
          _S::rotate(RIGHT, _S::newblock);
          _S::write_iter(RIGHT);
        } else {
          _S::write_gs(_S::gs,_S::iter,RIGHT);
        }
      } else {
//        build_target_states(0);
        _S::write_gs(_S::gs,_S::iter,RIGHT);
        _S::write_gs(_S::gs,_S::iter+1,RIGHT);
      }
      if (_S::signal_emit(SYSTEM_SIGNAL_MEASURE)) _S::measure();

      cout << "===========================================\n";
      cout << "CPU time (empty_sweep1): " << clock.LapTime() << endl;
      _S::signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
  }

  _S::signal_emit(SYSTEM_SIGNAL_END_SWEEP);
  cout << "CPU time(empty_sweep): " << clock.TotalTime() << endl;
  _in_empty_sweep = false;
}

/////////////////////////////////////////////////////////////////////////////
// init_time_iteration:
/////////////////////////////////////////////////////////////////////////////
template<class T>
void
SystemThermo<T>::init_time_iteration(const Block<T>&b1, const Block<T>&b2, const Block<T>&b3, const Block<T>&b4, bool rotate_states, bool build_operators)
{
  int ls = _S::lattice().ls();
  CTimer clock;

  clock.Start();

  _S::_b1 = &b1;
  _S::_b2 = &b2;
  _S::_b3 = &b3;
  _S::_b4 = &b4;

  cout << "|" << _S::_b1->lattice().size() << "|-|" << _S::_b2->lattice().size() << "|-|" << _S::_b3->lattice().size() << "|-|" << _S::_b4->lattice().size() << "|" << endl;

  int first_state = _propagate_gs ? 1 : 0;

//#ifdef USE_RK
#ifndef USE_PRODUCT_DEFAULT
  if(build_operators){
//    _S::create_interactions();
    _S::hint(MASK_BLOCK1|MASK_BLOCK2);
    _S::hint(MASK_BLOCK3|MASK_BLOCK4);
    if(_use_rk){
//      _S::hint(MASK_BLOCK2|MASK_BLOCK3, false);
      if(!_S::use_product_default()){
        _S::build_composite_operators();
        init_terms_composite(*this, this->_target[first_state]);
      }
    }

    cout << "CPU time (build interations) : " << clock.TotalTime() << endl;
  }
#endif
//#endif

  if(rotate_states){

    BMatrix<T> rho1;
    Basis basis1;
    BMatrix<T> rho2;
    Basis basis2;
  
    if(_S::iter <= 1 || _S::iter >= ls-2){
      if(_S::iter == 0)
        read_gs(_S::gs, ls-2-_S::iter, _S::dir == LEFT2RIGHT ? RIGHT : LEFT);
      else
#ifdef USE_RK
        read_gs(_S::gs, ls-2-_S::iter, _S::dir == LEFT2RIGHT ? RIGHT : LEFT);
#else
        if(_in_empty_sweep)
          read_gs(_S::gs, ls-1-_S::iter, _S::dir == LEFT2RIGHT ? RIGHT : LEFT);
        else 
          read_gs(_S::gs, _S::iter-1, _S::dir == LEFT2RIGHT ? LEFT : RIGHT);
//        _S::_target[first_state] = _S::gs;
#endif
    } else {
      if(_S::dir == LEFT2RIGHT){
        read_rho(rho1, basis1, _S::iter, LEFT);
        read_rho(rho2, basis2, _S::lattice().ls()-_S::iter-2+1, RIGHT);
        read_gs(_S::gs, _S::iter-1, LEFT);
//        _S::gs.resize(_S::_grand_canonical);
        _S::seed.set_qn_mask(_S::gs.qn(), _S::_grand_canonical);
        _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
        new_seed(_S::gs, _S::seed, rho1, rho2, basis1, basis2, LEFT);
        _S::gs = _S::seed;
        for(int i = 0; i < _S::_target.size(); i++) {
          cout << "NEW TARGET " << i << endl;
          _S::seed.set_qn_mask(_S::_target[i].qn(), _S::_grand_canonical);
          _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
//          _S::_target[i].resize(_S::_grand_canonical);
          new_seed(_S::_target[i], _S::seed, rho1, rho2, basis1, basis2, LEFT);
          _S::_target[i] = _S::seed;
        }
      } else {
        read_rho(rho1, basis1, _S::iter, RIGHT);
        read_rho(rho2, basis2, _S::lattice().ls()-_S::iter-2+1, LEFT);
        read_gs(_S::gs, _S::iter-1, RIGHT);
//        _S::gs.resize(_S::_grand_canonical);
        _S::seed.set_qn_mask(_S::gs.qn(), _S::_grand_canonical);
        _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
        new_seed(_S::gs, _S::seed, rho1, rho2, basis1, basis2, RIGHT);
        _S::gs = _S::seed;
        for(int i = 0; i < _S::_target.size(); i++) {
          cout << "NEW TARGET " << i << endl;
          _S::seed.set_qn_mask(_S::_target[i].qn(), _S::_grand_canonical);
          _S::seed.resize(_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis());
//          _S::_target[i].resize(_S::_grand_canonical);
          new_seed(_S::_target[i], _S::seed, rho1, rho2, basis1, basis2, RIGHT);
          _S::_target[i] = _S::seed;
        }
      }
    }
  }
  cout << "CPU time (init_time_iteration) : " << clock.TotalTime() << endl;
}

template<class T>
void
SystemThermo<T>::build_target_states(int pos)
{
  if(!_propagate_gs){
//    System<T>::build_target_states(pos);
    _S::_target[0] = _S::gs;
    _S::_target_weight[0] = 1.0;
  } else {
    _S::_target_weight = 0.5;
  }
#ifdef USE_RK
  if(_use_rk && time_evolve) {
    evolve(false);
    size_t first_state = _propagate_gs ? 1 : 0;
    _S::gs = _S::_target[first_state];
  }
  if(time_evolve){
    if(_propagate_gs){
      if(!_evolve_gs){
        _S::_target_weight[0] = 0.5;
        _S::_target_weight[1] = 0.25;
        _S::_target_weight[2] = _S::_target_weight[3] = _S::_target_weight[4] = 0.25/3;
      } else {
        _S::_target_weight[0] = 0.25;
        _S::_target_weight[1] = 0.25;
        for(int i = 2; i < 8; i++) _S::_target_weight[i] = 0.5/6.;
      }
    } else {
      _S::_target_weight[0] = 0.5;
      _S::_target_weight[1] = _S::_target_weight[2] = _S::_target_weight[3] = 0.5/3;
    }
  }
#endif // USE_RK
  if(_S::dm_ops.size() > 0) {
    _S::_target_weight[0] = 0.5;
    _S::_target.resize(_S::dm_ops.size()+1);
    _S::_target_weight.resize(_S::dm_ops.size()+1);
    typename Hami<T>::const_iterator titer;
    int nt = 1;
    for(titer = _S::dm_ops.begin(); titer != _S::dm_ops.end(); titer++, nt++){
      Term<T> t = *titer;
      BasicOp<T> top = t[0];

      const BasicOp<T> *_op = _S::operator()(top);
      if(_op){
        cout << "New target state: " << top.description() << "|gs>\n";
        _S::_target[nt].set_qn_mask(_S::qn+top.dqn, _S::_grand_canonical);
        _S::_target[nt] = VectorState<T> (_S::_b1->basis(),_S::_b2->basis(),_S::_b3->basis(),_S::_b4->basis(),_S::gs.qn()+top.dqn);
        product(*_op, _S::gs, _S::_target[nt], _S::block(top.site()));
        _S::_target[nt] /= sqrt(product(_S::_target[nt],_S::_target[nt]));

        _S::_target_weight[nt] = .5/_S::dm_ops.size();
      }
    }
  }
}

} //namspace dmtk

#endif // __DMTK_SYSTEM_TIME_H__
