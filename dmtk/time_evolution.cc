#include <dmtk/build_st.h>

namespace dmtk
{

bool do_measure;
double xtime;

void
calc_coef(Vector<double> &coef, int nexp, double dt)
{
  coef.resize(nexp+1);
  if(nexp == 9){
/* Fourth order expansion with 9 exps */
    double phi = 0.1720865590295143;
    double lambda = -0.09156203075515678;
    double xi = -0.1616217622107222;
    coef[1] = coef[9] = phi*dt;
    coef[2] = coef[8] = 0.5*(1.-2.*lambda)*dt;
    coef[3] = coef[7] = xi*dt;
    coef[4] = coef[6] = lambda*dt;
    coef[5] = (1.-2.*(xi+phi))*dt;
  } else if(nexp == 7) {
/* Fourth order expansion with 7 exps */
    double phi = 1.3512071919596578;
    coef[1] = coef[7] = 0.5*phi*dt;
    coef[2] = coef[6] = phi*dt;
    coef[3] = coef[5] = 0.5*(1.-phi)*dt;
    coef[4] = (1-2.*phi)*dt;
  } else if (nexp == 3){
/* Second order expansion */
    coef[1] = 0.5*dt;
    coef[2] = 1*dt;
    coef[3] = 0.5*dt;
  } else if (nexp <= 2){
    coef[1] = dt;
    coef[2] = dt;
  }

  for(int i = 1; i < coef.size(); i++) cout << "COEF " << i << " = " << coef(i) << endl;
}


template<class T>
bool handler(System<T> &S, size_t signal_id, void *data)
{
  if(do_measure && signal_id == SYSTEM_SIGNAL_MEASURE){

#ifndef USE_RK
    if(S._target.size() > 1)
      S._target[0] *= exp(complex<double>(0.,-xtime*S.energy[0]));
#endif
    S.gs = S._target[0];
    if(S._target.size() > 1){
      S.measure_n(1, &S._target[1]);
      S.gs = S._target[1];
    } else {
      S.measure();
    }
  }

  return false;
}


template<class T>
void
init_time_evolution(SystemThermo<T> &S, BasicOp<T> *_op1 = NULL, BasicOp<T> *_op2 = NULL, bool propagate_gs = true, bool clear_corr = true)
{
  int lx = S.h.lattice().size();

  S.set_imaginary_time(false);
  if(propagate_gs && !(_op1 == NULL)) S.propagate_gs(S.gs,false);
  S.phase = S.energy[0];
  S.dm_ops.clear();
  S.set_store_products(false);

  const BasicOp<T> *real_op1 = NULL;
  if(_op1){
    BasicOp<T> op1 = *_op1;
    if(!S(op1)) cout << "ERROR: OPERATOR NOT FOUND\n";
    real_op1 = S(op1);
    cout << "APPLYING OPERATOR " << op1.description() << endl;

    S.measure();
// Signal handler for measurements
    S.signal_connect(handler, SYSTEM_SIGNAL_MEASURE, NULL);
    do_measure = false;
  }

  if(clear_corr) S.corr.clear();
  if(_op2){
    BasicOp<T> op2 = *_op2;
    for(int i = 0; i < lx; i++){
      S.corr += op2.set_site(i) * T(1.);
    }
  }

#ifdef USE_RK
#ifndef USE_EXP_RK
  S.set_use_product_default(true);
#endif
#endif

  if(real_op1)
    S.start_sweep(real_op1, S.get_truncation(), LEFT2RIGHT, lx/2, S.block(_op1->site()));
  else
    S.start_sweep(real_op1, S.get_truncation(), LEFT2RIGHT, lx/2, BLOCK2);
#ifdef USE_RK
  S.empty_sweep(LEFT2RIGHT, lx/2, lx-3);
#endif // USE_RK
}

template<class T>
void
init_time_evolution(SystemThermo<T> &S, BasicOp<T> *_op1 = NULL, Hami<T> *new_corr = NULL, bool propagate_gs = true, bool clear_corr = true)
{
  int lx = S.h.lattice().size();

  S.set_imaginary_time(false);
  if(propagate_gs && !(_op1 == NULL)) S.propagate_gs(S.gs,false);
  S.phase = S.energy[0];
  S.dm_ops.clear();
  S.set_store_products(false);

  const BasicOp<T> *real_op1 = NULL;
  if(_op1){
    BasicOp<T> op1 = *_op1;
//    if(!S(op1.set_site(lx/2))) cout << "ERROR: OPERATOR NOT FOUND\n";
    if(!S(op1)) cout << "ERROR: OPERATOR NOT FOUND\n";
    real_op1 = S(op1);
    cout << "APPLYING OPERATOR " << op1.description() << endl;

    S.measure();
// Signal handler for measurements
    S.signal_connect(handler, SYSTEM_SIGNAL_MEASURE, NULL);
    do_measure = false;
  }

  if(clear_corr) S.corr.clear();
  if(new_corr) S.corr = *new_corr;

#ifdef USE_RK
#ifndef USE_EXP_RK
  S.set_use_product_default(true);
#endif
#endif

  if(real_op1)
    S.start_sweep(real_op1, S.get_truncation(), LEFT2RIGHT, lx/2, S.block(_op1->site()));
  else
    S.start_sweep(real_op1, S.get_truncation(), LEFT2RIGHT, lx/2, BLOCK2);

#ifdef USE_RK
  S.empty_sweep(LEFT2RIGHT, lx/2, lx-3);
#endif // USE_RK
}

template<class T>
void
init_time_evolution(SystemThermo<T> &S, Hami<T> *h = NULL, Hami<T> *new_corr = NULL, bool propagate_gs = true, bool clear_corr = true)
{
  int lx = S.h.lattice().size();

  S.set_imaginary_time(false);
  if(propagate_gs && !(h == NULL)) S.propagate_gs(S.gs,false);
  S.phase = S.energy[0];
  S.dm_ops.clear();
  S.set_store_products(false);

  S.signal_connect(handler, SYSTEM_SIGNAL_MEASURE, NULL);
  do_measure = false;

  if(clear_corr) S.corr.clear();
  if(new_corr) S.corr = *new_corr;

#ifdef USE_RK
#ifndef USE_EXP_RK
  S.set_use_product_default(true);
#endif
#endif
  S.start_sweep(h, S.get_truncation(), LEFT2RIGHT, lx/2);

#ifdef USE_RK
  S.empty_sweep(LEFT2RIGHT, lx/2, lx-3);
#endif // USE_RK
}



template<class T>
void
time_evolution(SystemThermo<T> &S, const Hami<T> &st_hami, double tstep, double tend, int nexp)
{
  VectorState<T> aux;
  int lx = S.h.lattice().size();
  std::vector<BasicOp<T> > hij(lx-1);
  std::vector<BasicOp<T> > hij_identity(lx-1);
  Vector<double> coef;
  S.dt = tstep;

// Build default identity operators
  hij_identity = default_st_ops(st_hami);
#warning If you are including the h.c. operator in the Hamiltonian, make user st_hami is set use_hc = true 

/*
  hij = default_st_ops(st_hami);
  for(int i = 0; i < lx-1; i++) {
     S.hsites[i] = &hij[i];
  }
*/

// Calculate coefs. for the s-t expansion
  calc_coef(coef, nexp, tstep);

// Start time evolution
  xtime = 0;
  while(xtime < tend){

    xtime += tstep;
    if(S.imaginary_time())
      cout << "TIME SWEEP BETA= " << xtime << " " << tstep << " " << nexp << endl;
    else
      cout << "TIME SWEEP T= " << xtime << " " << tstep << " " << nexp << endl;
#ifndef USE_RK

// Suzuki-Trotter
    for(int term = 1; term <= nexp; term++){
      double xtstep = coef[term];
      hij = build_st_ops(st_hami, xtstep, S.imaginary_time());

      if(nexp > 1){
        int is0 = term%2;
        for(int i = is0; i < lx-1; i+=2) S.hsites[i] = &hij[i];
        for(int i = 1-is0; i < lx-1; i+=2) S.hsites[i] = &hij_identity[i];
      } else {
        for(int i = 0; i < lx-1; i++) S.hsites[i] = &hij[i];
      }

      if(S.get_dir() == RIGHT2LEFT){
          cout << "NEW TIME SWEEP: LEFT2RIGHT " << S.get_truncation() << endl;
          S.time_sweep(S.get_truncation(),LEFT2RIGHT,0,true);
      }
      else
      {
          cout << "NEW TIME SWEEP: RIGHT2LEFT " << S.get_truncation() << endl;
          S.time_sweep(S.get_truncation(),RIGHT2LEFT,0,true);
      }
    }
    S.my_time += tstep;
#else // USE_RK                                                                                                                                       
// Runge-Kutta 
    if(S.get_dir() == RIGHT2LEFT) {
      cout << "NEW TIME SWEEP: LEFT2RIGHT " << S.get_truncation() << endl;
      S.empty_sweep(LEFT2RIGHT, 1, lx-3);
    } else {
      cout << "NEW TIME SWEEP: RIGHT2LEFT " << S.get_truncation() << endl;
      S.empty_sweep(RIGHT2LEFT, 1, lx-3);
    }
    S.evolve();

#endif // USE_RK          

// Run an empty sweep and measure correlations
    do_measure = true;
    bool build_operators = true;
#ifdef USE_RK
    build_operators = true;
#endif
    S.empty_sweep((S.get_dir() == RIGHT2LEFT ? LEFT2RIGHT : RIGHT2LEFT), 1, lx-3, build_operators);
    do_measure = false;

    if(S.imaginary_time() && S.rotate_terms()){
      S.gs = S._target[0];
      aux = product_default(S, S.gs, &S.h);
      T xe = product(S.gs,aux);
      cout << "ENER BETA = " << xtime << " " << 0.5/xtime << " " << real(xe)/(real(product(S.gs,S.gs))) << " " << log(real(product(S.gs,S.gs))) << endl;
      S.measure();
    }

    typename Hami<T>::iterator corr_iter;
    for(corr_iter = S.corr.begin(); corr_iter != S.corr.end(); corr_iter++){
      cout << "RESULTS " << xtime << " " << corr_iter->description() << " " << corr_iter->value() << endl;
    }

  }
}

} //namespace dmtk
