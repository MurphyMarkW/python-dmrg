#ifndef __DMTK_SYSTEM_H__
#define __DMTK_SYSTEM_H__
#warning You are using dmtk-1.6-devel-20100514
//*******************************************************************
// System Class: Defines system properties and performs DMRG loops
// This version of the system class works only on symmetric systems
//*******************************************************************
 
#include <iostream>
#include <iosfwd>
#include <list>
#include <string>
#include <map>
#include <iomanip>
#include <math.h>
#include "enums.h"
#include "vector.h"
#include "matrix.h"
#include "state.h"
#include "operators.h"
#include "block.h"
#include "lattice.h"
#include "hami.h"
#include "ctimer.h"
#include "util.h"
#include "globals.h"
#include "product.h"
#include "signals.h"
#ifdef WITH_PTHREADS
#include <pthread.h>
#include <unistd.h>
#endif // WITH_PTHREADS
#ifdef WITH_BZIP2
#include "bzip2stream.h"
#endif
#include "virtualdisk.h"

using namespace std;

namespace dmtk
{

enum{
 SYSTEM_SIGNAL_GS,  
 SYSTEM_SIGNAL_TRUNCATE,  
 SYSTEM_SIGNAL_BUILD_DM,  
 SYSTEM_SIGNAL_DM_READY,  
 SYSTEM_SIGNAL_ROTATE,  
 SYSTEM_SIGNAL_START_ITER,  
 SYSTEM_SIGNAL_END_ITER,  
 SYSTEM_SIGNAL_MEASURE,  
 SYSTEM_SIGNAL_END_SWEEP,  
 SYSTEM_SIGNAL_LAST,
};

#ifndef NUM_TRHEADS
//#define NUM_THREADS (sysconf(_SC_NPROCESSORS_ONLN)*2)
#define NUM_THREADS 40
#endif

#define COUT_PRODUCT_DEFAULT(term) \
   cout << "PRODUCT " << term.description() << endl;

#define COUT_PRODUCT(term) \
 if(calc_hc) \
   cout << "PRODUCT " << term.description() << " + h.c." << endl; \
 else  \
   cout << "PRODUCT " << term.description() << endl;

template<class T> 
VectorState<T> product(System<T>& ss, const VectorState<T>& vv);

template<class T> 
VectorState<T> product_custom(System<T>& ss, const VectorState<T>& vv);

template<class T> 
VectorState<T> product_composite(System<T>& ss, const VectorState<T>& vv, int mask_hc = (MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC));

template<class T> 
VectorState<T> product_default(System<T>& ss, const VectorState<T>& vv, const Hami<T>* hami = NULL, bool only_local = false);

template<class T>
VectorState<T> product_additive(System<T>& ss, const VectorState<T>& vv, const Term<T>& term);

template<class T>
class AuxTerm{
  public:
    Term<T> t;
    Term<T> cterm[4]; // composite term
    BasicOp<T> top[4]; // reference operator
    const BasicOp<T> *op[4]; // original operator
    const BasicOp<T> *ref_op[4]; // complementary operator
    BasicOp<T> sum_op; // complementary operator
    size_t pos[4]; // original position in the term
    size_t b1, b2, b3;
    int mask;
    bool done;
    int nblocks;

    AuxTerm(): nblocks(0), mask(0), done(false), b1(BLOCK_NONE), b2(BLOCK_NONE), b3(BLOCK_NONE) { op[0] = op[1] = op[2] = op[3] = NULL; ref_op[0] = ref_op[1] = ref_op[2] = ref_op[3] = NULL; }
};

template<class T>
void product_term(const System<T>& ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res);
template<class T> 
void init_terms_composite(System<T>& ss, const VectorState<T>& vv);
template<class T>
void init_term_composite(System<T>& ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res);

#include "hmatrix.h"

VirtualDisk disk;

template<class T>
class System
{
  protected:
    typedef typename dmtk::Block<T> B;

    Lattice _lattice;

    CTimer _timer;

    bool _in_warmup;
    bool _warmup_growing;
    size_t dir;
    int _sweep;
    int iter;

    int _grand_canonical;

    int m;
    int m1;
    int m2;
    int m3;
    int m4;

    double _error;
    bool _use_error;
    int _error_max_size;
    bool _target_subspaces; // build density matrix using all subspaces
    int _nsub; // number states per subspace

    int _verbose;
    int _calc_gap;
    bool _project;
    bool _use_hamis;
    bool _use_hc;
    bool _is_hermitian;
    bool _use_k;
    bool _use_seed;
    bool _use_hamis_seed;
    bool _use_basic_seed;
    bool _use_composite;
    bool _use_energy_rg;
    bool _grow_symmetric;
    bool _grow_outward;
    bool _custom_qns;
    bool _store_products;
    bool _measure_symmetric;
    bool _apply_hami;
    bool _apply_extern;
    bool _full_sweep; // run sweep from 1 to ls-3/ls-4
    bool _use_coef_tol; 
    bool _use_single_site;
    bool _save_hmatrix;
    bool _rotate_terms;
    bool _rotate_corr;
    bool _use_product_default;
    double _coef_tol;
    double _random_coef;
    Matrix<size_t> _sweeps;
    size_t _numsweeps;
    bool _print_dm_spectrum;

    std::map<string,T> _param;

    char _name[255];
    Vector<char> _stream_buffer;

    double _lanczos_tol;
    int _lanczos_maxiter;
    bool _use_lanczos;

    virtual void init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed = false, bool create_composite = true);
    void rotate_hami(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);
    void rotate_terms(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);
    void rotate_corr(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Hami<T> *this_hami = NULL);
    void rotate_additive(int position, Block<T> &b, Basis &basis, Basis &rho_basis, const Term<T> &term);
    void real_warmup_loop(size_t t, const Vector<QN> &qns, size_t first_iter = 1, bool randomize = false);

    void hint(int _mask, bool add_local_h = true, const Hami<T> *_h = NULL);
    void build_composite_operators(const Hami<T> *_h = NULL);

    typename std::list<Signal<System<T> > > handlers;

  public:

    std::vector<ProductTerm<T> > product_terms;

    B rightblock;
    B leftblock;
    B newblock;

    const B* _b1;
    const B* _b2;
    const B* _b3;
    const B* _b4;

    QN qnt;
    QN qn;
    double precision;

    size_t _ntargets;
    Vector<VectorState<T> > _target;   
    Vector<double> _target_weight;
    Vector<VectorState<T> *> _project_states;

    size_t _nstates;  // excited states for diagonalization
    Vector<VectorState<T> > _state; // gs and excited states  
    Vector<VectorState<T> > _propagate_state; // states that we want to transform as we sweep (besides the ground state)

    BMatrix<T> rho;
    Basis rho_basis;
    BMatrix<T> rho_left;
    Vector<double> rho_w;
    VectorState<T> gs;
    VectorState<T> seed;
    Vector<double> energy;

    Hami<T> h;
    std::vector<Hami<T> > hs;
    H<T> h12;
    H<T> h23;
    H<T> h34;
    std::list<AuxTerm<T> > aux_terms;

    Hami<T> dm_ops;
    Hami<T> dm_boundary_ops;
    Hami<T> ops;
    Hami<T> corr;
    Hami<T> control_ops; // single operators

//  Constructor
  
    System(): _numsweeps(0), _grow_outward(false), _nsub(1), _target_subspaces(false), _use_energy_rg(false), _use_composite(true), _custom_qns(true), _in_warmup(true), _sweep(1), iter(1), _verbose(0), _store_products(true), _use_hamis(false), _project (false), _calc_gap(0), _use_hc(false), _is_hermitian(true), _grow_symmetric(true), _measure_symmetric(false), _full_sweep(false), m1(2), m2(2), m4(2), m(10), _use_lanczos(true), _lanczos_tol(1.e-7), _lanczos_maxiter(-1), _nstates(1), _ntargets(1), _use_k(false), _grand_canonical(QN::default_mask()), _use_error(false), _error_max_size(-1), _use_seed(true), _use_hamis_seed(true), _use_basic_seed(false), _apply_hami(true), _apply_extern(true), _use_coef_tol(false), _coef_tol(1.e-5), _use_single_site(false), _save_hmatrix(false), _print_dm_spectrum(false), _rotate_terms(true), _rotate_corr(true), _use_product_default(false)
    { 
      _sweeps.resize(2,10); 
      _target.resize(1); 
      _target_weight.resize(1); 
      set_name(""); 
      set_stream_buffer_size(1024);
      _random_coef = 0.000001;
    };

    System(const Hami<T> &_h, const Lattice& lattice, const char *the_name): _numsweeps(0), _grow_outward(false), _nsub(1), _target_subspaces(false), _use_energy_rg(false), _use_composite(true), h(_h), _lattice(lattice), _custom_qns(true), _in_warmup(true), _sweep(1), iter(1), _verbose(0), _store_products(true), _use_hamis(false), _project (false), _calc_gap(0), _use_hc(false), _is_hermitian(true), _grow_symmetric(true), _measure_symmetric(false), _full_sweep(false), _nstates(1), _ntargets(1), _use_lanczos(true), _lanczos_tol(1.e-7), _lanczos_maxiter(-1), qnt(_h.lattice().ls(), 0), _use_k(false), _grand_canonical(QN::default_mask()), _use_error(false), _error_max_size(-1), _use_seed(true), _use_hamis_seed(true), _use_basic_seed(false), _apply_hami(true), _apply_extern(true), _use_coef_tol(false), _coef_tol(1.e-5), _use_single_site(false), _save_hmatrix(false), _print_dm_spectrum(false), _rotate_terms(true), _rotate_corr(true), _use_product_default(false)
    { 
      set_name(the_name); 
      m1 = m2 = m3 = m4 = _h.get_site(0).dim(); 
      m = m1 * m2; _sweeps.resize(2,10); 
      _target.resize(1); 
      _target_weight.resize(1); 
      _target_weight[0] = double(1); 
      set_stream_buffer_size(1024);
      _random_coef = 0.000001;
    }

//  Methods
    void set_parameter(const std::string &s, T val)
      { _param[s] = val; }
    T get_parameter(const std::string &s) const
      { return _param[s]; }
    int nparam() const { return _param.size(); } 

    void start(bool reorder = true) 
      { 
        if(reorder)
          h.reorder_terms(true);
        write_block(h.get_site(0),1,LEFT); 
        write_block(h.get_site(_lattice.size()-1),1,RIGHT); 
      }

    void start(size_t numsweeps, const Matrix<size_t> &nstates)
      { 
        _numsweeps = numsweeps;
        _sweeps = nstates;
        start();
      }

    void resume(size_t resume_sweep, size_t resume_dir, size_t resume_iter)
      {
        for(int _sweep = resume_sweep; _sweep <= _numsweeps; _sweep++){
          cout << "NEW SWEEP " << _sweep << " " << _sweeps(0,_sweep-1) << " " << _sweeps(1,_sweep-1) << endl;
          sweep(_sweeps(0,_sweep-1),_sweeps(1,_sweep-1),(size_t)resume_dir,resume_iter);
        }
      }

    void resume_default()
      {
        read_status();

        int resume_sweep = _sweep;
        size_t resume_dir = dir;
        int resume_iter = iter;

        resume(resume_sweep, resume_dir, resume_iter+1);
      }

    void run(int nwarmup = 20)
      {
        if(nwarmup > 0) warmup_loop(nwarmup);
        resume(1, RIGHT2LEFT, 1);
      }

    void set_error(double err, int error_max_size = -1) { _error = err; _error_max_size = error_max_size; _use_error = true; }

    void set_truncation(size_t _m) { m = _m; }
    size_t get_truncation () const { return m; }
    bool in_warmup() const { return _in_warmup; }
    void set_stream_buffer_size(size_t m) 
      { 
        _stream_buffer.resize(m); 
//        setvbuf(stdout, _stream_buffer.array(), _IOFBF, m);
      }

    void thermo_field(size_t state = 0, Matrix<size_t> *states = NULL, Vector<double> *weights = NULL, size_t ntrunc = 0, Vector<T> *gs_coef = NULL);

    void set_name(const char *the_name) 
      { snprintf(_name, 255, "%s", the_name); }
    const char* name() const { return _name; }

    void create_interactions(bool create_composite = true);
    BMatrix<T> density_matrix(int position);
    void delta_rho(void);
    virtual void diagonalize(bool use_seed = true);
    virtual void build_target_states(int pos = 0);
    virtual void truncate(int block, int new_size, bool build_targets = true, bool diagonalize = true);
    virtual void rotate(int block, Block<T>& b); 
    virtual void warmup_loop(size_t t, size_t first_iter = 1, bool randomize = false);
    virtual void warmup_loop(size_t t, const Vector<QN> &qns, size_t first_iter = 1, bool randomize = false);
    virtual void sweep(size_t t1, size_t t2, size_t dir = RIGHT2LEFT, int _start = 1);
    virtual void main_loop(size_t nsweeps, size_t t);
    virtual void init_correlations();
    virtual void final_sweep(size_t t, size_t dir = RIGHT2LEFT, int _start = 1, bool _rotate = false);
    virtual void empty_sweep(size_t t, size_t dir = RIGHT2LEFT, int _start = 1, bool half_sweep = false);
    virtual void measure(const VectorState<T> *v = NULL); 
    virtual void measure_n(size_t n = 0, const VectorState<T> *v = NULL); 
    virtual T measure_operator(const BasicOp<T> &top, const VectorState<T> *v, bool &info);
    virtual VectorState<T> product_operator(const BasicOp<T> &top, const VectorState<T> *v = NULL);
    virtual VectorState<T> product_term(const Term<T> &t, const VectorState<T> *v = NULL);

    size_t block(int site) const;

    const BasicOp<T>* operator()(const BasicOp<T>& op) const;
    const Block<T>& operator[](size_t pos) const 
      { 
        const Block<T> *b;
        switch(pos){
          case BLOCK1: b = _b1; break;
          case BLOCK2: b = _b2; break;
          case BLOCK3: b = _b3; break;
          case BLOCK4: b = _b4; break;
        }
        return *b;
      }


    size_t size() const 
      { return _b1->lattice().size()+ _b2->lattice().size()+
               _b3->lattice().size()+ _b4->lattice().size(); }
    int site1() const 
      { if(dir == LEFT2RIGHT) return iter; else return h.lattice().size()-iter-2; }
    int site2() const { return site1()+1; }

    const Lattice& lattice() const { return _lattice; } 

    System<T>& set_ntargets(size_t ntargets)
      { _ntargets = ntargets >= 1 ? ntargets : 1; return *this; } 
    size_t ntargets() { return _ntargets; }

    System<T>& set_target_weights(const Vector<double>& w)
      { _target_weight = w; return *this; }
    

    System<T>& set_store_products(bool b) { _store_products = b; return *this; }
    bool store_products() const { return _store_products; }
    System<T>& set_target_subspaces(bool b, int nsub = 1) { _target_subspaces = b; _nsub = nsub; return *this; }
    bool target_subspaces() const { return _target_subspaces; }
    System<T>& set_verbose(int b) { _verbose = b; return *this; }
    System<T>& set_verbose(bool b) { _verbose = b ? 1 : 0; return *this; }
    int verbose() const { return _verbose; }
    System<T>& set_translations(bool b) { _use_k = b; return *this; }
    bool translations() const { return _use_k; }
    System<T>& set_use_hc(bool b) { _use_hc = b; h.set_use_hc(b); return *this; }
    bool use_hc() const { return _use_hc; }
    System<T>& set_is_hermitian(bool b) { _is_hermitian = b; return *this; }
    bool is_hermitian() const { return _is_hermitian; }
    System<T>& set_calc_gap(int b) 
      { 
        _calc_gap = b; 
        if(b > 0){
          _ntargets = b+1;
          _target.resize(b+1);
          _target_weight.resize(b+1);
          _target_weight[0] = 0.875;
          for(int i = 0; i < b; i++) _target_weight[i+1] = 0.125/b;
        } else {
          _ntargets = 1;
          _target.resize(1);
          _target_weight.resize(1);
          _target_weight[0] = 1.0;
        }
        return *this; 
      }
    int calc_gap() const { return _calc_gap; }
    bool project() const { return _project; }
    bool use_hamis() const { return _use_hamis; }
    System<T>& set_use_hamis(bool b) 
      { 
        _use_hamis = b; 
        if(b){
          _ntargets = hs.size()+1;
          _target.resize(_ntargets);
          _target_weight.resize(_ntargets);
          _target_weight[0] = .5;
          for(int i = 1; i < _ntargets; i++) _target_weight[i] = 0.5/hs.size();
        } else {
          _ntargets = 1;
          _target.resize(1);
          _target_weight.resize(1);
          _target_weight[0] = 1.0;
        }
        return *this; 
      }
    System<T>& set_use_seed(bool b) { _use_seed = b; return *this; }
    bool use_seed() const { return _use_seed; }
    System<T>& set_use_hamis_seed(bool b) { _use_hamis_seed = b; return *this; }
    bool use_hamis_seed() const { return _use_hamis_seed; }
    System<T>& set_use_basic_seed(bool b) { _use_basic_seed = b; return *this; }
    bool use_basic_seed() const { return _use_basic_seed; }
    System<T>& set_apply_hami(bool b) { _apply_hami = b; return *this; }
    bool apply_hami() const { return _apply_hami; }
    System<T>& set_apply_extern(bool b) { _apply_extern = b; return *this; }
    bool apply_extern() const { return _apply_extern; }
    System<T>& set_use_coef_tol(bool b, double tol=1.e-5) { _use_coef_tol = b; _coef_tol = tol; return *this; }
    bool use_coef_tol() const { return _use_coef_tol; }
    double coef_tol() const { return _coef_tol; }
    System<T>& set_random_coef(double v) { _random_coef = v; return *this; }
    System<T>& set_use_single_site(bool b) { _use_single_site = b; return *this; }
    bool use_single_site() const { return _use_single_site; }
    System<T>& set_save_hmatrix(bool b) { _save_hmatrix = b; return *this; }
    bool save_hmatrix() const { return _save_hmatrix; }
    System<T>& set_grow_symmetric(bool b) { _grow_symmetric = b; return *this; }
    bool grow_symmetric() const { return _grow_symmetric; }
    System<T>& set_grow_outward(bool b) { _grow_outward = b; return *this; }
    bool grow_outward() const { return _grow_outward; }
    System<T>& set_use_composite(bool b) { _use_composite = b; return *this; }
    bool use_composite() const { return _use_composite; }
    System<T>& set_use_energy_rg(bool b) { _use_energy_rg = b; return *this; }
    bool use_energy_rg() const { return _use_energy_rg; }
    System<T>& set_measure_symmetric(bool b) { _measure_symmetric = b; return *this; }
    bool measure_symmetric() const { return _measure_symmetric; }
    System<T>& set_grand_canonical(bool b) { _grand_canonical = b ? 0 : QN::default_mask(); return *this; }
    bool grand_canonical() const { return (_grand_canonical == 0); }
    System<T>& set_qn_mask(int m = QN::default_mask()) { _grand_canonical = m; return *this; }
    int qn_mask() const { return _grand_canonical; }
    System<T>& set_full_sweep(bool b) { _full_sweep = b; return *this; }
    bool full_sweep() const { return _full_sweep; }
    CTimer timer() const { return _timer; }
    size_t sweep_direction() const { return dir; }
    int iteration() const { return iter; }
    void set_iteration(size_t _dir, size_t _iter, bool _rotate = false, int new_size = 20, bool use_seed = true);

    void set_lanczos_tolerance(double tol) { _lanczos_tol = tol; }
    double lanczos_tolerance() { return _lanczos_tol; }
    void set_lanczos_maxiter(int n) { _lanczos_maxiter = n; }
    int lanczos_maxiter() { return _lanczos_maxiter; }
    void set_use_lanczos(bool use) { _use_lanczos = use; }
    bool use_lanczos() const { return _use_lanczos; }
    void set_print_dm_spectrum(bool use) { _print_dm_spectrum = use; }
    bool print_dm_spectrum() const { return _print_dm_spectrum; }
    void set_rotate_terms(bool b) { _rotate_terms = b; }
    bool rotate_terms() { return _rotate_terms; }
    void set_rotate_corr(bool b) { _rotate_corr = b; }
    bool rotate_corr() { return _rotate_corr; }
    void set_use_product_default(bool b) { _use_product_default = b; }
    bool use_product_default() { return _use_product_default; }

    size_t get_iter() const { return iter; }
    size_t get_dir() const { return dir; }

    // Streams

    virtual void write_status() const
    {
      char file[255];
      snprintf(file, 255, "system_%s.dat", _name);
#ifdef USE_VDISK
      std::ostringstream outputfile(disk.oopen(file),std::ios::out|std::ios::binary);
#else
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
#endif

      outputfile.write((const char *)&_numsweeps, sizeof(size_t));
      for(int i = 1; i <= _numsweeps; i++){
        int l = _sweeps(0,i);
        outputfile.write((const char *)&l, sizeof(int));
        l = _sweeps(1,i);
        outputfile.write((const char *)&l, sizeof(int));
      }
      outputfile.write((const char *)&_in_warmup, sizeof(bool));
      outputfile.write((const char *)&_sweep, sizeof(int));
      outputfile.write((const char *)&dir, sizeof(size_t));
      outputfile.write((const char *)&iter, sizeof(int));
//      outputfile.write((const char *)&_ntargets, sizeof(int));
      int n = _propagate_state.size();
      outputfile.write((const char *)&n, sizeof(int));
      n = _target.size();
      outputfile.write((const char *)&n, sizeof(int));
      qnt.write(outputfile);
      qn.write(outputfile);

#ifndef USE_VDISK
      outputfile.close();
#else
      outputfile << flush;
      disk.close(outputfile.str());
#endif
    }

    virtual void read_status() 
    {
      char file[255];
      snprintf(file, 255, "system_%s.dat", _name);
#ifdef USE_VDISK
      std::istringstream inputfile(disk.iopen(file),std::ios::in|std::ios::binary);
#else
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
#endif
      inputfile.read((char *)&_numsweeps, sizeof(size_t));
      _sweeps.resize(2,_numsweeps+1);
      for(int i = 1; i <= _numsweeps; i++){
        inputfile.read((char *)&_sweeps(0,i), sizeof(int));
        inputfile.read((char *)&_sweeps(1,i), sizeof(int));
      }
      inputfile.read((char *)&_in_warmup, sizeof(bool));
      inputfile.read((char *)&_sweep, sizeof(int));
      inputfile.read((char *)&dir, sizeof(size_t));
      inputfile.read((char *)&iter, sizeof(int));
//      inputfile.read((char *)&_ntargets, sizeof(int));
      inputfile.read((char *)&_nstates, sizeof(int));
      int n;
      inputfile.read((char *)&n, sizeof(int));
      _propagate_state.resize(n);
      inputfile.read((char *)&n, sizeof(int));
      _target.resize(n);
      _target_weight.resize(n);
      qnt.read(inputfile);
      qn.read(inputfile);

#ifndef USE_VDISK
      inputfile.close();
#endif

      read_gs(gs,iter,dir == LEFT2RIGHT ? LEFT : RIGHT);
    }

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
      rho.write(s);
      rho_basis.write(s);
//      if(!is_hermitian()) rho_left.write(s);
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
      rho.read(s);
      rho_basis.read(s);
//      if(!is_hermitian()) rho_left.read(s);
      rightblock.read(s);
      leftblock.read(s);
    }

    void write_gs(const VectorState<T> &_gs, int n, int position = RIGHT) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file,255,"gs_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"gs_%s_%i_r.dat",_name,n);

#ifdef USE_VDISK
      std::ostringstream outputfile(disk.oopen(file),std::ios::out|std::ios::binary);
#else
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
#endif
      cout << "WRITING GS " << file << endl;
      VectorState<T> aux(_gs);
// CHANGE HERE IF YOU DO NOT WANT TO NORMALIZE THE GS
      aux /= sqrt(product(_gs,_gs));
      aux.write(outputfile);

      int _n = _target.size();
      outputfile.write((const char *)&_n, sizeof(int));
      for(int i = 0; i < _n; i++){
        aux = _target[i];
        aux /= sqrt(product(aux,aux));
        aux.write(outputfile);
      }

      _n = _propagate_state.size();
      outputfile.write((const char *)&_n, sizeof(int));
      for(int i = 0; i < _propagate_state.size(); i++){
        aux = _propagate_state[i];
        aux /= sqrt(product(aux,aux));
        aux.write(outputfile);
      }
#ifndef USE_VDISK
      outputfile.close();
#else
      outputfile << flush;
      disk.close(outputfile.str());
#endif
    }

    void read_gs(VectorState<T> &_gs, int n, int position = RIGHT)
    {
      char file[255];
      if(position == LEFT)
        snprintf(file,255,"gs_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"gs_%s_%i_r.dat",_name,n);

#ifdef USE_VDISK
      std::istringstream inputfile(disk.iopen(file),std::ios::in|std::ios::binary);
#else
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
#endif
      cout << "READING GS " << file << endl;
      _gs.read(inputfile);
      gs = _gs;

      inputfile.read((char *)&n, sizeof(int));
      _target.resize(n);
      _target_weight.resize(n);
      for(int i = 0; i < n; i++){
        _target[i].read(inputfile);
      }

      int nprop;
      inputfile.read((char *)&nprop, sizeof(int));
      _propagate_state.resize(nprop);
      for(int i = 0; i < _propagate_state.size(); i++){
        _propagate_state[i].read(inputfile);
      }
#ifndef USE_VDISK
      inputfile.close();
#endif
    }

    void write_rho(const BMatrix<T> &m, const Basis& basis, int n, int position, bool is_rho_left = false) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file,255,"rho_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"rho_%s_%i_r.dat",_name,n);

      if(is_rho_left)
        if(position == LEFT)
          snprintf(file,255,"rho2_%s_%i_l.dat",_name,n);
        else
          snprintf(file,255,"rho2_%s_%i_r.dat",_name,n);

      cout << "SAVING RHO " << file << endl;
#ifdef USE_VDISK
      std::ostringstream outputfile(disk.oopen(file),std::ios::out|std::ios::binary);
#else
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
#endif
      basis.write(outputfile);
      m.write(outputfile);
#ifndef USE_VDISK
      outputfile.close();
#else
      outputfile << flush;
      disk.close(outputfile.str());
#endif
    }

    void read_rho(BMatrix<T> &m, Basis &basis, int n, int position, bool is_rho_left = false)
    {
      char file[255];
      if(position == LEFT)
        snprintf(file,255,"rho_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"rho_%s_%i_r.dat",_name,n);

      if(is_rho_left)
        if(position == LEFT)
          snprintf(file,255,"rho2_%s_%i_l.dat",_name,n);
        else
          snprintf(file,255,"rho2_%s_%i_r.dat",_name,n);


      cout << "READING RHO " << file << endl;
#ifdef USE_VDISK
      std::istringstream inputfile(disk.iopen(file),std::ios::in|std::ios::binary);
#else
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
#endif
      basis.read(inputfile);
      m.read(inputfile);
#ifndef USE_VDISK
      inputfile.close();
#endif
    }

    void write_block(const Block<T> &b, int n, int position) const
    {
      char file[255];
      if(position == LEFT)
        snprintf(file,255,"block_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"block_%s_%i_r.dat",_name,n);

      cout << "SAVING BLOCK " << file << endl;
#ifdef USE_VDISK
      std::ostringstream outputfile(disk.oopen(file),std::ios::out|std::ios::binary);
#else
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
#endif
      b.write(outputfile);
#ifndef USE_VDISK
      outputfile.close();
#else
      outputfile << flush;
      disk.close(outputfile.str());
#endif
    }

    void read_block(Block<T> &b, int n, int position)
    {
      if(n == 1){
        if(position == LEFT) 
          b = h.get_site(0);
        else
          b = h.get_site(lattice().size()-1);
        return; 
      }

      char file[255];
      if(position == LEFT)
        snprintf(file,255,"block_%s_%i_l.dat",_name,n);
      else
        snprintf(file,255,"block_%s_%i_r.dat",_name,n);

      cout << "READING BLOCK " << file << endl;
#ifdef USE_VDISK
      std::istringstream inputfile(disk.iopen(file),std::ios::in|std::ios::binary);
#else
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
#endif
      b.read(inputfile);
#ifndef USE_VDISK
      inputfile.close();
#endif
    }


    void write_iter(int position, bool symmetric = false) const
    {
      write_status();
      write_block(newblock,iter+1,position);
      write_rho(rho,rho_basis,iter+1,position);
//      if(!is_hermitian) write_rho(rho_left,rho_basis,iter+1,position,true);
      write_gs(gs,iter,position);
      if(symmetric){
        write_block(newblock,iter+1,1-position);
        write_rho(rho,rho_basis,iter+1,1-position);
        write_gs(gs,iter,1-position);
      }
    }

    void read_iter(int position, bool symmetric = false) 
    {
      read_block(newblock,iter+1,position);
      read_rho(rho,rho_basis,iter+1,position);
//      if(!is_hermitian) read_rho(rho_left,rho_basis,iter+1,position,true);
      read_gs(gs,iter,position);
      if(symmetric){
        read_block(newblock,iter+1,1-position);
        read_rho(rho,rho_basis,iter+1,1-position);
        read_gs(gs,iter,1-position);
      }
    }

    // SIGNALS
    
    System<T>& signal_connect(typename Signal<System<T> >::callback handler, size_t signal_id, void *data)
      {
         Signal<System<T> > signal(handler, this, signal_id, data);
         handlers.push_back(signal);
         return *this;
      }

    bool signal_emit(size_t signal_id)
      {
         typename std::list<Signal<System<T> > >::reverse_iterator iter;
         for(iter = handlers.rbegin(); iter != handlers.rend(); iter++){
           Signal<System<T> > &signal = *iter;
           if(signal.signal_id() == signal_id) return signal.emit();
         }
         return true;
      }

};

/////////////////////////////////////////////////////////////////////////
// warmup:
// warmup_loop: 
// First sweep through the lattice to build the blocks.
// Symmetric growing (1D): we first sweep to ls/2-1 and build right and
// left blocks at the same time. We then continue from ls/2-1 to ls-3 for
// the remaining left blocks. The system is then ready for the main loop.
// Non symmetric growing (2D): we first weep from left to right, using small
// 1-site/2-sites environment blocks on the right (to have an even total 
// number of sites all the time), until we reach the final size. We then sweep
// from right to left to build the right blocks. Finally, a sweep from left to 
// right to leave the system ready for the main loop. 
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::warmup_loop(size_t t, size_t first_iter, bool randomize)
{
  Vector<QN> qns(lattice().ls());
  qns = qnt;
  _custom_qns = false;
  real_warmup_loop(t, qns, first_iter, randomize);
}

template <class T>
void
System<T>::warmup_loop(size_t t, const Vector<QN> &qns, size_t first_iter, bool randomize)
{
  _custom_qns = true;
  real_warmup_loop(t, qns, first_iter, randomize);
}

template <class T>
void
System<T>::real_warmup_loop(size_t t, const Vector<QN> &qns, size_t first_iter, bool randomize)
{
  m = t;
  int ls = lattice().ls();
  int l1 = std::max(size_t(1),h.site.lattice().size());
  _timer.Start();
  bool save_use_seed = _use_seed;
  bool save_use_hamis_seed = _use_hamis_seed;
  _use_seed = false;
  _use_hamis_seed = false;
  _warmup_growing = true;

  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;

  start();
//****************************************************************

// LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)
  cout << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";
  outputfile << "LEFT-TO-RIGHT sweep to get B(1)B(2)...B(L/2-1)\n";

  dir = LEFT2RIGHT;
  ls /= l1;
  int real_ls = ls;
  if(use_single_site()) ls++;
  int sweep_max = full_sweep() ? ls-2 : ls-3;
//  int sweep_max = ls-2;
  int mid_size = (ls/2)*2 == ls ? ls/2 : ls/2+1;
  int max = _grow_symmetric ? mid_size : sweep_max;

  for(iter = first_iter; iter < max; iter++)
  {
//    qn = qns[_grow_symmetric?iter*2+2-1:(SGN(iter) == 1?iter+3:iter+2)];
    qn = qns[_grow_symmetric?iter*2+2-1:iter+3];

    if(_grow_symmetric)
      read_block(rightblock, 2*iter+2 <= real_ls ? iter : iter-1, RIGHT);
    else
      rightblock = h.get_site(iter+2);
/*
#ifdef GRANDCANONICAL
    read_block(rightblock, _grow_symmetric ? iter : 1, RIGHT);
#else // GRANDCANONICAL
    read_block(rightblock, _grow_symmetric ? iter : (SGN(iter) == 1 && _grand_canonical  == 0 ? 2 : 1), RIGHT);
#endif // GRANDCANONICAL
*/
    read_block(leftblock, iter, LEFT);
    const Block<T>* _site1 = _grow_symmetric && _grow_outward ? &h.get_site(ls/2) : &h.get_site(iter);
    const Block<T>* _site2 = _grow_symmetric ? (_grow_outward ? &h.get_site(ls/2+1) : &h.get_site(ls-iter-1)) : &h.get_site(iter+1);

//    if(use_single_site()) _site2 = empty_block<T>::get(); 
    const Block<T> &site1 = *_site1;
    const Block<T> &site2 = *_site2;

    _timer.Lap();
    cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    init_iteration(leftblock, site1, site2, rightblock, false);
    if(randomize){
      gs.randomize();
      gs /= sqrt(product(gs,gs));
    } else
      diagonalize(false);

    m1 = m4 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    if(_grow_outward){
      _b2 = &(h.get_site(ls/2-iter-2)); 
      _b1 = &leftblock; 
    }
    rotate(LEFT, newblock);
    write_iter(LEFT);

#ifndef GRANDCANONICAL 
    if((_grow_symmetric && iter < mid_size-1) || (_grand_canonical  == 0 && iter == 1)) { 
      dir = RIGHT2LEFT;
      truncate(RIGHT, m4);
      if(_grow_outward){
        _b3 = &(h.get_site(ls/2+iter+2)); 
        _b4 = &leftblock; 
      }
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
      dir = LEFT2RIGHT;
    }
#endif //GRANDCANONICAL

    cout << "===========================================\n";
    cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
    cout << "===========================================\n";
    outputfile << "===========================================\n";
    outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
    outputfile << "===========================================\n";
  }

  if(_grow_symmetric){

//  LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)
    cout << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";
    outputfile << "LEFT-TO-RIGHT sweep to get B(L/2)B(L/2+1)...B(L-3)\n";

    qn = qnt;
    for(iter = std::max(int(first_iter),int(mid_size)); iter < sweep_max; iter++)
    {
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>* _site2 = &h.get_site(iter+1);
      if(use_single_site()) _site2 = empty_block<T>::get();
      const Block<T>& site2 = *_site2;
  
      _timer.Lap();
      cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
      outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
  
      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
  
      init_iteration(leftblock, site1, site2, rightblock, false);
      if(randomize) {
        gs.randomize();
        gs /= sqrt(product(gs,gs));
      } else
        diagonalize(false);
  
      m1 = std::min(m1*m2,m);
  
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
  
      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  } else {
    cout << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    outputfile << "RIGHT-TO-LEFT sweep to get B(1)B(2)...B(L/2-1)\n";
    dir = RIGHT2LEFT;

    for(iter = 1; iter < sweep_max; iter++)
    {
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>* _site1 = &h.get_site(real_ls-2-iter);
      const Block<T>& site2 = h.get_site(real_ls-2-iter+1);
      if(use_single_site()) _site1 = empty_block<T>::get();
      const Block<T>& site1 = *_site1;

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.Lap();
      cout << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "WARMUP RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, false); 
      if(randomize) {
        gs.randomize();
        gs /= sqrt(product(gs,gs));
      } else
        diagonalize(false);

      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);

      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }

    dir = LEFT2RIGHT;
    _use_seed = save_use_seed;
    _use_hamis_seed = save_use_hamis_seed;
    _warmup_growing = false;
    _custom_qns = false;

    for(iter = 1; iter < sweep_max; iter++)
    {
      read_block(rightblock, ls-iter-2, RIGHT);
      read_block(leftblock, iter, LEFT);
      const Block<T>& site1 = h.get_site(iter);
      const Block<T>* _site2 = &h.get_site(iter+1);
      if(use_single_site()) _site2 = empty_block<T>::get();
      const Block<T>& site2 = *_site2;
  
      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
 
      _timer.Lap();
      cout << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
      outputfile << "WARMUP LEFT-TO-RIGHT ITERATION " << iter << endl;
 
//      if(iter < sweep_max/2) {
//        init_iteration(leftblock, site1, site2, rightblock, false); 
//        diagonalize(false);
//      }else{
        init_iteration(leftblock, site1, site2, rightblock, save_use_seed); 
        if(randomize) {
          gs.randomize();
          gs /= sqrt(product(gs,gs));
        } else
          diagonalize(_use_seed);
//      }
  
      m1 = std::min(m1*m2,m);
  
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
  
      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
    }
  }
    
  qn = qnt;
  outputfile.close();
  _in_warmup = false;
  _use_seed = save_use_seed;
  _use_hamis_seed = save_use_hamis_seed;
}

/////////////////////////////////////////////////////////////////////////
// thermo_field:
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::thermo_field(size_t state, Matrix<size_t> *states, Vector<double> *weights, size_t ntrunc, Vector<T> *gs_coef)
{
  int ls = lattice().ls();
  _timer.Start();
  size_t save_qn_mask = _grand_canonical;
  _grand_canonical = 0;
//  _grand_canonical ^= MASK_QN_KX;
//  _grand_canonical ^= MASK_QN_KY;
//  _in_warmup = true;
  bool save_save_hmatrix = _save_hmatrix;
  _save_hmatrix = false;
  bool save_use_composite = _use_composite;
  _use_composite = false;

  start(false);
  init_correlations();

  Matrix<size_t> _states(1,ls);
  _states = state;
  if(states) _states = *states;

  Vector<double> _weights(_states.size1());
  _weights = 1;
  if(weights) _weights = *weights; 

  QN qn_iter;

//****************************************************************
// We build the thermo field
//****************************************************************
  _target.resize(_states.size1());
  _target_weight.resize(_states.size1());
  _target_weight = _weights;
  double tot = 0;
  for(int i = 0; i < _weights.size(); i++) tot += _weights[i];
  _target_weight /= tot;
  _ntargets = _states.size1();
  bool save_use_single_site = use_single_site();
  _use_single_site = false;
  
  for(iter = 1; iter <= ls-3; iter++)
  {
    dir = LEFT2RIGHT;
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);
    rightblock = h.get_site(iter+2);

    cout << "THERMO-FIELD BUILDING ITERATION (LEFT-TO-RIGHT) " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, false);
    cout << m1 << " " << m2 << " " << m3 << " " << m4 << endl;

    cout << _b1->basis().size() << " " << _b2->basis().size() << " " << _b3->basis().size() << " " << _b4->basis().size() << endl;
    qn_iter = QN();
    if(iter == 1){
      for(int i = 0; i < _states.size1(); i++){
        QN qn1 = _b1->basis()[_states(i,0)].qn();
        QN qn2 = _b2->basis()[_states(i,1)].qn();
        QN qn3 = _b3->basis()[_states(i,2)].qn();
        QN qn4 = _b4->basis()[_states(i,3)].qn();
        qn_iter = qn1 + qn2 + qn3 + qn4;
        cout << _states(i,0) << " " << qn1 << endl;
        cout << _states(i,1) << " " << qn2 << endl;
        cout << _states(i,2) << " " << qn3 << endl;
        cout << _states(i,3) << " " << qn4 << endl;
        _target[i].set_qn_mask(qn_iter, _grand_canonical);
        _target[i].resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
        cout << "GS QN " << _target[i].qn() << endl;
        cout << "GS SIZE " << _target[i].size() << endl;
        _target[i] = T(0);
        _target[i](_states(i,0),_states(i,1),_states(i,2),_states(i,3)) = T(1.);
        gs = _target[i];
//        measure();
      }
    }else{
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      read_rho(rho1, basis1, iter, LEFT);
      read_gs(gs, iter-1, LEFT);

      for(int i = 0; i < _states.size1(); i++){
        QN qn1 = _b1->basis()[0].qn();
        QN qn2 = _b2->basis()[_states(i,iter)].qn();
        QN qn3 = _b3->basis()[_states(i,iter+1)].qn();
        QN qn4 = _b4->basis()[_states(i,iter+2)].qn();
        seed = _target[i];
        seed.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
        seed = T(0);

        new_seed(_target[i], seed, rho1, basis1, LEFT);
        qn_iter = _target[i].qn() + qn4;
        cout << "GS QN " << qn_iter << endl;
        _target[i].set_qn_mask(qn_iter, _grand_canonical);
        _target[i].resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
        _target[i] = T(0);
  
        typename VectorState<T>::iterator siter;

/*
cout << "HOLA QN " << gs.qn().n() << " " << seed.qn().n() << endl;
      for(siter = seed.subspace_begin(); siter != seed.subspace_end(); siter++){
        StateSpace s = *siter;
        StateSpace sgs = gs.get_qn_space(s[1].qn(),s[2].qn(),s[3].qn(),qn4);
        state_slice<T> ss = seed(s); 
        state_slice<T> ssgs = gs(sgs);
cout << "HOLA " << ss.size() << " " << ssgs.size() << endl;
cout << "HOLA " << s[1].qn().n() << endl;
cout << "HOLA " << s[2].qn().n() << endl;
cout << "HOLA " << s[3].qn().n() << endl;
cout << "HOLA " << s[4].qn().n() << endl;
        for(int i1 = 0; i1 < s[1].dim(); i1++)
          for(int i2 = 0; i2 < s[2].dim(); i2++)
            for(int i3 = 0; i3 < s[3].dim(); i3++)
              for(int i4 = 0; i4 < s[4].dim(); i4++)
                cout << "HOLA " << i1 << " " << i2 << " "<< i3 << " " << i4 << " " << ss(i1,i2,i3,i4) << endl;
        
      }
*/

        for(siter = seed.subspace_begin(); siter != seed.subspace_end(); siter++){
          StateSpace s = *siter;
          StateSpace sgs = _target[i].get_qn_space(s[1].qn(),s[3].qn(),s[4].qn(),qn4);
          if(sgs.start() == -1) continue;
          state_slice<T> ss = seed(s); 
          state_slice<T> ssgs = _target[i](sgs);
          for(int i1 = 0; i1 < s[1].dim(); i1++){
            for(int i2 = 0; i2 < s[3].dim(); i2++){
              for(int i3 = 0; i3 < s[4].dim(); i3++){
                ssgs(i1,i2,i3,_states(i,iter+2)-sgs[4].begin()) = ss(i1,0,i2,i3);
              }
            }
          }
        }
        gs = _target[i];
//        measure();
      }
    }

    truncate(LEFT, std::max(ntrunc,site1.basis().size()*site1.basis().size()), false, false);
    rotate(LEFT, newblock);
    write_iter(LEFT);
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }

//////////////////////////////////////////////////////////////////
  Vector<T> _gs_coef(_target.size());
  _gs_coef = T(1);
  if(!gs_coef) 
    gs_coef = &_gs_coef;

  gs = _target[0]*gs_coef->operator()(0);
  for(int i = 1; i < _target.size(); i++) gs += _target[i]*gs_coef->operator()(i);
  cout << "COMBINED GROUND STATE " << _grand_canonical << " " << sqrt(product(gs,gs)) << endl;
  measure();

////////////////////////////////////////////////////////////
  typename VectorState<T>::iterator siter;

/*
  for(siter = _target[0].subspace_begin(); siter != _target[0].subspace_end(); siter++){
    StateSpace s = *siter;
    state_slice<T> ss = _target[0](s);
    QN qn0 = s[1].qn() + s[2].qn() + s[3].qn() + s[4].qn();

//    if(qn0 == qnt){
//      cout << "HOLA SUBSPACE " << s[1].qn() << " " << s[2].qn() << " " << s[3].qn() << " " << s[4].qn() << endl;
//      for(int i1 = 0; i1 < s[1].dim(); i1++)
//        for(int i2 = 0; i2 < s[2].dim(); i2++)
//          for(int i3 = 0; i3 < s[3].dim(); i3++)
//            for(int i4 = 0; i4 < s[4].dim(); i4++){
//              cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << ss(i1,i2,i3,i4) << endl;
//            }
//    }

      cout << "HOLA SUBSPACE " << s[1].qn() << " " << s[2].qn() << " " << s[3].qn() << " " << s[4].qn() << " " << qn0 << endl;
      for(int i1 = 0; i1 < s[1].dim(); i1++)
        for(int i2 = 0; i2 < s[2].dim(); i2++)
          for(int i3 = 0; i3 < s[3].dim(); i3++)
            for(int i4 = 0; i4 < s[4].dim(); i4++){
              if(fabs(ss(i1,i2,i3,i4)) > 1.e-4) cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << ss(i1,i2,i3,i4) << endl;
            }

  }
*/
  cout << "==========================================================\n";
  VectorState<T> auxgs = gs;  
  _grand_canonical = save_qn_mask;
  gs.set_qn_mask(qnt, _grand_canonical);
  gs.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
  gs = T(0);

  for(siter = gs.subspace_begin(); siter != gs.subspace_end(); siter++){
    StateSpace s = *siter;
    StateSpace auxs = auxgs.get_qn_space(s[1].qn(),s[2].qn(),s[3].qn(),s[4].qn());
    state_slice<T> ss = gs(s);
    state_slice<T> auxss = auxgs(auxs);
    for(int i1 = 0; i1 < s[1].dim(); i1++)
      for(int i2 = 0; i2 < s[2].dim(); i2++)
        for(int i3 = 0; i3 < s[3].dim(); i3++)
          for(int i4 = 0; i4 < s[4].dim(); i4++){
            ss(i1,i2,i3,i4) = auxss(i1,i2,i3,i4);
          }
  }
  gs /= sqrt(product(gs,gs));

  measure();
  _target[0] = product(*this, gs);
  cout << "THERMO FIELD NORM = " << sqrt(product(gs,gs)) << endl;
  cout << "THERMO FIELD VARIATIONAL ENERGY = "  << product(gs,_target[0])/product(gs,gs) << endl;
  _target.resize(1);
  _target_weight.resize(1);
  _ntargets = 1;
  _target[0] = gs;
  _target_weight[0] = 1.;
  write_gs(gs,ls-3,LEFT);
  _in_warmup = false;
//////////////////////////////////////////////////////////////////

  for(iter = 1; iter <= ls-3; iter++)
  {
    dir = RIGHT2LEFT;
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(leftblock, ls-iter-2, LEFT);
    read_block(rightblock, iter, RIGHT);
    const Block<T>& site1 = h.get_site(ls-2-iter);
    const Block<T>& site2 = h.get_site(ls-2-iter+1);

    cout << "THERMO-FIELD RIGHT-TO-LEFT ITERATION " << iter << endl;
    init_iteration(leftblock, site1, site2, rightblock, false);
    cout << m1 << " " << m2 << " " << m3 << " " << m4 << endl;

    cout << _b1->basis().size() << " " << _b2->basis().size() << " " << _b3->basis().size() << " " << _b4->basis().size() << endl;

/*

cout << "HOLA QN " << gs.qn().n() << " " << seed.qn().n() << endl;
      typename VectorState<T>::iterator siter;
      for(siter = gs.subspace_begin(); siter != gs.subspace_end(); siter++){
        StateSpace s = *siter;
        state_slice<T> ssgs = gs(s);
cout << "HOLA " << s[1].qn().n() << endl;
cout << "HOLA " << s[2].qn().n() << endl;
cout << "HOLA " << s[3].qn().n() << endl;
cout << "HOLA " << s[4].qn().n() << endl;
        for(int i1 = 0; i1 < s[1].dim(); i1++)
          for(int i2 = 0; i2 < s[2].dim(); i2++)
            for(int i3 = 0; i3 < s[3].dim(); i3++)
              for(int i4 = 0; i4 < s[4].dim(); i4++)
                cout << "HOLA " << i1 << " " << i2 << " "<< i3 << " " << i4 << " " << ssgs(i1,i2,i3,i4) << endl;
        
      }
*/

    if(iter == 1){
      read_gs(gs, ls-3, LEFT);
      cout << "GS SIZE " << gs.size() << endl;
      cout << "QN GS " << gs.qn().n() << " " << gs.qn().sz() << endl;
      cout << "QN " << qnt.n() << " " << qnt.sz() << endl;
    } else {
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      read_rho(rho1, basis1, iter, RIGHT);
      read_rho(rho2, basis2, lattice().ls()-iter-2+1, LEFT);
      read_gs(gs, iter-1, RIGHT);
      cout << "GS SIZE " << gs.size() << endl;
      cout << "QN GS " << gs.qn().n() << " " << gs.qn().sz() << endl;
      cout << "QN " << qnt.n() << " " << qnt.sz() << endl;
      cout << "NEW SEED" << endl;
      seed.set_qn_mask(gs.qn(), _grand_canonical);
      seed.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
      new_seed(gs, seed, rho1, rho2, basis1, basis2, RIGHT);
      gs = seed;
    }
    measure();

    truncate(RIGHT, std::max(ntrunc,site2.basis().size()*site2.basis().size()), true, false);
    rotate(RIGHT, newblock);
    write_iter(RIGHT);
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }
  for(iter = 1; iter <= ls-3; iter++)
  {
    dir = LEFT2RIGHT;
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(leftblock, iter, LEFT);
    read_block(rightblock, ls-iter-2, RIGHT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>& site2 = h.get_site(iter+1);

    cout << "THERMO-FIELD LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, false);
    cout << m1 << " " << m2 << " " << m3 << " " << m4 << endl;

    cout << _b1->basis().size() << " " << _b2->basis().size() << " " << _b3->basis().size() << " " << _b4->basis().size() << endl;
    cout <<"GS SIZE " << gs.size() << endl;

    if(iter != 1){
      BMatrix<T> rho1;
      Basis basis1;
      BMatrix<T> rho2;
      Basis basis2;

      read_rho(rho1, basis1, iter, LEFT);
      read_rho(rho2, basis2, lattice().ls()-iter-2+1, RIGHT);
      read_gs(gs, iter-1, LEFT);
      gs.resize(_grand_canonical);
      cout << "NEW SEED" << endl;
      seed.set_qn_mask(gs.qn(), _grand_canonical);
      seed.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis());
      new_seed(gs, seed, rho1, rho2, basis1, basis2, LEFT);
      gs = seed;
    }
    measure();

    if(iter < ls-3 || save_use_single_site){
      truncate(LEFT,std::max(ntrunc,site1.basis().size()*site1.basis().size()), true, false);
      rotate(LEFT, newblock);
    }
    write_iter(LEFT);
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }
  _target[0] = product(*this, gs);
  cout << "THERMO FIELD VARIATIONAL ENERGY = "  << product(gs,_target[0])/product(gs,gs) << endl;
  _use_single_site = save_use_single_site;
  _save_hmatrix = save_save_hmatrix;
  _use_composite = save_use_composite;
}

////////////////////////////////////////////////////////////////////////
// set_iteration:
////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::set_iteration(size_t _dir, size_t _iter, bool _rotate, int new_size, bool use_seed)
{
  int ls = lattice().ls();
  int l1 = std::max(size_t(1),h.site.lattice().size());
  ls /= l1;
  int real_ls = ls;
  qn = qnt;

  dir = _dir;
  iter = _iter;

  signal_emit(SYSTEM_SIGNAL_START_ITER);
  if(dir == RIGHT2LEFT){
    cout << "RIGHT-TO-LEFT ITERATION " << iter << endl;
    read_block(leftblock, ls-iter-2, LEFT);
    read_block(rightblock, iter, RIGHT);
    const Block<T>* _site1 = &h.get_site(real_ls-2-iter);
    const Block<T>& site2 = h.get_site(real_ls-2-iter+1);
    if(use_single_site()) _site1 = empty_block<T>::get();
    const Block<T>& site1 = *_site1;

    _b1 = &leftblock;
    _b2 = &site1;
    _b3 = &site2;
    _b4 = &rightblock;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    if(use_seed){
      init_iteration(leftblock, site1, site2, rightblock, _use_seed, false);
      gs = seed;
    } else {
      cout << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
      cout << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
      cout << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
      cout << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 

      create_interactions(false);
      read_gs(gs, iter, RIGHT);
    }
    if(_rotate){
//      diagonalize(_use_seed); 
      truncate(RIGHT, std::min(m4*m3,new_size));
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
    }
    if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
  } else {
    cout << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>* _site2 = &h.get_site(iter+1);
    if(use_single_site()) _site2 = empty_block<T>::get();
    const Block<T>& site2 = *_site2;

    _b1 = &leftblock;
    _b2 = &site1;
    _b3 = &site2;
    _b4 = &rightblock;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    if(use_seed){
      init_iteration(leftblock, site1, site2, rightblock, _use_seed, false);
      gs = seed;
    } else {
      cout << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
      cout << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
      cout << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
      cout << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 

      create_interactions(false);
      read_gs(gs, iter, LEFT);
    }
    if(_rotate){
//      diagonalize(_use_seed); 
      truncate(LEFT, std::min(m1*m2,new_size));
      rotate(LEFT, newblock);
      write_iter(LEFT);
    }
    if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
  }
}


/////////////////////////////////////////////////////////////////////////
// sweep: 
// left-to-right and right-to-left sweeps, with t1(t2) number
// of states, respectively. We repeat this procedure several times
// in the main program with growing t1(t2) until we reach convergence
// with the final desired number of states.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::sweep(size_t t1, size_t t2, size_t _dir, int start)
{
  int ls = lattice().ls();
  int l1 = std::max(size_t(1),h.site.lattice().size());
  _in_warmup = false;
  _custom_qns = false;
  ls /= l1;
  int real_ls = ls;
  if(use_single_site()) ls++;
  int sweep_max = full_sweep() ? ls-2 : ls-3;
  qn = qnt;

  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;

  cout << "===========================================\n";
  cout << "NEW SWEEP\n";
  cout << "===========================================\n";
  outputfile << "===========================================\n";
  outputfile << "NEW SWEEP\n";
  outputfile << "===========================================\n";

// RIGHT-TO-LEFT sweep
  
  dir = _dir;

  if(dir == RIGHT2LEFT){
    for(iter = start; iter < sweep_max; iter++)
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(leftblock, ls-iter-2, LEFT);
      read_block(rightblock, iter, RIGHT);
      const Block<T>* _site1 = &h.get_site(real_ls-2-iter);
      const Block<T>& site2 = h.get_site(real_ls-2-iter+1);
      if(use_single_site()) _site1 = empty_block<T>::get();
      const Block<T>& site1 = *_site1;

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();

      _timer.Lap();
      cout << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;

      init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
      diagonalize(_use_seed); 

      m = t1;
      m4 = std::min(m4*m3,m);

      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }

    start = 1;
  }

// LEFT-TO-RIGHT sweep

  dir = LEFT2RIGHT;

  for(iter = start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>* _site2 = &h.get_site(iter+1);
    if(use_single_site()) _site2 = empty_block<T>::get();
    const Block<T>& site2 = *_site2;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.Lap();
    cout << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 
    diagonalize(_use_seed);

    m = t2;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    rotate(LEFT, newblock);
    write_iter(LEFT);
    if(!_store_products){
      if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
    }

    cout << "===========================================\n";
    cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
    cout << "===========================================\n";
    outputfile << "===========================================\n";
    outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
    outputfile << "===========================================\n";
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }

  signal_emit(SYSTEM_SIGNAL_END_SWEEP);
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::empty_sweep(size_t t, size_t _dir, int _start, bool half_sweep)
{
  int ls = lattice().ls();
  int l1 = std::max(size_t(1),h.site.lattice().size());
  dir = _dir;
  ls /= l1;
  int real_ls = ls;
  if(use_single_site()) ls++;
  int sweep_max = full_sweep() ? ls-2 : ls-3;

  if(dir == RIGHT2LEFT){
    for(iter = _start; iter < sweep_max; iter++) 
    {
      cout << "EMPTY SWEEP ITERATION\n";
      set_iteration(dir, iter, true, t, true);
    }
    _start = 1;
    if(half_sweep) return;
  } 
  dir = LEFT2RIGHT;
  sweep_max = full_sweep() ? ls-2 : ls/2-1;
  if(dir == LEFT2RIGHT){
    for(iter = _start; iter < sweep_max; iter++) 
    {
      cout << "EMPTY SWEEP ITERATION\n";
      set_iteration(dir, iter, true, t, true);
    }
  } 
  
  signal_emit(SYSTEM_SIGNAL_END_SWEEP);

}
/////////////////////////////////////////////////////////////////////////
// init_correlations:
// We store in site blocks the operators used for measuring correlations
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::init_correlations()
{

  corr.reorder_terms(true);
  return;
}

/////////////////////////////////////////////////////////////////////////
// final_sweep:
// Like a regular sweep but from 1 to ls/2-1 to build
// right and left symmetric blocks (we don't care about the others for the
// measurement). The main difference is that at each step we add the 
// measurement operators to the block (stored in the list "ops")  
/////////////////////////////////////////////////////////////////////////

template<class T>
void
System<T>::final_sweep(size_t t, size_t _dir, int _start, bool _rotate)
{
// Last iteration to get a symmetric block B(L/2-1)..B(L/2-1)

  int ls = lattice().ls();
  int l1 = std::max(size_t(1),h.site.lattice().size());
  _in_warmup = false;

  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;
  outputfile << "Last iteration to get a symmetric block B(L/2-1)..B(L/2-1)\n";

//****************************************************************
// We store in site blocks the operators used for measuring correlations
//****************************************************************
  init_correlations();
//****************************************************************

  dir = _dir;
  ls /= l1;
  int real_ls = ls;
  if(use_single_site()) ls++;
  int sweep_max = full_sweep() ? ls-2 : ls-3;

  if(dir == RIGHT2LEFT){
    for(iter = _start; iter < sweep_max; iter++) // ls/2-1; iter++);
    {
      signal_emit(SYSTEM_SIGNAL_START_ITER);
      read_block(rightblock, iter, RIGHT);
      read_block(leftblock, ls-iter-2, LEFT);
      const Block<T>* _site1 = &h.get_site(real_ls-2-iter);
      const Block<T>& site2 = h.get_site(real_ls-2-iter+1);
      if(use_single_site()) _site1 = empty_block<T>::get();
      const Block<T>& site1 = *_site1;

      m1 = leftblock.dim();
      m4 = rightblock.dim();
      m2 = site1.dim();
      m3 = site2.dim();
   
      _timer.Lap();
      cout << "FINAL SWEEP ITERATION " << endl;
      outputfile << "FINAL SWEEP ITERATION " << endl;
      cout << "RIGHT-TO-LEFT ITERATION " << iter << endl;
      outputfile << "RIGHT-TO-LEFT ITERATION " << iter << endl;
  
      init_iteration(leftblock, site1, site2, rightblock, _use_seed);
  
      diagonalize(_use_seed);
  
      m = t;
      m4 = std::min(m4*m3,m);
  
      truncate(RIGHT, m4);
      rotate(RIGHT, newblock);
      write_iter(RIGHT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }

      cout << "===========================================\n";
      cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
      cout << "===========================================\n";
      outputfile << "===========================================\n";
      outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
      outputfile << "===========================================\n";
      signal_emit(SYSTEM_SIGNAL_END_ITER);
    }
    _start = 1;
  }

  dir = LEFT2RIGHT;
  sweep_max = full_sweep() ? ls-2 : ls/2-1;

  for(iter = _start; iter < sweep_max; iter++)
  {
    signal_emit(SYSTEM_SIGNAL_START_ITER);
    read_block(rightblock, ls-iter-2, RIGHT);
    read_block(leftblock, iter, LEFT);
    const Block<T>& site1 = h.get_site(iter);
    const Block<T>*_site2 = &h.get_site(iter+1);
    if(use_single_site()) _site2 = empty_block<T>::get();
    const Block<T>& site2 = *_site2;

    m1 = leftblock.dim();
    m4 = rightblock.dim();
    m2 = site1.dim();
    m3 = site2.dim();

    _timer.Lap();
    cout << "FINAL SWEEP ITERATION " << endl;
    outputfile << "FINAL SWEEP ITERATION " << endl;
    cout << "LEFT-TO-RIGHT ITERATION " << iter << endl;
    outputfile << "LEFT-TO-RIGHT ITERATION " << iter << endl;

    init_iteration(leftblock, site1, site2, rightblock, _use_seed); 

    diagonalize(_use_seed);

    m = t;
    m1 = std::min(m1*m2,m);

    truncate(LEFT, m1);
    rotate(LEFT, newblock);
    write_iter(LEFT);
    if(!_store_products){
      if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
    }

    cout << "===========================================\n";
    cout << "Iteration time: " << _timer.LapTime().c_str() << endl;
    cout << "===========================================\n";
    outputfile << "===========================================\n";
    outputfile << "Iteration time: " << _timer.LapTime().c_str() << endl;
    outputfile << "===========================================\n";
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }

// Calculate the ground state in the symmetric system

  if(!full_sweep()){
    cout << "FINAL SWEEP - LAST ITERATION " << endl;
    outputfile << "FINAL SWEEP - LAST ITERATION " << endl;
    read_block(rightblock, IS_EVEN(ls) ? ls/2-1 : ls/2, RIGHT);
    read_block(leftblock, ls/2-1, LEFT);
    const Block<T>& site1 = h.get_site(ls/2-1);
    const Block<T>* _site2 = &h.get_site(ls/2);
    if(use_single_site()) _site2 = empty_block<T>::get();
    const Block<T>& site2 = *_site2;
    init_iteration(leftblock, site1, site2, rightblock, _use_seed);  
//  set_lanczos_tolerance(std::min(1.e-16,_lanczos_tol));
//    set_lanczos_maxiter(30);
    diagonalize(_use_seed);
    if(_rotate){
      m = t;
      m1 = std::min(m1*m2,m);
                                                                               
      truncate(LEFT, m1);
      rotate(LEFT, newblock);
      write_iter(LEFT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }
    } else {
      write_gs(gs, ls/2-1, LEFT);
      if(!_store_products){
        if (signal_emit(SYSTEM_SIGNAL_MEASURE)) measure();
      }
    }
    signal_emit(SYSTEM_SIGNAL_END_ITER);
  }
  
  signal_emit(SYSTEM_SIGNAL_END_SWEEP);
  outputfile.close();
}

/////////////////////////////////////////////////////////////////////////
// main_loop: 
// if we don't care about increasing the number of states 
// in the block we can just run main_loop with a given number of states.
/////////////////////////////////////////////////////////////////////////

template<class T>
void
System<T>::main_loop(size_t nsweeps, size_t t)
{
  m = t;
  for(int _sweep = 1; _sweep <= nsweeps; _sweep++) { sweep(m,m); }
}

/////////////////////////////////////////////////////////////////////////
// init_iteration: 
// initialize the blocks, quantum numbers, and transform the ground
// state of the previous iteration as a new seed for the new one.
/////////////////////////////////////////////////////////////////////////
template <class T>
void
System<T>::create_interactions(bool create_composite)
{
#ifndef USE_PRODUCT_DEFAULT
  cout << "Creating interaction terms" << endl;
  hint(MASK_BLOCK1|MASK_BLOCK2);
  hint(MASK_BLOCK3|MASK_BLOCK4);
  hint(MASK_BLOCK2|MASK_BLOCK3, false);
  if(create_composite)
    build_composite_operators();
#endif
}

template <class T>
void
System<T>::init_iteration(const B&b1, const B&b2, const B&b3, const B&b4, bool use_seed, bool create_composite)
{
  CTimer clock;
  clock.Start(); 

  _b1 = &b1;
  _b2 = &b2;
  _b3 = &b3;
  _b4 = &b4;

//  if(!use_seed) return;

  int ltotal = 0, lnow = 0;
  for(int i = 0; i < h.lattice().ls(); i++) 
    ltotal += h.get_site(i).n_orbitals();
  for(int i = 0; i < size(); i++) 
    lnow += h.get_site(i).n_orbitals();

/*
  int n, spin; 
  if(qnt.n() <= ltotal){ 
    int nholes = ltotal - qnt.n();
    n = nholes >= lnow ? lnow : lnow - nholes;
    int nd = (int)(lnow * ((double)qnt.n() / (double)ltotal));
//    if(!_grow_symmetric && size() < h.lattice().ls())
//      n = lnow <  qnt.n() ? lnow : qnt.n();
      n = nd;
    if(n == 0) n = lnow;
    n = std::min(n,lnow);
//    if(SGN(n) == -1) n += 1; 
    spin = qnt.sz();
    if(!IS_EVEN(n) && IS_EVEN(spin)) spin += 1;
    if(IS_EVEN(n) && !IS_EVEN(spin)) spin -= 1;
    if(spin > n) spin = n;
  } else {
    n = lnow < ltotal ? lnow : qnt.n(); 
    spin = qnt.sz();
//    if(SGN(n) == -1) n += 1; 
    if(!IS_EVEN(n) && IS_EVEN(spin)) spin += 1;
    if(IS_EVEN(n) && !IS_EVEN(spin)) spin -= 1;
    if(spin > n) spin = n;
  }
  if(qnt.n() == 0) { spin = n = 0; }

  if(_in_warmup && !_custom_qns) 
    qn = QN(n,spin);
  if(!_in_warmup)
    qn = qnt;
*/

  QN aux_qn(999);
  PackedBasis::const_iterator biter1, biter2, biter3, biter4;
  for(biter1 = b1.basis().subspace_begin(); biter1 != b1.basis().subspace_end(); biter1++){
    for(biter2 = b2.basis().subspace_begin(); biter2 != b2.basis().subspace_end(); biter2++){
      for(biter3 = b3.basis().subspace_begin(); biter3 != b3.basis().subspace_end(); biter3++){
        for(biter4 = b4.basis().subspace_begin(); biter4 != b4.basis().subspace_end(); biter4++){
          QN sqn = (*biter1).qn()+(*biter2).qn()+(*biter3).qn()+(*biter4).qn();
//          cout <<  (*biter1).qn() << " " << (*biter2).qn() << "  "<< (*biter3).qn() << "  "<< (*biter4).qn() << endl;
          double this_dist = 0;
          double dist = 0;
#ifdef GRANDCANONICAL
          aux_qn.n() = mod2(aux_qn.n());
          sqn.n() = mod2(sqn.n());
#endif // GRANDCANONICAL
          for(int j = 0; j < QN::QN_LAST; j++){
            if(QN::qn_norm(j) > 1) {
              qnt[j] = qnt[j]%QN::qn_norm(j);
              sqn[j] = sqn[j]%QN::qn_norm(j);
            }
            if((_grand_canonical & (1 << j)) == 0) continue;
            int nd = (int)(lnow * ((double)qnt[j] / (double)ltotal));
            if(QN::qn_norm(j) > 1) {
              nd = nd%QN::qn_norm(j);
            }
#ifdef GRANDCANONICAL
            nd = mod2(nd);
#endif // GRANDCANONICAL
            this_dist += (sqn[j]-nd)*(sqn[j]-nd);
            dist += (aux_qn[j]-nd)*(aux_qn[j]-nd);
          }
          if(this_dist < dist) aux_qn = sqn;
        }
      }
    }
  }

  if(_in_warmup && !_custom_qns)
    qn = aux_qn;
  if(!_in_warmup || lnow == ltotal)
    qn = qnt;

#ifdef GRANDCANONICAL
    qn.n() = mod2(qn.n());
    qn.sz() = mod2(qn.sz());
#endif //GRANDCANONICAL

//  if(verbose() > 0) 
    {
       cout << "TOTAL SIZE : " << size() << " " << lnow << endl;
       cout << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
       cout << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
       cout << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
       cout << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 
       if(_grand_canonical == 0)
         cout << "Ntotal = GRAND CANONICAL\n";
       else
         cout << "Ntotal = " << qn.n() << " " << qnt.n() << endl;

       cout << "SZtotal = " << qn.sz() << " " << qnt.sz() << endl;
#ifdef USE_K
       cout << "KXtotal = " << qn.kx1() << "/" << qn.kx2() << endl;
       cout << "KYtotal = " << qn.ky1() << "/" << qn.ky2() << endl;
#endif // USE_K
       cout << "QN = " << qn << " " << qnt << endl;

    }
  char file[255];
  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;
  outputfile << "TOTAL SIZE : " << size() << " " << lnow << endl;
  outputfile << "LEFT BLOCK SIZE : " << _b1->lattice().size() << " " << _b1->n_orbitals() << endl;
  outputfile << "RIGHT BLOCK SIZE : " << _b4->lattice().size() << " " << _b4->n_orbitals() << endl;
  outputfile << "|" << _b1->n_orbitals() << "|-|" << _b2->n_orbitals() << "|-|" << _b3->n_orbitals() << "|-|" << _b4->n_orbitals() << "|" << endl; 
  outputfile << "|" << _b1->lattice().size() << "|-|" << _b2->lattice().size() << "|-|" << _b3->lattice().size() << "|-|" << _b4->lattice().size() << "|" << endl; 

  if(_grand_canonical == 0)
    outputfile << "Ntotal = GRAND CANONICAL\n";
  else
    outputfile << "Ntotal = " << qn.n() << " " << qnt.n() << endl;

  outputfile << "SZtotal = " << qn.sz() << " " << qnt.sz() << endl;
#ifdef USE_K
  outputfile << "KXtotal = " << qn.kx1() << "/" << qn.kx2() << endl;
  outputfile << "KYtotal = " << qn.ky1() << "/" << qn.ky2() << endl;
#endif // USE_K

  if(use_seed){
     BMatrix<T> rho1;
     Basis basis1;
     BMatrix<T> rho2;
     Basis basis2;

     if(_use_basic_seed){
       seed.set_qn_mask(qn, _grand_canonical);   
       seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
       seed = T(1);

     } else if(iter == 1){

       if(lattice().ls() != 4) {
         seed.set_qn_mask(qn, _grand_canonical);   
         seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);

         if(dir == RIGHT2LEFT){
           if(!use_single_site())
             read_rho(rho1, basis1, lattice().ls()-3, LEFT);
           else
             read_rho(rho1, basis1, lattice().ls()-3+1, LEFT);
           read_rho(rho2, basis2, 2, RIGHT);
           if(!use_single_site())           
             read_gs(gs, lattice().ls()-4,LEFT);
           else
             read_gs(gs, lattice().ls()-4+1,LEFT);
           gs.resize(_grand_canonical);
           cout << "NEW SEED" << endl; 
/*
           if(_target.size() > 1) {
             gs *= T(_target_weight[0]);
             for(int i = 1; i < _target.size(); i++){
               gs += T(_target_weight[i])*_target[i];
             }
           }
*/
           clock.Lap();
           new_seed(gs, seed, rho1, rho2, basis1, basis2, LEFT, use_single_site());

           VectorState<T> aux;
           for(int i = 0; i < _propagate_state.size(); i++){
             aux = _propagate_state[i]; 
             aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
             new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, LEFT, use_single_site());
             _propagate_state[i] = aux;
           }
           cout << "Lap: " << clock.LapTime().c_str() << endl;
         } else {
           if(!use_single_site())
             read_rho(rho1, basis1, lattice().ls()-3, RIGHT);
           else
             read_rho(rho1, basis1, lattice().ls()-3+1, RIGHT);
           read_rho(rho2, basis2, 2, LEFT);
           if(!use_single_site())           
             read_gs(gs, lattice().ls()-4,RIGHT);
           else
             read_gs(gs, lattice().ls()-4+1,RIGHT);
           gs.resize(_grand_canonical);
           cout << "NEW SEED" << endl; 
/*
           if(_target.size() > 1){
             gs *= T(_target_weight[0]);
             for(int i = 1; i < _target.size(); i++){
               gs += T(_target_weight[i])*_target[i];
             }
           }
*/
           clock.Lap();
           new_seed(gs, seed, rho1, rho2, basis1, basis2, RIGHT, use_single_site());

           VectorState<T> aux;
           for(int i = 0; i < _propagate_state.size(); i++){
             aux = _propagate_state[i]; 
             aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
             new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, RIGHT, use_single_site());
             _propagate_state[i] = aux;
           }
           cout << "Lap: " << clock.LapTime().c_str() << endl;
         }
       } else {
         seed = gs;
       }
     } else {

       seed.set_qn_mask(qn, _grand_canonical);   
       seed.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);

       if(dir == LEFT2RIGHT){
         read_rho(rho1, basis1, iter, LEFT);
         if(!use_single_site())
           read_rho(rho2, basis2, lattice().ls()-iter-2+1, RIGHT);
         else
           read_rho(rho2, basis2, lattice().ls()-iter-2+2, RIGHT);

         read_gs(gs, iter-1, LEFT);
         gs.resize(_grand_canonical);
         cout << "NEW SEED" << endl; 
/*
         if(_target.size() > 1){
           gs *= T(_target_weight[0]);
           for(int i = 1; i < _target.size(); i++){
             gs += T(_target_weight[i])*_target[i];
           }
         }
*/
         clock.Lap();
         new_seed(gs, seed, rho1, rho2, basis1, basis2, LEFT, use_single_site());

         VectorState<T> aux;
         for(int i = 0; i < _propagate_state.size(); i++){
           aux = _propagate_state[i]; 
           aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
           new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, LEFT, use_single_site());
           _propagate_state[i] = aux;
         }
         cout << "Lap: " << clock.LapTime().c_str() << endl;
       } else {
         read_rho(rho1, basis1, iter, RIGHT);
         if(!use_single_site())
           read_rho(rho2, basis2, lattice().ls()-iter-2+1, LEFT);
         else
           read_rho(rho2, basis2, lattice().ls()-iter-2+2, LEFT);
         read_gs(gs, iter-1, RIGHT);
         gs.resize(_grand_canonical);
/*
         cout << "NEW SEED" << endl; 
         if(_target.size() > 1){
           gs *= T(_target_weight[0]);
           for(int i = 1; i < _target.size(); i++){
             gs += T(_target_weight[i])*_target[i];
           }
         }
*/
         clock.Lap();
         new_seed(gs, seed, rho1, rho2, basis1, basis2, RIGHT, use_single_site());

         VectorState<T> aux;
         for(int i = 0; i < _propagate_state.size(); i++){
           aux = _propagate_state[i]; 
           aux.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
           new_seed(_propagate_state[i], aux, rho1, rho2, basis1, basis2, RIGHT, use_single_site());
           _propagate_state[i] = aux;
         }
         cout << "Lap: " << clock.LapTime().c_str() << endl;
       }

     }
 
  }

  gs.set_qn_mask(qn, _grand_canonical);
  gs.resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
#ifndef USE_PRODUCT_DEFAULT
  if(create_composite && _rotate_terms){
    create_interactions(create_composite);
    init_terms_composite(*this, this->gs); 
  }

  if(save_hmatrix()) {
    clock.Lap();
    hmatrix_create(*this);
    cout << "Save hmatrix time: " << clock.LapTime().c_str() << endl;
  }
#endif
//  hmatrix_verif(*this);
  cout << "Initialization time: " << clock.TotalTime().c_str() << endl;

}

#include "lanczos.cc"

/////////////////////////////////////////////////////////////////////////
// diagonalize:
// run lanczos to calculate the ground state
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::diagonalize(bool use_seed)
{   
  double tol = _lanczos_tol;
  Vector<double> a(1000);
  Vector<double> b(1000);
  Vector<double> g(1000);
  Vector<double> e(_nstates);

  CTimer clock;
  clock.Start(); 

  _state.resize(_nstates);
  for(int i = 0; i < _nstates; i++){
     _state[i].set_qn_mask(qn, _grand_canonical);
     _state[i].resize(_b1->_basis,_b2->_basis,_b3->_basis,_b4->_basis);
  }

//  if(verbose() > 0) 
    cout << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.qn() << " " << gs.size() << endl;

  int n = _lanczos_maxiter;
  if(_target_subspaces && _in_warmup && lattice().size() != size()) n = 5;
//  if(_target_subspaces && _in_warmup && _warmup_growing) n = 0;

  char file[255];
  snprintf(file,255,"vectors_%s.dat",_name);

  cout << ">>>>>>> LANCZOS <<<<<<<" << endl;

  if(real(product(seed,seed)) < 0.1) use_seed = false;

  if(gs.size() > 1) {
//    verif_hamiltonian<T, System<T>, VectorState<T> >(*this, _state[0]);
    if(use_lanczos()){
      double ei;
      VectorState<T> gs_left(_state[0]);
      if(this->is_hermitian())
        lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, _nstates, a, b, n, tol, use_seed, true, file, true);
      else 
        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], gs_left, seed, e[0], ei, n, tol, use_seed, true);
//        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], seed, e[0], n, tol, use_seed, true);

      gs = _state[0];
      gs /= sqrt(product(gs,gs));
    } else {
      gs = seed;
//      gs += product(*this, seed);
//      for(int i = 0; i < 5; i++) { gs += product(*this, gs); }
//      gs /= sqrt(product(gs,gs));
      _state[0] = product(*this,gs);
      e[0] = real(product(_state[0],gs)); 
    }
//    if(n == 0) gs.randomize(); 

/*
    if(n == 0){
      for(int i = 0; i < gs.size()/2; i++){
        gs[i] = (i-gs.size()/2.)*(i-gs.size()/2.)/gs.size()/gs.size();
      }
    }
*/
  } else {
    gs = T(1);
    _state[0] = product(*this,gs);
    e[0] = real(product(_state[0],gs));
    cout << "E0 = " << e[0] << endl;
    cout << "ENER = " << e[0] << endl;
  }

  energy.resize(_calc_gap+1);
  energy[0] = e[0];

  if(use_seed){

    typename VectorState<T>::const_iterator siter;
    Vector<SubSpace> ss(5);

    if(verbose() == 3)
      for(siter = seed.subspace_begin(); siter != seed.subspace_end(); siter++){
        ss[1] = (*siter)[1];
        ss[2] = (*siter)[2];
        ss[3] = (*siter)[3];
        ss[4] = (*siter)[4];

        cout << "------------------------------------------\n";
        for(int i1 = ss[1].begin(); i1 <= ss[1].end(); i1++)
        for(int i2 = ss[2].begin(); i2 <= ss[2].end(); i2++)
        for(int i3 = ss[3].begin(); i3 <= ss[3].end(); i3++)
        for(int i4 = ss[4].begin(); i4 <= ss[4].end(); i4++)
          cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << ss[1].qn().n() << " " << ss[2].qn().n() << " " << ss[3].qn().n() << " " << ss[4].qn().n() << " " << ss[1].dim() << " " << ss[2].dim() << " " << ss[3].dim() << " " << ss[4].dim() << " " << seed(i1,i2,i3,i4) << " " << gs(i1,i2,i3,i4) << endl;
      }
    
    cout << "OVERLAP " << product(seed,gs) << endl;
  } 
  else
  {
    typename VectorState<T>::const_iterator siter;
    Vector<SubSpace> ss(5);

    if(verbose() == 3)
      for(siter = gs.subspace_begin(); siter != gs.subspace_end(); siter++){
        ss[1] = (*siter)[1];
        ss[2] = (*siter)[2];
        ss[3] = (*siter)[3];
        ss[4] = (*siter)[4];

        cout << "------------------------------------------\n";
        for(int i1 = ss[1].begin(); i1 <= ss[1].end(); i1++)
        for(int i2 = ss[2].begin(); i2 <= ss[2].end(); i2++)
        for(int i3 = ss[3].begin(); i3 <= ss[3].end(); i3++)
        for(int i4 = ss[4].begin(); i4 <= ss[4].end(); i4++)
          cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << ss[1].qn().n() << " " << ss[2].qn().n() << " " << ss[3].qn().n() << " " << ss[4].qn().n() << " " << ss[1].dim() << " " << ss[2].dim() << " " << ss[3].dim() << " " << ss[4].dim() << " " << gs(i1,i2,i3,i4) << endl;
      }
  }

  snprintf(file,255,"iter_%s.dat",_name);
  ofstream outputfile(file,std::ios::out|std::ios::app);
  if(!outputfile) cout << "*** ERROR: could not open " << file << endl;

  outputfile << _b1->lattice().size() << " " << _b2->lattice().size() << " " << _b3->lattice().size() << " " <<  _b4->lattice().size() << endl;
  outputfile << _b1->dim() << " " << _b2->dim() << " " << _b3->dim() << " " <<  _b4->dim() << " " << gs.size() << endl;

  outputfile.precision(10);
  outputfile << "LANCZOS ENER = " << m << " " << gs.size() << " " << e[0] << endl;
  outputfile << "Lanczos time: " << clock.TotalTime().c_str() << endl;
  cout << "-------------------------------------------\n";
  cout << "Lanczos time: " << clock.TotalTime().c_str() << endl;
  cout << "-------------------------------------------\n";
  outputfile.close();

  signal_emit(SYSTEM_SIGNAL_GS);
}

/////////////////////////////////////////////////////////////////////////
// build_target_states:
// build target states for the density matrix. 
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::build_target_states(int pos)
{
  cout << "-------------------------------------------\n";
  cout << "Creating target states for the Density Matrix\n";

  if(_calc_gap == 0 && _use_hamis == false && dm_ops.size() == 0){
    if(_ntargets > 1){
      _target.resize(std::min(_nstates,_ntargets)); 
      _target_weight.resize(std::min(_nstates,_ntargets)); 
      for(int i = 0; i < std::min(_nstates,_ntargets); i++){
        _target[i] = _state[i];
      }
    } else {
      _target[0] = gs;
      _target_weight[0] = double(1.);
    }
    return;
  }

  BMatrix<T> rho1;
  Basis basis1;
  BMatrix<T> rho2;
  Basis basis2;
  _target[0] = gs;

  energy.resize(_ntargets);

  typename std::vector<Hami<T> >::iterator hsiter; 
  hsiter = hs.begin();

  for(int i = 0; i < _ntargets-1; i++){
    int i1, i2;

    if(_use_hamis_seed){
      seed = gs;
      if(iter == 1){
        if(dir == RIGHT2LEFT){
          i1 = lattice().ls() - 3;
          i2 = 2;
          read_rho(rho1, basis1, i1, LEFT);
          read_rho(rho2, basis2, i2, RIGHT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, LEFT, use_single_site());
        } else {
          i1 = lattice().ls() - 3;
          i2 = 2;
          read_rho(rho1, basis1, i1, RIGHT);
          read_rho(rho2, basis2, i2, LEFT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, RIGHT, use_single_site());
        }
      } else {
        if(dir == LEFT2RIGHT){
          i1 = iter;
          i2 = lattice().ls() - iter - 2 + 1;
          read_rho(rho1, basis1, i1, LEFT);
          read_rho(rho2, basis2, i2, RIGHT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, LEFT, use_single_site());
        } else {
          i1 = iter;
          i2 = lattice().ls() - iter - 2 + 1;
          read_rho(rho1, basis1, i1, RIGHT);
          read_rho(rho2, basis2, i2, LEFT);
          new_seed(_target[1+i], seed, rho1, rho2, basis1, basis2, RIGHT, use_single_site());
        }
      }
    }
 
    if(gs.size() <= i+1){
      _target[i+1] = gs;
      continue;
    }
 
    if(_calc_gap > 0) {
      _project = true;
      double tol = _lanczos_tol;
      Vector<double> a(1000);
      Vector<double> b(1000);
      Vector<double> g(1000);
      Vector<double> e(1);
      char file[255];
      snprintf(file,255,"vectors_%s.dat",_name);
      int n = _lanczos_maxiter;
    
      _project_states.resize(i+1);
      for(int j = 0; j <= i; j++) _project_states[j] = &_target[j];
   
      cout << "<<<<<<< LANCZOS >>>>>>>\n";
      double ei;
      VectorState<T> gs_left(_state[0]);
      if(this->is_hermitian())
        lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, 1, a, b, n, tol, _use_seed, true, file, true);
      else
        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], gs_left, seed, e[0], ei, n, tol, _use_seed, true);
//        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], seed, e[0], n, tol, _use_seed, true);

      _target[i+1] = _state[0];
  
      _project = false;
      for(int j = 0; j <= i; j++)
        cout << "EXCITED OVERLAP(" << j << ") = " << product(_target[j],_target[1+i]) << endl;
      cout << "E(" << i+1 << ") = " << e[0] << endl;
      energy[i+1] = e[0];
      cout << "GAP(" << i+1 << ") = " << energy[0] << " " << e[0] << " " << e[0] - energy[0] << endl;
    } else {
      double tol = _lanczos_tol;
      Vector<double> a(1000);
      Vector<double> b(1000);
      Vector<double> g(1000);
      Vector<double> e(1);
      char file[255];
      snprintf(file,255,"vectors_%s.dat",_name);
      int n = _lanczos_maxiter;
   
      Hami<T> haux = h;
      h = *hsiter; 
#ifndef USE_PRODUCT_DEFAULT
      cout << "Creating interaction terms" << endl;
      hint(MASK_BLOCK1|MASK_BLOCK2);
      hint(MASK_BLOCK3|MASK_BLOCK4);
      hint(MASK_BLOCK2|MASK_BLOCK3, false);
      build_composite_operators();
#endif
      double ei;
      VectorState<T> gs_left(_state[0]);
      if(this->is_hermitian())
        lanczos<T, System<T>, VectorState<T> >(*this, _state, seed, e, 1, a, b, n, tol, _use_hamis_seed, true, file, true);
      else
        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], gs_left, seed, e[0], ei, n, tol, _use_hamis_seed, true);
//        lanczos_unsymmetric_modified<T, System<T>, VectorState<T> >(*this, _state[0], seed, e[0], n, tol, _use_hamis_seed, true);

      _target[i+1] = _state[0];
      h = haux;
      hsiter++;
      cout << "OVERLAP = " << product(_state[0],seed) << endl;
  
      for(int j = 0; j <= i; j++)
        cout << "EXCITED OVERLAP(" << j << ") = " << product(_target[j],_target[1+i]) << endl;
      cout << "E(" << i+1 << ") = " << e[0] << endl;
      energy[i+1] = e[0];
      cout << "GAP(" << i+1 << ") = " << energy[0] << " " << e[0] << " " << e[0] - energy[0] << endl;
    }
     
  } 
#ifdef USE_K
/*
  if(_in_warmup && ((_grand_canonical & MASK_QN_KX != 0) || (_grand_canonical & MASK_QN_KY != 0))){
    bool final_size = false;
    if(lattice().size() == size()) final_size = true;

    if(final_size){
      typename VectorState<T>::iterator iter;
      for(iter = gs.subspace_begin(); iter != gs.subspace_end(); iter++){
        StateSpace s = *iter;
        state_slice<T> s_slice = gs(s);
        QN sqn = s[1].qn() + s[2].qn() + s[3].qn() + s[4].qn();

        if (sqn != qnt){ 
          for(int i1 = 0; i1 < s[1].dim(); i1++){
            for(int i2 = 0; i2 < s[2].dim(); i2++){
              for(int i3 = 0; i3 < s[3].dim(); i3++){
                for(int i4 = 0; i4 < s[4].dim(); i4++){
                  s_slice(i1,i2,i3,i4) = T(0);
                }
              }
            }
          }
        }
      }
    }
  }
*/
#endif // USE_K

  if(dm_ops.size() > 0) {
    double tweight = 0.5f;
    _target_weight[0] = 0.5f;
    _target.resize(dm_ops.size()+1);
    _target_weight.resize(dm_ops.size()+1);
    typename Hami<T>::const_iterator titer;
    int nt = 1;
    for(titer = dm_ops.begin(); titer != dm_ops.end(); titer++, nt++){
      Term<T> t = *titer;
      BasicOp<T> top = t[0];
 
      const BasicOp<T> *_op = operator()(top);
      if(_op){
        cout << "New target state: " << top.description() << "|gs>\n";
        _target[nt].set_qn_mask(qn+top.dqn, _grand_canonical);
        _target[nt] = VectorState<T> (_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis(),gs.qn()+top.dqn);
        product(*_op, gs, _target[nt], block(top.site()));
        _target[nt] /= sqrt(product(_target[nt],_target[nt]));
 /* 
        int n;
        double e;
        Vector<double> a(1000), b(1000);
        VectorState<T> seed;
        seed = _target[nt];
        
        lanczos<T, System<T>, VectorState<T> >(*this, _target[nt], seed, e, a, b, n, 1.e-5, true, true, "aux.dat");
*/
        _target_weight[nt] = .5f / dm_ops.size();
        tweight += _target_weight[nt];
      }
    }
//    _target_weight /= tweight;
  }
}
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
void
System<T>::delta_rho()
{
  int ls = h.lattice().size();
  Hami<T> &ops = dm_boundary_ops;
  if(ops.size() == 0) return;

  typename Hami<T>::iterator hiter;
  T weight = 0;
  for(hiter = ops.begin(); hiter != ops.end(); hiter++){
    weight += (*hiter).coef();
  }

  typename BMatrix<T>::iterator riter;
  for(riter = rho.begin(); riter != rho.end(); riter++){
    SubMatrix<T> &sm = (*riter);
    sm *= (1.-weight);
  }
  for(hiter = ops.begin(); hiter != ops.end(); hiter++){
    Term<T> &t = *hiter;
    if(get_dir() == LEFT2RIGHT){
      if(t.size() == 2){
        for(int l = 0; l <= get_iter(); l++){
          const BasicOp<T> *op1 = operator()(t[0].set_site(l));
          const BasicOp<T> *op2 = operator()(t[1].set_site(l));
          if(op1 && op2){
            cout << "APPLYING DM OPS " << op1->description() << " " << op2->description() << endl;
            BasicOp<T> drho;
            BasicOp<T> new_op1, new_op2;
            Basis basis(operator[](1).basis(),operator[](2).basis());
            basis.reorder();
            new_op1.dqn = op1->dqn;
            new_op2.dqn = op2->dqn;
            new_op1.resize(basis);
            new_op2.resize(basis);
            new_op1 = 0.;
            new_op2 = 0.;
            int pos = (l == get_iter()) ? RIGHT : LEFT;
            new_operator(new_op1, *op1, rho, basis, pos, T(1), false);
            new_operator(new_op2, *op2, rho, basis, pos, T(1), false);
            BasicOp<T> aux_rho; aux_rho.resize(basis); aux_rho = rho;
            BasicOp<T> aux_rho2; aux_rho2.resize(basis); aux_rho2 = rho;
            aux_rho2 = product(aux_rho,new_op2);
            aux_rho = product(new_op1,aux_rho2);
            drho = aux_rho;
            typename BMatrix<T>::iterator iter;
            T tr = 0.;
            for(iter = drho.begin(); iter != drho.end(); iter++){
              SubMatrix<T> &sm = (*iter);
              SubMatrix<T> *_sm0 = rho.block(sm.qn());
              if(!_sm0) continue;
              SubMatrix<T> &sm0 = *_sm0;
    
              if(sm.rows() != sm.cols()) cout << "ERROR: sm.rows() != sm.cols() : " << sm.rows() << " != " << sm.cols() << "\n";
              if(sm.rows() != sm0.rows()) cout << "ERROR: sm.rows() != sm0.rows() : " << sm.rows() << " != " << sm0.rows() << "\n";
              if(sm.cols() != sm0.cols()) cout << "ERROR: sm.cols() != sm0.cols() : " << sm.cols() << " != " << sm0.cols() << "\n";
              for(int i = 0; i < sm.rows(); i++){
                tr += sm(i,i);
                for(int j = i; j < sm.cols(); j++){
                  sm(j,i) += sm(i,j);
                  sm(i,j) = sm(j,i);
                }
              }
              sm0 += sm*t.coef();
            }
            cout << "DRHO TRACE " << tr << endl;
          }
        }
      } else {  // t.size() == 4
        for(int l1 = 0; l1 < get_iter(); l1++){
          for(int l2 = 0; l2 < get_iter(); l2++){
            BasicOp<T> ref_op1 = BasicOp<T>(t[1].set_site(l2)*t[0].set_site(l1));
            BasicOp<T> ref_op2 = BasicOp<T>(t[3].set_site(l1)*t[2].set_site(l2));
            const BasicOp<T> *op1 = operator()(ref_op1);
            const BasicOp<T> *op2 = operator()(ref_op2);
//            cout << "APPLYING DM OPS: LOOKING FOR " << ref_op1.description() << " " << ref_op2.description() << endl;
            if(op1 && op2){
              cout << "APPLYING DM OPS " << op1->description() << " " << op2->description() << endl;
              BasicOp<T> drho;
              BasicOp<T> new_op1, new_op2;
              Basis basis(operator[](1).basis(),operator[](2).basis());
              basis.reorder();
              new_op1.dqn = op1->dqn;
              new_op2.dqn = op2->dqn;
              new_op1.resize(basis);
              new_op2.resize(basis);
              new_op1 = 0.;
              new_op2 = 0.;
              int pos = LEFT;
              new_operator(new_op1, *op1, rho, basis, pos, T(1), false);
              new_operator(new_op2, *op2, rho, basis, pos, T(1), false);
              BasicOp<T> aux_rho; aux_rho.resize(basis); aux_rho = rho;
              BasicOp<T> aux_rho2; aux_rho2.resize(basis); aux_rho2 = rho;
              aux_rho2 = product(aux_rho,new_op2);
              aux_rho = product(new_op1,aux_rho2);
              drho = aux_rho;
              typename BMatrix<T>::iterator iter;
              T tr = 0.;
              for(iter = drho.begin(); iter != drho.end(); iter++){
                SubMatrix<T> &sm = (*iter);
                SubMatrix<T> *_sm0 = rho.block(sm.qn());
                if(!_sm0) continue;
                SubMatrix<T> &sm0 = *_sm0;
      
                if(sm.rows() != sm.cols()) cout << "ERROR: sm.rows() != sm.cols() : " << sm.rows() << " != " << sm.cols() << "\n";
                if(sm.rows() != sm0.rows()) cout << "ERROR: sm.rows() != sm0.rows() : " << sm.rows() << " != " << sm0.rows() << "\n";
                if(sm.cols() != sm0.cols()) cout << "ERROR: sm.cols() != sm0.cols() : " << sm.cols() << " != " << sm0.cols() << "\n";
                for(int i = 0; i < sm.rows(); i++){
                  tr += sm(i,i);
                  for(int j = i; j < sm.cols(); j++){
                    sm(j,i) += sm(i,j);
                    sm(i,j) = sm(j,i);
                  }
                }
                sm0 += sm*t.coef();
              }
              cout << "DRHO TRACE " << tr << endl;
            }
          }
        }

      }
    } else {
      if(t.size() == 2){
        for(int l = 0; l <= get_iter(); l++){
          const BasicOp<T> *op1 = operator()(t[0].set_site(ls-1-l));
          const BasicOp<T> *op2 = operator()(t[1].set_site(ls-1-l));
          if(op1 && op2){
            cout << "APPLYING " << op1->description() << " " << op2->description() << endl;
            BasicOp<T> drho;
            BasicOp<T> new_op1, new_op2;
            Basis basis(operator[](3).basis(),operator[](4).basis());
            basis.reorder();
            new_op1.dqn = op1->dqn;
            new_op2.dqn = op2->dqn;
            new_op1.resize(basis);
            new_op2.resize(basis);
            new_op1 = 0.;
            new_op2 = 0.;
            int pos = (l == get_iter()) ? LEFT : RIGHT;
            new_operator(new_op1, *op1, rho, basis, pos, T(1), false);
            new_operator(new_op2, *op2, rho, basis, pos, T(1), false);
            BasicOp<T> aux_rho; aux_rho.resize(basis); aux_rho = rho;
            BasicOp<T> aux_rho2; aux_rho2.resize(basis); aux_rho2 = rho;
            aux_rho2 = product(aux_rho,new_op2);
            aux_rho = product(new_op1,aux_rho2);
            drho = aux_rho;
            T tr = 0.;
            typename BMatrix<T>::iterator iter;
            for(iter = drho.begin(); iter != drho.end(); iter++){
              SubMatrix<T> &sm = (*iter);
              SubMatrix<T> *_sm0 = rho.block(sm.qn());
              if(!_sm0) continue;
              SubMatrix<T> &sm0 = *_sm0;
              if(sm.rows() != sm.cols()) cout << "ERROR: sm.rows() != sm.cols() : " << sm.rows() << " != " << sm.cols() << "\n";
              if(sm.rows() != sm0.rows()) cout << "ERROR: sm.rows() != sm0.rows() : " << sm.rows() << " != " << sm0.rows() << "\n";
              if(sm.cols() != sm0.cols()) cout << "ERROR: sm.cols() != sm0.cols() : " << sm.cols() << " != " << sm0.cols() << "\n";
              for(int i = 0; i < sm.rows(); i++){
                tr += sm(i,i);
                for(int j = i; j < sm.cols(); j++){
                  sm(j,i) += sm(i,j);
                  sm(i,j) = sm(j,i);
                }
              }
              sm0 += sm*t.coef();
            }
            cout << "DRHO TRACE " << tr << endl;
          }
        }
      } else { // t.size() == 4
        for(int l1 = 0; l1 < get_iter(); l1++){
          for(int l2 = 0; l2 < get_iter(); l2++){
            BasicOp<T> ref_op1 = BasicOp<T>(t[1].set_site(ls-1-l2)*t[0].set_site(ls-1-l1));
            BasicOp<T> ref_op2 = BasicOp<T>(t[3].set_site(ls-1-l1)*t[2].set_site(ls-1-l2));
            const BasicOp<T> *op1 = operator()(ref_op1);
            const BasicOp<T> *op2 = operator()(ref_op2);
//            cout << "APPLYING DM OPS: LOOKING FOR " << ref_op1.description() << " " << ref_op2.description() << endl;
            if(op1 && op2){
              cout << "APPLYING DM OPS " << op1->description() << " " << op2->description() << endl;
              BasicOp<T> drho;
              BasicOp<T> new_op1, new_op2;
              Basis basis(operator[](3).basis(),operator[](4).basis());
              basis.reorder();
              new_op1.dqn = op1->dqn;
              new_op2.dqn = op2->dqn;
              new_op1.resize(basis);
              new_op2.resize(basis);
              new_op1 = 0.;
              new_op2 = 0.;
              int pos = RIGHT;
              new_operator(new_op1, *op1, rho, basis, pos, T(1), false);
              new_operator(new_op2, *op2, rho, basis, pos, T(1), false);
              BasicOp<T> aux_rho; aux_rho.resize(basis); aux_rho = rho;
              BasicOp<T> aux_rho2; aux_rho2.resize(basis); aux_rho2 = rho;
              aux_rho2 = product(aux_rho,new_op2);
              aux_rho = product(new_op1,aux_rho2);
              drho = aux_rho;
              typename BMatrix<T>::iterator iter;
              T tr = 0.;
              for(iter = drho.begin(); iter != drho.end(); iter++){
                SubMatrix<T> &sm = (*iter);
                SubMatrix<T> *_sm0 = rho.block(sm.qn());
                if(!_sm0) continue;
                SubMatrix<T> &sm0 = *_sm0;
      
                if(sm.rows() != sm.cols()) cout << "ERROR: sm.rows() != sm.cols() : " << sm.rows() << " != " << sm.cols() << "\n";
                if(sm.rows() != sm0.rows()) cout << "ERROR: sm.rows() != sm0.rows() : " << sm.rows() << " != " << sm0.rows() << "\n";
                if(sm.cols() != sm0.cols()) cout << "ERROR: sm.cols() != sm0.cols() : " << sm.cols() << " != " << sm0.cols() << "\n";
                for(int i = 0; i < sm.rows(); i++){
                  tr += sm(i,i);
                  for(int j = i; j < sm.cols(); j++){
                    sm(j,i) += sm(i,j);
                    sm(i,j) = sm(j,i);
                  }
                }
                sm0 += sm*t.coef();
              }
              cout << "DRHO TRACE " << tr << endl;
            }
          }
        }

      }

    }
  }
}

template<class T>
BMatrix<T>
System<T>::density_matrix(int position)
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

/*
  if(!is_hermitian() && _target.size() == 2){
    int mask = position == LEFT ? (MASK_BLOCK1|MASK_BLOCK2) : (MASK_BLOCK3|MASK_BLOCK4);
    dm = state_density_matrix(_target[0], _target[1], mask, true);
  } else {
*/
  {
    for(int nt = 0; nt < _target.size(); nt++){
      if(_target_weight[nt] != 0.0) {
        VectorState<T> &s = _target[nt]; 
#ifdef USE_K
        VectorState<T> s2 = s;
        s2.randomize();
//        s2 *= T(0.000001);
        s2 *= T(_random_coef);
        s += s2; 
        gs = s;
//      if(_in_warmup && lattice().size() != size()) s.randomize();
#endif // USE_K
        T norm = product(s,s);
        cout << "Target vector " << nt << " " << _target_weight[nt] << " " << norm << endl;
        s /= sqrt(norm);
//        if(s.qn_mask() != QN::get_qn_mask()){
//          VectorState<T> s2 = s;
//          s2.set_qn_mask(s2.qn(),QN::get_qn_mask());
//          s2.resize(_b1->basis(),_b2->basis(),_b3->basis(),_b4->basis()); 
//          dm += s2.density_matrix(position)*T(_target_weight[nt]);
//        } else 
        {
          dm += s.density_matrix(position)*T(_target_weight[nt]);
        }
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
System<T>::truncate(int position, int new_size, bool build_targets, bool diagonalize)
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
//      if(is_hermitian()){
        rho.diagonalize(w);
/*
      } else {
        BMatrix<T> u, v;
        rho.svd(w, u, v);
        biorthonormalize(u,v);
        rho_left = u;
        rho = v;
      }
*/
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
//      if(is_hermitian()){
        aux_rho.diagonalize(w);
/*      } else {
        BMatrix<T> u, v;
        rho.svd(w, u, v);
      }
*/
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
//    if(is_hermitian()){
      rho.diagonalize(w);
/*
    } else {
      BMatrix<T> u, v;
      rho.svd(w, u, v);
      biorthonormalize(u,v);
      rho_left = u;
      rho = v;
    }
*/
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
// hint:
// build new interaction Hamiltonian between two blocks.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::hint(int _mask, bool add_local_h, const Hami<T> *_h)
{
  if(!_rotate_terms) return;
#ifdef USE_PRODUCT_DEFAULT
  return;
#endif //USE_PRODUCT_DEFAULT
  if(!_h) _h = &h;
  const Hami<T> &this_hami = *_h; 

  const Block<T> *pb1;
  const Block<T> *pb2;
  BasicOp<T> *_hij;
  size_t pos1, pos2;
  switch(_mask){
    case (MASK_BLOCK1|MASK_BLOCK2):
      pb1 = _b1;
      pb2 = _b2;
      _hij = &h12;
      pos1 = BLOCK1;
      pos2 = BLOCK2;
      break;
    case (MASK_BLOCK2|MASK_BLOCK3):
      pb1 = _b2;
      pb2 = _b3;
      _hij = &h23;
      pos1 = BLOCK2;
      pos2 = BLOCK3;
      break;
    case (MASK_BLOCK3|MASK_BLOCK4):
      pb1 = _b3;
      pb2 = _b4;
      _hij = &h34;
      pos1 = BLOCK3;
      pos2 = BLOCK4;
      break;
    default:
      cout << "*** WARNING: System::hint can not recognize mask\n";
  }

  const Block<T> &b1(*pb1);
  const Block<T> &b2(*pb2);
  BasicOp<T> &hij(*_hij);
//  hij.dqn.kx2() = qnt.kx2();
//  hij.dqn.ky2() = qnt.ky2();

  Basis aux_basis(b1.basis(),b2.basis());
  aux_basis.reorder();
  hij.resize(aux_basis);
  hij = T(0);

// We use indentity matrix for the new operators since 
// we don't want to rotate or truncate them

  BMatrix<T> aux_rho(aux_basis);
  aux_rho.clear();
  typename PackedBasis::const_iterator biter;
  for(biter = aux_basis.subspaces().begin(); biter != aux_basis.subspaces().end(); biter++){
    SubMatrix<T> block(biter->qn(),*biter,*biter);
    block=(I<T>()); 
    aux_rho.push_back(block);
  }

  CTimer clock;
  clock.Start();

  if(add_local_h && apply_hami()){
    clock.Lap();
//    if(verbose() > 0) 
      cout << "OPERATOR H2 " << this_hami.name() << endl;
    const BasicOp<T> *hop = b2(this_hami);
    if(hop) new_operator(hij, *hop, aux_rho, aux_basis, RIGHT, T(1), false);
//    cout << "Lap: " << clock.LapTime().c_str() << endl;

    clock.Lap();
//    if(verbose() > 0) 
      cout << "OPERATOR H1 " << this_hami.name() << endl;
    hop = b1(this_hami);
    if(hop) new_operator(hij, *hop, aux_rho, aux_basis, LEFT, T(1), false);
//   cout << "Lap: " << clock.LapTime().c_str() << endl;
  }

  typename Hami<T>::const_iterator hiter;
  for(hiter = this_hami.begin(); hiter != this_hami.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;
    if(t.type() == TERM_EXTERN) continue;
    if(t.size() == 1 && t.type() != TERM_LOCAL && apply_hami() && t[0].name() != this_hami.name()) {
//      continue;
      const BasicOp<T>&top = t[0];
      size_t ib = block(top.site());
      if((ib == BLOCK2 && _mask == (MASK_BLOCK1|MASK_BLOCK2)) || (ib == BLOCK3 && _mask == (MASK_BLOCK3|MASK_BLOCK4)) || (ib == BLOCK1 && _mask == (MASK_BLOCK1|MASK_BLOCK2) && b1.lattice().size() == 1) || (ib == BLOCK4 && _mask == (MASK_BLOCK3|MASK_BLOCK4) && b2.lattice().size() == 1)){

        bool calc_hc = this_hami.use_hc();
 
        const BasicOp<T>* op= operator()(top);
        if(op){
          if (op->is_diagonal()) calc_hc = false;

//          if(verbose() > 0)
            if(!calc_hc)
              cout << "H12 TERM " << t.name(true) << endl;
            else
              cout << "H12 TERM " << t.name(true) << " + h.c." << endl;

          clock.Lap();

          size_t position = _mask;
          if(_mask == (MASK_BLOCK1|MASK_BLOCK2)) 
            position = ib == BLOCK2 ? RIGHT : LEFT;
          else
            position = ib == BLOCK4 ? RIGHT : LEFT;

          new_operator(hij, *op, aux_rho, aux_basis, position, T(t.coef()));
          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }
    if(t.size() >= 2 && ((t.type() == TERM_PRODUCT && apply_hami()) || (t.type() == TERM_EXTERN && apply_extern()))){
      bool found1 = false;
      bool found2 = false;
      bool found = true;
      Term<T> aux_term1, aux_term2;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        bmask |= mask(block(top.site()));
//        cout << top.description() << " " << pos1 << " " << pos2 << " " << block(top.site()) << endl;
        if(block(top.site()) == pos1){
          aux_term1 *= top;
          found1 = true;
        }
        if(block(top.site()) == pos2){
          aux_term2 *= top;
          found2 = true;
        }
        if(block(top.site()) == BLOCK_NONE){
          found = false;
          break;
        }
      } 
      if(found && (bmask & _mask) == bmask && found1 && found2){ // we have a new composite operator

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
          cout << "ERROR : Operator " << top1.description() << " in term " << t.description() << " not found\n";
          continue;
        }
        if(!op2) {
          cout << "ERROR : Operator " << top2.description() << " in term " << t.description() << " not found\n";
          continue;
        }

        bool calc_hc = this_hami.use_hc();
//        if (op1->is_diagonal() && op2->is_diagonal()) calc_hc = false;
        if (t.is_diagonal()) calc_hc = false;
//        if (t.size() > 2) {
////          calc_hc = false; // you need to include h.c. terms explicitly
//          bool is_diag = true;
//          for(int i = 0; i < t.size(); i++) 
//            if(t[i].is_diagonal()) { is_diag = false; break; }
//          if(is_diag) calc_hc = false;
//        }

        if(op1 && op2){
          clock.Lap();
  
          if(!calc_hc)
            cout << "H12 TERM " << t.description() << endl;
          else
            cout << "H12 TERM " << t.description() << " + h.c." << endl;
          new_operator(hij, *op1, *op2, pos1, pos2, aux_rho, aux_basis, T(t.coef()), calc_hc);

          cout << "Lap: " << clock.LapTime().c_str() << endl;
        }
      }
    }

  }

/* Debugging info

  if (_mask == (MASK_BLOCK3|MASK_BLOCK4)){
    typename BMatrix<T>::iterator ii = hij.begin();
    for(; ii != hij.end(); ii++){
      SubMatrix<T> &kk = *ii;
      cout << "--------------------------------\n";
      cout << kk.qn().n() << " " << kk.qn().sz() << " " << kk.rows() << endl;
      cout << "--------------------------------\n";
      for(int i = 0; i < kk.rows(); i++)
        for(int j = 0; j < kk.cols(); j++)
          cout << i << " " << j << " " << kk(i,j) << endl;
    }
  }

*/
}

/////////////////////////////////////////////////////////////////////////
// rotate:
// build new block and operators, truncating with the density matrix.
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::rotate(int position, Block<T>& b)
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


/////////////////////////////////////////////////////////////////////////
// measure:
// measure physical quantities (mean value of the observables)
/////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::measure_n(size_t n, const VectorState<T> *v)
{
  Lattice &l = h.lattice();

  using namespace std;

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0 = sqrt(product(gs,gs));

  typename Hami<T>::iterator titer;
  for(titer = corr.begin(); titer != corr.end(); titer++){
    Term<T>& t = *titer;

    if(n != 0 && t.size() != n) continue;

    aux = *pv;
    if(gs.qn() != aux.qn()+t.dqn()) continue;

    typename Term<T>::iterator iter;
    iter = t.end();
    iter--;
    bool done = true;
    while(true){
      BasicOp<T> top = (*iter).internals();

      size_t ib = block(top.site());
      if(n == 1 && !_store_products && !operator[](ib).single_site()) {
        done = false;
        break;
      }

      const BasicOp<T>* op = operator()(top);
      if(op){ 
        res.set_qn_mask(aux.qn()+op->dqn, _grand_canonical);
        res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
        res = T(0);

        product(*op, aux, res, ib);
        aux = res;
      } else {
        done = false;
        break;
      }

      if(iter == t.begin()) break;
      iter--;
    }
    if(!done) continue;

    T val = T(0);
    if(gs.qn().equal(res.qn(), _grand_canonical) || _grand_canonical == 0){
      val = t.coef()*product(gs,res)/x/x0;
    }
    cout << setprecision(6) << t.coef();
    for(iter = t.begin(); iter != t.end(); iter++){
      BasicOp<T> &top = *iter;
      if(h.lattice().type() == LATTICE_2D)
        cout << "*" << top.name().c_str() << "(" << l.x(top.site()) << "," << l.y(top.site()) << "," << top.internal_site() << ")";
      if(h.lattice().type() == LATTICE_1D)
        cout << "*" << top.name().c_str() << "(" << l.x(top.site()) << "," << top.internal_site() << ")";
    }
    cout << setprecision(10) << " = " << val << endl;

// cout << setprecision(6) << t.name(true) << " = " << val << endl;
    t.set_value(val);
  }

}

template<class T>
void
System<T>::measure(const VectorState<T> *v)
{
  Lattice &l = h.lattice();

  using namespace std;
  DMTKglobals<T> *globals = get_globals(T(0));

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  VectorState<T> aux_res(*pv);
  T x = product(*pv,*pv);
  T x0 = product(*pv,*pv);


  typename Hami<T>::iterator titer;
//////////////////////////////////////////////////////////////////////////
  for(titer = corr.begin(); titer != corr.end(); titer++){
    Term<T> t = *titer;
    Term<T> &real_t = *titer;

    if(gs.qn() != pv->qn()+t.dqn()) continue;

    if(t.size() == 1 && t[0].is_hami()) {
      res = T(0);
      res = product_default(*this, *pv, &t[0].hami());
      T val = product(*pv,res)/x;  
      cout << setprecision(6) << t[0].hami().name() << " = " << val << endl;
      real_t.set_value(val);
      continue;
    }
    if(t.size() == 1 && !t[0].is_hami()) {
      bool info;
      T val = measure_operator(t[0], v, info);
      if(info) real_t.set_value(val);
      continue;
    }

// We measure terms with 2 operators first 
    if(t.size() == 2){
      BasicOp<T> top2 = t[0].internals();
      BasicOp<T> top1 = t[1].internals();
  
      size_t ib1 = block(top1.site());
      size_t ib2 = block(top2.site());
  
      if(ib1 != BLOCK_NONE && ib2 != BLOCK_NONE){
        if(ib1 != ib2){
            
          const BasicOp<T>* op1 = operator()(top1);
          const BasicOp<T>* op2 = operator()(top2);
  
          if(op1 && op2){
            if(op1->name() != top1.name() || op2->name() != top2.name()) 
                 continue; 
  
            res = T(0);
            product(*op1, *op2, *pv, res, ib1, ib2, T(t.coef()));
            T val = product(*pv,res)/x;
  
            if(h.lattice().type() == LATTICE_1D)
              cout << setprecision(6) << t.coef() << "*" << top1.name().c_str() << "(" << l.x(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << l.x(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
            if(h.lattice().type() == LATTICE_2D)
              cout << setprecision(6) << t.coef() << "*" << top1.name().c_str() << "(" << l.x(top1.site()) << "," << l.y(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << l.x(top2.site()) << "," << l.y(top2.site()) << "," << top2.internal_site() << ") = " << val << endl;
  
//          cout << setprecision(6) << t.name(true) << " = " << val << endl;
  
            real_t.set_value(val);
          }
  
        } else {

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = operator()(aux_op);
  
          if(op){
            res = T(0);
            product(*op, *pv, res, ib1, T(t.coef()));
            T val = product(*pv,res)/x;
  
            if(h.lattice().type() == LATTICE_1D)
              cout << setprecision(6) << t.coef() << "*" << top1.name().c_str() << "(" << l.x(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << l.x(top2.site()) << "," << top2.internal_site() << ") = " << setprecision(10) << val << endl;
            if(h.lattice().type() == LATTICE_2D)
              cout << setprecision(6) << t.coef() << "*" << top1.name().c_str() << "(" << l.x(top1.site()) << "," << l.y(top1.site()) << "," << top1.internal_site() << ")" << top2.name().c_str() << "(" << l.x(top2.site()) << "," << h.lattice().y(top2.site()) << "," << top2.internal_site() << ") = " << setprecision(10) << val << endl;
  
//          cout << setprecision(6) << t.name(true) << " = " << val << endl;

            real_t.set_value(val);
          }
        }
      }
  
    } else {
// We measure terms with arbitrary number of operators
      t = t.reorder();
      Term<T> aux_term[4];
      res = T(0);
      for(int ib = 0; ib < 4; ib++) aux_term[ib].clear();
      bool found = true;
      int bmask = 0;
      for(int i = 0; i < t.size(); i++){
        BasicOp<T> top = t[i].internals();
        if(block(top.site()) == BLOCK_NONE){
          found = false;
          break;
        }
        for(int ib = 1; ib <= 4; ib++){
          if(block(top.site()) == (size_t)ib){
            aux_term[ib-1] *= top;
            bmask |= mask((size_t)ib);
          }
        }
      } 
  
      if(!found) continue;
 
      found = true; 
      const BasicOp<T> *aux_op[4];
      aux_op[0] = aux_op[1] = aux_op[2] = aux_op[3] = NULL;
      for(int ib = 3; ib >= 0 ; ib--){
        if(aux_term[ib].size() != 0){
          BasicOp<T> top;
          if(aux_term[ib].size() == 1)
            top = aux_term[ib][0].internals();
          else
            top = aux_term[ib];
  
          aux_op[ib] = operator()(top);
          if(!aux_op[ib]) {
            cout << "ERROR : Operator " << top.description() << " not found\n";
            found = false;
            break;
          }
        }
      }

      if(!found) continue;

      if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        cout << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *pv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), false, globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        cout << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[1]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), false, globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        cout << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[0]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), false, globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        cout << setprecision(10) << t.description() << " = " << val << endl;
      }
      else if(bmask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        product(*(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, globals);
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val);
        cout << setprecision(10) << t.description() << " = " << val << endl;
      }
      else
      {
        aux = aux_res = *pv;
        for(int ib = 3; ib >= 0 ; ib--){
          if(aux_op[ib]){
            res.set_qn_mask(aux.qn()+aux_op[ib]->dqn, _grand_canonical);
            res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
            res = T(0);
  
            product(*aux_op[ib], aux, res, (size_t)(ib+1));
            aux = res;
          }
        }
        T val = product(*pv,res)/sqrt(x)/sqrt(x0);
        real_t.set_value(val*t.coef());
        cout << setprecision(10) << t.description() << " = " << real_t.value() << endl;
      }
    }
  }
}

template<class T>
T
System<T>::measure_operator(const BasicOp<T> &top, const VectorState<T> *v, bool &info)
{
  using namespace std;
  info = false;

  Lattice &l = h.lattice();

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0= sqrt(product(gs,gs));

  size_t ib = block(top.site());

  cout << "MEASURE OPERATOR " << top.description() << endl;
  const BasicOp<T>* op = operator()(top);

  if(!op) return T(0);
  if(op && op->name() != top.name()) return T(0);

  info = true;
  res.set_qn_mask(aux.qn()+op->dqn, _grand_canonical);
  res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
  res = T(0);

  product(*op, aux, res, ib);

  T val = product(gs,res)/x/x0;

  if(h.lattice().type() == LATTICE_2D)
    cout << top.name().c_str() << "(" << l.x(top.site()) << "," << l.y(top.site()) << "," << top.internal_site() << ")";
  if(h.lattice().type() == LATTICE_1D)
    cout << top.name().c_str() << "(" << l.x(top.site()) << "," << top.internal_site() << ")";
  cout << setprecision(6) << " = " << val << endl;
  return val;
}

template<class T>
VectorState<T>
System<T>::product_operator(const BasicOp<T> &top, const VectorState<T> *v)
{
  using namespace std;

  Lattice &l = h.lattice();

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  T x = sqrt(product(*pv,*pv));
  T x0= sqrt(product(gs,gs));

  size_t ib = block(top.site());

  const BasicOp<T>* op = operator()(top);

  if(op->name() != top.name()) { VectorState<T> res(*pv); res = T(0); return res; }

  res.set_qn_mask(aux.qn()+op->dqn, _grand_canonical);
  res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
  res = T(0);

  if(op->name() != top.name()) return res;
  product(*op, aux, res, ib);

  return res;
}

template<class T>
VectorState<T>
System<T>::product_term(const Term<T> &_t, const VectorState<T> *v)
{
  Lattice &l = h.lattice();

  using namespace std;
  DMTKglobals<T> *globals = get_globals(T(0));

  const VectorState<T> *pv = NULL;
  if(!v)
   pv = &gs;
  else
   pv = v;

  VectorState<T> res(*pv);
  VectorState<T> aux(*pv);
  VectorState<T> aux_res(*pv);
  T x = product(*pv,*pv);
  T x0 = product(*pv,*pv);


  typename Hami<T>::iterator titer;
//////////////////////////////////////////////////////////////////////////
  Term<T> t = _t;
  Term<T> real_t = t;
  res = T(0);

  if(gs.qn() != pv->qn()+t.dqn()) return res;

  if(t.size() == 1 && t[0].is_hami()) {
    res = product_default(*this, *pv, &t[0].hami());
    res *= t.coef();
    return res;
  }
  if(t.size() == 1 && !t[0].is_hami()) {
    res = product_operator(t[0], v);
    res *= t.coef();
    return res;
  }

// We measure terms with 2 operators first 
  if(t.size() == 2){
    BasicOp<T> top2 = t[0].internals();
    BasicOp<T> top1 = t[1].internals();

    size_t ib1 = block(top1.site());
    size_t ib2 = block(top2.site());

    if(ib1 != BLOCK_NONE && ib2 != BLOCK_NONE){
      if(ib1 != ib2){
          
        const BasicOp<T>* op1 = operator()(top1);
        const BasicOp<T>* op2 = operator()(top2);
  
        if(op1 && op2){
          res = T(0);
          if(op1->name() != top1.name() || op2->name() != top2.name()) 
            return res; 

          product(*op1, *op2, *pv, res, ib1, ib2, T(t.coef()));
          return res;
        }
  
      } else {
        res = T(0);
        Term<T> auxt = t;
        auxt.coef() = T(1);
        BasicOp<T> aux_op(auxt);
        const BasicOp<T>* op = operator()(aux_op);

        if(op){
          product(*op, *pv, res, ib1, T(t.coef()));
        } 
        return res;
      }
    } 
  } else {
// We measure terms with arbitrary number of operators
    t = t.reorder();
    Term<T> aux_term[4];
    res = T(0);
    for(int ib = 0; ib < 4; ib++) aux_term[ib].clear();
    bool found = true;
    int bmask = 0;
    for(int i = 0; i < t.size(); i++){
      BasicOp<T> top = t[i].internals();
      if(block(top.site()) == BLOCK_NONE){
        found = false;
        break;
      }
      for(int ib = 1; ib <= 4; ib++){
        if(block(top.site()) == (size_t)ib){
          aux_term[ib-1] *= top;
          bmask |= mask((size_t)ib);
        }
      }
    } 
  
    if(!found) return res;
 
    found = true; 
    const BasicOp<T> *aux_op[4];
    aux_op[0] = aux_op[1] = aux_op[2] = aux_op[3] = NULL;
    for(int ib = 3; ib >= 0 ; ib--){
      if(aux_term[ib].size() != 0){
        BasicOp<T> top;
        if(aux_term[ib].size() == 1)
          top = aux_term[ib][0].internals();
        else
          top = aux_term[ib];

        aux_op[ib] = operator()(top);
        if(!aux_op[ib]) {
          cout << "ERROR : Operator " << top.description() << " not found\n";
          found = false;
          break;
        }
      }
    }

    if(!found) return res;

    if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
      product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, globals);
      return res;
    }
    else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
      product(*(aux_op[0]), *(aux_op[1]), *(aux_op[2]), *pv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), false, globals);
    }
    else if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
      product(*(aux_op[0]), *(aux_op[1]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), false, globals);
    }
    else if(bmask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
      product(*(aux_op[0]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), false, globals);
    }
    else if(bmask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
      product(*(aux_op[1]), *(aux_op[2]), *(aux_op[3]), *pv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), false, globals);
    }
    else
    {
      aux = aux_res = *pv;
      for(int ib = 3; ib >= 0 ; ib--){
        if(aux_op[ib]){
          res.set_qn_mask(aux.qn()+aux_op[ib]->dqn, _grand_canonical);
          res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
          res = T(0);

          product(*aux_op[ib], aux, res, (size_t)(ib+1));
          aux = res;
        }
      }
      res = aux*t.coef();
    }
  }
  return res;
}



//////////////////////////////////////////////////////////////////////
/*
template<class T>
VectorState<T>
product(const System<T>& ss, const VectorState<T>& vv)
{
  VectorState<T> res(vv);
  res = T(0);

  typename Hami<T>::const_iterator hiter;
  for(hiter = ss.h.begin(); hiter != ss.h.end(); hiter++){
    const Term<T>& t = (*hiter);
    if(t.type() == TERM_PRODUCT){
      const BasicOp<T>& top2 = t[0];
      const BasicOp<T>& top1 = t[1];

      size_t b1 = ss.block(top1.site());
      size_t b2 = ss.block(top2.site());
      if(b1 != BLOCK_NONE && b2  != BLOCK_NONE && b1 != b2){
        const BasicOp<T>* op1 = ss(top1);
        const BasicOp<T>* op2 = ss(top2);

        if(ss.verbose() > 0) COUT_PRODUCT_DEFAULT(t);

        if(op1 && op2){
          VectorState<T> aux(*ss._b1,*ss._b2,*ss._b3,*ss._b4,vv.qn()+op2->dqn);
          aux = T(0);
          product(*op2, vv, aux, b2);
          product(*op1, aux, res, b1, T(t.coef()));
        }
      }
    }
  }

  product(*ss._b4->operator()(H<T>()), vv, res, BLOCK4);
  product(*ss._b3->operator()(H<T>()), vv, res, BLOCK2);
  product(*ss._b2->operator()(H<T>()), vv, res, BLOCK3);
  product(*ss._b1->operator()(H<T>()), vv, res, BLOCK1);

  return res;
}
*/

/////////////////////////////////////////////////////////////////////
// product_default:
// This version of product uses the four blocks
/////////////////////////////////////////////////////////////////////
// Original version of product, using 4 blocks
/////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
product_default(System<T>& ss, const VectorState<T>& vv, const Hami<T>* hami, bool only_local)
{
  VectorState<T> res(vv);
  res = T(0);

  const Hami<T> *phami;
  if(hami == NULL)
    phami = &ss.h;
  else
    phami = hami;
  const Hami<T> &h = *phami;

// Test new aproach
// Blocks 1 and 2
/*
  VectorState<T> v12 = vv.condense(MASK_BLOCK4|MASK_BLOCK3);
  VectorState<T> res12(v12);
//  VectorState<T> kk = v12.decondense(MASK_BLOCK4|MASK_BLOCK3,vv);
// for(int i = 0; i < kk.size(); i++) cout << "" << i << " " << vv[i] << " " << kk[i] << endl;
  res12 = T(0);
  if(ss.verbose() > 0) cout << "PRODUCT H12\n";
  product(ss.h34, v12, res12, BLOCK4);
  res = res12.decondense(MASK_BLOCK4|MASK_BLOCK3, res);
*/
/////////////////////////////////////////////////////////////
  typename Hami<T>::const_iterator hiter;
  VectorState<T> aux(res), aux_res(res);
  for(hiter = h.begin(); hiter != h.end(); hiter++){
    const Term<T>& t = (*hiter);
//    if(T(t.coef()) == T(0)) continue;   
    if(t.type() == TERM_EXTERN && !ss.apply_extern()) continue;
    if(t.size() > 2){  // you have to include h.c. terms explicitly!
      bool save_use_composite = ss.use_composite();
      ss.set_use_composite(false);
      AuxTerm<T> new_term;
      new_term.t = t;

      if(ss.verbose() > 0)
        cout << "PRODUCT TERM " << t.description() << endl;

      int bmask = 0;
      bool found = true;
      for(int i = 0; i < t.size(); i++){
        const BasicOp<T> top = t[i].internals();
        if(ss.block(top.site()) == BLOCK_NONE){
          found = false;
//          cout << "WARNING: *** Operator " << top.description() << " " << top.site() << " not found\n";
          break;
        }
        for(int ib = 1; ib <= 4; ib++){
          if(ss.block(top.site()) == (size_t)ib){
            new_term.cterm[ib-1] *= top;
            new_term.pos[ib-1] = i;
            bmask |= mask((size_t)ib);
          }
        }
      } 
      new_term.mask = bmask;
      new_term.nblocks = 0; 
      new_term.nblocks += ((bmask & MASK_BLOCK1) == MASK_BLOCK1 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK2) == MASK_BLOCK2 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK3) == MASK_BLOCK3 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK4) == MASK_BLOCK4 ? 1 : 0);
 
      if(found && bmask != MASK_BLOCK1 && bmask != MASK_BLOCK2 && bmask != MASK_BLOCK3 && bmask != MASK_BLOCK4){
// && ((t.type() != TERM_EXTERN && bmask != (MASK_BLOCK1|MASK_BLOCK2) && bmask != (MASK_BLOCK2|MASK_BLOCK3) && bmask != (MASK_BLOCK3|MASK_BLOCK4)) || t.type() == TERM_EXTERN)){
        if(bmask == (MASK_BLOCK1|MASK_BLOCK2)) {
          new_term.b1 = BLOCK1; new_term.b2 = BLOCK2;
        }else if(bmask == (MASK_BLOCK2|MASK_BLOCK3)) {
          new_term.b1 = BLOCK2; new_term.b2 = BLOCK3;
        }else if(bmask == (MASK_BLOCK3|MASK_BLOCK4)) {
          new_term.b1 = BLOCK3; new_term.b2 = BLOCK4;
        }else if(bmask == (MASK_BLOCK1|MASK_BLOCK3)) {
          new_term.b1 = BLOCK1; new_term.b2 = BLOCK3;
        }else if(bmask == (MASK_BLOCK1|MASK_BLOCK4)) {
          new_term.b1 = BLOCK1; new_term.b2 = BLOCK4;
        }else if(bmask == (MASK_BLOCK2|MASK_BLOCK4)) {
          new_term.b1 = BLOCK2; new_term.b2 = BLOCK4;
        }

        bool found_op = true;
        for(int ib = 3; ib >= 0 ; ib--){
          if(new_term.cterm[ib].size() != 0){
            if(new_term.cterm[ib].size() == 1)
              new_term.top[ib] = new_term.cterm[ib][0];
            else
              new_term.top[ib] = new_term.cterm[ib];
 
            new_term.op[ib] = ss(new_term.top[ib]);
            if(!new_term.op[ib]) {
              cout << "ERROR : Operator " << new_term.top[ib].description() << " not found\n";
              found_op = false;
            }
          }
        }
        if(found_op){
          if(ss.verbose() > 0)
            cout << "APPLYING " << t.description() << endl;
 
          product_term(ss, new_term, t, vv, res);
        }
      }
      ss.set_use_composite(save_use_composite);

/*
      if(bmask != MASK_BLOCK1 && bmask != MASK_BLOCK2 && bmask != MASK_BLOCK3 && bmask != MASK_BLOCK4 && found){ 
        aux = aux_res = vv;
        for(int ib = 3; ib >= 0 ; ib--){
          if(aux_term[ib].size() != 0){
            BasicOp<T> top;
            if(aux_term[ib].size() == 1)
              top = aux_term[ib][0].internals();
            else
              top = aux_term[ib];
    
            const BasicOp<T>* _op = ss(top);
            if(!_op) {
              cout << "ERROR : Operator " << top.description() << " not found\n";
            }
    
            if(_op){
              if(ss.verbose() > 0)
                cout << "APPLYING " << top.description() << endl;
    
              aux_res.set_qn_mask(aux.qn()+_op->dqn, ss.qn_mask());
              aux_res.resize(aux.b1(),aux.b2(),aux.b3(),aux.b4());
              aux_res = T(0);
    
              product(*_op, aux, aux_res, (size_t)(ib+1));
              aux = aux_res;
            }
          }
        }
        res += aux_res * t.coef(); 
      }
*/
    }  
    if(t.size() == 2 && t.type() == TERM_PRODUCT){
      const BasicOp<T>& top1 = t[0];
      const BasicOp<T>& top2 = t[1];

      bool calc_hc = h.use_hc();
      if (top1.is_diagonal() && top2.is_diagonal()) calc_hc = false;

      size_t b1 = ss.block(top1.site());
      size_t b2 = ss.block(top2.site());

      if(only_local && mask(b1,b2) == (MASK_BLOCK3|MASK_BLOCK4)) continue;
      if(only_local && mask(b1,b2) == (MASK_BLOCK1|MASK_BLOCK2)) continue;

      if(b1 != BLOCK_NONE && b2 != BLOCK_NONE){
        if(b1 != b2){
          if(ss.verbose() > 0)
            cout << "PRODUCT TERM " << t.description() << endl;

          const BasicOp<T>* op1 = ss(top1);
          const BasicOp<T>* op2 = ss(top2);

          if(op1 && op2){
            if(ss.verbose() > 0)
              COUT_PRODUCT(t);
  
            product(*op2, *op1, vv, res, b2, b1, T(t.coef()), calc_hc);
          }

        } else if(t.type() == TERM_EXTERN && ss.apply_extern()) {

          if(only_local && (b1 == BLOCK1 || b1 == BLOCK4)) continue;

          BasicOp<T> aux_op(t);
          const BasicOp<T>* op = ss(aux_op);

          if(ss.verbose() > 0)
            if(calc_hc) 
              cout << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ") + h.c." << endl;
            else 
              cout << "PRODUCT EXTERN " << t.coef() << "*" << top1.name() << "(" << top1.site() << ")" << top2.name() << "(" << top2.site() << ")" << endl;
  
          if(op){
            product(*op, vv, res, b1, T(t.coef()));
            // obviously, op1 has to preserve the quantum numbers,
            // or we are working in the grand canonical 
            if(calc_hc) 
              product(*op, vv, res, b1, T(t.coef()), calc_hc);
          }

        }
      }
    }
    else if(t.size() == 1){

      const BasicOp<T>& top = t[0];
      size_t b1 = ss.block(top.site());
      if(b1 == BLOCK_NONE) continue;

      if(top.is_hami() && top.name() != h.name()){
        res += product_default(ss, vv, &top.hami());
        continue;
      }
      if(top.name() == "I"){
        VectorState<T> aux = vv;
        aux *= t.coef();
        res += aux;
        continue;
      }

//      if((b1 == BLOCK1 || b1 == BLOCK4) && t.type() != TERM_LOCAL && t.type() != TERM_EXTERN) continue;
      if(t.type() == TERM_LOCAL || t.type() == TERM_EXTERN || ss[b1].single_site()) { 

        if(b1 != BLOCK_NONE && top.name() != h.name()){
          const BasicOp<T>* op1 = ss(top);
          bool calc_hc = h.use_hc();
          if (op1 && op1->is_diagonal()) calc_hc = false;
  
          if(ss.verbose() > 0)
            if(calc_hc)
              cout << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ") + h.c." << endl;
            else
              cout << "PRODUCT LOCAL " << t.coef() << "*" << top.name() << "(" << top.site() << ")" << endl;
  
          if(op1){
            product(*op1, vv, res, b1, T(t.coef()));
            // obviously, op1 has to preserve the quantum numbers,
            // or we are working in the grand canonical
            if(calc_hc) {
              product(*op1, vv, res, b1, T(t.coef()), calc_hc);
            }
          }
        }
      }
    }
  }

  const char *name = h.name().c_str();
  const BasicOp<T>* local_op4 = ss._b4->operator()(H<T>().set_name(name));
  if(local_op4 && local_op4->name() == h.name()){ 
    if(ss.verbose() > 0) 
      cout << "PRODUCT " << name << "4\n";
    product(*local_op4, vv, res, BLOCK4);
  } else {
#ifdef WITH_WARNINGS
    if(h.name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
#endif
  }
  const BasicOp<T>* local_op3 = ss._b3->operator()(H<T>().set_name(name));
  if(local_op3 && local_op3->name() == h.name()){ 
    if(ss.verbose() > 0) 
      cout << "PRODUCT " << name << "3\n";
    product(*local_op3, vv, res, BLOCK3);
  } else {
#ifdef WITH_WARNINGS
    if(h.name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
#endif
  }
  const BasicOp<T>* local_op2 = ss._b2->operator()(H<T>().set_name(name));
  if(local_op2 && local_op2->name() == h.name()){ 
    if(ss.verbose() > 0) 
      cout << "PRODUCT " << name << "2\n";
    product(*local_op2, vv, res, BLOCK2);
  } else {
#ifdef WITH_WARNINGS
    if(h.name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
#endif
  }
  const BasicOp<T>* local_op1 = ss._b1->operator()(H<T>().set_name(name));
  if(local_op1 && local_op1->name() == h.name()){ 
    if(ss.verbose() > 0) 
      cout << "PRODUCT " << name << "1\n";
    product(*local_op1, vv, res, BLOCK1);
  } else {
#ifdef WITH_WARNINGS
    if(h.name() == "H")
      cout << "WARNING: Add local (site) Hamiltonian term\n";
//    else
//      cout << "You can ignore the warning\n";
#endif
  }

  return res;
}

/////////////////////////////////////////////////////////////////////
// product_additive:
// This producting additive operators in the hamiltonian
// like S+(0)(S-(1)+ S-(2)+ S-(3))
/////////////////////////////////////////////////////////////////////
template<class T>
VectorState<T>
product_additive(const System<T>& ss, const VectorState<T>& vv, const Term<T>& term)
{
  VectorState<T> res(vv);
  res = T(0);

  const BasicOp<T>& top1 = term[0];
  const BasicOp<T>& top2 = term[1];

  if(!top2.is_hami()) {
    cout << "WARNING: (product_additive) Operator is not additive\n";
    return res;
  }
  bool calc_hc = ss.h.use_hc();
  if (term.is_diagonal()) calc_hc = false;

  CTimer aux_timer;
  aux_timer.Start();
  if(ss.verbose() > 0){
    cout << "PRODUCT ADDITIVE ";
    COUT_PRODUCT(term);
  }
  
  const Hami<T> &h = top2.hami();
  size_t b1 = ss.block(top1.site());
  const BasicOp<T>* op1 = ss(top1);
  T coef = T(term.coef());
  if(!op1 && (b1 == BLOCK2 || b1 == BLOCK3)) {
    BasicOp<T> _op = top1.internals();
    op1 = ss[b1](_op.set_site(0));
  }

  if(ss.verbose()  > 0 && !op1)
    cout << "WARNING: " << top1.description() << " not found " << endl;
  
  if(op1){
    size_t mm = mask(b1); // auxiliary mask
    typename Hami<T>::const_iterator hiter;
    for(hiter = h.begin(); hiter != h.end(); hiter++){
      const Term<T>& t = (*hiter);
  
      BasicOp<T> one_op = t[0].internals();
      size_t b2 = ss.block(one_op.site());
  
      if(!(mm & mask(b2))){
        mm |= mask(b2);
        const BasicOp<T>* op2 = ss[b2](top2);

        if(!op2 && (b2 == BLOCK2 || b2 == BLOCK3 || ss[b2].lattice().size() == 1)){
          BasicOp<T> _op = one_op.internals();
          op2 = ss[b2](_op.set_site(0));
          coef *= t.coef();
        } 
        if(ss.verbose() > 0 && !op2)
          cout << "WARNING: " << one_op.description() << " not found " << endl;
        if(op2) {
          if(ss.verbose() > 0){
            cout << "BLOCKS " << b1 << " " << b2 << " ";
            COUT_PRODUCT(term);
          }
          aux_timer.Lap();
          product(*op1, *op2, vv, res, b1, b2, coef, calc_hc);
          if(ss.verbose() > 0)  
            cout << "Lap time: " << aux_timer.LapTime() << endl;
        }
      }    
      if(mm == 15 && ss.verbose() == 0) break;
    }
  }
  if(ss.verbose() > 0)  
    cout << "CPU Time (product_additive): " << aux_timer.TotalTime() << endl;
  return res;
}
/////////////////////////////////////////////////////////////////////
// product:
// This version of product uses the interaction blocks
/////////////////////////////////////////////////////////////////////
#ifdef WITH_PTHREADS

pthread_mutex_t mutex_res;

template<class T>
class ThreadArg
{
  public:
    const ProductTerm<T> *term;
    const VectorState<T> *vv;
    VectorState<T> *res;
    int index;

    ThreadArg() {}
}; 

template<class T>
void *my_thread_func(void *_arg)
{
  ThreadArg<T> *arg = (ThreadArg<T> *)_arg;
  const ProductTerm<T> &term = *arg->term;
  const VectorState<T> &vv = *arg->vv;
  VectorState<T> res = *arg->res;
  res = T(0);
  int nproc = sysconf(_SC_NPROCESSORS_ONLN);

  cout << "THREAD INDEX " << arg->index << endl;
//  
/*
  cpu_set_t mask;
  CPU_ZERO(&mask);
  CPU_SET(arg->index%nproc, &mask);
*/
/*
  if(pthread_setaffinity_np(pthread_self(), sizeof(mask), &mask) < 0){
    perror("pthread_setaffinity_np");
  }
*/

  product_term(term, vv, res, (MASK_PRODUCT_DEFAULT|MASK_PRODUCT_HC), globals);

//  pthread_mutex_lock(&mutex_res);
  VectorState<T> &real_res = *arg->res;
//cout << "HOLA " << arg->index << " " << real_res.size() << " " << res.size() << endl;
  real_res += res;
//  pthread_mutex_unlock(&mutex_res);

  pthread_exit(NULL);
  return NULL;
}
#endif // WITH_PTHREADS

#ifdef USE_PRODUCT_DEFAULT
template<class T>
VectorState<T>
product(System<T>& ss, const VectorState<T>& vv)
{
  VectorState<T> res(vv);
  res = T(0);
  res = product_default(ss, vv);

  if(ss.project()){
    for(int i = 0; i < ss._project_states.size(); i++){
      VectorState<T> aux = *ss._project_states[i];
      T p = product(aux,vv);
      res = res + aux * p * T(100000);
    }
  }

  return res;
}
#else // USE_PRODUCT_DEFAULT
////////////////////////////////////////////////////////////////////////
// build_composite:
////////////////////////////////////////////////////////////////////////
template<class T>
void
System<T>::build_composite_operators(const Hami<T> *this_hami)
{
  CTimer aux_timer;
  aux_timer.Start();

  if(!this_hami) aux_terms.clear();

  const Hami<T> *ph = &h;
  if(this_hami) ph = this_hami;
  const Hami<T> &_h = *ph;

  typename Hami<T>::const_iterator hiter;
  for(hiter = _h.begin(); hiter != _h.end(); hiter++){
    const Term<T> t = *hiter;
    if(_use_coef_tol && fabs(real(t.coef())) <= _coef_tol) continue;

    if(t.size() == 1 && t[0].is_hami()) {
      build_composite_operators(&t[0].hami());
    }

    if(t.size() >= 2 && ((apply_hami() && t.type() == TERM_PRODUCT) || (apply_extern() && t.type() == TERM_EXTERN))){  
      AuxTerm<T> new_term;
      new_term.t = (*hiter);

      int bmask = 0;
      bool found = true;
      for(int i = 0; i < t.size(); i++){
        const BasicOp<T> top = t[i].internals();
        if(block(top.site()) == BLOCK_NONE){
          found = false;
//          cout << "WARNING: *** Operator " << top.description() << " " << top.site() << " not found\n";
          break;
        }
        for(int ib = 1; ib <= 4; ib++){
          if(block(top.site()) == (size_t)ib){
            new_term.cterm[ib-1] *= top;
            new_term.pos[ib-1] = i;
            bmask |= mask((size_t)ib);
          }
        }
      } 
      new_term.mask = bmask;
      new_term.nblocks = 0; 
      new_term.nblocks += ((bmask & MASK_BLOCK1) == MASK_BLOCK1 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK2) == MASK_BLOCK2 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK3) == MASK_BLOCK3 ? 1 : 0);
      new_term.nblocks += ((bmask & MASK_BLOCK4) == MASK_BLOCK4 ? 1 : 0);
// WARNING: FIXME
if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4) || bmask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)) new_term.nblocks = 3;
if(bmask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3) || bmask == (MASK_BLOCK4|MASK_BLOCK3|MASK_BLOCK2)) new_term.nblocks = 3;

      bool found_op = true;
      if(found && bmask != MASK_BLOCK1 && bmask != MASK_BLOCK2 && bmask != MASK_BLOCK3 && bmask != MASK_BLOCK4 && ((t.type() != TERM_EXTERN && bmask != (MASK_BLOCK1|MASK_BLOCK2) && bmask != (MASK_BLOCK2|MASK_BLOCK3) && bmask != (MASK_BLOCK3|MASK_BLOCK4)) || t.type() == TERM_EXTERN || this_hami != NULL)){
        for(int ib = 3; ib >= 0 ; ib--){
          if(new_term.cterm[ib].size() != 0){
            if(new_term.cterm[ib].size() == 1)
              new_term.top[ib] = new_term.cterm[ib][0];
            else
              new_term.top[ib] = new_term.cterm[ib];
 
            new_term.op[ib] = operator()(new_term.top[ib]);
            if(!new_term.op[ib]) {
              cout << "ERROR : Operator " << new_term.top[ib].description() << " not found\n";
              found_op = false;
            }
          }
        }

        if(found_op) aux_terms.push_back(new_term);
      }
    }
  }

  cout << "COMPOSITE TERMS: " << aux_terms.size() << endl;

  typename std::list<AuxTerm<T> >::iterator titer;
  for(titer = aux_terms.begin(); titer != aux_terms.end(); titer++){
    AuxTerm<T>& auxt = (*titer);
    const Term<T>& t = (auxt.t);
//    T sign = T(1);

    if(!auxt.done){

      if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK4)){
        auxt.b1 = BLOCK4, auxt.b2 = BLOCK1; 
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b1 = BLOCK1; auxt.b2 = BLOCK4; 
        } 
        if(auxt.cterm[0].size() == 1) { auxt.b1 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[3].size() == 1) { auxt.b1 = BLOCK4; auxt.b2 = BLOCK1; }
      } else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK4; 
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK1; 
//------------------------------------------------------------------
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1; 
      } else if(auxt.mask == (MASK_BLOCK3|MASK_BLOCK4)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK4; 
      } else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK3; 
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1; auxt.b3 = BLOCK3; 
      } else if(auxt.mask == (MASK_BLOCK4|MASK_BLOCK2|MASK_BLOCK3)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK4; auxt.b3 = BLOCK3; 
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
        auxt.b1 = BLOCK2; auxt.b2 = BLOCK1; auxt.b3 = BLOCK4; 
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; 
        } 
        if(auxt.cterm[3].size() > 1) { auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[0].size() > 1) { auxt.b3 = BLOCK4; auxt.b2 = BLOCK1; }
      } else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
        auxt.b1 = BLOCK3; auxt.b2 = BLOCK1; auxt.b3 = BLOCK4; 
        if(operator[](BLOCK1).lattice().size() < operator[](BLOCK4).lattice().size()){
          auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; 
        } 
        if(auxt.cterm[3].size() > 1) { auxt.b3 = BLOCK1; auxt.b2 = BLOCK4; }
        if(auxt.cterm[0].size() > 1) { auxt.b3 = BLOCK4; auxt.b2 = BLOCK1; }
      } else
        continue;

      if(!use_composite()) continue;

      auxt.ref_op[0] = (auxt.op[int(auxt.b1)-1]);
      auxt.sum_op = *(auxt.op[int(auxt.b2)-1]);
      if(auxt.nblocks >= 3) auxt.ref_op[1] = (auxt.op[int(auxt.b3)-1]);

//      if(auxt.pos[auxt.b1-1] < auxt.pos[auxt.b2-1] && auxt.b1 > auxt.b2 && (auxt.ref_op[0]->fermion() || auxt.sum_op.fermion())) sign = T(-1);

      typename BMatrix<T>::iterator biter;
      for(biter = auxt.sum_op.begin(); biter != auxt.sum_op.end(); biter++){
        SubMatrix<T> &sm = (*biter);
        sm *= t.coef();
      }
      auxt.done = true;

//      if(verbose() > 0)
      {
        cout << "COMPOSITE TERM " << auxt.t.description() << " " << auxt.nblocks << " " << auxt.ref_op[0]->description();
        if(auxt.nblocks >=3 ) cout << " " << auxt.ref_op[1]->description();
        cout << endl;
      }   

      int nterm = 1;
      typename std::list<AuxTerm<T> >::iterator titer2;
      for(titer2 = titer; titer2 != aux_terms.end(); titer2++){
        AuxTerm<T>& auxt2 = (*titer2);
        const Term<T>& t2 = auxt2.t;
//        sign = T(1);
// FIXME
        if(!auxt2.done && auxt2.mask == auxt.mask && auxt.t.is_diagonal() == auxt2.t.is_diagonal() && auxt.t.dqn().equal(auxt2.t.dqn(),MASK_QN_ALL) && ( (auxt.nblocks == 2 && auxt2.op[int(auxt.b1)-1]->dqn.equal(auxt.ref_op[0]->dqn,MASK_QN_ALL) && auxt2.op[int(auxt.b1)-1]->internals() == auxt.ref_op[0]->internals()) || (auxt.nblocks == 3 && auxt2.op[int(auxt.b1)-1]->dqn.equal(auxt.ref_op[0]->dqn,MASK_QN_ALL) && auxt2.op[int(auxt.b3)-1]->dqn.equal(auxt.ref_op[1]->dqn,MASK_QN_ALL) && auxt2.op[int(auxt.b1)-1]->internals() == auxt.ref_op[0]->internals() && auxt2.op[int(auxt.b3)-1]->internals() == auxt.ref_op[1]->internals()))){
//        if(!auxt2.done && auxt2.mask == auxt.mask && auxt.t.dqn() == auxt2.t.dqn() && (auxt2.op[int(auxt.b1)-1]->internals() == auxt.ref_op[0]->internals() && auxt.nblocks == 2)) {
          auxt2.done = true;
          if(verbose() > 0)
            cout << "COMPOSITE TERM SUM " << auxt2.t.description() << " " << auxt2.t.dqn() << endl;
          auxt2.mask = 0;
//          if(auxt2.pos[auxt.b1-1] < auxt2.pos[auxt.b2-1] && auxt2.b1 > auxt2.b2 && (auxt.ref_op[0]->fermion() || auxt.sum_op.fermion())) sign = T(-1);
          typename BMatrix<T>::iterator biter1;
          const BMatrix<T>& bm2 = *auxt2.op[int(auxt.b2)-1];

//////////////////////////////////////////////
// If this is the second term in the sum, 
// initialize sum_op to include all the blocks
//////////////////////////////////////////////
          if(nterm == 1){
            BasicOp<T> aux_sum;
            aux_sum = auxt.sum_op.internals();
            const Block<T> *_block;
            valarray<const Block<T>* > b(4);
            b[0] = _b1;
            b[1] = _b2;
            b[2] = _b3;
            b[3] = _b4;
            _block = b[int(auxt.b2)-1];
            aux_sum.resize(_block->basis().subspaces());
            for(biter1 = auxt.sum_op.begin(); biter1 != auxt.sum_op.end(); biter1++){
              SubMatrix<T> &_sm1 = (*biter1);
              SubMatrix<T> *_sm2 = aux_sum.block(_sm1.qn());
              *_sm2 = _sm1;
            }
            auxt.sum_op = aux_sum;
          }
//////////////////////////////////////////////

          nterm++;
          typename BMatrix<T>::const_iterator biter2;
          for(biter2 = bm2.begin(); biter2 != bm2.end(); biter2++){
            const SubMatrix<T> &_sm2 = (*biter2);
            SubMatrix<T> *_sm1 = auxt.sum_op.block(_sm2.qn());
            if(_sm1){
              Matrix<T> sm2 = _sm2;
              Matrix<T> &sm1 = *_sm1;
              sm2 *= t2.coef();
              sm1 += sm2; 
            } else {
              cout << "WARNING: Block not found " << auxt.t.description() << " " << auxt2.t.description() << endl;
;
            }
          }
          titer2 = aux_terms.erase(titer2);
        }
      }
    } 
  }
  cout << "COMPOSITE TERMS: " << aux_terms.size() << endl;
  cout << "Build composite time: " << aux_timer.TotalTime().c_str() << endl;
}


template<class T>
VectorState<T>
product_composite(System<T>& ss, const VectorState<T>& vv, int mask_hc)
{
  VectorState<T> res(vv);
  res = T(0);
  CTimer aux_timer;
  DMTKglobals<T> *globals = get_globals(T(0));

  aux_timer.Start();
  VectorState<T> v12 = vv.condense(MASK_BLOCK1|MASK_BLOCK2);
  VectorState<T> v1234 = v12.condense(MASK_BLOCK3|MASK_BLOCK4);
  VectorState<T> res1234(v1234);
  if(ss.verbose() > 0)
    cout << "CONDENSE TIME= " << aux_timer.LapTime().c_str() << endl;
  aux_timer.Lap();
  res1234 = T(0);
  if(ss.verbose() > 0)
    cout << "PRODUCT H12\n";
  product(ss.h12, v1234, res1234, BLOCK1, T(1), false, globals);
  if(ss.verbose() > 0)
    cout << aux_timer.LapTime().c_str() << endl;
  aux_timer.Lap();
  if(ss.verbose() > 0)
    cout << "PRODUCT H34\n";
  product(ss.h34, v1234, res1234, BLOCK4, T(1), false, globals);
  if(ss.verbose() > 0)
    cout << aux_timer.LapTime().c_str() << endl;
  aux_timer.Lap();
  VectorState<T> res12 = res1234.decondense(MASK_BLOCK3|MASK_BLOCK4, v12);
  res = res12.decondense(MASK_BLOCK1|MASK_BLOCK2, vv);
  if(ss.verbose() > 0)
    cout << "DECONDENSE TIME= " << aux_timer.LapTime().c_str() << endl;

// Blocks 2 and 3

  if(ss.verbose() > 0)
    cout << "PRODUCT H23\n";
  VectorState<T> v23 = vv.condense(MASK_BLOCK2|MASK_BLOCK3);
  VectorState<T> res23(v23);
  res23 = T(0);
  product(ss.h23, v23, res23, BLOCK2, T(1), false, globals);
  res += res23.decondense(MASK_BLOCK2|MASK_BLOCK3, vv);

#ifdef WITH_PTHREADS
  pthread_t p[NUM_THREADS];
  pthread_attr_t attr;
  Vector<ThreadArg<T> > arg(NUM_THREADS);
  for(int i = 0; i < arg.size(); i++) arg[i].index = i;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  int nthread = 0;
  int nactive = 0;
  typename std::vector<ProductTerm<T> >::iterator term_iter;
  for(term_iter = ss.product_terms.begin(); term_iter != ss.product_terms.end(); term_iter++){
    int pindex = nthread % NUM_THREADS;
    
    ProductTerm<T> &pterm = *term_iter;
    arg[pindex].term = &pterm;
    arg[pindex].vv = &vv;
    arg[pindex].res = &res;

    nactive++;
    pthread_create(&p[pindex], &attr, my_thread_func<T>, (void *)&arg[pindex]);
    if(pindex == NUM_THREADS-1 || nthread == arg.size()-1) {
      for(int n = 0; n < nactive; n++) {
        pthread_join(p[n],NULL);
      }
      nactive = 0;
    }

  }
  pthread_attr_destroy(&attr);

#else
  typename std::vector<ProductTerm<T> >::iterator term_iter;
  for(term_iter = ss.product_terms.begin(); term_iter != ss.product_terms.end(); term_iter++){
    ProductTerm<T> &pterm = *term_iter;
    product_term(pterm, vv, res, mask_hc, globals);
  }
#endif

 if(ss.project()){
    for(int i = 0; i < ss._project_states.size(); i++){
      VectorState<T> aux = *ss._project_states[i];
      T p = product(aux,vv);
      res = res + aux * p*T(100000);
    }
  }

  if(ss.verbose() > 0)
    cout << "CPU Time (product): " << aux_timer.TotalTime() << endl;

  return res;
}


template<class T>
void
init_terms_composite(System<T>& ss, const VectorState<T>& vv)
{
  if(!ss.rotate_terms()) return;

  VectorState<T> res(vv);
  res = T(0);
  std::vector<ProductTerm<T> > new_terms;
  cout << "Initializing product terms\n";
  ss.product_terms.clear();
// Interactions between blocks (1-3,2-4,2-3), and external and local terms
  VectorState<T> aux(res), aux_res(res);
  typename Hami<T>::const_iterator hiter;
//////////////////////////////////////////////////////////////////////////
// product-4+  products of 4 or more operators
//////////////////////////////////////////////////////////////////////////
  int nterm = 0;
  typename std::list<AuxTerm<T> >::const_iterator titer;
  for(titer = ss.aux_terms.begin(); titer != ss.aux_terms.end(); titer++){
    const AuxTerm<T> &auxt = (*titer);

    Term<T> t = (auxt.t);

    init_term_composite(ss, auxt, t, vv, res);
  }
//////////////////////////////////////////////////////////////////////////
// product-hami
//////////////////////////////////////////////////////////////////////////
  for(hiter = ss.h.begin(); hiter != ss.h.end(); hiter++){
    const Term<T>& t = *hiter;
     
    if(t.size() == 1 && t[0].is_hami()){
      const Hami<T>& h = t[0].hami();
      const BasicOp<T> *op = ss._b1->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK1);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b2->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK2);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b3->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK3);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
      op = ss._b4->operator()(t[0]);
      if(op) {
        new_terms = get_product_terms(*op, vv, res, BLOCK4);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    }
  }
//////////////////////////////////////////////////////////////////////////
// product-1 (only for TERM_EXTERN || TERM_LOCAL || IDENTITY)
//////////////////////////////////////////////////////////////////////////
  for(hiter = ss.h.begin(); hiter != ss.h.end(); hiter++){
    const Term<T>& t = *hiter;
     
    if(ss.use_coef_tol() && fabs(real(t.coef())) <= ss.coef_tol()) continue;

    if(t.size() == 1 && t[0].name() == "I"){
      VectorState<T> aux = vv;
      aux *= t.coef();
      res += aux;
      continue;
    }

    if((t.type() == TERM_EXTERN && ss.apply_extern() && t.size() == 1 && t[0].name() != ss.h.name()) || (t.type() == TERM_LOCAL && t[0].name() != ss.h.name())){
      size_t bmask = 0;
      for(int i = 0; i < t.size(); i++)
        bmask |= mask(ss.block(t[i].site())); 

      if(bmask == MASK_BLOCK1 || bmask == MASK_BLOCK2 || bmask == MASK_BLOCK3 || bmask == MASK_BLOCK4){
        BasicOp<T> top = t[0].internals();
        const BasicOp<T>* _op = ss(top);
        if(!_op) {
          cout << "ERROR EXTERN: Operator " << top.description() << " not found\n";
        }
  
        if(_op){
          bool calc_hc = ss.h.use_hc();
          if (_op && _op->is_diagonal()) calc_hc = false;

          size_t pos = ss.block(t[0].site()); 
          if(_op){
            new_terms = get_product_terms(*_op, vv, res, pos, T(t.coef()), false);
            ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
            // obviously, op1 has to preserve the quantum numbers,
            // or we are working in the grand canonical 
            if(calc_hc) {
              new_terms = get_product_terms(*_op, vv, res, pos, T(t.coef()), calc_hc);
              ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
            }
          }
        }
      }
    }
  } 
}

#endif //USE_PRODUCT_DEFAULT
//////////////////////////////////////////////////////////////////////////
// product_term:
//////////////////////////////////////////////////////////////////////////
template<class T>
void
product_term(const System<T> &ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res)
{
  CTimer aux_timer;
#ifdef WITH_PTHREADS
  DMTKglobals<T> *globals = NULL;
#else
  DMTKglobals<T> *globals = get_globals(T(0));
#endif

  aux_timer.Start();

  bool calc_hc = ss.h.use_hc();
// TODO WARNING: This may cause problems ?
#warning If you are using use_hc and you have problems, make sure you are not including the h.c. operator in the Hamiltonian 
  int ndiag = 0;
  if(t.is_diagonal()) calc_hc = false;
/*
  for(int i = 0; i < t.size(); i++){
    if(t[i].is_diagonal()) ndiag++;
  }
  if(ndiag == t.size()) calc_hc = false;
*/
///////////////////////////////////////////////
  if(ss.verbose() > 0)
    if(!calc_hc)
      cout << "APPLYING " << t.description() << endl;
    else
      cout << "APPLYING " << t.description() << " + h.c." << endl;


  if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

    aux_timer.Lap(); 

    if(ss.use_composite()){
//      if(auxt.ref_op[0]->is_diagonal() && auxt.sum_op.is_diagonal()) calc_hc = false;
      if(auxt.b1 < auxt.b2)
        product(*(auxt.ref_op[0]), auxt.sum_op, vv, res, auxt.b1, auxt.b2, T(1), calc_hc, globals, false);
      else
        product(auxt.sum_op, *(auxt.ref_op[0]), vv, res, auxt.b2, auxt.b1, T(1), calc_hc, globals, false);
    } else {
//      if(auxt.op[size_t(auxt.b1)-1]->is_diagonal() && auxt.op[size_t(auxt.b2)-1]->is_diagonal()) calc_hc = false;
      if(auxt.b1 < auxt.b2)
        product(*auxt.op[size_t(auxt.b1)-1], *auxt.op[size_t(auxt.b2)-1], vv, res, auxt.b1, auxt.b2, t.coef(), calc_hc, globals, false);
      else
        product(*auxt.op[size_t(auxt.b2)-1], *auxt.op[size_t(auxt.b1)-1], vv, res, auxt.b2, auxt.b1, t.coef(), calc_hc, globals, false);
    }

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  } 
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.Lap(); 

//    if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
    product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  } 
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
    aux_timer.Lap(); 
    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(1), calc_hc, globals);
    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
    aux_timer.Lap(); 

    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      if(auxt.b2 < auxt.b3){
        product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc, globals);
      } else {
        product(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc, globals);
      }

    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[1]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.Lap(); 

    if(ss.use_composite()){
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      if(auxt.b2 < auxt.b3){
        product(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc, globals);
      } else {
        product(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc, globals);
      }
    } else {
//      if(auxt.op[0]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[0]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  }
  else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    aux_timer.Lap(); 

    if(ss.use_composite()) {
//      if(auxt.sum_op.is_diagonal() && auxt.ref_op[0]->is_diagonal() && auxt.ref_op[1]->is_diagonal()) calc_hc = false;
      product(*(auxt.ref_op[0]), *(auxt.ref_op[1]), auxt.sum_op, vv, res, BLOCK2, BLOCK3, BLOCK4, T(1), calc_hc, globals);
    } else {
//      if(auxt.op[1]->is_diagonal() && auxt.op[2]->is_diagonal() && auxt.op[3]->is_diagonal()) calc_hc = false;
      product(*(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc, globals);
    }

    if(ss.verbose() > 0)
      cout << "PRODUCT TIME= " << aux_timer.LapTime().c_str() << endl;
  }
}

#ifndef USE_PRODUCT_DEFAULT
#ifndef USE_CUSTOM_PRODUCT
template<class T>
VectorState<T>
product(System<T>& ss, const VectorState<T>& vv)
{
  if(!ss.rotate_terms()) {
    cout << "WARNING: product not possible (rotate_terms == FALSE)\n";
    return vv;
  }
  if(ss.use_product_default()){
    return product_default(ss, vv);
  } else {
    if(ss.is_hermitian())
      return product_composite(ss, vv);
    else
      return product_composite(ss, vv, MASK_PRODUCT_DEFAULT);
  }
}

template<class T>
VectorState<T>
product_hc(System<T>& ss, const VectorState<T>& vv)
{
  if(!ss.rotate_terms()) {
    cout << "WARNING: product not possible (rotate_terms == FALSE)\n";
    return vv;
  }
  return product_composite(ss, vv, MASK_PRODUCT_HC);
}
#else
template<class T>
VectorState<T>
product(System<T>& ss, const VectorState<T>& vv)
{
  return product_custom(ss, vv);
}
#endif //USE_CUSTOM_PRODUCT
#endif //USE_PRODUCT_DEFAULT

//////////////////////////////////////////////////////////////////////////
// product_term:
//////////////////////////////////////////////////////////////////////////
template<class T>
void
init_term_composite(System<T> &ss, const AuxTerm<T> &auxt, const Term<T> &t, const VectorState<T> &vv, VectorState<T> &res)
{
  if(!ss.rotate_terms()) return;

  std::vector<ProductTerm<T> > new_terms;
  bool calc_hc = ss.h.use_hc();
// TODO WARNING: This may cause problems ?
#warning If you are using use_hc and you have problems, make sure you are not including the h.c. operator in the Hamiltonian 
  int ndiag = 0;
  if(t.is_diagonal()) calc_hc = false;
///////////////////////////////////////////////
  if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

    if(ss.use_composite()){
      if(auxt.b1 < auxt.b2){
        new_terms = get_product_terms(*(auxt.ref_op[0]), auxt.sum_op, vv, res, auxt.b1, auxt.b2, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), vv, res, auxt.b2, auxt.b1, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    } else {
      if(auxt.b1 < auxt.b2){
        new_terms = get_product_terms(*auxt.op[size_t(auxt.b1)-1], *auxt.op[size_t(auxt.b2)-1], vv, res, auxt.b1, auxt.b2, t.coef(), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*auxt.op[size_t(auxt.b2)-1], *auxt.op[size_t(auxt.b1)-1], vv, res, auxt.b2, auxt.b1, t.coef(), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    }
  } 
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){

    new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
    ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
  } 
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
    if(ss.use_composite()){
      new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(1), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[2]), vv, res, BLOCK1, BLOCK2, BLOCK3, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
    if(ss.use_composite()){
      if(auxt.b2 < auxt.b3){
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK2, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }

    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[1]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK2, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){

    if(ss.use_composite()){
      if(auxt.b2 < auxt.b3){
        new_terms = get_product_terms(auxt.sum_op, *(auxt.ref_op[0]), *(auxt.ref_op[1]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      } else {
        new_terms = get_product_terms(*(auxt.ref_op[1]), *(auxt.ref_op[0]), auxt.sum_op, vv, res, BLOCK1, BLOCK3, BLOCK4, T(1), calc_hc);
        ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
      }
    } else {
      new_terms = get_product_terms(*(auxt.op[0]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK1, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
  else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
    if(ss.use_composite()) {
      new_terms = get_product_terms(*(auxt.ref_op[0]), *(auxt.ref_op[1]), auxt.sum_op, vv, res, BLOCK2, BLOCK3, BLOCK4, T(1), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    } else {
      new_terms = get_product_terms(*(auxt.op[1]), *(auxt.op[2]), *(auxt.op[3]), vv, res, BLOCK2, BLOCK3, BLOCK4, T(t.coef()), calc_hc);
      ss.product_terms.insert(ss.product_terms.end(),new_terms.begin(),new_terms.end());
    }
  }
}

template <class T>
size_t 
System<T>::block(int site) const
{
  typename Lattice::const_iterator iter;
  valarray<const Block<T>* > b(4);
  valarray<int> offset(4);
  b[0] = _b1;
  b[1] = _b2;
  b[2] = _b3;
  b[3] = _b4;
  offset[0] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? h.lattice().size()/2 - _b1->lattice().size() - _b2->lattice().size() : 0;
  offset[1] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? offset[0] + _b1->lattice().size() : _b1->lattice().size();
  offset[2] = offset[1] + _b2->lattice().size();
  offset[3] = offset[2] + _b3->lattice().size();

  if(_in_warmup && _grow_symmetric && size() < h.lattice().size() && !_grow_outward){
    offset[3] = h.lattice().size() - _b4->lattice().size();
    offset[2] = offset[3] - _b3->lattice().size();
  }

  for(size_t i = 0; i < 4; i++){
    const Lattice &l = b[i]->lattice();
    int n = 0;
    for(iter = l.begin(); iter != l.end(); iter++){
        if((n++)+offset[i] == site){
           return (size_t)(i+1); 
        }
    }
  }
  return BLOCK_NONE;
}

template<class T>
const BasicOp<T>* 
System<T>::operator()(const BasicOp<T>& op) const
{
  typename Lattice::const_iterator iter;
  valarray<const Block<T>* > b(4);
  valarray<int> offset(4);
  b[0] = _b1;
  b[1] = _b2;
  b[2] = _b3;
  b[3] = _b4;
  offset[0] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? h.lattice().size()/2 - _b1->lattice().size() - _b2->lattice().size() : 0;
  offset[1] = _in_warmup && _grow_symmetric && _grow_outward && size() < h.lattice().size() ? offset[0] + _b1->lattice().size() : _b1->lattice().size();
  offset[2] = offset[1] + _b2->lattice().size();
  offset[3] = offset[2] + _b3->lattice().size();

  if(_in_warmup && _grow_symmetric && size() < h.lattice().size() && !_grow_outward){
    offset[3] = h.lattice().size() - _b4->lattice().size();
    offset[2] = offset[3] - _b3->lattice().size();
  }
  for(int i = 0; i < 4; i++){
    const Lattice &l = b[i]->lattice();
    int n = 0;
    for(iter = l.begin(); iter != l.end(); iter++){
        if((n++)+offset[i] == op.site()){
           BasicOp<T> _op(op);
           if(b[i]->single_site()) _op.set_site(op.site() - offset[i]);
// cout << op.name() << " " << op.site() << " " << op.internal_site() << " " << i << " " <<  _op.site() << endl;
           return (b[i]->operator()(_op));
        }
    }
  }
  return 0;
}



//////////////////////////////////////////////////////////////////

} // namespace dmtk

#endif // __DMTK_SYSTEM_H__
