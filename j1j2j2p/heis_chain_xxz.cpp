//#define WITH_WARNINGS
#define USE_VDISK  // Virtual disk: we store everything in memory
#define DO_MEASURE_DIMER // Measure correlation of dimers.
#define DO_MEASURE_SS // Measure correlation of Sz*Sz and S+*S-
#define DO_MEASURE_CHIRAL // Measure correlation of triangles.
#define DO_MEASURE_POLAR // Measures polar correlation.
#define MYTYPE double // Used for declaration of Hamiltonian types in (spliced in) Adrian's code.
#ifdef DO_MEASURE_CHIRAL
  #define MYTYPE complex<double>
#endif
#include <dmtk/dmtk.h>
#include "heis_site.cpp"

using namespace dmtk;

int main()
{
  int ntrunc;
  int niter;
  int dir = RIGHT2LEFT;
  int start = 1;
  int m;
  int custom;
  int nsteps;   
  int nx, xsz;
  Matrix<int> nstates(2,20);

// Input parameters for the problem
  cout << "Number of sites: " << endl;
  cin >> nx;
  cout << nx << endl;

  cout << "Sz: " << endl;
  cin >> xsz;
  cout << xsz << endl;

  double jz;
  cout << "Jz =" << endl;
  cin >> jz;
  cout << jz << endl; 

  double jx;
  cout << "Jx =" << endl;
  cin >> jx;
  cout << jx << endl; 

  double j2z;
  cout << "J2z =" << endl;
  cin >> j2z;
  cout << j2z << endl; 

  double j2x;
  cout << "J2x =" << endl;
  cin >> j2x;
  cout << j2x << endl;

  double j2pz;
  cout << "J2pz =" << endl;
  cin >> j2pz;
  cout << j2pz << endl; 

  double j2px;
  cout << "J2px =" << endl;
  cin >> j2px;
  cout << j2px << endl;
   

// Input parameters for the simulation
  cout << "Number of states:" << endl;
  cin >> ntrunc;

  cout << "Customize dmrg? (1->YES // 0->NO) :" << endl;
  cin >> custom;

  if(custom == 0){
    m = ntrunc/6;
    nstates(0,0) = 20; // warmup
    nstates(0,0) = 2*m; nstates(1,0) = 3*m;
    nstates(0,1) = 2*m; nstates(1,1) = 3*m;
    nstates(0,2) = 4*m; nstates(1,2) = 5*m;
    nstates(0,3) = 6*m; nstates(1,3) = 6*m;
    nstates(0,4) = 6*m; nstates(1,4) = 6*m;
    nsteps = 4;
  } else {
    cout << "Number of dmrg iterations:" << endl;
    cin >> nsteps;

    for(int i = 1; i <= nsteps; i++){
      cout << "Iteration " << i << " :" << endl;
      cout << "Number of states L2R:" << endl;
      cin >> nstates(0,i);
      cout << "Number of states R2L:" << endl;
      cin >> nstates(1,i);
    }
  }

  cout << "Start from:" << endl;
  cout << "    0 - Warmup loop" << endl;

  for(int i = 1; i <= nsteps; i++){
    cout << "    " << i << " - Iteration " << i << " : (" << nstates(0,i) << "," << nstates(1,i) << ") states" << endl;
  }

  cout << "  999 - Final Iteration" << endl;
  cin >> niter;

  if(niter >= 1 && niter <= nsteps || niter == 999){
    cout << "Direction (0->LEFT2RIGHT // 1->RIGHT2LEFT):\n";
    cin >> dir;
    cout << "Iteration:\n";
    cin >> start;
  }

// Specify quantum numbers for the problem
  QN::add_qn_index("Sz",false);
  QN::set_qn_mask(QN::default_mask());

// Initialize the system
  Hami<MYTYPE> hami = heis_chain<MYTYPE>(nx, OBC, jz, jx, j2z, j2x, j2pz, j2px); // TODO
  Lattice l = hami.lattice();
  System<MYTYPE> S(hami, hami.lattice(), "heis");
  S.set_use_single_site(false);
  S.set_use_hc(true);
  S.set_grow_symmetric(false);
  S.set_error(1.e-7);

// Quantum numbers for the ground state

  S.set_qn_mask(QN::default_mask());
  S.qnt["Sz"]=xsz;

// We calculate the gap
//  S.set_calc_gap(1);
//  S.set_verbose(1);

// Start simulation
  if(niter == 0){
    S.warmup_loop(20);
    niter = 1;
    cout << "**********************************\n";
  }

  if(niter != 999){
    S.set_lanczos_tolerance(1.e-9);
//    S.set_lanczos_maxiter(100);
    for(int i = niter; i <= nsteps; i++){
      cout << "NEW SWEEP " << i << " " << nstates(0,i) << " " << nstates(1,i) << endl;
      S.sweep(nstates(0,i),nstates(1,i),(size_t)dir,start);
      dir = RIGHT2LEFT; // or LEFT2RIGHT
      start = 1;
    }
  }

// Define correlations to measure
  S.set_store_products(false);
  for(int i = 0; i < l.size()-1; i++) {
    S.corr += Sz<MYTYPE>(i)*1;
    S.corr += Splus<MYTYPE>(i)*1;
  }


  /*
  for(int i = 0; i < l.size()-1; i++)
    S.corr += Sz<MYTYPE>(i)*Sz<MYTYPE>(i+1);

  
  for(int i = 0; i < l.size()-1; i++)
    for(int j = i+1; j < l.size(); j++)
      S.corr += Sz<double>(i)*Sz<double>(j);

  for(int i = 0; i < l.size()-1; i++)
    for(int j = i+1; j < l.size(); j++)
      S.corr += Splus<double>(i)*Sminus<double>(j);
  */

// ======================= BEG ADRIAN STUFF ====================
int jini = l.size()/2;
int jend = jini+1;
                                             
#ifdef DO_MEASURE_DIMER
  for(int j = jini; j < jend; j++){ // Inner-most dimers. Some range.
    Hami<MYTYPE> dimer0;
    dimer0.clear();
    dimer0 += Splus<MYTYPE>(j)*Sminus<MYTYPE>(j+1)*MYTYPE(0.5);
    dimer0 += Sminus<MYTYPE>(j)*Splus<MYTYPE>(j+1)*MYTYPE(0.5);
    dimer0 += Sz<MYTYPE>(j)*Sz<MYTYPE>(j+1)*MYTYPE(1);

    for(int k = 0; k < l.size()-1; k++){ // All dimers.
      Hami<MYTYPE> dimer1;
      dimer1.clear();
      dimer1 += Splus<MYTYPE>(k)*Sminus<MYTYPE>(k+1)*MYTYPE(0.5);
      dimer1 += Sminus<MYTYPE>(k)*Splus<MYTYPE>(k+1)*MYTYPE(0.5);
      dimer1 += Sz<MYTYPE>(k)*Sz<MYTYPE>(k+1)*MYTYPE(1);

      if(abs(k-j) > 1){
        char text[200];
        sprintf(text,"Hdimer(%d,%d)",j,k);
        Hami<MYTYPE> hdd;
        //hdd = hdimer0*hdimer1;
        hdd = dimer0*dimer1; // From line above - typo?
        hdd.set_name(text);
        S.corr += BasicOp<MYTYPE>(hdd);
      }
    }
  }
#endif

#ifdef DO_MEASURE_SS
//////////////////////////////// S_i.S_j
  for(int i = 0; i < l.size()-1; i++) {
    for(int j = i+1; j < l.size(); j++) {
      S.corr += Sz<MYTYPE>(i)*Sz<MYTYPE>(j);
      S.corr += Splus<MYTYPE>(i)*Sminus<MYTYPE>(j);
    }
  }

  /* // Original
  for(int i = 0; i < lx; i++){
    for(int j = 0; j < lx; j++){
      if(i != j) S.corr += Sz<MYTYPE>(i,0)*Sz<MYTYPE>(j,0);
      if(i != j) S.corr += Splus<MYTYPE>(i,0)*Sminus<MYTYPE>(j,0);
    }
  }
  */
#endif

#ifdef DO_MEASURE_CHIRAL
  //int jini = l.size()/2-1;
  //int jend = jini+1;

// TODO Question: Do we want only one center triangle
// (which will have J1 and either J2 or J2P)
// or do we want two center triangles?

//////////////////////////////// SCALAR CHIRALITY
  for(int j = jini; j < jend; j++){
    // Find the middle triangle.
    Hami<MYTYPE> hc;
    hc += Splus<MYTYPE>(j)*Sminus<MYTYPE>(j+1)*Sz<MYTYPE>(j+2)*complex<double>(0,0.5);
    hc += Sminus<MYTYPE>(j)*Splus<MYTYPE>(j+1)*Sz<MYTYPE>(j+2)*complex<double>(0,-0.5);
   
    hc += Sminus<MYTYPE>(j)*Sz<MYTYPE>(j+1)*Splus<MYTYPE>(j+2)*complex<double>(0,0.5);
    hc += Splus<MYTYPE>(j)*Sz<MYTYPE>(j+1)*Sminus<MYTYPE>(j+2)*complex<double>(0,-0.5);

    hc += Sz<MYTYPE>(j)*Splus<MYTYPE>(j+1)*Sminus<MYTYPE>(j+2)*complex<double>(0,0.5);
    hc += Sz<MYTYPE>(j)*Sminus<MYTYPE>(j+1)*Splus<MYTYPE>(j+2)*complex<double>(0,-0.5);
    char text[200];
    sprintf(text,"Hscalarlocal(%d)",j);
    hc.set_name(text);
    S.corr += BasicOp<MYTYPE>(hc);

    for(int m = 0; m < l.size()-2; m++){
      // Find all other triangles.
      Hami<MYTYPE> hc2;
      hc2 += Splus<MYTYPE>(m)*Sminus<MYTYPE>(m+1)*Sz<MYTYPE>(m+2)*complex<double>(0,0.5);
      hc2 += Sminus<MYTYPE>(m)*Splus<MYTYPE>(m+1)*Sz<MYTYPE>(m+2)*complex<double>(0,-0.5);
     
      hc2 += Sminus<MYTYPE>(m)*Sz<MYTYPE>(m+1)*Splus<MYTYPE>(m+2)*complex<double>(0,0.5);
      hc2 += Splus<MYTYPE>(m)*Sz<MYTYPE>(m+1)*Sminus<MYTYPE>(m+2)*complex<double>(0,-0.5);
  
      hc2 += Sz<MYTYPE>(m)*Splus<MYTYPE>(m+1)*Sminus<MYTYPE>(m+2)*complex<double>(0,0.5);
      hc2 += Sz<MYTYPE>(m)*Sminus<MYTYPE>(m+1)*Splus<MYTYPE>(m+2)*complex<double>(0,-0.5);
 
      if(abs(j-m) > 2) { // Do not include the middle-middle comparison. 
        char text[200];
        Hami<MYTYPE> hcc;
        hcc = hc*hc2;
        sprintf(text,"Hscalar_scalar(%d,%d)",j,m);
        hcc.set_name(text);
        S.corr += BasicOp<MYTYPE>(hcc);
      }
    }
  }

//////////////////////////////// VECTOR CHIRALITY

  for(int j = jini; j < jend; j++){
    // Find middle triangle.
    Hami<MYTYPE> hc;
    hc += Splus<MYTYPE>(j)*Sminus<MYTYPE>(j+1)*complex<double>(0,0.5);
    hc += Sminus<MYTYPE>(j)*Splus<MYTYPE>(j+1)*complex<double>(0,-0.5);

    char text[200];
    sprintf(text,"Hvectorlocal(%d)",j);
    hc.set_name(text);
    S.corr += BasicOp<MYTYPE>(hc);
   
    for(int m = 0; m < l.size()-1; m++){
      // Find all other triangles.
      Hami<MYTYPE> hc2;
      hc2 += Splus<MYTYPE>(m)*Sminus<MYTYPE>(m+1)*complex<double>(0,0.5);
      hc2 += Sminus<MYTYPE>(m)*Splus<MYTYPE>(m+1)*complex<double>(0,-0.5);
     
      if(abs(j-m) > 1) { // Do not include the middle-middle comparison.
        char text[200];
        Hami<MYTYPE> hcc;
        hcc = hc*hc2;
        sprintf(text,"Hvector_vector(%d,%d)",j,m);
        hcc.set_name(text);
        S.corr += BasicOp<MYTYPE>(hcc);
      }
    }
  }
#endif //DO_MEASURE_CHIRAL



#ifdef DO_MEASURE_POLAR
////////////////////////////////// MULTIPOLAR
  for(int j = jini; j < jend; j++){
    Hami<MYTYPE> hc;
    hc += Sminus<MYTYPE>(j)*Sminus<MYTYPE>(j+1);
   
    for(int m = j+2; m < l.size()-1; m++){
      Hami<MYTYPE> hc2;
      hc2 += Splus<MYTYPE>(m)*Splus<MYTYPE>(m+1)*complex<double>(0,0.5);
     
      if(abs(j-m) > 2) { 
        char text[200];
        Hami<MYTYPE> hcc;
        hcc = hc*hc2;
        sprintf(text,"Hpolar(%d,%d)",j,m);
        hcc.set_name(text);
        S.corr += BasicOp<MYTYPE>(hcc);
      }
    }
  }
#endif

/*
  S.set_store_products(false);
  S.final_sweep(ntrunc,ntrunc,(size_t)dir,start,false);

  for(iter = S.corr.begin(); iter != S.corr.end(); iter++){
    Term<MYTYPE> &t = (*iter);
    cout << "CORR " << t.name() << " " << t.value() << endl;
  }
}
*/


// ======================= END ADRIAN STUFF ====================


//  for(int i = 0; i < l.size()-3; i++)
//  for(int i = 0; i < l.size()-3; i++)
//    S.corr += Sz<double>(i)*Sz<double>(i+1)*Sz<double>(i+2)*Sz<double>(i+3);

// Final sweep
  S.final_sweep(ntrunc);
// Measure
//  S.measure();


  Hami<MYTYPE>::iterator corr_iter;
  for(corr_iter = S.corr.begin(); corr_iter != S.corr.end(); corr_iter++){
    Term<MYTYPE> &t = (*corr_iter);
    if(t.size() == 1)
      cout << t.description() << " " << t[0].site() << " " << t.value() << endl;
    else if(t.size() == 2) 
      cout << t.description() << " " << t[0].site() << " " << t[1].size() << " " << t.value() << endl;
    else 
      cout << t.description() << t.value() << endl;
  }


// Output
  if(S.calc_gap())
    cout << "RESULTS " << setprecision(12) << S.energy[0] << " " << S.energy[1] << " " << S.energy[1]-S.energy[0] << endl;
  else
    cout << "RESULTS " <<  S.energy[0] << endl;
}
                                                                     
                                                                     
                                                                     
