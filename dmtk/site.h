#ifndef __DMTK_SITE_OPS_H__
#define __DMTK_SITE_OPS_H__

namespace dmtk
{

template<class T>
Block<T>
build_complex_site(const Hami<T> & h, Hami<T> corr = Hami<T>(), bool add_hami = true, bool add_operators = true)
{
  Lattice l= h.lattice();
  Block<T> bl = h.get_site(0);

  typename Hami<T>::const_iterator hiter;
  typename Block<T>::iterator biter;
  typename Block<T>::const_iterator cbiter;
  for(int is = 0; is < l.size()-1; is++){
    cout << "BUILDING COMPLEX SITE: ADDING SITE " << is << endl;
    const Block<T> &br = h.get_site(is+1);
    Basis aux_basis(bl.basis(),br.basis());
    aux_basis.reorder();

    Block<T> aux_b;
    aux_b.clear();
    aux_b.resize(aux_basis);
    aux_b.set_lattice(Lattice(1,OBC));
    aux_b.set_orbitals(is+1);

    BMatrix<T> aux_rho(aux_basis);
    PackedBasis::const_iterator pbiter;

////////////////////////////////////////////////////////////////////
// We need the identity matrix -- ignore
////////////////////////////////////////////////////////////////////
    for(pbiter = aux_basis.subspaces().begin(); pbiter != aux_basis.subspaces().end(); pbiter++){

      SubMatrix<T> block(pbiter->qn(),*pbiter,*pbiter);
      block = I<T>();
      block.qn() = (*pbiter).qn();
      aux_rho.push_back(block);
    }

    if(add_operators){
////////////////////////////////////////////////////////////////////
//  Transform and add operators from the left block
////////////////////////////////////////////////////////////////////
      for(biter = bl.begin(); biter != bl.end(); biter++){
        BasicOp<T> op = *biter;
        if(op.name() != h.name()){
          BasicOp<T> new_op = op.internals();
          new_op.resize(aux_basis);
          cout << "ADDING NEW OPERATOR LEFT " << new_op.description() << endl;
          new_operator(new_op, op, aux_rho, aux_basis, LEFT, T(1), false);
          aux_b.push_back(new_op);
        }
      }
////////////////////////////////////////////////////////////////////
//  Transform and add operators from the right block
////////////////////////////////////////////////////////////////////
      for(cbiter = br.begin(); cbiter != br.end(); cbiter++){
        const BasicOp<T> &op = *cbiter;
        if(op.name() != h.name()){
          BasicOp<T> new_op = op.internals();
          new_op.resize(aux_basis);
          new_op.set_site(0).set_internal_site(is+1);
          cout << "ADDING NEW OPERATOR RIGHT " << new_op.description() << endl;
          new_operator(new_op, op, aux_rho, aux_basis, RIGHT, T(1), false);
          aux_b.push_back(new_op);
        }
      }
    }
////////////////////////////////////////////////////////////////////
//  Transform and add the new Hamiltonian
////////////////////////////////////////////////////////////////////
    H<T> hij(aux_basis,0);
    hij.set_name(h.name().c_str());

    const BasicOp<T> *hop = bl(H<T>(0).set_name(h.name().c_str()));
    if(hop && (add_hami || is > 0)){
      const Block<T> &b = h.get_site(0);
      cout << "ADDING NEW HAMILTONIAN TERM HL " << h.name() << endl;
      new_operator(hij, *hop, aux_rho, aux_basis, LEFT, T(1), false);
    }
    hop = br(H<T>(0).set_name(h.name().c_str()));
    if(hop && add_hami){
      const Block<T> &b = h.get_site(0);
      cout << "ADDING NEW HAMILTONIAN TERM HR " << h.name() << endl;;
      new_operator(hij, *hop, aux_rho, aux_basis, RIGHT, T(1), false);
    }

    for(hiter = h.begin(); hiter != h.end(); hiter++){
      const Term<T> &t = *hiter;
      if(t.size() == 1){
        if(t[0].site() == is+1 || (t[0].site() == 0 && is == 0)){
          const Block<T> &b = h.get_site(t[0].site());
          BasicOp<T> ref_op = t[0].internals();
          const BasicOp<T> *op = b(ref_op.set_site(0));
          if(op){
            int pos = (t[0].site() == 0) ? LEFT : RIGHT;
            cout << "ADDING NEW HAMILTONIAN TERM " << hij.name() << " " << t.description() << " " << pos << endl;
            new_operator(hij, *op, aux_rho, aux_basis, pos, t.coef(), false);
          }
        }
      } 
      if(t.size() == 2){
        if(t[0].site() == is+1 && t[1].site() <= is){
          BasicOp<T> ref_op0 = t[0].internals();
          BasicOp<T> ref_op1 = t[1].internals();
          const BasicOp<T> *op0 = br(ref_op0.set_site(0));
          ref_op1.set_site(0);
          if(is > 0) ref_op1.set_internal_site(t[1].site());
          BasicOp<T> *op1 = bl(ref_op1);
          if(op0 && op1){
            bool calc_hc = true;
            if(t.is_diagonal() || !h.use_hc()){
              calc_hc = false;
              cout << "ADDING NEW HAMILTONIAN TERM " << hij.name() << " " << t.description() << " " << op1->description() << " " << op0->description() << endl;
            } else {
              cout << "ADDING NEW HAMILTONIAN TERM " << hij.name() << " " << t.description() << " + h.c. " << op1->description() << " " << op0->description() << endl;
            }
            new_operator(hij, *op1, *op0, BLOCK1, BLOCK2, aux_rho, aux_basis, t.coef(), calc_hc, false);
          }
        }
      }

    }
    aux_b.push_back(hij);

/*
    typename BMatrix<T>::const_iterator biter;
    for(biter = hij.begin(); biter != hij.end(); biter++){
      const SubMatrix<T> &block = (*biter);
      for(int l = 0; l < block.rows(); l++)
       for(int j = 0; j < block.rows(); j++) cout << l << " " << j << " " << block(l,j) << endl;
    }
*/
////////////////////////////////////////////////////////////////////
//  Transform and add the new operators for correlations 
////////////////////////////////////////////////////////////////////
    for(hiter = corr.begin(); hiter != corr.end(); hiter++){
      const Term<T> &t = *hiter;
      if(t.size() == 2) {
        if(t[0].site() == is+1 && t[1].site() <= is){
          BasicOp<T> new_op;
          BasicOp<T> aux_op(t);
          new_op.set_fermion(aux_op.fermion());
          new_op.dqn = t.dqn();
          new_op.set_name(t.description().c_str());
          new_op.resize(aux_basis);
          BasicOp<T> ref_op0 = t[0].internals();
          BasicOp<T> ref_op1 = t[1].internals();
          const BasicOp<T> *op0 = br(ref_op0.set_site(0));
          ref_op1.set_site(0);
          if(is > 0) ref_op1.set_internal_site(t[1].site());
          BasicOp<T> *op1 = bl(ref_op1);
          if(op0 && op1){
            bool calc_hc = true;
            if(t.is_diagonal() || !corr.use_hc()){
              calc_hc = false;
              cout << "ADDING NEW CORRELATION " << t.description() << " " << op1->description() << " " << op0->description() << endl;
            } else {
              cout << "ADDING NEW CORRELATION " << t.description() << " + h.c. " << op1->description() << " " << op0->description() << endl;
            }
            new_operator(new_op, *op1, *op0, BLOCK1, BLOCK2, aux_rho, aux_basis, t.coef(), calc_hc, false);
            aux_b.push_back(new_op);
          }
        }
      }
      if(t.size() > 2){
        Term<T> aux_term;
        BasicOp<T> ref_op1;
        bool found1 = false;
        bool found2 = false;
        for(int i = 0; i < t.size(); i++){
          BasicOp<T> ref_op = t[i].internals();
          if(ref_op.site() <= is){
            aux_term *= ref_op;
            found1 = true;
          }
          if(ref_op.site() == is+1){
            ref_op1 = ref_op.internals();
            found2 = true;
          }
        }

        if(found1 && found2){
          Term<T> new_term;
          for(int i = 0; i < t.size(); i++){
            BasicOp<T> ref_op = t[i].internals();
            if(ref_op.site() <= is+1) new_term *= ref_op;
          }
          BasicOp<T> ref_op0;
          if(aux_term.size() == 1) 
            ref_op0 = aux_term[0].internals();
          else 
            ref_op0.set_name(aux_term.name().c_str());
          const BasicOp<T> *op0 = bl(ref_op0.set_site(0));
          const BasicOp<T> *op1 = br(ref_op1.set_site(0));
          if(op0 && op1){
            BasicOp<T> new_op(new_term);
            new_op.set_site(0);
            new_op.dqn = new_term.dqn();
            new_op.resize(aux_basis);
            bool add = aux_b.contains(new_op);

            if(!add){
              bool calc_hc = true;
              if(t.is_diagonal() || !corr.use_hc()){
                calc_hc = false;
                cout << "ADDING NEW CORRELATION TERM " << t.description() << " " << new_op.name() << " " << op1->description() << " " << op0->description() << endl;
              } else {
                cout << "ADDING NEW CORRELATION TERM " << t.description() << " + h.c. " << new_op.name() << " " << op1->description() << " " << op0->description() << endl;
              }
              T coef(1);
              if(is == l.size()-2) {
                new_op.set_name(t.description().c_str());
                coef = t.coef();
              }
              new_operator(new_op, *op0, *op1, BLOCK1, BLOCK2, aux_rho, aux_basis, coef, calc_hc, false);
              aux_b.push_back(new_op);
// for(int i = 0; i < aux_basis.size(); i++)
// for(int j = 0; j < aux_basis.size(); j++) cout << i << " " << j << " " << new_op(i,j) << endl;
            }
          }
        }
      }


    }
////////////////////////////////////////////////////////////////////
//  Transform and add the new operators for correlations 
////////////////////////////////////////////////////////////////////
    bl = aux_b;
  }

/*
  for(biter = bl.begin(); biter != bl.end(); biter++) {
//    if((*biter).name() == "H") (*biter).set_name(h.name().c_str());
    cout << (*biter).description() << endl;

    typename BMatrix<T>::const_iterator oiter;
    for(oiter = (*biter).begin(); oiter != (*biter).end(); oiter++){
      const SubMatrix<T> &block = (*oiter);
      for(int l = 0; l < block.rows(); l++)
       for(int j = 0; j < block.rows(); j++) cout << l << " " << j << " " << block(l,j) << endl;
    }

  }
*/
  
  bl.set_orbitals(l.size());
  return bl; 
}

} // namespace dmtk

#endif // __DMTK_SITE_OPS_H__
