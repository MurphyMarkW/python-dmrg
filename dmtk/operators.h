#ifndef __DMTK_OPERATORS_H__
#define __DMTK_OPERATORS_H__

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <string>
#include <iosfwd>
#include "enums.h"
#include "conj.h"
#include "vector.h"
#include "qn.h"
#include "subspace.h"
#include "block_matrix.h"
#include "state.h"
#include "globals.h"

namespace dmtk
{

// Basic operators
template <class T> class Block;
template <class T> class Term;
template <class T> class Hami;

enum
{  
  OP_SYSTEM,
  OP_MEASURE,
  OP_ADDITIVE,
};

template <class T>
class BasicOp: public dmtk::BMatrix<T>
{
  protected:
    int _site;     // site where the operator is seated on
    int _label; // if we use more than one site per unit cell;
    int _internal_site;    // any arbitrary internal_site we want to use
                   // in general: relative position for a product operator
    bool _is_fermion;
    bool _is_diagonal;

    size_t _type;

    std::string _name;

    virtual void init() {};
    void init_blocks();
    void bmatrix_from_matrix(const Matrix<T> &);

    Term<T> _term; // if this is a composite operator originated in a product
    bool _is_term;
    Hami<T> _hami; // if this is a sum of terms
    bool _is_hami;

  public:
    typedef BMatrix<T> _BM;
    typedef std::vector<SubMatrix<T> > _BMvector;

    QN dqn; // (Sz(after) - Sz(before))
    	    // (Nt(after) - Nt(before))

    BasicOp(): _internal_site(0), _site(0), _label(-1), _type(OP_SYSTEM), _is_term(false), _is_hami(false), _is_fermion(false), _is_diagonal(true) {}

    BasicOp(const char *name, int site, int internal_site = 0, int label = -1): _name(name), _site(site), _internal_site(internal_site), _label(label), _type(OP_SYSTEM), _is_term(false), _is_hami(false), _is_fermion(false), _is_diagonal(true) {}
    BasicOp(const char *name, const Basis& b, int site, int internal_site = 0, int label = -1): 
      BMatrix<T>(b), _name(name), _internal_site(internal_site), _site(site), _label(label), _type(OP_SYSTEM), _is_term(false), _is_hami(false), _is_fermion(false), _is_diagonal(true) {}
    BasicOp(const char *name, const Matrix<T> &_m, const Basis& b, int site, int internal_site = 0, int label = -1):
       _name(name), BMatrix<T>(b), _site(site), _internal_site(internal_site), _label(label), _type(OP_SYSTEM), _is_term(false), _is_hami(false), _is_fermion(false) {}

    BasicOp(const BasicOp<T> &o): BMatrix<T>(o), _internal_site(o._internal_site), dqn(o.dqn), _name(o._name), _type(o._type), _site(o._site), _label(o._label), _is_term(o._is_term), _is_hami(o._is_hami), _is_fermion(o._is_fermion), _is_diagonal(o._is_diagonal), _term(o._term), _hami(o._hami) {} 


    BasicOp(const Term<T>& t): _type(OP_SYSTEM)
      {
        if(t.size() == 0 && t.coef() == T(1)) {
          this->operator=(t[0]);
          return;
        }
        _name = string(t.name());
        _site = t[0].site();
        _internal_site = t[0].internal_site();
        _term = t;
        _is_term = true;
        _is_fermion = false;
        _is_diagonal = t.is_diagonal();
        typename vector<BasicOp<T> >::const_iterator iter;
        int nfermions = 0;
        for(iter = t.begin(); iter != t.end(); iter++){
          dqn += (*iter).dqn;
          if(iter->fermion()) nfermions++;
        }
        if(SGN(nfermions) == -1) _is_fermion = true;
      }

    BasicOp(const Hami<T> &h): _type(OP_ADDITIVE)
      {
        _name = string(h.name());
        _site = 0;
        _internal_site = 0;
        _hami = h;
        _is_hami = true;
        _is_fermion = false;
        _is_term = false;
        if(h.size() > 0){
          BasicOp<T> first_op = *(h.begin());
          dqn = first_op.dqn;
          _is_fermion = first_op.fermion();
        }
        _is_diagonal = true;
      }

    BasicOp(const BasicOp<T> &op1, const BasicOp<T> &op2);

//  Asignment

    BasicOp& operator=(const T& v)
      {BMatrix<T>::operator=(v); return *this;}
    BasicOp& operator=(const BMatrix<T> &m)
      {BMatrix<T>::operator=(m); return *this;}
    BasicOp& operator=(const BasicOp<T> &m)
      {
        BMatrix<T>::operator=(m); 
        dqn = m.dqn; 
        _type = m._type;
        _site = m._site;
        _internal_site = m._internal_site;
        _label = m._label;
        _name = m._name;
        _is_term = m._is_term;
        _is_fermion = m._is_fermion;
        _term = m._term;
        _is_hami = m._is_hami;
        _hami = m._hami;
        _is_diagonal = m._is_diagonal;
        return *this;
      }

    BasicOp& operator=(const Term<T>& t)
      {
        _type = OP_SYSTEM;
        _name = string(t.name());
        _site = t[0].site();
        _internal_site = t[0].internal_site();
        _term = t;
        _is_term = true;
        _is_fermion = false;
        _is_diagonal = t.is_diagonal();
        typename vector<BasicOp<T> >::const_iterator iter;
        int nfermions = 0;
        for(iter = t.begin(); iter != t.end(); iter++){
          dqn += (*iter).dqn;
          if(iter->fermion()) nfermions++;
        }
        if(SGN(nfermions) == -1) _is_fermion = true;
        return *this;
      }
/*
    BasicOp& operator=(const Hami<T> &h)
      {
        _type = OP_ADDITIVE;
        _name = string(h.name());
        _site = 0;
        _internal_site = 0;
        _hami = h;
        _is_hami = true;
        _is_fermion = false;
        _is_term = false;
        BasicOp<T> first_op = *(h.begin());
        dqn = first_op.dqn;
        _is_fermion = first_op._is_fermion;
        return *this;
      }
*/

    BasicOp& operator=(const Matrix<T> &m)
      { bmatrix_from_matrix(m); return *this; }
    BasicOp& set_matrix(const Matrix<T> &m, const Basis &b, QN _dqn)
      { dqn = _dqn; _BM::repack(b); bmatrix_from_matrix(m); return *this; }

    BasicOp internals() const
      {
        BasicOp<T> ret_op;

        ret_op.dqn = dqn; 
        ret_op._type = _type;
        ret_op._site = _site;
        ret_op._internal_site = _internal_site;
        ret_op._label = _label;
        ret_op._name = _name;
        ret_op._is_term = _is_term;
        ret_op._is_fermion = _is_fermion;
        ret_op._term = _term;
        ret_op._is_hami = _is_hami;
        ret_op._hami = _hami;
        ret_op._is_diagonal = _is_diagonal;

        return ret_op;
      }

    Term<T> operator*(T v) const 
      { return Term<T>(*this,v); }
    Term<T> operator*(BasicOp<T> other) const 
      { 
        if(this->is_term()) return (this->term()*other); 
        if(other.is_term()) return ((*this)*other.term()); 
        return Term<T>(*this,other);
      }
    Term<T> operator*(const Term<T>& t) const
      { 
        Term<T> aux(t); 
        aux.push_back(*this); 
        return aux; 
      }

    bool operator==(const BasicOp<T>& op) const
      { return (op._name == _name && op._site == _site && op._internal_site == _internal_site && op._label == _label); } 
    bool operator==(const Term<T>& t) const
      { return (name() == t.name()); }

//  Methods

    BasicOp& resize(const Basis &b) 
      { BMatrix<T>::repack(b); _BMvector::clear(); init_blocks(); return *this; }
    BasicOp& resize(const PackedBasis &b) 
      { BMatrix<T>::repack(b); _BMvector::clear(); init_blocks(); return *this; }


    BasicOp reshape(int qn_mask) const;

    bool fermion() const 
      { return (dqn.fermion_sign() == -1); }
    bool apply_sign() const { return fermion(); }
    bool is_qn_diagonal() const 
      { return (dqn.equal(QN(),QN::default_mask())); }
    bool is_diagonal() const { return _is_diagonal; }

    std::string const& name() const { return _name; }
    int site() const { return _site; }
    int internal_site() const { return _internal_site; }
    int orbital() const { return _internal_site; }
    int label() const { return _label; }
    size_t type() const { return _type; }
    bool is_term() const { return _is_term; }
    bool is_hami() const { return _is_hami; }

    const Hami<T> &hami() const { return _hami; }
    Hami<T> &hami() { return _hami; }
    const Term<T> &term() const { return _term; }
    Term<T> &term() { return _term; }

    BasicOp& set_type(size_t type) { _type = type; return *this; }
    BasicOp& set_internal_site(int internal_site) { _internal_site = internal_site; return *this; }
    BasicOp& set_fermion(bool fermion) { _is_fermion = fermion; return *this; }
    BasicOp& set_diagonal(bool diagonal) { _is_diagonal = diagonal; return *this; }
    BasicOp& set_site(int site) { _site = site; return *this; }
    BasicOp& set_name(const char* name) { _name = std::string(name); return *this; }
    BasicOp& set_term(const Term<T> *t) { _is_term= t == NULL ? false : true; if(t) _term = *t; return *this; }
    BasicOp& set_hami(const Hami<T> *h) { _is_hami = h == NULL ? false : true; if(h) _hami = *h; return *this; }
    QN &get_dqn() { return dqn; }
    QN get_dqn() const { return dqn; }

    friend class Block<T>;
    friend class Hami<T>;
    friend class Term<T>;

    string description() const 
      {
        ostringstream bf;
        if(_is_term)
          return _term.description();
        else if(_is_hami)
          bf << "(" << _hami.description() << ")";
        else
          bf << name() << "(" << site() << "," << internal_site() << ")";
           
        return bf.str();
      }


    // Streams

    void write(std::ostream &s) const
    {
      s.write((const char *)&_site, sizeof(int));
      s.write((const char *)&_internal_site, sizeof(int));
      s.write((const char *)&_label, sizeof(int));
      s.write((const char *)&_type, sizeof(size_t));
      s.write((const char *)&_is_term, sizeof(bool));
      s.write((const char *)&_is_fermion, sizeof(bool));
      s.write((const char *)&_is_diagonal, sizeof(bool));
      if(_is_term) _term.write(s);
      dqn.write(s);
      size_t n = _name.size()+1;
      char *t = new char[n];
      snprintf(t, n, _name.c_str());
      size_t l = strlen(t);
      s.write((const char *)&l, sizeof(size_t));
      s.write((const char *)t, l*sizeof(char));
      BMatrix<T>::write(s);
      delete(t);
    }

    void read(std::istream &s)
    {
      s.read((char *)&_site, sizeof(int));
      s.read((char *)&_internal_site, sizeof(int));
      s.read((char *)&_label, sizeof(int));
      s.read((char *)&_type, sizeof(size_t));
      s.read((char *)&_is_term, sizeof(bool));
      s.read((char *)&_is_fermion, sizeof(bool));
      s.read((char *)&_is_diagonal, sizeof(bool));
      if(_is_term) _term.read(s);
      dqn.read(s);
      size_t l;
      s.read((char *)&l, sizeof(size_t));
      int l1 = l+1;
      char *t = new char[l1];
      s.read((char *)t, l*sizeof(char));
      t[l] = '\0';
      _name = string(t);
      BMatrix<T>::read(s);
      delete(t);
    }

};
template <class T>
class NullOp: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    NullOp(): BasicOp<T>("NULL",0) 
      { init(); }
    NullOp(int site, int internal_site = 0): BasicOp<T>("NULL",site,internal_site)  
      { init(); }
    NullOp(const Basis& b, int site, int internal_site = 0): BasicOp<T>("NULL",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    NullOp(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("NULL",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
 
template <class T>
class Identity: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    Identity(): BasicOp<T>("I",0) 
      { init(); }
    Identity(int site, int internal_site = 0): BasicOp<T>("I",site,internal_site)  
      { init(); }
    Identity(const Basis& b, int site, int internal_site = 0): BasicOp<T>("I",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Identity(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("I",m,b,site,internal_site) 
      { init(); bmatrix_from_matrix(m); }
};
    
template <class T>
class H: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init() 
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = 0; _O::_type = OP_ADDITIVE; }
  public:
    H(): BasicOp<T>("H",0) 
      { init(); }
    H(int site, int internal_site = 0): BasicOp<T>("H",site,internal_site)  
      { init(); }
    H(const Basis& b, int site, int internal_site = 0): BasicOp<T>("H",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    H(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("H",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Sz: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    Sz(): BasicOp<T>("Sz",0) 
      { init(); }
    Sz(int site, int internal_site = 0): BasicOp<T>("Sz",site,internal_site)  
      { init(); }
    Sz(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Sz",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Sz(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Sz",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class N: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    N(): BasicOp<T>("N",0) 
      { init(); }
    N(int site, int internal_site = 0): BasicOp<T>("N",site,internal_site)  
      { init(); }
    N(const Basis& b, int site, int internal_site = 0): BasicOp<T>("N",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    N(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("N",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Splus: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn["Sz"] = 2; _O::_is_diagonal = false; }
  public:
    Splus(): BasicOp<T>("S+",0) 
      { init(); }
    Splus(int site, int internal_site = 0): BasicOp<T>("S+",site,internal_site)  
      { init(); }
    Splus(const Basis& b, int site, int internal_site = 0): BasicOp<T>("S+",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Splus(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("S+",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Sminus: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = false; _O::dqn["Sz"] = -2; _O::_is_diagonal = false; }
  public:
    Sminus(): BasicOp<T>("S-",0) 
      { init(); }
    Sminus(int site, int internal_site = 0): BasicOp<T>("S-",site,internal_site)  
      { init(); }
    Sminus(const Basis& b, int site, int internal_site = 0): BasicOp<T>("S-",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Sminus(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("S-",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};

template <class T>
class C: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.n() = -1; _O::_is_diagonal = false; }
  public:
    C(): BasicOp<T>("C",0) 
      { init(); }
    C(int site, int internal_site = 0): BasicOp<T>("C",site,internal_site)  
      { init(); }
    C(const Basis& b, int site, int internal_site = 0): BasicOp<T>("C",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    C(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("C",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Cd: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()  
      { _O::_is_fermion = true; _O::dqn.n() = 1; _O::_is_diagonal = false; }
  public:
    Cd(): BasicOp<T>("Cd",0) 
      { init(); }
    Cd(int site, int internal_site = 0): BasicOp<T>("Cd",site,internal_site)  
      { init(); }
    Cd(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cd",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cd(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cd",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    

template <class T>
class Cup: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = -1; _O::dqn.n() = -1; _O::_is_diagonal = false; }
  public:
    Cup(): BasicOp<T>("Cup",0) 
      { init(); }
    Cup(int site, int internal_site = 0): BasicOp<T>("Cup",site,internal_site)  
      { init(); }
    Cup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cup",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cup",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Cdup: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()  
      { _O::_is_fermion = true; _O::dqn.sz() = 1; _O::dqn.n() = 1; _O::_is_diagonal = false; }
  public:
    Cdup(): BasicOp<T>("Cdup",0) 
      { init(); }
    Cdup(int site, int internal_site = 0): BasicOp<T>("Cdup",site,internal_site)  
      { init(); }
    Cdup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdup",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cdup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdup",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Cdn: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = 1; _O::dqn.n() = -1; _O::_is_diagonal = false; }
  public:
    Cdn(): BasicOp<T>("Cdn",0) 
      { init(); }
    Cdn(int site, int internal_site = 0): BasicOp<T>("Cdn",site,internal_site)  
      { init(); }
    Cdn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdn",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cdn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cdn",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Cddn: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      { _O::_is_fermion = true; _O::dqn.sz() = -1; _O::dqn.n() = 1; _O::_is_diagonal = false; }
  public:
    Cddn(): BasicOp<T>("Cddn",0) 
      { init(); }
    Cddn(int site, int internal_site = 0): BasicOp<T>("Cddn",site,internal_site)  
      { init(); }
    Cddn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cddn",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Cddn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Cddn",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
// Boson creation and anihilation operators    
template <class T>
class B: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
     {
        _O::_is_fermion = false;
        _O::dqn.sz() = 0;
        _O::dqn.n() = 0;
        if(QN::get_qn_index("Nb") == QN_MAX_SIZE)
          _O::dqn.n() = -1;
        else
          _O::dqn["Nb"] = -1;
        _O::_is_diagonal = false;
      }
  public:
    B(): BasicOp<T>("B",0) 
      { init(); }
    B(int site, int internal_site = 0): BasicOp<T>("B",site,internal_site)  
      { init(); }
    B(const Basis& b, int site, int internal_site = 0): BasicOp<T>("B",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    B(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("B",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};
    
template <class T>
class Bd: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()   
      {
        _O::_is_fermion = false;
        _O::dqn.sz() = 0;
        _O::dqn.n() = 0;
        if(QN::get_qn_index("Nb") == QN_MAX_SIZE)
          _O::dqn.n() = 1;
        else
          _O::dqn["Nb"] = 1;
        _O::_is_diagonal = false;
      }
  public:
    Bd(): BasicOp<T>("Bd",0) 
      { init(); }
    Bd(int site, int internal_site = 0): BasicOp<T>("Bd",site,internal_site)  
      { init(); }
    Bd(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Bd",b,site,internal_site) 
      { init(); _O::init_blocks(); }
    Bd(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Bd",m,b,site,internal_site) 
      { init(); this->bmatrix_from_matrix(m); }
};

template <class T>
class Pair: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = -2; _O::_is_diagonal = false; }
  public:
    Pair(): BasicOp<T>("Pair",0)
      { init(); }
    Pair(int site, int internal_site = 0): BasicOp<T>("Pair",site,internal_site)
      { init(); }
    Pair(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Pair",b,site,internal_site)
      { init(); _O::init_blocks(); }
    Pair(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Pair",m,b,site,internal_site)
      { init(); this->bmatrix_from_matrix(m); }
;
};

template <class T>
class Paird: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()
      { _O::_is_fermion = false; _O::dqn.sz() = 0; _O::dqn.n() = 2; _O::_is_diagonal = false; }
  public:
    Paird(): BasicOp<T>("Paird",0)
      { init(); }
    Paird(int site, int internal_site = 0): BasicOp<T>("Paird",site,internal_site)
      { init(); }
    Paird(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Paird",b,site,internal_site)
      { init(); _O::init_blocks(); }
    Paird(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Paird",m,b,site,internal_site)
      { init(); this->bmatrix_from_matrix(m); }
};

template <class T>
class Double: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()
      { _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public: 
    Double(): BasicOp<T>("Double",0)
      { init(); }
    Double(int site, int internal_site = 0): BasicOp<T>("Double",site,internal_site)
      { init(); }
    Double(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Double",b,site,internal_site)
      { init(); _O::init_blocks(); }
    Double(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Double",m,b,site,internal_site)
      { init(); this->bmatrix_from_matrix(m); }
};

template <class T>
class Nup: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()
      { _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    Nup(): BasicOp<T>("Nup",0)
      { init(); }
    Nup(int site, int internal_site = 0): BasicOp<T>("Nup",site,internal_site)
      { init(); }
    Nup(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Nup",b,site,internal_site)
      { init(); _O::init_blocks(); }
    Nup(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Nup",m,b,site,internal_site)
      { init(); this->bmatrix_from_matrix(m); }
};

template <class T>
class Ndn: public dmtk::BasicOp<T>
{
  private:
    typedef BasicOp<T> _O;
    virtual void init()
      { _O::dqn.sz() = 0; _O::dqn.n() = 0; }
  public:
    Ndn(): BasicOp<T>("Ndn",0)
      { init(); }
    Ndn(int site, int internal_site = 0): BasicOp<T>("Ndn",site,internal_site)
      { init(); }
    Ndn(const Basis& b, int site, int internal_site = 0): BasicOp<T>("Ndn",b,site,internal_site)
      { init(); _O::init_blocks(); }
    Ndn(const Matrix<T> &m, const Basis& b, int site, int internal_site = 0): BasicOp<T>("Ndn",m,b,site,internal_site)
      { init(); this->bmatrix_from_matrix(m); }
};

template<class T>
inline void
BasicOp<T>::init_blocks()
{
  typename PackedBasis::iterator iter;

  typedef BMatrix<T> bm;

  for(iter = bm::subspace_begin(); iter != bm::subspace_end(); iter++){
    SubSpace col_range = *iter;
    SubSpace row_range = bm::subspace(col_range.qn()+dqn);
    if(row_range.begin() == DMTK_ERROR) continue;
    SubMatrix<T> b(col_range.qn(),col_range,row_range);
    this->push_back(b);
  }
}

template<class T>
void
BasicOp<T>::bmatrix_from_matrix(const Matrix<T>& m)
{
  typename PackedBasis::iterator iter;

  typedef BMatrix<T> bm;
  _BMvector::clear();

/*
  for(iter = bm::subspace_begin(); iter != bm::subspace_end(); iter++){
    SubSpace col_range = *iter;
    cout << "HOLA BMATRIX SUBSPACE " << col_range.qn() << "  "<< col_range.begin() << " " << col_range.end() << endl;
  }
*/

  for(iter = bm::subspace_begin(); iter != bm::subspace_end(); iter++){
    SubSpace col_range = *iter;
    SubSpace row_range = bm::subspace(col_range.qn()+dqn);
//    cout << "HOLA BMATRIX ADDING " << description() << " " << dqn << " " << col_range.qn() << " " << (col_range.qn()+dqn) << endl;
    if(row_range.begin() == DMTK_ERROR) continue;
//    cout << "HOLA ADDED\n";
    SubMatrix<T> b(col_range.qn(),col_range,row_range);
    b=(m(col_range,row_range));
    if(col_range.dim() == 1 && row_range.dim() == 1 && b(0,0) == T(0)) continue;
    this->push_back(b);
  }
}

/////////////////////////////////////////////////////////////////////
// Term: Product of two or more operators and a coefficient
/////////////////////////////////////////////////////////////////////
typedef enum
{
  TERM_EMPTY,
  TERM_LOCAL,
  TERM_PRODUCT,
  TERM_EXTERN,
  TERM_MEASURE,
} TermType;

template <class T>
class Term: public std::vector<BasicOp<T> >
{
  private:
    T _coef;
    TermType _term_type;
    T _value;  // mean value, after measurement

  public:
    typedef typename std::vector<BasicOp<T> > _V;
    typedef typename _V::const_iterator const_iterator;
    typedef typename _V::iterator iterator;

    Term(): _coef(1), _term_type(TERM_EMPTY), _value(T(0)) {}

    Term(const BasicOp<T>& local_op): _coef(1), _term_type(TERM_PRODUCT), _value(T(0)) 
      { 
        if(local_op.is_term())
          operator=(local_op.term());
        else
          { _V::clear(); _V::push_back(local_op); }
      }

    Term(const BasicOp<T>& local_op, T coef): _coef(coef), _term_type(TERM_PRODUCT), _value(T(0)) 
      { _V::clear(); _V::push_back(local_op); }

    Term(const BasicOp<T>& op1, const BasicOp<T>& op2): _coef(1), _value(T(0)), _term_type(TERM_PRODUCT) 
      { _V::clear(); _V::push_back(op2); _V::push_back(op1); }

    Term(const BasicOp<T>& _op1, const BasicOp<T>& _op2, T coef): _coef(coef), _value(T(0)), _term_type(TERM_PRODUCT)
      { _V::clear(); _V::push_back(_op2); _V::push_back(_op1); }

    Term(const Term& t): _V(t), _coef(t.coef()), _term_type(t.type()), _value(t.value()) {}

    Term& operator=(const Term& t) { _V::operator=(t); _coef = t.coef(); _term_type = t.type(); _value = t.value(); return *this; }

    Term operator*(const T& v) const 
      { Term aux(*this); aux._coef *= v; return aux; }

    Term operator*(const BasicOp<T>& op) const
      { 
        Term aux(op); 
//        if(this->size() == 1 && this->operator[](0).name() == "I") return aux;
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) aux.push_back(*iter);
        aux._term_type = TERM_PRODUCT; 
        aux.coef() = coef();
        return aux; 
      }
      
    Term operator*(const Term& t) const
      {
        Term aux = t;
//        if(this->size() == 1 && this->operator[](0).name() == "I") aux.clear();
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) aux.push_back(*iter);
        aux._coef *= _coef;
        aux._term_type = TERM_PRODUCT; 
        return aux; 
      }

    Term& operator*=(const T& v) { _coef *= v; return *this; }
    Term& operator*=(const BasicOp<T>& op) 
      { 
//        if(this->size() == 1 && this->operator[](0).name() == "I") this->clear();
        _V::push_back(op); 
        _term_type = TERM_PRODUCT; 
        return *this; 
      }
    Term& operator*=(const Term& t) 
      { 
//        if(this->size() == 1 && this->operator[](0).name() == "I") this->clear();
        const_iterator iter;
        for(iter = t.begin(); iter != t.end(); iter++) push_back(*iter);
        _coef *= t.coef();
        _term_type = TERM_PRODUCT; 
        return *this; 
      }

    QN dqn() const
    {
      QN _dqn;
      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++) {
         _dqn += (*iter).dqn;
      }
      return _dqn; 
    }

    Term reorder(bool use_sign = true) const
    {
      Term new_term;
      if(_V::size() == 1) {
        new_term = *this;
        BasicOp<T> &op = new_term[0];
        if(op.is_hami()){
          op.hami().reorder_terms();
          return new_term;
        } else if(op.is_term()){
          op.term() = op.term().reorder();
          return new_term;
        } else {
          return new_term;
        }
      }

      Term t = *this;
      Vector<int> sites(t.size());
      Vector<size_t> indx(t.size());
      Vector<bool> is_fermion(t.size());
      bool fermion = false;
      for(int i = 0; i < t.size(); i++){
        const BasicOp<T> &op = t[i];
        sites(i) = op.site();
        is_fermion(i) = op.fermion();
        if(op.fermion()) fermion = true;
      }
      bool changed = true;
      int fermion_sign = 1;
      new_term = t;
//cout << t.description() << endl;
      while(changed){
        changed = false;
        for(int i = 0; i < t.size()-1; i++){
          if(sites(i) < sites(i+1)){
            if(is_fermion(i) && is_fermion(i+1)) fermion_sign = -fermion_sign;
            int si = sites(i);
            bool fi = is_fermion(i);
            sites(i) = sites(i+1);
            is_fermion(i) = is_fermion(i+1);
            sites(i+1) = si;
            is_fermion(i+1) = fi;
            BasicOp<T> aux = new_term[i];
            new_term[i] = new_term[i+1];
            new_term[i+1] = aux; 
            changed = true;
//cout << "HOLA " << i << " " << new_term[i].description() << " " << new_term[i+1].description() << " " << sites(i) << " " << sites(i+1) << " " << is_fermion(i) << " " << is_fermion(i+1) << " " << fermion_sign << endl;
          }
        }
      }
     
      new_term.coef() = fermion && use_sign ? t.coef()*T(fermion_sign) : t.coef();
      new_term._term_type = _term_type;
      return new_term;
    }

    bool is_diagonal() const 
    { 
        bool val = true; 
        for(int i = 0; i < this->size(); i++)
           if(!_V::operator[](i).is_diagonal()) { val = false; break; }
        return val; 
    }

    Term& set_type(TermType type) { _term_type = type; return *this; }

    TermType type() const { return _term_type; }
    T coef() const { return _coef; }
    T& coef() { return _coef; }

    string name(bool use_coef = false) const 
      {
        ostringstream bf;
        if(use_coef) bf << coef() << "*";
        const_iterator iter;
        iter = _V::end();
        iter--;
        bf << (*iter).name() << "(" << (*iter).site() << "," << (*iter).internal_site() << ")"; 
        if(_V::size() > 1){
          while(1){
            iter--;
            const BasicOp<T> &op = *iter;
            bf << "*" << op.name() << "(" << op.site() << "," << op.internal_site() << ")"; 
            if(iter == _V::begin()) break;
          }
        }
        return bf.str();
      }

    string description() const 
      {
        ostringstream bf;
        const_iterator iter;
        iter = _V::end();
        iter--;
        bf << coef() << "*" << (*iter).description();
        if(_V::size() > 1){
          while(1){
            iter--;
//            bf << "*" << (*iter).description(); 
            bf << (*iter).description(); 
            if(iter == _V::begin()) break;
          }
        }
        return bf.str();
      }


    bool operator==(const BasicOp<T>& op)
      { return (name() == op.name()); }

    bool operator==(const Term<T>& t)
      {
        typename _V::iterator iter1;
        typename _V::const_iterator iter2;
        for(iter1 = _V::begin(), iter2 = t.begin(); iter1 != _V::end() && iter2 != t.end(); iter1++, iter2++){
          if(*iter1 != *iter2) return false;
        }

        return(_coef == t._coef && _V::size() == t.size());
      }

    T value() const { return _value; }
    Term<T>& set_value(T val) { _value = val; return *this; }

// Streams

    void write(std::ostream &s) const
    {
      size_t n = _V::size();
      s.write((const char *)&_coef, sizeof(T));
      s.write((const char *)&_term_type, sizeof(TermType));
      s.write((const char *)&_value, sizeof(T));
      s.write((const char *)&n, sizeof(size_t));
      const_iterator iter;
      for(iter = _V::begin(); iter != _V::end(); iter++) (*iter).write(s);
    }

    void read(std::istream &s)
    {
      size_t n;
      s.read((char *)&_coef, sizeof(T));
      s.read((char *)&_term_type, sizeof(TermType));
      s.read((char *)&_value, sizeof(T));
      s.read((char *)&n, sizeof(size_t));
      _V::clear();
      for(int i = 0; i < n; i++){
        BasicOp<T> op;
        op.read(s);
        this->push_back(op);
      }
    }
};

////////////////////////////////////////////////////////////////////////////

//////////// NEW

template<class T>
void
new_operator(BasicOp<T> &new_op, const BasicOp<T>&op,
             const BMatrix<T>& rho, const Basis &new_basis)
{
  typename BMatrix<T>::iterator biter;
  Matrix<T> aux;

  for(biter = new_op.begin(); biter != new_op.end(); biter++){

    SubMatrix<T> &block = (*biter);
    const SubMatrix<T> *_u1 = rho.block(block.qn());
    const SubMatrix<T> *_u2 = rho.block(block.qn()+new_op.dqn);

    if(!_u1 || !_u2) continue;

    const SubMatrix<T> &u1 = *_u1;
    const SubMatrix<T> &u2 = *_u2;

    SubSpace subspacei = new_basis(block.qn());
    SubSpace subspacej = new_basis(block.qn()+new_op.dqn);
    const SubMatrix<T> *m = op.block(block.qn());

    aux.reshape(u1.cols(),m->rows());
    matrix_matrix_product('N','N',*m,u1,aux);
    matrix_matrix_product('C','N',u2,aux,block,T(1),T(1));
  }

}

template<class T>
void
new_operator(BasicOp<T> &new_op, const BasicOp<T>&op,
             const BMatrix<T>& rho, const Basis &new_basis, 
             int position, T coef = T(1), bool rotate = true)
{
  QN qn(-999);
  bool start = true;
  typename BMatrix<T>::iterator biter;
  const SubMatrix<T> *_sub;
  Matrix<T> m, aux;
  CTimer clock;
  clock.Start();

//  cout << "HOLA NEW OP " << new_op.description() << " " << new_op.dqn << " " << op.description() << " " <<  op.dqn << endl;

  for(biter = new_op.begin(); biter != new_op.end(); biter++){

    SubMatrix<T> &block = (*biter);
    const SubMatrix<T> *_u1 = rho.block(block.qn());
    const SubMatrix<T> *_u2 = rho.block(block.qn()+new_op.dqn);

    if(!_u1 || !_u2) continue;

//  cout << "HOLA BLOCK " << _u1->qn() << " " << _u2->qn() << endl;

    const SubMatrix<T> &u1 = *_u1;
    const SubMatrix<T> &u2 = *_u2;

    SubSpace subspacei = new_basis(block.qn());
    SubSpace subspacej = new_basis(block.qn()+new_op.dqn);
    m.reshape(subspacei.dim(),subspacej.dim()); 
    m = T(0);
    int sign = 1;

//    cout << subspacei.begin() << " " << subspacei.end() << endl;
//    cout << subspacej.begin() << " " << subspacej.end() << endl;
//    cout << "SubMatrix " << block.qn().n() << " " << block.qn().sz() << endl;
//    clock.Lap();

    if(position == RIGHT){

      for(int col = subspacei.begin(); col <= subspacei.end(); col++){
        State si = new_basis(col);

        if((!start && si.qn2() != qn) || start){
          qn = si.qn2();
          _sub = op.block(qn);
        }

        start = false;
        if(!_sub) continue;
        const SubMatrix<T> &sub = *_sub;

        if(op.fermion()) sign = si.qn1().fermion_sign();

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){
          State sj = new_basis(row);

          if(sj.qn2() == qn+op.dqn && sj.i1 == si.i1){
            int i2 = si.i2 - sub.col_range().begin();
            int j2 = sj.i2 - sub.row_range().begin();

            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rcol,rrow) = T(sign)*coef*sub(i2,j2);
          }

        }
      }

    } else if(position == LEFT) { 

      for(int col = subspacei.begin(); col <= subspacei.end(); col++){
        State si = new_basis(col);
        if((!start && si.qn1() != qn) || start){
          qn = si.qn1();
          _sub = op.block(qn);
        }

        start = false;

        if(!_sub) continue;
        const SubMatrix<T> &sub = *_sub;

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){
          State sj = new_basis(row);

          if(sj.qn1() == qn+op.dqn && sj.i2 == si.i2){ 
            int i1 = si.i1 - sub.col_range().begin();
            int j1 = sj.i1 - sub.row_range().begin();

            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rcol,rrow) = coef*sub(i1,j1);
          }

        }
      }

    }


//    cout << "Lap: " << clock.LapTime() << endl;
//    clock.Lap();

    if(rotate){
      aux.reshape(u1.cols(),m.rows());
//      aux = product(m,u1);
//      Matrix<T> tu2(u2.rows(),u2.cols());
//      tu2 = ctranspose(u2);
//      block += product(tu2,aux);
      matrix_matrix_product('N','N',m,u1,aux);
      matrix_matrix_product('C','N',u2,aux,block,T(1),T(1));
    } else {
      block += m;
    }

/*
    for(int j = 0; j < block.col_range().dim(); j++){
      for(int i = 0; i < block.row_range().dim(); i++){
        T &val = block(j,i);
//        val = 0;
        for(int l = 0; l < u1.row_range().dim(); l++){
          for(int k = 0; k < u2.row_range().dim(); k++){
            val += std::conj(u2(i,k))*m(l,k)*u1(j,l);
          }
        }
      }
    }
*/
//    cout << "Lap: " << clock.LapTime() << endl;
  }

/*
  cout << "HOLA NEW OP " << op.description() << endl;

  typename BMatrix<T>::iterator iter;
  for(iter = new_op.begin(); iter != new_op.end(); iter++){
   cout << "BLOCK " << (*iter).qn() << endl;
   for(int i = 0; i < (*iter).row_range().size(); i++)
     for(int j = 0; j < (*iter).col_range().size(); j++)
       cout << i << " " << j << " " << (*iter)(j,i) << endl;
  }
*/
/*
  cout << "============================\n";
  for(int i = 0; i < 4; i++)
   for(int j = 0; j < 4; j++)
     cout << i << " " << j << " " << new_op(i,j) << endl;
*/
}


template<class T>
void
new_operator(BasicOp<T> &new_op, const BasicOp<T>&op1, const BasicOp<T>&op2,
             size_t m1, size_t m2,  
             const BMatrix<T>& rho, const Basis& new_basis, 
             T coef = T(1), bool use_hc = false, bool rotate = true)
{
  QN qn1(-999), qn2(-999);
  bool start = true;
  typename BMatrix<T>::iterator biter;
  const SubMatrix<T> *_sub1, *_sub2;
  Matrix<T> m, aux;
  int sign = 1;

  const BasicOp<T>* pop1 = &op1;
  const BasicOp<T>* pop2 = &op2;

  bool fermion1 = false;
  bool fermion2 = false;

  if(m2 < m1){
    pop2 = &op1;
    pop1 = &op2;
    if(op1.fermion()) fermion2 = true;        
  } else if(m2 > m1) {
    if(op2.fermion()) fermion1 = true;        
  }


  const BasicOp<T> &_op1 = *pop1;
  const BasicOp<T> &_op2 = *pop2;

//  cout << "HOLA NEW OP " << new_op.description() << " " << new_op.dqn << " " << _op1.description() << " " <<  _op1.dqn << " " << _op2.description() << _op2.dqn << endl;

  for(biter = new_op.begin(); biter != new_op.end(); biter++){

    SubMatrix<T> &block = (*biter);
    const SubMatrix<T> *_u1 = rho.block(block.qn());
    const SubMatrix<T> *_u2 = rho.block(block.qn()+new_op.dqn);

    if(!_u1 || !_u2) continue;

    const SubMatrix<T> &u1 = *_u1;
    const SubMatrix<T> &u2 = *_u2;


    SubSpace subspacei = new_basis(block.qn());
    SubSpace subspacej = new_basis(block.qn()+new_op.dqn);
    m.reshape(subspacei.dim(),subspacej.dim()); 
    m = T(0);
////////////////////////////////////
/*
    cout << "HOLA OP1 " << _op1.dqn.kx1() << endl;
    cout << "HOLA OP2 " << _op2.dqn.kx1() << endl;
    cout << "HOLA SUBSPACE I " << endl;
    for(int row = subspacei.begin(); row <= subspacei.end(); row++){ 
      State sj = new_basis(row);
cout << "HOLA I " << row << " " << sj.qn().n() << " " << sj.qn().sz() << " " << sj.qn().kx1() << " " << sj.qn1().kx1() << " " << sj.qn2().kx1() << endl;
    }
    cout << "HOLA SUBSPACE J " << endl;
    for(int row = subspacej.begin(); row <= subspacej.end(); row++){ 
      State sj = new_basis(row);
cout << "HOLA J " << row << " " << sj.qn().n() << " " << sj.qn().sz() << " " << sj.qn().kx1() << " " << sj.qn1().kx1() << " " << sj.qn2().kx1() << endl;
    }
*/
////////////////////////////////////

    for(int col = subspacei.begin(); col <= subspacei.end(); col++){ 
      State si = new_basis(col);
      if((!start && si.qn1() != qn1) || start){
        qn1 = si.qn1();
        _sub1 = _op1.block(qn1);
      }
      if(!_sub1) continue;

      if((!start && si.qn2() != qn2) || start){
        qn2 = si.qn2();
        _sub2 = _op2.block(qn2);
      }
      start = false;
      if(!_sub2 || (m1 == m2 && qn1 != qn2 + op2.dqn)) continue;
      const SubMatrix<T> &sub1 = *_sub1;
      const SubMatrix<T> &sub2 = *_sub2;

      int i1 = si.i1 - sub1.col_range().begin();
      int i2 = si.i2 - sub2.col_range().begin();

      if(fermion1) sign = qn1.fermion_sign();
      else if(fermion2) sign = qn1.fermion_sign()*_op2.dqn.fermion_sign();

      for(int row = subspacej.begin(); row <= subspacej.end(); row++){ 
        State sj = new_basis(row);

        if(sj.qn2() == qn2+_op2.dqn && sj.qn1() == qn1+_op1.dqn){ 
          int j1 = sj.i1 - sub1.row_range().begin();
          int j2 = sj.i2 - sub2.row_range().begin();

          int rcol = col - subspacei.begin();
          int rrow = row - subspacej.begin();
          m(rcol,rrow) += T(sign)*coef*sub1(i1,j1)*sub2(i2,j2);
        }

      }

      if(use_hc){

        for(int row = subspacej.begin(); row <= subspacej.end(); row++){ 
          State sj = new_basis(row);

          if(sj.qn2() == qn2+_op2.dqn && sj.qn1() == qn1+_op1.dqn){ 
            int j1 = sj.i1 - sub1.row_range().begin();
            int j2 = sj.i2 - sub2.row_range().begin();
      
            int rcol = col - subspacei.begin();
            int rrow = row - subspacej.begin();
            m(rrow,rcol) += T(sign)*std::conj(coef)*std::conj(sub1(i1,j1))*std::conj(sub2(i2,j2));
          }

        }
      }

    }

    if(rotate){   
      aux.reshape(u1.cols(),m.rows());
//      aux = product(m,u1);
//      Matrix<T> tu2(u2.rows(),u2.cols());
//      tu2 = ctranspose(u2);
//      block += product(tu2,aux);
      matrix_matrix_product('N','N',m,u1,aux);
      matrix_matrix_product('C','N',u2,aux,block,T(1),T(1));
    } else {
      block += m;
    }


/*
    for(int j = 0; j < block.col_range().dim(); j++){
      for(int i = 0; i < block.row_range().dim(); i++){
        T &val = block(j,i);
//        val = 0;
        for(int l = 0; l < u1.row_range().dim(); l++){
          for(int k = 0; k < u2.row_range().dim(); k++){
            val += std::conj(u2(i,k))*m(l,k)*u1(j,l);
          }
        }
      }
    }
*/

  }


}

////////////////////////////////////////////////////////////////////
// product:
////////////////////////////////////////////////////////////////////
template<class T>
void
product(const BasicOp<T> &the_op,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  if(the_op.name() == "I") {
    VectorState<T> aux(v);
    aux *= coef;
    res += aux;
  }

  Matrix<T> *_maux1, *_maux2;
  if(globals){
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
  } else {
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;
//------------------------------------------------------
  bool doreshape = false;
  const BasicOp<T> *_op;
  _op = &the_op;
/*
  BasicOp<T> new_op;
  if(v.qn_mask() == QN::get_qn_mask()){
  } else {
    int nmask1 = 0, nmask2 = 0;
    for(int i = 0; i < QN::QN_LAST; i++) {
      nmask1 += IBITS(v.qn_mask(),i);
      nmask2 += IBITS(QN::get_qn_mask(),i);
    }
    if(nmask1 > nmask2) {
      _op = &the_op;
    } else{
      new_op = the_op.reshape(v.qn_mask()); 
      _op = &new_op;
      doreshape = true;
    }
  } 
*/

  const BasicOp<T> &op = *_op;
//------------------------------------------------------
//
  Vector<QN> dqn(5);
  if(!hc) 
    dqn(m) += op.dqn; 
  else 
    dqn(m) -= op.dqn;

  char do_hc = (!hc) ? 'N' : 'C';
  T real_coef = (!hc) ? coef : std::conj(coef);

  const SubMatrix<T> *_block = NULL;

  typename VectorState<T>::const_iterator siter;

/*
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;
cout << "HOLA QN1 " << ss[1].qn().n() << " " << ss[1].dim() << endl;
cout << "HOLA QN2 " << ss[2].qn().n() << " " << ss[2].dim() << endl;
cout << "HOLA QN3 " << ss[3].qn().n() << " " << ss[3].dim() << endl;
cout << "HOLA QN4 " << ss[4].qn().n() << " " << ss[4].dim() << endl;
  }
*/

  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;


    if(!hc){
        _block = op.block(ss[m].qn());
    }else{
        _block = op.block(ss[m].qn()-op.dqn);
    }

    if(!_block) continue;

    const SubMatrix<T> &block(*_block);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    int sign = 1;
    if(op.fermion()){
      for(int ib = m-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }

    switch(m){
      case BLOCK1:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux2.reshape(subres.size2(),subres.size1());
            maux1 = subv;
            matrix_matrix_product(do_hc,'T',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres.transpose() += maux2.array();
          }
        }
        break;
      case BLOCK2:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size2(),subres.size1());
            matrix_matrix_product(do_hc,'T',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres.transpose() += maux2.array();
          }
        }
/* This is much slower !!!
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),slice(0,res_slice.size3(),1),i4);
            maux1 = subv;
            maux2.reshape(subres.size2(),subres.size1());
            matrix_matrix_product(do_hc,'T',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres.transpose() += maux2.array();
          }
        }
*/
        break;
      case BLOCK3:
        for(int i1 = 0; i1 < ss[1].dim(); i1++){
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),slice(0,v_slice.size3(),1),i4);
            gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),slice(0,res_slice.size3(),1),i4);
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres += maux2.array();
          }
        }
        break;
      case BLOCK4:
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,slice(0,v_slice.size4(),1));
            gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,slice(0,res_slice.size4(),1));
            maux1 = subv;
            maux2.reshape(subres.size1(),subres.size2());
            matrix_matrix_product(do_hc,'N',static_cast<Matrix<T> >(block),maux1,maux2,T(sign)*real_coef);
            subres += maux2.array();
          }
        }
        break;
    }
  }
  if(!globals){
    delete(_maux1);
    delete(_maux2);
  }
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL, bool use_condensed = false)
{
  int _mask = mask(m1,m2);
  int _m1, _m2;
  Vector<QN> dqn(5);
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<T> *_vaux1, *_vaux2;
  Matrix<T> *_maux1, *_maux2, *_maux3;
  if(globals){
    _vaux1 = &globals->v1;
    _vaux2 = &globals->v2;
    _maux1 = &globals->m1;
    _maux2 = &globals->m2;
    _maux3 = &globals->m3;
  } else {
    _vaux1 = new Vector<T>(MIN_VECTOR_SIZE);
    _vaux2 = new Vector<T>(MIN_VECTOR_SIZE);
    _maux1 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux2 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
    _maux3 = new Matrix<T>(MIN_MATRIX_SIZE,MIN_MATRIX_SIZE);
  }
  Vector<T> &vaux1 = *_vaux1;
  Vector<T> &vaux2 = *_vaux2;
  Matrix<T> &maux1 = *_maux1;
  Matrix<T> &maux2 = *_maux2;
  Matrix<T> &maux3 = *_maux3;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  int sign0 = 1;
  if(m1 <= m2){
    _op1 = &op1;
    _op2 = &op2;
    _m1 = m1;
    _m2 = m2;
  } else {
    _op1 = &op2;
    _op2 = &op1;
    _m1 = m2;
    _m2 = m1;
    if(_op1->fermion() && _op2->fermion()) sign0 = -1;
  }

/*
  bool doreshape = false;
  BasicOp<T> new_op1, new_op2;
  if(v.qn_mask() != QN::get_qn_mask()){
    int nmask1 = 0, nmask2 = 0;
    for(int i = 0; i < QN::QN_LAST; i++) {
      nmask1 += IBITS(v.qn_mask(),i);
      nmask2 += IBITS(QN::get_qn_mask(),i);
    }
    if(nmask1 < nmask2) {
      new_op1 = _op1->reshape(v.qn_mask()); 
      new_op2 = _op2->reshape(v.qn_mask()); 
      _op1 = &new_op1;
      _op2 = &new_op2;
      for(int i = 0; i < QN::QN_LAST; i++)
        if(IBITS(v.qn_mask(),i) == 0) 
          { 
            dqn[1][i] = 0; dqn[2][i] = 0; dqn[3][i] = 0; dqn[4][i] = 0; 
            dqn2[1][i] = 0; dqn2[2][i] = 0; dqn2[3][i] = 0; dqn2[4][i] = 0; 
          }
      doreshape = true;
    }
  } 
*/

/*
  cout << "HOLA " << _op1->description() << endl;
  for(int i = 0; i < _op1->size(); i++){
    const SubMatrix<T> *s = _op1->operator[](i);
    if(!s) cout << "ERROR\n";
    cout << s->qn()[0] << " " << s->qn()[1] << " " << s->qn()[2] << " " << s->qn()[3] << endl;
  }
  cout << "HOLA " << _op2->description() << endl;
  for(int i = 0; i < _op2->size(); i++){
    const SubMatrix<T> *s = _op2->operator[](i);
    if(!s) cout << "ERROR\n";
    cout << s->qn()[0] << " " << s->qn()[1] << " " << s->qn()[2] << " " << s->qn()[3] << endl;
  }
*/
/* Test 
  if(m1 == m2){
    VectorState<T> aux(v.b1(),v.b2(),v.b3(),v.b4(),v.qn()+op2.dqn, v.qn_mask());
    product(op2,v,aux,m2,T(1),false);
    product(op1,aux,res,m1,coef,false);
    if(hc){
      aux.set_qn_mask(v.qn()-op1.dqn, v.qn_mask());
      aux = T(0);
      product(op1,v,aux,m1,T(1),true);
      product(op2,aux,res,m2,coef,true);
    }
    return;
  }
*/

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;

  typename VectorState<T>::const_iterator siter;


/*
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;
cout << "HOLA QN1 " << ss[1].qn()[0] << " " << ss[1].qn()[1] << " " << ss[1].qn()[2] << " " << ss[1].qn()[3] << endl; 
cout << "HOLA QN2 " << ss[2].qn()[0] << " " << ss[2].qn()[1] << " " << ss[2].qn()[2] << " " << ss[2].qn()[3] << endl; 
cout << "HOLA QN3 " << ss[3].qn()[0] << " " << ss[3].qn()[1] << " " << ss[3].qn()[2] << " " << ss[3].qn()[3] << endl; 
cout << "HOLA QN4 " << ss[4].qn()[0] << " " << ss[4].qn()[1] << " " << ss[4].qn()[2] << " " << ss[4].qn()[3] << endl; 
//    cout << "HOLA SUBSPACE QN " << ss[1].qn()[0] << " " << ss[2].qn()[0] << " " << ss[3].qn()[0] << " " << ss[4].qn()[0] << " " << ss.dim() << endl;
  }
*/
/*
  for(siter = res.subspace_begin(); siter != res.subspace_end(); siter++){
    StateSpace ss = *siter;
    cout << "HOLA RES SUBSPACE QN " << ss[1].qn()[0] << " " << ss[2].qn()[0] << " " << ss[3].qn()[0] << " " << ss[4].qn()[0] << " " << ss.dim() << endl;
  }
*/
/*
  for(int i = 0; i < _op1->size(); i++){
    const SubMatrix<T> *sub = _op1->operator[](i);
    cout << "HOLA SUB 1 " << sub->qn()[0] << " " << sub->qn()[1] << " " << sub->qn()[2] << " " << sub->qn()[3] << endl;
  }
  for(int i = 0; i < _op2->size(); i++){
    const SubMatrix<T> *sub = _op2->operator[](i);
    cout << "HOLA SUB 2 " << sub->qn()[0] << " " << sub->qn()[1] << " " << sub->qn()[2] << " " << sub->qn()[3] << endl;
  }
*/

  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

     _block2 = _op2->block(ss[_m2].qn());
    if(m1 != m2){
        _block1 = _op1->block(ss[_m1].qn());
    } else {
        _block1 = _op1->block(ss[_m1].qn() + _op2->dqn);
    }
    if(!_block1 || !_block2) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss); 
    if(v_slice_hc.size() == 0 || res_slice_hc.size() == 0) continue; 

    int sign = sign0;
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign();
    }


    switch(_mask){
      case (MASK_BLOCK1|MASK_BLOCK2):
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),
                                         slice(0,v_slice.size2(),1),i3,i4);
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),
                                            slice(0,res_slice.size2(),1),i3,i4);
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),
                                              slice(0,v_slice_hc.size2(),1),i3,i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),
                                                 slice(0,res_slice_hc.size2(),1),i3,i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK1|MASK_BLOCK3):
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,
                                         slice(0,v_slice.size3(),1),i4);
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,
                                            slice(0,res_slice.size3(),1),i4);
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,
                                              slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,
                                                 slice(0,res_slice_hc.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK1|MASK_BLOCK4):
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cgslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,
                                            slice(0,res_slice.size4(),1));
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK2|MASK_BLOCK3):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),
                                         slice(0,v_slice.size3(),1),i4);
          gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),
                                            slice(0,res_slice.size3(),1),i4);
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),
                                           slice(0,v_slice_hc.size3(),1),i4);
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),
                                              slice(0,res_slice_hc.size3(),1),i4);
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK2|MASK_BLOCK4):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cgslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,
                                            slice(0,res_slice.size4(),1));
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cgslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case (MASK_BLOCK3|MASK_BLOCK4):
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++){
          cgslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),
                                         slice(0,v_slice.size4(),1));
          gslice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),
                                            slice(0,res_slice.size4(),1));
          if(use_condensed){
            product(block1,block2,subv,subres,maux3,coef*T(sign),T(1));
          } else {
            maux1 = subv;
            maux2 = subres;
            product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1));
            subres = maux2.array();
          }
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++){
            cgslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),
                                              slice(0,v_slice_hc.size4(),1));
            gslice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),
                                                 slice(0,res_slice_hc.size4(),1));
            if(use_condensed){
              product(block1,block2,subv,subres,maux3,coef*T(sign),T(1),true);
            } else {
              maux1 = subv;
              maux2 = subres;
              product(block1,block2,maux1,maux2,maux3,coef*T(sign),T(1),true);
              subres = maux2.array();
            }
          }
        }
        break;
      case MASK_BLOCK1:
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(slice(0,v_slice.size1(),1),i2,i3,i4);
          slice_iter<T> subres = res_slice(slice(0,res_slice.size1(),1),i2,i3,i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(slice(0,v_slice_hc.size1(),1),i2,i3,i4);
            slice_iter<T> subres = res_slice_hc(slice(0,res_slice_hc.size1(),1),i2,i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK2:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(i1,slice(0,v_slice.size2(),1),i3,i4);
          slice_iter<T> subres = res_slice(i1,slice(0,res_slice.size2(),1),i3,i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,slice(0,v_slice_hc.size2(),1),i3,i4);
            slice_iter<T> subres = res_slice_hc(i1,slice(0,res_slice_hc.size2(),1),i3,i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK3:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i4 = 0; i4 < ss[4].dim(); i4++){
          cslice_iter<T> subv = v_slice(i1,i2,slice(0,v_slice.size3(),1),i4);
          slice_iter<T> subres = res_slice(i1,i2,slice(0,res_slice.size3(),1),i4);
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,slice(0,v_slice_hc.size3(),1),i4);
            slice_iter<T> subres = res_slice_hc(i1,i2,slice(0,res_slice_hc.size3(),1),i4);
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
      case MASK_BLOCK4:
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++){
          cslice_iter<T> subv = v_slice(i1,i2,i3,slice(0,v_slice.size4(),1));
          slice_iter<T> subres = res_slice(i1,i2,i3,slice(0,res_slice.size4(),1));
          vaux1 = subv;
          vaux2 = subres;
          product(block1,block2,vaux1,vaux2,coef,T(1));
          subres = vaux2.array();
        }
        if(hc){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++){
            cslice_iter<T> subv = v_slice_hc(i1,i2,i3,slice(0,v_slice_hc.size4(),1));
            slice_iter<T> subres = res_slice_hc(i1,i2,i3,slice(0,res_slice_hc.size4(),1));
            vaux1 = subv;
            vaux2 = subres;
            product(block1,block2,vaux1,vaux2,coef,T(1),true);
            subres = vaux2.array();
          }
        }
        break;
    }
  }
  if(!globals){
    delete(_vaux1);
    delete(_vaux2);
    delete(_maux1);
    delete(_maux2);
    delete(_maux3);
  }
}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, const BasicOp<T>& op4,
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3, size_t m4,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3, _m4;
  dqn(m4) += op4.dqn;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;
  Vector<QN> dqn4(5);
  dqn4(m4) += op4.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;
  const BasicOp<T> *_op4;

  Vector<size_t> m(4), indx(4);
  Vector<const BasicOp<T>* > ops(4);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  ops(3) = &op4;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  m(3) = int(m4);
  indexx<size_t,Vector<size_t> >(4, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _m4 = (size_t)m(indx(3));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));
  _op4 = ops(indx(3));

/*
  bool doreshape = false;
  BasicOp<T> new_op1, new_op2, new_op3, new_op4;
  if(v.qn_mask() != QN::get_qn_mask()){
    int nmask1 = 0, nmask2 = 0;
    for(int i = 0; i < QN::QN_LAST; i++) {
      nmask1 += IBITS(v.qn_mask(),i);
      nmask2 += IBITS(QN::get_qn_mask(),i);
    }
    if(nmask1 < nmask2) {
      new_op1 = _op1->reshape(v.qn_mask());
      new_op2 = _op2->reshape(v.qn_mask());
      new_op3 = _op3->reshape(v.qn_mask());
      new_op4 = _op4->reshape(v.qn_mask());
      _op1 = &new_op1;
      _op2 = &new_op2;
      _op3 = &new_op3;
      _op4 = &new_op4;
      for(int i = 0; i < QN::QN_LAST; i++)
        if(IBITS(v.qn_mask(),i) == 0) 
          { 
            dqn[1][i] = 0; dqn[2][i] = 0; dqn[3][i] = 0; dqn[4][i] = 0; 
            dqn2[1][i] = 0; dqn2[2][i] = 0; dqn2[3][i] = 0; dqn2[4][i] = 0; 
            dqn3[1][i] = 0; dqn3[2][i] = 0; dqn3[3][i] = 0; dqn3[4][i] = 0; 
            dqn4[1][i] = 0; dqn4[2][i] = 0; dqn4[3][i] = 0; dqn4[4][i] = 0; 
          }
      doreshape = true;
    }
  }
*/

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;
  const SubMatrix<T> *_block3 = NULL;
  const SubMatrix<T> *_block4 = NULL;

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

     _block4 = _op4->block(ss[_m4].qn());
     _block3 = _op3->block(ss[_m3].qn());
     _block2 = _op2->block(ss[_m2].qn());
     _block1 = _op1->block(ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3 || !_block4) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);
    const SubMatrix<T> &block3(*_block3);
    const SubMatrix<T> &block4(*_block4);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(hc && (v_slice_hc.size() == 0 || res_slice_hc.size() == 0)) continue; 

    int sign = 1;
    if(op4.fermion()){
      for(int ib = m4-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign()*dqn4[ib].fermion_sign();
    }

    for(int i1 = 0; i1 < ss[1].dim(); i1++)
      for(int i2 = 0; i2 < ss[2].dim(); i2++)
        for(int i3 = 0; i3 < ss[3].dim(); i3++)
          for(int i4 = 0; i4 < ss[4].dim(); i4++)
            for(int j1 = 0; j1 < res_slice.size1(); j1++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3)*block4(i4,j4);
    if(hc){
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      res_slice_hc(i1,i2,i3,i4) += std::conj(coef)*T(sign)*v_slice_hc(j1,j2,j3,j4)*std::conj(block1(i1,j1)*block2(i2,j2)*block3(i3,j3)*block4(i4,j4));
    }

  }

}

template<class T>
void
product(const BasicOp<T> &op1, const BasicOp<T>& op2,
        const BasicOp<T> &op3, 
        const VectorState<T> &v, VectorState<T> &res,
        size_t m1, size_t m2, size_t m3,
        T coef = T(1), bool hc = false, DMTKglobals<T> *globals = NULL)
{
  Vector<QN> dqn(5);
  size_t _m1, _m2, _m3;
  dqn(m3) += op3.dqn;
  dqn(m2) += op2.dqn;
  dqn(m1) += op1.dqn;
  Vector<QN> dqn2(5);
  dqn2(m2) += op2.dqn;
  Vector<QN> dqn3(5);
  dqn3(m3) += op3.dqn;

  const BasicOp<T> *_op1;
  const BasicOp<T> *_op2;
  const BasicOp<T> *_op3;

  Vector<size_t> m(3), indx(3);
  Vector<const BasicOp<T>* > ops(3);
  ops(0) = &op1;
  ops(1) = &op2;
  ops(2) = &op3;
  m(0) = int(m1);
  m(1) = int(m2);
  m(2) = int(m3);
  indexx<size_t,Vector<size_t> >(3, m, indx);
  _m1 = (size_t)m(indx(0));
  _m2 = (size_t)m(indx(1));
  _m3 = (size_t)m(indx(2));
  _op1 = ops(indx(0));
  _op2 = ops(indx(1));
  _op3 = ops(indx(2));

/*
  bool doreshape = false;
  BasicOp<T> new_op1, new_op2, new_op3;
  if(v.qn_mask() != QN::get_qn_mask()){
    int nmask1 = 0, nmask2 = 0;
    for(int i = 0; i < QN::QN_LAST; i++) {
      nmask1 += IBITS(v.qn_mask(),i);
      nmask2 += IBITS(QN::get_qn_mask(),i);
    }
    if(nmask1 < nmask2) {
      new_op1 = _op1->reshape(v.qn_mask());
      new_op2 = _op2->reshape(v.qn_mask());
      new_op3 = _op3->reshape(v.qn_mask());
      _op1 = &new_op1;
      _op2 = &new_op2;
      _op3 = &new_op3;
      for(int i = 0; i < QN::QN_LAST; i++)
        if(IBITS(v.qn_mask(),i) == 0) 
          { 
            dqn[1][i] = 0; dqn[2][i] = 0; dqn[3][i] = 0; dqn[4][i] = 0; 
            dqn2[1][i] = 0; dqn2[2][i] = 0; dqn2[3][i] = 0; dqn2[4][i] = 0; 
            dqn3[1][i] = 0; dqn3[2][i] = 0; dqn3[3][i] = 0; dqn3[4][i] = 0; 
          }
      doreshape = true;
    }
  }
*/

  const SubMatrix<T> *_block1 = NULL;
  const SubMatrix<T> *_block2 = NULL;
  const SubMatrix<T> *_block3 = NULL;

  typename VectorState<T>::const_iterator siter;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++){
    StateSpace ss = *siter;

    _block3 = _op3->block(ss[_m3].qn());
    _block2 = _op2->block(ss[_m2].qn());
    _block1 = _op1->block(ss[_m1].qn());

    if(!_block1 || !_block2 || !_block3) continue;

    const SubMatrix<T> &block1(*_block1);
    const SubMatrix<T> &block2(*_block2);
    const SubMatrix<T> &block3(*_block3);

    cstate_slice<T> v_slice = v(ss);
    state_slice<T> res_slice = res(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    if(res_slice.size() == 0) continue; 

    cstate_slice<T> v_slice_hc = v(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                   ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
    state_slice<T> res_slice_hc = res(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
    if(hc && (v_slice_hc.size() == 0 || res_slice_hc.size() == 0)) continue; 

    int sign = 1;
    if(op3.fermion()){
      for(int ib = m3-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign();
    }
    if(op2.fermion()){
      for(int ib = m2-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn3[ib].fermion_sign();
    }
    if(op1.fermion()){
      for(int ib = m1-1; ib >= 1; ib--) sign *= ss[ib].qn().fermion_sign()*dqn2[ib].fermion_sign()*dqn3[ib].fermion_sign();
    }

    if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK3){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      res_slice(j1,j2,j3,i4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i3,j3);
      if(hc){
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i4 = 0; i4 < ss[4].dim(); i4++)
                for(int j1 = 0; j1 < res_slice.size1(); j1++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                        res_slice_hc(i1,i2,i3,i4) += std::conj(coef)*T(sign)*v_slice_hc(j1,j2,j3,i4)*std::conj(block1(i1,j1)*block2(i2,j2)*block3(i3,j3));
      }
    } else if(_m1 == BLOCK1 && _m2 == BLOCK2 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,j2,i3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i2,j2)*block3(i4,j4);
      if(hc){
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i4 = 0; i4 < ss[4].dim(); i4++)
                for(int j1 = 0; j1 < res_slice.size1(); j1++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      res_slice_hc(i1,i2,i3,i4) += std::conj(coef)*T(sign)*v_slice_hc(j1,j2,i3,j4)*std::conj(block1(i1,j1)*block2(i2,j2)*block3(i4,j4));
      }
    } else if(_m1 == BLOCK1 && _m2 == BLOCK3 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j1 = 0; j1 < res_slice.size1(); j1++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(j1,i2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i1,j1)*block2(i3,j3)*block3(i4,j4);
      if(hc){
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i4 = 0; i4 < ss[4].dim(); i4++)
                for(int j1 = 0; j1 < res_slice.size1(); j1++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      res_slice_hc(i1,i2,i3,i4) += std::conj(coef)*T(sign)*v_slice_hc(j1,i2,j3,j4)*std::conj(block1(i1,j1)*block2(i3,j3)*block3(i4,j4));
      }
    } else if(_m1 == BLOCK2 && _m2 == BLOCK3 && _m3 == BLOCK4){ 
      for(int i1 = 0; i1 < ss[1].dim(); i1++)
        for(int i2 = 0; i2 < ss[2].dim(); i2++)
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int j2 = 0; j2 < res_slice.size2(); j2++)
                for(int j3 = 0; j3 < res_slice.size3(); j3++)
                  for(int j4 = 0; j4 < res_slice.size4(); j4++)
                    res_slice(i1,j2,j3,j4) += coef*T(sign)*v_slice(i1,i2,i3,i4)*block1(i2,j2)*block2(i3,j3)*block3(i4,j4);
      if(hc){
        for(int i1 = 0; i1 < ss[1].dim(); i1++)
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i4 = 0; i4 < ss[4].dim(); i4++)
                for(int j2 = 0; j2 < res_slice.size2(); j2++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      res_slice_hc(i1,i2,i3,i4) += std::conj(coef)*T(sign)*v_slice_hc(i1,j2,j3,j4)*std::conj(block1(i2,j2)*block2(i3,j3)*block3(i4,j4));
      }
    }
  }

}

//////////////////////////////////////////////////////////////////////////
template<class T>
BasicOp<T>
product(const BasicOp<T> &op2, const BasicOp<T> &op1)
{
  BasicOp<T> res;
  QN dqn1 = op1.dqn;
  QN dqn2 = op2.dqn;
  res.dqn = dqn1+dqn2;
  res.repack(op2.subspaces());

  typename BMatrix<T>::const_iterator iter;
  for(iter = op1.begin(); iter != op1.end(); iter++){
    const SubMatrix<T> &sm1 = (*iter);
    QN qn1 = sm1.qn();
    const SubMatrix<T> *_sm2 = op2.block(qn1+dqn1);
    if(_sm2){
      const SubMatrix<T> &sm2 = *_sm2;
      if(sm1.rows() == sm2.cols()){
        SubMatrix<T> sm(qn1,sm1.col_range(),sm2.row_range());
        sm=(product(sm2,sm1));
        res.push_back(sm); 
      } else {
        cout << "ERROR : op1 and op2 are inconsistent\n";
      }
    }  
  }
  return res;
}
//////////////////////////////////////////////////////////////////////////
template<class T>
BasicOp<T> 
BasicOp<T>::reshape(int qn_mask) const
{
  BasicOp<T> new_op;
  new_op = this->internals(); 
  for(int i = 0; i < QN::QN_LAST; i++)
    if(IBITS(qn_mask,i) == 0) new_op.dqn[i] = 0;

  PackedBasis::const_iterator biter;
  PackedBasis new_subspaces;

  new_subspaces = this->subspaces().reshape(qn_mask);
  new_op.resize(new_subspaces);

/*
  cout << "=========================\n";
  for(biter = new_subspaces.begin(); biter != new_subspaces.end(); biter++)
   cout << "HOLA BEGIN " << (*biter).start() << " END " << (*biter).end() << " QN " << (*biter).qn().n() << " " << (*biter).dim() << endl;
  cout << "=========================\n";
*/

  typename BMatrix<T>::iterator iter;
  typename BMatrix<T>::const_iterator citer;
  for(iter = new_op.begin(); iter != new_op.end(); iter++){
    SubMatrix<T> &sub = (*iter);
    QN new_qn = sub.qn();

    for(citer = this->begin(); citer != this->end(); citer++){
      const SubMatrix<T> &m = (*citer);
      if(new_qn.equal(m.qn(), qn_mask)){
        for(int col = 0; col < m.cols(); col++){
          int new_col = col + m.col_range().begin()-sub.col_range().begin();
          for(int row = 0; row < m.rows(); row++){
            int new_row = row + m.row_range().begin()-sub.row_range().begin();
            sub(new_col,new_row) = m(col,row);
          }
        }
      } 
    }
  } 

  return new_op;
}
//////////////////////////////////////////////////////////////////////////

} // dmtk

#endif // __DMTK_OPERATORS_H__
