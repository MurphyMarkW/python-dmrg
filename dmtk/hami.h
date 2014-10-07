#ifndef __DMTK_HAMI_H__
#define __DMTK_HAMI_H__

#include "operators.h"
#include "block.h"

namespace dmtk{

template <class T>
class Hami : public std::vector<Term<T> > 
{
  private:
    Lattice _lattice;
    std::string _name;
    T _val;
    bool _use_hc;
  public:
    typedef typename std::vector<Term<T> > _V;
    typedef typename std::vector<Term<T> >::const_iterator const_iterator;
    typedef typename std::vector<Term<T> >::iterator iterator;

    Block<T> site;
    Vector<Block<T>* > sites;

    Hami(const char *__name = NULL) : _val(0), _use_hc(false)
      { if(__name) _name = std::string(__name); else _name = std::string("H"); };
    Hami(const Lattice& lattice, const char *__name = NULL): _lattice(lattice), _val(0), _use_hc(false) 
      { 
        sites.resize(lattice.size()); 
        if(__name) _name = std::string(__name); else _name = std::string("H"); 
      };
    Hami(const Hami& h): _V(h), _lattice(h._lattice), site(h.site), sites(h.sites), _name(h._name), _val(h._val), _use_hc(h._use_hc)
      { sites.resize(_lattice.size()); } 

    Hami& operator=(const Hami<T> &h)
      {
        _V::operator=(h);
        _lattice = h._lattice;
        site = h.site;
        sites = h.sites;
        _name = h._name;
        _val = h._val;
        _use_hc = h._use_hc;
        return *this; 
      }
/*
    ~Hami() 
      { 
        for(int i = 0; i < _lattice.size(); i++)
          { if(sites[i]) delete sites[i]; sites[i] = NULL; }
      }
*/
    Hami& operator+=(const Term<T>& t) 
      { this->push_back(t); return *this; }
    Hami& operator-=(const Term<T>& t) 
      { t.coef() = - t.coef(); this->push_back(t); return *this; }
    Hami& operator+=(const BasicOp<T>& t) 
      { this->push_back(Term<T>(t)); return *this; }
    Hami& operator-=(const BasicOp<T>& t) 
      { this->push_back(Term<T>(t,T(-1))); return *this; }
    Hami& operator+=(const Hami& h) 
      {
        const_iterator iter;
        for(iter = h.begin(); iter != h.end(); iter++)
          this->push_back(*iter);
        return *this;
      }
    Hami& operator*=(const T& v) 
      {
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) (*iter) *= v;
        return *this;
      }

    Hami& operator*=(const BasicOp<T>& op) 
      {
        iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++) (*iter) *= op;
        return *this;
      }

    Hami operator*(const Hami& h) const
      {
        Hami<T> res(*this);
        res.clear();
        const_iterator iter1;
        const_iterator iter2;
        for(iter1 = _V::begin(); iter1 != _V::end(); iter1++)
          for(iter2 = h.begin(); iter2 != h.end(); iter2++)
            res += ((*iter1)*(*iter2));
        return res;
      }
    Hami operator*(const T& v) const
      {
        Hami<T> res(*this);
        res.clear();
        const_iterator iter;
        for(iter = _V::begin(); iter != _V::end(); iter++)
          res += (*iter)*v;
        return res;
      }
    Hami operator+(const Hami& h) const
      {
        Hami<T> res(*this);
        const_iterator iter;
        for(iter = h.begin(); iter != h.end(); iter++)
          res += (*iter);
        return res;
      }
 
    const Lattice &lattice() const { return _lattice; }
    Lattice &lattice() { return _lattice; }
    Block<T>& get_site(int pos) 
      { 
        if(pos >= _lattice.size()) { return *(empty_block<T>::get()); };
        if(this->sites[pos]) {
           return *this->sites[pos]; 
        } else {
           return site;
        } 
      }
    const Block<T>& get_site(int pos) const
      { 
        if(pos >= _lattice.size()) { return *(empty_block<T>::get()); };
        if(this->sites[pos]) {
           return *this->sites[pos]; 
        } else { 
           return site;
        } 
      }

    std::string const& name() const { return _name; }; 
    Hami& set_name(const char *name) 
      { 
         std::string old_name = _name;
         _name = std::string(name); 
/*
         typename Vector<Block<T> *>::iterator iter;
         typename Block<T>::iterator biter;
         for(iter = this->sites.begin(); iter != this->sites.end(); iter++){
           Block<T> &block = *(*iter);
           for(biter = block.begin(); biter != block.end(); biter++){
             BasicOp<T> &op = *biter;
             if(op.name() == old_name) op.set_name(_name.c_str());
           }
         }
         for(biter = this->site.begin(); biter != this->site.end(); biter++){
           BasicOp<T> &op = *biter;
           if(op.name() == old_name) op.set_name(_name.c_str());
         }
*/
         return *this; 
      }; 
    T value() const { return _val; }
    Hami& set_value(T val) { _val = val; return *this; }

    Hami& set_use_hc(bool use) { _use_hc = use; return *this; }
    bool use_hc() const { return _use_hc; }

    Hami& reorder_terms(bool use_sign = true); 

    Hami find_common_factors() const;
    std::string description() const;
};

template<class T>
Hami<T> 
Hami<T>::find_common_factors() const
{
  Hami<T> new_hami;
  const_iterator iter1;
  const_iterator iter2;

  bool found[_V::size()];
  for(int i = 0; i < _V::size(); i ++) found[i] = false;

  cout << "H = " << endl;
  int n1 = 0;
  for(iter1 = _V::begin(); iter1 != _V::end(); iter1++){
    const Term<T> &t1 = *iter1;
  
    if(t1.size() != 2) continue;

    const BasicOp<T> &op11 = t1[0];
    const BasicOp<T> &op12 = t1[1];

    if(found[n1] == false){
      Hami<T> terms1;
      terms1.clear();
      int n2 = 0;
      for(iter2 = _V::begin(); iter2 != _V::end(); iter2++){
        const Term<T> &t2 = *iter2;
    
        if(t2.size() != 2) continue;
  
        const BasicOp<T> &op21 = t2[0];
        const BasicOp<T> &op22 = t2[1];
  
        if(op11 == op21 && found[n2] == false) { 
          Term<T> aux = op22 * t2.coef();
          terms1.push_back(aux);
          found[n2] = true;
        }
        n2++;
      }


// OUTPUT: Show terms
      if(terms1.size() == 1){
        cout << " + " << t1.name(true) << endl;    
        new_hami += t1;
      } else {
        BasicOp<T> add_op;
        add_op.set_name(t1[1].name().c_str());
        add_op.set_type(OP_ADDITIVE);
        add_op.set_site(0);
        add_op.set_internal_site(0);
        add_op._is_hami = true;
        add_op.dqn = -op11.dqn(); //QN(-op11.dqn.n(),-op11.dqn.sz());
        add_op._hami = terms1;

        Term<T> new_term;
        new_term = op11*add_op;
        new_hami += new_term;
/*
        cout << " + " << op11.name() << "(" << op11.site() << "," << op11.internal_site() << ") * (";
        for(int i = 0; i < terms1.size(); i++){
          cout << " + " << terms1[i].name(true); 
        }
        cout <<")" << endl;
*/
//        cout << " + " << new_term.name(true) << endl;
//        cout << " + " << add_op.description() << endl;
//        cout << " + " << new_term.description() << endl;
      }
//

    }
    n1++;
  }

  cout << new_hami.description() << endl;

  return new_hami;
}

template<class T>
Hami<T> &
Hami<T>::reorder_terms(bool use_sign)
{
  typename Hami<T>::iterator iter;
  for(iter = _V::begin(); iter != _V::end(); iter++){
    Term<T> &t = *iter;
    t = t.reorder(use_sign);
/*
    Vector<int> sites(t.size());
    Vector<size_t> indx(t.size());
    bool fermion = false;
    for(int i = 0; i < t.size(); i++){
      BasicOp<T> &op = t[i];
      sites(i) = op.site();
      if(op.fermion()) fermion = true;
    }
    int nex = indexx2<int, Vector<int> >(t.size(), sites, indx, false);
    if(nex > 0){
      Term<T> new_term;
      for(int i = 0; i < t.size(); i++){
        new_term *= t[indx(i)];
      }
      new_term.coef() = fermion && use_sign ? t.coef()*T(SGN(nex)) : t.coef();
      t.clear();
      t = new_term;
    }
*/
  }

  return *this;
}

template<class T>
Hami<T>
operator*(const BasicOp<T> &op, const Hami<T> & h)
{
  Hami<T> aux;

  typename Hami<T>::const_iterator iter;
  for(iter = h.begin(); iter != h.end(); iter++){
    aux += op*(*iter);
  }  
  return aux;
}

template<class T>
Hami<T>
operator*(const Hami<T> h, const BasicOp<T> &op)
{
  Hami<T> aux;

  typename Hami<T>::const_iterator iter;
  for(iter = h.begin(); iter != h.end(); iter++){
    aux += (*iter)*op;
  }  
  return aux;
}

template<class T>
std::string
Hami<T>::description() const
{
  ostringstream bf;

  if(_V::size() == 0) return std::string("");

  bf << _V::operator[](0).description(); 
  const_iterator iter;
  iter = _V::begin();
  iter++;
  for(; iter != _V::end(); iter++){
    bf << " + " << (*iter).description();
  }  

  return bf.str();
}

Hami<complex<double> >
jspin(int i, int j)
{
  Hami<complex<double> > J;
  char name[2000];
  sprintf(name,"jspin(%i,%i)",i,j);
  J.set_name(name);
  J += Splus<complex<double> >(i) * Sminus<complex<double> >(j);
  J += Sminus<complex<double> >(i) * Splus<complex<double> >(j) * complex<double> (-1.);

  J *= complex<double>(0.,0.5);
  return J;
}

Hami<complex<double> >
Jspin(int i, int j, int k, double d = 1.)
{
  Hami<complex<double> > J;
  char name[2000];
  sprintf(name,"Jspin(%i,%i,%i)",i,j,k);
  J.set_name(name);
  J += Sz<complex<double> >(i)*Splus<complex<double> >(j)*Sminus<complex<double> >(k) * complex<double>(d);
  J += Sz<complex<double> >(i)*Sminus<complex<double> >(j)*Splus<complex<double> >(k)* complex<double>(-d);
  J += Sminus<complex<double> >(i)*Sz<complex<double> >(j)*Splus<complex<double> >(k);
  J += Splus<complex<double> >(i)*Sz<complex<double> >(j)*Sminus<complex<double> >(k) * complex<double>(-1);
  J += Splus<complex<double> >(i)*Sminus<complex<double> >(j)*Sz<complex<double> >(k) * complex<double>(d);
  J += Sminus<complex<double> >(i)*Splus<complex<double> >(j)*Sz<complex<double> >(k) * complex<double>(-d);
  J *= complex<double>(0.,0.5);
  return J;
}

} // namespace dmtk

#endif // __DMTK_HAMI_H__
