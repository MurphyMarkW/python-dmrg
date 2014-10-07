#ifndef __DMTK_QN_H__
#define __DMTK_QN_H__

#include <iosfwd>
#include "constants.h"
#include "bits.h"
#include "meta.h"

namespace dmtk 
{

// For backward compatibility with dmtk-1
enum
{
  MASK_QN_N     = 1 << 0,
  MASK_QN_SZ    = 1 << 1,
  MASK_QN_KX    = 1 << 2,
  MASK_QN_KY    = 1 << 3,
};

#define MASK_QN_BASIC MASK_QN_N|MASK_QN_SZ
#define MASK_QN_ALL MASK_QN_N|MASK_QN_SZ|MASK_QN_KX|MASK_QN_KY

#ifndef QN_MAX_SIZE
#define QN_MAX_SIZE 2
#endif

class QN
{
  private:
    static Vector<std::string> _qn_name;
    static Vector<int> _qn_norm;
    static int _qn_fermion;
    static int _qn_mask;
    static int _qn_use_bz;

    int _qn0;    
    int _qn1;    
    int _qn2;    
    int _qn3;    
    int _qn4;    
    int _qn5;    

  public:
    static size_t QN_LAST;
    QN(): _qn0(0),_qn1(0),_qn2(0),_qn3(0),_qn4(0),_qn5(0) {}
    QN(int val): _qn0(val),_qn1(val),_qn2(val),_qn3(val),_qn4(val),_qn5(val) {}
    QN(const QN &_qn):
      _qn0(_qn._qn0),_qn1(_qn._qn1),_qn2(_qn._qn2),_qn3(_qn._qn3),_qn4(_qn._qn4),_qn5(_qn._qn5) {}
    // For backward compatibility
    QN(int _n, int _sz) 
      { 
        int idx = get_qn_index("N");
        if(idx != QN_MAX_SIZE) operator[](idx) = _n; 
        idx = get_qn_index("Sz");
        if(idx != QN_MAX_SIZE) operator[](idx) = _sz; 
      }
    QN(int _n, int _sz, int _kx1, int _kx2 = 1) 
      { 
        int idx = get_qn_index("N");
        if(idx != QN_MAX_SIZE) operator[](idx) = _n; 
        idx = get_qn_index("Sz");
        if(idx != QN_MAX_SIZE) operator[](idx) = _sz; 
        idx = get_qn_index("Kx");
        if(idx != QN_MAX_SIZE) {
          operator[](idx) = _kx1; 
          _qn_norm(idx) = _kx2;
        }
      }
//    QN(const dmtk::Vector<int> &v) : dmtk::Vector<int >(v) {}
//    QN(const QN &qn) : dmtk::Vector<int>(qn) {}

//    QN& operator=(const dmtk::Vector<int>&v) { this->operator=(v); return *this; }
    QN& operator=(const QN &_qn) 
      { _qn0 = _qn._qn0; _qn1 = _qn._qn1; _qn2 = _qn._qn2; _qn3 = _qn._qn3; _qn4 = _qn._qn4; _qn5 = _qn._qn5; return *this; }
    QN& operator=(int  val) { _qn0=_qn1=_qn2=_qn3=_qn4=_qn5=val; return *this; }

    int operator[](const std::string& str) const
      { return operator[](get_qn_index(str)); }
    int& operator[](const std::string& str)
      { return operator[](get_qn_index(str)); }

    bool equal(const QN &qn, int mask = QN::default_mask()) const;
    bool operator==(const QN &qn) const;
    bool operator!=(const QN &qn) const;
    bool operator>(const QN &qn) const;
    bool operator<(const QN &qn) const;
    bool operator>=(const QN &qn) const;
    bool operator<=(const QN &qn) const;

    QN& operator+=(const QN &v);
    QN& operator-=(const QN &v);

    int sz() const 
      { 
        size_t idx = get_qn_index("Sz"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("SZ"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("sz"); 
        return (operator[](idx)); 
      }
    int& sz() 
      { 
        size_t idx = get_qn_index("Sz"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("SZ"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("sz"); 
        return (operator[](idx)); 
      }
    int n() const 
      { 
        size_t idx = get_qn_index("N"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("n"); 
        return (operator[](idx)); 
      }
    int& n() 
      { 
        size_t idx = get_qn_index("N"); 
        if(idx == QN_MAX_SIZE) idx = get_qn_index("n"); 
        return (operator[](idx)); 
      }
    int kx1() const 
      { 
        size_t idx = get_qn_index("Kx"); 
        return (operator[](idx)); 
      }
    int& kx1() 
      { 
        size_t idx = get_qn_index("Kx"); 
        return (operator[](idx)); 
      }
    int ky1() const 
      { 
        size_t idx = get_qn_index("Ky"); 
        return (operator[](idx)); 
      }
    int kx2() const 
      { 
        size_t idx = get_qn_index("Kx"); 
        return (_qn_norm(idx)); 
      }
    int& kx2() 
      { 
        size_t idx = get_qn_index("Kx"); 
        return (_qn_norm(idx)); 
      }
    int ky2() const 
      { 
        size_t idx = get_qn_index("Ky"); 
        return (_qn_norm(idx)); 
      }
    int& ky2()  
      { 
        size_t idx = get_qn_index("Ky"); 
        return (_qn_norm(idx)); 
      }

    int fermion_sign() const;
    static size_t max_index() { return QN::QN_LAST; }
    static const std::string& qn_name(size_t i) { return _qn_name[i]; }
    static int qn_norm(size_t i) { return _qn_norm[i]; }

#ifdef USE_QN_NORM        
    QN & normalize() 
      {
        QN &qn = *this;
        for(int i = 0; i < QN_LAST; i++){ 
          int &i1 = qn[i]; 
          int norm = _qn_norm[i];
          if(qn_index_use_bz(i) && norm > 1) {
            i1 = (i1+10000*norm) % norm;
          } else if(norm > 1) {
            int _i1 = std::abs(i1);
            _i1 = i1;
            i1 = _i1%norm; 
            if(i1 < 0) i1 = norm+i1;
          }
        }
        return *this;
      }
#else
    QN & normalize() { return *this; }
#endif

    static size_t get_qn_index(const std::string& str) 
      {
        Vector<std::string>::iterator iter;
        size_t idx = 0;
        for(iter = _qn_name.begin(); iter != _qn_name.end(); iter++, idx++){
          if(*iter == str) return idx;
        }
        return QN_MAX_SIZE;
      };

    static void set_qn_index(size_t idx, const std::string& str)
      { _qn_name(idx) = str; }

    static size_t add_qn_index(const std::string& str, bool _fermion = false, int norm = 1, bool use_bz = false)
      {
         size_t idx = get_qn_index(str);
         if(idx == QN_MAX_SIZE) {
           idx = QN_LAST;
           _qn_name.resize(QN_LAST+1);
           _qn_norm.resize(QN_LAST+1);
           set_qn_index(QN_LAST, str);
           _qn_norm[QN_LAST] = norm;
           if(_fermion) _qn_fermion |= (1 << QN_LAST);
           if(use_bz) _qn_use_bz |= (1 << QN_LAST);
           QN_LAST++;
         }
         return idx;
      } 

    static int qn_index_fermion(size_t idx) { return (IBITS(_qn_fermion,idx)); } 
    static int qn_index_use_bz(size_t idx) { return (IBITS(_qn_use_bz,idx)); } 

    static int default_mask() { return ((1 << QN_LAST)-1); }
    static int get_qn_mask() { return (_qn_mask); }
    static int mask(const std::string &str) { return (1 << get_qn_index(str)); }
    static void init() { _qn_name = std::string(""); _qn_fermion = 0; _qn_use_bz = 0; QN_LAST = 0; }

    static void set_qn_mask(int qn_mask)
    {
      // For backward compatibility with dmtk-1
      if(_qn_mask == 0 && QN_LAST == 0){
        init();
        add_qn_index("N",true);
        add_qn_index("Sz");
        add_qn_index("Kx");
        add_qn_index("Ky");
      }
//      if(qn_mask & MASK_QN_N) add_qn_index("N",true);
//      if(qn_mask & MASK_QN_SZ) add_qn_index("Sz");
//      if(qn_mask & MASK_QN_KX) add_qn_index("Kx");
//      if(qn_mask & MASK_QN_KY) add_qn_index("Ky");
      _qn_mask = qn_mask;
    }
    static int old_qn_mask(int qn_mask)
    {
      int m = 0;
      if(qn_mask & MASK_QN_N) m |= mask("N");
      if(qn_mask & MASK_QN_SZ) m |= mask("Sz");
      if(qn_mask & MASK_QN_KX) m |= mask("Kx");
      if(qn_mask & MASK_QN_KY) m |= mask("Ky");
      return m;
    }

    int operator[](int i) const {
      switch(i){
        case 0: 
          return _qn0;
        case 1: 
          return _qn1;
        case 2: 
          return _qn2;
        case 3: 
          return _qn3;
        case 4: 
          return _qn4;
        case 5: 
          return _qn5;
        default:
          return DMTK_ERROR;
       }
    }
    int & operator[](int i) {
      switch(i){
        case 0: 
          return _qn0;
        case 1: 
          return _qn1;
        case 2: 
          return _qn2;
        case 3: 
          return _qn3;
        case 4: 
          return _qn4;
        case 5: 
          return _qn5;
       }
       return _qn0;
    }


    // Streams

    void write(std::ostream &s) const
    {
      for(int i = 0; i < QN_LAST; i++){
        int n = this->operator[](i);
        s.write((const char *)&n, sizeof(int));
      }
    }

    void read(std::istream &s)
    {
      for(int i = 0; i < QN_LAST; i++){
        int n;
        s.read((char *)&n, sizeof(int));
        operator[](i) = n;
      }
    }

};

size_t QN::QN_LAST = 0;
dmtk::Vector<std::string> QN::_qn_name = dmtk::Vector<std::string>(QN_MAX_SIZE+1);
dmtk::Vector<int> QN::_qn_norm = dmtk::Vector<int>(QN_MAX_SIZE+1);
int QN::_qn_fermion = 0;
int QN::_qn_mask = 0;
int QN::_qn_use_bz = 0;

inline bool
QN::operator==(const QN &qn) const
{
  return equal(qn, QN::get_qn_mask());
}

inline bool
QN::operator!=(const QN &qn) const
{
  return (!equal(qn, QN::get_qn_mask()));
}

inline bool
QN::equal(const QN &qn, int mask) const
{
  QN qn1 = *this; 
  QN qn2 = qn; 
  qn1.normalize(); 
  qn2.normalize(); 

  for(int i = 0; i < QN_LAST; i++){
    if(mask & (1 << i)){
      int i1 = qn1[i];
      int i2 = qn2[i];
      if(i1 != i2) return false;
    }
  }
  return true;
}

//      if(_qn_norm[i] > 1) \
//          { i1 = std::abs(i1%_qn_norm[i]); i2 = std::abs(i2%_qn_norm[i]); } \

#define OP_EXCLUSIVE(op,ap) \
inline bool \
op(const QN &qn) const \
{ \
  QN qn1 = *this; \
  QN qn2 = qn; \
  qn1.normalize(); \
  qn2.normalize(); \
  for(int i = 0; i < QN_LAST; i++){ \
    if(_qn_mask & (1 << i)){ \
      int i1 = qn1[i]; \
      int i2 = qn2[i]; \
      if(i1 == i2)  \
        continue; \
      else if(i1 ap i2)  \
        return true; \
      else; \
        return false; \
    } \
  } \
  \
  return false; \
} 

OP_EXCLUSIVE(QN::operator>,>)
OP_EXCLUSIVE(QN::operator<,<)
#undef OP_EXCLUSIVE

#define OP_INCLUSIVE(op,ap) \
inline bool \
op(const QN &qn) const \
{ \
  QN qn1 = *this; \
  QN qn2 = qn; \
  qn1.normalize(); \
  qn2.normalize(); \
  for(int i = 0; i < QN_LAST; i++){ \
    if(_qn_mask & (1 << i)){ \
      int i1 = qn1[i]; \
      int i2 = qn2[i]; \
      if(i1 == i2)  \
        continue; \
      else if(i1 ap i2)  \
        return true; \
      else; \
        return false; \
    } \
  } \
  \
  return true; \
} 

OP_INCLUSIVE(QN::operator>=,>=)
OP_INCLUSIVE(QN::operator<=,<=)
#undef OP_INCLUSIVE

template <int I, class BinOp> 
struct meta_op{
  static void op(QN& a, const QN& b, const QN &c)
    { a[I] = BinOp::apply(b[I],c[I]); meta_op<I-1,BinOp>::op(a,b,c); }
};

template<class BinOp>
struct meta_op<0,BinOp>{
  static void op(QN& a, const QN& b, const QN &c)
    { a[0] = BinOp::apply(b[0],c[0]); }
};

/*
    if(QN::get_qn_mask() & (1 << i)){  \
      int i1 = qn1[i];  \
      int i2 = qn2[i];  \
      qnres[i] = i1 ap i2; \
*/

#define OP_BIN(op,ap) \
QN \
op(const QN &qn1, const QN &qn2) \
{ \
  QN qnres; \
  for(int i = 0; i < QN::QN_LAST; i++){  \
    int i1 = qn1[i];  \
    int i2 = qn2[i];  \
    qnres[i] = i1 ap i2; \
  }  \
  return qnres; \
} 

OP_BIN(operator+,+)
OP_BIN(operator-,-)
#undef OP_BIN 

/*
    if(_qn_mask & (1 << i)){  \
      qn1[i] ap qn2[i]; \
    } \
*/

#define OP_SELF(op,ap) \
QN& \
op(const QN &qn2) \
{ \
  QN &qn1 = *this;  \
  for(int i = 0; i < QN::QN_LAST; i++){  \
    qn1[i] ap qn2[i]; \
  }  \
  return *this; \
} 

OP_SELF(QN::operator+=,+=)
OP_SELF(QN::operator-=,-=)
#undef OP_SELF

/*  
QN
operator+(const QN &a, const QN &b)
{
  QN c;
  meta_op<QN_MAX_SIZE-1,DMApAdd0<int > >::op(c, a, b);
  return c;
}

QN
operator-(const QN &a, const QN &b)
{
  QN c;
  meta_op<QN_MAX_SIZE-1,DMApSubs0<int > >::op(c, a, b);
  return c;
}

QN&
QN::operator+=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE-1,DMApAdd0<int > >::op(*this,*this,v); 
  return *this; 
}

QN&
QN::operator-=(const QN &v)
{ 
  meta_op<QN_MAX_SIZE-1,DMApSubs0<int > >::op(*this,*this,v); 
  return *this; 
}


template <int I> 
struct meta_and{
  static int op(const QN& a)
    { return (QN::qn_index_fermion(I)*a[I] + meta_and<I-1>::op(a)); }
};

template<>
struct meta_and<0>{
  static int op(const QN& a)
    { return QN::qn_index_fermion(0)*a[0]; }
};


inline int 
QN::fermion_sign() const
{
//cout << "HOLA SIGN " << *this << " " << SGN(meta_and<QN_MAX_SIZE-1>::op(*this)) << endl;
//for(int i = 0; i < QN_MAX_SIZE-1; i++) cout << "HOLA " << i << " " << this->operator[](i) << " " << qn_index_fermion(i) << endl;
  return SGN(meta_and<QN_MAX_SIZE-1>::op(*this));  
}
*/


inline int 
QN::fermion_sign() const
{
  int sign = 1;
  for(int i = 0; i < QN_LAST; i++){
//    if((_qn_mask & (1 << i)) && QN::qn_index_fermion(i) == 1) 
    if(QN::qn_index_fermion(i) == 1) 
      sign *= SGN(this->operator[](i));
  }
//cout << "HOLA FERMION " << *this << sign << endl;
  return sign;
}


std::ostream& operator << (std::ostream& s, const QN& q)
{
  s << "(" << q[0];
  for(int i = 1; i < QN::QN_LAST; i++)
    s << "," << q[i];
  s << ")";

  return s;
}

} // namespace dmtk

#endif // __DMTK_QN_H__
