#ifndef __DMTK_QN_H__
#define __DMTK_QN_H__

#include <iosfwd>
#include "constants.h"
#include "bits.h"

namespace dmtk 
{


// NOTE: Although the k code has been implemented, it is not actually been used
// This is because the operators mix subspaces with different k's and I haven't
// figured out a way to deal with it, yet

enum
{
  QN_N,
  QN_SZ,
  QN_NB,
  QN_KX1,
  QN_KY1,
  QN_KX2,
  QN_KY2,
};

enum
{
  MASK_QN_N 	= 1 << 0,
  MASK_QN_SZ 	= 1 << 1,
  MASK_QN_NB 	= 1 << 2,
  MASK_QN_KX 	= 1 << 3,
  MASK_QN_KY 	= 1 << 4,
};

#define MASK_QN_BASIC MASK_QN_N|MASK_QN_SZ
#define MASK_QN_ALL MASK_QN_N|MASK_QN_NB|MASK_QN_SZ|MASK_QN_KX|MASK_QN_KY

#ifdef USE_K
#define QN_LAST 5
#elif USE_NB
#define QN_LAST 3
#else
#define QN_LAST 2
#endif

class QN
{
  private:
    int _sz;
    int _n;
    int _nb;
    int _kx1, _kx2; // numerator and denominator of 2*PI*kx1/kx2
    int _ky1, _ky2; // numerator and denominator of 2*PI*ky1/ky2

    static int _qn_mask;

    QN& bz()
      {
#ifdef USE_K
#ifdef USE_BZ
        _kx1 = (_kx1+10000*_kx2) % _kx2;
        _ky1 = (_ky1+10000*_ky2) % _ky2;
#endif // USE_BZ
#endif // USE_K
        return *this;
      }
    void equalize(QN &qn1, QN&qn2) const
      {

#ifdef GRANDCANONICAL
        qn1._n = mod2(qn1._n); //*(qn1._n < 0 ? -1 : 1);
        qn2._n = mod2(qn2._n); //*(qn2._n < 0 ? -1 : 1);
//        qn1._sz = mod2(qn1._sz); //*(qn1._n < 0 ? -1 : 1);
//        qn2._sz = mod2(qn2._sz); //*(qn2._n < 0 ? -1 : 1);
#endif
#ifdef GRANDCANONICAL_N
        qn1._n = 0;
        qn2._n = 0;
#endif
#ifdef GRANDCANONICAL_SZ
        qn1._sz = 0;
        qn2._sz = 0;
#endif
#ifdef USE_K
        if(qn1._kx2 != qn2._kx2){
          qn1._kx1 *= qn2._kx2;
          qn2._kx1 *= qn1._kx2;
          qn1._kx2 = qn1._kx2 * qn2._kx2;
          qn2._kx2 = qn1._kx2;
        }  
        if(qn1._ky2 != qn2._ky2){
          qn1._ky1 *= qn2._ky2;
          qn2._ky1 *= qn1._ky2;
          qn1._ky2 = qn1._ky2 * qn2._ky2;
          qn2._ky2 = qn1._ky2;
        }  
#ifdef USE_BZ
        qn1.bz();
        qn2.bz();
#endif // USE_BZ
#endif // USE_K
      } 

  public:
    QN(): _n(0), _sz(0), _nb(0), _kx1(0), _kx2(1), _ky1(0), _ky2(1) {}
    explicit QN(int v): _sz(v), _n(v), _nb(0), _kx1(0), _kx2(1), _ky1(0), _ky2(1) {} 
    QN(int n, int sz, int kx1 = 0, int kx2 = 1, int ky1 = 0, int ky2 = 1): _sz(sz), _n(n), _nb(0), _kx1(kx1), _kx2(kx2), _ky1(ky1), _ky2(ky2) {} 
    QN(const QN& qn): _sz(qn._sz), _n(qn._n), _nb(qn._nb), _kx1(qn._kx1), _kx2(qn._kx2), _ky1(qn._ky1), _ky2(qn._ky2) { bz(); } 

    QN& operator=(int v) 
      { _sz = v; _n = v; return *this; }
    QN& operator=(const QN& qn) 
      { _sz = qn._sz; _n = qn._n; _nb = qn._nb; _kx1 = qn._kx1; _kx2 = qn._kx2; _ky1 = qn._ky1, _ky2 = qn._ky2 ; bz(); return *this; }

    int operator[](int t) const 
      { 
        switch(t){
          case QN_SZ:
            return _sz;
          case QN_N:
            return _n;
          case QN_NB:
            return _nb;
          case QN_KX1:
            return _kx1;
          case QN_KX2:
            return _kx2;
          case QN_KY1:
            return _ky1;
          case QN_KY2:
            return _ky2;
          default:
            cout << "*** WARNING QN: index out of bounds\n";
            return 0;
        }
      }
    int& operator[](int t) 
      {
        switch(t){
          case QN_SZ:
            return _sz;
          case QN_N:
            return _n;
          case QN_NB:
            return _nb;
          case QN_KX1:
            return _kx1;
          case QN_KX2:
            return _kx2;
          case QN_KY1:
            return _ky1;
          case QN_KY2:
            return _ky2;
          default:
            cout << "*** WARNING QN: index out of bounds\n";
            return _sz;
        }
      }

    QN operator+(const QN&qn) const 
      {
        QN aux;
        aux._sz = _sz + qn._sz;
        aux._n = _n + qn._n;
#ifdef USE_NB 
        aux._nb = _nb + qn._nb;
#endif
#ifdef USE_K
        if(_kx2 == qn._kx2){
          aux._kx1 = _kx1 + qn._kx1;
          aux._kx2 = _kx2;
        } else {
          aux._kx1 = _kx1*qn._kx2 + qn._kx1*_kx2;
          aux._kx2 = qn._kx2 * _kx2;
        }  
        if(_ky2 == qn._ky2){
          aux._ky1 = _ky1 + qn._ky1;
          aux._ky2 = _ky2;
        } else {
          aux._ky1 = _ky1*qn._ky2 + qn._ky1*_ky2;
          aux._ky2 = qn._ky2 * _ky2;
        }  
        aux.bz();
#endif // USE_K
        return aux;
      }
    QN operator-(const QN&qn) const 
      {
        QN aux;
        aux._sz = _sz - qn._sz;
        aux._n = _n - qn._n;
#ifdef USE_NB 
        aux._nb = _nb - qn._nb;
#endif
#ifdef USE_K
        if(_kx2 == qn._kx2){
          aux._kx1 = _kx1 - qn._kx1;
          aux._kx2 = _kx2;
        } else {
          aux._kx1 = _kx1*qn._kx2 - qn._kx1*_kx2;
          aux._kx2 = qn._kx2 * _kx2;
        }  
        if(_ky2 == qn._ky2){
          aux._ky1 = _ky1 - qn._ky1;
          aux._ky2 = _ky2;
        } else {
          aux._ky1 = _ky1*qn._ky2 - qn._ky1*_ky2;
          aux._ky2 = qn._ky2 * _ky2;
        }  
        aux.bz();
#endif // USE_K
        return aux;
      }
    QN& operator+=(const QN&qn) 
      { 
        _sz += qn._sz; _n += qn._n; 
#ifdef USE_NB
        _nb += qn._nb;
#endif
#ifdef USE_K
        if(_kx2 == qn._kx2){
          _kx1 = _kx1 + qn._kx1;
          _kx2 = _kx2;
        } else {
          _kx1 = _kx1*qn._kx2 + qn._kx1*_kx2;
          _kx2 = qn._kx2 * _kx2;
        }  
        if(_ky2 == qn._ky2){
          _ky1 = _ky1 + qn._ky1;
          _ky2 = _ky2;
        } else {
          _ky1 = _ky1*qn._ky2 + qn._ky1*_ky2;
          _ky2 = qn._ky2 * _ky2;
        }  
        bz();
#endif // USE_K
        return *this;
      }
    QN& operator-=(const QN&qn)  
      { 
        _sz -= qn._sz; _n -= qn._n; 
#ifdef USE_NB
        _nb -= qn._nb;
#endif
#ifdef USE_K
        if(_kx2 == qn._kx2){
          _kx1 = _kx1 - qn._kx1;
          _kx2 = _kx2;
        } else {
          _kx1 = _kx1*qn._kx2 - qn._kx1*_kx2;
          _kx2 = qn._kx2 * _kx2;
        }  
        if(_ky2 == qn._ky2){
          _ky1 = _ky1 - qn._ky1;
          _ky2 = _ky2;
        } else {
          _ky1 = _ky1*qn._ky2 - qn._ky1*_ky2;
          _ky2 = qn._ky2 * _ky2;
        }  
        bz();
#endif // USE_K
        return *this;
      }

    bool operator==(const QN &qn) const;
    bool operator!=(const QN &qn) const;
    bool operator>(const QN &qn) const;
    bool operator<(const QN &qn) const;
    bool operator>=(const QN &qn) const;
    bool operator<=(const QN &qn) const;

    bool equal(const QN &qn, int mask = MASK_QN_BASIC) const;

    int sz() const { return _sz; }
    int& sz() { return _sz; }
    int n() const { return _n; }
    int& n() { return _n; }
    int nb() const { return _nb; }
    int& nb() { return _nb; }
    int kx1() const { return _kx1; }
    int& kx1() { return _kx1; }
    int kx2() const { return _kx2; }
    int& kx2() { return _kx2; }
    int ky1() const { return _ky1; }
    int& ky1() { return _ky1; }
    int ky2() const { return _ky2; }
    int& ky2() { return _ky2; }
    double kx() const { return 2*double(_kx1)/double(_kx2); }
    double ky() const { return 2*double(_ky1)/double(_ky2); }
    double vkx() const { return DMTK_TWOPI*double(_kx1)/double(_kx2); }
    double vky() const { return DMTK_TWOPI*double(_ky1)/double(_ky2); }
    static void set_qn_mask(int qn_mask) { _qn_mask = qn_mask; }
    static int get_qn_mask() { return _qn_mask; }

    // Streams

    void write(std::ostream &s) const
    {
      for(int i = 0; i < 7; i++){
        int n = operator[](i);
        s.write((const char *)&n, sizeof(int));
      }
    }

    void read(std::istream &s)
    {
      for(int i = 0; i < 7; i++){
        s.read((char *)&operator[](i), sizeof(int));
      }
    }

};

int QN::_qn_mask = MASK_QN_BASIC;

inline bool
QN::operator==(const QN &qn) const
{
  return equal(qn, _qn_mask);
}

inline bool
QN::operator!=(const QN &qn) const
{
  QN qn1 = *this;
  QN qn2 = qn;
  equalize(qn1,qn2);
  for(int i = 0; i < QN_LAST; i++){
    if((_qn_mask & (1 << i)) && qn1[i] != qn2[i]) return true;
  }
  return false;
}

inline bool
QN::equal(const QN &qn, int mask) const
{
  QN qn1 = *this;
  QN qn2 = qn;
  equalize(qn1,qn2);
  for(int i = 0; i < QN_LAST; i++){
    if((mask & (1 << i)) && qn1[i] != qn2[i]) return false;
  }
  return true;
}

#define OP_EXCLUSIVE(op,ap) \
inline bool \
op(const QN &qn) const \
{ \
  QN qn1 = *this; \
  QN qn2 = qn; \
  equalize(qn1,qn2); \
  for(int i = 0; i < QN_LAST; i++){ \
    if((_qn_mask & (1 << i))) { \
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
  equalize(qn1,qn2); \
  for(int i = 0; i < QN_LAST; i++){ \
    if((_qn_mask & (1 << i))) { \
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

} // namespace dmtk

#endif // __DMTK_QN_H__
