#ifndef __DMTK_COUNT_REF_H__
#define __DMTK_COUNT_REF_H__

#include <stdio.h>
#include <iostream>
#include <iosfwd>

using namespace std;

namespace dmtk
{

#define ASSERT(a) { if(!a){ cerr << "*** ERROR !a \n" << endl; exit(0); } }

template <class T> class CountedRef;

template<class T>
class Counted
{
    friend class CountedRef<T>;

    static int _tot_counted;

  private:
    Counted (T* pt): _pt(pt), _count(0) { ASSERT(_pt != 0); _tot_counted++; }

    ~Counted() { ASSERT(_count == 0); delete _pt; _tot_counted--; }

    size_t get_ref() { return ++_count; }
    size_t unref() { ASSERT(_count != 0); return --_count; }
    size_t ref_count() const { return _count; }

    T* _pt;
    size_t _count;

  public:
    static int tot_counted() { return _tot_counted; } 
};
  
template<class T>
class CountedRef
{
  private:
    static int _tot_counted;

    Counted<T>* _counted;
    
    void unbind()
      {
        if(!null() && _counted->unref() == 0) { delete _counted; }
        _counted = NULL;
      }
    
  public:
    CountedRef(): _counted(NULL) { _tot_counted++; }
    CountedRef(T* pt) { _counted = new Counted<T>(pt); _counted->get_ref(); _tot_counted++; }
    CountedRef(const T& pt) { T* new_counted = new T(pt); _counted = new Counted<T>(new_counted); _counted->get_ref(); _tot_counted++; }

    ~CountedRef() { unbind(); _tot_counted--; }

    CountedRef(const CountedRef<T>& other)
      {
        _counted = other._counted;
        if(!null()) _counted->get_ref();
        _tot_counted++;
      }

    CountedRef& operator=(const CountedRef<T>& other)
      {
        if(!other.null()) other._counted->get_ref();

        unbind();
        _counted = other._counted;
        return *this;
      }

    CountedRef& operator=(T* other)
      {
        unbind();
        _counted = new Counted<T>(other); 
        _counted->get_ref(); 
        return *this;
      }

    CountedRef& operator=(const T* other)
      {
        unbind();
        if(null()){
          T* new_data = new T(*other);
          _counted->_pt = new_data;
          _counted->get_ref();
        } else {
          _counted->_pt->operator=(*other);
        }
        return *this;
      }

    CountedRef& operator<<=(const T& other)
      {
        if(null()){
          T* new_data = new T(other);
          _counted->_pt = new_data;
          _counted->get_ref();
        } else {
          _counted->_pt->operator=(other);
        }
        return *this;
      }

    CountedRef& operator<<=(const CountedRef<T>& other)
      {
        if(null()){
          T* new_data = new T(*(other._counted->_pt));
          _counted->_pt = new_data;
          _counted->get_ref();
        } else {
          _counted->_pt->operator=(*(other._counted->_pt));
        }
        return *this;
      }

    T* operator->() 
      { 
        if(null()) { cerr << "*** ERROR: NULL reference\n"; exit(0); }
        return _counted->_pt;
      }

    const T* operator->() const 
      { 
        if(null()) { cerr << "*** ERROR: NULL reference\n"; exit(0); }
        return _counted->_pt;
      }

    T* ref() 
      { 
        if(null()) { cerr << "*** ERROR: NULL reference\n"; exit(0); }
        return _counted->_pt;
      }

    const T* ref() const 
      { 
        if(null()) { cerr << "*** ERROR: NULL reference\n"; exit(0); }
        return _counted->_pt;
      }

    bool operator==(const CountedRef<T>& other)
      { return _counted == other._counted; }

    bool null() const { return _counted == 0; } 
    size_t ref_count() const 
      { 
        if(!null()) 
          return _counted->ref_count(); 
        else
          return 0;
      }

    void  set_null() { unbind(); } 

    static int tot_counted() { return _tot_counted; } 
};

template<class T>
int CountedRef<T>::_tot_counted = 0;

template<class T>
int Counted<T>::_tot_counted = 0;

} // namespace dmtk

#endif // __DMTK_COUNT_REF_H__
