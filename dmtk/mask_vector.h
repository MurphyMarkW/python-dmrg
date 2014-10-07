#ifndef __DMTK_MASK_VECTOR_H__
#define __DMTK_MASK_VECTOR_H__

namespace dmtk 
{

template <class T>
class MVector
{
  private:
    dmtk::Vector<T> _v;
    dmtk::Vector<bool> _mask;

    T& ref(size_t i)
      {
        Vector<T>::iterator iter_v;
        Vector<bool>::iterator iter_mask;

        size_t j;
        for(iter_mask = _mask.begin(); iter_mask != _mask.end() && j < i; iter_mask++, iter_v++)
          if(*iter_mask) j++;

        return *iter_v;
      }
    T ref(size_t i) const
      {
        Vector<T>::const_iterator iter_v;
        Vector<bool>::const_iterator iter_mask;

        size_t j;
        for(iter_mask = _mask.begin(); iter_mask != _mask.end() && j < i; iter_mask++, iter_v++)
          if(*iter_mask) j++;

        return *iter_v;
      }
       
  public:
    MVector(const Vector<T> &v):
      _v(v.v), _mask(Vector<bool>(v.size(),true)) {}
    MVector(const MVector<T> &v):
      _v(v.v), _mask(v.mask) {}
    MVector(const Vector<T> &v, const Vector<bool> &mask):
      _v(v), _mask(mask) {}

    T& operator()(size_t i) {return ref(i);}
    T operator()(size_t i) const {return ref(i);}
    T& operator[](size_t i) {return ref(i);}
    T operator[](size_t i) const {return ref(i);}
};


} // namsepace dmtk

#endif // __DMTK_MASK_VECTOR_H__
