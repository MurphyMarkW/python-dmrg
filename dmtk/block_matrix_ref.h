#ifndef __DMTK_BLOCK_MATRIX_REF_H__
#define __DMTK_BLOCK_MATRIX_REF_H__

#include "counted.h"
#include "block_matrix.h"

namespace dmtk
{

template <class T>
class BMatrixRef: public CountedRef<BMatrix<T> >
{
  public:
    BMatrixRef(): {}
    BMatrixRef(const BMatrix<T>* m): CountedRef<BMatrix<T> >(m) {}
    BMatrixRef(const BMatrixRef<T>& other ): CountedRef<BMatrix<T> >(other) {}

    BMatrixRef& repack(const Basis &basis)
      { return _counted->_pt->repack(basis); }

    BMatrixRef& operator=(const T& v)
      { return _counted->_pt->operator=(v); }

    const SubMatrix<T>* block(const QN &qn) const
      { return _counted->_pt->block(qn); }

    SubMatrix<T>* block(const QN &qn) 
      { return _counted->_pt->block(qn); }

    SubSpace subspace(const QN &qn) const
      { return _counted->_pt->subspace(qn); }

    SubSpace subspace(size_t index) const
      { return _counted->_pt->subspace(index); }

    const PackedBasis& subspaces() const 
      { return _counted->_pt->subspaces(); }

    T operator() (size_t col, size_t row) const
      { return _counted->_pt->operator()(col,row); }

    PackedBasis::iterator subspace_begin() 
      { return _counted->_pt->subspace_begin(); }
    PackedBasis::const_iterator subspace_begin() const 
      { return _counted->_pt->subspace_begin(); }
    PackedBasis::iterator subspace_end() 
      { return _counted->_pt->subspace_end(); }
    PackedBasis::const_iterator subspace_end() const
      { return _counted->_pt->subspace_end(); }

    BMatrixRef& diagonalize(Vector<double>& ev)
      { _counted->_pt->diagonalize(ev); } 

    size_t dim() const
      { return (_counted->_pt->dim()); }

    BMatrixRef& resize(const Basis& b) 
      { _counted->_pt->resize(b); return *this; }
    BMatrixRef& resize(const Basis& b1, const Basis &b2) 
      { _counted->_pt->resize(b1,b2); return *this; }

//  Streams

    void read(istream& s)
      { _counted->_pt->read(s); } 

    void write(ostream& s) const
      { _counted->_pt->write(s); } 

};

} // namespace dmtk

#endif // __DMTK_BLOCK_MATRIX_REF_H__
