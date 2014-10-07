#ifndef __DMTK_HMATRIX__
#define __DMTK_HMATRIX__

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
#include "system.h"

//namespace dmtk
//{

template<class T>
void
hmatrix_write_column(const BMatrix<T> &column, size_t size1, size_t size2, ofstream &file)
{
  Vector<T> data(size2);
  Vector<size_t> index(size2);

  typename BMatrix<T>::const_iterator biter;
  for(int i = 0; i < size1; i++){
    int ncol = 1;
    data(0) = T(0);
    for(biter = column.begin(); biter != column.end(); biter++){
      const SubMatrix<T> &block = (*biter);
      int jstart = 0;
      if(block.row_range().qn().n() == block.col_range().qn().n()){
        jstart = i+1;
        data(0) = block(i,i);
      }
      for(int j = jstart; j < block.rows(); j++){
        T val = block(i,j);
        if(fabs(val) > 1.e-10) {
          index(ncol) = j+block.row_range().begin(); 
          data(ncol) = val;
          ncol++;
        }
      }      
    }
    index(0) = ncol;
    data.resize(ncol);
    index.resize(ncol);
    index.write(file);
    data.write(file);
  }
}

//////////////////////////////////////////////////////////////////////
template<class T>
void
hmatrix_create(const System<T>& s)
{
  CTimer clock;
  clock.Start();
  int total_size = 0;

// Saving file
  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", s.name());
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
//////////////////////////////////////////////////////////////////////

  VectorState<T> v(s.seed);

  int index_i = 0;
  BMatrix<T> column;
  vector<bool> used(v.subspaces().size());
  typename VectorState<T>::const_iterator siter, siter2;
  for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++, index_i++){
    StateSpace ss = *siter;

//    int index_i = v.get_qn_space_index(ss[1].qn(),ss[2].qn(),ss[3].qn(),ss[4].qn());
/*
 cout << "HOLA COLUMN " << index_i << endl;
cout << "HOLA QN " << ss[1].qn().n() << " " << ss[2].qn().n() << " " << ss[3].qn().n() << " " << ss[4].qn().n() << endl;
cout << "HOLA QN " << ss[1].qn().kx() << " " << ss[2].qn().kx() << " " << ss[3].qn().kx() << " " << ss[4].qn().kx() << endl;
QN qn = ss[1].qn() + ss[2].qn() + ss[3].qn() + ss[4].qn();
cout << "HOLA QNT " << qn.n() << " " << qn.kx() << endl;
*/
// We build a block matrix for the column
    column.clear();
    int istart = index_i;
    siter2 = siter;
    if(s.h.use_hc()){
      istart = 0;
      siter2 = v.subspace_begin();
    }
    for(; siter2 != v.subspace_end(); siter2++){
      StateSpace ss2 = *siter2;
      int index_j = v.get_qn_space_index(ss2[1].qn(),ss2[2].qn(),ss2[3].qn(),ss2[4].qn());
      SubSpace col_range(QN(index_i),ss.start(),ss.start()+ss.dim()-1);
      SubSpace row_range(QN(index_j),ss2.start(),ss2.start()+ss2.dim()-1);
// cout << "HOLA RANGE " << index_j << " " << ss.dim() << " " << ss2.dim() << endl;
      SubMatrix<T> aux(QN(index_j),col_range,row_range);
      column.push_back(aux);
// cout << "HOLA ROW " << index_j << " " << col_range.size() << " " << row_range.size() << endl;
    }
    used.resize(column.size());
    for(int i = 0; i < used.size(); i++) used[i] = false;
///////////////////

    clock.Lap();
    typename std::list<AuxTerm<T> >::const_iterator titer;
    for(titer = s.aux_terms.begin(); titer != s.aux_terms.end(); titer++){
      const AuxTerm<T> &auxt = (*titer);
      Vector<const BasicOp<T>*> ops(5);
      Vector<int> bindex(5);
      int nops = 0;

      Term<T> t = (auxt.t);
      T coef = s.use_composite() ? T(1) : T(t.coef());

/*
      bool calc_hc = s.h.use_hc();
      BasicOp<T> _op(auxt.t);
      if(_op.is_diagonal()) calc_hc = false;
*/

      if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

        if(s.use_composite()){
          if(auxt.b1 < auxt.b2){ 
            ops[0] = auxt.ref_op[0];
            ops[1] = &auxt.sum_op;
            bindex[0] = auxt.b1;
            bindex[1] = auxt.b2;
          } else {
            ops[1] = auxt.ref_op[0];
            ops[0] = &auxt.sum_op;
            bindex[0] = auxt.b2;
            bindex[1] = auxt.b1;
          }
        } else {
          if(auxt.b1 < auxt.b2){ 
            ops[0] = auxt.op[size_t(auxt.b1)-1];
            ops[1] = auxt.op[size_t(auxt.b2)-1];
            bindex[0] = auxt.b1;
            bindex[1] = auxt.b2;
          } else {
            ops[0] = auxt.op[size_t(auxt.b2)-1];
            ops[1] = auxt.op[size_t(auxt.b1)-1];
            bindex[0] = auxt.b2;
            bindex[1] = auxt.b1;
          }
        }
        nops = 2;
      }
      else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        ops[0] = auxt.op[0];
        ops[1] = auxt.op[1];
        ops[2] = auxt.op[2];
        ops[3] = auxt.op[3];
        bindex[0] = BLOCK1;
        bindex[1] = BLOCK2;
        bindex[2] = BLOCK3;
        bindex[3] = BLOCK4;
        nops = 4;
      }
      else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
        if(s.use_composite()){
          ops[0] = &auxt.sum_op;
          ops[1] = auxt.ref_op[0];
          ops[2] = auxt.ref_op[1];
        } else {
          ops[0] = auxt.op[0];
          ops[1] = auxt.op[1];
          ops[2] = auxt.op[2];
        }
        bindex[0] = BLOCK1;
        bindex[1] = BLOCK2;
        bindex[2] = BLOCK3;
        nops = 3;
      }
      else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
        if(s.use_composite()){
          if(auxt.b2 < auxt.b3){
            ops[0] = &auxt.sum_op;
            ops[1] = auxt.ref_op[0];
            ops[2] = auxt.ref_op[1];
          } else {
            ops[0] = auxt.ref_op[1];
            ops[1] = auxt.ref_op[0];
            ops[2] = &auxt.sum_op;
          }
        } else {
          ops[0] = auxt.op[0];
          ops[1] = auxt.op[1];
          ops[2] = auxt.op[3];
        }
        bindex[0] = BLOCK1;
        bindex[1] = BLOCK2;
        bindex[2] = BLOCK4;
        nops = 3;
      }
      else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
        if(s.use_composite()){
          if(auxt.b2 < auxt.b3){
            ops[0] = &auxt.sum_op;
            ops[1] = auxt.ref_op[0];
            ops[2] = auxt.ref_op[1];
          } else {
            ops[0] = auxt.ref_op[1];
            ops[1] = auxt.ref_op[0];
            ops[2] = &auxt.sum_op;
          }
        } else {
          ops[0] = auxt.op[0];
          ops[1] = auxt.op[2];
          ops[2] = auxt.op[3];
        }
        bindex[0] = BLOCK1;
        bindex[1] = BLOCK3;
        bindex[2] = BLOCK4;
        nops = 3;
      }
      else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
        if(s.use_composite()){
          ops[0] = auxt.ref_op[0];
          ops[1] = auxt.ref_op[1];
          ops[2] = &auxt.sum_op;
        } else {
          ops[0] = auxt.op[1];
          ops[1] = auxt.op[2];
          ops[2] = auxt.op[3];
        }
        bindex[0] = BLOCK2;
        bindex[1] = BLOCK3;
        bindex[2] = BLOCK4;
        nops = 3;
      }
      int sign = 1;
      Vector<QN> dqn(5);
      dqn = QN(0);
      for(int iop = nops-1; iop >=0; iop--){
        if(ops[iop]->fermion()){
          int saux = 0;
          for(int ib = bindex[iop]-1; ib >= 1; ib--) saux += ss[ib].qn().n();
          sign *= SGN(saux);
        }
        dqn[bindex[iop]] = ops[iop]->dqn;
      }
      state_slice<T> v_slice = v(ss);
      int index_j = v.get_qn_space_index(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                         ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
      if(index_j == -1) continue;
      StateSpace ss2 = v.subspaces()[index_j];
      state_slice<T> res_slice = v(ss2);
      if(res_slice.size() == 0) continue;

      if((!s.use_hc() && index_j >= index_i) || s.use_hc()){
        Matrix<T> m(v_slice.size(),res_slice.size());
        m = T(0);
        const SubMatrix<T> *_sm1, *_sm2, *_sm3, *_sm4;
        if(ops[0]) _sm1 = ops[0]->block(ss[bindex[0]].qn());
        if(ops[1]) _sm2 = ops[1]->block(ss[bindex[1]].qn());
        if(ops[2]) _sm3 = ops[2]->block(ss[bindex[2]].qn());
        if(ops[3]) _sm4 = ops[3]->block(ss[bindex[3]].qn());

        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i2 = 0; i2 < ss[2].dim(); i2++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++){
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,j2,i3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2);
                    }
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i3 = 0; i3 < ss[3].dim(); i3++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,i2,j3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i3,j3);
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,i2,i3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i2 = 0; i2 < ss[2].dim(); i2++)
                for(int i3 = 0; i3 < ss[3].dim(); i3++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,j2,j3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i3,j3);
        }
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK4) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i2 = 0; i2 < ss[2].dim(); i2++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,j2,i3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,i2,j3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i3,j3)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3 && _sm4){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j3 = 0; j3 < res_slice.size3(); j3++)
                        for(int j4 = 0; j4 < res_slice.size4(); j4++){
                          m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                            res_slice.index(j1,j2,j3,j4)-res_slice.start1()) = 
                             t.coef()*T(sign)*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i3,j3)*_sm4->operator()(i4,j4);
                        }
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3) && _sm1 && _sm2 && _sm3){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j3 = 0; j3 < res_slice.size3(); j3++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,j2,j3,i4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i3,j3);
        } 
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,j2,i3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i4,j4);
        } 
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,i2,j3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i3,j3)*_sm3->operator()(i4,j4);
        } 
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          used[index_j-istart] = true;
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(i1,j2,j3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i3,j3)*_sm3->operator()(i4,j4);
        } 
        SubMatrix<T> *_block = column[index_j-istart];
        if(!_block) cout << "*** ERROR: hmatrix_create \n";
        SubMatrix<T> &block = *_block;
        block += m;
      }

    }

    cout << "Hmatrix column Lap: " << index_i << " " << clock.LapTime() << endl;
    clock.Lap();
  
    typename BMatrix<T>::iterator biter;
    vector<bool>::iterator uiter;
    uiter = used.begin(); biter = column.begin();
    while(biter != column.end()){
      if(!*uiter) {
        biter = column.erase(biter);
        uiter = used.erase(uiter);
      } else {
        SubMatrix<T> &block = (*biter);

//cout << "HMATRIX BLOCK " << block.col_range().qn().n() << " " << block.row_range().qn().n() << endl;
//cout << "HMATRIX BLOCK " << block.cols() << " " << block.rows() << endl;
/*
for(int i = 0; i < block.cols(); i++)
for(int j = 0; j < block.rows(); j++) cout << "HOLA COL ROW " << i << " " << j << " " << block(i,j) << endl;
*/
        total_size += block.cols()*block.rows();
        biter++; uiter++;
      }
    }
    cout << "Hmatrix column cleanup Lap: " << index_i << " " << clock.LapTime() << endl;
    clock.Lap();
#ifdef USE_SPARSE
    hmatrix_write_column(column,ss.dim(),v.size(),outputfile);
#else
    column.write(outputfile);
#endif
    cout << "Hmatrix column write Lap: " << index_i << " " << clock.LapTime() << endl;
    clock.Lap();
  }
  cout << "HMATRIX SIZE " << total_size << " " << v.size()*v.size() << " " << double(total_size)/double(v.size()*v.size()) << endl;
  cout << "Hmatrix time: " << clock.TotalTime() << endl;

  outputfile.close();
}


//////////////////////////////////////////////////////////////////////
template<class T>
void
hmatrix_create_new(const System<T>& s)
{
  CTimer clock;
  clock.Start();
  int total_size = 0;

// Saving file
  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", s.name());
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
//////////////////////////////////////////////////////////////////////

  VectorState<T> v(s.seed);

///////////////////

  BMatrix<T> column;
  typename std::list<AuxTerm<T> >::const_iterator titer;
  for(titer = s.aux_terms.begin(); titer != s.aux_terms.end(); titer++){
    const AuxTerm<T> &auxt = (*titer);
    Vector<const BasicOp<T>*> ops(5);
    Vector<int> bindex(5);
    int nops = 0;

    clock.Lap();
    column.clear();
    Term<T> t = (auxt.t);
    T coef = s.use_composite() ? T(1) : T(t.coef());

    if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) || auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) || auxt.mask == (MASK_BLOCK2|MASK_BLOCK4)){

      if(s.use_composite()){
        if(auxt.b1 < auxt.b2){ 
          ops[0] = auxt.ref_op[0];
          ops[1] = &auxt.sum_op;
          bindex[0] = auxt.b1;
          bindex[1] = auxt.b2;
        } else {
          ops[1] = auxt.ref_op[0];
          ops[0] = &auxt.sum_op;
          bindex[0] = auxt.b2;
          bindex[1] = auxt.b1;
        }
      } else {
        if(auxt.b1 < auxt.b2){ 
          ops[0] = auxt.op[size_t(auxt.b1)-1];
          ops[1] = auxt.op[size_t(auxt.b2)-1];
          bindex[0] = auxt.b1;
          bindex[1] = auxt.b2;
        } else {
          ops[0] = auxt.op[size_t(auxt.b2)-1];
          ops[1] = auxt.op[size_t(auxt.b1)-1];
          bindex[0] = auxt.b2;
          bindex[1] = auxt.b1;
        }
      }
      nops = 2;
    }
    else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
      ops[0] = auxt.op[0];
      ops[1] = auxt.op[1];
      ops[2] = auxt.op[2];
      ops[3] = auxt.op[3];
      bindex[0] = BLOCK1;
      bindex[1] = BLOCK2;
      bindex[2] = BLOCK3;
      bindex[3] = BLOCK4;
      nops = 4;
    }
    else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3)){
      if(s.use_composite()){
        ops[0] = &auxt.sum_op;
        ops[1] = auxt.ref_op[0];
        ops[2] = auxt.ref_op[1];
      } else {
        ops[0] = auxt.op[0];
        ops[1] = auxt.op[1];
        ops[2] = auxt.op[2];
      }
      bindex[0] = BLOCK1;
      bindex[1] = BLOCK2;
      bindex[2] = BLOCK3;
      nops = 3;
    }
    else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4)){
      if(s.use_composite()){
        if(auxt.b2 < auxt.b3){
          ops[0] = &auxt.sum_op;
          ops[1] = auxt.ref_op[0];
          ops[2] = auxt.ref_op[1];
        } else {
          ops[0] = auxt.ref_op[1];
          ops[1] = auxt.ref_op[0];
          ops[2] = &auxt.sum_op;
        }
      } else {
        ops[0] = auxt.op[0];
        ops[1] = auxt.op[1];
        ops[2] = auxt.op[3];
      }
      bindex[0] = BLOCK1;
      bindex[1] = BLOCK2;
      bindex[2] = BLOCK4;
      nops = 3;
    }
    else if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4)){
      if(s.use_composite()){
        if(auxt.b2 < auxt.b3){
          ops[0] = &auxt.sum_op;
          ops[1] = auxt.ref_op[0];
          ops[2] = auxt.ref_op[1];
        } else {
          ops[0] = auxt.ref_op[1];
          ops[1] = auxt.ref_op[0];
          ops[2] = &auxt.sum_op;
        }
      } else {
        ops[0] = auxt.op[0];
        ops[1] = auxt.op[2];
        ops[2] = auxt.op[3];
      }
      bindex[0] = BLOCK1;
      bindex[1] = BLOCK3;
      bindex[2] = BLOCK4;
      nops = 3;
    }
    else if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4)){
      if(s.use_composite()){
        ops[0] = auxt.ref_op[0];
        ops[1] = auxt.ref_op[1];
        ops[2] = &auxt.sum_op;
      } else {
        ops[0] = auxt.op[1];
        ops[1] = auxt.op[2];
        ops[2] = auxt.op[3];
      }
      bindex[0] = BLOCK2;
      bindex[1] = BLOCK3;
      bindex[2] = BLOCK4;
      nops = 3;
    }

    Vector<QN> dqn(5);
    dqn = QN(0);
    for(int iop = nops-1; iop >=0; iop--){
      dqn[bindex[iop]] = ops[iop]->dqn;
    }

    typename VectorState<T>::const_iterator siter, siter2;
    int index_i = 0;
    for(siter = v.subspace_begin(); siter != v.subspace_end(); siter++, index_i++){
      StateSpace ss = *siter;

      state_slice<T> v_slice = v(ss);
      int index_j = v.get_qn_space_index(ss[1].qn()+dqn[1],ss[2].qn()+dqn[2],
                                         ss[3].qn()+dqn[3],ss[4].qn()+dqn[4]);
      if(index_j == -1) continue;
      StateSpace ss2 = v.subspaces()[index_j];
      state_slice<T> res_slice = v(ss2);
      if(res_slice.size() == 0) continue;
      SubSpace col_range(QN(index_i),ss.start(),ss.start()+ss.dim()-1);
      SubSpace row_range(QN(index_j),ss2.start(),ss2.start()+ss2.dim()-1);

      if((!s.use_hc() && index_j >= index_i) || s.use_hc()){
        SubMatrix<T> m(QN(index_j),col_range,row_range);
        m = T(0);

        const SubMatrix<T> *_sm1, *_sm2, *_sm3, *_sm4;
        if(ops[0]) _sm1 = ops[0]->block(ss[bindex[0]].qn());
        if(ops[1]) _sm2 = ops[1]->block(ss[bindex[1]].qn());
        if(ops[2]) _sm3 = ops[2]->block(ss[bindex[2]].qn());
        if(ops[3]) _sm4 = ops[3]->block(ss[bindex[3]].qn());

        int sign = 1;
        for(int iop = nops-1; iop >=0; iop--){
          if(ops[iop]->fermion()){
            int saux = 0;
            for(int ib = bindex[iop]-1; ib >= 1; ib--) saux += ss[ib].qn().n();
            sign *= SGN(saux);
          }
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2) && _sm1 && _sm2){
          for(int i3 = 0; i3 < ss[3].dim(); i3++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i2 = 0; i2 < ss[2].dim(); i2++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++){
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,j2,i3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2);
                    }
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3) && _sm1 && _sm2){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i3 = 0; i3 < ss[3].dim(); i3++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,i2,j3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i3,j3);
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK4) && _sm1 && _sm2){
          for(int i2 = 0; i2 < ss[2].dim(); i2++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i1 = 0; i1 < ss[1].dim(); i1++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(j1,i2,i3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3) && _sm1 && _sm2){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i4 = 0; i4 < ss[4].dim(); i4++)
              for(int i2 = 0; i2 < ss[2].dim(); i2++)
                for(int i3 = 0; i3 < ss[3].dim(); i3++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,j2,j3,i4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i3,j3);
        }
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK4) && _sm1 && _sm2){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i3 = 0; i3 < ss[3].dim(); i3++)
              for(int i2 = 0; i2 < ss[2].dim(); i2++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,j2,i3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j3 = 0; j3 < res_slice.size3(); j3++)
                    for(int j4 = 0; j4 < res_slice.size4(); j4++)
                      m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                        res_slice.index(i1,i2,j3,j4)-res_slice.start1()) = 
                         T(sign)*coef*_sm1->operator()(i3,j3)*_sm2->operator()(i4,j4);
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3 && _sm4){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j3 = 0; j3 < res_slice.size3(); j3++)
                        for(int j4 = 0; j4 < res_slice.size4(); j4++){
                          m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                            res_slice.index(j1,j2,j3,j4)-res_slice.start1()) = 
                             t.coef()*T(sign)*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i3,j3)*_sm4->operator()(i4,j4);
                        }
        }
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK3) && _sm1 && _sm2 && _sm3){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j3 = 0; j3 < res_slice.size3(); j3++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,j2,j3,i4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i3,j3);
        } 
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK2|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j2 = 0; j2 < res_slice.size2(); j2++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,j2,i3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i2,j2)*_sm3->operator()(i4,j4);
        } 
        if(auxt.mask == (MASK_BLOCK1|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j1 = 0; j1 < res_slice.size1(); j1++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(j1,i2,j3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i1,j1)*_sm2->operator()(i3,j3)*_sm3->operator()(i4,j4);
        } 
        if(auxt.mask == (MASK_BLOCK2|MASK_BLOCK3|MASK_BLOCK4) && _sm1 && _sm2 && _sm3){
          for(int i1 = 0; i1 < ss[1].dim(); i1++)
            for(int i2 = 0; i2 < ss[2].dim(); i2++)
              for(int i3 = 0; i3 < ss[3].dim(); i3++)
                for(int i4 = 0; i4 < ss[4].dim(); i4++)
                  for(int j2 = 0; j2 < res_slice.size2(); j2++)
                    for(int j3 = 0; j3 < res_slice.size3(); j3++)
                      for(int j4 = 0; j4 < res_slice.size4(); j4++)
                        m(v_slice.index(i1,i2,i3,i4)-v_slice.start1(),
                          res_slice.index(i1,j2,j3,j4)-res_slice.start1()) = 
                           T(sign)*coef*_sm1->operator()(i2,j2)*_sm2->operator()(i3,j3)*_sm3->operator()(i4,j4);
        } 
        column.push_back(m);
      }

    }

    cout << "Hmatrix column Lap: " << " " << clock.LapTime() << endl;
    clock.Lap();
    column.write(outputfile);
  }

//  cout << "HMATRIX SIZE " << total_size << " " << v.size()*v.size() << " " << double(total_size)/double(v.size()*v.size()) << endl;
  cout << "Hmatrix time: " << clock.TotalTime() << endl;

  outputfile.close();
}


template<class T>
void
hmatrix_product(const System<T>& s, const VectorState<T> &v, VectorState<T> &res)
{
  const Vector<T> &aux_v = v;
  Vector<T> &aux_res = res;

// Reading file
  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", s.name());
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

#ifdef USE_SPARSE
  Vector<T> data;
  Vector<size_t> index;
  for(int i = 0; i < v.size(); i++){
    index.read(inputfile);
    data.read(inputfile);
    typename Vector<T>::iterator data_iter;
    Vector<size_t>::iterator index_iter;
    data_iter = data.begin();
    index_iter = index.begin();
    data_iter++;
    index_iter++;
    aux_res(i) += data(0)*aux_v(i);
//cout << "HOLA " << i << " " << index(0) << endl;
//cout << "HOLA ---------------\n";
    for(; index_iter != index.end(); index_iter++, data_iter++){
      int j = *index_iter;
      T val = *data_iter;
//cout << "HOLA " << j << endl;
      aux_res(j) += val*aux_v(i);
      aux_res(i) += std::conj(val)*aux_v(j);
    }
  }
#else
  int ntotal = 0, nztotal = 0;
  for(int index_i = 0; index_i < v.subspaces().size(); index_i++){
    BMatrix<T> column;
    column.read(inputfile);
    
    StateSpace sub_i = v.subspaces()[index_i];
    cslice_iter<T> v_slice(&aux_v,slice(sub_i.start(),sub_i.dim(),1));

    typename BMatrix<T>::iterator biter;
    for(biter = column.begin(); biter != column.end(); biter++){
      SubMatrix<T> &block = (*biter);

/*
      for(int i = 0; i < block.size1(); i++)
        for(int j = 0; j < block.size2(); j++){
          ntotal++;
          if(fabs(block(i,j)) > 1.e-10) nztotal++;
        }
*/
             
      int index_j = block.qn().n();
      StateSpace sub_j = v.subspaces()[index_j];
      slice_iter<T> res_slice(&aux_res,slice(sub_j.start(),sub_j.dim(),1));
      const T* v_array = v_slice.get_pointer(0);
      T* res_array = &(*res_slice.begin());
      matrix_vector_product('N',block.rows(),block.cols(),block.array(),v_array,res_array,T(1),T(1));
    }

    slice_iter<T> res_slice(&aux_res,slice(sub_i.start(),sub_i.dim(),1));

    for(biter = column.begin(); biter != column.end(); biter++){
      SubMatrix<T> &block = (*biter);

      int index_j = block.qn().n();
      
      if(index_i != index_j){
        StateSpace sub_j = v.subspaces()[index_j];
        v_slice = cslice_iter<T>(&aux_v,slice(sub_j.start(),sub_j.dim(),1));
        const T* v_array = v_slice.get_pointer(0);
        T* res_array = &(*res_slice.begin());
        matrix_vector_product('C',block.rows(),block.cols(),block.array(),v_array,res_array,T(1),T(1));
      }
    }
  }
//  cout << "HOLA NZTOTAL = " << ntotal << " " << nztotal << " " << double(nztotal)/double(ntotal) << endl;

#endif // USE_SPARSE

  inputfile.close();
}

template<class T>
void
hmatrix_product_new(const System<T>& s, const VectorState<T> &v, VectorState<T> &res)
{
  const Vector<T> &aux_v = v;
  Vector<T> &aux_res = res;

// Reading file
  char file[255];
  snprintf(file, 255, "hmatrix_%s.dat", s.name());
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

  int ntotal = 0, nztotal = 0;
  BMatrix<T> column;
  for(int i = 0; i < s.aux_terms.size(); i++){
    column.read(inputfile);

    typename BMatrix<T>::iterator biter;
    for(biter = column.begin(); biter != column.end(); biter++){
      SubMatrix<T> &block = (*biter);

      int index_i = block.col_range().qn().n();
      int index_j = block.row_range().qn().n();
      StateSpace sub_i = v.subspaces()[index_i];
      StateSpace sub_j = v.subspaces()[index_j];

      cslice_iter<T> v_slice(&aux_v,slice(sub_i.start(),sub_i.dim(),1));
      slice_iter<T> res_slice(&aux_res,slice(sub_j.start(),sub_j.dim(),1));
      const T* v_array = v_slice.get_pointer(0);
      T* res_array = &(*res_slice.begin());
      matrix_vector_product('N',block.rows(),block.cols(),block.array(),v_array,res_array,T(1),T(1));

      if(index_i != index_j){
        res_slice = slice_iter<T>(&aux_res,slice(sub_i.start(),sub_i.dim(),1));
        v_slice = cslice_iter<T>(&aux_v,slice(sub_j.start(),sub_j.dim(),1));
        v_array = v_slice.get_pointer(0);
        res_array = &(*res_slice.begin());
        matrix_vector_product('C',block.rows(),block.cols(),block.array(),v_array,res_array,T(1),T(1));
      }
    }
  }
//  cout << "HOLA NZTOTAL = " << ntotal << " " << nztotal << " " << double(nztotal)/double(ntotal) << endl;

  inputfile.close();
}


template<class T>
void
hmatrix_verif(const System<T>& s)
{
  VectorState<T> aux = s.gs;
  VectorState<T> aux2 = s.gs;
  Vector<T> &auxv = aux;
  for(int i = 0; i < aux.size(); i++){
    auxv = T(0);
    auxv(i) = T(1);
    aux2 = T(0);
    aux2 = product(s,aux);
    for(int j = 0; j < auxv.size(); j++){
      cout << "HOLA VERIF " << i << " " << j << " " << aux2(j) << endl;
    } 
  }

}

//} //namepsace dmtk

#endif
