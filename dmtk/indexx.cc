/* modifications by ZA 28/03/97:
 * - changed floats to doubles everywhere
 */

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

template<class T, class A>
void indexx(size_t n, const A& arr_orig, Vector<size_t>& indx_orig, bool ascending = true)
{
  size_t i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0;
  T a;
  Vector<size_t> istack(NSTACK);

  Vector<size_t> indx(n+1);
  A arr(n+1);
  for(i = 1; i <= n; i++) arr[i] = arr_orig[i-1];

  for (j=1;j<=n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        indxt=indx[j];
        a=arr[indxt];
        for (i=j-1;i>=1;i--) {
          if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
        }
        indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP(indx[l+1],indx[ir])
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
        SWAP(indx[l],indx[ir])
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
        SWAP(indx[l+1],indx[l])
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
        do i++; while (arr[indx[i]] < a);
        do j--; while (arr[indx[j]] > a);
        if (j < i) break;
        SWAP(indx[i],indx[j])
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) cerr << "NSTACK too small in indexx." << endl;
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }

  for(j = 0; j < n; j++) indx_orig[j] = indx[j+1]-1;
  if(!ascending)
    for(j = 0; j <= n/2-1; j++)
      SWAP(indx_orig[j],indx_orig[n-j-1]);
}

#undef M
#undef NSTACK
