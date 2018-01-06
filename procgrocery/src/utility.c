#include<R.h>
#include<stdlib.h>
#define SQ(x) ((x)*(x))

void lasttruegap(const int *ivec, const int *len, int *ovec)
{
  int length = *len;
  int i, cnt=-length-1;
  for(i=0; i<length; i++, cnt++)
  {
    if(ivec[i]) cnt = 0;
    ovec[i] = cnt;
  }
}
void lagsum(const double *ivec, const int *len,
	    const int *lag, double *ovec)
{
  int length = *len;
  int nlag = *lag;
  int i, j;
  double sum = 0.0;
  for(i = 0; i < nlag; ++i)
    sum += ivec[i];
  for(j=0 ;i<length; ++i, ++j)
  {
    ovec[i] = sum;
    sum -= ivec[j];
    sum += ivec[i];
  }
}
/* lag must > 1 */
void lagvar(const double *ivec, const int *len,
	    const int *lag, double *ovec)
{
  int length = *len;
  int nlag = *lag;
  int i, j;
  double sx = 0.0, sxsq = 0.0;
  double idn = 1.0/((double) nlag);
  double idnm1 = 1.0/((double)(nlag-1));
  for(i = 0; i < nlag; ++i)
  {
    sx += ivec[i];
    sxsq += SQ(ivec[i]);
  }
  for(j = 0; i < length; ++i, ++j)
  {
    ovec[i] = (sxsq - SQ(sx)*idn)*idnm1;
    sx -= ivec[j];
    sxsq -= SQ(ivec[j]);
    sx += ivec[i];
    sxsq += SQ(ivec[i]);
  }
}
void periodsum(const double *ivec, const int *len_,
	       const int *start_, const int *end_,
	       const int *by_, int *rend_, int *slen_,
	       double * ovec)
{
  int len = *len_, start = *start_;
  int end = *end_, by = *by_;
  if(start<0 || end < start || by <0)
    error("incorrect sequence!");
  int slen = (end-start)/by+1;
  int rend = start + (slen-1)*by;
  int sp = (slen-1)*by;
  int sm = slen*by;
  int i, j, k, idx;
  double *sum;
  if(len < rend)
  {
    *rend_ = len;
    *slen_ = slen;
    return;
  }
  sum = (double*) malloc(sizeof(double) * by);
  for(i=0; i < by; ++i) sum[i] = 0.0;
  for(i=0; i<sp; ++i)
  {
    idx = i%by;
    sum[idx] += ivec[i];
  }
  for(j = rend; i < sm; ++i, ++j)
  {
    idx = i%by;
    sum[idx] += ivec[i];
    ovec[j] = sum[idx];
  }
  for(k=0;j<len;++i,++j,++k)
  {
    idx = i%by;
    sum[idx] -= ivec[k];
    sum[idx] += ivec[i];
    ovec[j] = sum[idx];
  }
  *rend_ = rend;
  *slen_ = slen;
  free(sum);
}

void periodvar(const double *ivec, const int *len_,
	       const int *start_, const int *end_,
	       const int *by_, int *rend_, int *slen_,
	       double * ovec)
{
  int len = *len_, start = *start_;
  int end = *end_, by = *by_;
  if(start<0 || end < start || by <0)
    error("incorrect sequence!");
  int slen = (end-start)/by+1;
  int rend = start + (slen-1)*by;
  int sp = (slen-1)*by;
  int sm = slen*by;
  int i, j, k, idx;
  double *sx, *sxsq;
  double idn = 1.0/((double) slen);
  double idnm1 = 1.0/((double)(slen-1));
  if(len < rend)
  {
    *rend_ = len;
    *slen_ = slen;
    return;
  }
  sx = (double*) malloc(sizeof(double) * by);
  sxsq = (double*) malloc(sizeof(double) * by);
  for(i=0; i < by; ++i)
  {
    sx[i] = 0.0;
    sxsq[i] = 0.0;
  }
  for(i=0; i<sp; ++i)
  {
    idx = i%by;
    sx[idx] += ivec[i];
    sxsq[idx] += SQ(ivec[i]);
  }
  for(j = rend; i < sm; ++i, ++j)
  {
    idx = i%by;
    sx[idx] += ivec[i];
    sxsq[idx] += SQ(ivec[i]);
    ovec[j] = (sxsq[idx]-SQ(sx[idx])*idn)*idnm1;
  }
  for(k=0;j<len;++i,++j,++k)
  {
    idx = i%by;
    sx[idx] -= ivec[k];
    sx[idx] += ivec[i];
    sxsq[idx] -= SQ(ivec[k]);
    sxsq[idx] += SQ(ivec[i]);
    ovec[j] = (sxsq[idx]-SQ(sx[idx])*idn)*idnm1;
  }
  *rend_ = rend;
  *slen_ = slen;
  free(sx);
  free(sxsq);
}
