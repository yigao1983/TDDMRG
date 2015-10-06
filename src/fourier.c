#include "dbg.h"
#include "sys.h"
#include "fftw3.h"
#include "fourier.h"

const int fourier_forward  = FFTW_FORWARD;
const int fourier_backward = FFTW_BACKWARD;

void fourier_dcmplx_1d(int nfft, dcmplx *fr, dcmplx *ff, int isign)
{
  fftw_plan plan;
  
  check(isign == fourier_forward || isign == fourier_backward, "Invalid isign");
  
  plan = fftw_plan_dft_1d(nfft, fr, ff, isign, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  return;
  
 error:
  abort();
}

void fourier_dcmplx_3d(int nfftx, int nffty, int nfftz, dcmplx *fr, dcmplx *ff, int isign)
{
  fftw_plan plan;
  
  check(isign == fourier_forward || isign == fourier_backward, "Invalid isign");
  
  plan = fftw_plan_dft_3d(nfftx, nffty, nfftz, fr, ff, isign, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  return;
  
 error:
  abort();
}
