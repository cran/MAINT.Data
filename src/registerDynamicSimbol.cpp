#include "tle.h"
#include "msnCP_dev.h"
#include "CEMGauss.h"
// #include "tmp.h"

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

static const R_CallMethodDef CallEntries[] = {
  {"Cfasttle",       (DL_FUNC) &Cfasttle,       12},
  {"Cfulltle",       (DL_FUNC) &Cfulltle,        6},
  {"msnCP_dev",      (DL_FUNC) &msnCP_dev,      14},
  {"msnCP_dev_grad", (DL_FUNC) &msnCP_dev_grad, 13},
  {"CEMGauss",       (DL_FUNC) &CEMGauss, 13},
//  {"tmp",            (DL_FUNC) &tmp, 6},
  {NULL, NULL, 0}
};

void R_init_MAINT_Data(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
