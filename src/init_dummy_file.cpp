// Dummy file required so that useDynLib(Rceattle, .registration=TRUE) doesn't fail on empty 'src'

#if defined(__clang__)
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wunknown-warning-option"
#endif
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#if defined(__clang__)
# pragma clang diagnostic pop
#endif

void attribute_visible R_init_Rceattle(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
