#ifndef UTIL_H
#define UTIL_H

#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mod_poly.h>
#include <iostream>
extern "C" {
#include <relic/relic.h>
}


fmpz** shamir_share(int t, int n, fmpz_t message, fmpz_mod_ctx_t ctx, fmpz_mod_poly_t f);


void shamir_recon(fmpz_t result, int t, fmpz** share, fmpz** point, fmpz_mod_ctx_t ctx);

void view_data(fmpz*** data, int size1, int size2);

void fmpz2bn(fmpz_t in, bn_t out);


void bn2fmpz(fmpz_t out, bn_t in);

#endif