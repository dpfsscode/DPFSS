#ifndef DPFSS_H
#define DPFSS_H

#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mpoly.h>
#include <iostream>
#include <cmath>
#include <assert.h>
#include "util.h"
extern "C" {
#include <relic/relic.h>
}


class dpfss
{
public:
    fmpz_t p;
    int lambda;
    fmpz_t g;
    fmpz_t h;
    int t = 20;      // Degree of sharing polynomial
    int n = 100;      // Number of participants
    int N = 10;      // N, where k \in [N]
    int m = 100;      // Number of secrete
    fmpz_mod_ctx_t ctx;

    fmpz_mod_mpoly_t* fk;
    fmpz_mod_mpoly_ctx_t mpctx;

    ep_t g1_gen;
    ep2_t g2_gen;
    fp12_t gT_gen;

    

    dpfss(/* args */);
    dpfss(int input_N, int input_n, int input_t, int input_m);
    ~dpfss();


    fmpz*** si;    // An m*n matrix containing all shares
    fmpz*** com;   // An m*(t+1) matrix containing m commitments
    fmpz*** M;     // Vandermonde Matrix of size N*n

    fmpz*** com_r1;     // N*(t+1) Matrix returned in process prepare
    fmpz*** com_r2;     // N*(t+1) Matrix returned in process prepare
    fmpz*** r1;         // N*n Matrix output by prepare
    fmpz*** r2;         // N*n Matrix output by prepare

    fmpz** u;  // n RCS used in prepare

    fmpz*** y;          // N*n Matrix returned by refresh

    fmpz*** y_tilde;    // N*n Matrix return refresh
    fmpz*** com_yt;     // N*(t+1) Matrix returned by refresh

    fmpz** recon_s;     // N result finally


    void share(fmpz** s);
    fmpz** fc_commit(fmpz_mod_poly_t f);
    int fc_verify(fmpz** com, fmpz_t i, fmpz_t y);
    void prepare();
    void refresh();
    void reconstruct();



    // Functions for debugging
    void test_fc();
};



#endif