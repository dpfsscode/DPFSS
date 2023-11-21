#include "dpfss.h"


dpfss::dpfss(){
    fmpz_init(p);
    fmpz_set_str(p,"57896044618658097711785492504343953926634992332820282019728792003956564820063",10);
    fmpz_init(g);
    fmpz_set_ui(g,3);
    fmpz_init(h);
    fmpz_set_ui(h,3);
    fmpz_mod_ctx_init(ctx,p);

    // Initialize the N degree-2 m-variate polynomial;
    flint_rand_t rnd;
    flint_randinit(rnd);
    fmpz_mod_mpoly_ctx_init(mpctx, m, ORD_LEX, p);
    fk= FLINT_ARRAY_ALLOC(N, fmpz_mod_mpoly_t);
    for (int i = 0; i < N; i++){
        fmpz_mod_mpoly_init(fk[i], mpctx);
        fmpz_mod_mpoly_randtest_bound(fk[i], rnd, m*m, 2, mpctx);
    }

    core_init();
    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(EP_MTYPE);
    ep_curve_get_gen(g1_gen);
    ep2_curve_get_gen(g2_gen);
    
    pp_map_oatep_k12(gT_gen, g1_gen, g2_gen);

}

dpfss::dpfss(int input_N, int input_n, int input_t, int input_m){
    N = input_N;
    n = input_n;
    t = input_t;
    m = input_m;

    fmpz_init(p);
    fmpz_set_str(p,"57896044618658097711785492504343953926634992332820282019728792003956564820063",10);
    fmpz_init(g);
    fmpz_set_ui(g,3);
    fmpz_init(h);
    fmpz_set_ui(h,3);
    fmpz_mod_ctx_init(ctx,p);
   

    // Initialize the N degree-2 m-variate polynomial;
    flint_rand_t rnd;
    flint_randinit(rnd);
    fmpz_mod_mpoly_ctx_init(mpctx, m, ORD_LEX, p);
    fk= FLINT_ARRAY_ALLOC(N, fmpz_mod_mpoly_t);
    for (int i = 0; i < N; i++){
        fmpz_mod_mpoly_init(fk[i], mpctx);
        fmpz_mod_mpoly_randtest_bound(fk[i], rnd, m*m, 3, mpctx);
    }

    core_init();
    ep_curve_init();
    ep_param_set(B12_P381);
    ep2_curve_init();
    ep2_curve_set_twist(EP_MTYPE);
    ep_curve_get_gen(g1_gen);
    ep2_curve_get_gen(g2_gen);
    
    pp_map_oatep_k12(gT_gen, g1_gen, g2_gen);

   
}

dpfss::~dpfss(){}

fmpz** dpfss::fc_commit(fmpz_mod_poly_t f){
    fmpz** result = FLINT_ARRAY_ALLOC(t+1, fmpz*);
    for (int i = 0; i <= t; i++){
        result[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(result[i]);
        fmpz_set_ui(result[i], 1);
        fmpz_t ai;
        fmpz_init(ai);
        fmpz_mod_poly_get_coeff_fmpz(ai, f, i, ctx);
        fmpz_mod_pow_fmpz(result[i], g, ai, ctx);
    }

    return result;
}

void dpfss::share(fmpz** s){
    si = FLINT_ARRAY_ALLOC(m, fmpz**);
    com = FLINT_ARRAY_ALLOC(m, fmpz**);
    for (int i = 0; i < m; i++){
        fmpz_mod_poly_t f;
        si[i] = shamir_share(t, n, s[i], ctx, f);
        com[i] = fc_commit(f);
    }
    //std::cout << "si,si[j],com[i]: " << std::endl;
   // std::cout << sizeof(s[3]) << " " ;
   // std::cout << sizeof(si[3]) << " " ;
   // std::cout << sizeof(com[3]) << " " << std::endl;
}

void dpfss::prepare(){
    u = FLINT_ARRAY_ALLOC(n, fmpz*);
    flint_rand_t state;
    flint_randinit(state);
    for (int i = 0; i < n; i++){
        u[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(u[i]);
        fmpz_randm(u[i], state, p);
    }

    fmpz*** s1 = FLINT_ARRAY_ALLOC(n, fmpz**);
    fmpz*** com1 = FLINT_ARRAY_ALLOC(n, fmpz**);
    fmpz*** s2 = FLINT_ARRAY_ALLOC(n, fmpz**);
    fmpz*** com2 = FLINT_ARRAY_ALLOC(n, fmpz**);
    for (int i = 0; i < n; i++){
        fmpz_mod_poly_t f;
        s1[i] = shamir_share(t, n, u[i], ctx, f);
        com1[i] = fc_commit(f);
        s2[i] = shamir_share(t, n, u[i], ctx, f);
        com2[i] = fc_commit(f);
    }

    // Step 3 in protocol 2
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            fmpz_t temp_j;
            fmpz_init(temp_j);
            fmpz_set_si(temp_j, j);
            int v = fc_verify(com1[i], temp_j, s1[i][j]);
            // assert(v == 1);
        }
    }

    // Step 4 in Protocol 2
    for (int i = 0; i < n; i++){
        int v = fmpz_equal(com1[i][0], com2[i][0]);
        // assert(v == 1);
        
    }

    // Step 5 in Protocol 2
        // Generate the Vandermonde Matrix of size N*n
    M = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        M[i] = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (int j = 0; j < n; j++){
            M[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(M[i][j]);
            int t = (int) pow(j+1,i);
            fmpz_set_si(M[i][j], t);
        }
    }
        // Generate the shares, two r matrix N*n
    r1 = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        r1[i] = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (int j = 0; j < n; j++){
            r1[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(r1[i][j]);
            fmpz_set_si(r1[i][j],0);
            for (int z = 0; z < n; z++){
                //M[i][z]  s1[z][j]
                fmpz_t t;
                fmpz_init(t);
                fmpz_mod_mul(t, M[i][z], s1[z][j], ctx);
                fmpz_mod_add(r1[i][j], r1[i][j], t, ctx);
            }
        }
    }
    r2 = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        r2[i] = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (int j = 0; j < n; j++){
            r2[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(r2[i][j]);
            fmpz_set_si(r2[i][j],0);
            for (int z = 0; z < n; z++){
                fmpz_t t;
                fmpz_init(t);
                fmpz_mod_mul(t, M[i][z], s2[z][j], ctx);
                fmpz_mod_add(r2[i][j], r2[i][j], t, ctx);
            }
        }
    }
        // compute com_r1 and com_r2, which are two N*(t+1) array.
    com_r1 = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        com_r1[i] = FLINT_ARRAY_ALLOC(t+1, fmpz*);
        for (int j = 0; j < t+1; j++){
            com_r1[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(com1[i][j]);
            fmpz_set_ui(com1[i][j], 1);
            for (int z = 0; z < n; z++){
                fmpz_t t;
                fmpz_init(t);
                // com1[z][j]  M[i][z]
                fmpz_mod_pow_fmpz(t, com1[z][j], M[i][z], ctx);
                fmpz_mod_mul(com1[i][j], com1[i][j], t, ctx);
            }
        }
    }
    com_r2 = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        com_r2[i] = FLINT_ARRAY_ALLOC(t+1, fmpz*);
        for (int j = 0; j < t+1; j++){
            com_r2[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(com2[i][j]);
            fmpz_set_ui(com2[i][j], 1);
            for (int z = 0; z < n; z++){
                fmpz_t t;
                fmpz_init(t);
                // com1[z][j]  M[i][z]
                fmpz_mod_pow_fmpz(t, com2[z][j], M[i][z], ctx);
                fmpz_mod_mul(com2[i][j], com2[i][j], t, ctx);
            }
        }
    }

    //std::cout << "rk[i],comrk[i]: " << std::endl;
    //std::cout << sizeof(r2[3]) << " " ;
    //std::cout << sizeof(com_r2) << " " << std::endl;
}




int dpfss::fc_verify(fmpz** com, fmpz_t i, fmpz_t y){
    fmpz_t w;
    fmpz_init(w);
    fmpz_set_ui(w,1);

    for (int z = 0; z <= t; z++){
        fmpz_t exp;
        fmpz_init(exp);
        fmpz_pow_ui(exp, i, z);
        fmpz_mod_pow_fmpz(exp, com[z], exp, ctx);
        fmpz_mod_mul(w, w, exp, ctx);
    }
    fmpz_t gy;
    fmpz_init(gy);
    fmpz_mod_pow_fmpz(gy, g, y, ctx);
    return fmpz_equal(gy, w);
}


void dpfss::test_fc(){
    fmpz_t m;
    fmpz_init(m);
    fmpz_set_ui(m, 2333);
    fmpz_mod_poly_t f;
    fmpz** shares = shamir_share(t, n, m, ctx, f);
    fmpz** commitment = fc_commit(f);

    fmpz_t y, i;
    fmpz_init(y);
    fmpz_init(i);
    fmpz_set_ui(i,3);
    fmpz_mod_poly_evaluate_fmpz(y, f, i, ctx);
    int result = fc_verify(commitment, i, y);
}


// Input: si, com, r1, r2, com_r1, com_r2, paring parameter
// Output: Y_i com_yi
void dpfss::refresh(){
    // 1. Calculate y_{k,j} = r_{k,j}+F_k(S_j). y: N*n
    y = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        y[i] = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (int j = 0; j < n; j++){
            y[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(y[i][j]);
        }
    }
    for (int i = 0; i < N; i++){                // for k \in [N]
        for (int j = 0; j < n; j++){            // for j \in [n] -> P_j
            fmpz_t result;
            fmpz_init(result);
            fmpz** Sj = FLINT_ARRAY_ALLOC(m, fmpz*);
            for (int k = 0; k < m; k++){        
                Sj[k] = FLINT_ARRAY_ALLOC(1, fmpz);
                fmpz_init_set(Sj[k], si[k][j]);
            }
            fmpz_mod_mpoly_evaluate_all_fmpz(result, fk[i], Sj, mpctx);
            fmpz_mod_add_fmpz(result, result, r1[i][j], ctx);
            fmpz_set(y[i][j], result);
        }
    }


    // 2. Check the equation correctness of y_{k,j}
    for (int k = 0; k < N; k++){
        for (int j = 0; j < n; j++){
            fp12_t lhs, rhs;
            fp12_new(lhs);
            fp12_new(rhs);
            bn_t e1;
            fmpz2bn(y[k][j], e1);
            fp12_exp(lhs, gT_gen, e1);

            // the first term in rhs
            fp12_t temp;
            fp12_new(temp);
            ep_t g1;
            ep_new(g1);
            fmpz2bn(r1[k][j], e1);
            ep_mul_gen(g1, e1);
            core_init();
            ep_curve_init();
            ep_param_set(B12_P381);
            ep2_curve_init();
            ep2_curve_set_twist(EP_MTYPE);
            pp_map_oatep_k12(temp, g1, g2_gen);
            fp12_mul_basic(rhs, rhs, temp);

            // the second term in rhs
            // compute a_0^k firstly
            fmpz_t coef;
            fmpz_init(coef);
            ulong* exp = new ulong[m];
            for(int t = 0; t < m; t++){
                exp[t] = 0;
            }
            fmpz_mod_mpoly_get_coeff_fmpz_ui(coef, fk[k], exp, mpctx);
            fmpz2bn(coef, e1);
            fp12_exp(temp, gT_gen, e1);
            fp12_mul_basic(rhs, rhs, temp);


            // the third term in the right side
            // which is a product of m terms
            for (int l1 = 0; l1 < m; l1++){
                fmpz2bn(si[l1][k], e1);
                ep_mul_gen(g1, e1);
                pp_map_oatep_k12(temp, g1, g2_gen);
                exp[l1] = 1;
                fmpz_mod_mpoly_get_coeff_fmpz_ui(coef, fk[k], exp, mpctx);
                exp[l1] = 0;
                fmpz2bn(coef, e1);
                fp12_exp(temp, temp, e1);
                fp12_mul_basic(rhs, rhs, temp);
            }


            // the fourth term in the right
            // the product of m^2 terms
            fmpz_t exponent;
            fmpz_init(exponent);
            for (int l1 = 0; l1 < m; l1++){
                exp[l1] ++;
                for (int l2 = 0; l2 < m; l2++){
                    fmpz_mod_mul_fmpz(exponent, si[l1][j], si[l2][j], ctx);
                    exp[l2] ++;
                    fmpz_mod_mpoly_get_coeff_fmpz_ui(coef, fk[k], exp, mpctx);
                    fmpz_mod_mul_fmpz(exponent, exponent, coef, ctx);
                    fmpz2bn(exponent, e1);
                    fp12_exp(temp, gT_gen, e1);
                    fp12_mul_basic(rhs, rhs, temp);
                }
            }

            // Finally verify the equation. 
            fp12_cmp(lhs, rhs);
        }
    }
    fmpz** pts = FLINT_ARRAY_ALLOC(2*t+1, fmpz*);
    for (int i = 0; i < 2*t+1; i++){
        pts[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(pts[i]);
        fmpz_set_ui(pts[i], i+1);
    }
    fmpz** yk = FLINT_ARRAY_ALLOC(N, fmpz*);
    for (int i = 0; i < N; i++){
        yk[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(yk[i]);
        shamir_recon(yk[i], t, y[i], pts, ctx);
    }

    fmpz*** com_yk= FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        com_yk[i] = FLINT_ARRAY_ALLOC(t+1, fmpz*);
        for (int j = 0; j < t+1; j++){
            com_yk[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init_set_ui(com_yk[i][j], 1);
        }
        fmpz_mod_pow_fmpz(com_yk[i][0], g, yk[i], ctx);
    }

    // Step 3 in refresh
    y_tilde = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        y_tilde[i] = FLINT_ARRAY_ALLOC(n, fmpz*);
        for (int j = 0; j < n; j++){
            y_tilde[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(y_tilde[i][j]);
            fmpz_mod_sub_fmpz(y_tilde[i][j], yk[i], r2[i][j], ctx);
        }
    }

    com_yt = FLINT_ARRAY_ALLOC(N, fmpz**);
    for (int i = 0; i < N; i++){
        com_yt[i] = FLINT_ARRAY_ALLOC(t+1, fmpz*);
        for (int j = 0; j < t+1; j++){
            com_yt[i][j] = FLINT_ARRAY_ALLOC(1, fmpz);
            fmpz_init(com_yt[i][j]);
            fmpz_mod_divides(com_yt[i][j], com_yk[i][j], com_r2[i][j], ctx);
        }
    }
}

void dpfss::reconstruct(){
    // Step 1 in Protocol 6
    for ( int i = 0; i < N; i++){
        for (int j = 0; j < n; j++){
            fmpz_t temp_j;
            fmpz_init(temp_j);
            fmpz_set_si(temp_j, j+1);

            //verification for ~Pi in ~C
            int v1 = fc_verify(com_yt[i], temp_j, y_tilde[i][j]);
             //assert(v1 == 1);

        }
    }


    // interpolating point
    fmpz** point = FLINT_ARRAY_ALLOC(n, fmpz*);
    for (int i = 0; i < n; i++){
        point[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(point[i]);
        fmpz_set_ui(point[i], i+1);
    }

    // Step 2 in Protocol 6
    recon_s = FLINT_ARRAY_ALLOC(N, fmpz*);
    for (int i = 0; i < N; i++){
        recon_s[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(recon_s[i]);
        shamir_recon(recon_s[i], t, y_tilde[i], point, ctx);
    }
}
