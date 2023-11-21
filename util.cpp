# include "util.h"


fmpz** shamir_share(int t, int n, fmpz_t message, fmpz_mod_ctx_t ctx, fmpz_mod_poly_t f){
    fmpz** result = FLINT_ARRAY_ALLOC(n, fmpz*);
    for (size_t i = 0; i < n; i++){
        result[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(result[i]);
    }

    // fmpz_mod_poly_t f;
    fmpz_mod_poly_init(f, ctx);
    flint_rand_t state;
    flint_randinit(state);
    fmpz_mod_poly_randtest_monic(f, state, t+1, ctx);
    fmpz_mod_poly_set_coeff_fmpz(f,0,message, ctx); 

    fmpz_t a;
    fmpz_init(a);
    for (size_t i = 0; i < n; i++){
        fmpz_set_ui(a, i+1);
        fmpz_mod_poly_evaluate_fmpz(result[i], f,  a, ctx);
    }
    return result;
}

void view_data(fmpz*** data, int size1, int size2){
    for (int i = 0; i < size1; i++){
        for (int j = 0; j < size2; j++){
            fmpz_print(data[i][j]);
            std::cout << " ";
        }
        std::cout << "\n";
    }
    return;
}

void fmpz2bn(fmpz_t in, bn_t out){
    bn_new(out);
    size_t size = fmpz_sizeinbase(in, 10);
    char* buffer = new char[size+2];
    freopen("/dev/null", "a", stdout);
    setbuf(stdout, buffer);
    fmpz_print(in);
    freopen ("/dev/tty", "a", stdout);
    bn_read_str(out, buffer, size+2, 10);
    fcloseall();
}


void bn2fmpz(fmpz_t out, bn_t in){
    int size = bn_size_str(in, 10);
    char *in_char = new char[size];
    bn_write_str(in_char, size, in, 10);
    fmpz_set_str(out, in_char, 10);
}


void shamir_recon(fmpz_t result, int t, fmpz** shares, fmpz** point , fmpz_mod_ctx_t ctx ){
    // fmpz_t result;
    fmpz_init(result);
    fmpz_set_ui(result,0);

    // P1,...,P_t+1 give t+1 shares; 
    // point= {1,...,t+1}
    for (int i = 0; i < t+1; i++ ){
        fmpz_t ell_i;   //Lagrange coefficient 
        fmpz_init(ell_i);
        fmpz_set_ui(ell_i,1);
        fmpz_t temp1,temp2,temp3;
        fmpz_init(temp1); 
        fmpz_init(temp2);
        fmpz_init(temp3);
        fmpz_set_ui(temp1,1); 
        fmpz_set_ui(temp2,1);
        fmpz_set_ui(temp3,1);
        for (int z = 0; z < t+1 ; z++ ){
            if (z != i){
                fmpz_mul(temp1, temp1, point[z]);
                fmpz_sub(temp2, point[z], point[i]);   // tepm2=point[z]-point[i]
                fmpz_mul(temp3, temp3, temp2);
            }  
        }
        fmpz_mod_inv(temp3, temp3, ctx);   // (temp3)^-1
        fmpz_mod_mul(ell_i, temp1, temp3, ctx);
        
        //Interpolate
        fmpz_mod_mul(shares[i], shares[i], ell_i, ctx);
        fmpz_mod_add(result, result, shares[i], ctx);
    }
    // fmpz_print(message);

}