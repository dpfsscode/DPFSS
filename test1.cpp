# include "dpfss.h"
# include "util.h"
# include <iostream>

void test();


int main(){
    test();
    std::cout << "test finished" << std::endl;
    return 0;
}


void test(){
    dpfss scheme;
    int m = scheme.m; 
    fmpz** s = FLINT_ARRAY_ALLOC(m, fmpz*);
    for (int i = 0; i < m; i++){
        s[i] = FLINT_ARRAY_ALLOC(1, fmpz);
        fmpz_init(s[i]);
        fmpz_set_ui(s[i], 3456);
    }
    scheme.share(s);
    scheme.prepare();
    scheme.refresh();
    scheme.reconstruct();
}