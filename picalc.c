#include <gmp.h>
#include <iostream>

void compute_pi(unsigned long digits) {
    mpf_set_default_prec(digits * 3.32193); // log2(10) ≈ 3.32193
    mpf_t pi;
    mpf_init(pi);

    mpf_t sum, term;
    mpf_init(sum);
    mpf_set_ui(sum, 0);
    mpf_init(term);

    mpz_t kfact, sixk, threek, k3fact;
    mpz_init(kfact);
    mpz_init(sixk);
    mpz_init(threek);
    mpz_init(k3fact);

    const unsigned long ITER = digits / 14 + 10; // roughly enough iterations

    for (unsigned long k = 0; k < ITER; ++k) {
        mpz_fac_ui(kfact, k);                   // k!
        mpz_ui_pow_ui(sixk, 6*k, 3);           // (6k)^3
        mpz_mul_ui(sixk, sixk, 13591409 + 545140134*k); // multiply numerator
        mpz_fac_ui(k3fact, 3*k);                // (3k)!
        
        mpf_set_z(term, sixk);
        mpf_div_z(term, term, k3fact);          // divide by (3k)!
        mpf_add(sum, sum, term);
    }

    mpf_mul_ui(pi, sum, 12);
    mpf_div_ui(pi, pi, 640320*640320*640320);

    gmp_printf("Pi ≈ %.Ff\n", pi);

    mpf_clear(pi);
    mpf_clear(sum);
    mpf_clear(term);
    mpz_clear(kfact);
    mpz_clear(sixk);
    mpz_clear(threek);
    mpz_clear(k3fact);
}

int main() {
    unsigned long digits;
    std::cout << "Enter number of digits: ";
    std::cin >> digits;

    compute_pi(digits);
    return 0;
}
