#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <omp.h>

// Result struct for binary splitting
typedef struct {
    mpz_t P;  // product part
    mpz_t Q;  // denominator
    mpz_t T;  // numerator
} bs_result_t;

// Initialize a binary-splitting result
static inline void bs_result_init(bs_result_t *r) {
    mpz_inits(r->P, r->Q, r->T, NULL);
}

// Clear a binary-splitting result
static inline void bs_result_clear(bs_result_t *r) {
    mpz_clears(r->P, r->Q, r->T, NULL);
}

// Binary splitting recursion for Chudnovsky terms
static void binary_split(bs_result_t *res, unsigned long a, unsigned long b) {
    if (b - a == 1) {
        unsigned long k = a;
        mpz_t Pk, Qk, Tk, tmp;

        mpz_inits(Pk, Qk, Tk, tmp, NULL);

        // Pk = (6k)! / ((3k)! (k!)^3)
        mpz_fac_ui(Pk, 6 * k);
        mpz_fac_ui(Qk, 3 * k);
        mpz_divexact(Pk, Pk, Qk);

        mpz_fac_ui(Qk, k);
        mpz_pow_ui(Qk, Qk, 3);
        mpz_divexact(Pk, Pk, Qk);

        // Qk = Pk
        mpz_set(Qk, Pk);

        // Tk = Pk * (545140134k + 13591409) * (-1)^k
        mpz_mul_ui(Tk, Pk, 545140134);
        mpz_mul_ui(Tk, Tk, k);
        mpz_add_ui(Tk, Tk, 13591409);
        mpz_mul(Tk, Tk, Pk);
        if (k % 2) mpz_neg(Tk, Tk);

        mpz_set(res->P, Pk);
        mpz_set(res->Q, Qk);
        mpz_set(res->T, Tk);

        mpz_clears(Pk, Qk, Tk, tmp, NULL);
        return;
    }

    unsigned long m = (a + b) / 2;

    bs_result_t left, right;
    bs_result_init(&left);
    bs_result_init(&right);

    // Parallelize the two halves
#pragma omp parallel sections
    {
#pragma omp section
        binary_split(&left, a, m);
#pragma omp section
        binary_split(&right, m, b);
    }

    // Combine results
    mpz_mul(res->P, left.P, right.P);
    mpz_mul(res->Q, left.Q, right.Q);

    // T = right.T * left.Q + left.T * right.P
    mpz_t tmp1, tmp2;
    mpz_inits(tmp1, tmp2, NULL);
    mpz_mul(tmp1, right.T, left.Q);
    mpz_mul(tmp2, left.T, right.P);
    mpz_add(res->T, tmp1, tmp2);
    mpz_clears(tmp1, tmp2, NULL);

    bs_result_clear(&left);
    bs_result_clear(&right);
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <thousands of digits>\n", argv[0]);
        return 1;
    }

    unsigned long digits = strtoul(argv[1], NULL, 10) * 1000;
    unsigned long terms = digits / 14 + 2; // ~14 digits per Chudnovsky term

    // Precision: log2(10) ≈ 3.32193 bits per decimal digit
    mpf_set_default_prec(digits * 3.322);

    bs_result_t res;
    bs_result_init(&res);

    binary_split(&res, 0, terms);

    // Compute pi = (426880 * sqrt(10005)) / res.T
    mpf_t pi, invpi, sqrtC;
    mpf_inits(pi, invpi, sqrtC, NULL);

    mpf_sqrt_ui(sqrtC, 10005);
    mpf_mul_ui(sqrtC, sqrtC, 426880);

    mpf_set_z(invpi, res.T);
    mpf_div(invpi, sqrtC, invpi);
    mpf_ui_div(pi, 1, invpi);

    gmp_printf("%.Ff\n", pi);

    bs_result_clear(&res);
    mpf_clears(pi, invpi, sqrtC, NULL);
    return 0;
}
