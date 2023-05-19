#include <gmp.h>
#include <iostream>

/*
 *
 * Requires gmp, link with gmp using -lgmp -lgmpxx and -lstdc++
 *
 * */

void euclid(mpz_t, mpz_t, mpz_t);
void euclidean(mpz_t, mpz_t, mpz_t);
void binary(mpz_t, mpz_t, mpz_t);
void lehmer(mpz_t, mpz_t, mpz_t);

int main() {
    mpz_t integer;
    mpz_init(integer);

    mpz_t left, right, result, result_euclidean, result_binary, result_lehmer;
    mpz_inits(left, right, result, result_euclidean, result_binary, result_lehmer, nullptr);

    mpz_set_ui(left, 14);
    mpz_set_ui(right, 56);

    euclid(result, left, right);
    euclidean(result_euclidean, left, right);
    binary(result_binary, left, right);
    lehmer(result_lehmer, left, right);

    
    std::cout << "gcd(14, 56) = " << result << std::endl;
    std::cout << "gcd(14, 56) = " << result_euclidean << std::endl;
    std::cout << "gcd(14, 56) = " << result_binary << std::endl;
    std::cout << "gcd(14, 56) = " << result_lehmer << std::endl;
    std::cout << "gcd(" << left << "," << right << ") = " << result_lehmer << std::endl;
}

void lehmer(mpz_t result, mpz_t left, mpz_t right) {
    mpz_t smaller, larger;
    mpz_t smaller_prime, larger_prime;
    mpz_inits(smaller, larger, smaller_prime, larger_prime, nullptr);

    /* under- and over-estimates of each number*/
    mp_limb_t smaller_u, smaller_o, larger_u, larger_o;
    /* Corresponding quotients and remainders */
    mp_limb_t q_u, r_u, q_o, r_o;

    mpz_set(smaller, left);
    mpz_set(larger, right);


    if (mpz_cmp(smaller, larger) > 0) {
        mpz_swap(smaller, larger);
    }
    /* Begin as with the Euclidean algorithm */
    while (mpz_cmp_ui(smaller, 0)) {
        if (mpz_cmp(smaller, larger) > 0) {
            mpz_swap(smaller, larger);
        }
        if (!mpz_cmp_ui(smaller, 0)) {
            break;
        }
        /* Get the most significant limb from each number */
        smaller_u = mpz_getlimbn(smaller, mpz_size(larger) - 1);
        larger_u = mpz_getlimbn(larger, mpz_size(larger) - 1);

        smaller_o = smaller_u + 1;
        larger_o = larger_u + 1;
        
        /*
         * The larger number is much larger than the smaller number
         *  We are forced to compute a multi-precision division
         *  This happens extremely rarely
         * */
        if (!smaller_u) {
            mpz_mod(larger, larger, smaller);
            continue;
        }

        // After some steps of the euclidean algorithm
        //  the two numbers become some linear combination
        //  of the original values
        // Use a, b to denote the coefficients
        //  for the larger and smaller number
        mp_limb_t
          a_l = 0,
          b_l = 1,
          a_s = 1,
          b_s = 0;
        while (1) {
            if (smaller_o && smaller_u) {
                q_u = larger_u / smaller_o;
                q_o = larger_o / smaller_u;
            }

            if (larger_u * smaller_u != larger_o * smaller_o || q_u != q_o) {
                if(a_l == 0 && b_l == 1 && a_s == 1 && b_s == 0) {
                    mpz_mod(larger, larger, smaller);
                    break;
                }
                mpz_mul_ui(smaller_prime, smaller, a_s);
                mpz_addmul_ui(smaller_prime, larger, b_s);

                mpz_mul_ui(larger_prime, smaller, a_l);
                mpz_addmul_ui(larger_prime, larger, b_l);
                break;
            } 
            // Perform an iteration of Lehmer's algorithm
            r_u = larger_u - q_u * smaller_o;
            r_o = larger_o - q_o * smaller_u;

            // This remainder is less than both of these numbers
            // The new largest number is the original smallest number
            // The new smaller number is the remainder
            larger_u = smaller_o;
            larger_o = smaller_u;
            
            smaller_u = r_u;
            smaller_o = r_o;

            // We have 
            // remainder = larger - q * smaller
            // larger' = smaller
            // smaller' = remainder
            // if larger = a_l * smaller_0 + b_l * larger_0
            //    smaller = a_s * smaller_0 + b_s * larger_0
            // then
            //    smaller' = (a_l - q * a_s) * smaller_0 + (b_l - q * b_s) * larger_0
            //    larger' = a_s * smaller_0 + b_s * larger_0
            // thus
            //
            // a_s' = a_l - q * a_s
            // b_s' = b_l - q * b_s
            //
            // a_l' = a_s
            // b_l' = b_s

            mp_limb_t
              _a_l = a_s,
              _b_l = b_s;
            a_s = a_l - q_u * a_s;
            b_s = b_l - q_o * b_s;
            a_l = _a_l;
            b_l = _b_l; 
        }
    }
    mpz_set(result, larger);
}

void binary(mpz_t result, mpz_t left, mpz_t right) {
    mpz_t smaller, larger;
    int k = 0;
    mpz_inits(smaller, larger, nullptr);
    mpz_set(smaller, left);
    mpz_set(larger, right);
    if (mpz_cmp(smaller, larger) > 0) {
        mpz_swap(smaller, larger);
    }
    while (mpz_divisible_ui_p(smaller, 2) && mpz_divisible_ui_p(larger, 2)) {
        ++k;
        mpz_tdiv_q_ui(smaller, smaller, 2);
        mpz_tdiv_q_ui(larger, larger, 2);
    }
    while (mpz_cmp_ui(smaller, 0)) {
        while (mpz_divisible_ui_p(larger, 2)) {
            mpz_tdiv_q_ui(larger, larger, 2);
        }
        while (mpz_divisible_ui_p(smaller, 2)) {
            mpz_tdiv_q_ui(smaller, smaller, 2);
        }
        mpz_sub(larger, larger, smaller);
        if (mpz_cmp(smaller, larger) > 0) {
            mpz_swap(smaller, larger);
        }
    }
    mpz_mul_2exp(larger, larger, k);
    mpz_set(result, larger);
}

void euclidean(mpz_t result, mpz_t left, mpz_t right) {
    mpz_t smaller, larger;
    if (mpz_cmp(left, right) < 0) {
        mpz_set(smaller, left);
        mpz_set(larger, right);
    }
    else {
        mpz_set(smaller, right);
        mpz_set(larger, right);
    }

    while (mpz_cmp_ui(smaller, 0)) {
        mpz_mod(larger, larger, smaller);
        if (mpz_cmp(smaller, larger) > 0) {
            mpz_swap(smaller, larger);
        }
    }
    mpz_set(result, larger);
}

void euclid(mpz_t result, mpz_t left, mpz_t right) {
    mpz_t smaller, larger;

    mpz_inits(smaller, larger, nullptr);
    if (mpz_cmp(left, right) < 0) {
        mpz_set(smaller, left);
        mpz_set(larger, right);
    }
    else {
        mpz_set(smaller, right);
        mpz_set(larger, left);
    }

    while (mpz_cmp_ui(smaller, 0)) {
        mpz_sub(larger, larger, smaller);
        if(mpz_cmp(smaller, larger) > 0) {
            mpz_swap(smaller, larger);
        }
    }

    mpz_set(result, larger);
}
