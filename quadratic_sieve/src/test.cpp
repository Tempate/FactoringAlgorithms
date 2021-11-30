#include <cassert>
#include <array>
#include <cmath>

#include "auxiliary.h"
#include "quadratic_sieve.h"


void test_prime_generation();
void test_is_factorizable();
void test_factorize();
void test_gcd();


int main() {
    test_prime_generation();
    test_is_factorizable();
    test_factorize();
    test_gcd();
}


void test_prime_generation() {

    for (const auto prime : primes_below_bound(100)) {
        // Ensure that the primes generated don't have a divisor
        for (int number = 2; number < prime; number++) {
            assert(prime % number != 0);
        }
    }

    std::cout << "[+] Primes are being generated correctly" << std::endl;
}


void test_is_factorizable() {
    const auto factorbase = primes_below_bound(100);

    std::array<int, 3> factorizables = {
        27 * 12 * 33,
        121 * 125,
        23
    };

    for (const auto number : factorizables) {
        assert(is_factorizable(number, factorbase) == true);
    }

    std::array<int, 3> not_factorizables = {
        101, 127, 137
    };

    for (const auto number : not_factorizables) {
        assert(is_factorizable(number, factorbase) == false);
    }

    std::cout << "[+] Factorizable numbers are being identified correctly" << std::endl;
}


void test_factorize() {
    const auto factorbase = primes_below_bound(100);
    
    std::array<int, 3> numbers = {
        8 * 19 * 29 * 29 * 121,
        2 * 27 * 121 * 125, 
        23
    };

    for (const auto number : numbers) {
        const auto decomp = factorize(number, factorbase);

        int recomp = 1;

        // Ensure numbers can be recomposed
        for (int i = 0; i < factorbase.size(); i++) {
            recomp *= (int) pow(factorbase[i], decomp[i]);
        }

        assert(number == recomp);
    }

    std::cout << "[+] Numbers are being factorized correctly" << std::endl;
}


void test_gcd() {
    std::array<std::array<int, 3>, 3> tests = {{
        {2 * 5, 3 * 5, 5},
        {2 * 127, 3 * 127, 127},
        {15 * 44, 14 * 44, 44}
    }};

    for (const auto test : tests) {
        assert(gcd(test[0], test[1]) == test[2]);
    }

    std::cout << "[+] The GCD is being calculated correctly" << std::endl;
}