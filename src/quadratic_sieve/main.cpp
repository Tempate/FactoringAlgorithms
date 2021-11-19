#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <math.h>

#include "main.h"

#include "GaussBin.h"


const int N = 32749 * 15173;
const int B = 70;

const int PROB_SMOOTH_NUMBERS = 200000;

char const* PARITIES_FILE  = "parities.txt";
char const* NULLSPACE_FILE = "nullspace.txt";


int main() {
    const auto factors = quadratic_sieve();
    std::cout << factors.first << " " << factors.second << std::endl;
}

// We are going to find pairs such that x² = y² mod N. Then we can 
// try to factorize N thanks to: x² - y² = (x - y)(x + y) = 0 mod N.
// To do this we start by finding numbers such that if x² = y mod N,
// then y can be decomposed with the factorbase (y is B-smooth).
pair_type quadratic_sieve() {

    const auto factorbase = primes_below_bound(B);
    const int factorbase_size = factorbase.size();

    auto pairs = numbers_with_a_smooth_square(factorbase);

    // We remove pairs that are too similar.
    // We compute the y factorization and its parity.
    std::vector<std::vector<int>> parities;
    const auto ys_decomps = factorize(pairs, parities, factorbase);

    std::cout << pairs.size() << " valid pairs" << std::endl;

    // The nullspace gives us a way to multiply numbers
    // such that x² = y² mod N
    for (const auto vector : nullspace(parities, pairs.size())) {

        // We compute the product indicated by the solution
        int x = 1;

        decomp_type y_decomp(factorbase_size, 0);

        for (int i = 0; i < factorbase_size; i++) {
            if (vector[i] == 1) {
                x *= pairs[i].first;

                for (int j = 0; j < factorbase_size; j++) {
                    y_decomp[j] += ys_decomps[i][j];
                }
            }
        }
        
        int y = 1;

        for (int i = 0; i < factorbase_size; i++) {
            y *= pow(factorbase[i], y_decomp[i] / 2);
        }

        // The odds of this process giving a factorization are 
        // about one half
        const int factor = std::gcd(y - x, N);
        
        if (factor > 1) {
            return pair_type(factor, N / factor);
        }
    }

    return pair_type(1, N);
}


base_type primes_below_bound(const int bound) {
    std::vector<int> numbers;

    for (int i = 2; i < bound; i++) {
        numbers.push_back(i);
    }

    std::vector<int> primes;

    while (numbers.size() > 0) {
        const int new_prime = numbers[0];

        for (int i = numbers.size() - 1; i >= 1; i--) {
            if (numbers[i] % new_prime == 0) {
                numbers.erase(numbers.begin() + i);
            }
        }

        numbers.erase(numbers.begin());

        primes.push_back(new_prime);
    }

    return primes;
}


// Find numbers such that if x² = y mod N, then y can be 
// decomposed with the factorbase
std::vector<pair_type> numbers_with_a_smooth_square(const base_type &factorbase) {
    std::vector<pair_type> pairs;

    for (const int x : probably_smooth_numbers(PROB_SMOOTH_NUMBERS)) {
        const int y = (x * x) % N;

        if (is_factorizable(y, factorbase)) {
            pairs.emplace_back(x, y);
        }
    }

    return pairs;
}


std::vector<int> probably_smooth_numbers(const int min_amount) {
    std::vector<int> numbers;

    for (int factor = 2; numbers.size() < min_amount; factor++) {
        const int base = sqrt(factor * N);

        for (int shift = 2; shift <= factor; shift++) {
            numbers.push_back(base + shift);
        }
    }

    return numbers;
}


bool is_factorizable(const int number, const base_type &factorbase) {
    if (number == 0) {
        return false;
    }

    int cocient = number;

    for (const int prime : factorbase) {
        while (cocient % prime == 0) {
            cocient /= prime;

            if (cocient == 1) {
                return true;
            }
        }
    }

    return cocient == 1;
}


std::vector<decomp_type> factorize(std::vector<pair_type> &pairs, std::vector<decomp_type> &parities, const base_type &factorbase) {
    std::vector<decomp_type> y_decomps;

    for (int i = pairs.size() - 1; i >= 0; i--) {

        const auto y_decomp = decompose(pairs[i].second, factorbase);
        const auto parity = calc_parity<decomp_type>(y_decomp);

        // Only keep one factorization for each exponent parity
        if (contains<decomp_type>(parities, parity)) {
            pairs.erase(pairs.begin() + i);
            continue;
        }

        y_decomps.push_back(y_decomp);
        parities.push_back(parity);
    }

    return y_decomps;
}


decomp_type decompose(const int number, const base_type &factorbase) {
    int cocient = number;

    decomp_type exponents;

    for (const int prime : factorbase) {
        int exponent = 0;

        while (cocient % prime == 0) {
            cocient /= prime;
            exponent++;
        }

        exponents.push_back(exponent);
    }

    return exponents;
}


std::vector<std::vector<int>> nullspace(const std::vector<decomp_type> &parities, const int pair_count) {

    const int factorbase_size = parities[0].size();

    // Write the parity matrix into a file
    std::ofstream ofile;
    ofile.open(PARITIES_FILE);

    ofile << parities.size() << " " << factorbase_size << std::endl;

    for (const auto parity : parities) {
        for (const auto value : parity) {
            ofile << value << " ";
        }

        ofile << std::endl;
    }

    ofile.close();

    // Find the xs such that xA = 0
    GaussBin_Elimination(PARITIES_FILE, NULLSPACE_FILE);

    // Read the nullspace from the file
    std::ifstream ifile; 
    ifile.open(NULLSPACE_FILE);

    int length;
    ifile >> length;

    std::vector<std::vector<int>> solutions;
    solutions.reserve(length);

    for (int i = 0; i < length; i++) {
        std::vector<int> solution;
        
        for (int j = 0; j < pair_count; j++) {
            int value;
            ifile >> value;
            solution.push_back(value);
        }

        solutions.push_back(solution);
    }

    ifile.close();
    
    return solutions;
}