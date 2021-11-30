#include <fstream>
#include <cassert>

#include "auxiliary.h"
#include "quadratic_sieve.h"
#include "GaussBin.h"


LongInt N(1);

const int B = 250;
const int L = 10;

char const* PARITIES_FILE  = "parities.txt";
char const* NULLSPACE_FILE = "nullspace.txt";


// We are going to find pairs such that x² = y² mod N. Then we can 
// try to factorize N thanks to: x² - y² = (x - y)(x + y) = 0 mod N.
// To do this we start by finding numbers such that if x² = y mod N,
// then y can be decomposed with the factorbase (y is B-smooth).
pair_type quadratic_sieve(LongInt N) {

    const auto factorbase = primes_below_bound(B);
    const int factorbase_size = factorbase.size();

    std::vector<std::vector<int>> parities;
    std::vector<base_type> ys_decomps;

    auto pairs = gen_pairs(N, factorbase, parities, ys_decomps);

    std::cout << pairs.size() << " valid pairs" << std::endl;

    // The nullspace gives us a way to multiply numbers
    // such that x² = y² mod N
    for (const auto vector : nullspace(parities, pairs.size())) {
        
        // We compute the product indicated by the solution
        LongInt x = 1;

        decomp_type y_decomp(factorbase_size, 0);

        for (int i = 0; i < vector.size(); i++) {
            if (vector[i] == 1) {
                x *= pairs[i].first;

                for (int j = 0; j < factorbase_size; j++) {
                    y_decomp[j] += ys_decomps[i][j];
                }
            }
        }

        LongInt y = 1;

        for (int i = 0; i < factorbase_size; i++) {
            assert(factorbase[i] > 1);
            assert(y_decomp[i] % 2 == 0);

            LongInt number = big_pow(factorbase[i], y_decomp[i] / 2);
            y *= number;
        }
        
        // The odds of this process giving a factorization are 
        // about one half
        LongInt aux;

        if (x > y) {
            aux = x - y;
        } else {
            aux = y - x;
        }
        
        LongInt factor = gcd(aux, N);

        if (factor > 1 && factor < N) {
            return pair_type(factor, N / factor);
        }
    }

    return pair_type(1, N);
}


std::vector<pair_type> gen_pairs(LongInt N, const base_type &factorbase, 
                                 std::vector<decomp_type> &parities, 
                                 std::vector<decomp_type> &ys_decomps) {

    std::vector<pair_type> pairs;

    int factor = 1;

    while (pairs.size() < factorbase.size() + L) {

        std::vector<pair_type> new_pairs;

        // We generate numbers with a smooth square
        for (LongInt x : probably_smooth_numbers(N, factor)) {
            LongInt y = (x * x) % N;

            if (is_factorizable(y, factorbase)) {
                new_pairs.emplace_back(x, y);
            }
        }
        
        // We decompose each pair and add it if it's sufficiently different
        for (const auto pair : new_pairs) {

            const auto y_decomp = factorize(pair.second, factorbase);
            const auto parity = calc_parity(y_decomp);

            // Only keep one factorization for each exponent parity
            if (!contains(parities, parity)) {
                pairs.push_back(pair);
                ys_decomps.push_back(y_decomp);
                parities.push_back(parity);
            }
        }

        factor++;
    }

    return pairs;
}


std::vector<LongInt> probably_smooth_numbers(LongInt N, const int factor) {
    std::vector<LongInt> numbers;

    LongInt base = (N * factor).powfn(0.5);

    for (int shift = 1; shift <= factor; shift++) {
        LongInt shift_ = shift;
        numbers.push_back(base + shift_);
    }

    return numbers;
}


bool is_factorizable(LongInt number, const base_type &factorbase) {
    if (number == 0) {
        return false;
    }

    LongInt cocient = number;

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


decomp_type factorize(LongInt number, const base_type &factorbase) {
    LongInt cocient = number;

    decomp_type exponents;
    exponents.reserve(factorbase.size());

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

    std::vector<std::vector<int>> solutions(length);

    for (int i = 0; i < length; i++) {
        std::vector<int> solution(pair_count);
     
        for (int j = 0; j < pair_count; j++) {
            ifile >> solution[j];
        }

        solutions[i] = solution;
    }

    ifile.close();
    
    return solutions;
}
