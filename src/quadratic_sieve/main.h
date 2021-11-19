#ifndef MAIN_HPP_
#define MAIN_HPP_


#include <utility>
#include <vector>
#include <array>

using pair_type = std::pair<int, int>;
using decomp_type = std::vector<int>;
using base_type = std::vector<int>;

pair_type quadratic_sieve();

base_type primes_below_bound(const int bound);

std::vector<pair_type> numbers_with_a_smooth_square(const base_type &factorbase);
std::vector<int> probably_smooth_numbers(const int min_amount);
bool is_factorizable(const int number, const base_type &factorbase);

std::vector<decomp_type> factorize(std::vector<pair_type> &pairs, std::vector<decomp_type> &parities, const base_type &factorbase);

decomp_type decompose(const int number, const base_type &factorbase);
decomp_type calc_parity(const decomp_type &decomp_type);

std::vector<std::vector<int>> nullspace(const std::vector<decomp_type> &parities, const int pair_count);


template <class T>
T calc_parity(const T &vector) {
    T parity;
        
    for (const int value : vector) {
        parity.push_back(value % 2);
    }

    return parity;
}


template <class T>
bool contains(const std::vector<T> &vector, const T element) {
    return find(vector.begin(), vector.end(), element) != vector.end();
}


#endif // #ifndef MAIN_HPP_