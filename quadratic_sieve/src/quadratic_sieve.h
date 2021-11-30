#ifndef MAIN_HPP_
#define MAIN_HPP_

#include <iostream>
#include <vector>

#include "LongInt.h"

using pair_type = std::pair<LongInt, LongInt>;
using decomp_type = std::vector<int>;
using base_type = std::vector<int>;


pair_type quadratic_sieve(LongInt N);

std::vector<LongInt> probably_smooth_numbers(LongInt N, const int min_amount);

std::vector<pair_type> gen_pairs(LongInt N, 
                                 const base_type &factorbase, 
                                 std::vector<decomp_type> &parities, 
                                 std::vector<decomp_type> &ys_decomps);

bool is_factorizable(LongInt number, const base_type &factorbase);
decomp_type factorize(LongInt number, const base_type &factorbase);

std::vector<std::vector<int>> nullspace(const std::vector<decomp_type> &parities, 
                                        const int pair_count);


#endif // #ifndef MAIN_HPP_