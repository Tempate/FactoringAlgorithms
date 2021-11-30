#ifndef AUXILIARY_H_
#define AUXILIARY_H_

#include <iostream>
#include <vector>

#include "LongInt.h"


std::vector<int> primes_below_bound(const int bound);
LongInt gcd(LongInt m, LongInt n);

LongInt big_pow(const int base, const int power);

std::vector<int> calc_parity(const std::vector<int> &vector);

bool contains(const std::vector<std::vector<int>> &vectors, const std::vector<int> v);


#endif // #ifndef AUXILIARY_H_