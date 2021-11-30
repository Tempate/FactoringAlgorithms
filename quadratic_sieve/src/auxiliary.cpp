#include "auxiliary.h"


std::vector<int> primes_below_bound(const int bound) {
    std::vector<int> numbers;
    numbers.reserve(bound - 2);

    for (int i = 2; i < bound; i++) {
        numbers.push_back(i);
    }

    std::vector<int> primes;

    while (numbers.size() > 0) {
        const int new_prime = numbers[0];

        // Remove all the multiples of the new prime
        for (int i = numbers.size() - 1; i >= 1; i--) {
            if (numbers[i] % new_prime == 0) {
                numbers.erase(numbers.begin() + i);
            }
        }

        // Remove and save the new prime 
        numbers.erase(numbers.begin());
        primes.push_back(new_prime);
    }

    return primes;
}


LongInt gcd(LongInt m, LongInt n) {
    LongInt p = m;
    LongInt q = n;
    
    while (p != 0 && q != 0) {
        LongInt aux = p;

        p = q;
        q = aux % q;
    }

    if (p == 0) {
        return q;
    }
    
    if (q == 0) {
        return p;
    }
}


LongInt big_pow(const int base, const int power) {
    LongInt result(1);

    for (int i = 0; i < power; i++) {
        result *= base;
    }

    return result;
}


std::vector<int> calc_parity(const std::vector<int> &vector) {
    std::vector<int> parity;
    parity.reserve(vector.size());
        
    for (const int value : vector) {
        parity.push_back(value % 2);
    }

    return parity;
}


bool contains(const std::vector<std::vector<int>> &vectors, const std::vector<int> v) {
    for (auto vector : vectors) {
        int count = 0;

        for (int i = 0; i < v.size(); i++) {
            count += v[i] == vector[i];
        }

        if (count == v.size()) {
            return true;
        }
    }

    return false;
}