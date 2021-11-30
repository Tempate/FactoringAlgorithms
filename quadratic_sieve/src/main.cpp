#include "quadratic_sieve.h"


int main() {
    LongInt N(1);
    N = "98183149570452781423651";

    auto factors = quadratic_sieve(N);

    factors.first.DecOutput("factor = ");
    factors.second.DecOutput("factor = ");

    return 0;
}
