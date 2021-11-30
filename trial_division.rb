def trial_division1(n) 
    n_sqrt = Math.sqrt(n).to_i

    2.upto(n_sqrt).each{ |number|
        if n % number == 0
            return false
        end
    }

    true
end

def trial_division2(n)
    is_prime(n, $primes)
end

def primes_below_bound(bound)
    primes = Array.new

    # Iterate avoiding multiples of 2 and 3
    1.upto((bound - 1) / 6).each{ |k|
        if is_prime(6*k - 1, primes)
            primes.push(6*k - 1)
        end

        if is_prime(6*k + 1, primes)
            primes.push(6*k + 1)
        end
    }

    if is_prime(bound - 1, primes)
        primes.push(bound - 1)
    end

    primes.unshift(2,3)
end

def is_prime(n, primes)
    n_sqrt = Math.sqrt(n).to_i

    # A number is prime if there isn't a prime
    # below its square root that divides it 
    primes.each{ |prime|
        if n % prime == 0
            return false
        end

        if prime >= n_sqrt
            return true
        end
    }

    true
end

N = 13881737 * 27453061

$primes = primes_below_bound(Math.sqrt(N))

s = Time.new
puts trial_division1(N)
puts Time.new - s

s = Time.new
puts trial_division2(N)
puts Time.new - s