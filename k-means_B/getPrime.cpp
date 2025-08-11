#include "getPrime.h"

using namespace std;

void getPrime(mpz_t p, int bytes)
{
    int op;
    mpz_t min;
    mpz_init_set_ui(min,2);
    mpz_pow_ui(min, min, bytes-1);
    do
    {
        getRandomB(p, bytes);
        op = mpz_probab_prime_p(p, 20);
        
    }while( op <1/*!= 2*/ || !(mpz_cmp(p, min)>0)); //mpz_cmp_ui(p,2) == 0 || ); //dont want 2 as prime and op=2 means definitely prime
    
    /*Determine whether n is prime. Return 2 if n is definitely prime, return 1 if n is probably prime (without being certain), or return 0 if n is definitely non-prime.
     This function performs some trial divisions, then reps Miller-Rabin probabilistic primality tests*/
    
    mpz_clear(min);
}


void getRandomB(mpz_t rand, int size)//get random no of "bytes" bytes
{
    char* buf;
    mpz_t seed;
    gmp_randstate_t state;
    
    buf = new char [size];
    getSeed(buf, size);
    
    gmp_randinit_default(state);
    mpz_init(seed);
    mpz_import(seed, size, 1, 1, 0, 0, buf);
    gmp_randseed(state, seed);
    
    mpz_clear(seed);
    delete[] (buf);
    
    mpz_urandomb(rand, state, size);
    
    gmp_randclear(state);
}
void getRandomN(mpz_t rand, int size, mpz_t range)//get random no between 0 to range-1
{
    //  long R;
    char* buf;
    mpz_t seed;
    gmp_randstate_t state;
    
    //R = mpz_get_ui(range);
    
    buf = new char [size];
    getSeed(buf, size);
    
    gmp_randinit_default(state);
    mpz_init(seed);
    mpz_import(seed, size, 1, 1, 0, 0, buf);
    gmp_randseed(state, seed);
    
    mpz_clear(seed);
    delete [](buf);
    
    mpz_urandomm(rand, state, range);
    
    gmp_randclear(state);
}

void getRandomR(mpz_t rand, mpz_t range)//get random no between 0 to range-1
{
    gmp_randstate_t state;
    gmp_randinit_default(state);
    mpz_urandomm(rand, state, range);
    gmp_randclear(state);
}

void getSeed(char* buffer, int len)
{
    std::ifstream seed ("/dev/urandom", std::ifstream::binary);//("/dev/random", std::ifstream::binary);//("/dev/urandom", std::ifstream::binary);
    seed.read (buffer,len);
    seed.close();
}

/*
 using namespace std;
 
 int main()
 {
 mpz_t prime;// = new mpz_t();
 mpz_init(prime);
 int K=8;
 
 getPrime(prime, K);
 
 gmp_printf ("Random var %Zd\n", prime);//comment out #include <stdio.h> too
 
 return 1;
 }
 
 //Command to compile: g++ --std=c++11  getPrime.cpp -lgmp -o getPrime
 
 */
