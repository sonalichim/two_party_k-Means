#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <gmp.h>
#include "getPrime.h"

//int K = 8; //sec parameter

class PaillierSystem
{
protected:
    
    long int K, primeSize;
    bool isNum;
    mpz_t p, q, N, phi, Nsqr, invPhi;
    
    void getRandomFromGroup(mpz_t r);
    void genModulus(void);
    void genInvPhi(void);
    
public:
    
    PaillierSystem();
    PaillierSystem ( long int x);
    PaillierSystem(long int x, mpz_t pk);
    PaillierSystem(long int x, mpz_t pk, mpz_t sk);///Only for code test purpose
    
    void getThreeFactorPaillier( int * newSize, mpz_t newN, PaillierSystem **pObj);
    void encrypt(mpz_t c, mpz_t m);
    void decrypt(mpz_t m, mpz_t c);
    
    void getPublicKey(mpz_t);
    void getPhi(mpz_t privKey);
    void addCiphers(mpz_t sumCiph, mpz_t c1, mpz_t c2);
    
    ~PaillierSystem()
    {
        mpz_clear(p);
        mpz_clear(q);
        mpz_clear(N);
        mpz_clear(phi);
        mpz_clear(Nsqr);
    }
    
};





