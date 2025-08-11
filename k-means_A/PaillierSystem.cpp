#include "PaillierSystem.h"
#include <stdio.h>

using namespace std;
using namespace std::chrono;

PaillierSystem::PaillierSystem ()
{
    K =0;
    primeSize=0;
    
    mpz_init(p);
    mpz_set_ui(p,0);
    
    mpz_init(q);
    mpz_set_ui(q,0);
    
    mpz_init(N);
    mpz_set_ui(N,0);
    
    mpz_init(Nsqr);
    mpz_set_ui(Nsqr,0);
    
    mpz_init(phi);
    mpz_set_ui(phi,0);
    
}
PaillierSystem::PaillierSystem (long int x)
{
    K = 2*x;
    primeSize=x;
    
    mpz_init(p);
    mpz_set_ui(p,0);
    
    mpz_init(q);
    mpz_set_ui(q,0);
    
    mpz_init(N);
    mpz_set_ui(N,0);
    
    mpz_init(Nsqr);
    mpz_set_ui(Nsqr,0);
    
    mpz_init(phi);
    mpz_set_ui(phi,0);
    
    this->genModulus();
}

PaillierSystem::PaillierSystem(long int x, mpz_t pk)
{
    K =2*x;
    primeSize=x;
    
    mpz_init(p);
    mpz_set_ui(p,0);
    
    mpz_init(q);
    mpz_set_ui(q,0);
    
    mpz_init(N);
    mpz_set(N, pk);
    
    mpz_init(Nsqr);
    mpz_mul(Nsqr, N, N);
    
    mpz_init(phi);
    mpz_set_ui(phi,0);
    
}

PaillierSystem::PaillierSystem(long int x, mpz_t pk, mpz_t sk)
{
    K =2*x;
    primeSize=x;
    
    mpz_init(p);
    mpz_set_ui(p,0);
    
    mpz_init(q);
    mpz_set_ui(q,0);
    
    mpz_init(N);
    mpz_set(N, pk);
    
    mpz_init(Nsqr);
    mpz_mul(Nsqr, N, N);
    
    mpz_init(phi);
    mpz_set(phi, sk);
    
    genInvPhi();
}


void PaillierSystem::getRandomFromGroup(mpz_t r)
{
    mpz_t rInverse;
    mpz_init(rInverse);
    do
    {
        getRandomN(r, K, N);
    }while (!mpz_invert (rInverse, r, N));
    
    mpz_clear(rInverse);
}

void PaillierSystem::genModulus()
{
    mpz_t t;
    mpz_init(t);
    
    getPrime(p, primeSize);
    do{
        getPrime(q, primeSize);
        if(mpz_cmp(p,q)>0) //if( p>q) t= (p-1)/q;
            mpz_div(t, p, q);
        else
            mpz_div(t, q, p);
        
    }while(mpz_cmp(p,q) == 0 || mpz_cmp_ui(t,2) >= 0);
    mpz_clear(t);
    
    //gmp_printf ("Random var p %Zd\n", p);//comment out #include <stdio.h> too
    
    //gmp_printf ("Random var q %Zd\n", q);//comment out #include <stdio.h> too
    
    mpz_init(N);
    mpz_mul(N, p, q);
    gmp_printf ("\nModulus N =%Zd ", N);
    
    mpz_init(Nsqr);
    mpz_mul(Nsqr, N, N);
    //gmp_printf ("Nsqr =%Zd\n", Nsqr);
    
    mpz_t temp1, temp2;
    mpz_init(temp1);
    mpz_init(temp2);
    
    mpz_sub_ui(temp1, p, 1);
    mpz_sub_ui(temp2, q, 1);
    
    mpz_init(phi);
    mpz_mul(phi, temp1, temp2);
    //gmp_printf ("Phi =%Zd\n", phi);
    
    genInvPhi();
    
    mpz_clear(temp1);
    mpz_clear(temp2);
    
}
void PaillierSystem::genInvPhi()
{
    mpz_init(invPhi);
    if(!mpz_invert (invPhi, phi, N))
    {
        std::cout<<std::endl<<"####ERROR####"<<std::endl;
        exit(0);
    }
}


void PaillierSystem::encrypt(mpz_t c, mpz_t m)
{
    //cout<<endl<<"Inside encrypt";
    
    mpz_t c1, r;
    mpz_init(c1);
    mpz_init(r);
    
    //mpz_set_ui(m, stoi(oss.str()));
    
    mpz_add_ui(c, N, 1);
    mpz_powm(c, c, m, Nsqr);
    
    getRandomFromGroup(r); //get a uniformly random no from group Z^*_N
    
    //gmp_printf ("\nGroup Random  =%Zd\n", r);
    
    mpz_powm(c1, r, N, Nsqr);
    mpz_set_ui(r, 0);
    mpz_mul(r, c, c1);
    
    //c = c % Nsqr;
    mpz_mod(c, r, Nsqr);
    
    mpz_clear(c1);
    mpz_clear(r);
    
    
}


void PaillierSystem::decrypt(mpz_t m, mpz_t c)
{
    //cout<<endl<<"Inside Decrypt";
    
    if(K==0)
    {
        cout<<endl<<"Can not decrypt. No secret key available."<<endl;
        exit(0);
    }
    
    mpz_t c1;//, invPhi;
    
    mpz_set_ui(m, 0);
    mpz_init(c1);
    //mpz_init(invPhi);
    
    mpz_powm(c1, c, phi,  Nsqr);
    mpz_sub_ui(c1, c1, 1);
    mpz_cdiv_q (c1, c1, N);
    
    /*if(!mpz_invert (invPhi, phi, N))
     {
     std::cout<<std::endl<<"####ERROR####"<<std::endl;
     exit(0);
     }*/
    
    mpz_mul(m, c1, invPhi); //check the function arguments
    mpz_mod(m, m, N);
    
    mpz_clear(c1);
}

void PaillierSystem::getPublicKey(mpz_t pubKey)
{
    mpz_set(pubKey, N);
}

void PaillierSystem::getPhi(mpz_t privKey)
{
    mpz_set(privKey, phi);
}

void PaillierSystem::addCiphers(mpz_t sumCiph, mpz_t c1, mpz_t c2)
{
    mpz_mul(sumCiph, c1, c2);
    mpz_mod(sumCiph, sumCiph, Nsqr);
}


/*
 int main()
 {
 //Generating PK, SK
 PaillierSystem* PS = new PaillierSystem(16);
 PS->genModulus();
 
 mpz_t encMsg, decMsg, c;
 mpz_init(encMsg);
 mpz_init(decMsg);
 mpz_init(c);
 
 mpz_set_ui(encMsg, 123);
 gmp_printf ("\nMsg= %Zd", encMsg);
 
 time_point<high_resolution_clock> start = high_resolution_clock::now();
 
 PS->encrypt(encMsg, c);
 
 time_point<high_resolution_clock> stop = high_resolution_clock::now();
 auto duration = duration_cast<microseconds>(stop-start);
 
 gmp_printf ("\nCiphertext= %Zd", c);
 
 cout<<endl<<"time of encryption :"<<duration.count()<<" microseconds"<<endl;
 
 start = high_resolution_clock::now();
 
 PS->decrypt(c, decMsg);
 
 stop = high_resolution_clock::now();
 duration = duration_cast<microseconds>(stop-start);
 
 gmp_printf ("\nDecrypted= %Zd\n", decMsg);
 
 cout<<endl<<"time of decryption :"<<duration.count()<<" microseconds"<<endl;
 
 //Object only for encrypting plaintext using given pk
 mpz_t modulus, c2, encMsg2, decMsg2;
 mpz_init(modulus);
 mpz_init(c2);
 mpz_init(encMsg2);
 mpz_init(decMsg2);
 mpz_set_ui(encMsg2, 321);
 
 //mpz_set_ui(modulus, 1407915989);
 PS->getPublicKey(modulus);
 PaillierSystem *PSEnc = new PaillierSystem(modulus);
 
 gmp_printf ("\nMsg2= %Zd", encMsg2);
 
 PSEnc->encrypt(encMsg2, c2);
 
 gmp_printf ("\nCiphertext= %Zd", c2);
 
 //mpz_add(c, c, c2);
 PS->decrypt(c2, decMsg2);
 
 gmp_printf ("\nDecrypted Msg2 =%Zd\n", decMsg2);
 
 mpz_t sumMsg, sumCiph;
 mpz_init(sumMsg);
 mpz_init(sumCiph);
 
 PSEnc->addCiphers(sumCiph, c, c2);
 
 PS->decrypt(sumCiph, sumMsg);
 gmp_printf ("\nDecrypted sumMsg =%Zd\n", sumMsg);
 
 mpz_clear(c);
 mpz_clear(c2);
 
 mpz_clear(encMsg);
 mpz_clear(encMsg2);
 mpz_clear(decMsg);
 mpz_clear(decMsg2);
 mpz_clear(modulus);
 
 delete PS;
 //delete PSEnc;
 return 0;
 }*/

