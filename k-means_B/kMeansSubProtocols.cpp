//
//  subProtocols.cpp
//  k-means
//
//  Created by sonali chim on 23/08/18.
//  Copyright © 2018 personal. All rights reserved.
//

#include "kMeansProtocol.h"

using namespace std;
using namespace std::chrono;

void swap (  long   int *a,  long   int *b)
{
    long   int temp = *a;
    *a = *b;
    *b = temp;
}  
void permute_numbers (long   int arr[], long   int n)// Fisher–Yates shuffle Algorithm
{
    // Use a different seed value so that
    // we don't get same result each time
    // we run this program
    srand (time(NULL));
  
    // Start from the last element and swap
    // one by one. We don't need to run for
    // the first element that's why i > 0
    for (long   int i = n - 1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        long   int j = rand() % (i + 1);
  
        // Swap arr[i] with the element
        // at random index
        swap(&arr[i], &arr[j]);
    }
}

//Supporting -NO write
void kMeansProtocol::prod_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB)
{
    mpz_t mul;
    mpz_init(mul);
    mpz_set_ui(mul,0);

    mpz_mul(mul, XB, YB);
    othPaillierObj->encrypt(E_op, mul); //XBYB
    mpz_mul(E_op, E_op, E_XAYA);//XBYB+XAYA
    mpz_mod(E_op, E_op, Msqr);
    mpz_powm(mul, E_XA, YB, Msqr);//XAYB
    mpz_mul(E_op, E_op, mul);//XBYB+XAYA+XAYB
    mpz_mod(E_op, E_op, Msqr);
    mpz_powm(mul, E_YA, XB, Msqr);//XBYA
    mpz_mul(E_op, E_op, mul);//XBYB+XAYA+XAYB+XBYA
    mpz_mod(E_op, E_op, Msqr);

    mpz_clear(mul);
}

void kMeansProtocol::scalar_prod(mpz_t E_op, mpz_t E_X, mpz_t Y)
{
    mpz_powm(E_op, E_X, Y, Msqr);
}


void kMeansProtocol::neg_shared_num(mpz_t E_op, mpz_t E_op_A, mpz_t op_B, mpz_t E_XA, mpz_t XB)
{
    mpz_neg(op_B, XB);
    mpz_add_ui(op_B, op_B, 1); //1-XB
    mpz_mod(op_B,op_B, M);
    
    mpz_invert(E_op_A, E_XA, Msqr);//-XA
    mpz_mod(E_op_A, E_op_A, Msqr);
    
    mpz_t E_temp;
    mpz_init(E_temp);
    othPaillierObj->encrypt(E_temp, op_B);
    mpz_mul(E_op, E_temp, E_XA); //1-(XA+XB)
    mpz_mod(E_op, E_op, Msqr);
}

void kMeansProtocol::XOR_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB)
{
    //XOR=a+b-2ab
    mpz_t E_temp;
    mpz_init(E_temp);
    
    this->prod_shared_num(E_temp, E_XA, E_YA, E_XAYA, XB, YB);//XY
#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_temp);
#endif
    
    mpz_powm_ui(E_temp, E_temp, 2, Msqr);//2XY
#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_temp);
#endif
    mpz_invert(E_temp, E_temp, Msqr);//-2XY
#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_temp);
#endif
    mpz_mul(E_op, E_YA, E_XA);//XA+YA
#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_op);
#endif
    mpz_mul(E_op, E_op, E_temp);//XA+YA-2XY
#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_op);
#endif
    mpz_add(E_temp, XB, YB);//XB+YB
    othPaillierObj->encrypt(E_XAYA, E_temp);
    mpz_mul(E_op, E_op, E_XAYA);//X+Y-2XY
    mpz_mod(E_op, E_op, Msqr);

#ifdef xor_test
        othPaillierObj->decrypt(E_XAYA, E_op);
#endif
    
    mpz_clear(E_temp);

}
void kMeansProtocol::XOR(mpz_t E_op, mpz_t E_X, mpz_t Y)
{
    //XOR=a+b-2ab
    mpz_t E_temp;
    mpz_init(E_temp);

    mpz_powm(E_op, E_X,Y, Msqr);//XY
    mpz_powm_ui(E_op, E_op,2, Msqr);//2XY
    mpz_invert(E_op, E_op, Msqr);//-2XY
    mpz_mul(E_op, E_op, E_X);//X-2XY
    othPaillierObj->encrypt(E_temp, Y);
    mpz_mul(E_op, E_op, E_temp);//X+Y-2XY
    mpz_mod(E_op, E_op, Msqr);
    
    mpz_clear(E_temp);

}

void kMeansProtocol::OR_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB)
{
    //OR=a+b-ab
    mpz_t E_temp;
    mpz_init(E_temp);
    
    this->prod_shared_num(E_temp, E_XA, E_YA, E_XAYA, XB, YB);//XY
    mpz_invert(E_temp, E_temp, Msqr);//-XY
    mpz_mul(E_op, E_YA, E_XA);//XA+YA
    mpz_mul(E_op, E_op, E_temp);//XA+YA-XY
    mpz_add(E_temp, XB, YB);//XB+YB
    othPaillierObj->encrypt(E_temp, E_temp);
    mpz_mul(E_op, E_op,E_temp);//X+Y-XY
    mpz_mod(E_op, E_op, Msqr);

    mpz_clear(E_temp);
}


void kMeansProtocol::OR(mpz_t E_op, mpz_t E_X, mpz_t Y)
{
    //OR=a+b-ab
    mpz_t E_temp;
    mpz_init(E_temp);

    mpz_powm(E_op, E_X,Y, Msqr);//XY
    mpz_invert(E_op, E_op, Msqr);//-XY
    mpz_mul(E_op, E_op, E_X);//X-XY
    othPaillierObj->encrypt(E_temp, Y);
    mpz_mul(E_op, E_op, E_temp);//X+Y-XY
    mpz_mod(E_op, E_op, Msqr);
    
    mpz_clear(E_temp);

}

void kMeansProtocol::share(mpz_t E_op_A, mpz_t op_B, mpz_t E_op)
{
    mpz_t E_temp, temp;
    mpz_init(E_temp);
    mpz_init(temp);

    //getRandomN(op_B, L, M); // random binary no
    mpz_set_ui(op_B, 123);
    
    othPaillierObj->encrypt(E_temp, op_B);
    //othPaillierObj->decrypt(temp, E_temp);
    
    mpz_invert(E_temp, E_temp, Msqr);// = enc(-SB)
    //othPaillierObj->decrypt(temp, E_temp);

    mpz_mod(E_op, E_op, Msqr);
    //othPaillierObj->decrypt(temp, E_op);

    //othPaillierObj->decrypt(temp, E_op_A);
    mpz_mul(E_op_A, E_op, E_temp); //= enc( S -SB )
    mpz_mod(E_op_A, E_op_A, Msqr);
    //othPaillierObj->decrypt(temp, E_op_A);

    mpz_clear(E_temp);
    mpz_clear(temp);

}

void kMeansProtocol::write(mpz_t E_value_A)
{
    char * strValue = mpz_get_str(NULL, 10, E_value_A);
    connWithAlice->writeDataWithlength(&strValue);
    free(strValue);
}
void kMeansProtocol::share_write(mpz_t E_op_A, mpz_t op_B, mpz_t E_op)
{
    this->share(E_op_A, op_B, E_op);
    this->write(E_op_A);
}


void kMeansProtocol::read(mpz_t E_value_A)
{
    char * strValue;
    connWithAlice->readDataWithLength(&strValue);
    mpz_set_str(E_value_A, strValue, 10);
    free(strValue);
}

//sub-proto -write
void kMeansProtocol::MOTN(mpz_t E_cmp, mpz_t E_cmp_A, mpz_t cmp_B, mpz_t E_XA, mpz_t YA, mpz_t XB, mpz_t YB)
{
    mpz_t E_temp, temp, E_X, E_Y, E_P, E_Q, r_1, r_2, coin, bit;
    mpz_init(E_temp);
    mpz_init(temp);
    mpz_init(E_X);
    mpz_init(E_Y);
    mpz_init(E_P);
    mpz_init(E_Q);
    mpz_init(r_1);
    mpz_init(r_2);
    mpz_init(coin);
    mpz_init(bit);

    //XA+XB
    othPaillierObj->encrypt(E_temp, XB);
    mpz_mul(E_X,E_XA, E_temp);
    mpz_mod(E_X, E_X, Msqr);
    othPaillierObj->encrypt(E_temp, YB);
    mpz_mul(E_Y,E_Y, E_temp);
    mpz_mod(E_Y, E_Y, Msqr);

    //rX+r
    //getRandomN(r_1, K, N); // random binary no
    //getRandomN(r_2, K, N); // random binary no
    mpz_set_ui(r_1, 123);
    mpz_set_ui(r_2, 234);

    mpz_powm(E_P, E_X, r_1, Msqr);
    othPaillierObj->encrypt(E_temp, r_2);
    mpz_mul(E_P, E_P, E_temp);
    mpz_mod(E_P, E_P, Msqr);
    
    mpz_powm(E_Q, E_Y, r_1, Msqr);
    othPaillierObj->encrypt(E_temp, r_2);
    mpz_mul(E_Q, E_Q, E_temp);
    mpz_mod(E_Q, E_Q, Msqr);

    //flip a coin
    int size=1;
    mpz_set_ui(temp, 2); //range
    //getRandomN(coin, size, temp);
    mpz_set_ui(coin, 1);

    //send P and Q
    if(mpz_cmp_ui(coin, 1)==0)//if coin=1
    {
        write(E_Q);
        write(E_P);
    }
    else
    {
        write(E_P);
        write(E_Q);
    }

    
    read(bit);
    
#ifdef test
    othPaillierObj->decrypt(temp, bit);
#endif
    
    XOR(E_cmp, bit, coin);
    
    share_write(E_cmp_A, cmp_B, E_cmp);
    
#ifdef test
    othPaillierObj->decrypt(temp, E_cmp);
#endif
    
    mpz_clear(E_temp);
    mpz_clear(temp);
    mpz_clear(E_X);
    mpz_clear(E_Y);
    mpz_clear(E_P);
    mpz_clear(E_Q);
    mpz_clear(r_1);
    mpz_clear(r_2);
    mpz_clear(coin);
    mpz_clear(bit);
}

void kMeansProtocol::MOTN_value(mpz_t E_min, mpz_t E_min_A, mpz_t min_B, mpz_t E_XA, mpz_t E_YA, mpz_t XB, mpz_t YB)
{
    //cout<<"\nIn MOTN_VAL"<<endl;
    mpz_t E_cmp, E_cmp_A, cmp_B, temp, E_neg_cmp_A, neg_cmp_B, E_neg_cmp, E_mul1, E_mul2;
    mpz_init(E_cmp);
    mpz_init(E_cmp_A);
    mpz_init(cmp_B);
    mpz_init(temp);
    mpz_init(E_neg_cmp_A);
    mpz_init(neg_cmp_B);
    mpz_init(E_neg_cmp);
    mpz_init(E_mul1);
    mpz_init(E_mul2);
    

    
    MOTN(E_cmp, E_cmp_A, cmp_B, E_XA, E_YA, XB, YB);

#ifdef motn_val_test
    othPaillierObj->decrypt(temp, E_XA);
    mpz_add(temp, temp, XB);
    mpz_mod(temp, temp, M);
    gmp_printf("\nMin of :%Zd, ",temp);
    othPaillierObj->decrypt(temp, E_YA);
    mpz_add(temp, temp, YB);
    mpz_mod(temp, temp, M);
    gmp_printf("%Zd = ",temp);
    othPaillierObj->decrypt(temp, E_cmp);
    gmp_printf("%Zd",temp);
#endif
    
    neg_shared_num(E_neg_cmp, E_neg_cmp_A, neg_cmp_B, E_cmp_A, cmp_B);
    
    read(E_mul1);
    
    //compute E_min = neg_cmp X + cmp Y
    prod_shared_num(E_min, E_neg_cmp_A, E_XA, E_mul1, neg_cmp_B, XB);
    //othPaillierObj->decrypt(temp, E_min);
    
    read(E_mul2);

#ifdef motn_val_test
    othPaillierObj->decrypt(temp, E_cmp_A);
    gmp_printf("  %Zd + %Zd =",temp, cmp_B);
    mpz_add(temp, temp, cmp_B);
    mpz_mod(temp, temp, M);
    gmp_printf("%Zd ",temp);
    othPaillierObj->decrypt(temp, E_YA);
    gmp_printf("  %Zd + %Zd =",temp, YB);
    mpz_add(temp, temp, YB);
    mpz_mod(temp, temp, M);
    gmp_printf("%Zd ",temp);
    othPaillierObj->decrypt(E_mul1, E_mul2);
    gmp_printf("mul %Zd ",E_mul1);
#endif
    
    prod_shared_num(temp, E_cmp_A, E_YA, E_mul2, cmp_B, YB);
    //othPaillierObj->decrypt(E_mul1, temp);

    mpz_mul(E_min, E_min, temp);
    mpz_mod(E_min, E_min, Msqr);
    
    share_write(E_min_A, min_B, E_min);

#ifdef motn_val_test
    othPaillierObj->decrypt(temp, E_min);
    gmp_printf(" Value %Zd",temp);
#endif
    
    mpz_clear(E_cmp);
    mpz_clear(E_cmp_A);
    mpz_clear(cmp_B);
    mpz_clear(temp);
    mpz_clear(E_neg_cmp_A);
    mpz_clear(neg_cmp_B);
    mpz_clear(E_neg_cmp);
    mpz_clear(E_mul1);
    mpz_clear(E_mul2);
}

void kMeansProtocol::MOTN_termination(bool *cmp, mpz_t E_XA, mpz_t E_YA, mpz_t XB, mpz_t YB)
{
    mpz_t E_temp, temp, E_X, E_Y, E_P, E_Q, r_1, r_2, bit;
    mpz_init(E_temp);
    mpz_init(temp);
    mpz_init(E_X);
    mpz_init(E_Y);
    mpz_init(E_P);
    mpz_init(E_Q);
    mpz_init(r_1);
    mpz_init(r_2);
    mpz_init(bit);
    
    //XA+XB
    othPaillierObj->encrypt(E_temp, XB);
    mpz_mul(E_X, E_XA, E_temp);
    mpz_mod(E_X, E_X, Msqr);
    othPaillierObj->encrypt(E_temp, YB);
    mpz_mul(E_Y, E_YA, E_temp);
    mpz_mod(E_Y, E_Y, Msqr);

    
    /*othPaillierObj->decrypt(r_1, E_X);
    gmp_printf("\n\n%Zd ",r_1);
    othPaillierObj->decrypt(r_1, E_Y);
    gmp_printf("%Zd ",r_1);*/
    
    //rX+r
    //getRandomN(r_1, K, N); // random binary no
    //getRandomN(r_2, K, N); // random binary no
    
    mpz_set_ui(r_1, 1);
    mpz_set_ui(r_2, 2);
    
    mpz_powm(E_P, E_X, r_1, Msqr);
    othPaillierObj->encrypt(E_temp, r_2);
    mpz_mul(E_P, E_P, E_temp);
    mpz_mod(E_P, E_P, Msqr);
    
    mpz_powm(E_Q, E_Y, r_1, Msqr);
    othPaillierObj->encrypt(E_temp, r_2);
    mpz_mul(E_Q, E_Q, E_temp);
    mpz_mod(E_Q, E_Q, Msqr);

        write(E_P);
        write(E_Q);
 
    
    read(bit);
    *cmp= (bool)mpz_get_ui(bit);
    
    mpz_clear(E_temp);
    mpz_clear(temp);
    mpz_clear(E_X);
    mpz_clear(E_Y);
    mpz_clear(E_P);
    mpz_clear(E_Q);
    mpz_clear(r_1);
    mpz_clear(r_2);
    mpz_clear(bit);
}

void kMeansProtocol::MOMN_n(mpz_t E_cmp[NUMDATAPOINT], mpz_t E_cmp_A[NUMDATAPOINT], mpz_t cmp_B[NUMDATAPOINT], mpz_t E_X[NUMDATAPOINT], mpz_t E_XA[NUMDATAPOINT], mpz_t XB[NUMDATAPOINT])
{
    long   int permutation[NUMDATAPOINT];
    mpz_t temp, E_P, E_Q, r_1, r_2, coin, bit, E_W_A, W_B, E_neg_A, neg_B, E_W_min_A, W_min_B, E_W_min;
    mpz_init(temp);
    mpz_init(E_P);
    mpz_init(E_Q);
    mpz_init(r_1);
    mpz_init(r_2);
    mpz_init(E_W_A);
    mpz_init(W_B);
    mpz_init(E_neg_A);
    mpz_init(neg_B);
    mpz_init(E_W_min_A);
    mpz_init(W_min_B);
    mpz_init(E_W_min);
    mpz_init(coin);
    mpz_init(bit);


    //permute_numbers(permutation, NUMDATAPOINT);
    
    /*for(int i =0; i<NUMDATAPOINT; ++i)
    {
        mpz_init(W[i]);
        mpz_set(W[i], E_X[permutation[i]]);
    }*/
    
    //mpz_set(E_W_min, E_X[permutation[0]]);
    mpz_set(E_W_min, E_X[0]);

    mpz_set_ui(temp, 0);
    othPaillierObj->encrypt(E_cmp_A[0], temp);
    mpz_set_ui(cmp_B[0], 1);
    othPaillierObj->encrypt(E_cmp[0], cmp_B[0]);
    
#ifdef  momn_test
    cout<<"\nMOMN input :";
    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        othPaillierObj->decrypt(temp, E_XA[i]);
        gmp_printf("\n%Zd+", temp);
        gmp_printf("%Zd=", XB[i]);

        othPaillierObj->decrypt(temp, E_X[i]);
        gmp_printf("%Zd", temp);

    }
#endif
    
    
    for(long   int i =1; i<NUMDATAPOINT; ++i)
    {
        //rX+r
        //getRandomN(r_1, K, N); // random binary no
        //getRandomN(r_2, K, N); // random binary no
        mpz_set_ui(r_1, 1);
        mpz_set_ui(r_2, 2);
        
        mpz_powm(E_P, E_W_min, r_1, Msqr);
        othPaillierObj->encrypt(temp, r_2);
        mpz_mul(E_P, E_P, temp);
        mpz_mod(E_P, E_P, Msqr);
        
        //mpz_powm(E_Q, E_X[permutation[i]], r_1, Msqr);
        mpz_powm(E_Q, E_X[i], r_1, Msqr);

        othPaillierObj->encrypt(temp, r_2);
        mpz_mul(E_Q, E_Q, temp);
        mpz_mod(E_Q, E_Q, Msqr);

#ifdef  momn_test
        othPaillierObj->decrypt(temp, E_W_min);
        gmp_printf("\n Com: %Zd", temp);
        othPaillierObj->decrypt(temp, E_X[i]);
        gmp_printf(", %Zd", temp);
#endif
    
        
        //flip a coin
        int size=1;
        mpz_set_ui(temp, 2); //range
        //getRandomN(coin, size, temp);
        mpz_set_ui(coin, 0);
        
        //send P and Q
        if(mpz_cmp_ui(coin, 1)==0)//if coin=1
        {
            write(E_Q);
            write(E_P);
        }
        else
        {
            write(E_P);
            write(E_Q);
        }
        
        //party A does processing and sends back the result bit
        
        read(bit);
        
#ifdef  momn_test
        othPaillierObj->decrypt(temp, bit);
#endif
        //XOR(E_cmp[permutation[i]], bit, coin); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector
        XOR(E_cmp[i], bit, coin); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector

        if(i < NUMDATAPOINT-1)
        {
            share_write(E_W_min_A, W_min_B, E_W_min);
            //share_write(E_W_A, W_B, E_X[permutation[0]]);
            share_write(E_W_A, W_B, E_X[i]);

            //share_write(E_cmp_A[permutation[i]], cmp_B[permutation[i]], E_cmp[permutation[i]]);
            share_write(E_cmp_A[i], cmp_B[i], E_cmp[i]);
            
#ifdef  momn_test
            othPaillierObj->decrypt(temp, E_cmp_A[i]);
            gmp_printf("\n\n%Zd ", temp);
#endif
            //neg_shared_num(temp, E_neg_A, neg_B, E_cmp_A[permutation[i]], cmp_B[permutation[i]]);
            neg_shared_num(temp, E_neg_A, neg_B, E_cmp_A[i], cmp_B[i]);

            read(E_P); //neg_cmp_A*W_min_A
            read(E_Q);//cmp_A * W_A
            
            prod_shared_num(E_W_min, E_neg_A, E_W_min_A, E_P, neg_B, W_min_B);
#ifdef  momn_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
            
//            prod_shared_num(temp, E_cmp_A[permutation[i]], E_W_A, E_Q, cmp_B[permutation[i]], W_B);
  //          prod_shared_num(temp, E_cmp_A[permutation[i]], E_W_A, E_Q, cmp_B[permutation[i]], W_B);
            prod_shared_num(temp, E_cmp_A[i], E_W_A, E_Q, cmp_B[i], W_B);
#ifdef  momn_test
            othPaillierObj->decrypt(E_P, temp);
#endif
            mpz_mul(E_W_min, E_W_min, temp);
#ifdef  momn_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
        
            mpz_mod(E_W_min, E_W_min, Msqr);
            
#ifdef  momn_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
            
        }
        else
        {
            //share_write(E_cmp_A[permutation[i]], cmp_B[permutation[i]], E_cmp[permutation[i]]);
            share_write(E_cmp_A[i], cmp_B[i], E_cmp[i]);
            
#ifdef  momn_test
            othPaillierObj->decrypt(temp, E_cmp_A[i]);
            gmp_printf("\n\n%Zd ", temp);
#endif
        }
        
#ifdef  momn_test
        othPaillierObj->decrypt(temp, E_cmp_A[i]);
        gmp_printf("\n%d cmp: %Zd+%Zd ", i, temp, cmp_B[i]);
        othPaillierObj->decrypt(temp, E_cmp[i]);
        gmp_printf("%Zd ", temp);
        
        cout<<"\n\n :";
        for(long   int i =0; i<NUMDATAPOINT; ++i)
        {
            gmp_printf("\n%Zd+ ", cmp_B[i]);

            othPaillierObj->decrypt(temp, E_cmp_A[i]);
            gmp_printf("%Zd ", temp);
            othPaillierObj->decrypt(temp, E_cmp[i]);
            gmp_printf("%Zd ", temp);

        }
        cout<<"\n";
#endif

    }
    
#ifdef  momn_test
            cout<<"\nMOMN output1 :";
            for(long   int i =0; i<NUMDATAPOINT; ++i)
            {
                gmp_printf("\n%Zd+ ", cmp_B[i]);

                othPaillierObj->decrypt(temp, E_cmp_A[i]);
                gmp_printf("%Zd ", temp);
                othPaillierObj->decrypt(temp, E_cmp[i]);
                gmp_printf("%Zd ", temp);

            }
#endif
    
    //get only one 1 in CMp vector...using w as y variable from protocol
    mpz_set_ui(W_B, 0);
    othPaillierObj->encrypt(E_W_A, W_B);
    
    for(long int i=NUMDATAPOINT-2; i>=0; --i)
    {

        
        read(E_P);
        
//        XOR_shared_num(bit, E_W_A, E_cmp_A[permutation[i]], E_P, W_B, cmp_B[permutation[i]]);
    
        mpz_add(W_B, W_B, cmp_B[i+1]);
        mpz_mul(E_W_A, E_W_A, E_cmp_A[i+1]);
        
#ifdef  momn_test
        othPaillierObj->decrypt(temp, E_W_A);
        gmp_printf("\nW: %Zd+%Zd ", temp, W_B);

        othPaillierObj->decrypt(temp, E_cmp_A[i]);
        gmp_printf("cmp: %Zd+%Zd ", temp, cmp_B[i]);

        othPaillierObj->decrypt(temp, E_P);
        gmp_printf("cmp: %Zd ", temp, E_P);

#endif
#ifdef  momn_test
       othPaillierObj->decrypt(temp, E_cmp_A[i]);
#endif
        
        XOR_shared_num(bit, E_W_A, E_cmp_A[i], E_P, W_B, cmp_B[i]);
#ifdef  momn_test
        othPaillierObj->decrypt(temp, E_cmp_A[i]);
#endif

/*#ifdef  momn_test
        othPaillierObj->decrypt(temp, bit);
        gmp_printf("XOR=%Zd ", temp);

#endif*/
        
        share_write(E_W_min_A, W_min_B, bit);
        
        neg_shared_num(temp, E_neg_A, neg_B, E_W_A, W_B);
        
#ifdef  momn_test
        othPaillierObj->decrypt(temp, E_neg_A);
        othPaillierObj->decrypt(temp, E_cmp_A[i]);
#endif
        
        read(E_P);
        
//        prod_shared_num(E_cmp[permutation[i]], E_W_min_A, E_neg_A, E_W_min_A, W_min_B, neg_B);
        prod_shared_num(E_cmp[i], E_W_min_A, E_neg_A, E_P, W_min_B, neg_B);
        
        share_write(E_cmp_A[i], cmp_B[i], E_cmp[i]);

        
#ifdef  momn_test
        gmp_printf("\n%Zd+ ", cmp_B[i]);

        othPaillierObj->decrypt(temp, E_cmp_A[i]);
        gmp_printf("%Zd ", temp);
        
        othPaillierObj->decrypt(temp, E_cmp[i]);
        gmp_printf("=%Zd ", temp);

#endif

    }
    
    //unpermute cmp[permutation[i]] which is cmp[i]
    
    //share cmp[i]
#ifdef momn_test

    cout<<"\nMOMN output :";

    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        //share_write(E_cmp_A[i], cmp_B[i], E_cmp[i]);
        othPaillierObj->decrypt(temp, E_cmp[i]);
        gmp_printf("%Zd ", temp);

    }
#endif

    //share_write(E_cmp_A[i], cmp_B[i], E_cmp[permutation[i]]);
    
    mpz_clear(temp);
    mpz_clear(E_P);
    mpz_clear(E_Q);
    mpz_clear(r_1);
    mpz_clear(r_2);
    mpz_clear(E_W_A);
    mpz_clear(W_B);
    mpz_clear(E_neg_A);
    mpz_clear(neg_B);
    mpz_clear(E_W_min_A);
    mpz_clear(W_min_B);
    mpz_clear(E_W_min);
    mpz_clear(coin);
    mpz_clear(bit);
    
}

void kMeansProtocol::MOMN_k(mpz_t E_X[NUMCENTERS], mpz_t E_XA[NUMCENTERS], mpz_t XB[NUMCENTERS], long int dp_num)
{
    long   int permutation[NUMCENTERS];
    mpz_t temp, E_P, E_Q, r_1, r_2, coin, bit, E_W_A, W_B, E_neg_A, neg_B, E_W_min_A, W_min_B, E_W_min;
    mpz_init(temp);
    mpz_init(E_P);
    mpz_init(E_Q);
    mpz_init(r_1);
    mpz_init(r_2);
    mpz_init(E_W_A);
    mpz_init(W_B);
    mpz_init(E_neg_A);
    mpz_init(neg_B);
    mpz_init(E_W_min_A);
    mpz_init(W_min_B);
    mpz_init(E_W_min);
    mpz_init(coin);
    mpz_init(bit);


    //permute_numbers(permutation, NUMCENTERS);
    
    /*for(int i =0; i<NUMCENTERS; ++i)
    {
        mpz_init(W[i]);
        mpz_set(W[i], E_X[permutation[i]]);
    }*/
    
    //mpz_set(E_W_min, E_X[permutation[0]]);
    mpz_set(E_W_min, E_X[0]);

    mpz_set_ui(temp, 0);
    othPaillierObj->encrypt(E_closest_center_A[dp_num][0], temp);
    mpz_set_ui(closest_center_B[dp_num][0], 1);
    othPaillierObj->encrypt(E_closest_center[dp_num][0], closest_center_B[dp_num][0]);
    
#ifdef  momnk_test
    cout<<"\nMOMN input :";
    for(long   int i =0; i<NUMCENTERS; ++i)
    {
        othPaillierObj->decrypt(temp, E_XA[i]);
        gmp_printf("%Zd+", temp);
        gmp_printf("%Zd=", XB[i]);

        othPaillierObj->decrypt(temp, E_X[i]);
        gmp_printf("%Zd ", temp);

    }
#endif
    
    
    for(long   int i =1; i<NUMCENTERS; ++i)
    {
        //rX+r
        //getRandomN(r_1, K, N); // random binary no
        //getRandomN(r_2, K, N); // random binary no
        mpz_set_ui(r_1, 1);
        mpz_set_ui(r_2, 2);
        
        mpz_powm(E_P, E_W_min, r_1, Msqr);
        othPaillierObj->encrypt(temp, r_2);
        mpz_mul(E_P, E_P, temp);
        mpz_mod(E_P, E_P, Msqr);
        
        //mpz_powm(E_Q, E_X[permutation[i]], r_1, Msqr);
        mpz_powm(E_Q, E_X[i], r_1, Msqr);

        othPaillierObj->encrypt(temp, r_2);
        mpz_mul(E_Q, E_Q, temp);
        mpz_mod(E_Q, E_Q, Msqr);

#ifdef  momnk_test
        othPaillierObj->decrypt(temp, E_W_min);
        gmp_printf("\n Com: %Zd", temp);
        othPaillierObj->decrypt(temp, E_X[i]);
        gmp_printf(", %Zd", temp);
#endif
    
        
        //flip a coin
        int size=1;
        mpz_set_ui(temp, 2); //range
        //getRandomN(coin, size, temp);
        mpz_set_ui(coin, 0);
        
        //send P and Q
        if(mpz_cmp_ui(coin, 1)==0)//if coin=1
        {
            write(E_Q);
            write(E_P);
        }
        else
        {
            write(E_P);
            write(E_Q);
        }
        
        //party A does processing and sends back the result bit
        
        read(bit);
        
#ifdef  momnk_test
        othPaillierObj->decrypt(temp, bit);
#endif
        //XOR(E_cmp[permutation[i]], bit, coin); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector
        XOR(E_closest_center[dp_num][i], bit, coin); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector

        if(i < NUMCENTERS-1)
        {
            share_write(E_W_min_A, W_min_B, E_W_min);
            //share_write(E_W_A, W_B, E_X[permutation[0]]);
            share_write(E_W_A, W_B, E_X[i]);

            //share_write(E_cmp_A[permutation[i]], cmp_B[permutation[i]], E_cmp[permutation[i]]);
            share_write(E_closest_center_A[dp_num][i], closest_center_B[dp_num][i], E_closest_center[dp_num][i]);

            //neg_shared_num(temp, E_neg_A, neg_B, E_cmp_A[permutation[i]], cmp_B[permutation[i]]);
            neg_shared_num(temp, E_neg_A, neg_B, E_closest_center_A[dp_num][i], closest_center_B[dp_num][i]);

            read(E_P); //neg_cmp_A*W_min_A
            read(E_Q);//cmp_A * W_A
            
            prod_shared_num(E_W_min, E_neg_A, E_W_min_A, E_P, neg_B, W_min_B);
#ifdef  momnk_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
            
//            prod_shared_num(temp, E_cmp_A[permutation[i]], E_W_A, E_Q, cmp_B[permutation[i]], W_B);
  //          prod_shared_num(temp, E_cmp_A[permutation[i]], E_W_A, E_Q, cmp_B[permutation[i]], W_B);
            prod_shared_num(temp, E_closest_center_A[dp_num][i], E_W_A, E_Q, closest_center_B[dp_num][i], W_B);
#ifdef  momnk_test
            othPaillierObj->decrypt(E_P, temp);
#endif
            mpz_mul(E_W_min, E_W_min, temp);
#ifdef  momnk_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
        
            mpz_mod(E_W_min, E_W_min, Msqr);
            
#ifdef  momnk_test
            othPaillierObj->decrypt(temp, E_W_min);
#endif
            
        }
        else
        {
            //share_write(E_cmp_A[permutation[i]], cmp_B[permutation[i]], E_cmp[permutation[i]]);
            share_write(E_closest_center_A[dp_num][i], closest_center_B[dp_num][i], E_closest_center[dp_num][i]);
        }
        
#ifdef  momnk_test
        othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
        gmp_printf("\n%d cmp: %Zd+%Zd ", i, temp, closest_center_B[dp_num][i]);
        othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
        gmp_printf("%Zd ", temp);
        
        cout<<"\n\n :";
        for(long   int i =0; i<NUMCENTERS; ++i)
        {
            gmp_printf("\n%Zd+ ", closest_center_B[dp_num][i]);

            othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
            gmp_printf("%Zd ", temp);
            othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
            gmp_printf("%Zd ", temp);

        }
        cout<<"\n";
#endif

    }
    
#ifdef  momnk_test
            cout<<"\nMOMN output1 :";
            for(long   int i =0; i<NUMCENTERS; ++i)
            {
                gmp_printf("\n%Zd+ ", closest_center_B[dp_num][i]);

                othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
                gmp_printf("%Zd ", temp);
                othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
                gmp_printf("%Zd ", temp);

            }
#endif
    
    //get only one 1 in CMp vector...using w as y variable from protocol
    mpz_set_ui(W_B, 0);
    othPaillierObj->encrypt(E_W_A, W_B);
    
    for(long int i=NUMCENTERS-2; i>=0; --i)
    {

//        XOR_shared_num(bit, E_W_A, E_cmp_A[permutation[i]], E_P, W_B, cmp_B[permutation[i]]);
    
        mpz_add(W_B, W_B, closest_center_B[dp_num][i+1]);
        
        mpz_mul(E_W_A, E_W_A, E_closest_center_A[dp_num][i+1]);
        
/*#ifdef  momnk_test
        othPaillierObj->decrypt(temp, E_W_A);
        gmp_printf("\nW: %Zd+%Zd ", temp, W_B);

        othPaillierObj->decrypt(temp, E_cmp_A[i]);
        gmp_printf("cmp: %Zd+%Zd ", temp, cmp_B[i]);

        othPaillierObj->decrypt(temp, E_P);
        gmp_printf("cmp: %Zd ", temp, E_P);

#endif*/

        read(E_P);
        
        XOR_shared_num(bit, E_W_A, E_closest_center_A[dp_num][i], E_P, W_B, closest_center_B[dp_num][i]);
        
#ifdef  momnk_test
        othPaillierObj->decrypt(temp, E_W_A);
        gmp_printf(" %Zd + %Zd *", temp, W_B);


        othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
        gmp_printf(" %Zd + %Zd = ", temp, closest_center_B[dp_num][i]);
        othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
        gmp_printf("%Zd ", temp);
      
        othPaillierObj->decrypt(temp, E_P);
        gmp_printf(" mul=%Zd ", temp);

#endif


#ifdef  momnk_test
        othPaillierObj->decrypt(temp, bit);
        gmp_printf(" XOR=%Zd ", temp);

#endif
        
        share_write(E_W_min_A, W_min_B, bit);
        
        neg_shared_num(temp, E_neg_A, neg_B, E_W_A, W_B);
        
#ifdef  momnk_test
        othPaillierObj->decrypt(temp, E_neg_A);
        othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
#endif
        
        read(E_P);
        
//        prod_shared_num(E_cmp[permutation[i]], E_W_min_A, E_neg_A, E_W_min_A, W_min_B, neg_B);
        prod_shared_num(E_closest_center[dp_num][i], E_W_min_A, E_neg_A, E_P, W_min_B, neg_B);
        
        share_write(E_closest_center_A[dp_num][i], closest_center_B[dp_num][i], E_closest_center[dp_num][i]);

        
#ifdef  momnk_test
        gmp_printf("\n%Zd+ ", closest_center_B[dp_num][i]);

        othPaillierObj->decrypt(temp, E_closest_center_A[dp_num][i]);
        gmp_printf("%Zd ", temp);
        
        othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
        gmp_printf("=%Zd ", temp);

#endif

    }
    
    //unpermute cmp[permutation[i]] which is cmp[i]
    
    //share cmp[i]
#ifdef  momnk_test

    cout<<"\nMOMN output :";

    for(long   int i =0; i<NUMCENTERS; ++i)
    {
        //share_write(E_cmp_A[i], cmp_B[i], E_cmp[i]);
        othPaillierObj->decrypt(temp, E_closest_center[dp_num][i]);
        gmp_printf("%Zd ", temp);

    }
#endif

    //share_write(E_cmp_A[i], cmp_B[i], E_cmp[permutation[i]]);
    
    mpz_clear(temp);
    mpz_clear(E_P);
    mpz_clear(E_Q);
    mpz_clear(r_1);
    mpz_clear(r_2);
    mpz_clear(E_W_A);
    mpz_clear(W_B);
    mpz_clear(E_neg_A);
    mpz_clear(neg_B);
    mpz_clear(E_W_min_A);
    mpz_clear(W_min_B);
    mpz_clear(E_W_min);
    mpz_clear(coin);
    mpz_clear(bit);
    
}


void kMeansProtocol::RVS(mpz_t E_rand, mpz_t E_rand_A, mpz_t rand_B, mpz_t E_QA, mpz_t QB)
{
#ifdef test_RVS
    cout<<"\nIn RVS"<<endl;
#endif
    
    mpz_t E_q[LValue], E_qA[LValue], qB[LValue], S, S_prev, E_cmp, E_cmp_A, cmp_B, temp, E_temp, E_w, E_wA, wB, E_vA, E_v, vB, aB, aB_prev;
    mpz_t E_eA[LValue], E_e[LValue], eB[LValue];
    long   int perm[LValue];
    
    mpz_init(S);
    mpz_init(S_prev);
    mpz_init(E_cmp);
    mpz_init(E_cmp_A);
    mpz_init(cmp_B);
    mpz_init(temp);
    mpz_init(E_temp);
    mpz_init(E_w);
    mpz_init(E_wA);
    mpz_init(wB);
    mpz_init(E_vA);
    mpz_init(E_v);
    mpz_init(vB);
    mpz_init(aB);
    mpz_init(aB_prev);



    for(long   int i =0; i<L; ++i)
    {
        mpz_init(E_q[i]);
        mpz_init(E_qA[i]);
        mpz_init(qB[i]);
        
        mpz_init(E_eA[i]);
        mpz_init(E_e[i]);
        mpz_init(eB[i]);
        perm[i]=i;
    }
    
     TBP(E_q, E_qA, qB, E_QA, QB);
    
    //reordering_L(perm);

#ifdef test_RVS
    cout<<endl<<"q :";
    for(long   int i =1; i<L; ++i)
    {
        othPaillierObj->decrypt(temp, E_q[i]);
        gmp_printf(" %Zd ", temp);
    }
#endif

    mpz_set(eB[perm[0]],qB[perm[0]]);
    mpz_set(E_eA[perm[0]],E_qA[perm[0]]);
    mpz_set(E_e[perm[0]], E_q[perm[0]]);

    mpz_set_ui(cmp_B,0);
    othPaillierObj->encrypt(E_cmp, cmp_B);
    othPaillierObj->encrypt(E_cmp_A, cmp_B);
    
    mpz_set_ui(temp, 0);
    othPaillierObj->encrypt(E_rand, temp);
 
#ifdef test_RVS
    cout<<endl<<"e: ";
#endif
    
    for(long   int i =1; i<L; ++i)
    {
        mpz_add(cmp_B, cmp_B, eB[perm[i-1]]);
        mpz_mul(E_cmp_A, E_cmp_A, E_eA[perm[i-1]]);
        
        read(temp);
        
        prod_shared_num(E_w, E_qA[perm[i]], E_cmp_A, temp, qB[perm[i]], cmp_B);
        mpz_invert(temp, E_w, Msqr);
        
        mpz_mul(E_e[i], E_q[i], temp);
        mpz_mod(E_e[i], E_e[i], Msqr);
        
        share_write(E_eA[i], eB[i], E_e[i]);
    
#ifdef test_RVS
        othPaillierObj->decrypt(temp, E_e[i]);
        gmp_printf(" %Zd ", temp);
#endif
    }

    //cout<<endl;
   
    //choose Si
    mpz_set_ui(S,0);
     mpz_set_ui(S_prev, 1);
    mpz_set_ui(aB, 0);
    mpz_set_ui(wB, 0);

    for(long   int i =0; i<L; ++i)
    {
        if(i>0)
        {

            mpz_mul_ui(S, S_prev, 2);
            mpz_set(S_prev, S);

            getRandomN(wB, i, S);

            mpz_mul(aB, aB_prev, qB[i-1]);

            if(i>1)
                mpz_mul_ui(aB, aB, 2);

        }

        
        read(E_wA);
        mpz_set_ui(temp, 0);
        othPaillierObj->encrypt(E_temp, temp);
        
        MOTN(E_cmp, E_cmp_A, cmp_B, E_wA, E_temp, wB, S);
            
        othPaillierObj->encrypt(E_temp, wB);
        mpz_mul(E_w, E_wA, E_temp);
        mpz_mod(E_w, E_w, Msqr);
        
        mpz_powm(E_temp, E_cmp, S, Msqr);
        mpz_invert(temp, E_temp, Msqr);
        mpz_mul(E_v, temp, E_w);
        mpz_mod(E_v,E_v,Msqr);
        
        share_write(E_vA, vB, E_v);
        
        mpz_set(aB_prev, aB);
        
        //comp final
        read(E_vA); //vA= vA+aA
        mpz_add(vB, vB, aB);
        mpz_mod(vB, vB, M);
    
        read(E_wA);//eA(vA+aA)
        
        prod_shared_num(temp, E_vA, E_eA[i], E_wA, vB, eB[i]);
        
        //othPaillierObj->decrypt(E_vA, temp);
        
        mpz_mul(E_rand, E_rand, temp);
        mpz_mod(E_rand, E_rand, Msqr);
        //othPaillierObj->decrypt(temp, E_rand);
        mpz_mod(E_rand, E_rand, Msqr);

    }
    
    share_write(E_rand_A, rand_B, E_rand);

    //othPaillierObj->decrypt(temp, E_rand);


    mpz_clear(S);
    mpz_clear(S_prev);
    mpz_clear(E_cmp);
    mpz_clear(E_cmp_A);
    mpz_clear(cmp_B);
    mpz_clear(temp);
    mpz_clear(E_temp);
    mpz_clear(E_w);
    mpz_clear(E_wA);
    mpz_clear(wB);
    mpz_clear(E_vA);
    mpz_clear(E_v);
    mpz_clear(vB);
    mpz_clear(aB);
    mpz_clear(aB_prev);

    for(long   int i =0; i<L; ++i)
    {
        mpz_clear(E_q[i]);
        mpz_clear(E_qA[i]);
        mpz_clear(qB[i]);
        
        mpz_clear(E_eA[i]);
        mpz_clear(E_e[i]);
        mpz_clear(eB[i]);
    }


}

//other sub-protocol -write
/*void kMeansProtocol::comp_distance(mpz_t E_dist, mpz_t E_dist_A, mpz_t dist_B, long   int dp_num, long   int cent_num) //correct but no needed
{
    mpz_t temp, E_x, E_y, E_z, z1_B, sq_sum_dp_B, sq_sum_cent_B;
    mpz_init(temp);
    mpz_init(E_x);
    mpz_init(E_y);
    mpz_init(E_z);
    mpz_init(z1_B);
    mpz_init(sq_sum_dp_B);
    mpz_init(sq_sum_cent_B);

    mpz_set_ui(z1_B,0);
    
    othPaillierObj->encrypt(E_x, z1_B);//CA.CB
    othPaillierObj->encrypt(E_y, z1_B);//dB.CA
    othPaillierObj->encrypt(E_z, z1_B);//dA.CB
    
    mpz_set_ui(sq_sum_dp_B, 0);//d^2
    mpz_set_ui(sq_sum_cent_B, 0);//CB^2


    for(long   int j=0; j<NUMATTRIBUTE; ++j)
    {
        mpz_addmul(sq_sum_dp_B, data_pt_B[dp_num][j], data_pt_B[dp_num][j]);//mpz_addmul: a=a+ b*b
        mpz_mod(sq_sum_dp_B, sq_sum_dp_B, M); //dB^2
        
        othPaillierObj->decrypt(temp, sq_sum_dp_B);

        mpz_addmul(sq_sum_cent_B, center_B[cent_num][j], center_B[cent_num][j]);
        mpz_mod(sq_sum_cent_B, sq_sum_cent_B, M); //CB^2
        
        othPaillierObj->decrypt(temp, sq_sum_cent_B);

        mpz_powm(temp, E_center_A[cent_num][j], center_B[cent_num][j], Msqr);
        mpz_mul(E_x, E_x, temp);
        mpz_mod(E_x, E_x, Msqr);//CA.CB
        othPaillierObj->decrypt(temp, E_x);

        mpz_powm(temp, E_center_A[cent_num][j], data_pt_B[dp_num][j], Msqr);
        mpz_mul(E_y, E_y, temp);
        mpz_mod(E_y, E_y, Msqr);//CA.dB
        othPaillierObj->decrypt(temp, E_y);

        mpz_powm(temp, E_data_pt_A[dp_num][j], center_B[cent_num][j], Msqr);
        mpz_mul(E_z, E_z, temp);
        mpz_mod(E_z, E_z, Msqr);//dA.CB
        othPaillierObj->decrypt(temp, E_z);

        mpz_addmul(z1_B, data_pt_B[dp_num][j], center_B[cent_num][j]);
        mpz_mod(z1_B, z1_B, M);//dB.CB

    }


    mpz_powm_ui(E_x, E_x, 2, Msqr); //2 CA.CB
    //share_write(E_x_A, x_B, E_x);
    othPaillierObj->decrypt(temp, E_x);

    mpz_powm_ui(E_y, E_y, 2, Msqr);
    mpz_invert(E_y, E_y, Msqr);//-2CA.dB
    //share_write(E_y_A, y_B, E_y);
    othPaillierObj->decrypt(temp, E_y);

    mpz_mul(E_dist, E_x, E_y);//2 CA.CB-2CA.dB
    othPaillierObj->decrypt(temp, E_dist);
    
    mpz_powm_ui(E_z, E_z, 2, Msqr);
    mpz_invert(E_z, E_z, Msqr);//-2dA.CB
    //share_write(E_z_A, z_B, E_z);
    othPaillierObj->decrypt(temp, E_z);

    mpz_mul(E_dist, E_dist, E_z);//2 CA.CB-2CA.dB-2dA.CB
    mpz_mod(E_dist, E_dist, Msqr);
    othPaillierObj->decrypt(temp, E_dist);

    mpz_mul_ui(z1_B, z1_B, -2);
    mpz_mod(z1_B, z1_B, M);//-2dB.CB
    
    mpz_add(dist_B, sq_sum_cent_B, sq_sum_dp_B);//dB^2+CB^2
    mpz_add(dist_B, dist_B, z1_B);
    mpz_mod(dist_B, dist_B, M); //dB^2+CB^2-2dB.CB
    
    read(E_dist_A);//dA^2+CA^2-2CA.DA
    mpz_mul(E_dist, E_dist, E_dist_A); //dA^2+CA^2-2CA.DA +2 CA.CB-2CA.dB-2dA.CB
    mpz_mod(E_dist, E_dist, Msqr);//dA^2+CA^2-2CA.DA +2 CA.CB-2CA.dB-2dA.CB
    othPaillierObj->decrypt(temp, E_dist);

    othPaillierObj->encrypt(temp, dist_B);
    mpz_mul(E_dist, E_dist, temp);//dA^2+CA^2-2CA.DA +dB^2+CB^2-2CB.DB  +2 CA.CB-2CA.dB-2dA.CB
    mpz_mod(E_dist, E_dist, Msqr);//dA^2+dB^2 +CA^2+CB^2+2 CA.CB -2CA.DA-2CA.DB -2CB.DA-2CB.DB
    
    share_write(E_dist_A, dist_B, E_dist);
    
    othPaillierObj->decrypt(temp, E_dist);
    gmp_printf("\nDist between %d-%d:%Zd ",dp_num, cent_num, temp);
    
    mpz_clear(temp);
    mpz_clear(E_x);
    mpz_clear(E_y);
    mpz_clear(E_z);
    mpz_clear(z1_B);
    mpz_clear(sq_sum_dp_B);
    mpz_clear(sq_sum_cent_B);
    
}*/

void kMeansProtocol::comp_distance_short(mpz_t E_sc_dist, mpz_t E_sc_dist_A, mpz_t sc_dist_B, long   int dp_num, long   int cent_num)
//dist scaled by n^2
{
    //cout<<"\nIn Comp dist"<<endl;

    mpz_t temp, E_x, E_x_A, x_B, E_xAsqr, E_y, E_z, sum;
    mpz_init(temp);
    mpz_init(E_x);
    mpz_init(E_y);
    mpz_init(E_z);
    mpz_init(E_x_A);
    mpz_init(x_B);
    mpz_init(E_xAsqr);
    mpz_init(sum);

    mpz_set_ui(temp,0);
    
    othPaillierObj->encrypt(E_sc_dist, temp);


    for(long   int j=0; j<NUMATTRIBUTE; ++j)
    {
        mpz_invert(E_y, E_center[cent_num][j], Msqr);//-c
        mpz_mul(E_x, E_data_pt[dp_num][j], E_y);//d-c
        mpz_mod(E_x, E_x, Msqr);
        mpz_powm_ui(E_x, E_x, NUMDATAPOINT, Msqr);//n(d-c)
        
        share_write(E_x_A, x_B, E_x);
        
        read(E_xAsqr);
        
        prod_shared_num(E_z, E_x_A, E_x_A, E_xAsqr, x_B, x_B);//n^2(d-c)^2
        mpz_mul(E_sc_dist, E_sc_dist, E_z);
        mpz_mod(E_sc_dist, E_sc_dist, Msqr);
    }
   
    share_write(E_sc_dist_A, sc_dist_B, E_sc_dist);
    
#ifdef test_dist
    othPaillierObj->decrypt(temp, E_sc_dist);
    gmp_printf("\nDist between %d-%d:%Zd ",dp_num, cent_num, temp);
#endif
    
    mpz_clear(temp);
    mpz_clear(E_x);
    mpz_clear(E_y);
    mpz_clear(E_z);
    mpz_clear(E_x_A);
    mpz_clear(x_B);
    mpz_clear(E_xAsqr);
    mpz_clear(sum);
    
}


void kMeansProtocol::TBP(mpz_t E_q[LValue], mpz_t E_qA[LValue], mpz_t qB[LValue], mpz_t E_QA, mpz_t QB)
{
#ifdef test_TBP
    cout<<"\nIn TBP"<<endl;
#endif
    mpz_t S, E_mul, E_cmp, E_cmp_A, cmp_B, E_neg_cmp, E_neg_cmp_A, neg_cmp_B, d, E_XOR1,E_XOR1_A, XOR1_B, E_XOR2,E_XOR2_A, XOR2_B;
    mpz_t E_carry_prev, E_carry_prev_A, carry_prev_B, E_carry1, E_carry1_A, carry1_B, E_carry2, E_carry2_A, carry2_B, E_carry3, E_carry3_A, carry3_B, temp;
    mpz_t E_Ex2carry_prev, E_Ex2carry_prev_A, Ex2carry_prev_B, E_Ex2carry1, E_Ex2carry1_A, Ex2carry1_B, E_Ex2carry2, E_Ex2carry2_A, Ex2carry2_B, E_Ex2carry3, E_Ex2carry3_A, Ex2carry3_B;
    //mpz_t /*E_QA_bin[LValue], QB_bin[LValue],*/ E_XOR1[LValue],E_XOR1_A[LValue], XOR1_B[LValue], E_XOR2[LValue],E_XOR2_A[LValue], XOR2_B[LValue] /*, d[LValue]*/;
    mpz_init(d);
    mpz_init(S);
    mpz_init(E_carry1);
    mpz_init(E_carry1_A);
    mpz_init(carry1_B);
    mpz_init(E_carry2);
    mpz_init(E_carry2_A);
    mpz_init(carry2_B);
    mpz_init(E_carry3);
    mpz_init(E_carry3_A);
    mpz_init(carry3_B);
    mpz_init(E_carry_prev);
    mpz_init(E_carry_prev_A);
    mpz_init(carry_prev_B);
    mpz_init(temp);
    
    mpz_init(E_Ex2carry1);
    mpz_init(E_Ex2carry1_A);
    mpz_init(Ex2carry1_B);
    mpz_init(E_Ex2carry2);
    mpz_init(E_Ex2carry2_A);
    mpz_init(Ex2carry2_B);
    mpz_init(E_Ex2carry3);
    mpz_init(E_Ex2carry3_A);
    mpz_init(Ex2carry3_B);
    mpz_init(E_Ex2carry_prev);
    mpz_init(E_Ex2carry_prev_A);
    mpz_init(Ex2carry_prev_B);
    
    mpz_init(E_mul);
    mpz_init(E_cmp);
    mpz_init(E_cmp_A);
    mpz_init(cmp_B);
    mpz_init(E_neg_cmp);
    mpz_init(E_neg_cmp_A);
    mpz_init(neg_cmp_B);

    mpz_init(E_XOR1);
    mpz_init(E_XOR1_A);
    mpz_init(XOR1_B);
    mpz_init(E_XOR2);
    mpz_init(E_XOR2_A);
    mpz_init(XOR2_B);
    
    mpz_set_ui(Ex2carry_prev_B, 0);
    othPaillierObj->encrypt(E_Ex2carry_prev, Ex2carry_prev_B);
    mpz_set(E_Ex2carry_prev_A, E_Ex2carry_prev);
    
    mpz_set_ui(carry_prev_B, 0);
    othPaillierObj->encrypt(E_carry_prev, carry_prev_B);
    mpz_set(E_carry_prev_A, E_carry_prev);
    
    //binary of QB
    mpz_set(S, QB);
    for(long   int j=0; j<LValue; j++)
    {
     //mpz_init(E_QA_bin[j]);
        //mpz_init(QB_bin[j]);
       //mpz_init(E_XOR1[j]);
       // mpz_init(E_XOR1_A[j]);
        //mpz_init(XOR1_B[j]);
        //mpz_init(E_XOR2[j]);
        //mpz_init(E_XOR2_A[j]);
        //mpz_init(XOR2_B[j]);

        //mpz_init(d[j]);

        mpz_mod_ui(qB[j], S, 2);
        mpz_div_ui(S, S, 2);
    }
    
    //run FM2NP
    mpz_sub(S, M, QB);
    mpz_set_ui(carry1_B, 0);
    othPaillierObj->encrypt(E_carry1_A, carry1_B);
    MOTN(E_cmp, E_cmp_A, cmp_B, E_QA, E_carry1_A, carry1_B, S);
    neg_shared_num(E_neg_cmp, E_neg_cmp_A, neg_cmp_B, E_cmp_A, cmp_B);
    //cmp=0 if XA+XB<M
    
    //d = 2^L+1 -M
    //2^L-M
    mpz_set_ui(S, 2);
    mpz_pow_ui(S, S, L);
    mpz_sub(S, S, M);
    
    /*binary of D
    for(long   int j=0; j<L; j++)
    {
        mpz_mod_ui(d[j], temp, 2);
        mpz_div_ui(temp, temp, 2);
    }*/
    

    for(long   int i=0; i<LValue; i++)
    {
        read(E_qA[i]);
        
#ifdef test_TBP
        othPaillierObj->decrypt(temp, E_qA[i]);
        gmp_printf("\n%Zd x %Zd", temp, qB[i]);
#endif
        
        XOR(E_XOR1, E_qA[i], qB[i]);
        share_write(E_XOR1_A, XOR1_B, E_XOR1);

        //othPaillierObj->decrypt(temp, E_XOR1);
        
        read(E_mul);//XOR1_A[i]*E_carry_prev_A
        //xor with prev carry = a+b
        XOR_shared_num(E_XOR1, E_XOR1_A, E_carry_prev_A, E_mul, XOR1_B, carry_prev_B);
        share_write(E_XOR1_A, XOR1_B, E_XOR1);
        //othPaillierObj->decrypt(temp, E_XOR1);

        //curr carry =a.b or b.c or c.a = a.b or c(a or b)
        scalar_prod(E_carry1, E_qA[i], qB[i]); //a.b
        share_write(E_carry1_A, carry1_B, E_carry1);
        //othPaillierObj->decrypt(temp, E_carry1);

        OR(E_carry2, E_qA[i], qB[i]); //a or b
        share_write(E_carry2_A, carry2_B, E_carry2);
        //othPaillierObj->decrypt(temp, E_carry2);

        read(E_mul);
        prod_shared_num(E_carry3, E_carry_prev_A, E_carry2_A, E_mul, carry_prev_B, carry2_B); //c(a or b)
        share_write(E_carry3_A, carry3_B, E_carry3);
        //othPaillierObj->decrypt(temp, E_carry3);

        read(E_mul);
        OR_shared_num(E_carry_prev, E_carry1_A, E_carry3_A, E_mul, carry1_B, carry3_B);
        share_write(E_carry_prev_A, carry_prev_B, E_carry_prev);
        //othPaillierObj->decrypt(temp, E_carry_prev);


    
    //Now compute a + b + d =xor+d

        mpz_mod_ui(d, S, 2);
        mpz_div_ui(S, S, 2);
        
        //a ex b ex d
        XOR(E_XOR2, E_XOR1, d);
        share_write(E_XOR2_A, XOR2_B, E_XOR2);
        read(E_mul);
        XOR_shared_num(E_XOR2, E_XOR2_A, E_Ex2carry_prev_A, E_mul, XOR2_B, Ex2carry_prev_B);
        share_write(E_XOR2_A, XOR2_B, E_XOR2);
                
        
        //compute carry to add in next iteration
        // compute (a ex b)*d
        //curr carry =a.b or b.c or c.a = a.b or c(a or b)
        scalar_prod(E_Ex2carry1, E_XOR1, d); //a.b
        share_write(E_Ex2carry1_A, Ex2carry1_B, E_Ex2carry1);
        OR(E_Ex2carry2, E_XOR1, d); //a or b
        share_write(E_Ex2carry2_A, Ex2carry2_B, E_Ex2carry2);
        read(E_mul);
        prod_shared_num(E_Ex2carry3, E_Ex2carry_prev_A, E_Ex2carry2_A, E_mul, Ex2carry_prev_B, Ex2carry2_B); //c(a or b)
        share_write(E_Ex2carry3_A, Ex2carry3_B, E_Ex2carry3);
        read(E_mul);
        OR_shared_num(E_Ex2carry_prev, E_Ex2carry1_A, E_Ex2carry3_A, E_mul, Ex2carry1_B, Ex2carry3_B);
        share_write(E_Ex2carry_prev_A, Ex2carry_prev_B, E_Ex2carry_prev);

    
    //now we have a+b and a+b+d

        read(E_mul);
        prod_shared_num(E_carry1, E_neg_cmp_A, E_XOR1_A, E_mul, neg_cmp_B, XOR1_B);
        
        read(E_mul);
        prod_shared_num(E_carry2, E_cmp_A, E_XOR2_A, E_mul, cmp_B, XOR2_B);
        
        mpz_mul(E_q[i], E_carry1, E_carry2);
        mpz_mod(E_q[i], E_q[i], Msqr);
        
        share_write(E_qA[i], qB[i], E_q[i]);
        
    }
    
    
    mpz_clear(d);
    mpz_clear(S);
    mpz_clear(E_carry1);
    mpz_clear(E_carry1_A);
    mpz_clear(carry1_B);
    mpz_clear(E_carry2);
    mpz_clear(E_carry2_A);
    mpz_clear(carry2_B);
    mpz_clear(E_carry3);
    mpz_clear(E_carry3_A);
    mpz_clear(carry3_B);
    mpz_clear(E_carry_prev);
    mpz_clear(E_carry_prev_A);
    mpz_clear(carry_prev_B);
    
    mpz_clear(E_Ex2carry1);
    mpz_clear(E_Ex2carry1_A);
    mpz_clear(Ex2carry1_B);
    mpz_clear(E_Ex2carry2);
    mpz_clear(E_Ex2carry2_A);
    mpz_clear(Ex2carry2_B);
    mpz_clear(E_Ex2carry3);
    mpz_clear(E_Ex2carry3_A);
    mpz_clear(Ex2carry3_B);
    mpz_clear(E_Ex2carry_prev);
    mpz_clear(E_Ex2carry_prev_A);
    mpz_clear(Ex2carry_prev_B);
    
    mpz_clear(E_mul);
    mpz_clear(E_cmp);
    mpz_clear(E_cmp_A);
    mpz_clear(cmp_B);
    mpz_clear(E_neg_cmp);
    mpz_clear(E_neg_cmp_A);
    mpz_clear(neg_cmp_B);

    mpz_clear(E_XOR1);
    mpz_clear(E_XOR1_A);
    mpz_clear(XOR1_B);
    mpz_clear(E_XOR2);
    mpz_clear(E_XOR2_A);
    mpz_clear(XOR2_B);
    
}


void kMeansProtocol::reordering_L( long    int perm[LValue])
{
        cout<<endl<<"In re-order:";
   /*
        mpz_t U, temp1, r, a[LValue], two, min, max;
        mpz_init(U);
        mpz_init(temp1);
        mpz_init(r);
        mpz_init_set_ui(two, 2);
        mpz_init_set_ui(min, 0);
        mpz_init_set_ui(max, 0);

        long int u[LValue];
        bool first =true;
        
        for(int i =0; i<L; ++i)
            mpz_init_set_ui(a[i],0);
        
        for(int i =0; i<L; ++i)
        {
            first =true;
            if(i==0)
            {
                mpz_pow_ui(U, two, K);
                mpz_sub_ui(U, U, 1);
            }
            else
            {
                mpz_pow_ui(temp1, two, perm[i-1]);
                mpz_sub(U, U, temp1);
            }
            
            //binary of U
            mpz_set(temp1, U);
            gmp_printf("\n%Zd :", U);
            
            for(int j=0; j<L; j++)
            {
                mpz_mod_ui(r, temp1, 2);
                u[j] = (int) mpz_get_ui(r);
                mpz_div_ui(temp1, temp1, 2);
                cout<<u[j];
            }
            
            getRandomN(r, L, U);
            gmp_printf("\nr=%Zd", r);
            for(int j =0; j<L; ++j)
            {
                if(j>0)
                {
                    mpz_pow_ui(temp1, two, (j-1));
                    mpz_mul_ui(temp1, temp1, u[j-1]);
                    mpz_add(a[j], a[j-1], temp1);
                }
                else
                    mpz_set_ui(a[j], 0);
                gmp_printf("\na[j]=%Zd", a[j]);
                
                if(u[j] == 1)
                {
                    
                    mpz_set(min, a[j]);
                    mpz_pow_ui(temp1, two, j);
                    mpz_add(max, min, temp1);

                    if( first == 1)
                    {
                        mpz_set_ui(min, 0);
                        first =0;
                    }
                    else
                        mpz_add_ui(min, min, 1);
                    
                    gmp_printf(" %Zd:%Zd", min, max);
                    
                    if( mpz_cmp(r, min) >= 0 && mpz_cmp(r, max) <= 0)
                    {
                        perm[i] = j;
                        break;
                    }
                    
                }
                
                //if Uj ==0, then that intervals is already taken
            }
            //cout<<" "<<i<<":"<<sigma[i];
        }
        cout<<endl;
        for(int j =0; j<L; ++j)
            cout<<" "<<j<<":"<<perm[j];
        
        mpz_clear(U);
        mpz_clear(temp1);
        mpz_clear(r);
        mpz_clear(min);
        mpz_clear(max);
        for(int i =0; i<L; ++i)
            mpz_clear(a[i]);*/
        

}


void kMeansProtocol::choose_center( mpz_t E_randA, mpz_t randB, long   int center_num)
{
    mpz_t E_Z[NUMDATAPOINT], E_ZA[NUMDATAPOINT], ZB[NUMDATAPOINT],E_cmp[NUMDATAPOINT], E_cmp_A[NUMDATAPOINT], cmp_B[NUMDATAPOINT], temp, enc, temp2, mul, sum;
    mpz_init(enc);
    mpz_init(temp2);
    mpz_init_set_ui(temp,0);
    mpz_init(mul);
    mpz_init(sum);

    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        mpz_init(E_Z[i]);
        mpz_init(E_ZA[i]);
        mpz_init(ZB[i]);
        mpz_init(E_cmp[i]);
        mpz_init(E_cmp_A[i]);
        mpz_init(cmp_B[i]);

    }
    
    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        //write(temp);
        read(E_ZA[i]);
    }
    
    if(center_num==0)
    {
        mpz_mul_ui(temp, min_sc_dist_B[0], NUMDATAPOINT);
        mpz_add(ZB[0], avg_min_sc_dist_B, temp);
        mpz_mod(ZB[0], ZB[0], M);
        mpz_sub(ZB[0], ZB[0], randB);
        mpz_mod(ZB[0], ZB[0], M);
        

        for(long   int i =1; i<NUMDATAPOINT; ++i)
        {
            //ZB=avg_min_sc_dist_B + nsc_dist_B[i]
            mpz_mul_ui(temp, min_sc_dist_B[0], NUMDATAPOINT);
            mpz_add(temp, temp, avg_min_sc_dist_B);
            mpz_mod(ZB[i], ZB[i], M);
            mpz_add(ZB[i], ZB[i-1], temp);
            mpz_sub(ZB[i], ZB[i], randB);
            mpz_mod(ZB[i], ZB[i], M);
        }
    
    }
    else
    {
        mpz_sub(ZB[0], min_sc_dist_B[0], randB);
        mpz_mod(ZB[0], ZB[0], M);

        for(long   int i =1; i<NUMDATAPOINT; ++i)
        {
            //ZB= sc_dist_B[i]
            mpz_add(ZB[i], ZB[i-1], min_sc_dist_B[i]);
            mpz_mod(ZB[i], ZB[i], M);
        }
            
    }
    
    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        othPaillierObj->encrypt(temp, ZB[i]);
        mpz_mul(E_Z[i], E_ZA[i], temp);
        mpz_mod(E_Z[i],E_Z[i],Msqr);
    }
    

    MOMN_n(E_cmp, E_cmp_A, cmp_B, E_Z, E_ZA, ZB);
    
    cout<<"\ndatapt chosen as center :";
    for(long   int j=0; j < NUMDATAPOINT; ++j)
    {
            othPaillierObj->decrypt(temp, E_cmp[j]);
            gmp_printf("%Zd ",temp);
    }
    
    
    for(long   int j =0; j<NUMATTRIBUTE; ++j)
    {
        //cout<<"\nAttribute "<<j<<endl;
        mpz_set_ui(temp,  0);
        othPaillierObj->encrypt(sum, temp);
        
        for(long   int i =0; i<NUMDATAPOINT; ++i)
        {
            //othPaillierObj->decrypt(mul, enc);

            //mpz_mul(temp2, cmp_B[i], data_pt_B[i][j]);
            //mpz_add(mul, mul, temp2);
            //mpz_mod(mul, mul, M);
            
            mpz_mul(temp, cmp_B[i], data_pt_B[i][j]);
            othPaillierObj->encrypt(temp2, temp);
            mpz_mul(sum, sum, temp2);
            mpz_mod(sum, sum, Msqr);
            //othPaillierObj->decrypt(temp, enc);

            
            mpz_powm(temp2, E_data_pt_A[i][j], cmp_B[i], Msqr);
            
            //othPaillierObj->decrypt(temp, temp2);
            //mpz_add(mul, mul, temp);
            //mpz_mod(mul, mul, M);
            
            mpz_mul(sum, sum, temp2);
            mpz_mod(sum, sum, Msqr);
            //othPaillierObj->decrypt(temp, sum);
            //othPaillierObj->encrypt(sum, mul);

            mpz_powm(temp2, E_cmp_A[i], data_pt_B[i][j], Msqr);
            
            //othPaillierObj->decrypt(temp, temp2);
            //mpz_add(mul, mul, temp);
            //mpz_mod(mul, mul, M);

            mpz_mul(sum, sum, temp2);
            mpz_mod(sum, sum, Msqr);
            //othPaillierObj->decrypt(temp, sum);
            
            //mpz_add(sum, sum, mul);
            //gmp_printf("%Zd ",temp);

        }
  
        read(E_center[center_num][j]);
        
        mpz_mul(E_center[center_num][j], E_center[center_num][j], sum);
        mpz_mod(E_center[center_num][j], E_center[center_num][j], Msqr);
        
        share_write(E_center_A[center_num][j], center_B[center_num][j], E_center[center_num][j]);
        
#ifdef ccp_test
        othPaillierObj->decrypt(temp, E_center_A[center_num][j]);
        gmp_printf(" \n %Zd +%Zd=",temp, center_B[center_num][j]);
            othPaillierObj->decrypt(temp, E_center[center_num][j]);
            gmp_printf("%Zd ",temp);
#endif
        
    }
    
    mpz_clear(enc);

    mpz_clear(temp);
    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        mpz_clear(E_Z[i]);
        mpz_clear(E_ZA[i]);
        mpz_clear(ZB[i]);
        mpz_clear(E_cmp[i]);
        mpz_clear(E_cmp_A[i]);
        mpz_clear(cmp_B[i]);

    }
    
}

//Scaling -write
void kMeansProtocol::scaling_fact()
{
    mpz_t temp, E_num_prod_A[NUMCENTERS], num_prod_B[NUMCENTERS], E_num_prod[NUMCENTERS];
    mpz_init(temp);
    
    for(long   int j =0; j<NUMCENTERS; ++j)
    {
        mpz_init(E_num_prod_A[j]);
        mpz_init(num_prod_B[j]);
        mpz_init(E_num_prod[j]);

        mpz_set_ui(alpha_B[j], 1);
        mpz_set_ui(alpha_B[j], 1);

        mpz_set_ui(temp, 0);
        othPaillierObj->encrypt(E_alpha_A[j], temp);
        othPaillierObj->encrypt(E_beta_A[j], temp);

        read(E_num_prod_A[j]);
        
    }
    for(long   int j =0; j<NUMCENTERS; ++j)
    {
        prod_shared_num(E_num_prod[j], E_num_dp_new_A[j], E_num_dp_A[j], E_num_prod_A[j], num_dp_new_B[j], num_dp_B[j]);
        
        /*
        othPaillierObj->decrypt(temp, E_num_dp_new_A[j]);
        gmp_printf("\n%Zd+%Zd",temp, num_dp_new_B[j]);
        othPaillierObj->decrypt(temp, E_num_dp_new[j]);
        gmp_printf("=%Zd ",temp);
        
        othPaillierObj->decrypt(temp, E_num_dp_A[j]);
        gmp_printf("\n%Zd+%Zd",temp, num_dp_B[j]);
        othPaillierObj->decrypt(temp, E_num_dp[j]);
        gmp_printf("=%Zd ",temp);
        */
        
        //othPaillierObj->decrypt(temp, E_num_prod[j]);

        share_write(E_num_prod_A[j], num_prod_B[j], E_num_prod[j]);
    }
    
    for(long   int j =0; j<NUMCENTERS; ++j)
    {
        for(long   int i =0; i<NUMCENTERS && i!=j ; ++i)
        {
            read(temp);
            prod_shared_num(E_alpha[j], E_alpha_A[j], E_num_dp_new_A[i], temp, alpha_B[j], num_dp_new_B[i]);
            
            //othPaillierObj->decrypt(temp, E_alpha[j]);

            read(temp);
            prod_shared_num(E_beta[j], E_beta_A[j], E_num_prod_A[i], temp, beta_B[j], num_prod_B[i]);
            
            //othPaillierObj->decrypt(temp, E_beta[j]);

            share_write(E_alpha_A[j], alpha_B[j], E_alpha[j]);
            share_write(E_beta_A[j], beta_B[j], E_beta[j]);
        }
    }
    
    read(temp);
    prod_shared_num(E_G_alpha, E_alpha_A[0], E_num_dp_new_A[0], temp, alpha_B[0], num_dp_new_B[0]);
    //othPaillierObj->decrypt(temp, E_G_alpha);

    read(temp);
    prod_shared_num(E_G_beta, E_beta_A[0], E_num_prod_A[0], temp, beta_B[0], num_prod_B[0]);
    //othPaillierObj->decrypt(temp, E_G_alpha);

    share_write(E_G_alpha_A, G_alpha_B, E_G_alpha);
    share_write(E_G_alpha_A, G_alpha_B, E_G_alpha);
    
    mpz_clear(temp);
    for(long   int j =0; j<NUMCENTERS; ++j)
    {
        mpz_clear(E_num_prod_A[j]);
        mpz_clear(num_prod_B[j]);
        mpz_clear(E_num_prod[j]);
    }
    
}



void kMeansProtocol::comp_scaled_dist(long   int dp_num, long   int cent_num )
{
    
    mpz_t E_mul1, E_mul2, E_w, E_w_A, w_B, x, y, E_wsqr, temp1, temp2;
    mpz_init(E_mul1);
    mpz_init(E_mul2);
    mpz_init(E_w);
    mpz_init(E_w_A);
    mpz_init(w_B);
    mpz_init(x);
    mpz_init(y);
    mpz_init(E_wsqr);
    mpz_init(temp1);
    mpz_init(temp2);


    for(long   int i =0; i<NUMATTRIBUTE; ++i)
    {
        read(E_mul1);
        read(E_mul2);
#ifdef dist_test
    othPaillierObj->decrypt(temp1, E_alpha_A[cent_num]);
    //gmp_printf("\n%Zd+%Zd",temp1,alpha_B[cent_num]);
        othPaillierObj->decrypt(temp1, E_center_A[cent_num][i]);
        //gmp_printf("\n%Zd+%Zd",temp1,center_B[cent_num][i]);
#endif
        prod_shared_num(x, E_alpha_A[cent_num], E_center_A[cent_num][i], E_mul1, alpha_B[cent_num], center_B[cent_num][i]);
        prod_shared_num(y, E_G_alpha_A, E_data_pt_A[dp_num][i], E_mul2, G_alpha_B, data_pt_B[dp_num][i]);
        
#ifdef dist_test
        othPaillierObj->decrypt(temp1, x);
        othPaillierObj->decrypt(temp1, y);

#endif
        mpz_invert(y, y, Msqr);
        mpz_mul(E_w, x, y);
        mpz_mod(E_w, E_w, Msqr);
        
        share_write(E_w_A, w_B, E_w);
        read(E_mul1);
        
        prod_shared_num(E_wsqr, E_w_A, E_w_A, E_mul1, w_B, w_B);
        
#ifdef dist_test
        othPaillierObj->decrypt(temp1, E_wsqr);
#endif
        
        mpz_mul(E_sc_dist[dp_num][cent_num], E_sc_dist[dp_num][cent_num], E_wsqr);
        mpz_mod(E_sc_dist[dp_num][cent_num],E_sc_dist[dp_num][cent_num],Msqr);
        
#ifdef dist_test
        othPaillierObj->decrypt(temp1, E_sc_dist[dp_num][cent_num]);
#endif
    }
    
    share_write( E_sc_dist_A[dp_num][cent_num],  sc_dist_B[dp_num][cent_num], E_sc_dist[dp_num][cent_num]);
#ifdef dist_test
    othPaillierObj->decrypt(E_mul1, E_sc_dist_A[dp_num][cent_num]);
    gmp_printf("\ndist %d-%d: %Zd+%Zd=",dp_num, cent_num, E_mul1, sc_dist_B[dp_num][cent_num]);
    othPaillierObj->decrypt(E_mul1, E_sc_dist[dp_num][cent_num]);
    gmp_printf("\n%Zd",E_mul1);
#endif
    
    mpz_clear(E_mul1);
    mpz_clear(E_mul2);
    mpz_clear(E_w);
    mpz_clear(E_w_A);
    mpz_clear(w_B);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(E_wsqr);
}


void kMeansProtocol::termination_condition(bool *b) //b=0 terminate: L<=R i.e. dist <= epsilon
{
    mpz_t E_mul1, E_mul2, E_w, E_w_A, w_B, x, y, E_wsqr, E_v, E_v_A, v_B, E_L, E_L_A, L_B, E_R, E_R_A, R_B, E_betasqr, E_betasqr_A, betasqr_B, temp;
    mpz_init(E_mul1);
    mpz_init(E_mul2);
    mpz_init(E_w);
    mpz_init(E_w_A);
    mpz_init(w_B);
    mpz_init(x);
    mpz_init(y);
    mpz_init(E_wsqr);
    mpz_init(E_v);
    mpz_init(E_v_A);
    mpz_init(v_B);
    mpz_init(E_L);
    mpz_init(E_L_A);
    mpz_init(L_B);
    mpz_init(E_R);
    mpz_init(R_B);
    mpz_init(E_R_A);
    mpz_init(E_betasqr);
    mpz_init(E_betasqr_A);
    mpz_init(betasqr_B);
    mpz_init_set_ui(temp,0);
    
    othPaillierObj->encrypt(E_v, temp);
    mpz_set_ui(temp, 1);
    othPaillierObj->encrypt(E_L, temp);
    othPaillierObj->encrypt(E_R, temp);

    for(long   int j =0; j<NUMCENTERS; ++j)
    {
        for(long   int i =0; i<NUMATTRIBUTE; ++i)
        {
            
            read(E_mul1);
            read(E_mul2);
            
           /* othPaillierObj->decrypt(temp, E_num_dp_new_A[j]);
            gmp_printf("\n%Zd+%Zd",temp, num_dp_new_B[j]);
            othPaillierObj->decrypt(temp, E_num_dp_new[j]);
            gmp_printf("=%Zd ",temp);
     
            othPaillierObj->decrypt(temp, E_center_A[j][i]);
            gmp_printf("* %Zd+%Zd",temp, center_B[j][i]);
            othPaillierObj->decrypt(temp, E_center[j][i]);
            gmp_printf("=%Zd ",temp);
            
            othPaillierObj->decrypt(temp, E_mul1);
            gmp_printf(" mul=%Zd ",temp);*/
            
            prod_shared_num(x, E_num_dp_new_A[j], E_center_A[j][i], E_mul1, num_dp_new_B[j], center_B[j][i]);
            othPaillierObj->decrypt(temp, x);

            
  /*          othPaillierObj->decrypt(temp, E_num_dp_A[j]);
            gmp_printf("\n%Zd+%Zd",temp, num_dp_B[j]);
            othPaillierObj->decrypt(temp, E_num_dp[j]);
            gmp_printf("=%Zd ",temp);
     
            othPaillierObj->decrypt(temp, E_center_new_A[j][i]);
            gmp_printf("* %Zd+%Zd",temp, center_new_B[j][i]);
            othPaillierObj->decrypt(temp, E_center_new[j][i]);
            gmp_printf("=%Zd ",temp);*/
            
            
            prod_shared_num(y, E_num_dp_A[j], E_center_new_A[j][i], E_mul2, num_dp_B[j], center_new_B[j][i]);
            
            
            //othPaillierObj->decrypt(temp, y);

            mpz_invert(y, y, Msqr);
            mpz_mul(E_w, x, y);
            mpz_mod(E_w, E_w, Msqr);
            
            //othPaillierObj->decrypt(temp, E_w);

            share_write(E_w_A, w_B, E_w);
        
            read(E_mul1);
            
            prod_shared_num(E_wsqr, E_w_A, E_w_A, E_mul1, w_B, w_B);
            mpz_mul(E_v, E_v, E_wsqr);
            mpz_mod(E_v, E_v,Msqr);
            
            //othPaillierObj->decrypt(temp, E_v);

     
        }
        share_write(E_v_A, v_B, E_v);
    
        read(E_mul1);
        
        prod_shared_num(E_betasqr, E_beta_A[j], E_beta_A[j], E_mul1, beta_B[j], beta_B[j]);
        othPaillierObj->decrypt(temp, E_betasqr);

        share_write(E_betasqr_A, betasqr_B, E_betasqr);
        
        read(E_mul1);
        prod_shared_num(E_w, E_betasqr_A, E_v_A, E_mul1, betasqr_B, v_B);
        
        //othPaillierObj->decrypt(temp, E_w);

        mpz_mul(E_L, E_L, E_w);
        mpz_mod(E_L, E_L, Msqr);
    }
    
    //othPaillierObj->decrypt(temp, E_L);

    share_write(E_L_A, L_B, E_L);
    
    read(E_mul1);
    prod_shared_num(E_betasqr, E_G_beta_A, E_G_beta_A, E_mul1, G_beta_B, G_beta_B);

    mpz_powm_ui(E_R, E_betasqr, EPSILON, Msqr);
    
    share_write(E_R_A, R_B, E_R);
    

    MOTN_termination(b, E_L_A, E_R_A, L_B, R_B);
    cout<<"\nTerm cond: "<<*b<<endl;

#ifdef test_Termination
    othPaillierObj->decrypt(E_mul1, E_L);
    othPaillierObj->decrypt(E_mul2, E_R);

    gmp_printf(" -- %Zd Vs %Zd",E_mul1,E_mul2);
#endif
    
    if(*b==1)
    {
        
        for(long   int j =0; j<NUMCENTERS; ++j)
        {
            for(long   int i =0; i<NUMATTRIBUTE; ++i)
            {
                mpz_set(center_B[j][i], center_new_B[j][i]);
                mpz_set(E_center_A[j][i], E_center_new_A[j][i]);
                mpz_set(E_center[j][i], E_center_new[j][i]);
                
                mpz_set(sum_dp_B[j][i], sum_dp_new_B[j][i]);
                mpz_set(E_sum_dp_A[j][i], E_sum_dp_new_A[j][i]);
                mpz_set(E_sum_dp[j][i], E_sum_dp_new[j][i]);

            }
            mpz_set(num_dp_B[j], num_dp_new_B[j]);
            mpz_set(E_num_dp_A[j], E_num_dp_new_A[j]);
            mpz_set(E_num_dp[j], E_num_dp_new[j]);
            
        }
    }
    
    mpz_clear(E_mul1);
    mpz_clear(E_mul2);
    mpz_clear(E_w);
    mpz_clear(E_w_A);
    mpz_clear(w_B);
    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(E_wsqr);
    mpz_clear(E_v);
    mpz_clear(E_v_A);
    mpz_clear(v_B);
    mpz_clear(E_L);
    mpz_clear(E_L_A);
    mpz_clear(L_B);
    mpz_clear(E_R);
    mpz_clear(R_B);
    mpz_clear(E_R_A);
    mpz_clear(E_betasqr);
    mpz_clear(E_betasqr_A);
    mpz_clear(betasqr_B);
    
    
}
