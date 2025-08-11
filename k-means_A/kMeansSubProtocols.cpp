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

//static long int numOfEnc=0;

//subprotocols

//Supporting -NO write

void kMeansProtocol::encrypt_write(mpz_t value_A)
{
    char * strValue;
    mpz_t ciph;
    mpz_init(ciph);
    
    mpz_mod(value_A, value_A, M);

    //static time_point<high_resolution_clock> start = high_resolution_clock::now();
    //numOfEnc++;
    
    myPaillierObj->encrypt(ciph, value_A);
    
    //time_point<high_resolution_clock> stop = high_resolution_clock::now();
    //auto duration_sec = duration_cast<microseconds>(stop-start);
    //cout<<"\ntime for "<< numOfEnc<<" enc: "<<duration_sec.count()<<endl;
    
    strValue =  mpz_get_str(NULL, 10, ciph);
    
    connWithBob->writeDataWithLength(&strValue);
    free(strValue);
    
    mpz_clear(ciph);
}


void kMeansProtocol::read_decrypt(mpz_t value_A)
{
    char * strValue;
    mpz_t ciph;
    mpz_init(ciph);

    connWithBob->readDataWithLength(&strValue);
    mpz_set_str(ciph, strValue, 10);
    free(strValue);
    
    myPaillierObj->decrypt(value_A, ciph);
    
    mpz_clear(ciph);
}

void kMeansProtocol::neg_num_A(mpz_t op, mpz_t XA)
{
    /*
     mpz_neg(op_B, XB);
     mpz_add_ui(op_B, op_B, 1); //1-XB
     mpz_mod(op_B,op_B, M);
     
     mpz_invert(E_op_A, E_XA, Msqr);//-XA
     mpz_mod(E_op_A, E_op_A, Msqr);
     */
    
    mpz_sub(op, M, XA);//M-XA =XA
    mpz_mod(op, op, M);
}


//sub-proto -write
void kMeansProtocol::MOTN(mpz_t cmp_A)
{
    mpz_t P,Q, bit;
    mpz_init(P);
    mpz_init(Q);
    mpz_init(bit);

    read_decrypt(P);
    read_decrypt(Q);
    
    short  int cmp = mpz_cmp(P,Q); //Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
    if(cmp >0) // P>Q
        mpz_set_ui(bit, 1);
    else mpz_set_ui(bit, 0); //P<=Q
    
    encrypt_write(bit);
        
    mpz_clear(P);
    mpz_clear(Q);
    mpz_clear(bit);
    
    read_decrypt(cmp_A);
    
}

void kMeansProtocol::MOTN_value(mpz_t min_A, mpz_t XA, mpz_t YA)
{
    mpz_t cmp_A, neg_cmp_A, temp;
    
    mpz_init(cmp_A);
    mpz_init(neg_cmp_A);
    mpz_init(temp);

    MOTN(cmp_A);
    
    neg_num_A(neg_cmp_A, cmp_A);

    //compute Epk[neg_Cmp_A ·XA ], Epk[Cmp_A ·Y_A ]
    mpz_mul(temp, neg_cmp_A, XA);
    mpz_mod(temp, temp, M);
    encrypt_write(temp);

    mpz_mul(temp, cmp_A, YA);
    mpz_mod(temp, temp, M);
    encrypt_write(temp);
    
    read_decrypt(min_A);
    
    mpz_clear(cmp_A);
    mpz_clear(neg_cmp_A);
    mpz_clear(temp);
    
}

void kMeansProtocol::MOTN_termination(bool * cmp)
{
    mpz_t P,Q, bit;
    mpz_init(P);
    mpz_init(Q);
    mpz_init(bit);

    read_decrypt(P);
    read_decrypt(Q);
    
    short  int op = mpz_cmp(P,Q); //Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
    if(op >0) // P>Q
        mpz_set_ui(bit, 1);
    else mpz_set_ui(bit, 0);//P<=Q terminate: L<=R i.e. dist <= epsilon

    char * strValue = mpz_get_str(NULL, 10, bit);
    connWithBob->writeDataWithLength(&strValue);
    free(strValue);
        
    *cmp = (bool)mpz_get_ui(bit);
    
    mpz_clear(P);
    mpz_clear(Q);
    mpz_clear(bit);

}



void kMeansProtocol::MOMN_n( mpz_t cmp_A[NUMDATAPOINT])
{
    mpz_t P,Q, bit, W_min_A, W_A, neg_A, temp;
    mpz_init(P);
    mpz_init(Q);
    mpz_init(bit);
    mpz_init(W_min_A);
    mpz_init(W_A);
    mpz_init(neg_A);
    mpz_init(temp);

    mpz_set_ui(cmp_A[0],0);
    
    for(long  int i =1; i<NUMDATAPOINT; ++i)
    {
        //MOTN
        read_decrypt(P);
        read_decrypt(Q);
        
        int cmp = mpz_cmp(P,Q); //Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
        if(cmp >0) // P>Q
            mpz_set_ui(bit, 1);
        else mpz_set_ui(bit, 0); //P<=Q
        
        encrypt_write(bit);

        if(i != NUMDATAPOINT-1)
        {
        read_decrypt(W_min_A);
        read_decrypt(W_A);
        read_decrypt(cmp_A[i]); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector
            //gmp_printf("\ncmp :%Zd", cmp_A[i] );

        neg_num_A(neg_A, cmp_A[i]);

        //compute Epk[neg_Cmp_A ·XA ], Epk[Cmp_A ·Y_A ]
        mpz_mul(temp, neg_A, W_min_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);

        mpz_mul(temp, cmp_A[i], W_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);
            
        }
        else
        {
            read_decrypt(cmp_A[i]);
            //gmp_printf("\ncmp :%Zd", cmp_A[i] );

        }

    }
    
    mpz_set_ui(W_A, 0);
    
    for(long int i=NUMDATAPOINT-2; i>=0; --i)
    {
        //gmp_printf("\ncmp i:%Zd, i+1:%Zd", cmp_A[i], cmp_A[i+1] );
        mpz_add(W_A, W_A, cmp_A[i+1]);
        
        mpz_mul(temp, W_A, cmp_A[i]);
        encrypt_write(temp);
        
        read_decrypt(W_min_A);
        
        neg_num_A(neg_A, W_A);
        
        mpz_mul(temp, neg_A, W_min_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);
        
        read_decrypt(cmp_A[i]);

    }
    
 /*   for(long  int i =0; i<NUMDATAPOINT; ++i)
        read_decrypt(cmp_A[i]);*/
    
    
    mpz_clear(P);
    mpz_clear(Q);
    mpz_clear(bit);
    mpz_clear(W_min_A);
    mpz_clear(W_A);
    mpz_clear(neg_A);
    mpz_clear(temp);

}
    
    
void kMeansProtocol::MOMN_k(long int dp_num)
{
    mpz_t P,Q, bit, W_min_A, W_A, neg_A, temp;
    mpz_init(P);
    mpz_init(Q);
    mpz_init(bit);
    mpz_init(W_min_A);
    mpz_init(W_A);
    mpz_init(neg_A);
    mpz_init(temp);

    mpz_set_ui(closest_center_A[dp_num][0], 0);

    for(long  int i =1; i<NUMCENTERS; ++i)
    {
        //MOTN
        read_decrypt(P);
        read_decrypt(Q);
        
        int cmp = mpz_cmp(P,Q); //Return a positive value if op1 > op2, zero if op1 = op2, or a negative value if op1 < op2.
        if(cmp >0) // P>Q
            mpz_set_ui(bit, 1);
        else mpz_set_ui(bit, 0); //P<=Q
        
        encrypt_write(bit);

        if(i != NUMCENTERS-1)
        {
        read_decrypt(W_min_A);
        read_decrypt(W_A);
        read_decrypt(closest_center_A[dp_num][i]); //this is not for i th location...but to save memory we r storing it in this..this is result of permuted vector
        
        neg_num_A(neg_A, closest_center_A[dp_num][i]);

        //compute Epk[neg_Cmp_A ·XA ], Epk[Cmp_A ·Y_A ]
        mpz_mul(temp, neg_A, W_min_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);

        mpz_mul(temp, closest_center_A[dp_num][i], W_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);
            
        }
        else
        {
            read_decrypt(closest_center_A[dp_num][i]);
        }
        
    }
    
    mpz_set_ui(W_A, 0);
    
    for(long int i=NUMCENTERS-2; i>=0; --i)
    {
        mpz_add(W_A, W_A, closest_center_A[dp_num][i+1]);
        
        mpz_mul(temp, W_A, closest_center_A[dp_num][i]);
        encrypt_write(temp);
        
        read_decrypt(W_min_A);
        
        neg_num_A(neg_A, W_A);
        
        mpz_mul(temp, neg_A, W_min_A);
        mpz_mod(temp, temp, M);
        encrypt_write(temp);
        
        read_decrypt(closest_center_A[dp_num][i]);

    }
    
 /*   for(long  int i =0; i<NUMCENTERS; ++i)
        read_decrypt(cmp_A[i]);*/
    
    
    mpz_clear(P);
    mpz_clear(Q);
    mpz_clear(bit);
    mpz_clear(W_min_A);
    mpz_clear(W_A);
    mpz_clear(neg_A);
    mpz_clear(temp);

}


void kMeansProtocol::RVS(mpz_t rand_A, mpz_t QA)
{
    
    mpz_t qA[LValue], S, S_prev, w, aA, aA_prev, cmp_A, vA, eA[LValue];
    long  int  perm[LValue];
    
    mpz_init(S);
    mpz_init(w);
    mpz_init(cmp_A);
    mpz_init(aA);
    mpz_init(vA);
    mpz_init(S_prev);
    mpz_init(aA_prev);
    
    for(long  int i =0; i<L; ++i)
    {
        mpz_init(qA[i]);

        mpz_init(eA[i]);
        
        perm[i]=i;

    }
    
    
    
    
    TBP(qA, QA);
    
    //reordering_L(perm);
    
    mpz_set(eA[perm[0]], qA[perm[0]]);
    mpz_set_ui(cmp_A,0);
    
    for(long  int i =1; i<L; ++i)
    {
        
        mpz_add(cmp_A, cmp_A, eA[perm[i-1]]);
        mpz_mul(w, cmp_A, qA[perm[i]]);
        
        encrypt_write(w);
        
        read_decrypt(eA[i]);
        
    }

    mpz_set_ui(S, 1);
    mpz_set_ui(aA_prev, 0);

    for(long  int i =0; i<L; ++i)
    {
        if(i>0)
        {
            mpz_mul_ui(S, S_prev, 2);
            mpz_mul(aA, aA_prev, qA[i-1]);

            if(i>1)
                mpz_mul_ui(aA, aA, 2);
            
        }
        mpz_sub_ui(S, S, 1);
        
        if(i==0)
            mpz_set_ui(w, 0);
        else
            getRandomN(w, i, S);
 
        encrypt_write(w);
        
        MOTN(cmp_A);
        
        read_decrypt(vA);
    
        mpz_set(aA_prev, aA);
        mpz_set(S_prev, S);
        
        //comp final
        mpz_add(w, vA, aA);//vA+aA
        encrypt_write(w);
        
#ifdef sync
        //sleep(SLEEPTIME);
        read_decrypt(w);
#endif
        
        mpz_mul(w, w, eA[i]);//eA(vA+aA)
        encrypt_write(w);
        
    }

    read_decrypt(rand_A);

    mpz_clear(w);
    mpz_clear(cmp_A);
    mpz_clear(S);
    mpz_clear(aA);
    mpz_clear(vA);
    mpz_clear(S_prev);
    mpz_clear(aA_prev);
    
    for(long  int i =0; i<L; ++i)
    {
        mpz_clear(qA[i]);

        mpz_clear(eA[i]);
    }
}

//other sub-protocol -write
void kMeansProtocol::comp_distance( mpz_t dist_A, long  int dp_num, long  int cent_num)
{
    mpz_t y1_A, sq_sum_dp_A, sq_sum_cent_A;
    mpz_init(y1_A);
    mpz_init(sq_sum_dp_A);
    mpz_init(sq_sum_cent_A);

    mpz_set_ui(sq_sum_dp_A, 0);//d^2
    mpz_set_ui(sq_sum_cent_A, 0);//cA^2
    
    for(long  int j=0; j<NUMATTRIBUTE; ++j)
    {
        mpz_addmul(sq_sum_dp_A, data_pt_A[dp_num][j], data_pt_A[dp_num][j]);
        mpz_mod(sq_sum_dp_A, sq_sum_dp_A, M);//dA^2

        mpz_addmul(sq_sum_cent_A, center_A[cent_num][j], center_A[cent_num][j]);
        mpz_mod(sq_sum_cent_A, sq_sum_cent_A, M);//cA^2


        mpz_addmul(y1_A, center_A[cent_num][j], data_pt_A[dp_num][j]);
        mpz_mod(y1_A, y1_A, M);//CA.DA
    }
    
    mpz_add(dist_A, sq_sum_dp_A, sq_sum_cent_A);//dA^2+CA^2
    
    mpz_mul_ui(y1_A, y1_A,2);
    mpz_mod(y1_A, y1_A, M);//2CA.DA
    mpz_invert(y1_A, y1_A, M);//-2CA.DA
    
    mpz_add(dist_A, dist_A, y1_A);//dA^2+CA^2-2CA.DA
    mpz_mod(dist_A, dist_A, M);

    encrypt_write(dist_A);
    
    read_decrypt(dist_A);
    
    mpz_clear(y1_A);
    mpz_clear(sq_sum_dp_A);
    mpz_clear(sq_sum_cent_A);
    
}


void kMeansProtocol::comp_distance_short( mpz_t dist_A, long  int dp_num, long  int cent_num)
{
    mpz_t x_A, xAsqr;
    mpz_init(x_A);
    mpz_init(xAsqr);
    
    for(long  int j=0; j<NUMATTRIBUTE; ++j)
    {
        read_decrypt(x_A);
        mpz_mul(xAsqr, x_A, x_A);
        mpz_mod(xAsqr, xAsqr, M);
        
        encrypt_write(xAsqr);
    }
    
    
    read_decrypt(dist_A);
    
    mpz_clear(x_A);
    mpz_clear(xAsqr);
    
}




void kMeansProtocol::TBP(mpz_t qA[LValue], mpz_t QA)
{
    mpz_t S, XOR1_A, carry1_A, carry2_A, carry3_A, carry_prev_A, Ex2carry1_A, Ex2carry2_A, Ex2carry3_A, Ex2carry_prev_A, mul, XOR2_A, cmp_A, neg_cmp_A;
    mpz_init(S);
    mpz_init(carry1_A);
    mpz_init(carry2_A);
    mpz_init(carry3_A);
    mpz_init(carry_prev_A);
    mpz_init(Ex2carry1_A);
    mpz_init(Ex2carry2_A);
    mpz_init(Ex2carry3_A);
    mpz_init(Ex2carry_prev_A);
    mpz_init(mul);
    mpz_init(cmp_A);
    mpz_init(neg_cmp_A);
    mpz_init(XOR1_A);
    mpz_init(XOR2_A);
    
    //binary of QB
    mpz_set(S, QA);
    for(long  int j=0; j<L; j++)
    {
        //mpz_init(QA_bin[j]);
        //mpz_init(XOR1_A[j]);
        //mpz_init(XOR2_A[j]);

        mpz_mod_ui(qA[j], S, 2);
        mpz_div_ui(S, S, 2);
    }
    
    mpz_set_ui(carry_prev_A, 0);
    mpz_set_ui(Ex2carry_prev_A, 0);

    //run FM2NP
    MOTN(cmp_A);
    neg_num_A(neg_cmp_A, cmp_A);
    //cmp=0 if XA+XB<M
    
    for(long  int j=0; j<L; j++)
    {
        encrypt_write(qA[j]);
        read_decrypt(XOR1_A);
        
        mpz_mul(mul, XOR1_A, carry_prev_A);//XOR1_A[j]*E_carry_prev_A
        encrypt_write(mul);
        
        //xor with prev carry = a+b
        read_decrypt(XOR1_A);
        
        //curr carry
        read_decrypt(carry1_A);
        read_decrypt(carry2_A);
        
        mpz_mul(mul, carry_prev_A, carry2_A);
        encrypt_write(mul);
        
        read_decrypt(carry3_A);
        mpz_mul(mul, carry1_A, carry3_A);
        encrypt_write(mul);
        
        read_decrypt(carry_prev_A);
        

    //Now compute a + b + d =xor+d
        //a ex b ex d
        read_decrypt(XOR2_A);
        mpz_mul(mul, XOR2_A, Ex2carry_prev_A);
        encrypt_write(mul);
        read_decrypt(XOR2_A);
                        
        //compute carry to add in next iteration
        // compute (a ex b)*d
        //curr carry =a.b or b.c or c.a = a.b or c(a or b)
        read_decrypt(Ex2carry1_A);
        read_decrypt(Ex2carry2_A);
        
        mpz_mul(mul, Ex2carry_prev_A, Ex2carry2_A);
        encrypt_write(mul);
        read_decrypt(Ex2carry3_A);
        mpz_mul(mul, Ex2carry1_A, Ex2carry3_A);
        encrypt_write(mul);
        read_decrypt(Ex2carry_prev_A);
    
    
    //now we have a+b and a+b+d

        mpz_mul(mul, neg_cmp_A, XOR1_A);
        encrypt_write(mul);
            
        mpz_mul(mul, cmp_A, XOR2_A);
        encrypt_write(mul);

        read_decrypt(qA[j]);
    }
    
    mpz_clear(S);
    mpz_clear(carry1_A);
    mpz_clear(carry2_A);
    mpz_clear(carry3_A);
    mpz_clear(carry_prev_A);
    mpz_clear(Ex2carry1_A);
    mpz_clear(Ex2carry2_A);
    mpz_clear(Ex2carry3_A);
    mpz_clear(Ex2carry_prev_A);
    mpz_clear(mul);
    mpz_clear(cmp_A);
    mpz_clear(neg_cmp_A);
    mpz_clear(XOR1_A);
    mpz_clear(XOR2_A);
    
}

void kMeansProtocol::reordering_L( long  int  perm[LValue])
{
    
}
void kMeansProtocol::choose_center(mpz_t randA, long  int center_num)
{
    mpz_t ZA[NUMDATAPOINT], cmp_A[NUMDATAPOINT], temp;
    
    mpz_init(temp);
    for(long  int i =0; i<NUMDATAPOINT; ++i)
    {
        mpz_init(ZA[i]);
        mpz_init(cmp_A[i]);
    }
    
    if(center_num==0)
    {
        mpz_mul_ui(temp, min_sc_dist_A[0], NUMDATAPOINT);
        mpz_add(ZA[0], avg_min_sc_dist_A, temp);
        mpz_mod(ZA[0], ZA[0], M);
        mpz_sub(ZA[0], ZA[0], randA);
        mpz_mod(ZA[0], ZA[0], M);
        encrypt_write(ZA[0]);

        for(long  int i =1; i<NUMDATAPOINT; ++i)
        {
            //ZA=avg_min_sc_dist_A + nsc_dist_A[i]
            mpz_mul_ui(temp, min_sc_dist_A[0], NUMDATAPOINT);
            mpz_add(temp, temp, avg_min_sc_dist_A);
            mpz_mod(ZA[i], ZA[i], M);
            mpz_add(ZA[i], ZA[i-1], temp);
            mpz_sub(ZA[i], ZA[i], randA);
            mpz_mod(ZA[i], ZA[i], M);
            
            encrypt_write(ZA[i]);
        }
    
    }
    else
    {
        mpz_sub(ZA[0], min_sc_dist_A[0], randA);
        mpz_mod(ZA[0], ZA[0], M);
        encrypt_write(ZA[0]);

        for(long  int i =1; i<NUMDATAPOINT; ++i)
        {
            //ZA= sc_dist_A[i]
            mpz_add(ZA[i], ZA[i-1], min_sc_dist_A[i]);
            mpz_mod(ZA[i], ZA[i], M);
            
            encrypt_write(ZA[i]);
        }
    
    }

    MOMN_n(cmp_A);
    
    for(long  int j =0; j<NUMATTRIBUTE; ++j)
    {
        mpz_set_ui(center_A[center_num][j],  0);
        //cout<<"\nAttribute "<<j<<endl;

        for(long  int i =0; i<NUMDATAPOINT; ++i)
        {
            /*center_A[j] = center_A[j] + cmp_A[i] * data_pt_A[i][j];
*/
            
            mpz_mul(temp, cmp_A[i], data_pt_A[i][j]);
            mpz_add(center_A[center_num][j], center_A[center_num][j], temp);
            mpz_mod(center_A[center_num][j], center_A[center_num][j], M);
            
            //gmp_printf("%Zd ",center_A[center_num][j]);

        }
        
        encrypt_write(center_A[center_num][j]);

        read_decrypt(center_A[center_num][j]);
        //gmp_printf(" %Zd ",center_A[center_num][j]);

    }
    
    mpz_clear(temp);
    for(long  int i =0; i<NUMDATAPOINT; ++i)
    {
        mpz_clear(ZA[i]);
        mpz_clear(cmp_A[i]);
    }
    
    
}

//Scaling -write
void kMeansProtocol::scaling_fact()
{
    mpz_t temp, num_prod_A[NUMCENTERS];
    mpz_init(temp);
            
    for(long  int j =0; j<NUMCENTERS; ++j)
    {
        mpz_init(num_prod_A[j]);

        mpz_set_ui(alpha_A[j], 1);
        mpz_set_ui(alpha_A[j], 1);

        mpz_mul(temp, num_dp_A[j], num_dp_new_A[j]);
        encrypt_write(temp);
        
    }
    
    for(long  int j =0; j<NUMCENTERS; ++j)
    read_decrypt(num_prod_A[j]);
    
    for(long  int j =0; j<NUMCENTERS; ++j)
    {
        for(long  int i =0; i<NUMCENTERS && i!=j ; ++i)
        {
            mpz_mul(temp, alpha_A[j], num_dp_new_A[i]);
            encrypt_write(temp);
            mpz_mul(temp, beta_A[j], num_prod_A[i]);
            encrypt_write(temp);
            
            read_decrypt(alpha_A[j]);
            read_decrypt(beta_A[j]);
        }
    }

    mpz_mul(temp, alpha_A[1], num_dp_new_A[1]);
    encrypt_write(temp);
    mpz_mul(temp, beta_A[1], num_prod_A[1]);
    encrypt_write(temp);
    
    read_decrypt(G_alpha_A);
    read_decrypt(G_beta_A);
    
    mpz_clear(temp);
            
    for(long  int j =0; j<NUMCENTERS; ++j)
    {
        mpz_clear(num_prod_A[j]);
    }
    
    
}


void kMeansProtocol::comp_scaled_dist( long  int dp_num, long  int cent_num )
{
    mpz_t mul, w;
    mpz_init(mul);
    mpz_init(w);
    
    for(long  int i =0; i<NUMATTRIBUTE; ++i)
    {
        mpz_mul(mul, alpha_A[cent_num], center_A[cent_num][i]);
        encrypt_write(mul);
        
        mpz_mul(mul, G_alpha_A, data_pt_A[dp_num][i]);
        encrypt_write(mul);
        
        read_decrypt(w);
        
        mpz_mul(mul,w, w);
        encrypt_write(mul);
    }
    
    read_decrypt(sc_dist_A[dp_num][cent_num]);

    mpz_clear(mul);
    mpz_clear(w);
                
}


void kMeansProtocol::termination_condition(bool *b)//b=0 terminate: L<=R i.e. dist <= epsilon
{
    mpz_t mul, w, v_A, betasqr_A, L_A, R_A;
    mpz_init(mul);
    mpz_init(w);
    mpz_init(v_A);
    mpz_init(betasqr_A);
    mpz_init(L_A);
    mpz_init(R_A);

        for(long  int j =0; j<NUMCENTERS; ++j)
    {
        
        for(long  int i =0; i<NUMATTRIBUTE; ++i)
        {
            mpz_mul(mul, num_dp_new_A[j], center_A[j][i]);
            encrypt_write(mul);
            
            mpz_mul(mul, num_dp_A[j], center_new_A[j][i]);
            encrypt_write(mul);
            
            read_decrypt(w);
            
            mpz_mul(mul,w, w);
            encrypt_write(mul);
        }
        
        read_decrypt(v_A);
        
        mpz_mul(mul, beta_A[j], beta_A[j]);
        encrypt_write(mul);
        
        read_decrypt(betasqr_A);
        
        mpz_mul(mul, betasqr_A, v_A);
        encrypt_write(mul);
        
    }
    
    read_decrypt(L_A);
    
    mpz_mul(mul, G_beta_A, G_beta_A);
    encrypt_write(mul);
    
    read_decrypt(R_A);
    
    MOTN_termination(b);
    
    if(*b==1)
    {
        
        for(long   int j =0; j<NUMCENTERS; ++j)
        {
            for(long   int i =0; i<NUMATTRIBUTE; ++i)
            {
                mpz_set(center_A[j][i], center_new_A[j][i]);
                
                mpz_set(sum_dp_A[j][i], sum_dp_new_A[j][i]);
            }
            mpz_set(num_dp_A[j], num_dp_new_A[j]);
            
        }
    }


    mpz_clear(mul);
    mpz_clear(w);
    mpz_clear(v_A);
    mpz_clear(betasqr_A);
    mpz_clear(L_A);
    mpz_clear(R_A);
    
}
