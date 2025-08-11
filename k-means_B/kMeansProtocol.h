//
//  kMeans.hpp
//  k-means
//
//  Created by sonali chim on 22/08/18.
//  Copyright © 2018 personal. All rights reserved.
//

#include <thread>
#include <cmath>
#include <vector>
#include<array>

#include "PaillierSystem.h"
#include "server.h"

#define LValue 14//2048//3072//32//1024//512//
#define KValue 6//1000//1500//10//500//250//
#define NUMDATAPOINT 10//569
#define NUMMYDATAPOINT 5//B:300 A:269
#define NUMATTRIBUTE 2//30
#define NUMCENTERS 2
#define EPSILON 14807104//14807104 //sqr root 3848
#define STARTSKIP 2
#define TAILSKIP 28
#define SCALING 1
#define SLEEPTIME 0

//#define dist_test
//#define test
//#define ccp_test
//#define momn_test
//#define momnk_test
//#define pickRandom
//#define xor_test

extern connection * connWithAlice;


class kMeansProtocol
{
    
    //PaillierSystem *myPaillierObj;
   // PaillierSystem *othPaillierObj;
    
    long   int L; //K=security parameter
    mpz_t M; //modulus N chosen by Alice = pub key of A in case of pailliers
    long   int K; //L > 2K+1
    mpz_t N; //range of size K
    mpz_t Msqr;
    
    
    mpz_t data_pt_B[NUMDATAPOINT][NUMATTRIBUTE]; //two dim array of int
    mpz_t E_data_pt_A[NUMDATAPOINT][NUMATTRIBUTE]; //two dim array of int
    mpz_t E_data_pt[NUMDATAPOINT][NUMATTRIBUTE]; //two dim array of int

    //centers
    mpz_t center_B[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_center_A[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_center[NUMCENTERS][NUMATTRIBUTE];
    mpz_t fact_cent[NUMCENTERS];

    mpz_t center_new_B[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_center_new_A[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_center_new[NUMCENTERS][NUMATTRIBUTE];
    mpz_t fact_cent_new[NUMCENTERS];

    mpz_t cog_B[NUMATTRIBUTE];
    mpz_t E_cog_A[NUMATTRIBUTE];
    mpz_t E_cog[NUMATTRIBUTE];

    //distances
    /*mpz_t dist_B[NUMDATAPOINT][NUMCENTERS]; //k+1 for cog
    mpz_t E_dist_A[NUMDATAPOINT][NUMCENTERS];
    mpz_t E_dist[NUMDATAPOINT][NUMCENTERS];*/

    mpz_t COGdist_B[NUMDATAPOINT];
    mpz_t E_COGdist_A[NUMDATAPOINT];
    mpz_t E_COGdist[NUMDATAPOINT];

    /*mpz_t min_dist_B[NUMDATAPOINT];
    mpz_t E_min_dist_A[NUMDATAPOINT];
    mpz_t E_min_dist[NUMDATAPOINT];*/

    mpz_t avg_min_dist_B;
    mpz_t E_avg_min_dist_A;
    mpz_t E_avg_min_dist;

    mpz_t sc_dist_B[NUMDATAPOINT][NUMCENTERS];
    mpz_t E_sc_dist_A[NUMDATAPOINT][NUMCENTERS];
    mpz_t E_sc_dist[NUMDATAPOINT][NUMCENTERS];

    mpz_t min_sc_dist_B[NUMDATAPOINT];
    mpz_t E_min_sc_dist_A[NUMDATAPOINT];
    mpz_t E_min_sc_dist[NUMDATAPOINT];

    mpz_t avg_min_sc_dist_B;
    mpz_t E_avg_min_sc_dist_A;
    mpz_t E_avg_min_sc_dist;

    //closest center
    mpz_t closest_center_B[NUMDATAPOINT][NUMCENTERS];//closest center details at the current iteration for each datapoint
    mpz_t E_closest_center_A[NUMDATAPOINT][NUMCENTERS];
    mpz_t E_closest_center[NUMDATAPOINT][NUMCENTERS];


    mpz_t num_dp_new_B[NUMCENTERS];
    mpz_t E_num_dp_new_A[NUMCENTERS];
    mpz_t E_num_dp_new[NUMCENTERS];
        
    mpz_t num_dp_B[NUMCENTERS];
    mpz_t E_num_dp_A[NUMCENTERS];
    mpz_t E_num_dp[NUMCENTERS];
 
    mpz_t sum_dp_new_B[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_sum_dp_new_A[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_sum_dp_new[NUMCENTERS][NUMATTRIBUTE];

    mpz_t sum_dp_B[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_sum_dp_A[NUMCENTERS][NUMATTRIBUTE];
    mpz_t E_sum_dp[NUMCENTERS][NUMATTRIBUTE];
    
    //scaling
    mpz_t G_alpha_B;
    mpz_t E_G_alpha_A;
    mpz_t E_G_alpha;

    mpz_t G_beta_B;
    mpz_t E_G_beta_A;
    mpz_t E_G_beta;

    mpz_t alpha_B[NUMCENTERS];
    mpz_t E_alpha_A[NUMCENTERS];
    mpz_t E_alpha[NUMCENTERS];

    mpz_t beta_B[NUMCENTERS];
    mpz_t E_beta_A[NUMCENTERS];
    mpz_t E_beta[NUMCENTERS];

    PaillierSystem *othPaillierObj;

    public:

    //SubProtocols
    
    //Supporting -NO write
    void prod_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB);
    void scalar_prod(mpz_t E_op, mpz_t E_X, mpz_t Y);
    void neg_shared_num(mpz_t E_op, mpz_t E_op_A, mpz_t op_B, mpz_t E_XA, mpz_t XB);
    void XOR_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB);
    void XOR(mpz_t E_op, mpz_t E_X, mpz_t Y);
    void OR_shared_num(mpz_t E_op, mpz_t E_XA, mpz_t E_YA, mpz_t E_XAYA, mpz_t XB, mpz_t YB);
    void OR(mpz_t E_op, mpz_t E_X, mpz_t Y);

    void share_write(mpz_t E_op_A, mpz_t op_B, mpz_t E_op);
    void share(mpz_t E_op_A, mpz_t op_B, mpz_t E_op);
    void write(mpz_t E_value_A);
    void read (mpz_t E_value_A);

    //sub-proto -write
    void MOTN(mpz_t E_cmp, mpz_t E_cmp_A, mpz_t cmp_B, mpz_t E_XA, mpz_t YA, mpz_t XB, mpz_t YB);
    void MOTN_value(mpz_t E_min, mpz_t E_min_A, mpz_t min_B, mpz_t E_XA, mpz_t YA, mpz_t XB, mpz_t YB);
    void MOTN_termination(bool * cmp, mpz_t E_XA, mpz_t YA, mpz_t XB, mpz_t YB);
    void MOMN_n(mpz_t E_cmp_arr[NUMDATAPOINT], mpz_t E_cmp_arr_A[NUMDATAPOINT], mpz_t cmp_arr_B[NUMDATAPOINT], mpz_t E_X[NUMDATAPOINT], mpz_t E_XA[NUMDATAPOINT], mpz_t XB[NUMDATAPOINT]);
    void MOMN_k(mpz_t E_X[NUMDATAPOINT], mpz_t E_XA[NUMCENTERS], mpz_t XB[NUMCENTERS], long int dp_num);
    void RVS(mpz_t E_rand, mpz_t E_rand_A, mpz_t rand_B, mpz_t E_QA, mpz_t QB);
    
    //other sub-protocol -write
    //void comp_distance(mpz_t E_dist, mpz_t E_dist_A, mpz_t dist_B, long   int dp_num, long   int cent_num);
    void comp_distance_short(mpz_t E_sc_dist, mpz_t E_sc_dist_A, mpz_t sc_dist_B, long   int dp_num, long   int cent_num);
    void TBP(mpz_t E_q[LValue], mpz_t E_qA[LValue], mpz_t qB[LValue], mpz_t E_QA, mpz_t QB);
    void reordering_L(long    int perm[LValue]);
    void choose_center(mpz_t E_randA, mpz_t randB, long   int center_num);
    
    //Scaling -write
    void scaling_fact();
            /*(mpz_t E_GAlpha, mpz_t E_GAlpha_A, mpz_t GAlpha_B,
                      mpz_t E_Alpha[NUMCENTERS], mpz_t E_Alpha_A[NUMCENTERS], mpz_t Alpha_B[NUMCENTERS],
                      mpz_t E_GBeta, mpz_t E_GBeta_A, mpz_t GBlpha_B,
                      mpz_t E_Beta[NUMCENTERS], mpz_t E_Beta_A[NUMCENTERS], mpz_t Blpha_B[NUMCENTERS]);*/
    void comp_scaled_dist(long   int dp_num, long   int cent_num );
    void termination_condition(bool *b); //b=0 terminate: L<=R i.e. dist <= epsilon
    
    //Steps for the protocol by Bob
    //init
    void CalcCOG();
    void CalcDistToCOG_AvgSqrdDist();
    //void CalcAvgSqrdDistanceToCOG( );
    void CalcAvgSqrdDistance();
    void CalcDistToCent(long   int centerNum);
    void PickCenter(long   int centerNum);

    //Lloyd's
    void findClosestCenter();
    void calcNewCenter();
    
   
    //public:

    
    //init
    void getModulus(mpz_t);
    void readOwnData();
    void initDataPoints();
    void readEncDataPts();

    //pick initial cluster centers
    void InitAllClusterCenters();
    
    void runLloydsStep(bool * b);
 
    kMeansProtocol();
    kMeansProtocol( long   int x, mpz_t modulus);
    kMeansProtocol( long   int x, mpz_t modulus, mpz_t othSK);//using this for test purpose
    kMeansProtocol( long   int x, PaillierSystem * pObj);
    ~kMeansProtocol();
};
