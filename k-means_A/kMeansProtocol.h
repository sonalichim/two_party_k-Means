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
#include "client.h"

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


//#define pickRandom

extern connection * connWithBob;

class kMeansProtocol
{
    long  int L, K; //L=security parameter
    mpz_t M; //modulus M chosen by Alice = pub key of M in case of pailliers
    mpz_t Msqr;
    
    mpz_t data_pt_A[NUMDATAPOINT][NUMATTRIBUTE]; //data point belonging to A

    //centers
    mpz_t center_A[NUMCENTERS][NUMATTRIBUTE]; //center points share last iteration

    mpz_t center_new_A[NUMCENTERS][NUMATTRIBUTE]; //center points share scurr iteration

    mpz_t cog_A[NUMATTRIBUTE]; //center of gravity
    mpz_t fact_cog;

    //distances
    mpz_t dist_A[NUMDATAPOINT][NUMCENTERS]; //dist of a dp from a center
    mpz_t COGdist_A[NUMDATAPOINT]; //dist of a dp frm COG

    mpz_t min_dist_A[NUMDATAPOINT]; //min of the dist of a dp from all centers+COg

    mpz_t avg_min_dist_A; //avg of the abv min dist of all dp

    mpz_t sc_dist_A[NUMDATAPOINT][NUMCENTERS]; //scaled dist of a dp from a center

    mpz_t min_sc_dist_A[NUMDATAPOINT]; //min of the scaled dist of a dp from all centers+COg

    mpz_t avg_min_sc_dist_A; //avg of the abv min scaled dist of all dp

    //closest center
    mpz_t closest_center_A[NUMDATAPOINT][NUMCENTERS];//closest center details at the current iteration for each datapoint
    
    mpz_t num_dp_A[NUMCENTERS];//no of dp in a center in prev iter
    mpz_t num_dp_new_A[NUMCENTERS]; //no of dp in a center in curr iter

    mpz_t sum_dp_new_A[NUMCENTERS][NUMATTRIBUTE]; //sum of dp in a center in prev iter
    mpz_t sum_dp_A[NUMCENTERS][NUMATTRIBUTE]; //sum of dp in a center in prev iter

    //scaling
    mpz_t G_alpha_A; //

    mpz_t G_beta_A; //

    mpz_t alpha_A[NUMCENTERS]; //

    mpz_t beta_A[NUMCENTERS]; //
   
public:
    //SubProtocols
    
    //Supporting -NO write

    void encrypt_write(mpz_t value_A);
    void read_decrypt(mpz_t E_value_A);
    void neg_num_A(mpz_t E_op, mpz_t XA);

    //sub-proto -write
    void MOTN(mpz_t cmp_A);
    void MOTN_value(mpz_t min_A, mpz_t XA, mpz_t YA);
    void MOTN_termination(bool * cmp);
    void MOMN_n(mpz_t cmp_A[NUMDATAPOINT]);
    void MOMN_k(long int dp_num);
    void RVS(mpz_t rand_A, mpz_t QA);
    
    //other sub-protocol -write
    void comp_distance(mpz_t dist_A, long  int dp_num, long  int cent_num);
    void comp_distance_short( mpz_t dist_A, long  int dp_num, long  int cent_num);
    void TBP(mpz_t qA[LValue], mpz_t QA);
    void reordering_L( long  int perm[LValue]);
    void choose_center(mpz_t randA, long  int center_num);
    
    //Scaling -write
    void scaling_fact();
    void comp_scaled_dist(long  int dp_num, long  int cent_num );
    void termination_condition(bool *b);//b=0 terminate: L<=R i.e. dist <= epsilon
    
    //Steps for the protocol by Bob
    //init
    void CalcCOG();
    void CalcDistToCOG_AvgSqrdDist();
    //void CalcAvgSqrdDistanceToCOG( );
    void CalcAvgSqrdDistance();
    void CalcDistToCent(long  int centerNum);
    void PickCenter(long  int centerNum);

    //Lloyd's
    void findClosestCenter();
    void calcNewCenter();

    
    //change it to protected after done testing subprotocols
    PaillierSystem *myPaillierObj, *othPaillierObj;
    
    
   //   public:
    //Init
    void initDataPoints();
    void readOwnData();
    void shareDataPoints();
    void getMyModulus(mpz_t);
    void getPrivKey(mpz_t);
    void setOthModulus(mpz_t, long  int x);
    
    void InitAllClusterCenters();
    void runLloydsStep(bool * b);
    
    kMeansProtocol();
    kMeansProtocol( long  int x);
    kMeansProtocol( long  int x, PaillierSystem * pObj);
    ~kMeansProtocol();
};
