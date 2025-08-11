//
//  kMeans.cpp
//  k-means
//
//  Created by sonali chim on 22/08/18.
//  Copyright © 2018 personal. All rights reserved.
//

#include "kMeansProtocol.h"
using namespace std;
using namespace std::chrono;

void kMeansProtocol::readOwnData()
{
    cout<<"\nReading own data"<<endl;

    string word;
    double var;
    
    fstream file;
    file.open("/Users/sonalichim/Dropbox/Code/Sonali/k-Means/k-means_B/k-means_Bob/wdbc.data");
    //file.open ("./k-means_Bob/wdbc.data");
    //file.open ("/home1/chimeknath/new/Bob/wdbc.data");
    
    if(! file.is_open())
        cout<<"%%%%%FILE OPEN ERROR";
    
    for(long   int i=NUMDATAPOINT-NUMMYDATAPOINT; i < NUMDATAPOINT; ++i)
    {
        //cout<<endl<<i<<")";
        for(long   int u=0; u < STARTSKIP; ++u)
        file >> word;
        
        for(long   int j=0; j < NUMATTRIBUTE; ++j)
        {
            file >> word;
            var = atof(word.c_str());
            var=var*SCALING;//add scaling
            mpz_set_ui(data_pt_B[i][j], var);
            
            othPaillierObj->encrypt(E_data_pt[i][j], data_pt_B[i][j]);
            //gmp_printf(" %Zd ",data_pt_B[i][j]);

        }
        for(long   int u=NUMATTRIBUTE; u < NUMATTRIBUTE+TAILSKIP; ++u)
        file >> word;
    }
    
    //closing the file
    file.close();
}



void kMeansProtocol::initDataPoints()
{
    cout<<"\nInit dp"<<endl;

    mpz_t temp, enc;
    mpz_init_set_ui(temp,0);
    mpz_init(enc);
    othPaillierObj->encrypt(enc, temp);
    
    for(long   int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long   int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_init_set_ui(data_pt_B[i][j],0);
            mpz_init_set(E_data_pt[i][j], enc);
            mpz_init_set(E_data_pt_A[i][j], enc);
            
        }
    }
    
    readOwnData();
}

void kMeansProtocol::readEncDataPts()
{
    cout<<"\n Read shared data"<<endl;

    mpz_t ciph;
    mpz_init(ciph);
    char * strValue;
    
    long int count =NUMDATAPOINT*NUMATTRIBUTE;
    
    //connWithAlice->readnData(count);
    
    for(long   int i=0; i<NUMDATAPOINT; ++i)
    {
        cout<<endl<<i<<")";
        for(long   int j=0; j<NUMATTRIBUTE; ++j, count++)
        {
            /*mpz_init(E_data_pt_A[i][j]);
            connWithAlice->readDataWithLength(&strValue);
            mpz_set_str(E_data_pt_A[i][j], strValue, 10);
            free(strValue);
            */
            //connWithAlice->readDataWithLength(&strValue);
            //free(strValue);
            
            read(E_data_pt_A[i][j]);
            mpz_mul(E_data_pt[i][j], E_data_pt[i][j], E_data_pt_A[i][j]);
            mpz_mod(E_data_pt[i][j], E_data_pt[i][j], Msqr);
            
            othPaillierObj->decrypt(ciph, E_data_pt[i][j]);
           gmp_printf(" %Zd ",ciph);
        }
    }
    mpz_clear(ciph);
}

void kMeansProtocol::getModulus(mpz_t modulus)
{
    mpz_set(modulus, M);
}



kMeansProtocol::kMeansProtocol()
{
    
}


//use this constructor
kMeansProtocol::kMeansProtocol(long   int x, mpz_t modulus, mpz_t othSK)
{
    L = x;
    othPaillierObj = new PaillierSystem(x/2, modulus, othSK);
    mpz_init(M);
    mpz_init(Msqr);
    mpz_set(M, modulus);
    mpz_mul(Msqr, M, M);
    mpz_init(N);
    mpz_ui_pow_ui(N,2, KValue);
}





kMeansProtocol::~kMeansProtocol()
{
    delete connWithAlice;
    //delete myPaillierObj;
    delete othPaillierObj;
    mpz_clear(M);
    mpz_clear(Msqr);
    
    for(long   int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long   int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_clear(data_pt_B[i][j]);
            mpz_clear(E_data_pt[i][j]);
        }
    }
    
    for(long   int i=0; i < NUMCENTERS; ++i)
    {
        for(long   int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_clear(center_new_B[i][j]);
            mpz_clear(E_center_new[i][j]);
        }
    }
}

void kMeansProtocol::CalcCOG()
{
    mpz_t temp;
    mpz_init(temp);
    
    cout<<"\nCalculating COG"<<endl;
    
    for(long   int j=0; j <  NUMATTRIBUTE; ++j)
    {
        mpz_set_ui(cog_B[j],0);
        
        for(long   int i=0; i < NUMDATAPOINT; ++i)
        {
            //othPaillierObj->decrypt(temp2, E_data_pt[i][j]);
            if(mpz_cmp_ui(data_pt_B[i][j],0) !=0)
                mpz_add(cog_B[j], cog_B[j], data_pt_B[i][j]);
            
            //othPaillierObj->decrypt(temp2, sum);
        }
        
        //gmp_printf("\nSum CoG_B(%d)=%Zd", j,w);
        
        //gmp_printf("  B Sum=%Zd", cog_B[j]);
        
        read(E_cog_A[j]);
        
        othPaillierObj->encrypt(E_cog[j], cog_B[j]);
        mpz_mul(E_cog[j],E_cog[j],E_cog_A[j]);
        mpz_mod(E_cog[j], E_cog[j], Msqr);
        
        share_write(E_cog_A[j], cog_B[j], E_cog[j]);
        
        //othPaillierObj->decrypt(temp, E_cog[j]);
        //gmp_printf(" %Zd ", temp);

    }
    
    //mpz_set_ui(centerOfGravity[0],51); for mod 143
    //mpz_set_ui(centerOfGravity[1],81);
    
    mpz_clear(temp);
}



void kMeansProtocol::CalcDistToCOG_AvgSqrdDist()
{
    cout<<"\nCalculating dist to COG";

    mpz_t E_mul1, E_w, E_w_A, w_B, x, y, E_wsqr, temp;
    mpz_init_set_ui(E_mul1,0);
    mpz_init(E_w);
    mpz_init(E_w_A);
    mpz_init(w_B);
    mpz_init(x);
    mpz_init(y);
    mpz_init(E_wsqr);
    mpz_init(temp);
    
    mpz_set_ui(avg_min_sc_dist_B, 0);
    othPaillierObj->encrypt(E_avg_min_sc_dist_A, avg_min_sc_dist_B);
    othPaillierObj->encrypt(E_avg_min_sc_dist, avg_min_sc_dist_B);

    mpz_set_ui(w_B, 0);
    othPaillierObj->encrypt(x, w_B);
    
    for(long   int i =0; i<NUMDATAPOINT; ++i)
    {
        for(long   int j =0; j<NUMATTRIBUTE; ++j)
        {
            mpz_powm_ui(y, E_data_pt[i][j], NUMDATAPOINT, Msqr);//n.d[i,j]
            //othPaillierObj->decrypt(temp, E_data_pt[i][j]);
            //othPaillierObj->decrypt(temp, y);

            mpz_invert(y, y, Msqr);//-n.d[i,j]
            mpz_mul(E_w, E_cog[j], y);//c[0][j]-n.d[i,j]
            mpz_mod(E_w, E_w, Msqr);
            
            //othPaillierObj->decrypt(temp, E_w);
            
            share_write(E_w_A, w_B, E_w);
            
            read(E_mul1);
            
            prod_shared_num(E_wsqr, E_w_A, E_w_A, E_mul1, w_B, w_B);//(c[0][j]-n.d[i,j])^2
            
            mpz_mul(E_COGdist[i], E_COGdist[i], E_wsqr);
            mpz_mod(E_COGdist[i], E_COGdist[i], Msqr);
            

        }

        
        share_write( E_COGdist_A[i],  COGdist_B[i], E_COGdist[i]);
        mpz_set(min_sc_dist_B[i], COGdist_B[i]); //scaled by n^2
        mpz_set(E_min_sc_dist_A[i], E_COGdist_A[i]);
        mpz_set(E_min_sc_dist[i], E_COGdist[i]);
        
        mpz_add(avg_min_sc_dist_B, avg_min_sc_dist_B, min_sc_dist_B[i]);//scaled by n^3 since we r not diving the dist
        mpz_mod(avg_min_sc_dist_B, avg_min_sc_dist_B, M);
        
        mpz_mul(E_avg_min_sc_dist_A, E_avg_min_sc_dist_A, E_min_sc_dist_A[i]);
        mpz_mul(E_avg_min_sc_dist, E_avg_min_sc_dist, E_min_sc_dist[i]);

        mpz_mod(E_avg_min_sc_dist, E_avg_min_sc_dist, Msqr);
        
#ifdef cog_test
        othPaillierObj->decrypt(E_mul1, E_COGdist[i]);
        gmp_printf("\nE_COGdist[i]:%Zd, ",E_mul1);
#endif
    }
    
    /*othPaillierObj->encrypt(E_w, avg_min_sc_dist_B);
    mpz_mul(E_avg_min_sc_dist, E_avg_min_sc_dist_A, E_w);
    mpz_mod(E_avg_min_sc_dist, E_avg_min_sc_dist, Msqr);*/

#ifdef cog_test
    othPaillierObj->decrypt(E_mul1, E_avg_min_sc_dist);
    gmp_printf("\nAvg Dist COG:%Zd, ",E_mul1);
#endif
    
    mpz_clear(E_mul1);
    mpz_clear(E_w);
    mpz_clear(E_w_A);
    mpz_clear(w_B);
    mpz_clear(y);
    mpz_clear(E_wsqr);
    mpz_clear(temp);
}



void kMeansProtocol::PickCenter(long   int centerNum)
{
    cout<<"\n\nPick center "<<centerNum<<endl;

    mpz_t E_rand, E_rand_A, rand_B,  E_QA, QB;
    mpz_init(E_rand);
    mpz_init(E_rand_A);
    mpz_init(rand_B);
    mpz_init(E_QA);
    mpz_init(QB);

    if(centerNum == 0)
    {
        long   int mul = NUMDATAPOINT*2;
        mpz_mul_ui(QB, avg_min_sc_dist_B, mul);
        mpz_mod(QB, QB, M);
        mpz_powm_ui(E_QA, E_avg_min_sc_dist, mul, Msqr);
    }
    else
    {
        mpz_set(QB, avg_min_sc_dist_B);
        mpz_set(E_QA, E_avg_min_sc_dist_A);
    }
 
#ifdef pickRandom
    othPaillierObj->decrypt(E_rand, E_avg_min_sc_dist);
    gmp_printf("\n Q:%Zd", E_rand);
#endif
    
    RVS(E_rand, E_rand_A, rand_B, E_QA, QB);
    
    choose_center(E_rand_A, rand_B, centerNum);
    
    /*long int dpNum;
    mpz_set_ui(E_QA, NUMDATAPOINT);
    
    getRandomR(E_rand, E_QA);
    dpNum = mpz_get_ui(E_rand);
    for(long   int j =0; j<NUMATTRIBUTE; ++j)
    {
        mpz_set(E_center[centerNum][j], E_data_pt[dpNum][j]);
        mpz_set(E_center_A[centerNum][j], E_data_pt_A[dpNum][j]);
        mpz_set(center_B[centerNum][j], data_pt_B[dpNum][j]);
        
        write(E_center_A[centerNum][j]);
    }*/
    

    cout<<"\nCenter "<<centerNum<<" :";

    for(long   int j=0; j < NUMATTRIBUTE; ++j)
    {
        /*othPaillierObj->decrypt(E_rand, E_center_A[centerNum][j]);
        gmp_printf(" %Zd+%Zd=",E_rand, center_B[centerNum][j]);*/
        othPaillierObj->decrypt(E_rand, E_center[centerNum][j]);
        gmp_printf("%Zd ",E_rand);
    }
    
    mpz_clear(E_rand);
    mpz_clear(E_rand_A);
    mpz_clear(rand_B);
    mpz_clear(E_QA);
    mpz_clear(QB);
    
}



//Main protocols steps

void kMeansProtocol::InitAllClusterCenters()
{
    cout<<"\n\nInitialize cluster centers"<<endl;

    for(long   int i=0; i < NUMCENTERS; ++i)
    {
        mpz_init(fact_cent[i]);
        mpz_init(fact_cent_new[i]);
        
        for(long   int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_init(center_B[i][j]);
            mpz_init(E_center_A[i][j]);
            mpz_init(E_center[i][j]);
            
            mpz_init(center_new_B[i][j]);
            mpz_init(E_center_new_A[i][j]);
            mpz_init(E_center_new[i][j]);
        }
    }
    for(long   int j=0; j < NUMATTRIBUTE; ++j)
    {
        mpz_init(cog_B[j]);
        mpz_init(E_cog_A[j]);
        mpz_init(E_cog[j]);
    }
    
    for(long   int i=0; i < NUMDATAPOINT; ++i)
    {
        mpz_init(min_sc_dist_B[i]);
        mpz_init(E_min_sc_dist_A[i]);
        mpz_init(E_min_sc_dist[i]);
        
        mpz_init_set_ui(COGdist_B[i],0);
        mpz_init(E_COGdist_A[i]);
        mpz_init(E_COGdist[i]);
        othPaillierObj->encrypt(E_COGdist_A[i], COGdist_B[i]);
        othPaillierObj->encrypt(E_COGdist[i], COGdist_B[i]);

        /*mpz_init(min_dist_B[j]);
        mpz_init(E_min_dist_A[j]);
        mpz_init(E_min_dist[j]);*/
        
        for(long   int j=0; j < NUMCENTERS; ++j)
        {
            //distances
            /*mpz_init_set_ui(dist_B[i][j],0); //k+1 for cog
            mpz_init(E_dist_A[i][j]);
            mpz_init(E_dist[i][j]);
            
            othPaillierObj->encrypt(E_dist_A[i][j], dist_B[i][j]);
            othPaillierObj->encrypt(E_dist[i][j], dist_B[i][j]);*/

            mpz_init_set_ui(sc_dist_B[i][j],0);
            mpz_init(E_sc_dist_A[i][j]);
            mpz_init(E_sc_dist[i][j]);
            
            othPaillierObj->encrypt(E_sc_dist_A[i][j], sc_dist_B[i][j]);
            othPaillierObj->encrypt(E_sc_dist[i][j], sc_dist_B[i][j]);
        }
    }
    
    mpz_init(avg_min_sc_dist_B);
    mpz_init(E_avg_min_sc_dist_A);
    mpz_init(E_avg_min_sc_dist);
    mpz_init(avg_min_dist_B);
    mpz_init(E_avg_min_dist_A);
    mpz_init(E_avg_min_dist);

    
    
    
    CalcCOG();
    CalcDistToCOG_AvgSqrdDist();
    PickCenter(0);

    long   int mul = NUMDATAPOINT*NUMDATAPOINT;
    mpz_t E_min, E_min_A, min_B, temp;
    mpz_init(E_min);
    mpz_init(E_min_A);
    mpz_init(min_B);
    mpz_init(temp);
    
    //cout<<"\navg_min_sc_dist :"<<endl;
    for(long   int j = 1; j<NUMCENTERS; ++j)
    {
        mpz_set_ui(avg_min_sc_dist_B, 0);
        othPaillierObj->encrypt(E_avg_min_sc_dist_A, avg_min_sc_dist_B);
        othPaillierObj->encrypt(E_avg_min_sc_dist, avg_min_sc_dist_B);

        for(long   int i= 0; i<NUMDATAPOINT; ++i)
        {
            if (i>350)
                cout<<endl;
            
            comp_distance_short(E_sc_dist[i][j-1], E_sc_dist_A[i][j-1], sc_dist_B[i][j-1], i, j-1); //dist already scaled by n^2
            
            /*mpz_mul_ui(sc_dist_B[i][j-1], dist_B[i][j-1], mul);//scaled by n^2
            mpz_mod(sc_dist_B[i][j-1], sc_dist_B[i][j-1], M);
            
            mpz_pow_ui(E_sc_dist_A[i][j-1], E_dist_A[i][j-1], mul);//scaled by n^2
            mpz_mod(E_sc_dist_A[i][j-1], E_sc_dist_A[i][j-1], Msqr);
            
#ifdef initcent_test
            othPaillierObj->decrypt(temp, E_dist[i][j-1]);
            
#endif
            mpz_pow_ui(E_sc_dist[i][j-1], E_dist[i][j-1], mul);//scaled by n^2
            mpz_mod(E_sc_dist[i][j-1], E_sc_dist[i][j-1], Msqr);
            othPaillierObj->decrypt(temp, E_sc_dist[i][j-1]);*/
            
            MOTN_value(E_min, E_min_A, min_B, E_sc_dist_A[i][j-1], E_min_sc_dist_A[i], sc_dist_B[i][j-1], min_sc_dist_B[i]);
            mpz_set(min_sc_dist_B[i], min_B);
            mpz_set(E_min_sc_dist_A[i], E_min_A);
            mpz_set(E_min_sc_dist[i], E_min);
            
            /*mpz_set(min_sc_dist_B[i], sc_dist_B[i][j]);
            othPaillierObj->encrypt(E_min_sc_dist[i], min_sc_dist_B[i]);
            mpz_set(E_min_sc_dist_A[i], E_dist_A[i][j]);
            mpz_mul(E_min_sc_dist[i], E_min_sc_dist[i], E_min_sc_dist_A[i]);
            mpz_mod(E_min_sc_dist[i], E_min_sc_dist[i], Msqr);*/
            
            mpz_add(avg_min_sc_dist_B, avg_min_sc_dist_B, min_sc_dist_B[i]);//scaled by n^3 since we r not diving the dist
            mpz_mod(avg_min_sc_dist_B, avg_min_sc_dist_B, M);
            
            mpz_mul(E_avg_min_sc_dist_A, E_avg_min_sc_dist_A, E_min_sc_dist_A[i]);
            mpz_mod(E_avg_min_sc_dist_A, E_avg_min_sc_dist_A, Msqr);
                    
        }
        
        othPaillierObj->encrypt(E_min, avg_min_sc_dist_B);
        mpz_mul(E_avg_min_sc_dist, E_avg_min_sc_dist_A, E_min);
        mpz_mod(E_avg_min_sc_dist, E_avg_min_sc_dist, Msqr);

#ifdef test_init
        othPaillierObj->decrypt(temp, E_avg_min_sc_dist_A);
        gmp_printf("Before picking center avg_min_sc_dist: %Zd+%Zd=",temp,avg_min_sc_dist_B);
        othPaillierObj->decrypt(temp, E_avg_min_sc_dist);
        gmp_printf("%Zd\n",temp);
#endif
        PickCenter(j);
    
    }
    
    mpz_clear(E_min);
    mpz_clear(E_min_A);
    mpz_clear(min_B);
    mpz_clear(temp);
    
    for(long   int i=0; i < NUMCENTERS; ++i)
    {
        mpz_clear(fact_cent[i]);
        mpz_clear(fact_cent_new[i]);
    }
    for(long   int j=0; j < NUMATTRIBUTE; ++j)
    {
        mpz_clear(cog_B[j]);
        mpz_clear(E_cog_A[j]);
        mpz_clear(E_cog[j]);
    }
    for(long   int j=0; j < NUMDATAPOINT; ++j)
    {
        mpz_clear(min_sc_dist_B[j]);
        mpz_clear(E_min_sc_dist_A[j]);
        mpz_clear(E_min_sc_dist[j]);
        mpz_clear(COGdist_B[j]);
        mpz_clear(E_COGdist_A[j]);
        mpz_clear(E_COGdist[j]);
        
        /*mpz_clear(min_dist_B[j]);
        mpz_clear(E_min_dist_A[j]);
        mpz_clear(E_min_dist[j]);*/
    }
    
    mpz_clear(avg_min_sc_dist_B);
    mpz_clear(E_avg_min_sc_dist_A);
    mpz_clear(E_avg_min_sc_dist);
    mpz_clear(avg_min_dist_B);
    mpz_clear(E_avg_min_dist_A);
    mpz_clear(E_avg_min_dist);
    
    
    
    for(long   int j=0; j < NUMATTRIBUTE; ++j)
    {
        for(long   int i=0; i < NUMCENTERS; ++i)
        {
            mpz_init(sum_dp_new_B[i][j]);
            mpz_init(E_sum_dp_new_A[i][j]);
            mpz_init(E_sum_dp_new[i][j]);
            
            mpz_init(sum_dp_B[i][j]);
            mpz_init(E_sum_dp_A[i][j]);
            mpz_init(E_sum_dp[i][j]);
        }

        
    }
    
    for(long   int i=0; i < NUMCENTERS; ++i)
    {
        
        for(long   int j=0; j < NUMDATAPOINT; ++j)
        {
            
            mpz_init(closest_center_B[j][i]);//closest center details at the current iteration for each datapoint
            mpz_init(E_closest_center_A[j][i]);
            mpz_init(E_closest_center[j][i]);
            
        }

        
        mpz_init_set_ui(num_dp_new_B[i], 0);
        mpz_init(E_num_dp_new_A[i]);
        othPaillierObj->encrypt(E_num_dp_new_A[i], num_dp_new_B[i]);
        mpz_set_ui(num_dp_new_B[i], 1);
        mpz_init(E_num_dp_new[i]);
        othPaillierObj->encrypt(E_num_dp_new[i], num_dp_new_B[i]);
        
        mpz_init_set_ui(num_dp_B[i], 0);
        mpz_init(E_num_dp_A[i]);
        othPaillierObj->encrypt(E_num_dp_A[i], num_dp_B[i]);
        mpz_set_ui(num_dp_B[i], 1);
        mpz_init(E_num_dp[i]);
        othPaillierObj->encrypt(E_num_dp[i], num_dp_B[i]);

        mpz_init_set_ui(alpha_B[i],0);
        mpz_init(E_alpha_A[i]);
        othPaillierObj->encrypt(E_alpha_A[i], alpha_B[i]);
        mpz_set_ui(alpha_B[i],1);
        mpz_init(E_alpha[i]);
        othPaillierObj->encrypt(E_alpha[i], alpha_B[i]);
        
        mpz_init_set_ui(beta_B[i],0);
        mpz_init(E_beta_A[i]);
        othPaillierObj->encrypt(E_beta_A[i], beta_B[i]);
        mpz_set_ui(beta_B[i],1);
        mpz_init(E_beta[i]);
        othPaillierObj->encrypt(E_beta[i], beta_B[i]);
        
    }
    
    
    //scaling
    mpz_init(G_alpha_B);
    mpz_init(E_G_alpha_A);
    mpz_init(E_G_alpha);
    
    mpz_init_set_ui(G_alpha_B,0);
    mpz_init(E_G_alpha_A);
    othPaillierObj->encrypt(E_G_alpha_A, G_alpha_B);
    mpz_set_ui(G_alpha_B,1);
    mpz_init(E_G_alpha);
    othPaillierObj->encrypt(E_G_alpha, G_alpha_B);
    
    mpz_init(G_beta_B);
    mpz_init(E_G_beta_A);
    mpz_init(E_G_beta);

    
    
}



void kMeansProtocol::runLloydsStep(bool * b)
{
    cout<<"\nRun Lloyd's step"<<endl;

    mpz_t temp;
    mpz_init(temp);
    
    cout<<"\nIn Lloyd's step"<<endl;
    findClosestCenter();
    
    calcNewCenter();

    scaling_fact();
    
    termination_condition(b);
    
}

void kMeansProtocol::findClosestCenter()
{
    cout<<"\nFind closest center"<<endl;

    mpz_t temp;
    mpz_init(temp);
    
    for(long   int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long   int j=0; j < NUMCENTERS; ++j)
        {
            comp_scaled_dist(i, j);
            
        }
        
        MOMN_k( E_sc_dist[i], E_sc_dist_A[i], sc_dist_B[i], i);

        
    }
    
    for(long   int i=0; i < NUMDATAPOINT; ++i)
    {
    cout<<"\nClosest for "<<i<<" :";
    for(long   int j=0; j < NUMCENTERS; ++j)
    {
#ifdef closest_test
        othPaillierObj->decrypt(temp, E_closest_center_A[i][j]);
            gmp_printf(" %Zd+",temp);
        gmp_printf("%Zd=",closest_center_B[i][j]);
#endif
            othPaillierObj->decrypt(temp, E_closest_center[i][j]);
            gmp_printf("%Zd ",temp);
    }
        
    }
        
    mpz_clear(temp);
    
}

void kMeansProtocol::calcNewCenter()
{
    cout<<"\nCalculate new center"<<endl;

    mpz_t mul1, mul2, temp;
    mpz_init(mul1);
    mpz_init(mul2);
    mpz_init(temp);
    
    for(long   int j=0; j < NUMCENTERS; ++j)
    {
        cout<<"\nNew Center "<<j<<" :";

        for(long   int l=0; l < NUMATTRIBUTE; ++l)
        {
            mpz_set_ui(sum_dp_new_B[j][l], 0);
            othPaillierObj->encrypt(E_sum_dp_new[j][l], sum_dp_new_B[j][l]);

            for(long   int i=0; i < NUMDATAPOINT; ++i)
            {
                mpz_powm(mul1, E_closest_center[i][j], data_pt_B[i][l], Msqr);
                //othPaillierObj->decrypt(temp, mul1);

                mpz_powm(mul2, E_data_pt_A[i][l], closest_center_B[i][j], Msqr);
                //othPaillierObj->decrypt(temp, mul2);

                mpz_mul(E_sum_dp_new[j][l], E_sum_dp_new[j][l], mul1);
                mpz_mod(E_sum_dp_new[j][l], E_sum_dp_new[j][l], Msqr);

                //othPaillierObj->decrypt(temp, E_sum_dp_new[j][l]);

                mpz_mul(E_sum_dp_new[j][l], E_sum_dp_new[j][l], mul2);
                mpz_mod(E_sum_dp_new[j][l], E_sum_dp_new[j][l], Msqr);

                //othPaillierObj->decrypt(temp, E_sum_dp_new[j][l]);

                
               /* read(E_sum_dp_new_A[j][l]);
                mpz_mul(mul1, E_sum_dp_new_A[j][l], E_sum_dp_new[j][l] );
                mpz_mod(mul1, mul1, Msqr);

                othPaillierObj->decrypt(temp, mul1);
                gmp_printf("%Zd ",temp);*/

            }
                       
            read(E_sum_dp_new_A[j][l]);
            mpz_mul(E_sum_dp_new[j][l], E_sum_dp_new_A[j][l], E_sum_dp_new[j][l]);
            
            share_write(E_sum_dp_new_A[j][l], sum_dp_new_B[j][l], E_sum_dp_new[j][l]);

            //othPaillierObj->decrypt(mul1, E_sum_dp_new[j][l]);

            mpz_set(center_new_B[j][l] , sum_dp_new_B[j][l]);
            mpz_set(E_center_new_A[j][l] , E_sum_dp_new_A[j][l]);
            mpz_set(E_center_new[j][l] , E_sum_dp_new[j][l]);
            
                    othPaillierObj->decrypt(mul1, E_center_new[j][l]);
                    gmp_printf("%Zd ",mul1);

        }
    }
    
    
    for(long   int j=0; j < NUMCENTERS; ++j)
    {
        cout<<"\n no of data points in Center "<<j<<" :";

           mpz_set_ui(num_dp_new_B[j], 0);
            othPaillierObj->encrypt(E_num_dp_new[j], num_dp_new_B[j]);

            for(long   int i=0; i < NUMDATAPOINT; ++i)
        {
            mpz_mul(E_num_dp_new[j], E_num_dp_new[j], E_closest_center[i][j]);
            mpz_mod(E_num_dp_new[j], E_num_dp_new[j],Msqr);

            /*othPaillierObj->decrypt(mul1, E_closest_center[i][j]);

            othPaillierObj->decrypt(mul1, E_num_dp_new[j]);
            gmp_printf("%Zd ",mul1);*/

        }
        
        mpz_mod(E_num_dp_new[j], E_num_dp_new[j],Msqr);
        
        share_write(E_num_dp_new_A[j], num_dp_new_B[j], E_num_dp_new[j]);
        
        /*othPaillierObj->decrypt(mul1, E_num_dp_new_A[j]);

        gmp_printf(" %Zd+%Zd=",mul1, num_dp_new_B[j]);*/
         othPaillierObj->decrypt(mul1, E_num_dp_new[j]);
         gmp_printf("%Zd ",mul1);

    }


    mpz_clear(mul1);
    mpz_clear(mul2);
}
