
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

//Init
void kMeansProtocol::readOwnData()
{
    cout<<"\nReading own data"<<endl;
    string word;
    double var;
    
    fstream file;
    file.open ("/Users/sonalichim/Dropbox/Code/Sonali/k-Means/k-means_A/k-means_Alice/wdbc.data");
    //file.open ("./k-means_Alice/wdbc.data");
    
    if(! file.is_open())
        cout<<"%%%%%FILE OPEN ERROR";
    
    for(long  int i=0; i < NUMMYDATAPOINT; ++i)
    {
        //cout<<endl<<i<<")";

        for(long  int u=0; u < STARTSKIP; ++u)
        file >> word;
        
        for(long  int j=0; j < NUMATTRIBUTE; ++j)
        {
            file >> word;
            var = atof(word.c_str());
            var=var*SCALING;//add scaling
            mpz_set_ui(data_pt_A[i][j], var);
            
            //gmp_printf(" %Zd ",data_pt_A[i][j]);

        }
        for(long  int u=NUMATTRIBUTE; u< NUMATTRIBUTE+TAILSKIP; ++u)
        file >> word;
        
    }
    
    //closing the file
    file.close();
}


void kMeansProtocol::shareDataPoints()
{
    cout<<"\nSharing own data"<<endl;

    mpz_t ciph;
    mpz_init(ciph);
    char * strValue;
    char * strValue_t;

    strValue_t = (char *) malloc(DATASIZE);
    memset(strValue_t, '\0' ,DATASIZE);

    writeCount=0;
    
    long int count =0;
    long int temp;
    sleep(SLEEPTIME);
    for(long  int i=0; i<NUMDATAPOINT; ++i)
    {
        cout<<endl<<i;

        //multiply each attribute of the datapoint by corresponding LA share
        for(long  int j=0; j<NUMATTRIBUTE; ++j)
        {
            //temp=(i+1)*(j+1);

            //share with Bob
            
            /*myPaillierObj->encrypt(ciph, data_pt_A[i][j]);
            strValue =  mpz_get_str(NULL, 10, ciph);
            connWithBob->writeDataWithLength(&strValue);
            free(strValue);*/
            
            encrypt_write(data_pt_A[i][j]);
            gmp_printf(" %Zd ",data_pt_A[i][j]);
            
            //count++;
            /*memset(strValue_t, '\0' ,DATASIZE);
            strcpy(strValue_t, (to_string(count)).c_str());
            connWithBob->writeDataWithLength(&strValue_t);*/
        
           // count++;

        }
    }
    free(strValue_t);

    mpz_clear(ciph);
}

void kMeansProtocol::initDataPoints()
{
    cout<<"\nInit dp"<<endl;

    for(long  int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long  int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_init_set_ui(data_pt_A[i][j],0);
        }
    }
    
    /*mpz_set_ui(data_pt_A[0][0], 1);
     mpz_set_ui(data_pt_A[0][1], 2);
     mpz_set_ui(data_pt_A[1][0], 3);
     mpz_set_ui(data_pt_A[1][1], 2);
     mpz_set_ui(data_pt_A[4][0], 6);
     mpz_set_ui(data_pt_A[4][1], 6);
     mpz_set_ui(data_pt_A[5][0], 6);
     mpz_set_ui(data_pt_A[5][1], 4);*/
    
    readOwnData();
    
    /*  for(long  int i=0; i < NUMCENTERS; ++i)
     {
     for(long  int j=0; j < NUMATTRIBUTE; ++j)
     {
     mpz_set_ui(center_new_A[i][j],0);
     mpz_set_ui(center_A[i][j], 0);
     }
     }
     
     for(long  int j=0; j < NUMATTRIBUTE; ++j)
     {
     mpz_set_ui(cog_A[j],0);
     }
     */
    /*for(long  int i = 0; i<L; ++i)
     mpz_init(GGA[i]);*/
    
    //mpz_init(CBar);
}

void kMeansProtocol::getMyModulus(mpz_t modulus)
{
    mpz_set(modulus, M);
}

void kMeansProtocol::getPrivKey(mpz_t priv)
{
    myPaillierObj->getPhi(priv);
}

kMeansProtocol::kMeansProtocol( long  int x)
{
    L = x;
    myPaillierObj = new PaillierSystem(x/2);
    mpz_init(M);
    mpz_init(Msqr);
    myPaillierObj->getPublicKey(M);
    mpz_mul(Msqr, M, M);
}

kMeansProtocol::kMeansProtocol()
{
    
}

kMeansProtocol::~kMeansProtocol()
{
    delete connWithBob;
    delete myPaillierObj;
    //delete othPaillierObj;
    mpz_clear(M);
    mpz_clear(Msqr);
    
    for(long  int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long  int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_clear(data_pt_A[i][j]);
        }
    }
    for(long  int i=0; i < NUMCENTERS; ++i)
    {
        for(long  int j=0; j < NUMATTRIBUTE; ++j)
        {
            mpz_clear(center_new_A[i][j]);
        }
    }
    
    //for(int i = 0; i<M; ++i)
    //  mpz_clear(GGA[i]);
}




//Main protocol steps
void kMeansProtocol::InitAllClusterCenters()
{
    cout<<"\n\nInit cluster centers"<<endl;

    for(long  int j =0; j<NUMCENTERS; ++j)
    {
        for(long  int i =0; i<NUMATTRIBUTE; ++i)
        {
            //centers
            mpz_init(center_A[j][i]); //center points share last iteration
            
            mpz_init(center_new_A[j][i]); //center points share scurr iteration
            
            mpz_init(sum_dp_new_A[j][i]); //sum of dp in a center in prev iter
            mpz_init(sum_dp_A[j][i]); //sum of dp in a center in prev iter
        }
        

        for(long  int i =0; i<NUMDATAPOINT; ++i)
        {
            mpz_init(dist_A[i][j]); //dist of a dp from a center
            mpz_init(sc_dist_A[i][j]); //scaled dist of a dp from a center
            //closest center
            mpz_init(closest_center_A[i][j]);//closest center details at the current iteration for each datapoint
            
        }
        
        mpz_init(num_dp_A[j]);//no of dp in a center in prev iter
        mpz_init(num_dp_new_A[j]); //no of dp in a center in curr iter
        
        mpz_init_set_ui(alpha_A[j],0); //
        
        mpz_init_set_ui(beta_A[j],0); //
                
    }
    
    for(long  int i =0; i<NUMDATAPOINT; ++i)
    {
        //distances
        mpz_init(COGdist_A[i]); //dist of a dp frm COG
        
        mpz_init(min_dist_A[i]); //min of the dist of a dp from all centers+COg
        
        mpz_init(min_sc_dist_A[i]); //min of the scaled dist of a dp from all centers+COg
        
    }

    for(long  int i =0; i<NUMATTRIBUTE; ++i)
        mpz_init(cog_A[i]); //center of gravity
    
    mpz_init(avg_min_dist_A); //avg of the abv min dist of all dp
    
    mpz_init(avg_min_sc_dist_A); //avg of the abv min scaled dist of all dp
    
    //scaling
    mpz_init_set_ui(G_alpha_A,0); //
    
    mpz_init_set_ui(G_beta_A,0); //
    
    

    CalcCOG();
    CalcDistToCOG_AvgSqrdDist();
    
    PickCenter(0);
    
    
    long  int mul = NUMDATAPOINT*NUMDATAPOINT;
    mpz_set_ui(avg_min_sc_dist_A,0);

    for(long  int j = 1; j<NUMCENTERS; ++j)
    {
        mpz_set_ui(avg_min_sc_dist_A, 0);

        for(long  int i= 0; i<NUMDATAPOINT; ++i)
        {
            comp_distance_short(dist_A[i][j-1], i, j-1);
            mpz_mul_ui(sc_dist_A[i][j-1], dist_A[i][j-1], mul);//scaled by n^2
            mpz_mod(sc_dist_A[i][j-1], sc_dist_A[i][j-1], M);
            
            MOTN_value(min_sc_dist_A[i], sc_dist_A[i][j-1], min_sc_dist_A[i]);
            
            mpz_add(avg_min_sc_dist_A, avg_min_sc_dist_A, min_sc_dist_A[i]);
            mpz_mod(avg_min_sc_dist_A, avg_min_sc_dist_A, M);
            
        }
        PickCenter(j);

    }
    
    

    for(long  int j=0; j < NUMATTRIBUTE; ++j)
    {
        mpz_clear(cog_A[j]);
    }
    for(long  int j=0; j < NUMDATAPOINT; ++j)
    {
        mpz_clear(min_sc_dist_A[j]);
        mpz_clear(COGdist_A[j]);
        
        /*mpz_clear(min_dist_B[j]);
        mpz_clear(E_min_dist_A[j]);
        mpz_clear(E_min_dist[j]);*/
    }
    
    mpz_clear(avg_min_sc_dist_A);
    mpz_clear(avg_min_dist_A);

}



void kMeansProtocol::CalcCOG()
{
    cout<<"\nCalc COG"<<endl;

    for(long  int j=0; j <  NUMATTRIBUTE; ++j)
    {
        mpz_set_ui(cog_A[j],0);
        
        for(long  int i=0; i < NUMDATAPOINT; ++i)
        {
            //othPaillierObj->decrypt(temp2, E_data_pt[i][j]);
            if(mpz_cmp_ui(data_pt_A[i][j],0) !=0)
                mpz_add(cog_A[j], cog_A[j], data_pt_A[i][j]);
            
            //othPaillierObj->decrypt(temp2, sum);
        }
        
        //gmp_printf("\nSum CoG_B(%d)=%Zd", j,w);
        
        //gmp_printf("  A Sum=%Zd", cog_A[j]);
        
        encrypt_write(cog_A[j]);
        
        read_decrypt(cog_A[j]);
    }
    
    //mpz_set_ui(centerOfGravity[0],51); for mod 143
    //mpz_set_ui(centerOfGravity[1],81);
}

void kMeansProtocol::CalcDistToCOG_AvgSqrdDist()
{
    cout<<"\nCalc dist to COG"<<endl;

    mpz_t mul, w;
    mpz_init(mul);
    mpz_init(w);
    
    mpz_set_ui(avg_min_sc_dist_A,0);
    
    for(long  int i =0; i<NUMDATAPOINT; ++i)
    {
        for(long  int j =0; j<NUMATTRIBUTE; ++j)
        {
            read_decrypt(w);
            mpz_mul(mul, w, w);
            encrypt_write(mul);
        }
        
        read_decrypt(COGdist_A[i]);
        mpz_set(min_sc_dist_A[i], COGdist_A[i]);
        
        mpz_add(avg_min_sc_dist_A, avg_min_sc_dist_A, min_sc_dist_A[i]);
        mpz_mod(avg_min_sc_dist_A, avg_min_sc_dist_A, M);
    }
    
    mpz_clear(mul);
    mpz_clear(w);
}






void kMeansProtocol::PickCenter(long  int centerNum)
{
    cout<<"\nPick center "<<centerNum<<endl;

    mpz_t rand_A,  QA;
    mpz_init(rand_A);
    mpz_init(QA);
    
    if(centerNum == 0)
    {
        long  int mul = NUMDATAPOINT*2;
        mpz_mul_ui(QA, avg_min_sc_dist_A, mul);
        mpz_mod(QA, QA, M);

    }
    else
    {
        mpz_set(QA, avg_min_sc_dist_A);
    }
 


    RVS(rand_A, QA);
    
    choose_center(rand_A, centerNum);
    /*for(long   int j =0; j<NUMATTRIBUTE; ++j)
    {
        read_decrypt(center_A[centerNum][j]);
    }*/
    
    
    
    mpz_clear(rand_A);
    mpz_clear(QA);
    


}



void kMeansProtocol::runLloydsStep(bool * b)
{
    cout<<"\nIn Lloyd's step"<<endl;
    
    findClosestCenter();
    
    calcNewCenter();
    
    scaling_fact();
    
    termination_condition(b);
    
}

void kMeansProtocol::findClosestCenter()
{
    cout<<"\nFind closest center"<<endl;

    for(long  int i=0; i < NUMDATAPOINT; ++i)
    {
        for(long  int j=0; j < NUMCENTERS; ++j)
        {
            comp_scaled_dist(i, j);
            
        }
        
        MOMN_k(i);
        
    }
    
}

void kMeansProtocol::calcNewCenter()
{
    cout<<"\nCalc new center"<<endl;

    mpz_t mul;
    mpz_init(mul);
    
    for(long  int j=0; j < NUMCENTERS; ++j)
    {
        
    for(long  int l=0; l < NUMATTRIBUTE; ++l)
    {
        mpz_set_ui(sum_dp_new_A[j][l], 0);
        for(long  int i=0; i < NUMDATAPOINT; ++i)
        {
            mpz_mul(mul, closest_center_A[i][j], data_pt_A[i][l]);
            mpz_mod(mul, mul, M);

            mpz_add(sum_dp_new_A[j][l], sum_dp_new_A[j][l], mul);
            mpz_mod(sum_dp_new_A[j][l], sum_dp_new_A[j][l], M);
            
            //encrypt_write(sum_dp_new_A[j][l]);
        
        }
        
        encrypt_write(sum_dp_new_A[j][l]);
        
        read_decrypt(sum_dp_new_A[j][l]);
        mpz_set(center_new_A[j][l] , sum_dp_new_A[j][l]);
    }
        

    }
    
    for(long   int j=0; j < NUMCENTERS; ++j)
            read_decrypt(num_dp_new_A[j]);

    
    
    mpz_clear(mul);
}
