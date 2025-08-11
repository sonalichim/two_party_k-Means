//
//  main.cpp
//  k-means
//
//  Created by sonali chim on 22/08/18.
//  Copyright © 2018 personal. All rights reserved.
//


#include "kMeansProtocol.h"

using namespace std;
using namespace std::chrono;

connection * connWithBob;
int main()
{
    mpz_t Modulus, pk;
    mpz_init(Modulus);
    mpz_init(pk);
    //string strValue;
    char * strValue;
    
    kMeansProtocol *kMeansObj = new kMeansProtocol(LValue);
    kMeansObj->getMyModulus(Modulus);
    
    connWithBob =  new connection();
    connWithBob->connectToServer();
    strValue = mpz_get_str(NULL, 10, Modulus);
    connWithBob->writeDataWithLength(&strValue);
    free(strValue);
    
    kMeansObj->getPrivKey(pk);
    strValue = mpz_get_str(NULL, 10, pk); //memory allocated
    connWithBob->writeDataWithLength(&strValue);
    free(strValue);
      
/*    mpz_t PB[4], DA, QA, temp;
     bool b[KValue];
    
    mpz_init_set_ui(PB[0], 2);
    mpz_init_set_ui(PB[1], 2);
    mpz_init_set_ui(PB[2], 2);
    mpz_init_set_ui(PB[3], 2);
    
    mpz_init_set_ui(QA, 8);
    
    kMeansObj->RVS(temp, QA);*/

    
    kMeansObj->initDataPoints();
    kMeansObj->shareDataPoints();
    bool op=true;
    int i = 0;
    time_point<high_resolution_clock> start = high_resolution_clock::now();

    kMeansObj->InitAllClusterCenters();
    
    time_point<high_resolution_clock> stop_init = high_resolution_clock::now();
    auto duration_sec_init = duration_cast<microseconds>(stop_init-start);
    cout<<endl<<"Time of init :"<<duration_sec_init.count()<<" microseconds"<<endl;
    
    while(i<2)
    {
        time_point<high_resolution_clock> start_iter = high_resolution_clock::now();

        cout<<endl<<"Iteration "<<i++<<endl;
        kMeansObj->runLloydsStep(&op);
        
        time_point<high_resolution_clock> stop = high_resolution_clock::now();
        auto duration_sec = duration_cast<microseconds>(stop-start_iter);
        cout<<endl<<"Time of iteration :"<<duration_sec.count()<<" microseconds"<<endl;
    }
    
    time_point<high_resolution_clock> stop = high_resolution_clock::now();
    auto duration_sec = duration_cast<microseconds>(stop-start);
    cout<<endl<<"Time of total clusering(microsec) :"<<duration_sec.count()<<" microseconds"<<endl;

    cout<<endl;
    mpz_clear(Modulus);
    //delete kMeansObj;
    delete connWithBob;
    return 0;
    
}

