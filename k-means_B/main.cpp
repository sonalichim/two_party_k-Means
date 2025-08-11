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

connection * connWithAlice;

void divP(mpz_t QB, mpz_t PB, mpz_t DB);
void compGamma_int(mpz_t * GB, mpz_t DB);


int main()
{
    mpz_t othModulus, othSK;
    mpz_init(othModulus);
    mpz_init(othSK);
    char * strValue;
    
    /*struct timespec st;
    clock_gettime(CLOCK_REALTIME, &st);
    cout<<"\n"<<st.tv_nsec<<endl;*/
    
    connWithAlice = new connection();
    connWithAlice->getReady();
    connWithAlice->readDataWithLength(&strValue);
    mpz_set_str(othModulus, strValue, 10);
    free(strValue);
    gmp_printf ("\nAlice's PubKey =%Zd \n", othModulus);
    
    connWithAlice->readDataWithLength(&strValue);
    mpz_set_str(othSK, strValue, 10);
    
    kMeansProtocol *kMeansObj = new kMeansProtocol(LValue, othModulus, othSK);
        
    
/*     mpz_t PB[LValue], E_R,E_RA, RB, PA[LValue], temp, QA, QB;
     bool b[LValue];
        
    mpz_init_set_ui(PB[0], 1);
    mpz_init_set_ui(PB[1], 0);
    mpz_init_set_ui(PB[2], 0);
    mpz_init_set_ui(PB[3], 1);
    mpz_init(PA[0]);
    mpz_init(PA[1]);
    mpz_init(PA[2]);
    mpz_init(E_RA);
    mpz_init(E_R);
    mpz_init(RB);
    mpz_init(QB);
    
    mpz_init_set_ui(temp, 8);
    kMeansObj->othPaillierObj->encrypt(QA, temp);
    mpz_init_set_ui(QB, 2);

    
    kMeansObj->RVS(E_R, E_RA, RB, QA, QB);*/
    
    bool op=true;
    int i = 0;
    /*  long sigma[LValue];
    kMeansObj->reordering(sigma);*/

    kMeansObj->initDataPoints();
    kMeansObj->readEncDataPts();
    
    time_point<high_resolution_clock> start = high_resolution_clock::now();

    kMeansObj->InitAllClusterCenters();
    
    time_point<high_resolution_clock> stop_init = high_resolution_clock::now();
    auto duration_sec_init = duration_cast<microseconds>(stop_init-start);
    cout<<endl<<"Time for init :"<<duration_sec_init.count()<<" microseconds"<<endl;
    
    while(i<2)
    {
        time_point<high_resolution_clock> start_iter = high_resolution_clock::now();

        cout<<endl<<"Iteration :"<<i++<<endl;
        kMeansObj->runLloydsStep(&op);
        
        time_point<high_resolution_clock> stop = high_resolution_clock::now();
        auto duration_sec = duration_cast<microseconds>(stop-start_iter);
        cout<<endl<<"Time of this iteration :"<<duration_sec.count()<<" microseconds"<<endl;
    }
    time_point<high_resolution_clock> stop = high_resolution_clock::now();
    auto duration_sec = duration_cast<microseconds>(stop-start);
    cout<<endl<<"Time of total clusering :"<<duration_sec.count()<<" microseconds"<<endl;

    cout<<endl;
    mpz_clear(othModulus);
    //delete kMeansObj;
    delete connWithAlice;
    return 0;
}


