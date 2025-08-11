//
//  client.h
//  k-means_Alice
//
//  Created by sonali chim on 14/12/18.
//  Copyright © 2018 personal. All rights reserved.
//

#ifndef client_h
#define client_h

#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>

#include <iostream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h>
#include <time.h>
#include <vector>

//#define test_conn_deep

#define DATASIZE 6144//2*3072//4096 //40KB
#define MSGLENGHTSIZE 5 //1024 bits 309 digits //The largest value a 4096-bit integer can represent is 2^4096, which happens to have 1234 digits in decimal

using namespace std;
using std::vector;

static int writeCount=-2;
static int readCount=-1;
static bool isMsgComplete=1;//1 yes 0 no
static string * RemMsg = NULL;

/*void error(const char* msg)
{
    perror(msg);
    exit(0);
}*/

class connection
{
    int af, type, protocol, port;
    struct sockaddr_in ipOfServer;
    vector<string*> readMsgPool;
    
public:
    int ServerSocket;
    connection();
    void connectToServer();
    
    void readDataWithLength(char ** msg);
    void writeDataWithLength(char ** msg);

    void readData(char ** msg);
    void readnData(int n);
    void writeData(char ** msg);
    
    ~connection();
    
};


#endif /* client_h */
