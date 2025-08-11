//
//  server.h
//  k-means_Bob
//
//  Created by sonali chim on 14/12/18.
//  Copyright © 2018 personal. All rights reserved.
//

#ifndef server_h
#define server_h

#include <iostream> // standard input and output library
#include <stdlib.h> // this includes functions regarding memory allocation
#include <string> // contains string functions
#include <string.h>
#include <errno.h> //It defines macros for reporting and retrieving error conditions through error codes
#include <time.h> //contains various functions for manipulating date and time
#include <unistd.h> //contains various constants
#include <sys/types.h> //contains a number of basic derived types that should be used whenever appropriate
#include <arpa/inet.h> // defines in_addr structure
#include <sys/socket.h> // for socket creation
#include <netinet/in.h> //contains constants and structures needed for internet domain addresses
#include <vector>

//#define test_conn_deep

#define BACKLOG = 1 //since only one party
#define DATASIZE 6144 //2*3072//4096 //40KB
#define MSGLENGHTSIZE 5 //1024 bits 309 digits //The largest value a 4096-bit integer can represent is 2^4096, which happens to have 1234 digits in decimal//2^6144: 1850 decimal digits requires 4 digits

/*void error(const char* msg)
{
    perror(msg);
    exit(0);
}*/

using namespace std;
using std::vector;

static int writeCount=-1;
static int readCount=-2;
static bool isMsgComplete=1;//1 yes 0 no
static string * RemMsg = NULL;
static bool atLeastOneFullMsgReceived = 0;

class connection
{
    int af, type, protocol, clientConn, port;
    struct sockaddr_in addressOfServer;
    //vector<string*> readMsgPool;
    vector<string*> readMsgPool;
    
    public:
        int mySocket;
        connection();
        void getReady();
        void readData(char ** msg);
        void readnData(long int num);
        void writeData(char ** msg);

        void writeDataWithlength(char ** msg);
    void readDataWithLength(char ** msg);
        ~connection();
};
#endif /* server_h */
