//
//  server.c
//  k-means_Bob
//
//  Created by sonali chim on 14/12/18.
//  Copyright © 2018 personal. All rights reserved.
//

#include "server.h"


connection::connection()
{
    af= AF_INET; type =  SOCK_STREAM; protocol = 0; port = 2220;
    clientConn = 0;
    mySocket = 0;
    memset(&addressOfServer, '0', sizeof(addressOfServer));
}

void connection::getReady()
{
    socklen_t adressLen = sizeof(addressOfServer);
    
    // Creating socket file descriptor
    if ((mySocket = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    {
        perror("socket failed");
        exit(EXIT_FAILURE);
    }
    
    //    int opt = 1;
    //    // Forcefully attaching socket to the port 8080
    //    if (setsockopt(mySocket, SOL_SOCKET, SO_REUSEADDR | SO_REUSEPORT, &opt, sizeof(opt)))
    //    {
    //        perror("setsockopt");
    //        exit(EXIT_FAILURE);
    //    }
    
    addressOfServer.sin_family = af;
    addressOfServer.sin_addr.s_addr = htonl(INADDR_ANY);
    addressOfServer.sin_port = htons(port);         // this is the port number of running server
    
    if (bind(mySocket, (const struct sockaddr *)&addressOfServer, adressLen)<0)
    {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }
    
    
    //socket created by socket function, bind assigns address to it, it assigns adress in 2nd argument to socket created which is 1st argument, 3rd arg is size of adress in bytes,
    //returns 0 if success else -1
    
    std::cout<<endl<<"Server listening.."<<std::endl;
    
    if (listen(mySocket, MSGLENGHTSIZE) < 0)
    {
        perror("listen");
        exit(EXIT_FAILURE);
    }
    
    //1st arg is sock provided by socket function, 2nd arg is no of connections a system can handle at a time
    if ((clientConn = accept(mySocket, (struct sockaddr *)&addressOfServer, &adressLen))<0)
    {
        perror("accept");
        exit(EXIT_FAILURE);
    }
    
    //cout << "Server stopped listening" << endl;
}

//void connection::getReady()
//{
//    int adressLen = sizeof(addressOfServer);
//    mySocket = socket(af, type, protocol); // creating socket, returns file descriptor
//
//    bind(mySocket, (struct sockaddr*)&addressOfServer , adressLen);
//    //socket created by socket function, bind assigns address to it, it assigns adress in 2nd argument to socket created which is 1st argument, 3rd arg is size of adress in bytes,
//    //returns 0 if success else -1
//
//    std::cout<<endl<<"Server listening.."<<std::endl;
//
//    while(listen(mySocket , 5) < 0); //BACKLOG = 5
//    //1st arg is sock provided by socket function, 2nd arg is no of connections a system can handle at a time
//
//    do
//    {
//        // whenever a request from client came. It will be processed here.
//        clientConn = accept(mySocket, (struct sockaddr*) &addressOfServer, (socklen_t *) adressLen);
//        //waits for connection function run by client, returns a socket file descriptor
//    }while(clientConn < 0);
//}

void connection::readData(char** msg)
{
    
    if (readMsgPool.size() >0)
    {
        string * temp = readMsgPool.at(0);
        *msg = (char *) malloc(temp->length() + 1);
        strcpy(*msg, temp->c_str());
        
        delete readMsgPool.at(0);
        readMsgPool.erase(readMsgPool.begin());
    }
    else
    {
        char dataReceived[DATASIZE];
        memset(dataReceived, '\0' ,DATASIZE);
        int n = 0;
        n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        }
        
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        }
        if(n <=0)
        {
            printf("\nError on read. Stopping now..\n");
            exit(0);
        }
        
        string temp;
        temp.append((char *) dataReceived); //append rcvd data here to temp
        
        *msg = (char *) malloc(temp.length() + 1); //assign  mem to output
        memset(*msg, '\0' , temp.length()+1); //fill with null char
        
        strcpy(*msg, temp.c_str());//copy the content
        printf("\n%d:%s",readCount++, temp.c_str());
        
    }
    
}

void connection::readnData(long int num) //allocated memory using new....please use dleete after use
{
    char dataReceived[DATASIZE];
    memset(dataReceived, '\0' ,DATASIZE);
    int n = 0;
    
    for(int i=0; i<num; ++i)
    {
        memset(dataReceived, '\0' ,DATASIZE);
        
        n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        }
        
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
        }
        if(n <=0)
        {
            printf("\n readnData:Error on read. Stopping now..\n");
            exit(0);
        }
        
        string temp;
        temp.append((char *) dataReceived);
        
        string * msg1 = new string(temp.c_str());
        readMsgPool.push_back(msg1);
        printf("\n%d:%s",i, msg1->c_str());
    }
}


void connection::readDataWithLength(char ** msg)
{
#ifdef test_conn
    cout<<endl<<"readCOunt="<<readCount++<<endl;
#endif
    
    if (readMsgPool.size() >0 && !(readMsgPool.size()==1 && isMsgComplete==0))
    {
        
        string * temp = readMsgPool.at(0);
        *msg = (char *) malloc(temp->length() + 1);
        strcpy(*msg, temp->c_str());
        
        
        
#ifdef test_conn_deep
        printf("\n%d:%s",readCount++, temp->c_str());
#endif
        delete readMsgPool.at(0);
        readMsgPool.erase(readMsgPool.begin());
    }
    else
    {
        char dataReceived[DATASIZE];
        memset(dataReceived, '\0' ,DATASIZE);
        int n = 0;
        atLeastOneFullMsgReceived=0;
        
        //while(atLeastOneFullMsgReceived ==0 )
        //{
            n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
            if(n <=0)
            {
                printf("\nError on read. Trying again..\n");
                n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
            }
            
            if(n <=0)
            {
                printf("\nError on read. Trying again..\n");
                n = (int)read(clientConn, dataReceived, sizeof(dataReceived)-1);
            }
            if(n <=0)
            {
                printf("\nError on read. Stopping now..\n");
                exit(0);
            }
            
#ifdef test_conn_deep
            printf("\n\n%s\n\n", dataReceived);
#endif
            string temp;
            
            if(readMsgPool.size()==1 && isMsgComplete==0)//chk for last msg is complete or not
            {
                temp.append((char*)RemMsg->c_str());
                
                isMsgComplete=1;
                delete RemMsg;
#ifdef test_conn_deep
                printf("\n\n%s\n\n", temp.c_str());
#endif
            }
            
            temp.append((char *) dataReceived);
            
#ifdef test_conn_deep
            printf("\n\n%s\n\n", temp.c_str());
#endif
            while (temp.size() >= MSGLENGHTSIZE)
            {
                int msgLen =  stoi(temp.substr(0, MSGLENGHTSIZE));
                string * msg1 = new string(temp.substr(MSGLENGHTSIZE, msgLen));
                
                if(msgLen > msg1->size())
                {
                    isMsgComplete=0;
                    RemMsg= new string(temp.c_str());
                }
                else
                {
                    readMsgPool.push_back(msg1);
                }
                /*msgLen.clear();
                 msgLen.append(msg->substr( 5+msg1->size(), msg->size()));*/
                
                temp.erase(0, MSGLENGHTSIZE+msg1->size());
            }
            
            if(temp.size() !=0)
            {
                isMsgComplete=0;
                RemMsg= new string(temp.c_str());
            }
            
            if (readMsgPool.size() >0)
            {
                atLeastOneFullMsgReceived = 1;
                
                temp = (readMsgPool.at(0))->c_str();
                *msg = (char *) malloc(temp.length() + 1);
                strcpy(*msg, temp.c_str());
                
#ifdef test_conn_deep
                
                printf("\n%d:%s",readCount++, temp.c_str());
#endif
                
                delete readMsgPool.at(0);
                readMsgPool.erase(readMsgPool.begin());
                
                //std::cout<<endl<<"Data received from Bob"<<msg->c_str()<<std::endl;

            }
        else
        {
            *msg = (char *) malloc(2);
            strcpy(*msg, "0");

        }
        }
    //}
    
}


void connection::writeDataWithlength(char** msg)
{
#ifdef test_conn
    cout<<endl<<"writeCount="<<writeCount++<<endl;
#endif
    //char dataSending[DATASIZE]; // Actually this is called packet in Network Communication, which contain data and send through.
    //time_t clock;
    //clock = time(NULL);
    //snprintf(dataSending, sizeof(dataSending), "%.24s\r\n", ctime(&clock)); // Printing successful message
    
    string lenPlusMsg, temp;
    temp.append((char*)*msg) ;
    lenPlusMsg.append(to_string(temp.size()));
    int len = (int)lenPlusMsg.size(); //lenght of msg lenght
    lenPlusMsg.clear();
    
    while((MSGLENGHTSIZE - len) >0)
    {
        lenPlusMsg.append("0");
        ++len;
    }
    lenPlusMsg.append(to_string(temp.size()));
    lenPlusMsg.append(temp.c_str());
    
    write(clientConn, lenPlusMsg.c_str(), lenPlusMsg.size());
    
    
    //std::cout<<"Data sent to Alice :"<<lenPlusMsg.c_str()<<std::endl;
}

void connection::writeData(char** msg)
{
    
    //char dataSending[DATASIZE]; // Actually this is called packet in Network Communication, which contain data and send through.
    //time_t clock;
    //clock = time(NULL);
    //snprintf(dataSending, sizeof(dataSending), "%.24s\r\n", ctime(&clock)); // Printing successful message
    
    string Message;
    Message.append((char*)*msg) ;
    
    write(clientConn, Message.c_str(), Message.size());
    //std::cout<<"Data sent to Alice :"<<lenPlusMsg.c_str()<<std::endl;
}


connection::~connection()
{
    close(clientConn);
}
