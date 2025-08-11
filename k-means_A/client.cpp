//
//  client.c
//  k-means_Alice
//
//  Created by sonali chim on 14/12/18.
//  Copyright © 2018 personal. All rights reserved.
//

#include "client.h"
connection::connection()
{
    af= AF_INET; type =  SOCK_STREAM; protocol = 0; port = 2220;
    ServerSocket = 0;
    
    memset(&ipOfServer, '0', sizeof(ipOfServer));
    
    ipOfServer.sin_family = af;
    ipOfServer.sin_port = htons(port);         // this is the port number of running server
    ipOfServer.sin_addr.s_addr = inet_addr("127.0.0.1");//inet_addr("10.192.13.11"); //inet_addr("127.0.0.1");
    
}

void connection::connectToServer()
{    
    if((ServerSocket = socket(af, type, protocol))< 0) // creating socket, returns file descriptor, 1st argument specifies the communication doamin, adress family, AF_INET belongs to IPv4 family, 2nd arg: type, SOCK_STREAM=TCP protocol, SOC_DGRAM=UDP, 3rd argument: protocol, 0= default for TCP
    {
        printf("Socket not created \n");
    }
    int n = connect(ServerSocket, (struct sockaddr *)&ipOfServer, sizeof(ipOfServer));
    while(n <0)
    {
        printf("Connection failed due to port and ip problems\n");
        //fprintf(stderr, "%s\n", xn_getlasterror(ServerSocket, (struct sockaddr *)&ipOfServer, sizeof(ipOfServer)));
        n = connect(ServerSocket, (struct sockaddr *)&ipOfServer, sizeof(ipOfServer));
    }
    
    std::cout<<"\nConnection successful"<<std::endl;
    
}

void connection::readDataWithLength(char ** msg) //allocated memory using new....please use dleete after use
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
        n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        }
        
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
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
        while (temp.size() >MSGLENGHTSIZE)
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
}

void connection::readData(char ** msg) //allocated memory using new....please use dleete after use
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
        n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        }
        
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
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

    }
}


void connection::readnData(int num) //allocated memory using new....please use dleete after use
{
    char dataReceived[DATASIZE];
    memset(dataReceived, '\0' ,DATASIZE);
    int n = 0;
        
    for(int i=0; i<num; ++i)
    {
        n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        }
        
        if(n <=0)
        {
            printf("\nError on read. Trying again..\n");
            n = (int)read(ServerSocket, dataReceived, sizeof(dataReceived)-1);
        }
        if(n <=0)
        {
            printf("\nreadnData: Error on read. Stopping now..\n");
            exit(0);
        }
        
        string temp;
        temp.append((char *) dataReceived);
        
        string * msg1 = new string(temp.c_str());
        readMsgPool.push_back(msg1);
           
     }
}


void connection::writeDataWithLength(char ** msg)
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
    int len = (int)lenPlusMsg.size(); ///lenght of msg lenght
    lenPlusMsg.clear();
    while((MSGLENGHTSIZE - len) >0)
    {
        lenPlusMsg.append("0");
        ++len;
    }
    lenPlusMsg.append(to_string(temp.size()));
    lenPlusMsg.append(temp.c_str());
    
    write(ServerSocket, lenPlusMsg.c_str(), lenPlusMsg.size());
    
#ifdef test_conn_deep
    cout<<writeCount++<<":"<<lenPlusMsg.c_str()<<endl;
#endif
    //cout<<writeCount++<<endl;

    //delete [] *msg; we dont know if mem is allocated with malloc or new
    //std::cout<<endl<<" Sending to Bob :"<<lenPlusMsg.c_str()<<std::endl;
}

void connection::writeData(char ** msg)
{


    //char dataSending[DATASIZE]; // Actually this is called packet in Network Communication, which contain data and send through.
    //time_t clock;
    //clock = time(NULL);
    //snprintf(dataSending, sizeof(dataSending), "%.24s\r\n", ctime(&clock)); // Printing successful message
    
    string Message;
    Message.append((char*)*msg) ;
    
    write(ServerSocket, Message.c_str(), Message.size());
    
    cout<<writeCount++<<":"<<Message.c_str()<<endl;

}

connection::~connection()
{
    close(ServerSocket);
}
