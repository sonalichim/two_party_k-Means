# M.Tech Cryptography Project

This repository contains the implementation of my **M.Tech Research Project** in cryptography **Efficient and secure two-party k-means clustering protocol**
The project explores the design and implementation of cryptographic algorithms and demonstrates how they can be executed via the command line.


## 🔑 Features
- Implementation of **Paillier Cryptosystem**
- Command-line interface for execution
- Modular code for experimenting with cryptographic primitives

## ⚙️ Requirements
- C++ compiler (e.g., `g++`)  
- GMP Library 
- Works on Linux / macOS command prompt

## 🚀 How to Build & Run
- 
- Command to build: 
  1. In directory of party B or Bob: g++ --std=c++11 -g getPrime.cpp PaillierSystem.cpp server.cpp kMeansSubProtocols.cpp kMeansProtocol.cpp main.cpp -lgmp -o kMeansBob
  2. 1. In directory of party A or Alice: g++ --std=c++11 -g getPrime.cpp PaillierSystem.cpp client.cpp kMeansSubProtocols.cpp kMeansProtocol.cpp main.cpp -lgmp -o kMeansAlice

## Repository Structure 
two_party_k-Means/
│── k-means_A/                 # Party A (generates public-private key pair)
    |── k-means_Alice/         # This directory contains data of party A
    |── client.cpp             # TCP connection establishment with party B
    |── PaillierSystem.cpp     # Paillier function
    |── kMeansProtocol.cpp     # k-Means protocol for party A
    |── kMeansSubProtocols.cpp # k-Means sub-protocol
    |── main.cpp               # main entry point
    
│── k-means_B/                 # Party B
    |── k-means_Bob/           # This directory contains data of party B
    |── server.cpp             # TCP connection establishment with party A
    |── PaillierSystem.cpp     # Paillier function
    |── kMeansProtocol.cpp     # k-Means protocol for party B
    |── kMeansSubProtocols.cpp # k-Means sub-protocol
    |── main.cpp               # main entry point









