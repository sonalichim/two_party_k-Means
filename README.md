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
- Command to build: 
  1. In directory of party B or Bob: g++ --std=c++11 -g getPrime.cpp PaillierSystem.cpp server.cpp kMeansSubProtocols.cpp kMeansProtocol.cpp main.cpp -lgmp -o kMeansBob
  2. In directory of party A or Alice: g++ --std=c++11 -g getPrime.cpp PaillierSystem.cpp client.cpp kMeansSubProtocols.cpp kMeansProtocol.cpp main.cpp -lgmp -o kMeansAlice

## 📂 Repository Structure
```
two_party_k-Means/
│── k-means_A/                  # Party A (generates public-private key pair)
│   ├── k-means_Alice/          # Data directory for Party A
│   ├── client.cpp              # TCP connection establishment with Party B
│   ├── PaillierSystem.cpp      # Paillier cryptosystem implementation
│   ├── kMeansProtocol.cpp      # k-Means protocol for Party A
│   ├── kMeansSubProtocols.cpp  # Sub-protocols for k-Means
│   └── main.cpp                # Entry point for Party A
│
│── k-means_B/                  # Party B
│   ├── k-means_Bob/            # Data directory for Party B
│   ├── server.cpp              # TCP connection establishment with Party A
│   ├── PaillierSystem.cpp      # Paillier cryptosystem implementation
│   ├── kMeansProtocol.cpp      # k-Means protocol for Party B
│   ├── kMeansSubProtocols.cpp  # Sub-protocols for k-Means
│   └── main.cpp                # Entry point for Party B
```

## 👩‍💻 Author
**Sonali Chim**  
M.Tech (Research), Indian Institute of Science (IISc)  
Specialization:Cryptography  
[LinkedIn](https://www.linkedin.com/in/sonali-chim) | [GitHub](https://github.com/sonalichim)

---

## 📝 License
This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.


















