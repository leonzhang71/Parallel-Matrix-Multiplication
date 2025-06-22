# 🔢 Parallel Matrix Multiplication

A high-performance implementation of matrix multiplication using **Strassen’s algorithm** and **OpenMP-based parallelization**. Developed as part of a Parallel Computing course to explore recursive algorithms and multithreading performance gains.

---

## 🚀 Objective

Optimize matrix multiplication for large-scale square matrices by reducing computational complexity and leveraging parallel threads for improved runtime performance.

---

## 🧠 Algorithm

- Implemented a recursive version of **Strassen’s matrix multiplication** algorithm.
- Applied **OpenMP** directives to parallelize recursive steps.
- Used **cache-aware block allocation** to enhance memory performance.

---

## 📊 Results

- Achieved up to **3x speedup** over naive matrix multiplication for 1024x1024 matrices.
- Benchmarked against standard C++ implementations using high-resolution timers.

---

## 🔧 Tech Stack

- **Language:** C++  
- **Parallelization:** OpenMP  
- **Tools:** g++, Unix `time` command, custom timers

---

## 🧪 How to Run

1. Clone this repo  
2. Go throught Major-Project.pdf, it shows how to run steps by steps

---

---

## 📌 Acknowledgments

- Developed as part of coursework in Parallel Computing
- Special thanks to course instructors for algorithmic guidance
