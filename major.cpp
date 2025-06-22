#include <omp.h>
#include <bits/stdc++.h>
#include <stdio.h>      
#include <stdlib.h>  
#include <math.h>
#include <time.h>   

using namespace std;

#define MAX_THREADS 65536
#define MAX_MATRIX_SIZE INT_MAX

int matrix_size, num_threads, base;

void print(int size, unsigned long long** matrix){
    for (int row = 0; row < size; row++){
        for (int col = 0; col < size; col++)
            std::cout << matrix[row][col] << ' ';
        std::cout << '\n';
    }
    std::cout << '\n';
}

unsigned long long** createSquareMatrix(int dimension) {
    unsigned long long* values = (unsigned long long*)malloc(dimension * dimension * sizeof(unsigned long long));
    unsigned long long** matrixRows = (unsigned long long**)malloc(dimension * sizeof(unsigned long long*));

    for (int row = 0; row < dimension; row++)
        matrixRows[row] = &values[dimension * row];

    return matrixRows;
}

void populateMatrix(int size, unsigned long long**& matrix){
    srand48(0);
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            matrix[row][col] = lrand48() % INT_MAX;
        }
    }
}


unsigned long long** parallel_omp(int dimension, unsigned long long** matrixA, unsigned long long** matrixB) {
    unsigned long long** matrixProduct = createSquareMatrix(dimension);

    #pragma omp parallel for collapse(2)
    for (int row = 0; row < dimension; row++) {
        for (int col = 0; col < dimension; col++) {
            matrixProduct[row][col] = 0;
            for (int inner = 0; inner < dimension; inner++)
                matrixProduct[row][col] += matrixA[row][inner] * matrixB[inner][col];
        }
    }
    return matrixProduct;
}

unsigned long long** extractSubMatrix(int size, unsigned long long** originalMatrix, int startRow, int startCol) {
    int subSize = size / 2;
    unsigned long long** subMatrix = createSquareMatrix(subSize);
    for (int row = 0; row < subSize; row++) {
        for (int col = 0; col < subSize; col++) {
            subMatrix[row][col] = originalMatrix[startRow + row][startCol + col];
        }
    }
    return subMatrix;
}


unsigned long long** addMatrices(int n, unsigned long long** matrixA, unsigned long long** matrixB, bool op){
    unsigned long long** result = createSquareMatrix(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            result[i][j] = op ? matrixA[i][j] + matrixB[i][j] : matrixA[i][j] - matrixB[i][j];
        }
    }
    return result;
}

unsigned long long** mergeSubMatrix(int subSize, unsigned long long** topLeft, unsigned long long** topRight, unsigned long long** bottomLeft, unsigned long long** bottomRight) {
    int fullSize = 2 * subSize;
    unsigned long long** mergedMatrix = createSquareMatrix(fullSize);

    for (int row = 0; row < fullSize; row++) {
        for (int col = 0; col < fullSize; col++) {
            if (row < subSize) {
                mergedMatrix[row][col] = (col < subSize) ? topLeft[row][col] : topRight[row][col - subSize];
            } else {
                mergedMatrix[row][col] = (col < subSize) ? bottomLeft[row - subSize][col] : bottomRight[row - subSize][col - subSize];
            }
        }
    }
    return mergedMatrix;
}


unsigned long long** strassen(int n, unsigned long long** matrixA, unsigned long long** matrixB){
    if (n <= base) return parallel_omp(n, matrixA, matrixB);

    int m = n/2;
    bool add = true, sub = false;

    unsigned long long** a11 = extractSubMatrix(n, matrixA, 0, 0);
    unsigned long long** a12 = extractSubMatrix(n, matrixA, 0, m);
    unsigned long long** a21 = extractSubMatrix(n, matrixA, m, 0);
    unsigned long long** a22 = extractSubMatrix(n, matrixA, m, m);
    unsigned long long** b11 = extractSubMatrix(n, matrixB, 0, 0);
    unsigned long long** b12 = extractSubMatrix(n, matrixB, 0, m);
    unsigned long long** b21 = extractSubMatrix(n, matrixB, m, 0);
    unsigned long long** b22 = extractSubMatrix(n, matrixB, m, m);

    unsigned long long** m1, **m2, **m3, **m4, **m5, **m6, **m7;

    #pragma omp parallel sections
    {
        //Compute M1: M1 = strassen(A11 + A22, B11 + B22)
        #pragma omp section
        {
            unsigned long long** temp1 = addMatrices(m, a11, a22, true);// true for addition
            unsigned long long** temp2 = addMatrices(m, b11, b22, true);// true for addition
            m1 = strassen(m, temp1, temp2);
            free(temp1);
            free(temp2);
        }
        //Compute M2: M2 = strassen(A21 + A22, B11)
        #pragma omp section
        {
            unsigned long long** temp = addMatrices(m, a21, a22, true);// true for addition
            m2 = strassen(m, temp, b11);
            free(temp);
        }
        #pragma omp section
        {
        // Compute M3: M3 = strassen(a11, b12 - b22)
        unsigned long long** temp = addMatrices(m, b12, b22, false); // false for subtraction
        m3 = strassen(m, a11, temp);
        free(temp);
        }
    #pragma omp section
        {
        // Compute M4: M4 = strassen(a22, b21 - b11)
        unsigned long long** temp = addMatrices(m, b21, b11, false); // false for subtraction
        m4 = strassen(m, a22, temp);
        free(temp);
        }
    #pragma omp section
        {
        // Compute M5: M5 = strassen(a11 + a12, b22)
        unsigned long long** temp = addMatrices(m, a11, a12, true); // true for addition
        m5 = strassen(m, temp, b22);
        free(temp);
        }
    #pragma omp section
        {
        // Compute M6: M6 = strassen(a21 - a11, b11 + b12)
        unsigned long long** temp1 = addMatrices(m, a21, a11, false); // false for subtraction
        unsigned long long** temp2 = addMatrices(m, b11, b12, true);  // true for addition
        m6 = strassen(m, temp1, temp2);
        free(temp1);
        free(temp2);
        }
    #pragma omp section
        {
        // Compute M7: M7 = strassen(a12 - a22, b21 + b22)
        unsigned long long** temp1 = addMatrices(m, a12, a22, false); // false for subtraction
        unsigned long long** temp2 = addMatrices(m, b21, b22, true);  // true for addition
        m7 = strassen(m, temp1, temp2);
        free(temp1);
        free(temp2);
        }
    }

    unsigned long long** c11 = addMatrices(m, m1, m4, add);
    unsigned long long** c12 = addMatrices(m, m3, m5, add);
    unsigned long long** c21 = addMatrices(m, m2, m4, add);
    unsigned long long** c22 = addMatrices(m, m1, addMatrices(m, m3, m6, sub), sub);

    unsigned long long** prod = mergeSubMatrix(m, c11, c12, c21, c22);

    // Free allocated memory
    free(c11); free(c12); free(c21); free(c22);

    return prod;
}

int main(int argc, char *argv[]){
    int k, q, n;
    struct timespec start, stop;
    double total_time;
    
    if (argc != 4){
        printf("Please Enter the Size of Matrix and Number of Threads!\n");
        exit(0);
    }
    else{
        base = atoi(argv[argc-3]);
        k = atoi(argv[argc-2]);
        n = 1 << k;
        q = atoi(argv[argc-1]);
        num_threads = 1 << q;
    } 

    matrix_size = n;

    unsigned long long** a = createSquareMatrix(n);
    unsigned long long** b = createSquareMatrix(n);
    unsigned long long** prod;

    populateMatrix(n, a);
    populateMatrix(n, b);

    omp_set_dynamic(1);
    omp_set_num_threads(num_threads);

    clock_gettime(CLOCK_REALTIME, &start);
    prod = strassen(n, a, b);
    clock_gettime(CLOCK_REALTIME, &stop);

    total_time = (stop.tv_sec-start.tv_sec) + 0.000000001*(stop.tv_nsec-start.tv_nsec);

    printf("k': %d, Matrix Size: %d x %d, Threads: %d, Time: %8.5f sec\n", base, n, n, num_threads, total_time);

    // Free memory
    free(a);
    free(b);
    free(prod);

    return 0;
}

