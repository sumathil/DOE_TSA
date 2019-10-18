#pragma once

#include <stdio.h>
#include <iostream>
#include <complex.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <sys/time.h>
#include <cuda_runtime.h>

using namespace std;

#define Efd 2.4296
#define Tm 0.9
#define x_d 1.93
#define x_dp 0.23
#define x_q 1.77
#define x_qp 0.50
#define T_dp 5.2
#define T_qp 0.81
#define H 3.74
#define Pt 0.9
#define Qt 0.36
#define Vt 1.0928
#define Eb 1.0
#define pi 3.14159265358979323846

#define CHECK(call) \
{                                                                        \
        const cudaError_t error = call;                                       \
        if (error != cudaSuccess)                                             \
        {                                                                     \
                printf("Error: %s:%d, ", __FILE__, __LINE__);                      \
                printf("code:%d, reason: %s\n", error, cudaGetErrorString(error)); \
                exit(1);                                                           \
        }                                                                     \
}

