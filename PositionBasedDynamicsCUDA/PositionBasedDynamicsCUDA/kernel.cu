//#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
//
//#include <stdio.h>
//#include <glad/glad.h> 
//#include <GLFW/glfw3.h>
//#include<shader/shader_m.h>
//#include <iostream>
//#include<fstream>
//#include<vector>
//#include <string>
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
//#include"Storage.cuh"
//# include <thrust/device_vector.h>

//using namespace std;
//__global__ void Plus(float A[], float B[], float C[], int n)
//{
//    int i = blockDim.x * blockIdx.x + threadIdx.x;
//    C[i] = A[i] + B[i];
//}

//__global__ void Advect(vector<vec3f> velBuffer,
//    vector<vec3f> positionBuffer,
//    vector<vec3f> prdPBuffer,
//    vec3f acceleration,
//    float dt,
//    double dampParam)
//{
//    int positionIndex = threadIdx.x + blockIdx.x * 512.0f;
//    velBuffer[positionIndex] = velBuffer[positionIndex] + acceleration * dt;
//    velBuffer[positionIndex] = velBuffer[positionIndex] * powf(dampParam, dt);
//
//    prdPBuffer[positionIndex] = positionBuffer[positionIndex] + velBuffer[positionIndex] * dt;
//}

//__global__ void Advect(thrust::device_vector <float3>  velBuffer,
//    thrust::device_vector <float3 > positionBuffer,
//    thrust::device_vector <float3 > prdPBuffer,
//    vec3f acceleration,
//    float dt,
//    double dampParam)
//{
//    int positionIndex = threadIdx.x + blockIdx.x * 512.0f;
//    velBuffer[positionIndex] = velBuffer[positionIndex] + acceleration * dt;
//    velBuffer[positionIndex] = velBuffer[positionIndex] * powf(dampParam, dt);
//
//    prdPBuffer[positionIndex] = positionBuffer[positionIndex] + velBuffer[positionIndex] * dt;
//}