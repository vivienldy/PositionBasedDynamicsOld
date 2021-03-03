#include <glad/glad.h> 
#include <GLFW/glfw3.h>
#include<shader/shader_m.h>
#include <iostream>
#include<fstream>
#include<vector>
#include <string>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include"Storage.h"
#include "cuda_runtime.h"
# include <thrust/host_vector.h>

#include"PBD_Basic.cuh"

using namespace std;

int resY = 32;
int resX = 32;
float dampingRate = 0.9f;
float sizeX = 48.0f;
float sizeY = 64.0f;
float3 gravity = make_float3(0.0, -10.0, 0.0);
int startFrame = 1;
int endFrame = 200;
int substep = 4;
int iteration = 10;
HardwareType ht = CPU;
SolverType st = GAUSSSEIDEL;


/*
    CPU Run
    
    EdgeColoring
        Houdini Visual 

    SortColor SPMD
        一个点 拉特别远


*/
/*
int* color
int* sortedColor
int* sortedPrim

int main
{
sort(color, sortedColor, sortedPrim)
    //color: 11220001110222
    //sortedColor:000111222
    //sortedPrim: 
int2 workSets = Eval(SortedColor) 
    //workSets.x start
    //workSets.y num
for(auto workset: workSets)
{
    int start  = workSets.x
    int num = workSets.y
    kernel project<<<numBlock, numThread>>>(num, sortedPrim+start, prdPbuffer, restlength)
    {
        if(index > = num)
            return;
        int2 prim = sortedPrimId[index]
        int i0 = primList[prim.x]
        int i1 = primList[prim.x + 1]

        float3 p0 = prdPbuffer[i0]
        float3 p1 = prdPbuffer[i1]

    }
}
}

*/


void EdgeColoring()
{

}


int main()
{ 
    SolverPBD solver;
   // PBDObject pbdObj();
    //pbdObj.dampingRate = dampingRate;
    PBDObject pbdObj(dampingRate, gravity, resX, resY, sizeX, sizeY);
    solver.SetTarget(&pbdObj);
   // solver.SetTarget(PBDObject::Create("E:/top.cache","E:/cc.cache"));

    int FPS = 25;
    float dt = 1.0 / FPS / (float)substep;

    solver.Advect(dt, ht);
    solver.ProjectConstraint(ht, st, iteration);
    solver.Integration(dt, ht);
    //for (size_t i = startFrame; i < endFrame; i++)
    //{
    //    for (size_t s = 0; s < substep; s++)
    //    {
    //        solver.Advect(dt, ht);
    //        solver.ProjectConstraint(ht, st, iteration);
    //        solver.Integration(dt, ht);
    //    }
    //}      
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void  CreateSphereVertiecsIndices(vector<float>& sphereVertices, vector<int>& sphereIndices, vec3f gSphereOrigin, double sphereRadius);
vec2i  PrepareSphereBuffer();
vec2i  PrepareClothBuffer(vector<vec3f>& positionBuffer, vector<int> indices);
void CreateOpenGLIndices(vector<int>& openGLIndices, int rows, int cols);
void CreatePosition(vector<vec3f>& gPositionBuffer, vec2f cord, double length, double height, int rows, int cols);
double Distance(vec3f p1, vec3f p2);
void CreateIndices(vector<int>& indices, int rows, int cols);
void Init(
    vector<vec3f>& positionBuffer,
    vec2f cord,
    double length,
    double height,
    int rows,
    int cols,
    vector<int>& indices,
    vector<vec2i>& primList,
    vector<double>& restLength,
    vector<vec3f>& prdPBuffer,
    vector<vec3f>& velBuffer,
    vector<double>& massBuffer,
    vector<double>& stiffnessBuffer,
    vector<vec3f>& restPosBuffer);
void IntegrateVelocitys(
    std::vector<vec3f>& velBuffer,
    std::vector<vec3f>& positionBuffer,
    std::vector<vec3f>& prdPBuffer,
    float dt);
void Projection(
    int iterations,
    std::vector<vec2i>& primList,
    std::vector<int>& indices,
    std::vector<vec3f>& prdPBuffer,
    std::vector<double> massBuffer,
    std::vector<double> restLengthBuffer,
    std::vector<double> stiffnessBuffer,
    vec3f sphereOrigin,
    double sphereRadius,
    std::vector<vec3f> restPosBuffer,
    int cols,
    vec3f groundCenter);
vec3f GenerateMoveVectorSphere(
    vec3f sphereOrigin,
    double sphereRadius,
    vec3f p,
    int i);
void Advect(
    std::vector<vec3f>& velBuffer,
    vec3f acceleration,
    float dt,
    std::vector<vec3f>& prdPBuffer,
    double dampParam,
    std::vector<vec3f> positionBuffer);

// settings 窗口大小
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

vector<float> gSphereVertices;
vector<int> gSphereIndices;


double gSphereRadius = 5.0;
vec3f gGroundCenter(25.0, -40.0, 35.0);
vec2f gCord(0.0, 0.0);

vec3f gSphereOrigin(10, -10, 40);

std::vector<int> gIndices;
std::vector<int> gOpenGLIndices;
std::vector<double>gMassBuffer;
std::vector<double> gRestLengthBuffer;
std::vector<double> gStiffnessBuffer;
std::vector<vec2i> gPrimList;
std::vector<vec3f> gPositionBuffer;
std::vector<vec3f> gVelBuffer;
std::vector<vec3f> gPrdPBuffer;
std::vector<vec3f> gRestPostionBuffer;


// camera
glm::vec3 cameraPos = glm::vec3(-35.0f, -15.0f, 33.0f);
glm::vec3 cameraFront = glm::vec3(1.0f, 0.0f, 0.0f);
glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);

//int main()
//{
//
//
//    Init(gPositionBuffer,
//        gCord,
//        gLength,
//        gHeight,
//        gRows,
//        gCols,
//        gIndices,
//        gPrimList,
//        gRestLengthBuffer,
//        gPrdPBuffer,
//        gVelBuffer,
//        gMassBuffer,
//        gStiffnessBuffer,
//        gRestPostionBuffer);
//
//    vector<vec3f>* positionBufferd;
//    vector<vec3f>* velBufferd;
//    vector<vec3f>* prdPBufferd;
//    int size = gPositionBuffer.size() * sizeof(gPositionBuffer[0]);
//
//
//    // GPU端分配内存
//    cudaMalloc((void**)&positionBufferd, size);
//    cudaMalloc((void**)&velBufferd, size);
//    cudaMalloc((void**)&prdPBufferd, size);
//
//    // CPU的数据拷贝到GPU端
//    cudaMemcpy(velBufferd, &gVelBuffer[0], size, cudaMemcpyHostToDevice);
//    cudaMemcpy(positionBufferd, &gPositionBuffer[0], size, cudaMemcpyHostToDevice);
//    cudaMemcpy(prdPBufferd, &gPrdPBuffer[0], size, cudaMemcpyHostToDevice);
//
//    dim3 DimGrid(4);  //4个block
//    dim3 DimBlock(512);
//
//    //loop
//    int FPS = 25;
//    float dt = 1.0 / FPS / (float)gSubstep;
//
//    // 执行kernel
//    Advect << <DimGrid, DimBlock >> > (velBufferd, positionBufferd, prdPBufferd, gAcceleration, dt, gDampParam);
//
//    // 将在GPU端计算好的结果拷贝回CPU端
//    cudaMemcpy(&gVelBuffer[0], velBufferd, size, cudaMemcpyDeviceToHost);
//    cudaMemcpy(&gPositionBuffer[0], positionBufferd, size, cudaMemcpyDeviceToHost);
//    cudaMemcpy(&gPrdPBuffer[0], prdPBufferd, size, cudaMemcpyDeviceToHost);



    //CreateOpenGLIndices(gOpenGLIndices, gRows, gCols);

    //glfwInit();
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    ////创建窗口
    //GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    //if (window == NULL)
    //{
    //    //std::cout << "Failed to create GLFW window" << std::endl;
    //    glfwTerminate();
    //    return -1;
    //}
    //glfwMakeContextCurrent(window);
    //glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    //// glad: load all OpenGL function pointers
    //if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    //{
    //    //std::cout << "Failed to initialize GLAD" << std::endl;
    //    return -1;
    //}

    ////load shader
    //Shader ourShader("D:\\OpenGL\\shaders\\vertexshader.vs", "D:\\OpenGL\\shaders\\fragmentshader.fs");

    ////prepare sphere buffers
    //CreateSphereVertiecsIndices(gSphereVertices, gSphereIndices, gSphereOrigin, gSphereRadius);
    //unsigned int VAO = PrepareSphereBuffer().x;
    //unsigned int VBO = PrepareSphereBuffer().y;


    ////loop
    //int FPS = 25;
    //float dt = 1.0 / FPS / (float)gSubstep;
    //int frameIndex = 1;
    //for (size_t i = gStartFrame; i < gEndFrame; i++)
    //{
    //    //std::cout << "ENTRY FRAME:" << i << std::endl;
    //    for (size_t s = 0; s < gSubstep; s++)
    //    {
    //        //std::cout << "      ENTRY Step Sub : " << s << std::endl;
    //        Advect(gVelBuffer,
    //            gAcceleration,
    //            dt,
    //            gPrdPBuffer,
    //            gDampParam,
    //            gPositionBuffer);

    //        Projection(
    //            gIteration,
    //            gPrimList,
    //            gIndices,
    //            gPrdPBuffer,
    //            gMassBuffer,
    //            gRestLengthBuffer,
    //            gStiffnessBuffer,
    //            gSphereOrigin,
    //            gSphereRadius,
    //            gRestPostionBuffer,
    //            gCols,
    //            gGroundCenter);

    //        IntegrateVelocitys(
    //            gVelBuffer,
    //            gPositionBuffer,
    //            gPrdPBuffer,
    //            dt);

    //        // std::cout << "  " << std::endl;
    //        frameIndex++;
    //    }
    //    // std::cout << "----END FRAME-----" << std::endl;

    //     //prepare cloth buffer
    //    unsigned int VAOTri = PrepareClothBuffer(gPositionBuffer, gOpenGLIndices).x;
    //    unsigned int VBOTri = PrepareClothBuffer(gPositionBuffer, gOpenGLIndices).y;

    //    //input
    //    processInput(window);
    //    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
    //    glClear(GL_COLOR_BUFFER_BIT);

    //    // activate shader
    //    ourShader.use();

    //    //prare projection and view matrix
    //    glm::mat4 projection = glm::perspective(glm::radians(80.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    //    ourShader.setMat4("projection", projection);

    //    glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    //    ourShader.setMat4("view", view);

    //    // rendering commands
    //   /* glEnable(GL_CULL_FACE);
    //    glCullFace(GL_BACK);*/
    //    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    //    //--------------------------------------------draw cloth--------------------------------------------
    //    glm::mat4 transform = glm::mat4(1.0f);
    //    transform = glm::translate(transform, glm::vec3(0.5f, -0.5f, 0.0f));
    //    ourShader.setMat4("model", transform);

    //    glBindVertexArray(VAOTri);
    //    glDrawElements(GL_TRIANGLES, gOpenGLIndices.size(), GL_UNSIGNED_INT, 0);

    //    //--------------------------------------------draw sphere--------------------------------------------
    //    glm::mat4 transformSphere = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    //    transformSphere = glm::translate(transformSphere, glm::vec3(gSphereOrigin.x, gSphereOrigin.y, gSphereOrigin.z));
    //    transformSphere = glm::scale(transformSphere, glm::vec3(gSphereRadius, gSphereRadius, gSphereRadius));
    //    ourShader.setMat4("model", transformSphere);

    //    glBindVertexArray(VAO);
    //    glDrawElements(GL_TRIANGLES, 20 * 20 * 6, GL_UNSIGNED_INT, 0);

    //    // check and call events and swap the buffers
    //    glfwSwapBuffers(window);
    //    glfwPollEvents();
    //}
    //glfwTerminate();
//    return 0;
//}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void  CreateSphereVertiecsIndices(vector<float>& sphereVertices, vector<int>& sphereIndices, vec3f gSphereOrigin, double sphereRadius)
{
    int Y_SEGMENTS = 20;
    int X_SEGMENTS = 20;
    const float PI = 3.14159265358979323846f;
    for (int y = 0; y <= Y_SEGMENTS; y++)
    {
        for (int x = 0; x <= X_SEGMENTS; x++)
        {
            float xSegment = (float)x / (float)X_SEGMENTS;
            float ySegment = (float)y / (float)Y_SEGMENTS;
            float xPos = std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI);
            float zPos = std::cos(ySegment * PI);
            float yPos = std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI);
            sphereVertices.push_back(xPos);
            sphereVertices.push_back(yPos);
            sphereVertices.push_back(zPos);
        }
    }
    for (int i = 0; i < Y_SEGMENTS; i++)
    {
        for (int j = 0; j < X_SEGMENTS; j++)
        {
            sphereIndices.push_back(i * (X_SEGMENTS + 1) + j);
            sphereIndices.push_back((i + 1) * (X_SEGMENTS + 1) + j);
            sphereIndices.push_back((i + 1) * (X_SEGMENTS + 1) + j + 1);
            sphereIndices.push_back(i * (X_SEGMENTS + 1) + j);
            sphereIndices.push_back((i + 1) * (X_SEGMENTS + 1) + j + 1);
            sphereIndices.push_back(i * (X_SEGMENTS + 1) + j + 1);
        }
    }
}

vec2i PrepareSphereBuffer()
{
    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    //生成并绑定球体的VAO和VBO
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    //将顶点数据绑定至当前默认的缓冲中
    glBufferData(GL_ARRAY_BUFFER, gSphereVertices.size() * sizeof(float), &gSphereVertices[0], GL_STATIC_DRAW);

    GLuint element_buffer_object;//EBO
    glGenBuffers(1, &element_buffer_object);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, element_buffer_object);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, gSphereIndices.size() * sizeof(int), &gSphereIndices[0], GL_STATIC_DRAW);

    //设置顶点属性指针
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    //解绑VAO和VBO
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    vec2i bufferList = vec2i(VAO, VBO);
    return bufferList;
}

vec2i  PrepareClothBuffer(vector<vec3f>& positionBuffer, vector<int> indices)
{
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(positionBuffer[0]) * positionBuffer.size(), &positionBuffer[0], GL_DYNAMIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0]) * indices.size(), &indices[0], GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(positionBuffer[0]), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    vec2i bufferList = vec2i(VAO, VBO);
    return bufferList;
}

vec2i  PrepareGroundBuffer(vector<vec3f>& positionBuffer, vector<int> indices)
{
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(positionBuffer[0]) * positionBuffer.size(), &positionBuffer[0], GL_DYNAMIC_DRAW);
    //glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    //glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0]) * indices.size(), &indices[0], GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(positionBuffer[0]), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    vec2i bufferList = vec2i(VAO, VBO);
    return bufferList;
}

double Distance(vec3f p1, vec3f p2)
{
    return pow(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2) + pow((p1.z - p2.z), 2), 0.5);
}

void CreatePosition(vector<vec3f>& gPositionBuffer, vec2f cord, double length, double height, int rows, int cols)
{
    double lengthInterval = length / (cols - 1);
    double heightInterval = height / (rows - 1);
    int num = rows * cols;
    int index = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            vec3f p;
            p.x = cord.x + j * lengthInterval;
            p.y = 0;
            p.z = cord.y + i * heightInterval;
            gPositionBuffer.push_back(p);
            index++;
        }
    }
}

void CreateOpenGLIndices(vector<int>& openGLIndices, int rows, int cols)
{
    int num = rows * cols;
    for (int i = 0; i < num - cols; i++)
    {
        if (i % cols == cols - 1)
            continue;
        openGLIndices.push_back(i);
        openGLIndices.push_back(i + cols + 1);
        openGLIndices.push_back(i + cols);
        openGLIndices.push_back(i);
        openGLIndices.push_back(i + 1);
        openGLIndices.push_back(i + cols + 1);

    }
}

void CreateIndices(vector<int>& indices, int rows, int cols)
{
    int num = rows * cols;
    for (int i = 0; i < num - cols; i++)
    {
        if (i % cols == 0)
        {
            indices.push_back(i);
            indices.push_back(i + 1);
            indices.push_back(i);
            indices.push_back(i + cols);
        }
        else if (i % cols == cols - 1)
        {
            indices.push_back(i);
            indices.push_back(i + cols);
            indices.push_back(i);
            indices.push_back(i + cols - 1);
        }
        else
        {
            indices.push_back(i);
            indices.push_back(i + 1);
            indices.push_back(i);
            indices.push_back(i + cols);
            indices.push_back(i);
            indices.push_back(i + cols - 1);
        }
    }
}

void CreatePrimList(vector<vec2i>& primList, vector<int>& indices, int count)
{
    for (int i = 0; i < indices.size() / count; i++)
    {
        vec2i p;
        p.x = i * count;
        p.y = count;;
        primList.push_back(p);
    }
}

void CreateRestLength(vector<double>& restLength, vector<vec2i>& primList, vector<int>& indices, vector<vec3f>& positionBuffer)
{
    for (int i = 0; i < primList.size(); i++)
    {
        int i0 = indices[primList[i].x];
        int i1 = indices[primList[i].x + 1];
        double d = Distance(positionBuffer[i0], positionBuffer[i1]);
        restLength.push_back(d);
    }
}

void Init(
    vector<vec3f>& positionBuffer,
    vec2f cord,
    double length,
    double height,
    int rows,
    int cols,
    vector<int>& indices,
    vector<vec2i>& primList,
    vector<double>& restLength,
    vector<vec3f>& prdPBuffer,
    vector<vec3f>& velBuffer,
    vector<double>& massBuffer,
    vector<double>& stiffnessBuffer,
    vector<vec3f>& restPosBuffer)
{
    CreatePosition(positionBuffer, cord, length, height, rows, cols);
    CreateIndices(indices, rows, cols);
    CreatePrimList(primList, indices, 2);
    CreateRestLength(restLength, primList, indices, positionBuffer);
    for (int i = 0; i < rows * cols; i++)
    {
        prdPBuffer.push_back(positionBuffer[i]);
        restPosBuffer.push_back(positionBuffer[i]);
        velBuffer.push_back(vec3f(0.0, 0.0, 0.0));

        massBuffer.push_back(1.0);
    }
    for (int j = 0; j < primList.size(); j++)
    {
        stiffnessBuffer.push_back(1.0);
    }
    //std::cout << "Init" << std::endl;
}

void DampVelocities(
    std::vector<vec3f>& velBuffer,
    float dt,
    double dampVelParam)
{
    for (int i = 0; i < velBuffer.size(); i++)
    {
        velBuffer[i] *= pow(dampVelParam, dt);
    }
}

void Advect(
    std::vector<vec3f>& velBuffer,
    vec3f acceleration,
    float dt,
    std::vector<vec3f>& prdPBuffer,
    double dampParam,
    std::vector<vec3f> positionBuffer)
{
    for (int i = 0; i < velBuffer.size(); i++)
    {
        velBuffer[i] += acceleration * dt;
    }
    // std::cout << "       V+=A*dt:" << dt << std::endl;

    DampVelocities(velBuffer, dt, dampParam);
    // std::cout << "       DampVelocities" << std::endl;

    for (int j = 0; j < prdPBuffer.size(); j++)
    {
        prdPBuffer[j] = positionBuffer[j] + velBuffer[j] * dt;
    }
    // std::cout << "       PrdP+=V*dt:" << dt << std::endl;
}

bool ColliderSphere(vec3f pointPos, vec3f sphereOrigin, double r, int i)
{
    double d = Distance(pointPos, sphereOrigin);

    if (d - r > 0.001)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool CollideGround(vec3f pointPos, vec3f groundCenter)
{
    if (pointPos.y - groundCenter.y < 0.001)
    {
        return true;
    }
    else
    {
        return false;
    }
}

vec3f GenerateMoveVectorSphere(
    vec3f sphereOrigin,
    double sphereRadius,
    vec3f p,
    int i)
{
    double moveDistance = sphereRadius - Distance(sphereOrigin, p);
    vec3f moveDirection = (p - sphereOrigin) / Distance(sphereOrigin, p);
    vec3f moveLength = moveDirection * moveDistance;
    return moveLength;
}

void Projection(
    int iterations,
    std::vector<vec2i>& primList,
    std::vector<int>& indices,
    std::vector<vec3f>& prdPBuffer,
    std::vector<double> massBuffer,
    std::vector<double> restLengthBuffer,
    std::vector<double> stiffnessBuffer,
    vec3f sphereOrigin,
    double sphereRadius,
    std::vector<vec3f> restPosBuffer,
    int cols,
    vec3f groundCenter)
{
    // std::cout << "       projection: " << iterations << std::endl;
     //约束
     //预测位置
    for (size_t ii = 0; ii < iterations; ii++)
    {
        for (size_t i = 0; i < primList.size(); i++)
        {
            if (primList[i].y != 2)
                continue;
            int i0 = indices[primList[i].x];
            int i1 = indices[primList[i].x + 1];
            vec3f dp1;
            vec3f dp2;
            double d = Distance(prdPBuffer[i0], prdPBuffer[i1]);
            vec3f v;
            v = prdPBuffer[i0] - prdPBuffer[i1];
            dp1.x = -massBuffer[i0] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.x / d;
            dp1.y = -massBuffer[i0] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.y / d;
            dp1.z = -massBuffer[i0] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.z / d;
            dp2.x = massBuffer[i1] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.x / d;
            dp2.y = massBuffer[i1] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.y / d;
            dp2.z = massBuffer[i1] / (massBuffer[i0] + massBuffer[i1]) * (d - restLengthBuffer[i]) * v.z / d;
            double k = 1 - pow(1 - stiffnessBuffer[i], 1.0 / (ii + 1));
            dp1 *= k;
            dp2 *= k;
            prdPBuffer[i0] += dp1;
            prdPBuffer[i1] += dp2;
        }

        for (size_t j = 0; j < prdPBuffer.size(); j++)
        {
            //attach points
            if (j == 0 || j == cols - 1)
            {
                prdPBuffer[j] = restPosBuffer[j];
            }

            //point collide with sphere
            bool isCollideSphere = ColliderSphere(prdPBuffer[j], sphereOrigin, sphereRadius, j);
            if (isCollideSphere) //move the point to the point which intersect with sphere
            {
                vec3f moveVector = GenerateMoveVectorSphere(sphereOrigin, sphereRadius, prdPBuffer[j], j);
                prdPBuffer[j] += moveVector;
            }
            //point collide with ground
            bool isCollideGoround = CollideGround(prdPBuffer[j], groundCenter);
            if (isCollideGoround)
            {
                prdPBuffer[j].y = groundCenter.y;
            }
        }
    }
}

void IntegrateVelocitys(
    std::vector<vec3f>& velBuffer,
    std::vector<vec3f>& positionBuffer,
    std::vector<vec3f>& prdPBuffer,
    float dt)
{
    for (size_t i = 0; i < positionBuffer.size(); i++)
    {
        velBuffer[i] = (prdPBuffer[i] - positionBuffer[i]) / dt;
        positionBuffer[i] = prdPBuffer[i];
    }
    // std::cout << "       IntegrateVelocitys" << std::endl;
}

