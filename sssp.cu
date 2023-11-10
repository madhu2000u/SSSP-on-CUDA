#include <iostream>
#include "sssp_common.h"


__device__ bool d_done = false;

float* computeMatrixMult(matElement*);

void setUpArrays(float *d_c, int *vertex, int *edges, bool *threadMask, float* cost, float* intermediateCost, matElement* minElement)
{   
    int edgeIndex = 0;
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {   
            threadMask[i * MATRIX_SIZE + j] = false;
            cost[i * MATRIX_SIZE + j] = __FLT_MAX__;

            vertex[i * MATRIX_SIZE + j] = edgeIndex;
            if((j + 1) < MATRIX_SIZE)
            {
                intermediateCost[edgeIndex] = __FLT_MAX__;            //intemediateCost array size is same as that of edges so do it in this step itself so that you can fill all the spots in the array just like edges array instead of making another loop for it
                edges[edgeIndex++] = i * MATRIX_SIZE + (j + 1);
            }
            if((i + 1) < MATRIX_SIZE)
            {
                intermediateCost[edgeIndex] = __FLT_MAX__;
                edges[edgeIndex++] = (i + 1) * MATRIX_SIZE + j;
            }
            if((j - 1) >= 0)
            {
                intermediateCost[edgeIndex] = __FLT_MAX__;
                edges[edgeIndex++] = i * MATRIX_SIZE + (j - 1);
            }
            if((i - 1) >= 0)
            {
                intermediateCost[edgeIndex] = __FLT_MAX__;
                edges[edgeIndex++] = (i - 1) * MATRIX_SIZE + j;
            }
        }

        threadMask[minElement[0].row * MATRIX_SIZE + minElement[0].col] = true;             //Make the thread of source vertex executable initially since that is the starting point
        cost[minElement[0].row * MATRIX_SIZE + minElement[0].col] = 0.0f;                   //Cost from source to source is 0
        intermediateCost[minElement[0].row * MATRIX_SIZE + minElement[0].col] = 0.0f;       
        
    }
    
}

// void printMatrixInRange(int rowStart, int rowEnd, float *d_c)
// {
//     for (int i = rowStart; i < MATRIX_SIZE; i++)
//     {
//         for (int j = rowEnd; j < MATRIX_SIZE; j++)
//         {
//             std::cout<<d_c[rowStart * MATRIX_SIZE + rowEnd]<<"       ";
//         }
//         std::cout<<"\n";
        
//     }
    
// }

void printNeighbors(int index, float *d_c, int *vertex, int* edges)
{
    for (int i = vertex[index]; i < vertex[index + 1]; i++)
    {
        std::cout<<i<<std::endl;
        std::cout<<d_c[edges[i]]<<std::endl;
        std::cout<<"\n";
    }
    
}

__global__ void computeIntermediates()
{
    
}

int main()
{   
    int *vertex, *edges;
    float *d_c, *cost, *intermediateCost;
    bool *threadMask;

    bool h_done = false;

    matElement minElement[2];
    size_t size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    //Note: Diagnonal neighbours are not considered
    int numEdges = (4 * 2                                       /*Each corner values in matrix has 2 neighbours*/ 
                    + ((MATRIX_SIZE - 2) * 3) * 4               /*Each element of the 4 boundary sides excluding the 2 corner elements for each boundary side has 3 neighbours*/ 
                    + (MATRIX_SIZE - 2) * 4 * (MATRIX_SIZE - 2) /*Each element not on the boundary has 4 neighbours*/);
    
    //Compute and get the pointer to the result matrix of the matrix muliplications
    d_c = computeMatrixMult(minElement);


    vertex = (int*)malloc(((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int));                            // + 1 because we need a location at the end of the vertex that stores the ending index of the edge
    threadMask = new bool[MATRIX_SIZE * MATRIX_SIZE];
    edges = (int*)malloc( numEdges * sizeof(int));

    cost = (float*)malloc(size);
    intermediateCost = (float*)malloc(numEdges * sizeof(float));

    // float test_data[16] = {1, 5, 3, 2, 9, 4, 10, 7, 6, 8, 12, 11, 14, 15, 16, 17};

    setUpArrays(d_c, vertex, edges, threadMask, cost, intermediateCost, minElement);
    vertex[MATRIX_SIZE * MATRIX_SIZE] = numEdges;                                                     //last value in vertex is total numEdges so that we can use the starting and ending index when getting the neighbors


    //Setup CUDA device memories for the data
    CHECK(cudaMallocHost(&vertex, ((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int)));
    CHECK(cudaMallocHost(&edges, numEdges * sizeof(int)));
    CHECK(cudaMallocHost(&threadMask, MATRIX_SIZE * MATRIX_SIZE * sizeof(bool)));
    CHECK(cudaMallocHost(&cost, size));
    CHECK(cudaMallocHost(&intermediateCost, numEdges * sizeof(float)));

    //Start computing SSSP
    while(!h_done)
    {
        h_done = true;
        //call kernel 1
        //call kernel 2
        //memcpy d_done to h_done
    }


    // int a = 0;
    // while(a != 100)
    // {
    //     std::cout<<"Enter the index - ";
    //     std::cin>>a;
    //     // int i = a/MATRIX_SIZE;
    //     // printMatrixInRange(i - 1, i + 1, d_c);

    //     printNeighbors(a, d_c, vertex, edges);
    // }
}