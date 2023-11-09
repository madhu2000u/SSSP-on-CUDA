#include <iostream>
#include "sssp_common.h"


float* computeMatrixMult();

void setUpArrays(float *d_c, int *vertex, int *edges)
{   
    int edgeIndex = 0;
    int validNeighborCount = 0;
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            vertex[i * MATRIX_SIZE + j] = edgeIndex;
            if((j + 1) < MATRIX_SIZE)
            {
                edges[edgeIndex++] = i * MATRIX_SIZE + (j + 1);
            }
            if((i + 1) < MATRIX_SIZE)
            {
                edges[edgeIndex++] = (i + 1) * MATRIX_SIZE + j;
            }
            if((j - 1) >= 0)
            {
                edges[edgeIndex++] = i * MATRIX_SIZE + (j - 1);
            }
            if((i - 1) >= 0)
            {
                edges[edgeIndex++] = (i - 1) * MATRIX_SIZE + j;
            }
        }
        
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

int main()
{   
    int *vertex, *edges;
    float *d_c;
    bool *threadMask;
    int numEdges = (4 * 2 + ((MATRIX_SIZE - 2) * 3) * 4 + (MATRIX_SIZE - 2) * 4 * (MATRIX_SIZE - 2));
    
    d_c = computeMatrixMult();

    vertex = (int*)malloc(((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int));                            // + 1 because we need a location at the end of the vertex that stores the ending index of the edge
    threadMask = new bool[MATRIX_SIZE * MATRIX_SIZE];
    edges = (int*)malloc( numEdges * sizeof(int));

    float test_data[16] = {1, 5, 3, 2, 9, 4, 10, 7, 6, 8, 12, 11, 14, 15, 16, 17};

    setUpArrays(d_c, vertex, edges);
    vertex[MATRIX_SIZE * MATRIX_SIZE] = numEdges;                                                     //last value in vertex is total numEdges so that we can use the starting and ending index when getting the neighbors
    
    int a = 0;
    while(a != 100)
    {
        std::cout<<"Enter the index - ";
        std::cin>>a;
        // int i = a/MATRIX_SIZE;
        // printMatrixInRange(i - 1, i + 1, d_c);

        printNeighbors(a, d_c, vertex, edges);
    }
}