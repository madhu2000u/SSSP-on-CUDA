#include <iostream>
#include "sssp_common.h"


// __device__ bool d_done = false;

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
            intermediateCost[i * MATRIX_SIZE + j] = __FLT_MAX__;

            vertex[i * MATRIX_SIZE + j] = edgeIndex;
            if((j + 1) < MATRIX_SIZE)
            {
                //intermediateCost[edgeIndex] = __FLT_MAX__;            //intemediateCost array size is same as that of edges so do it in this step itself so that you can fill all the spots in the array just like edges array instead of making another loop for it
                edges[edgeIndex++] = i * MATRIX_SIZE + (j + 1);
            }
            if((i + 1) < MATRIX_SIZE)
            {
                //intermediateCost[edgeIndex] = __FLT_MAX__;
                edges[edgeIndex++] = (i + 1) * MATRIX_SIZE + j;
            }
            if((j - 1) >= 0)
            {
                //intermediateCost[edgeIndex] = __FLT_MAX__;
                edges[edgeIndex++] = i * MATRIX_SIZE + (j - 1);
            }
            if((i - 1) >= 0)
            {
                //intermediateCost[edgeIndex] = __FLT_MAX__;
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

__device__ __forceinline__ float atomicMin(float *addr, float value)
{
    int current = __float_as_int(*addr);                                            //reinterpret to int since atomicCAS() requires int
    while (value < __int_as_float(current))                                         
    {
        int old = current;
        current = atomicCAS((int*)addr, old, __float_as_int(value));                //if *addr == old then it puts value into addr and returns old else it does nothing and just retunrs whatever was there in addr
        if(current == old) break;                                                   //if value was successfully put into addr then the current thread was successful in it's atomic operation else it has to re-run with the new "current value" from addr(that might have been changed by another thread's atomic operation) and do the swapping again
    }
    return __int_as_float(current);
    
}

__global__ void computeIntermediates(float *d_c, int *vertex, int *edges, bool *threadMask, float *cost, float *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;                              //this is the global thread index to index into the entire matrix and not just for threads within each block
    // if(threadIndex == (407 * MATRIX_SIZE + 976))printf("threadindex - %d\n", threadMask[threadIndex]);
    if(threadMask[threadIndex])
    {//printf("threadIndex - %d\n", threadIndex);
        threadMask[threadIndex] = false;
        for (int i = vertex[threadIndex]; i < vertex[threadIndex + 1]; i++)
        {   
            // printf("\nintermedate cost before swap - %f\n", intermediateCost[edges[i]]);
            atomicMin(&intermediateCost[edges[i]], (cost[threadIndex] + d_c[edges[i]]));
            // printf("intermedate cost after swap - %f\n", intermediateCost[edges[i]]);
            // printf("computing for each neighbour");
        }
        
    }
    // while(true)if(threadIndex==0){printf("looping");}

    // __syncthreads();
    // if(threadIndex == 0)printf("done with intermediates\n");
}

__global__ void computeFinalCostsAndPath(bool *d_done, int *vertex, int *edges, bool *threadMask, float *cost, float *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;
    if(threadIndex == 0)printf("started 2nd kernel from cuda\n");
    if(intermediateCost[threadIndex] < cost[threadIndex])
    {
        cost[threadIndex] = intermediateCost[threadIndex];
        threadMask[threadIndex] = true;                                                 //since cost of this vertex changed, make it executable again to update it's neigbours
        *d_done = false;                                                                 //no atomicity required as all threads write false value only
        // printf("inside here");
    }

    intermediateCost[threadIndex] = cost[threadIndex];

    // __syncthreads();
}

int main()
{   /*backup - int *h_vertex, *d_vertex, *h_edges, *d_vertex;
    float *d_c, *h_cost, *d_cost, *h_intermediateCost, *d_intermediateCost;
    bool *h_threadMask,;*/

    int *vertex, *edges;
    float *d_c;
    float *cost, *intermediateCost;
    bool *threadMask;

    bool h_done = false;
    bool *d_done_ptr;

    matElement minElement[2];
    size_t size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
    //Note: Diagnonal neighbours are not considered
    int numEdges = (4 * 2                                       /*Each corner values in matrix has 2 neighbours*/ 
                    + ((MATRIX_SIZE - 2) * 3) * 4               /*Each element of the 4 boundary sides excluding the 2 corner elements for each boundary side has 3 neighbours*/ 
                    + (MATRIX_SIZE - 2) * 4 * (MATRIX_SIZE - 2) /*Each element not on the boundary has 4 neighbours*/);
    
    //Compute and get the pointer to the result matrix of the matrix muliplications
    d_c = computeMatrixMult(minElement);

    //Use test data for sssp checking
    #if(TEST)
    float test_data[16] = {1, 5, 3, 2, 9, 4, 10, 7, 6, 8, 12, 11, 14, 15, 16, 17};
    CHECK(cudaMemcpy(d_c, &test_data, size, cudaMemcpyHostToDevice));

    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            std::cout<<d_c[i * MATRIX_SIZE + j]<<"  ";
        }
        printf("\n");
        
    }
    //Test values for source and target
    minElement[0].row = 2; minElement[0].col = 3;
    minElement[1].row = 1; minElement[1].col = 1;
    #endif

    vertex = (int*)malloc(((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int));                            // + 1 because we need a location at the end of the vertex that stores the ending index of the edge
    threadMask = new bool[MATRIX_SIZE * MATRIX_SIZE];
    edges = (int*)malloc( numEdges * sizeof(int));

    cost = (float*)malloc(size);
    intermediateCost = (float*)malloc(size);                                                          //each neighbor need not have it's own cost location because the intermediate cost for a vertex is the same memory location updated by all neighbouring threads.

   

    

    //Setup CUDA device memories for the data
    CHECK(cudaMallocHost(&vertex, ((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int)));
    CHECK(cudaMallocHost(&edges, numEdges * sizeof(int)));
    CHECK(cudaMallocHost(&threadMask, MATRIX_SIZE * MATRIX_SIZE * sizeof(bool)));
    CHECK(cudaMallocHost(&cost, size));
    CHECK(cudaMallocHost(&intermediateCost, numEdges * sizeof(float)));
    CHECK(cudaMallocHost(&d_done_ptr, sizeof(bool)));

    setUpArrays(d_c, vertex, edges, threadMask, cost, intermediateCost, minElement);
    vertex[MATRIX_SIZE * MATRIX_SIZE] = numEdges;                                                     //last value in vertex is total numEdges so that we can use the starting and ending index when getting the neighbors

    std::cout<<"threadmask from host - "<<threadMask[minElement[0].row * MATRIX_SIZE + minElement[0].col]<<std::endl;

    dim3 blockPerGrid(MATRIX_SIZE / BLOCK_DIM , MATRIX_SIZE / BLOCK_DIM);
    dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM);

    //Start computing SSSP
    while(!h_done)
    {
        h_done = true;
        //memcpy h_done to d_done
        CHECK(cudaMemcpy(d_done_ptr, &h_done, sizeof(bool), cudaMemcpyHostToDevice));
        //call kernel 1
        computeIntermediates<<<blockPerGrid, threadsPerBlock>>>(d_c, vertex, edges, threadMask, cost, intermediateCost);
        //call kernel 2
        std::cout<<"after 1st kernel"<<std::endl;
        std::cout<<"intermedeate cost after first kernel - "<<intermediateCost[minElement[1].row * MATRIX_SIZE + minElement[1].col + 1]<<std::endl;
        cudaDeviceSynchronize();
        std::cout<<"synch done - "<<std::endl;
        // std::cout<<"intermedate from host - "<<intermediateCost[10200]<<std::endl;
        computeFinalCostsAndPath<<<blockPerGrid, threadsPerBlock>>>(d_done_ptr, vertex, edges, threadMask, cost, intermediateCost);
        //memcpy d_done to h_done
        cudaDeviceSynchronize();
        std::cout<<"after 2nd kernel"<<std::endl;
        CHECK(cudaMemcpy(&h_done, d_done_ptr, sizeof(bool), cudaMemcpyDeviceToHost));

        std::cout<<"h_done - "<<h_done<<std::endl;
    }

    std::cout<<"cost of target - "<<cost[minElement[1].row * MATRIX_SIZE + minElement[1].col] - d_c[minElement[1].row * MATRIX_SIZE + minElement[1].col]<<std::endl;                        //exclude target's weight


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