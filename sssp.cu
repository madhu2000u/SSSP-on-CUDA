#include <iostream>
#include "sssp_common.h"




float* computeMatrixMult(matElement*);

void setUpArrays(float *d_c, int *vertex, int *edges, bool *threadMask, float* cost, float* intermediateCost, int* path, matElement* minElement)
{   
    int edgeIndex = 0;
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {   
            threadMask[i * MATRIX_SIZE + j] = false;
            cost[i * MATRIX_SIZE + j] = __FLT_MAX__;
            intermediateCost[i * MATRIX_SIZE + j] = __FLT_MAX__;
            path[i * MATRIX_SIZE + j] = -1;

            vertex[i * MATRIX_SIZE + j] = edgeIndex;
            if((j + 1) < MATRIX_SIZE) edges[edgeIndex++] = i * MATRIX_SIZE + (j + 1);
        
            if((i + 1) < MATRIX_SIZE) edges[edgeIndex++] = (i + 1) * MATRIX_SIZE + j;

            if((j - 1) >= 0) edges[edgeIndex++] = i * MATRIX_SIZE + (j - 1);

            if((i - 1) >= 0) edges[edgeIndex++] = (i - 1) * MATRIX_SIZE + j;

        }

        threadMask[minElement[0].row * MATRIX_SIZE + minElement[0].col] = true;             //Make the thread of source vertex executable initially since that is the starting point
        cost[minElement[0].row * MATRIX_SIZE + minElement[0].col] = 0.0f;                   //Cost from source to source is 0
        intermediateCost[minElement[0].row * MATRIX_SIZE + minElement[0].col] = 0.0f;       
        
    }
    
}

void printNeighbors(int index, float *d_c, int *vertex, int* edges)
{
    for (int i = vertex[index]; i < vertex[index + 1]; i++)
    {
        std::cout<<i<<std::endl;
        std::cout<<d_c[edges[i]]<<std::endl;
        std::cout<<"\n";
    }
    
}

void printPath(int *path)
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("(%d, %d),  ", i * MATRIX_SIZE + j, path[i * MATRIX_SIZE + j]);
        }
        
        
    }
}

typedef struct pathElement
{
    float value;
    int pathIndex;
};


unsigned long long int __pathElement_as_ulli(pathElement *pathElement)
{
    unsigned long long int ulli = (unsigned long long int) pathElement;
    return ulli;

}

__device__ __forceinline__ float atomicMinPath(float *intermediateAddr, int *pathAddr, float value, int indexOfCallingThreadToPointToForCorrespondingAddrToGivePath)
{
    int currentIntermediateCost = __float_as_int(*intermediateAddr);                                            //reinterpret to int since atomicCAS() requires int
    int currentPath = *pathAddr;
    while (value < __int_as_float(currentIntermediateCost))                                         
    {//atomicCAS_block()
        int old = currentPath;
        currentPath = atomicCAS(pathAddr, old, indexOfCallingThreadToPointToForCorrespondingAddrToGivePath);                //if *addr == old then it puts value into addr and returns old else it does nothing and just retunrs whatever was there in addr
        if(currentPath == old) break;                                                   //if value was successfully put into addr then the current thread was successful in it's atomic operation else it has to re-run with the new "current value" from addr(that might have been changed by another thread's atomic operation) and do the swapping again
    }
    return currentPath;
    
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

__global__ void computePaths(float *d_c, int *vertex, int *edges, bool *threadMask, float *cost, float *intermediateCost, int* path)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;                              //this is the global thread index to index into the entire matrix and not just for threads within each block
    if(threadMask[threadIndex])
    {
        //threadMask[threadIndex] = false;
        for (int i = vertex[threadIndex]; i < vertex[threadIndex + 1]; i++)
        {   
            atomicMinPath(&intermediateCost[edges[i]], &path[edges[i]], (cost[threadIndex] + d_c[edges[i]]), threadIndex);
        }
        
    }
}

__global__ void computeIntermediates(float *d_c, int *vertex, int *edges, bool *threadMask, float *cost, float *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;                              //this is the global thread index to index into the entire matrix and not just for threads within each block
    if(threadMask[threadIndex])
    {
        threadMask[threadIndex] = false;
        for (int i = vertex[threadIndex]; i < vertex[threadIndex + 1]; i++)
        {   
            atomicMin(&intermediateCost[edges[i]], (cost[threadIndex] + d_c[edges[i]]));
        }
        
    }
}

__global__ void computeFinalCostsAndPath(bool *d_done, int *vertex, int *edges, bool *threadMask, float *cost, float *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;
    if(intermediateCost[threadIndex] < cost[threadIndex])
    {
        cost[threadIndex] = intermediateCost[threadIndex];
        threadMask[threadIndex] = true;                                                 //since cost of this vertex changed, make it executable again to update it's neigbours
        *d_done = false;                                                                 //no atomicity required as all threads write false value only
    }

    intermediateCost[threadIndex] = cost[threadIndex];

}

int main()
{

    int *vertex, *edges, *path;
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
    float test_data[16] = {1.2, 5.4, 1.0f, 1.0f, 9.7, 4.9, 1.0f, 7.6, 6.1, 8.4, 1.0f, 11.5, 14.3, 15.8, 1.0f, 17.7};
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
    minElement[0].row = 0; minElement[0].col = 3;
    minElement[1].row = 3; minElement[1].col = 3;
    #endif

    // vertex = (int*)malloc(((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int));                            
    // threadMask = new bool[MATRIX_SIZE * MATRIX_SIZE];
    // edges = (int*)malloc( numEdges * sizeof(int));

    // cost = (float*)malloc(size);
    // intermediateCost = (float*)malloc(size);                                                          

   

    

    //Setup CUDA device memories for the data
    CHECK(cudaMallocHost(&vertex, ((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int)));                     // + 1 because we need a location at the end of the vertex that stores the ending index of the edge
    CHECK(cudaMallocHost(&edges, numEdges * sizeof(int)));
    CHECK(cudaMallocHost(&threadMask, MATRIX_SIZE * MATRIX_SIZE * sizeof(bool)));
    CHECK(cudaMallocHost(&cost, size));
    CHECK(cudaMallocHost(&intermediateCost, numEdges * sizeof(float)));                                  //each neighbor need not have it's own cost location because the intermediate cost for a vertex is the same memory location updated by all neighbouring threads.
    CHECK(cudaMallocHost(&path, size));
    CHECK(cudaMallocHost(&d_done_ptr, sizeof(bool)));

    setUpArrays(d_c, vertex, edges, threadMask, cost, intermediateCost, path, minElement);
    vertex[MATRIX_SIZE * MATRIX_SIZE] = numEdges;                                                     //last value in vertex is total numEdges so that we can use the starting and ending index when getting the neighbors


    dim3 blockPerGrid(MATRIX_SIZE / BLOCK_DIM , MATRIX_SIZE / BLOCK_DIM);
    dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM);

    //Start computing SSSP
    while(!h_done)
    {
        h_done = true;
        //memcpy h_done to d_done
        CHECK(cudaMemcpy(d_done_ptr, &h_done, sizeof(bool), cudaMemcpyHostToDevice));

        computePaths<<<blockPerGrid, threadsPerBlock>>>(d_c, vertex, edges, threadMask, cost, intermediateCost, path);
        cudaDeviceSynchronize();

        //call kernel 1
        computeIntermediates<<<blockPerGrid, threadsPerBlock>>>(d_c, vertex, edges, threadMask, cost, intermediateCost);
        cudaDeviceSynchronize();

        //call kernel 2
        computeFinalCostsAndPath<<<blockPerGrid, threadsPerBlock>>>(d_done_ptr, vertex, edges, threadMask, cost, intermediateCost);
        cudaDeviceSynchronize();

        //memcpy d_done to h_done
        CHECK(cudaMemcpy(&h_done, d_done_ptr, sizeof(bool), cudaMemcpyDeviceToHost));
        
        printPath(path);
        printf("\n");

    }
    
    
    printf("\ncost of target - %f\n", d_c[minElement[0].row * MATRIX_SIZE + minElement[0].col] + cost[minElement[1].row * MATRIX_SIZE + minElement[1].col] - d_c[minElement[1].row * MATRIX_SIZE + minElement[1].col]);       //include source's weight and exclude target's weight

}