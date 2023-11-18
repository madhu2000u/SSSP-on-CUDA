#include <iostream>
#include "sssp_common.h"




float* computeMatrixMult(matElement*);
void showSSSP(float*, pathElement*);

void setUpArrays(float *d_c, int *vertex, int *edges, bool *threadMask, float* cost, pathElement* intermediateCost, int* path, matElement* minElement)
{   
    int edgeIndex = 0;
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {   
            threadMask[i * MATRIX_SIZE + j] = false;
            cost[i * MATRIX_SIZE + j] = __FLT_MAX__;
            intermediateCost[i * MATRIX_SIZE + j].value = __FLT_MAX__;
            intermediateCost[i * MATRIX_SIZE + j].pathIndex = -1;
            path[i * MATRIX_SIZE + j] = -1;

            vertex[i * MATRIX_SIZE + j] = edgeIndex;
            if((j + 1) < MATRIX_SIZE) edges[edgeIndex++] = i * MATRIX_SIZE + (j + 1);
        
            if((i + 1) < MATRIX_SIZE) edges[edgeIndex++] = (i + 1) * MATRIX_SIZE + j;

            if((j - 1) >= 0) edges[edgeIndex++] = i * MATRIX_SIZE + (j - 1);

            if((i - 1) >= 0) edges[edgeIndex++] = (i - 1) * MATRIX_SIZE + j;

        }

        threadMask[minElement[0].row * MATRIX_SIZE + minElement[0].col] = true;             //Make the thread of source vertex executable initially since that is the starting point
        cost[minElement[0].row * MATRIX_SIZE + minElement[0].col] = 0.0f;                   //Cost from source to source is 0
        intermediateCost[minElement[0].row * MATRIX_SIZE + minElement[0].col].value = 0.0f;       
        
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

void printPath(pathElement *path)
{
    for (int i = 0; i < MATRIX_SIZE; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            printf("(%d, %d),  ", i * MATRIX_SIZE + j, path[i * MATRIX_SIZE + j].pathIndex);
        }
        
        
    }
}

void findPathBoundaries(pathElement *path)
{

}

/*
    function:   reinterpret the 8 bytes pathElement(4 int + 4 float) into a value of 
                type unsigned long long int which is also of 8 bytes
*/
__device__ unsigned long long int __pathElement_as_ulli(pathElement *pathElement)
{
    unsigned long long int *ulli = reinterpret_cast<unsigned long long int*>(pathElement);
    return *ulli;

}


/*
    function:   reinterpret the 8 bytes unsigned long long int back into a value of 
                type pathElement(4 int + 4 float) which is also of 8 bytes
*/
__device__ pathElement* __ulli_as_pathElement(unsigned long long int *ulli)
{
    pathElement *element = reinterpret_cast<pathElement*>(ulli);
    return element;
}


__device__ __forceinline__ pathElement* atomicMin(pathElement *addr, pathElement* pathElement)
{
    unsigned long long int currentPathElement = __pathElement_as_ulli(addr);                                            //reinterpret to unsigned long long int since atomicCAS() supports it and a few others only
    while (pathElement->value < __ulli_as_pathElement(&currentPathElement)->value)                                         
    {
        unsigned long long int old = currentPathElement;
        currentPathElement = atomicCAS((unsigned long long int*)addr, old, __pathElement_as_ulli(pathElement));          //do atomicCAS on the reinterpreted value of ulli, if *addr == old then it puts value into addr and returns old else it does nothing and just retunrs whatever was there in addr
        if(currentPathElement == old) break;                                                                             //if value was successfully put into addr then the current thread was successful in it's atomic operation else it has to re-run with the new "current value" from addr(that might have been changed by another thread's atomic operation) and do the swapping again
    }
    return __ulli_as_pathElement(&currentPathElement);
    
}


__global__ void computeIntermediatesAndPath(float *d_c, int *vertex, int *edges, bool *threadMask, float *cost, pathElement *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;                                              //this is the global thread index to index into the entire matrix and not just for threads within each block
    if(threadMask[threadIndex])
    {
        threadMask[threadIndex] = false;
        for (int i = vertex[threadIndex]; i < vertex[threadIndex + 1]; i++)
        {   
            pathElement costPlusWeightOfCurrentThread;                                      //package cost + weight into pathElement for atomicCAS
            costPlusWeightOfCurrentThread.value = cost[threadIndex] + d_c[edges[i]];
            costPlusWeightOfCurrentThread.pathIndex = threadIndex;
            atomicMin(&intermediateCost[edges[i]], &costPlusWeightOfCurrentThread);
        }
        
    }
}

__global__ void computeFinalCosts(bool *d_done, int *vertex, int *edges, bool *threadMask, float *cost, pathElement *intermediateCost)
{
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    int threadIndex = row * MATRIX_SIZE + col;
    if(intermediateCost[threadIndex].value < cost[threadIndex])
    {
        cost[threadIndex] = intermediateCost[threadIndex].value;
        threadMask[threadIndex] = true;                                                     //since cost of this vertex changed, make it executable again to update it's neigbours
        *d_done = false;                                                                    //no atomicity required as all threads write false value only
    }

    intermediateCost[threadIndex].value = cost[threadIndex];

}

int main()
{

    int *vertex, *edges, *path;
    float *d_c;
    float *cost;
    pathElement *intermediateCostAndPath;
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
    float test_data[16] = {1.2, 5.4, 1.0f, 1.0f, 9.7, 4.9, 1.0f, 7.6, 4.0, 8.4, 1.0f, 11.5, 14.3, 2, 30.0f, 17.7};
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
    minElement[0].row = 3; minElement[0].col = 1;
    minElement[1].row = 3; minElement[1].col = 3;
    #endif   


    //Setup CUDA device memories for the data
    CHECK(cudaMallocHost(&vertex, ((MATRIX_SIZE * MATRIX_SIZE) + 1) * sizeof(int)));                    // + 1 because we need a location at the end of the vertex that stores the ending index of the edge
    CHECK(cudaMallocHost(&edges, numEdges * sizeof(int)));
    CHECK(cudaMallocHost(&threadMask, MATRIX_SIZE * MATRIX_SIZE * sizeof(bool)));
    CHECK(cudaMallocHost(&cost, size));
    CHECK(cudaMallocHost(&intermediateCostAndPath, MATRIX_SIZE * MATRIX_SIZE * sizeof(pathElement)));          //each neighbor need not have it's own cost location because the intermediate cost for a vertex is the same memory location updated by all neighbouring threads.
    CHECK(cudaMallocHost(&path, size));
    CHECK(cudaMallocHost(&d_done_ptr, sizeof(bool)));

    setUpArrays(d_c, vertex, edges, threadMask, cost, intermediateCostAndPath, path, minElement);
    vertex[MATRIX_SIZE * MATRIX_SIZE] = numEdges;                                                       //last value in vertex is total numEdges so that we can use the starting and ending index when getting the neighbors


    dim3 blockPerGrid(MATRIX_SIZE / BLOCK_DIM , MATRIX_SIZE / BLOCK_DIM);
    dim3 threadsPerBlock(BLOCK_DIM, BLOCK_DIM);

    //Start computing SSSP
    while(!h_done)
    {
        h_done = true;

        //memcpy h_done to d_done
        CHECK(cudaMemcpy(d_done_ptr, &h_done, sizeof(bool), cudaMemcpyHostToDevice));

        //call kernel 1
        computeIntermediatesAndPath<<<blockPerGrid, threadsPerBlock>>>(d_c, vertex, edges, threadMask, cost, intermediateCostAndPath);
        cudaDeviceSynchronize();

        //call kernel 2
        computeFinalCosts<<<blockPerGrid, threadsPerBlock>>>(d_done_ptr, vertex, edges, threadMask, cost, intermediateCostAndPath);
        cudaDeviceSynchronize();

        //memcpy d_done to h_done
        CHECK(cudaMemcpy(&h_done, d_done_ptr, sizeof(bool), cudaMemcpyDeviceToHost));
        

    }

    #if(TEST)
        printPath(intermediateCost);
    #endif

    showSSSP(d_c, intermediateCostAndPath);
    
    // pathElement x;
    // x.value = 1.2f;
    // x.pathIndex = 5;

    // std::cout << std::bitset<32>(*(reinterpret_cast<int*>(&x.value))) << std::endl;
    // std::cout << std::bitset<32>(x.pathIndex) <<std::endl;

    // unsigned long long int y = __pathElement_as_ulli(&x);
    // std::cout << std::bitset<64>(y) << std::endl;

    // float xval = __ulli_as_pathElement(&y)->value;
    // int xpathIndex = (__ulli_as_pathElement(&y))->pathIndex;

    // std::cout << std::bitset<32>(*(reinterpret_cast<int*>(&xval))) << std::endl;
    // std::cout << std::bitset<32>(xpathIndex) << std::endl;
    
    
    printf("\ncost of target - %f\n", d_c[minElement[0].row * MATRIX_SIZE + minElement[0].col] + cost[minElement[1].row * MATRIX_SIZE + minElement[1].col] - d_c[minElement[1].row * MATRIX_SIZE + minElement[1].col]);       //include source's weight and exclude target's weight

}