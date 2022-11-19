/*
    Felipe Cassiano Naranjo             RA: 142489
    Pedro Tanajura Freire Meira Lima    RA: 140651
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <omp.h>
#include <sys/time.h>
#include "mpi.h"

#define SIZE 2048
#define ITERATIONS 2000

#define THREADS 8

struct timeval startSum, endSum, start, end;

void initializeGrid(bool **grid, int splitedSize);
void freeGridMemory(bool **grid);
void createGlider(bool **grid);
void createRPentonimo(bool **grid);
int getNeighbors(bool **grid, int i, int j);
void calculateNewGridState(bool **grid1, bool **grid2, int startIdx, int endIdx);
int aliveCounter(bool **grid);
void executeGame(bool **grid1, bool **grid2);
float convertTime(struct timeval start, struct timeval end);
bool isInteger(float N);

void initializeGrid(bool **grid, int splitedSize)
{
    int i, j;

    for (i = 0; i<SIZE; i++)
        grid[i] = malloc(splitedSize * sizeof (bool)) ;

    for (i = 0; i<SIZE; i++)
        for (j = 0; j<splitedSize; j++)
            grid[i][j] = 0;
}

void freeGridMemory(bool **grid)
{
    int i;

    for (i = 0; i<SIZE; i++)
        free(grid[i]);
    free(grid);
}

void createGlider(bool **grid)
{
    int lin = 1, col = 1;

    grid[lin  ][col+1] = 1;
    grid[lin+1][col+2] = 1;
    grid[lin+2][col  ] = 1;
    grid[lin+2][col+1] = 1;
    grid[lin+2][col+2] = 1;
}

void createRPentonimo(bool **grid)
{
    int lin = 10, col = 30;

    grid[lin  ][col+1] = 1;
    grid[lin  ][col+2] = 1;
    grid[lin+1][col  ] = 1;
    grid[lin+1][col+1] = 1;
    grid[lin+2][col+1] = 1;
}

int getNeighbors(bool **grid, int i, int j)
{
    if(i != 0 && j != 0 && i != SIZE-1 && j != SIZE-1)
    {
        return  grid[i-1][j-1] +    //top left
                grid[i-1][j] +      //top middle
                grid[i-1][j+1] +    //top right
                grid[i][j-1] +      //middle left
                grid[i][j+1] +      //middle right
                grid[i+1][j-1] +    //bottom left
                grid[i+1][j] +      //bottom middle
                grid[i+1][j+1];     //bottom right
    }

    return  grid[(i-1+SIZE)%SIZE][(j-1+SIZE)%SIZE] +    //top left
            grid[(i-1+SIZE)%SIZE][j] +                  //top middle
            grid[(i-1+SIZE)%SIZE][(j+1)%SIZE] +         //top right
            grid[i][(j-1+SIZE)%SIZE] +                  //middle left
            grid[i][(j+1)%SIZE] +                       //middle right
            grid[(i+1)%SIZE][(j-1+SIZE)%SIZE] +         //bottom left
            grid[(i+1)%SIZE][j] +                       //bottom middle
            grid[(i+1)%SIZE][(j+1)%SIZE];               //bottom right
}

void calculateNewGridState(bool **grid1, bool **grid2, int startIdx, int endIdx)
{
    int i, j;
    int neighbors;

    for(i = startIdx; i<endIdx; i++)
    {
        for(j = 0; j<SIZE; j++)
        {
            neighbors = getNeighbors(grid1, i, j);

            if(grid1[i][j])
            {
                if(neighbors == 2 || neighbors == 3)
                {
                    grid2[i][j] = 1;
                    continue;
                }

                grid2[i][j] = 0;
                continue;
            }

            if(neighbors == 3)
            {
                grid2[i][j] = 1;
                continue;
            }

            grid2[i][j] = 0;
        }
    }
}

int aliveCounter(bool **grid)
{
    int i, j, localSum = 0;

    #pragma omp parallel for private(i, j) reduction(+:localSum)
        for(i = 0; i<SIZE; i++)
            for(j = 0; j<SIZE; j++)
                localSum += grid[i][j];
    return localSum;
}

void executeGame(bool **grid1, bool **grid2)
{
    int i, startIdx, endIdx;
    bool **aux;
    int threadInterval = SIZE/THREADS;

    #pragma omp parallel shared(threadInterval, grid1, grid2, aux) private(i, startIdx, endIdx)
    {
        startIdx = omp_get_thread_num() * threadInterval;
        endIdx = startIdx + threadInterval;

        for(i = 0; i<ITERATIONS; i++)
        {
            calculateNewGridState(grid1, grid2, startIdx, endIdx);
            #pragma omp barrier

            #pragma omp single
            {
                aux = grid1;
                grid1 = grid2;
                grid2 = aux;
            }
            #pragma omp barrier
        }
    }
}

float convertTime(struct timeval start, struct timeval end)
{
    return (end.tv_sec - start.tv_sec) + 1e-6*(end.tv_usec - start.tv_usec);
}

bool isInteger(float N)
{
 
    // Convert float value
    // of N to integer
    int X = N;
 
    float aux = N - X;
 
    // If N is not equivalent
    // to any integer
    if (aux > 0) {
        return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
    int processId, noProcesses;
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses); 
    MPI_Comm_rank(MPI_COMM_WORLD, &processId); 

    bool **grid1, **grid2;
    int i, totalSumReduction, aux, totalRest, totalSum = 0;
    float rest, totalSumTimeReduction, totalTimeReduction, linesPerProcess;

    grid1 = malloc(SIZE * sizeof (bool*));
    linesPerProcess = SIZE/noProcesses;
    // Se a divisão está correta, divide igualmente
    // Se não, divide igualmente em todas e adiciona o resto no primeiro processo
    if(isInteger(linesPerProcess)){
        initializeGrid(grid1, linesPerProcess);
    } else{
        aux = linesPerProcess;
        rest = linesPerProcess - aux;
        if(processId == 0){
            totalRest = rest*noProcesses
            initializeGrid(grid1, (aux + totalRest));
        }else initializeGrid(grid1, aux);
    }
    createGlider(grid1);
    createRPentonimo(grid1);

    grid2 = malloc(SIZE * sizeof (bool*));
    initializeGrid(grid2);

    gettimeofday(&start, NULL);
    executeGame(grid1, grid2);
    gettimeofday(&end, NULL);

    gettimeofday(&startSum, NULL);
    totalSum = aliveCounter(grid1);
    gettimeofday(&endSum, NULL);

    float totalTime = convertTime(start, end);;
    float totalSumTime = convertTime(startSum, endSum);

    MPI_Reduce(&totalSum, &totalSumReduction, 1, MPI_FLOAT, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce(&totalTime, &totalTimeReduction, 1, MPI_FLOAT, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce(&totalSumTime, &totalSumTimeReduction, 1, MPI_FLOAT, MPI_SUM, 0,
                MPI_COMM_WORLD);

    if (processId == 0){
        printf("Celulas vivas: %d\n", totalSum);
        printf("Tempo de execucao das iteracoes: %.6fs\n", totalTimeReduction);
        printf("Tempo de execucao da soma: %.6fs\n", totalSumTimeReduction);
        printf("Tempo total: %.6fs\n", totalTimeReduction + totalSumTimeReduction);
    }

    freeGridMemory(grid1);
    freeGridMemory(grid2);

    MPI_Finalize();
    return 0;
}
