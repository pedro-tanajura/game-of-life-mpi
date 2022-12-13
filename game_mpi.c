%%writefile game-mpi.c

/*
    Felipe Cassiano Naranjo             RA: 142489
    Pedro Tanajura Freire Meira Lima    RA: 140651
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include "mpi.h"

#define SIZE 2048
#define ITERATIONS 2000

struct timeval startSum, endSum, start, end;

void initializeGrid(bool **grid, int linesPerProcess);
void freeGridMemory(bool **grid, int linesPerProcess);
void createGlider(bool **grid);
void createRPentonimo(bool **grid);
int getNeighbors(bool **grid, int i, int j, bool *bot_border, bool *top_border, int linesPerProcess);
void calculateNewGridState(bool **grid1, bool **grid2, int linesPerProcess);
int aliveCounter(bool **grid, int linesPerProcess);
void executeGame(bool **grid1, bool **grid2, int linesPerProcess);
float convertTime(struct timeval start, struct timeval end);
void print_grid(bool **grid, int linesPerProcess);

void print_grid(bool **grid, int linesPerProcess)
{
    int i, j;
    
    for(i = 0; i<linesPerProcess; i++){
        for(j = 0; j<SIZE; j++)
            printf("%d ", grid[i][j]);
        printf("\n");
    }
}

void initializeGrid(bool **grid, int linesPerProcess)
{
    int i, j;

    for (i = 0; i<linesPerProcess; i++)
        grid[i] = malloc(SIZE * sizeof (bool)) ;

    for (i = 0; i<linesPerProcess; i++)
        for (j = 0; j<SIZE; j++)
            grid[i][j] = 0;
}

void freeGridMemory(bool **grid, int linesPerProcess)
{
    int i;

    for (i = 0; i<linesPerProcess; i++)
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

int getNeighbors(bool **grid, int i, int j, bool *bot_border, bool *top_border, int linesPerProcess)
{
    if(i != 0 && j != 0 && i != linesPerProcess-1 && j != SIZE-1) // Pega parte interna do grid
    {
        return  grid[i-1][j-1] +                            //top left
                grid[i-1][j] +                              //top middle
                grid[i-1][j+1] +                            //top right
                grid[i][j-1] +                              //middle left
                grid[i][j+1] +                              //middle right
                grid[i+1][j-1] +                            //bottom left
                grid[i+1][j] +                              //bottom middle
                grid[i+1][j+1];                             //bottom right
    } 
    int soma = 0;
    if(i == 0){
        return  top_border[(j-1+SIZE)%SIZE] +               //top left
                top_border[j] +                             //top middle
                top_border[(j+1)%SIZE] +                    //top right
                grid[i][(j-1+SIZE)%SIZE] +                  //middle left
                grid[i][(j+1)%SIZE] +                       //middle right
                grid[i+1][j-1] +                            //bottom left
                grid[i+1][j] +                              //bottom middle
                grid[i+1][j+1];                             //bottom right
    } if (i == linesPerProcess-1){
        return  grid[i-1][j-1] +                            //top left
                grid[i-1][j] +                              //top middle
                grid[i-1][j+1] +                            //top right
                grid[i][(j-1+SIZE)%SIZE] +                  //middle left
                grid[i][(j+1)%SIZE] +                       //middle right
                bot_border[(j-1+SIZE)%SIZE] +               //bottom left
                bot_border[j] +                             //bottom middle
                bot_border[(j+1)%SIZE];                     //bottom right
    }

        return  grid[i-1][j-1] +                            //top left
                grid[i-1][j] +                              //top middle
                grid[i-1][j+1] +                            //top right
                grid[i][(j-1+SIZE)%SIZE] +                  //middle left
                grid[i][(j+1)%SIZE] +                       //middle right
                grid[i+1][j-1] +                            //bottom left
                grid[i+1][j] +                              //bottom middle
                grid[i+1][j+1];                             //bottom right

}

void calculateNewGridState(bool **grid1, bool **grid2, int linesPerProcess)
{
    int i, j;
    int neighbors;
    int processId, noProcesses;
    MPI_Status status;
    bool *top_border = malloc(SIZE * sizeof (bool));
    bool *bot_border = malloc(SIZE * sizeof (bool));

    MPI_Comm_size(MPI_COMM_WORLD, &noProcesses); 
    MPI_Comm_rank(MPI_COMM_WORLD, &processId); 

    MPI_Sendrecv(grid1[0], SIZE, MPI_C_BOOL, ((processId+1)%noProcesses), 10,
	            bot_border, SIZE,  MPI_C_BOOL, ((processId-1+noProcesses)%noProcesses), 10,
	            MPI_COMM_WORLD, &status);

    MPI_Sendrecv(grid1[linesPerProcess-1], SIZE, MPI_C_BOOL, ((processId-1+noProcesses)%noProcesses), 11,
	            top_border, SIZE,  MPI_C_BOOL, ((processId+1)%noProcesses), 11,
	            MPI_COMM_WORLD, &status);

    for(i = 0; i<linesPerProcess; i++)
    {
        for(j = 0; j<SIZE; j++)
        {
            neighbors = getNeighbors(grid1, i, j, bot_border, top_border, linesPerProcess);

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

int aliveCounter(bool **grid, int linesPerProcess)
{
    int i, j, localSum = 0;
    
    for(i = 0; i<linesPerProcess; i++)
        for(j = 0; j<SIZE; j++)
            localSum += grid[i][j];
    return localSum;
}

void executeGame(bool **grid1, bool **grid2, int linesPerProcess)
{
    int i;
    bool **aux;

    for(i = 0; i<ITERATIONS; i++)
    {
        printf("%d\n", i);
        calculateNewGridState(grid1, grid2, linesPerProcess);
        aux = grid1;
        grid1 = grid2;
        grid2 = aux;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

float convertTime(struct timeval start, struct timeval end)
{
    return (end.tv_sec - start.tv_sec) + 1e-6*(end.tv_usec - start.tv_usec);
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

    linesPerProcess = SIZE/noProcesses;

    grid1 = malloc(linesPerProcess * sizeof (bool*));
    initializeGrid(grid1, linesPerProcess);

    grid2 = malloc(linesPerProcess * sizeof (bool*));
    initializeGrid(grid2, linesPerProcess);
    
    if(processId == 0){
        createGlider(grid1);
        createRPentonimo(grid1);
    }

    gettimeofday(&start, NULL);
    executeGame(grid1, grid2, linesPerProcess);
    gettimeofday(&end, NULL);

    gettimeofday(&startSum, NULL);
    totalSum = aliveCounter(grid1, linesPerProcess);
    gettimeofday(&endSum, NULL);

    float totalTime = convertTime(start, end);;
    float totalSumTime = convertTime(startSum, endSum);

    MPI_Reduce(&totalSum, &totalSumReduction, 1, MPI_INT, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce(&totalTime, &totalTimeReduction, 1, MPI_FLOAT, MPI_SUM, 0,
                MPI_COMM_WORLD);
    MPI_Reduce(&totalSumTime, &totalSumTimeReduction, 1, MPI_FLOAT, MPI_SUM, 0,
                MPI_COMM_WORLD);

    if (processId == 0){
        printf("Celulas vivas: %d\n", totalSumReduction);
        printf("Tempo de execucao das iteracoes: %.6fs\n", totalTimeReduction);
        printf("Tempo de execucao da soma: %.6fs\n", totalSumTimeReduction);
        printf("Tempo total: %.6fs\n", totalTimeReduction + totalSumTimeReduction);
    }

    freeGridMemory(grid1, linesPerProcess);
    freeGridMemory(grid2, linesPerProcess);

    MPI_Finalize();
    return 0;
}