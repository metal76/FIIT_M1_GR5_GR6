#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <time.h>
using namespace std;
void Rec_Mult(int *C, const int *A, const int *B, int n, int rowsize)
{
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			for (size_t k = 0; k < n; k++)
			{
				C[i*n+j] += A[i*n+k] * B[k*n+j];
			}
		}
	}
}


#define ROW_COUNT 2880
int Block_Size = 2;

void printtrMatrix(const char *name, const int *mat,int block)
{
	printf("%s:\n", name);

	for (int i = 0; i < ((ROW_COUNT / block)*((ROW_COUNT / block) + 1) / 2)*block*block; ++i)
	{
			printf("%4d ", mat[i]);
			if (1%10==0)
				printf("\n");
	}
	printf("\n");
}

void printMatrix(const char *name, const int *mat)
{
	printf("%s:\n", name);

	for (int i = 0; i < ROW_COUNT; ++i)
	{
		for (int j = 0; j < ROW_COUNT; ++j)
		{
			printf("%4d", mat[i * ROW_COUNT + j]);
		}
		printf("\n");
	}
	printf("\n");
}

int *  CreateA()
{
	int *  res = new int[ROW_COUNT * ROW_COUNT];
	for (size_t i = 0; i < ROW_COUNT; i++)
	{
		for (size_t j = 0; j < ROW_COUNT; j++)
		{
			res[i*ROW_COUNT + j] = (j >= i) ? 1 : 0;
		}
	}
	return res;
}
int *  CreateB()
{
	int *  res = new int[ROW_COUNT * ROW_COUNT];
	for (size_t i = 0; i < ROW_COUNT; i++)
	{
		for (size_t j = 0; j < ROW_COUNT; j++)
		{
			int value = (i >= j) ? 1 : 0;
			res[i*ROW_COUNT + j] = (i >= j) ? 1 : 0;
		}
	}
	return res;
}

int *  CreateAtr(const int * A, int block)
{
	int length = ((ROW_COUNT / block+1)*((ROW_COUNT / block) + 2) / 2)*block*block;
	int *res = new int[((ROW_COUNT/block + 1)*((ROW_COUNT / block)+2) / 2)*block*block];
	int num = 0;
	for (size_t j = 0; j < ROW_COUNT; j+= block)
	{
		for (size_t i = 0; i <= j; i += block)
		{
			for (size_t n = 0; n < block; n++)
			{
				for (size_t m = 0; m < block; m++)
				{
					res[num++]=A[(i+n)*ROW_COUNT+j+m];
				}
			}
		}
	}

	return res;
}

int *  CreateBtr(const int * B, int block)
{
	int length = ((ROW_COUNT / block)*((ROW_COUNT / block) + 1) / 2)*block*block;
	int *  res = new int[((ROW_COUNT / block)*((ROW_COUNT / block) + 1) / 2)*block*block];
	int num = 0;
	for (size_t j = 0; j < ROW_COUNT; j += block)
	{
		for (size_t i = j; i < ROW_COUNT; i += block)
		{
			for (size_t n = 0; n < block; n++)
			{
				for (size_t m = 0; m < block; m++)
				{
					res[num++] = B[(i + n)*ROW_COUNT + j + m];
				}
			}
		}
	}
	return res;
}

int main()
{
	int* block_Sizes = new int[10] {4,6,8,12,16 ,20, 24, 36, 48, 72};
	for (size_t bi = 0; bi < 10; bi++)
	{
		int Block_count = ROW_COUNT / block_Sizes[bi];
		printf("Block size %4d ", block_Sizes[bi]);
		Block_Size = block_Sizes[bi];
		int Block_Capasity = Block_Size * Block_Size;
		const int * matA;
		matA = CreateA();
		int * matAtr = CreateAtr(matA, Block_Size);
		//printtrMatrix("Matrix Atr", matAtr, Block_Size);

		const int * matB;
		matB = CreateB();
		int * matBtr = CreateBtr(matB, Block_Size);
		int * matC = new int[ROW_COUNT * ROW_COUNT];
		//int matD[ROW_COUNT * ROW_COUNT];
		//memset(matD, 0, sizeof(matD));

		//printtrMatrix("Matrix Btr", matBtr, Block_Size);
		//Rec_Mult(matD, matA, matB, ROW_COUNT, ROW_COUNT);
		clock_t time;

		time = clock();
		for (size_t i = 0; i < Block_count; i++)
		{
			for (int j = 0; j < i; j++)
			{
				for (int k = i; k < Block_count; k++)
				{
					int algsum = ((k)*(k + 1) / 2 + i);
#pragma omp parallel for
					for (int row = 0; row < Block_Size; row++) {
						for (int col = 0; col < Block_Size; col++) {
							for (int inner = 0; inner < Block_Size; inner++) {								
								matC[(Block_Size*i + row)*ROW_COUNT + Block_Size * j + col] +=
									matAtr[algsum*Block_Capasity + row * Block_Size + inner] *
									matBtr[((2 * (Block_count) - j + 1)*j / 2 + k - j)*Block_Capasity + Block_Size * inner + col];
							}
						}
					}
				}
			}
			for (int j = i; j < Block_count; j++)
			{
				for (int k = j; k < Block_count; k++)
				{
					int algsum = ((k)*(k + 1) / 2 + i);
#pragma omp parallel for
					for (int row = 0; row < Block_Size; row++) {
						for (int col = 0; col < Block_Size; col++) {
							for (int inner = 0; inner < Block_Size; inner++) {
								matC[(Block_Size*i + row)*ROW_COUNT + Block_Size * j + col] += 
									matAtr[algsum*Block_Capasity + row * Block_Size + inner] *
									matBtr[((2 * Block_count - j + 1)*j / 2 + k - j)*Block_Capasity + Block_Size * inner + col];
							}
						}
					}
				}
			}

		}
		time = clock() - time;
		printf(" time %f \n", (double)time / CLOCKS_PER_SEC); //время выполнения "каких-то действий"

	}
	//printMatrix("Matrix A", matA);
	//printMatrix("Matrix B", matB);
	//printMatrix("Multiply My", matC);
	//printMatrix("Multiply Correct", matD);
	
	system("pause");
	return 0;
}