#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "strassens.h"

static bool ONLY_STATS = false;
static int MAX_NUM = 100;
static int MIN_NUM = 10;

void caller(matrix* (funct)(const matrix*, const matrix*), char* text, matrix* m1, matrix* m2)
{
	int i, j;
	clock_t begin;
	clock_t end;
	matrix* prod;
	double time_taken;

	begin = clock();
	prod = funct(m1, m2);
	end = clock();
	time_taken = (double)(end - begin) / CLOCKS_PER_SEC;
	if(ONLY_STATS)
		printf("%s,%lf\n",text, time_taken);
	else
	{
		char path[50];
		char snum[6];
		itoa(m1->row_size, snum, 10);
		strcpy(path, "Output\\");
		strcat(path, text);
		strcat(path, "_");
		strcat(path, snum);
		strcat(path, ".txt");
		printf("Time taken by %s: %f seconds\n", text, time_taken);
		//output_to_file(prod, path);
		printf("\n");
	}
	free_matrix(prod);
}

int main(int argc, char** argv)
{
	int i, j;
	int dimensions[] = {16, 32, 64, 128, 256, 512, 1024, 2048, 4096};

	switch(argc)
	{
		case 4: MAX_NUM = atoi(argv[3]);
		case 3: MIN_NUM = atoi(argv[2]);
		case 2: ONLY_STATS = (strcmpi(argv[1], "true") == 0);
	}

	if(MIN_NUM > MAX_NUM)
		exit(-1);

	srand(time(NULL));

	for(i=0; i<9; i++)
	{
		matrix* m1;
		matrix* m2;

		m1 = generate_random_matrix(dimensions[i], dimensions[i], MIN_NUM, MAX_NUM);
		m2 = generate_random_matrix(dimensions[i], dimensions[i], MIN_NUM, MAX_NUM);

		//display_matrix(m1);
		//display_matrix(m2);

		if(ONLY_STATS)
			printf("%d\n", dimensions[i]);
		else
			printf("Dimension: %dx%d\n", dimensions[i], dimensions[i]);

		caller(brute_force_multiplication, "brute_force_multiplication", m1, m2);
		caller(divide_and_conquer, "divide_and_conquer", m1, m2);
		caller(divide_and_conquer_parallel, "divide_and_conquer_parallel", m1, m2);
		caller(strassens, "strassens", m1, m2);
		caller(strassens_parallel, "strassens_parallel", m1, m2);

		free_matrix(m1);
		free_matrix(m2);

		if(!ONLY_STATS)
		{
			for(j=1; j<=30; j++)
				printf("- ");
			printf("\n\n");
		}
	}
	return 0;
}