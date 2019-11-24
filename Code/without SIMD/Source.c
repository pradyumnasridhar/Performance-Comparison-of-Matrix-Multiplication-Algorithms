#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <windows.h>
#include "strassens.h"

matrix* read_matrix(int row_size, int col_size)
{
	int i, j;
	matrix* m = (matrix*)malloc(sizeof(matrix));
	m->mat = (double**) malloc(sizeof(double*) * row_size);
	m->col_size = col_size;
	m->row_size = row_size;

	for (i = 0; i < m->col_size; i++)
	{
		m->mat[i] = (double*) malloc(sizeof(double) * col_size);
		for (j = 0; j < m->row_size; j++)
			scanf("%lf", &m->mat[i][j]);
	}

	return m;
}

void display_matrix(const matrix* m)
{
	int i, j;

	for (i = 0; i < m->row_size; i++)
	{
		for (j = 0; j < m->col_size; j++)
		{
			printf("%.3lf\t", m->mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void output_to_file(const matrix* m, const char* path)
{
	int i, j;
	FILE* f = fopen(path, "w");

	for (i = 0; i < m->row_size; i++)
	{
		for (j = 0; j < m->col_size; j++)
		{
			fprintf(f, "%.3lf\t", m->mat[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	printf("Output saved to %s\n", path);
}

void free_matrix(matrix* m)
{
	int i;

	for (i = 0; i < m->row_size; i++)
		free(m->mat[i]);
	free(m->mat);
	free(m);
}

matrix* generate_random_matrix(int row_size, int col_size, const int MIN, const int MAX)
{
	int i, j;

	if(MIN > MAX)
		return NULL;

	matrix* m = (matrix*)malloc(sizeof(matrix));
	m->col_size = col_size;
	m->row_size = row_size;
	m->mat = (double**)malloc(sizeof(double*) * m->row_size);

	for (i = 0; i < m->row_size; i++)
	{
		m->mat[i] = (double*)malloc(sizeof(double) * m->col_size);
		for (j = 0; j < m->col_size; j++)
			m->mat[i][j] = (double)rand() / RAND_MAX * (MAX - MIN) + MIN;
	}

	return m;
}

matrix* get_submatrix(const matrix* m, int x, int y, int dimension)
{
	int i, j, k;
	matrix* n = (matrix*)malloc(sizeof(matrix));
	n->row_size = n->col_size = dimension;
	n->mat = (double**)malloc(sizeof(double*)*dimension);

	for(i=0; i<dimension; i++)
	{
		n->mat[i] = (double*)malloc(sizeof(double)*dimension);
		for (j = 0; j < dimension; j++)
			n->mat[i][j] = m->mat[x + i][y + j];
	}
	return n;
}

matrix* brute_force_multiplication(const matrix* one, const matrix* two)
{
	int i, j, k;

	if (one->col_size != two->row_size)
		return NULL;

	matrix* prod = (matrix*)malloc(sizeof(matrix));
	prod->row_size = one->row_size;
	prod->col_size = two->col_size;
	prod->mat = (double**)malloc(sizeof(double*) * prod->row_size);

	for (i = 0; i < one->row_size; i++)
	{
		prod->mat[i] = (double*)malloc(sizeof(double)*two->col_size);
		for (k = 0; k < two->col_size; k++)
		{
			prod->mat[i][k] = 0;
			for (j = 0; j < one->col_size; j++)
			{
				prod->mat[i][k] += one->mat[i][j] * two->mat[j][k];
			}
		}
	}
	return prod;
}

matrix* divide_and_conquer_inner(const matrix* one, const matrix* two, int x1, int y1, int x2, int y2, int dimension)
{
	int i, j, k;
	matrix* prod = (matrix*)malloc(sizeof(matrix));
	prod->row_size = prod->col_size = dimension;
	prod->mat = (double**)malloc(sizeof(double*)*dimension);

	if (dimension <= 16)
	{
		double sum = 0;
		for (i = 0; i < dimension; i++)
		{
			prod->mat[i] = (double*)malloc(sizeof(double)*dimension);
			for (k = 0; k < dimension; k++)
			{
				sum = 0;
				for (j = 0; j < dimension; j++)
					sum += one->mat[x1 + i][y1 + j] * two->mat[x2 + j][y2 + k];
				prod->mat[i][k] = sum;
			}
		}
		return prod;
	}
	else
	{
		int d = dimension / 2;
		int free_index = 0;
		matrix** matrices = (matrix**)malloc(sizeof(matrix) * 8);

		matrices[0] = divide_and_conquer_inner(one, two, x1, y1, x2, y2, d);					//ae
		matrices[1] = divide_and_conquer_inner(one, two, x1, y1 + d, x2 + d, y2, d);			//bg
		matrices[2] = divide_and_conquer_inner(one, two, x1, y1, x2, y2 + d, d);				//af
		matrices[3] = divide_and_conquer_inner(one, two, x1, y1 + d, x2 + d, y2 + d, d);		//bh
		matrices[4] = divide_and_conquer_inner(one, two, x1 + d, y1, x2, y2, d);				//ce
		matrices[5] = divide_and_conquer_inner(one, two, x1 + d, y1 + d, x2 + d, y2, d);		//dg
		matrices[6] = divide_and_conquer_inner(one, two, x1 + d, y1, x2, y2 + d, d);			//cf
		matrices[7] = divide_and_conquer_inner(one, two, x1 + d, y1 + d, x2 + d, y2 + d, d);	//dh

		for (i = 0; i<dimension; i++)
			prod->mat[i] = (double*)malloc(sizeof(double)*dimension);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[i][j] = matrices[0]->mat[i][j] + matrices[1]->mat[i][j];				//ae + bg
			}
		}
		free_matrix(matrices[free_index++]);
		free_matrix(matrices[free_index++]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[i][d + j] = matrices[2]->mat[i][j] + matrices[3]->mat[i][j];			//af + bh
			}
		}
		free_matrix(matrices[free_index++]);
		free_matrix(matrices[free_index++]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[d + i][j] = matrices[4]->mat[i][j] + matrices[5]->mat[i][j];			//ce + dg
			}
		}
		free_matrix(matrices[free_index++]);
		free_matrix(matrices[free_index++]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[d + i][d + j] = matrices[6]->mat[i][j] + matrices[7]->mat[i][j];		//cf + dh
			}
		}
		free_matrix(matrices[free_index++]);
		free_matrix(matrices[free_index++]);

		free(matrices);

		return prod;
	}
}

matrix* divide_and_conquer(const matrix* one, const matrix* two)
{
	if (one->row_size != one->col_size || two->col_size != two->row_size)
		return NULL;
	if (one->row_size != two->row_size)
		return NULL;
	if (!(one->row_size>0 && ((one->row_size & (one->row_size - 1)) == 0)))
		return NULL;

	return divide_and_conquer_inner(one, two, 0, 0, 0, 0, one->row_size);
}

void* divide_and_conquer_parallel_inner(void* args)
{
	parallel_params* pp = (parallel_params*)args;
	return (void*)divide_and_conquer_inner(pp->one, pp->two, pp->x1, pp->y1, pp->x2, pp->y2, pp->dimension);
}

matrix* divide_and_conquer_parallel(const matrix* one, const matrix* two)
{
	int i, j;
	int free_index = 0;
	int dimension = one->row_size;
	int d = dimension / 2;
	int no_of_cores;
	SYSTEM_INFO sys_info;

	matrix* prod = (matrix*)malloc(sizeof(matrix));
	prod->row_size = prod->col_size = dimension;
	prod->mat = (double**)malloc(sizeof(double*)*dimension);

	matrix** matrices = (matrix**)malloc(sizeof(matrix) * 8);
	parallel_params* pp = (parallel_params*) malloc(sizeof(parallel_params) * 8);

	for(i=0; i<8; i++)
	{
		pp[i].one = one;
		pp[i].two = two;
		pp[i].dimension = d;
	}

	GetSystemInfo(&sys_info);
	no_of_cores = (int)sys_info.dwNumberOfProcessors > 8 ? 8 : (int)sys_info.dwNumberOfProcessors;
	pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t)*no_of_cores);
	void** temp_result = malloc(sizeof(void*)*8);

	pp[0].x1 = pp[0].y1 = pp[0].x2 = pp[0].y2 = 0;
	pp[1].x1 = 0, pp[1].y1 = d, pp[1].x2 = d, pp[1].y2 = 0;
	pp[2].x1 = 0, pp[2].y1 = 0, pp[2].x2 = 0, pp[2].y2 = d;
	pp[3].x1 = 0, pp[3].y1 = d, pp[3].x2 = d, pp[3].y2 = d;
	pp[4].x1 = d, pp[4].y1 = 0, pp[4].x2 = 0, pp[4].y2 = 0;
	pp[5].x1 = d, pp[5].y1 = d, pp[5].x2 = d, pp[5].y2 = 0;
	pp[6].x1 = d, pp[6].y1 = 0, pp[6].x2 = 0, pp[6].y2 = d;
	pp[7].x1 = pp[7].y1 = pp[7].x2 = pp[7].y2 = d;

	for(i=0; i<8/no_of_cores; i++)
	{
		for(j=0; j<no_of_cores; j++)
		{
			pthread_create(&threads[j], NULL, divide_and_conquer_parallel_inner, &pp[i*no_of_cores + j]);
		}
		for(j=0; j<no_of_cores; j++)
		{
			pthread_join(threads[j], &temp_result[i*no_of_cores + j]);
			matrices[i*no_of_cores + j] = (matrix*)temp_result[i*no_of_cores + j];
		}
	}

	free(pp);
	free(threads);

	for (i = 0; i<dimension; i++)
			prod->mat[i] = (double*)malloc(sizeof(double)*dimension);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[i][j] = matrices[0]->mat[i][j] + matrices[1]->mat[i][j];				//ae + bg
		}
	}

	free_matrix(matrices[free_index++]);
	free_matrix(matrices[free_index++]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[i][d + j] = matrices[2]->mat[i][j] + matrices[3]->mat[i][j];			//af + bh
		}
	}
	free_matrix(matrices[free_index++]);
	free_matrix(matrices[free_index++]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[d + i][j] = matrices[4]->mat[i][j] + matrices[5]->mat[i][j];			//ce + dg
		}
	}
	free_matrix(matrices[free_index++]);
	free_matrix(matrices[free_index++]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[d + i][d + j] = matrices[6]->mat[i][j] + matrices[7]->mat[i][j];		//cf + dh
		}
	}
	free_matrix(matrices[free_index++]);
	free_matrix(matrices[free_index++]);
	free(matrices);
	free(temp_result);

	return prod;
}

matrix* strassens_inner(const matrix* one, const matrix* two, int x1, int y1, int x2, int y2, int dimension)
{
	int i, j, k;
	matrix* prod = (matrix*)malloc(sizeof(matrix));
	prod->row_size = prod->col_size = dimension;
	prod->mat = (double**)malloc(sizeof(double*)*dimension);

	if (dimension <= 16)
	{
		double sum = 0;
		for (i = 0; i < dimension; i++)
		{
			prod->mat[i] = (double*)malloc(sizeof(double)*dimension);
			for (k = 0; k < dimension; k++)
			{
				sum = 0;
				for (j = 0; j < dimension; j++)
					sum += one->mat[x1 + i][y1 + j] * two->mat[x2 + j][y2 + k];
				prod->mat[i][k] = sum;
			}
		}
		return prod;
	}
	else
	{
		int d = dimension / 2;
		matrix** matrices = (matrix**)malloc(sizeof(matrix) * 7);

		matrix* temp1 = (matrix*)malloc(sizeof(matrix));
		temp1->row_size = temp1->col_size = dimension;
		temp1->mat = (double**)malloc(sizeof(double*)*dimension);

		matrix* temp2 = (matrix*)malloc(sizeof(matrix));
		temp2->row_size = temp2->col_size = dimension;
		temp2->mat = (double**)malloc(sizeof(double*)*dimension);


		for (i = 0; i < dimension; i++)
		{
			temp1->mat[i] = (double*)malloc(sizeof(double)*dimension);
			temp2->mat[i] = (double*)malloc(sizeof(double)*dimension);
		}

		//one:
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = one->mat[i + x1][j + y1] + one->mat[i + x1 + d][j + y1 + d];			//a + d
				temp2->mat[i][j] = two->mat[i + x2][j + y2] + two->mat[i + x2 + d][j + y2 + d];			//e + h
			}
		}
		matrices[0] = strassens_inner(temp1, temp2, 0, 0, 0, 0, d);										//(a + d)(e + h)

		//two:
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = one->mat[i + x1 + d][j + y1] + one->mat[i + x1 + d][j + y1 + d];		//c + d
			}
		}
		matrices[1] = strassens_inner(temp1, two, 0, 0, x2, y2, d);										//(c + d)e

		//three
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = two->mat[i + x2][j + y2 + d] - two->mat[i + x2 + d][j + y2 + d];		//f - h
			}
		}
		matrices[2] = strassens_inner(one, temp1, x1, y1, 0, 0, d);										//a(f - h)

		//four
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = two->mat[i + x2 + d][j + y2] - two->mat[i + x2][j + y2];				//g - e
			}
		}
		matrices[3] = strassens_inner(one, temp1, x1 + d, y1 + d, 0, 0, d);								//d(g - e)

		//five
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = one->mat[i + x1][j + y1] + one->mat[i + x1][j + y1 + d];				//a + b
			}
		}
		matrices[4] = strassens_inner(temp1, two, 0, 0, x2 + d, y2 + d, d);								//(a + b)h

		//six
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = one->mat[i + x1 + d][j + y1] - one->mat[i + x1][j + y1];				//c - a
				temp2->mat[i][j] = two->mat[i + x2][j + y2] + two->mat[i + x2][j + y2 + d];				//e + f
			}
		}
		matrices[5] = strassens_inner(temp1, temp2, 0, 0, 0, 0, d);										//(c - a)(e + f)

		//seven
		for(i = 0; i < d; i++)
		{
			for(j = 0; j < d; j++)
			{
				temp1->mat[i][j] = one->mat[i + x1][j + y1 + d] - one->mat[i + x1 + d][j + y1 + d];		//b - d
				temp2->mat[i][j] = two->mat[i + x2 + d][j + y2] + two->mat[i + x2 + d][j + y2 + d];		//g + h
			}
		}
		matrices[6] = strassens_inner(temp1, temp2, 0, 0, 0, 0, d);										//(b - d)(g + h)

		free_matrix(temp1);
		free_matrix(temp2);

		for (i = 0; i < dimension; i++)
			prod->mat[i] = (double*)malloc(sizeof(double)*dimension);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[i][j] = matrices[0]->mat[i][j] + matrices[3]->mat[i][j] - matrices[4]->mat[i][j] + matrices[6]->mat[i][j];			//0 + 3 - 4 + 6
			}
		}

		free_matrix(matrices[6]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[d + i][d + j] = matrices[0]->mat[i][j] - matrices[1]->mat[i][j] + matrices[2]->mat[i][j] + matrices[5]->mat[i][j];	//0 - 1 + 2 + 5
			}
		}
		free_matrix(matrices[0]);
		free_matrix(matrices[5]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[i][d + j] = matrices[2]->mat[i][j] + matrices[4]->mat[i][j];															//2 + 4
			}
		}
		free_matrix(matrices[4]);
		free_matrix(matrices[2]);

		for (i = 0; i < d; i++)
		{
			for (j = 0; j < d; j++)
			{
				prod->mat[d + i][j] = matrices[1]->mat[i][j] + matrices[3]->mat[i][j];															//1 + 3
			}
		}
		free_matrix(matrices[1]);
		free_matrix(matrices[3]);

		free(matrices);

		return prod;
	}
}

matrix* strassens(const matrix* one, const matrix* two)
{
	if (one->row_size != one->col_size || two->col_size != two->row_size)
		return NULL;
	if (one->row_size != two->row_size)
		return NULL;
	if (!(one->row_size>0 && ((one->row_size & (one->row_size - 1)) == 0)))
		return NULL;

	return strassens_inner(one, two, 0, 0, 0, 0, one->row_size);
}

void* strassens_parallel_inner(void* args)
{
	parallel_params* pp = (parallel_params*)args;
	return (void*)strassens_inner(pp->one, pp->two, pp->x1, pp->y1, pp->x2, pp->y2, pp->dimension);
}

matrix* strassens_parallel(matrix* one, matrix* two)
{
	int i, j;
	int dimension = one->row_size;
	int d = dimension / 2;

	matrix* prod = (matrix*)malloc(sizeof(matrix));
	prod->row_size = prod->col_size = dimension;
	prod->mat = (double**)malloc(sizeof(double*)*dimension);

	matrix** matrices = (matrix**)malloc(sizeof(matrix) * 7);
	parallel_params_mutable* pp = (parallel_params_mutable*) malloc(sizeof(parallel_params_mutable) * 7);

	for(i=0; i<7; i++)
	{
		pp[i].dimension = d;
		pp[i].one = (matrix*)malloc(sizeof(matrix));
		pp[i].two = (matrix*)malloc(sizeof(matrix));

		if(!(i == 2 || i == 3))
		{
			pp[i].one->row_size = pp[i].one->col_size = dimension;
			pp[i].one->mat = (double**)malloc(sizeof(double*)*dimension);
			for (j=0; j<dimension; j++)
			{
				pp[i].one->mat[j] = (double*)malloc(sizeof(double)*dimension);
			}
		}
		if(!(i == 1 || i == 4))
		{
			pp[i].two->row_size = pp[i].two->col_size = dimension;
			pp[i].two->mat = (double**)malloc(sizeof(double*)*dimension);
			for (j=0; j<dimension; j++)
			{
				pp[i].two->mat[j] = (double*)malloc(sizeof(double)*dimension);
			}
		}
	}

	//one:
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[0].one->mat[i][j] = one->mat[i][j] + one->mat[i + d][j + d];			//a + d
			pp[0].two->mat[i][j] = two->mat[i][j] + two->mat[i + d][j + d];			//e + h
		}
	}
	pp[0].x1 = pp[0].y1 = 0;
	pp[0].x2 = pp[0].y2 = 0;

	//two:
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[1].one->mat[i][j] = one->mat[i + d][j] + one->mat[i + d][j + d];		//c + d
		}
	}
	pp[1].two = two;
	pp[1].x1 = 0, pp[1].y1 = 0;
	pp[1].x2 = 0, pp[1].y2 = 0;

	//three
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[2].two->mat[i][j] = two->mat[i][j + d] - two->mat[i + d][j + d];		//f - h
		}
	}
	pp[2].one = one;
	pp[2].x1 = 0, pp[2].y1 = 0;
	pp[2].x2 = 0, pp[2].y2 = 0;

	//four
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[3].two->mat[i][j] = two->mat[i + d][j] - two->mat[i][j];				//g - e
		}
	}

	pp[3].one = one;
	pp[3].x1 = d, pp[3].y1 = d;
	pp[3].x2 = 0, pp[3].y2 = 0;

	//five
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[4].one->mat[i][j] = one->mat[i][j] + one->mat[i][j + d];				//a + b
		}
	}
	pp[4].two = two;
	pp[4].x1 = 0, pp[4].y1 = 0;
	pp[4].x2 = d, pp[4].y2 = d;

	//six
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[5].one->mat[i][j] = one->mat[i + d][j] - one->mat[i][j];				//c - a
			pp[5].two->mat[i][j] = two->mat[i][j] + two->mat[i][j + d];				//e + f
		}
	}
	pp[5].x1 = 0, pp[5].y1 = 0;
	pp[5].x2 = 0, pp[5].y2 = 0;

	//seven
	for(i = 0; i < d; i++)
	{
		for(j = 0; j < d; j++)
		{
			pp[6].one->mat[i][j] = one->mat[i][j + d] - one->mat[i + d][j + d];		//b - d
			pp[6].two->mat[i][j] = two->mat[i + d][j] + two->mat[i + d][j + d];		//g + h
		}
	}
	pp[6].x1 = 0, pp[6].y1 = 0;
	pp[6].x2 = 0, pp[6].y2 = 0;


	pthread_t* threads = (pthread_t*)malloc(sizeof(pthread_t) * 7);
	void** temp_result = malloc(sizeof(void*) * 7);

	for(j=0; j<7; j++)
	{
		pthread_create(&threads[j], NULL, strassens_parallel_inner, &pp[j]);
	}
	for(j=0; j<7; j++)
	{
		pthread_join(threads[j], &temp_result[j]);
		matrices[j] = (matrix*)temp_result[j];
	}

	for(i=0; i<7; i++)
	{
		if(!(i == 2 || i == 3))
			free_matrix(pp[i].one);
		if(!(i == 1 || i == 4))
			free_matrix(pp[i].two);
	}
	free(pp);

	for (i = 0; i < dimension; i++)
		prod->mat[i] = (double*)malloc(sizeof(double)*dimension);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[i][j] = matrices[0]->mat[i][j] + matrices[3]->mat[i][j] - matrices[4]->mat[i][j] + matrices[6]->mat[i][j];			//0 + 3 - 4 + 6
		}
	}

	free_matrix(matrices[6]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[d + i][d + j] = matrices[0]->mat[i][j] - matrices[1]->mat[i][j] + matrices[2]->mat[i][j] + matrices[5]->mat[i][j];	//0 - 1 + 2 + 5
		}
	}
	free_matrix(matrices[0]);
	free_matrix(matrices[5]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[i][d + j] = matrices[2]->mat[i][j] + matrices[4]->mat[i][j];															//2 + 4
		}
	}
	free_matrix(matrices[4]);
	free_matrix(matrices[2]);

	for (i = 0; i < d; i++)
	{
		for (j = 0; j < d; j++)
		{
			prod->mat[d + i][j] = matrices[1]->mat[i][j] + matrices[3]->mat[i][j];															//1 + 3
		}
	}
	free_matrix(matrices[1]);
	free_matrix(matrices[3]);

	free(matrices);

	return prod;
}