#pragma once

typedef enum bool{false, true} bool;

typedef struct matrix
{
	int row_size;
	int col_size;
	double** mat;
}matrix;

typedef struct parallel_params
{
	const matrix* one;
	const matrix* two;
	int x1;
	int y1;
	int x2;
	int y2;
	int dimension;
}parallel_params;

typedef struct parallel_params_mutable
{
	matrix* one;
	matrix* two;
	int x1;
	int y1;
	int x2;
	int y2;
	int dimension;
}parallel_params_mutable;

matrix* read_matrix(int row_size, int col_size);
void display_matrix(const matrix* m);
void output_to_file(const matrix* m, const char* path);
void free_matrix(matrix* m);
matrix* generate_random_matrix(int row_size, int col_size, const int MIN, const int MAX);
matrix* get_submatrix(const matrix* m, int x, int y, int dimension);
matrix* brute_force_multiplication(const matrix* one, const matrix* two);
matrix* brute_force_multiplication_inner(const matrix* a, const matrix* b, int x1, int y1, int x2, int y2, int dimension);
matrix* divide_and_conquer_inner(const matrix* one, const matrix* two, int x1, int y1, int x2, int y2, int dimension);
matrix* divide_and_conquer(const matrix* one, const matrix* two);
void* divide_and_conquer_parallel_inner(void* args);
matrix* divide_and_conquer_parallel(const matrix* one, const matrix* two);
matrix* strassens_inner(const matrix* one, const matrix* two, int x1, int y1, int x2, int y2, int dimension);
matrix* strassens(const matrix* one, const matrix* two);
void* strassens_parallel_inner(void* args);
matrix* strassens_parallel(matrix* one, matrix* two);
void caller(matrix* (funct)(const matrix*, const matrix*), char* text, matrix* m1, matrix* m2);


// matrix* brute_force_multiplication_with_simd(const matrix* a, const matrix* b);
