/*
    util.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"

char *strconcat(const char *s, const char *t){
    static char z[BUFSIZ];
    strcpy(z, s);
    strcat(z, t);
    return z;
}

void normalize_matrix_row(double **dst, double **src, int rows, int cols){
    int i, j;
    double z;
    
    for(i = 0;i < rows;i++){
        for(j = 0, z = 0;j < cols;j++)
            z += src[i][j];
        for(j = 0;j < cols;j++)
            dst[i][j] = src[i][j] / z;
    }
}

void normalize_matrix_col(double **dst, double **src, int rows, int cols){
    int i, j;
    double z;
    
    for(j = 0;j < cols;j++){
        for(i = 0, z = 0;i < rows;i++)
            z += src[i][j];
        for(i = 0;i < rows;i++)
            dst[i][j] = src[i][j] / z;
    }
}


double logsumexp(double *nums, size_t ct){
    double max_exp = nums[0], sum = 0.0;
    size_t i;
    for (i = 1 ; i < ct ; i++)
        if (nums[i] > max_exp)
            max_exp = nums[i];

    for (i = 0; i < ct ; i++)
        sum += exp(nums[i] - max_exp);

    return log(sum) + max_exp;
}