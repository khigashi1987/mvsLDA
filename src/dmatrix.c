/*
    dmatrix.c
    an implementation of double matrix
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dmatrix.h"
#define BUFFSIZE 65536

double **dmatrix(int rows, int cols){
    double **matrix;
    int i;
    
    matrix = (double **)calloc(rows, sizeof(double *));
    if(matrix == NULL)
        return NULL;
    for(i = 0;i < rows;i++){
        matrix[i] = (double *)calloc(cols, sizeof(double));
        if(matrix[i] == NULL)
            return NULL;
    }
    
    return matrix;
}

void free_dmatrix(double **matrix, int rows){
    int i;
    for(i = 0;i < rows;i++){
        free(matrix[i]);
    }
    free(matrix);
}

double **load_dmatrix(char *filename, int ndoc, int nresp){
    FILE *fp;
    char line[BUFFSIZE];
    char *tok;
    int i, j;
    double **resp;

    if((fp = fopen(filename, "r")) == NULL)
        return NULL;
    
    resp = dmatrix(ndoc, nresp);

    i = 0;
    while(fgets(line, sizeof(line), fp)){
        j = 0;
        char field[256];
        int offset = 0;
        while (sscanf(line + offset, "%255[^ \t\n]", field) == 1){
            if(j == nresp)
                break;
            resp[i][j] = atof(field);
            j++;
            offset += strlen(field);
            if (line[offset] != '\0')
                offset++;
        }
        i++;
        if(i > ndoc){
            fprintf(stderr, "dmatrix: n_responses must be identical to ndoc:\n");
            exit(1);
        }
    }

    fclose(fp);
    return resp;
}