/*
    writer.h
    a header file of matrix writer
*/
#ifndef WRITER_H
#define WRITER_H
#include <stdio.h>

extern void mvslda_write(FILE *pp, FILE *tp, FILE *ep,
        double **phi, double **theta, double **eta,
        int nclass, int nlex, int nresp, int ndoc);
void write_matrix(FILE *fp, double **matrix, int rows, int cols);
void write_imatrix(FILE *fp, int **matrix, int rows, int cols);
void write_imatrix_transpose(FILE *fp, int **matrix, int rows, int cols);
#endif
