#ifndef SEQ_MATRIX_MUL_H
#define SEQ_MATRIX_MUL_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <papi.h>
#include <cblas.h>

#include "papi_timer.h"
#define NUM_EXP (2)

#define NUM_EVENTS (3)


//#define VERBOSE


/*
 * calculate matrix multiplication using order: ijk
 * input:
 *  matrix C, A, B, all in flat format
 * output:
 *  matrix C
 */
void seq_cal_ijk(double *C, double *A, double *B, int n);

#endif
