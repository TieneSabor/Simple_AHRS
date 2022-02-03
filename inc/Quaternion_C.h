#ifndef Qtn_C_H
#define Qtn_C_H

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "basic_linear.h"

#define USE_BASIC_LINEAR

//#ifdef BASIC_LINEAR
//#endif

typedef struct Quaternion{
    TYPE w,x,y,z;
}Qtn;

//TYPE test;

void print_Qtn(Qtn* q);

void from_Rod(Qtn* result, TYPE theta, TYPE x, TYPE y, TYPE z);
TYPE get_theta(Qtn* q);
void add(Qtn* q1, Qtn* q2, Qtn* q3);
void sub(Qtn* q1, Qtn* q2, Qtn* q3);
void mul(Qtn* q1, Qtn* q2, Qtn* q3);
void scal(Qtn* q1, TYPE b, Qtn* q3);
void inv(Qtn* q, Qtn* qt);
TYPE norm_2(Qtn* q);
void normalize(Qtn* q, Qtn* qn);

#ifdef USE_BASIC_LINEAR

TYPE** Qtn_to_dmatrix(Qtn* q);

Qtn* dmatrix_to_Qtn(TYPE** d);

#endif

#endif
