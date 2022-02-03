#include "Quaternion_C.h"

void print_Qtn(Qtn* q){
    printf("w: %lf, x: %lf, y: %lf, z: %lf \r\n", q->w, q->x, q->y, q->z);
}

//(x,y,z) is the rotation axis and theta is the angle
void from_Rod(Qtn* result, TYPE theta, TYPE x, TYPE y, TYPE z){
    TYPE norm = sqrt(x*x+y*y+z*z);
    x = x/norm;
    y = y/norm;
    z = z/norm;
    result->w = cos(theta/2);
    result->x = sin(theta/2)*x;
    result->y = sin(theta/2)*y;
    result->z = sin(theta/2)*z;
}

//theta is ranged from -PI to PI
TYPE get_theta(Qtn* q){
    return acos(fabs(q->w))*2;
}

void add(Qtn* q1, Qtn* q2, Qtn* q3){
    q3->w = q1->w+q2->w;
    q3->x = q1->x+q2->x;
    q3->y = q1->y+q2->y;
    q3->z = q1->z+q2->z;
}

void sub(Qtn* q1, Qtn* q2, Qtn* res){
    res->w = q1->w-q2->w;
    res->x = q1->x-q2->x;
    res->y = q1->y-q2->y;
    res->z = q1->z-q2->z;
}

void mul(Qtn* q1, Qtn* q2, Qtn* res){
    res->w = q1->w*q2->w - q1->x*q2->x - q1->y*q2->y - q1->z*q2->z;
    res->x = q1->w*q2->x + q1->x*q2->w + q1->y*q2->z - q1->z*q2->y;
    res->y = q1->w*q2->y + q1->y*q2->w + q1->z*q2->x - q1->x*q2->z;
    res->z = q1->w*q2->z + q1->z*q2->w + q1->x*q2->y - q1->y*q2->x;
}

void scal(Qtn* q1, TYPE b, Qtn* res){
    res->w = q1->w*b;
    res->x = q1->x*b;
    res->y = q1->y*b;
    res->z = q1->z*b;
}

void inv(Qtn* q, Qtn* res){
    res->w = q->w;
    res->x = -q->x;
    res->y = -q->y;
    res->z = -q->z;
}

TYPE norm_2(Qtn* q){
    return sqrt(pow(q->w,2)+pow(q->x,2)+pow(q->y,2)+pow(q->z,2));
}

void normalize(Qtn* q, Qtn* res){
    TYPE norm = norm_2(q);
    res->w = q->w/norm;
    res->x = q->x/norm;
    res->y = q->y/norm;
    res->z = q->z/norm;
}

#ifdef USE_BASIC_LINEAR

TYPE** Qtn_to_dmatrix(Qtn* q){
    TYPE** res = dmatrix(4,1);
    res[0][0] = q->w;
    res[1][0] = q->x;
    res[2][0] = q->y;
    res[3][0] = q->z;
    return res;
}

Qtn* dmatrix_to_Qtn(TYPE** d){
    Qtn* res;
    res->w = d[0][0];
    res->x = d[1][0];
    res->y = d[2][0];
    res->z = d[3][0];
    return res;
}

#endif
