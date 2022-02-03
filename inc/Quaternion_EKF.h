#ifndef Qtn_EKF_H
#define Qtn_EKF_H

//#include "Quaternion_Gauss_Newton.h"
#include "Quaternion_C.h"
#include "basic_linear.h"

#define H_dim 7

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct Filter{
    TYPE** x;
    TYPE** w;
    TYPE** P;
    TYPE** F;
    TYPE** H;
    TYPE** Q;// Conti. Time Q
    TYPE** R;// Discrete Time R
    TYPE* tau;
}Ft;

/* w is an array of angular velocity: {wx, wy, wz}
 * q is an parameter array of covariance matrix of transition process: {Q1_1, Q2_2....Q10_10}
 * r is an parameter array of covariance matrix of observation process: {R1_1, R2_2...R10_10}
 * tau is an array of time constant of gyroscope dynamics: {tau_rx, tau_ry, tau_rz}
 */
void Initialize_EKF_state(Ft* K, TYPE* w, TYPE* q, TYPE* r, TYPE* tau);

/* w is an array of angular velocity: {wx, wy, wz}
 * dt is the second between prediction steps
 * The explicit integral (Eular integral) was used.
 * We can see if it is necessary to make it more accurate (maybe R-K)
 * Or we can use a inplicit solution: 
 * d/dt[X; Z] = [F, Q; 0, -F^T]*[X; Z]; P*Z = X; X0=P0; Z0=I
 * Then solve it with Laplace or sth...QQ
 * Please refer to " Introduction to Random Signals and Applied Kalman Filtering", pg290~
 */ 
void predict_EKF_state(Ft* K, TYPE* w, TYPE dt);

/* y is a 7*1 matrix of observed quaternion: {wx, wy, wz, qw, qx, qy, qz}
 * This step is also called "update"
 */
void estimate_EKF_state(Ft* K, TYPE** y);

#ifdef __cplusplus
 }
#endif

#endif
