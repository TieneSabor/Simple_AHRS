#include "Quaternion_Gauss_Newton.h"

void get_rpy(TYPE** rpy, TYPE** q){
    TYPE q0 = q[0][0];
    TYPE q1 = q[1][0];
    TYPE q2 = q[2][0];
    TYPE q3 = q[3][0];
    rpy[0][0] = atan2(2*(q0*q1+q2*q3),1-2*(q1*q1+q2*q2));
    rpy[1][0] = asin(2*(q0*q2-q3*q1));
    rpy[2][0] = atan2(2*(q0*q3+q1*q2),1-2*(q2*q2+q3*q3));
}

#ifdef NINE_AXIS

TYPE** pos_def;
TYPE** pos_def_inv;
TYPE** J;
TYPE** JT;
TYPE** search_dir;
TYPE** search_vec;
TYPE** search_vec_pad;
TYPE** err;
TYPE** I;
TYPE** L;
TYPE** LT;
TYPE** X;

void initialize_Gauss_Newton_9axis(void){ 
  pos_def = dmatrix(4,4);
  pos_def_inv = dmatrix(4,4);
  J = dmatrix(6,4);
  JT = dmatrix(4,6);
  search_dir = dmatrix(4,6);
  search_vec = dmatrix(4,1);
  search_vec_pad = dmatrix(6,1);
  err = dmatrix(6,1);
  I = dmatrix(4,4);
  L = dmatrix(4,4);
  LT = dmatrix(4,4);
  X = dmatrix(4,4);
	makeeyes(I, 4);
}

void get_Jacobian_9axis(TYPE** J, TYPE **qk, TYPE **bmk, TYPE **bak){
    TYPE q1 = qk[0][0];
    TYPE q2 = qk[1][0];
    TYPE q3 = qk[2][0];
    TYPE q4 = qk[3][0];

    TYPE bmx = bmk[0][0];
    TYPE bmy = bmk[1][0];
    TYPE bmz = bmk[2][0];
    TYPE bax = bak[0][0];
    TYPE bay = bak[1][0];
    TYPE baz = bak[2][0];

    J[0][0] = -2*(q1*bmx-q4*bmy+q3*bmz);
    J[1][0] = -2*(q4*bmx+q1*bmy-q2*bmz);
    J[2][0] = -2*(-q3*bmx+q2*bmy+q1*bmz);
    J[3][0] = -2*(q1*bax-q4*bay+q3*baz);
    J[4][0] = -2*(q4*bax+q1*bay-q2*baz);
    J[5][0] = -2*(-q3*bax+q2*bay+q1*baz);

    J[0][1] = -2*(q2*bmx+q3*bmy+q4*bmz);
    J[1][1] = -2*(q3*bmx-q2*bmy-q1*bmz);
    J[2][1] = -2*(q4*bmx+q1*bmy-q2*bmz);
    J[3][1] = -2*(q2*bax+q3*bay+q4*baz);
    J[4][1] = -2*(q3*bax-q2*bay-q1*baz);
    J[5][1] = -2*(q4*bax+q1*bay-q2*baz);

    J[0][2] = -2*(-q3*bmx+q2*bmy+q1*bmz);
    J[1][2] = -2*(q2*bmx+q3*bmy+q4*bmz);
    J[2][2] = -2*(-q1*bmx+q4*bmy-q3*bmz);
    J[3][2] = -2*(-q3*bax+q2*bay+q1*baz);
    J[4][2] = -2*(q2*bax+q3*bay+q4*baz);
    J[5][2] = -2*(-q1*bax+q4*bay-q3*baz);

    J[0][3] = -2*(-q4*bmx-q1*bmy+q2*bmz);
    J[1][3] = -2*(q1*bmx-q4*bmy+q3*bmz);
    J[2][3] = -2*(q2*bmx+q3*bmy+q4*bmz);
    J[3][3] = -2*(-q4*bax-q1*bay+q2*baz);
    J[4][3] = -2*(q1*bax-q4*bay+q3*baz);
    J[5][3] = -2*(q2*bax+q3*bay+q4*baz);
}

void get_error_9axis(TYPE** err, TYPE** qk, TYPE** eak, TYPE** bak, TYPE** emk, TYPE** bmk){
    TYPE bmx = bmk[0][0];
    TYPE bmy = bmk[1][0];
    TYPE bmz = bmk[2][0];
    TYPE bax = bak[0][0];
    TYPE bay = bak[1][0];
    TYPE baz = bak[2][0];

    TYPE emx = emk[0][0];
    TYPE emy = emk[1][0];
    TYPE emz = emk[2][0];
    TYPE eax = eak[0][0];
    TYPE eay = eak[1][0];
    TYPE eaz = eak[2][0];

    TYPE q1 = qk[0][0];
    TYPE q2 = qk[1][0];
    TYPE q3 = qk[2][0];
    TYPE q4 = qk[3][0];

    TYPE R11 = pow(q1,2)+pow(q2,2)-pow(q3,2)-pow(q4,2);
    TYPE R12 = 2*(q2*q3-q1*q4);
    TYPE R13 = 2*(q2*q4+q1*q3);
    TYPE R21 = 2*(q2*q3+q1*q4);
    TYPE R22 = pow(q1,2)+pow(q3,2)-pow(q2,2)-pow(q4,2);
    TYPE R23 = 2*(q3*q4-q1*q2);
    TYPE R31 = 2*(q2*q4-q1*q3);
    TYPE R32 = 2*(q3*q4+q1*q2);
    TYPE R33 = pow(q1,2)+pow(q4,2)-pow(q3,2)-pow(q2,2);

    err[0][0] = emx-R11*bmx-R12*bmy-R13*bmz;
    err[1][0] = emy-R21*bmx-R22*bmy-R23*bmz;
    err[2][0] = emz-R31*bmx-R32*bmy-R33*bmz;

    err[3][0] = eax-R11*bax-R12*bay-R13*baz;
    err[4][0] = eay-R21*bax-R22*bay-R23*baz;
    err[5][0] = eaz-R31*bax-R32*bay-R33*baz;
}

void update_quaternion_9axis(TYPE** newq, TYPE** qk, TYPE** eak, TYPE** bak, TYPE** emk, TYPE** bmk){
    get_Jacobian_9axis(J, qk, bmk, bak);
    transpose(J, JT, 6, 4);
    matmult(JT, J, pos_def, 4, 6, 4);
    Cholesky_Inv(pos_def, X, pos_def_inv, L, LT, I, 4);
    TYPE** sbi = dmatrix(4,4);
    matmult(pos_def, pos_def_inv, sbi, 4, 4, 4);
    matmult(pos_def_inv, JT, search_dir, 4, 4, 6);
    get_error_9axis(err, qk, eak, bak, emk, bmk);
    matmult(search_dir, err, search_vec, 4,6,1);
    matsub(qk, search_vec, newq, 4, 1);
		/*
		freedmat(pos_def,4,4);
		freedmat(J,6,4);
		freedmat(search_dir,4,6);
		freedmat(search_vec,4,1);
		freedmat(err,6,1);
		freedmat(JT,4,6);
		freedmat(pos_def_inv,4,4);
		freedmat(sbi,4,4);*/
}

#elif defined SIX_AXIS

TYPE** pos_def;
TYPE** pos_def_inv;
TYPE** J;
TYPE** JT;
TYPE** search_dir;
TYPE** search_vec;
TYPE** search_vec_pad;
TYPE** err;
TYPE** I;
TYPE** L;
TYPE** LT;
TYPE** X;

void initialize_Gauss_Newton_6axis(void){
		pos_def = dmatrix(3,3);
    pos_def_inv = dmatrix(3,3);
    J = dmatrix(3,3);
    JT = dmatrix(3,3);
    search_dir = dmatrix(3,3);
    search_vec = dmatrix(3,1);
    search_vec_pad = dmatrix(4,1);
    err = dmatrix(3,1);
    I = dmatrix(3,3);
    L = dmatrix(3,3);
    LT = dmatrix(3,3);
    X = dmatrix(3,3);
		makeeyes(I, 3);
}

void get_Jacobian_6axis(TYPE** J, TYPE **qk, TYPE **bak){
    TYPE q1 = qk[0][0];
    TYPE q2 = qk[1][0];
    TYPE q3 = qk[2][0];
    TYPE q4 = qk[3][0];
    
    TYPE bax = bak[0][0];
    TYPE bay = bak[1][0];
    TYPE baz = bak[2][0];
    
    J[0][0] = -2*(q1*bax-q4*bay+q3*baz);
    J[1][0] = -2*(q4*bax+q1*bay-q2*baz);
    J[2][0] = -2*(-q3*bax+q2*bay+q1*baz);
    
    J[0][1] = -2*(q2*bax+q3*bay+q4*baz);
    J[1][1] = -2*(q3*bax-q2*bay-q1*baz);
    J[2][1] = -2*(q4*bax+q1*bay-q2*baz);
    
    J[0][2] = -2*(-q3*bax+q2*bay+q1*baz);
    J[1][2] = -2*(q2*bax+q3*bay+q4*baz);
    J[2][2] = -2*(-q1*bax+q4*bay-q3*baz);
}

void get_error_6axis(TYPE** err, TYPE** qk, TYPE** eak, TYPE** bak){
    TYPE bax = bak[0][0];
    TYPE bay = bak[1][0];
    TYPE baz = bak[2][0];

    TYPE eax = eak[0][0];
    TYPE eay = eak[1][0];
    TYPE eaz = eak[2][0];

    TYPE q1 = qk[0][0];
    TYPE q2 = qk[1][0];
    TYPE q3 = qk[2][0];
    TYPE q4 = 0;//qk[3][0];

    TYPE R11 = pow(q1,2)+pow(q2,2)-pow(q3,2)-pow(q4,2);
    TYPE R12 = 2*(q2*q3-q1*q4);
    TYPE R13 = 2*(q2*q4+q1*q3);
    TYPE R21 = 2*(q2*q3+q1*q4);
    TYPE R22 = pow(q1,2)+pow(q3,2)-pow(q2,2)-pow(q4,2);
    TYPE R23 = 2*(q3*q4-q1*q2);
    TYPE R31 = 2*(q2*q4-q1*q3);
    TYPE R32 = 2*(q3*q4+q1*q2);
    TYPE R33 = pow(q1,2)+pow(q4,2)-pow(q3,2)-pow(q2,2);

    err[0][0] = eax-R11*bax-R12*bay-R13*baz;
    err[1][0] = eay-R21*bax-R22*bay-R23*baz;
    err[2][0] = eaz-R31*bax-R32*bay-R33*baz;
}

void update_quaternion_6axis(TYPE** newq, TYPE** qk, TYPE** eak, TYPE** bak){
    // Get the 3*4 Jacobian matrix of current (q, imu_measurement)
    get_Jacobian_6axis(J, qk, bak);
    // Get pos_def = J*JT
    transpose(J, JT, 3, 3);
    matmult(JT, J, pos_def, 3, 3, 3);
		// Get pos_def^-1
    Cholesky_Inv(pos_def, X, pos_def_inv, L, LT, I, 3);
    // search_vector = (J*JT)^(-1)*JT*err
    matmult(pos_def_inv, JT, search_dir, 3, 3, 3);
    get_error_6axis(err, qk, eak, bak);
    matmult(search_dir, err, search_vec, 3,3,1);

    search_vec_pad[0][0] = search_vec[0][0];
    search_vec_pad[1][0] = search_vec[1][0];
    search_vec_pad[2][0] = search_vec[2][0];
    search_vec_pad[3][0] = 0;
    // q_k+1 = q_k + search_vector
    matsub(qk, search_vec_pad, newq, 4, 1);
}

#endif
