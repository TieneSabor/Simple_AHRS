#include "Quaternion_EKF.h"
/*
TYPE x[10][1];
TYPE w[3][1];
TYPE P[10][10];
TYPE F[10][10];
TYPE H[H_dim][10];
TYPE Q[10][10];
TYPE R[4][4];

TYPE dx[10][1];
TYPE dP[10][10];
TYPE dP1[10][10];
TYPE dP2[10][10];
TYPE FT[10][10];

TYPE Kg[10][H_dim];
TYPE Kg1a[H_dim][H_dim];
TYPE Kg1aT[H_dim][H_dim];
TYPE Kg1b[H_dim][10];
TYPE Kg1c[H_dim][H_dim];
TYPE HT[H_dim][H_dim];

TYPE X[H_dim][H_dim];
TYPE L[H_dim][H_dim];
TYPE LT[H_dim][H_dim];
TYPE I[H_dim][H_dim];

TYPE error[H_dim][1];
TYPE fix[10][1];

TYPE P1[10][10];
TYPE P1T[10][10];
TYPE P2[10][10];
TYPE P3[10][10];
*/
TYPE** x;
TYPE** w;
TYPE** P;
TYPE** F;
TYPE** H;
TYPE** Q;
TYPE** R;

TYPE** dx;
TYPE** dP;
TYPE** dP1;
TYPE** dP2;
TYPE** FT;

TYPE** Kg;
TYPE** KgT;
TYPE** Kg1a;
TYPE** Kg1aT;
TYPE** Kg1b;
TYPE** Kg1c;
TYPE** PHT;
TYPE** HT;

TYPE** KMX;
TYPE** KML;
TYPE** KMLT;
TYPE** KMI;

TYPE** error;
TYPE** fix;

TYPE** P1;
TYPE** P1T;
TYPE** P2;
TYPE** P3;
TYPE** KgR;

void Initialize_EKF_state(Ft* K, TYPE* ini_w, TYPE* q, TYPE* r, TYPE* tau){
    // Initialize the tau
    K->tau = tau;

    // Initialize the state x
    x = dmatrix(10,1);
    K->x = x;
    K->x[0][0] = ini_w[0];
    K->x[1][0] = ini_w[1];
    K->x[2][0] = ini_w[2];
    K->x[8][0] = 1;

    // Initialize the angular velocity w
    w = dmatrix(3,1);
    K->w = w;
    K->w[0][0] = ini_w[0];
    K->w[1][0] = ini_w[1];
    K->w[2][0] = ini_w[2];

    // Initialize the state covariance matrix
    P = dmatrix(10,10);
    K->P = P;
    makeeyes(K->P, 10);

    // Initialize linearlize transition matrix F
    F = dmatrix(10,10);
    K->F = F;
    makezero(K->F, 10, 10);
    K->F[0][0] = -1/tau[0];
    K->F[1][0] = -1/tau[1];
    K->F[2][0] = -1/tau[2];
    K->F[0][3] = -1/tau[0];
    K->F[1][4] = -1/tau[1];
    K->F[2][5] = -1/tau[2];
    // q1 part
    //K->F[6][0] = -0.5*q2;
    //K->F[6][1] = -0.5*q3;
    //K->F[6][2] = -0.5*q4;
    K->F[6][7] = -0.5*ini_w[0];
    K->F[6][8] = -0.5*ini_w[1];
    K->F[6][9] = -0.5*ini_w[2];
    // q2 part
    K->F[7][0] = 0.5*1;//q1
    //K->F[7][1] = 0.5*q4;
    //K->F[7][2] = -0.5*q3;
    K->F[7][6] = 0.5*ini_w[0];
    K->F[7][8] = -0.5*ini_w[2];
    K->F[7][9] = 0.5*ini_w[1];
    // q3 part
    //K->F[8][0] = -0.5*q4;//q1
    K->F[8][1] = 0.5*1;
    //K->F[8][2] = 0.5*q2;
    K->F[8][6] = 0.5*ini_w[1];
    K->F[8][7] = 0.5*ini_w[2];
    K->F[8][9] = -0.5*ini_w[0];
    // q4 part
    //K->F[9][0] = 0.5*q3;//q1
    //K->F[9][1] = -0.5*q2;
    K->F[9][2] = 0.5*1;
    K->F[9][6] = 0.5*ini_w[2];
    K->F[9][7] = -0.5*ini_w[1];
    K->F[9][8] = 0.5*ini_w[0];

    // Initialize the linearlize observation matrix H
    H = dmatrix(H_dim,10);
    K->H = H;
    makezero(K->H, H_dim, 10);
    
    K->H[0][0] = 1;
    K->H[1][1] = 1;
    K->H[2][2] = 1;
    K->H[0][3] = 1;
    K->H[1][4] = 1;
    K->H[2][5] = 1;
    
    K->H[3][6] = 1;
    K->H[4][7] = 1;
    K->H[5][8] = 1;
    K->H[6][9] = 1;

    // Initialize the covariance matrix of transition process Q
    Q = dmatrix(10,10);
    K->Q = Q;
    makezero(K->Q, 10, 10);
    for(int i=0;i<10;i++){
        K->Q[i][i] = q[i];
    }

    // Initialize the covariance matrix of observation process R
    R = dmatrix(H_dim,H_dim);
    K->R = R;
    makezero(K->R, H_dim, H_dim);
    for(int i=0;i<H_dim;i++){
        K->R[i][i] = r[i];
    }

    dx = dmatrix(10,1);
    dP = dmatrix(10,10);
    dP1 = dmatrix(10,10);
    dP2 = dmatrix(10,10);
    FT = dmatrix(10,10);

    Kg = dmatrix(10,H_dim);
    KgT = dmatrix(H_dim,10);
    PHT = dmatrix(10,H_dim);
    Kg1a = dmatrix(H_dim,H_dim);
    Kg1aT = dmatrix(H_dim,H_dim);
    Kg1b = dmatrix(H_dim,10);
    Kg1c = dmatrix(H_dim,H_dim);
    HT = dmatrix(10,H_dim);

    KMX = dmatrix(H_dim,H_dim);
    KML = dmatrix(H_dim,H_dim);
    KMLT = dmatrix(H_dim,H_dim);
    KMI = dmatrix(H_dim,H_dim);

    error = dmatrix(H_dim,1);
    fix = dmatrix(10,1);

    P1 = dmatrix(10,10);
    P1T = dmatrix(10,10);
    P2 = dmatrix(10,10);
    P3 = dmatrix(10,10);
    KgR = dmatrix(10,10);
}

void predict_EKF_state(Ft* K, TYPE* w, TYPE dt){
    TYPE q1 = K->x[6][1];
    TYPE q2 = K->x[7][1];
    TYPE q3 = K->x[8][1];
    TYPE q4 = K->x[9][1];
    
    // Predict the state
    makezero(dx,10,1);
    TYPE* tau = K->tau;
    dx[0][0] = (-1/tau[0])*K->x[0][0]+(-1/tau[0])*K->x[3][0]+(1/tau[0])*w[0];
    dx[1][0] = (-1/tau[1])*K->x[1][0]+(-1/tau[1])*K->x[4][0]+(1/tau[1])*w[1];
    dx[2][0] = (-1/tau[2])*K->x[2][0]+(-1/tau[2])*K->x[5][0]+(1/tau[2])*w[2];

    dx[6][0] = 0.5*(-K->x[0][0]*K->x[7][0]-K->x[1][0]*K->x[8][0]-K->x[2][0]*K->x[9][0]);
    dx[7][0] = 0.5*(+K->x[0][0]*K->x[6][0]-K->x[2][0]*K->x[8][0]+K->x[1][0]*K->x[9][0]);
    dx[8][0] = 0.5*(+K->x[1][0]*K->x[6][0]+K->x[2][0]*K->x[7][0]-K->x[0][0]*K->x[9][0]);
    dx[9][0] = 0.5*(+K->x[2][0]*K->x[6][0]-K->x[1][0]*K->x[7][0]+K->x[0][0]*K->x[8][0]);

    //TYPE** dxdt = dmatrix(10,1);
    matscal(dx,dx,dt,10,1);
    matadd(K->x,dx,K->x,10,1);
    //printf("state: \r\n");
    //printdmat(K->x,10,1);

    // Predict the state error covariance matrix
    // First build the linearize transition matrix F
    // q1 part
    K->F[6][0] = -0.5*q2;
    K->F[6][1] = -0.5*q3;
    K->F[6][2] = -0.5*q4;
    K->F[6][7] = -0.5*w[0];
    K->F[6][8] = -0.5*w[1];
    K->F[6][9] = -0.5*w[2];
    // q2 part
    K->F[7][0] = 0.5*q1;
    K->F[7][1] = 0.5*q4;
    K->F[7][2] = -0.5*q3;
    K->F[7][6] = 0.5*w[0];
    K->F[7][8] = -0.5*w[2];
    K->F[7][9] = 0.5*w[1];
    // q3 part
    K->F[8][0] = -0.5*q4;
    K->F[8][1] = 0.5*q1;
    K->F[8][2] = 0.5*q2;
    K->F[8][6] = 0.5*w[1];
    K->F[8][7] = 0.5*w[2];
    K->F[8][9] = -0.5*w[0];
    // q4 part
    K->F[9][0] = 0.5*q3;
    K->F[9][1] = -0.5*q2;
    K->F[9][2] = 0.5*q1;
    K->F[9][6] = 0.5*w[2];
    K->F[9][7] = -0.5*w[1];
    K->F[9][8] = 0.5*w[0];
    //printf("F:\n");
    //printdmat(K->F,10,10);

    //printf("Pk-1:\n");
    //printdmat(K->P,10,10);
    // calculate dP = F*P+P*F^T+Q
    matmult(K->F,K->P,dP1,10,10,10);
    //printf("FP:\n");
    //printdmat(dP1,10,10);
		transpose(K->F,FT,10,10);
    matmult(K->P,FT,dP2,10,10,10);
    //printf("PFT:\n");
    //printdmat(dP2,10,10);
    matadd(dP2,K->Q,dP2,10,10);
    matadd(dP1,dP2,dP,10,10);
    //printf("dP:\n");
    //printdmat(dP,10,10);
    // calculate Pk+1 = Pk+dP*dt
    matscal(dP,dP,dt,10,10);
    //printdmat(K->P,10,10);
    matadd(K->P,dP,K->P,10,10);
    //printf("Pk|k-1:\n");
    //printdmat(K->P,10,10);
}

void estimate_EKF_state(Ft* K, TYPE** y){
    //calculate Kalman gain Kg = P*H^T*(H*P*H^T+R)^-1
		// Kg1b = H*P
    matmult(K->H,K->P,Kg1b,H_dim,10,10);
    // Kg1a = Kg1b*H^T = H*P*H^T
		transpose(K->H,HT,H_dim,10);
    matmult(Kg1b,HT,Kg1a,H_dim,10,H_dim);
    // Kg1a = Kg1a+R = (H*P*H^T+R)
    matadd(Kg1a,K->R,Kg1a,H_dim,H_dim);
    //printf("(H*P*H^T+R):\n");
    //printdmat(Kg1a,H_dim,H_dim);

    // Kg1a = Kg1a^-1 = (H*P*H^T+R)^-1
    // Safer Inverse: 
    // A*A^-1=I
    // A^T*A*A^-1=A^T
    // A^-1=(A^T*A)^-1*A^T
		makeeyes(KMI, H_dim);
		transpose(Kg1a,Kg1aT,H_dim,H_dim);
    matmult(Kg1a,Kg1aT,Kg1c,H_dim,H_dim,H_dim);
    Cholesky_Inv(Kg1c, KMX, Kg1c, KML, KMLT, KMI, H_dim);
    matmult(Kg1c,Kg1aT,Kg1a,H_dim,H_dim,H_dim);
    //printf("(H*P*H^T+R)^-1:\n");
    //printdmat(Kg1a,H_dim,H_dim);
    // PHT = P*H^T
    matmult(K->P,HT,PHT,10,10,H_dim);
    // Kg = PHT*Kg1a = P*H^T*(H*P*H^T+R)^-1
    matmult(PHT,Kg1a,Kg,10,H_dim,H_dim);
    //printf("K:\n");
    //printdmat(Kg,10,H_dim);

    // calculate estimation of state x = x+Kg*(q-H*x)
    matmult(K->H,K->x,error,H_dim,10,1);
    matsub(y, error, error, H_dim, 1);
    matmult(Kg,error,fix,10,H_dim,1);
    matadd(K->x,fix,K->x,10,1);
    
    // calculate estimation of error covariance matrix P 
    // P = (I-Kg*H)*P*(I-Kg*H)^T+Kg*R*Kg^T (Joseph Form)
    /*
    // P1 = Kg*H
    matmult(Kg,K->H,P1,10,4,10);
    // P1 = I-Kg*H
    matsub(deyes(10),P1,P1,10,10);
    //printf("I-Kg*H:\n");
    //printdmat(P1,10,10);
    //printf("Pk|k-1:\n");
    //printdmat(K->P,10,10);
    // P = (I-Kg*H)*P
    matmult(P1,K->P,K->P,10,10,10);
    //printf("Pk|k:\n");
    //printdmat(K->P,10,10);
    */
    
    //printf("Pk|k-1:\n");
    //printdmat(K->P,10,10);
    //printf("K:\n");
    //printdmat(Kg,10,H_dim);
    // P3 = Kg*R
    matmult(Kg,K->R,KgR,10,H_dim,H_dim);
    // P3 = P3*Kg^T = Kg*R*Kg^T
    transpose(Kg,KgT,10,H_dim);
    matmult(KgR,KgT,P3,10,H_dim,10);
    // P1 = Kg*H
    matmult(Kg,K->H,P1,10,H_dim,10);
    // P1 = I-Kg*H
    matsub(deyes(10),P1,P1,10,10);
    // P2 = (I-Kg*H)*P
    matmult(P1,K->P,P2,10,10,10);
    // P1 = P2*(I-Kg*H)^T = (I-Kg*H)*P*(I-Kg*H)^T
		transpose(P1,P1T,10,10);
    matmult(P2,P1T,P1,10,10,10);
    // P = P1 + P3
    matadd(P1,P3,K->P,10,10);
    //printf("Pk|k:\n");
    //printdmat(K->P,10,10);
}
