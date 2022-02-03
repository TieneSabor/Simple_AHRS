#include "Quaternion_C.h"
#include "basic_linear.h"
#include "Quaternion_Gauss_Newton.h"
#include "Quaternion_EKF.h"
#define PI 3.14159

int main(){
    // Below are Basic Linear Library test code
    /*
    double** A = deyes(2);
    double** B = deyes(2);
    double** C = dzero(2,2);
    A[0][1] = 3;
    A[1][0] = 2;
    B[0][1] = 6;
    B[1][0] = 0.2;
    
    matadd(A,B,C,2,2);
    printdmat(C,2,2);
    matscal(B,C,4,2,2);
    printdmat(C,2,2);
    matmult(A,B,C,2,2,2);
    printdmat(C,2,2);
    */
    // Below are quaternion test code
    /*
    Qtn* q1,* q2,* qt,* qi,* qa;
    q1 = (Qtn *)malloc(sizeof(Qtn));
    q2 = (Qtn *)malloc(sizeof(Qtn));
    qi = (Qtn *)malloc(sizeof(Qtn));
    qt = (Qtn *)malloc(sizeof(Qtn));
    qa = (Qtn *)malloc(sizeof(Qtn));
    q1->w=0;
    q1->x=0;
    q1->y=0;
    q1->z=1;
    print_Qtn(q1);
    from_Rod(q2, PI/2, 0, 1, 0);
    print_Qtn(q2);
    printf("The result quaternion: ");
    mul(q1, q2, qt);
    inv(q2,qi);
    mul(qi,qt,qa);
    print_Qtn(qa);
    */
    // Below are Cholesky/Basic Linear test code
    /*
    double** A = dmatrix(3,3);
    double** B = dmatrix(3,3);
    A[0][0] = 4;
    A[0][1] = 12;
    A[0][2] = -16;
    A[1][0] = 12;
    A[1][1] = 37;
    A[1][2] = -43;
    A[2][0] = -16;
    A[2][1] = -43;
    A[2][2] = 98;
    printdmat(A,3,3);
    Cholesky_getL(A, B, 3);
    printdmat(B,3,3);
    // solve Bx = b
    double** b = dmatrix(3,1);
    b[0][0] = 2;
    b[1][0] = 7;
    b[2][0] = 0;
    double** x = dmatrix(3,1);
    Solve_LowTri_Sys(B, x, b, 3, 0);
    printdmat(x,3,1);

    // solve BTx = c
    double** c = dmatrix(3,1);
    double** BT = dmatrix(3,3);
    double** x2 = dmatrix(3,1);
    c[0][0] = 0;
    c[1][0] = 6;
    c[2][0] = 3;
    transpose(B,BT,3,3);
    Solve_UpTri_Sys(BT, x2 ,c, 3, 0);
    printdmat(x2,3,1);
    
    // do the inverse
    double** Ainv = dmatrix(3,3);
    double** CX = dmatrix(3,3);
    double** L = dmatrix(3,3);
    double** LT = dmatrix(3,3);
    double** I = deyes(3);
    double **eye = dmatrix(3,3);
    Cholesky_Inv(A, CX, Ainv, L, LT, I, 3);
    printdmat(Ainv,3,3);
    matmult(A, Ainv, eye, 3, 3, 3);
    //A*Ainv should be eye. Lets see if it is true
    printdmat(eye, 3, 3);
    */
    // Below are Quaternion Gauss test code
    // 9-axis part
    /*
    // We want to find a transform R s.t. ea=R(ba) and em=R(bm)
    double** ea = dmatrix(3,1);
    double** ba = dmatrix(3,1);
    ea[0][0] = 0;
    ea[1][0] = 0;
    ea[2][0] = -1;
    // First produce a ground truth (tf) and ea,ba,em,bm...
    Qtn* qea,* qba,* tf,* qt1,* qt2,* qem,* qbm;
    qea = (Qtn *)malloc(sizeof(Qtn));
    qba = (Qtn *)malloc(sizeof(Qtn));
    qem = (Qtn *)malloc(sizeof(Qtn));
    qbm = (Qtn *)malloc(sizeof(Qtn));
    tf = (Qtn *)malloc(sizeof(Qtn));
    qt1 = (Qtn *)malloc(sizeof(Qtn));
    qt2 = (Qtn *)malloc(sizeof(Qtn));
    qea->w = 0;
    qea->x = ea[0][0];
    qea->y = ea[1][0];
    qea->z = ea[2][0];

    double** em = dmatrix(3,1);
    double** bm = dmatrix(3,1);
    em[0][0] = 0;
    em[1][0] = 1;
    em[2][0] = 0;
    qem->w = 0;
    qem->x = em[0][0];
    qem->y = em[1][0];
    qem->z = em[2][0];

    from_Rod(tf, PI/3, 1, 3, 2);
    inv(tf,qt1);
    mul(qea, qt1, qt2);
    mul(tf, qt2, qba);
    mul(qem, qt1, qt2);
    mul(tf, qt2, qbm);
    ba[0][0] = qba->x;
    ba[1][0] = qba->y;
    ba[2][0] = qba->z;
    bm[0][0] = qbm->x;
    bm[1][0] = qbm->y;
    bm[2][0] = qbm->z;

    printf("qba:\r\n");
    print_Qtn(qba);
    printf("qea:\r\n");
    print_Qtn(qea);
    printf("qbm:\r\n");
    print_Qtn(qbm);
    printf("qem:\r\n");
    print_Qtn(qem);
    printf("tf:\r\n");
    print_Qtn(tf);
    // Initial value
    double** qk = dmatrix(4,1);
    double** newq = dmatrix(4,1);
    qk[0][0] = 1.0;
    qk[1][0] = 0;
    qk[2][0] = 0;
    qk[3][0] = 0;
    // Initialize solver: the result should get closer to the tf
    initialize_Gauss_Newton_9axis();
    for(int i=0;i<5;i++){
        update_quaternion_9axis(newq, qk, ea, ba, em, bm);
        printf("i: %d, qk: \n", i);
        printdmat(newq, 4, 1);
        qk[0][0] = newq[0][0];
        qk[1][0] = newq[1][0];
        qk[2][0] = newq[2][0];
        qk[3][0] = newq[3][0];
    }
    */
    // Below are Quaternion Gauss test code
    // 6-axis part
    // setting earth frame vector
    /*
    double** ea = dmatrix(3,1);
    double** ba = dmatrix(3,1);
    ea[0][0] = 0;
    ea[1][0] = 0;
    ea[2][0] = -1;
    // declare quaternions
    Qtn* qea,* qba,* tf,* qt1,* qt2;
    qea = (Qtn *)malloc(sizeof(Qtn));
    qba = (Qtn *)malloc(sizeof(Qtn));
    tf = (Qtn *)malloc(sizeof(Qtn));
    qt1 = (Qtn *)malloc(sizeof(Qtn));
    qt2 = (Qtn *)malloc(sizeof(Qtn));

    qea->w = 0;
    qea->x = ea[0][0];
    qea->y = ea[1][0];
    qea->z = ea[2][0];
    // produce reference body frame quaternion from ground truth (tf)
    from_Rod(tf, PI/3, 1, 3, 0);
    inv(tf,qt1);
    mul(qea, qt1, qt2);
    mul(tf, qt2, qba);
    ba[0][0] = qba->x;
    ba[1][0] = qba->y;
    ba[2][0] = qba->z;

    printf("qba:\r\n");
    print_Qtn(qba);
    printf("qea:\r\n");
    print_Qtn(qea);
    printf("tf:\r\n");
    print_Qtn(tf);
    // initial value
    double** qk = dmatrix(4,1);
    double** newq = dmatrix(4,1);
    qk[0][0] = 1.0;
    qk[1][0] = 0;
    qk[2][0] = 0;
    qk[3][0] = 0;
    // Initialize the solver
    initialize_Gauss_Newton_6axis();
    // start G-N iteration
    for(int i=0;i<5;i++){
        update_quaternion_6axis(newq, qk, ea, ba);
        printf("i: %d, qk: \n", i);
        printdmat(newq, 4, 1);
        qk[0][0] = newq[0][0];
        qk[1][0] = newq[1][0];
        qk[2][0] = newq[2][0];
        qk[3][0] = newq[3][0];
    }
    */
    // Below are Quaternion EKF test code
    
    double w[3];
    w[0] = 0.05;
    w[1] = 0.02;
    w[2] = 0.05;
    double wb[3];
    wb[0] = w[0]+0.005;
    wb[1] = w[1]+0.005;
    wb[2] = w[2]+0.005;
    Qtn* qw;
    qw = (Qtn *)malloc(sizeof(Qtn));
    qw->w=0;
    qw->x=w[0];
    qw->y=w[1];
    qw->z=w[2];

    double q[10];
    for(int i=0;i<10;i++){
        q[i] = 5e-3;
    }
    double r[4];
    for(int i=0;i<4;i++){
        r[i] = 5e-5;
    }
    double tau[3];
    double tt = 1;
    tau[0] = tt;
    tau[1] = tt;
    tau[2] = tt;

    Ft* EKF;
    EKF =  (Ft *)malloc(sizeof(Ft));
    Initialize_EKF_state(EKF, wb, q, r, tau);

    //printf("t=0, state x=\n");
    //printdmat(EKF->x,10,1);
    
    Qtn* Ori;
    Ori = (Qtn *)malloc(sizeof(Qtn));
    Ori->w=0;
    Ori->x=0;
    Ori->y=1;
    Ori->z=0;
    double dt=0.01;
    for(int i=0;i<20000;i++){
        Qtn* dq, *qt,* qt2;
        dq = (Qtn *)malloc(sizeof(Qtn));
        qt = (Qtn *)malloc(sizeof(Qtn));
        qt2 = (Qtn *)malloc(sizeof(Qtn));

        /* (dq/dt)_now=0.5*(Omega(aka qw) "quaternion multiply" Ori)
         * The equation is from "quaternion kinemetic"
         * So Ori_new ~= Ori + dq
         */
        mul(qw, Ori, qt);
        scal(qt, 0.5*dt, dq);
        add(Ori, dq, qt2);
        normalize(qt2, Ori);
        //print_Qtn(Ori);

        // Do the predict step of the EKF
        predict_EKF_state(EKF, wb, dt);

        // Produce the observe state matrix
        //double** dmOri = Qtn_to_dmatrix(Ori);
        double** y = dmatrix(7,1);
        y[0][0] = wb[0];
        y[1][0] = wb[1];
        y[2][0] = wb[2] ;
        y[3][0] = Ori->w;
        y[4][0] = Ori->x;
        y[5][0] = Ori->y;
        y[6][0] = Ori->z;

        // Do the estimate step of the EKF
        //printf("observed y=\n");
        //printdmat(y,7,1);
        //printf("P=\n");
        //printdmat(EKF->P,10,10);
        estimate_EKF_state(EKF, y);
        //printf("t=%lf, state x=\n",dt*i+dt);
        //printdmat(EKF->x,10,1);
        printf("t:, %lf, wx:, %lf\n",dt*i+dt, EKF->x[0][0]);
    }
    
}