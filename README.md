# Basic_Linear_Library
## Intro
- It's an easy C library for linear algebra and control on mcus.
- Contain 4 parts:
    - basic_linear.h
    - Quaternion_C.h
    - Quaternion_Gauss_Newton.h
    - Quaternion_EKF.h
- It should be OK to use this in your C++ project

## Feature
- For basic_linear.h
    - 2d matrix initialization
    - add, substract, dot, scaling, transpose
    - matrix inverse
    - Cholesky Decomposition
    - Linear solver for Upper/Lower triangular matrix
- For Quaternion_C.h
    - Quaternion initialization
    - Quaternion add, substract, multiply, inverse, scaling, norm
    - Rodriguez Equation
- For Quaternion_Gauss_Newton.h
    - Calculate orientation in 3D space
    - from accelerometer only or from compass+accerelometer
- For Quaternion_EKF.h
    - Sensor fusion from gyroscope and orientation data
    - Gyro bias estimation
    - Potentially provide orientation with faster response
    - Orientation data can be obtained from Quaternion_Gauss_Newton.h

## Usage
- To test this library:
    - You can test it on your linux PC
    - make all to build code
    - make clean to wipe the binary
    - please check out main.c in src/ for examples. The comment should be clear.
        - There are several sections in main.c. Uncomment required part to run different examples
- Simply drag the file to inc/ and src/
    - Please check makefile if there is any linking problem.
    - We currently rely on math.h. Please check if your env. support this.
- Please refer to .h for how to use those functions

## Note and Problem
- basic_linear.h
    - Please specify **TYPE** macro for matrix type.
    - Chelosky Inverse is only applicable on Positive-Definite Matrix
        - However, you can do A^-1=(A^T*A)^-1*A^T, A can be any kind of matrix (but should be full rank). 
        - That is because A^T*A is Positive-Definite
    - no inplace multiplication for matrix, that is, you can't do $A = A \dot B$
    - You should create matrix only once because I use malloc in it.
        - Also it is not recommended to free the matrix, especially on mcus.
- Quaternion_C.h
    - rely on basic_linear.h
- Quaternion_Gauss_Newton.h
    - rely on basic_linear.h
    - Please specify **SIX_AXIS**(accerelometer, mpu6050) or **NINE_AXIS**(acc.+compass, mpu9250) macro for different IMU configuration.
    - Use global variables. Try not to use this library from 2 different place.
- Quaternion_EKF.h
    - rely on basic_linear.h
    - Use global variables. Try not to use this library from 2 different place.
    - When the noise covariance is very small (<1e-5), the EKF will go explode

## TODO
- basic_linear.h
    - Norm, 2-Norm, inf-Norm...
    - Eigen Value Decomposition
    - Singular Value Decomposition
    - Advance indicing
- Quaternion_EKF.h
    - Consider Velocity and Position data from i.e. GPS sensor

## Reference
- [Chelosky Decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition)
    - You can check out Cholesky–Banachiewicz and Cholesky–Crout algorithms for turning this into codes
- "An Extended Kalman Filter for Quaternion-Based Orientation Estimation Using MARG Sensors", JoZo Luis Marins et.el., 2001
