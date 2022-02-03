#ifndef Qtn_GN_H
#define Qtn_GN_H

#include "basic_linear.h"

//#define SIX_AXIS //for accerelometer only IMU like mpu6050
#define NINE_AXIS //for IMU with accerelometer + compass like mpu9250

#ifdef __cplusplus
 extern "C" {
#endif

/* the rpy: {roll,pitch,yall} in radian
 * the q: {w,x,y,z}
 * It is Body 3-2-1 or Lab 1-2-3 rotation
 */
void get_rpy(TYPE** rpy, TYPE** q);

#ifdef NINE_AXIS

/* Initialize all the static matrix
 */
void initialize_Gauss_Newton_9axis();

/* get Jacobian with current quaternion (qk)
 * amd current gravity(bak)/earth magnetic(bmk) 
 * field measured on the body frame, 
 * therefore we need 9-axis imu to get those data
 * Jacobian (J) must be a 6*4 TYPE matrix
 * bak/bmk must be a 3*1 TYPE matrix
 */
void get_Jacobian_9axis(TYPE** J, TYPE **qk, TYPE **bmk, TYPE **bak);

/* get error calculated by err = ey-M*by
 * ey = [em; ea], by = [bm; ba], M = [R, 0; 0, R] and R = R(qk)
 * err must be a 6*1 TYPE matrix
 * bak/bmk/eak/emk must be a 3*1 TYPE matrix
 */
void get_error_9axis(TYPE** err, TYPE** qk, TYPE** eak, TYPE** bak, TYPE** emk, TYPE** bmk);

/* update quaternion using Gauss-Newton method
 * it also require reference gravity(eak)/earth magnetic(emk) 
 * field measured on the earth frame
 * newq must be a 4*1 TYPE matrix
 * bak/bmk/eak/emk must be a 3*1 TYPE matrix
 */
void update_quaternion_9axis(TYPE** newq, TYPE** qk, TYPE** eak, TYPE** bak, TYPE** emk, TYPE** bmk);

#elif defined SIX_AXIS

/* Initialize all the static matrix
 */
void initialize_Gauss_Newton_6axis(void);

/* get Jacobian with current quaternion (qk)
 * amd current gravity(bak) field measured on the 
 * body frame, 
 * therefore we need 6-axis imu to get those data
 * Jacobian (J) must be a 3*3 TYPE matrix
 * bak/eak must be a 3*1 TYPE matrix
 */
void get_Jacobian_6axis(TYPE** J, TYPE **qk, TYPE **bak);

/* get error calculated by err = ey-M*by
 * ey = [ea], by = [ba], M = [R] and R = R(qk)
 * err must be a 3*1 TYPE matrix
 * bak/eak must be a 3*1 TYPE matrix
 */
void get_error_6axis(TYPE** err, TYPE** qk, TYPE** eak, TYPE** bak);

/* update quaternion using Gauss-Newton method
 * it also require reference gravity(eak)/earth magnetic(emk) 
 * field measured on the earth frame
 * newq must be a 4*1 TYPE matrix
 * bak/eak must be a 3*1 TYPE matrix
 */
void update_quaternion_6axis(TYPE** newq, TYPE** qk, TYPE** eak, TYPE** bak);

#endif

#ifdef __cplusplus
 }
#endif

/* 
 *
 */

#endif
