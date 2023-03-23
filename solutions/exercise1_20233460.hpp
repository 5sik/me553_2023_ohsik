//
// Created by Jemin Hwangbo on 2022/03/17.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_

#include <Eigen/Core>
#include "raisim/math.hpp"

/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////
    Eigen::Vector3d r00, r01, r12, r23, r3e,foot;

    r00<< 0,0,0.54; r01<< 0.2399,-0.051,0; r12<< 0,-0.083,0;
    r23<< 0,0,-0.25; r3e<< 0,0,-0.25;

    Eigen::Matrix3d R01, R11, R12, R22, R23, R33;
    double theta_hip = gc[7], theta_thigh = gc[8], theta_calf = gc[9];

    raisim::Mat<3,3> R00;
    raisim::Vec<4> v;

    for(int i=0;i<4;i++)
       v[i]=gc[i+3];
    raisim::quatToRotMat(v,R00);

    R01 <<  1,0,0,
            0,1,0,
            0,0,1;
    R11 <<  1,        0,         0,
            0, cos(theta_hip), -sin(theta_hip),
            0, sin(theta_hip),  cos(theta_hip);
    R12 <<  1,0,0,
            0,1,0,
            0,0,1;
    R22 <<  cos(theta_thigh), 0, sin(theta_thigh),
            0, 1,        0,
            -sin(theta_thigh), 0, cos(theta_thigh);
    R23 <<  1,0,0,
            0,1,0,
            0,0,1;
    R33 <<  cos(theta_calf), 0, sin(theta_calf),
            0, 1,        0,
            -sin(theta_calf), 0, cos(theta_calf);

    foot = r00+R00.e()*(r01+R01*R11*(r12+R12*R22*(r23+R23*R33*r3e)));

    return foot;
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
