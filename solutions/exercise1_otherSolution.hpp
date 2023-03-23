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

    Eigen::Vector3d position;
    raisim::Mat<3,3> rot;

    position << 0,0,-0.25; // RL_calf_joint to RL_foot_fixed

    raisim::angleAxisToRotMat(Eigen::Vector3d{0,1,0}, gc(9),rot);
    position<< rot.e() * position;
    position<< position + Eigen::Vector3d{0,0,-0.25}; // RL_thigh_fixed to RL_foot_fixed

    raisim::angleAxisToRotMat(Eigen::Vector3d{0, 1, 0}, gc(8), rot);
    position << rot.e() * position;
    position << position + Eigen::Vector3d{0, -0.083, 0}; // RL_hip_joint to RL_foot_fixed

    raisim::angleAxisToRotMat(Eigen::Vector3d{1, 0, 0}, gc(7), rot);
    position << rot.e() * position;
    position << position + Eigen::Vector3d{0.2399, -0.051, 0}; // base to RL_foot_fixed
    position << position + Eigen::Vector3d{gc[0],gc[1],gc[2]}; // world to RL_foot_fixed


    return position ; /// replace this
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
