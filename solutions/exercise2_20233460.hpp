#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"

/// do not change the name of the method

inline Eigen::Vector3d getEndEffectorPosition (const Eigen::VectorXd& gc){

    Eigen::Vector3d position;
    position << 0,0,-0.25; // RL_calf_joint to RL_foot_fixed
    raisim::Mat<3,3> rot;

    raisim::angleAxisToRotMat(Eigen::Vector3d{0,1,0},gc[9],rot);
    position << rot.e() * position;
    position << position + Eigen::Vector3d{0,0,-0.25}; // RL_thigh_fixed to RL_foot_fixed

    raisim::angleAxisToRotMat(Eigen::Vector3d{0,1,0},gc[8],rot);
    position << rot.e() * position;
    position << position + Eigen::Vector3d{0,-0.083,0}; // RL_hip_joint to RL_thigh_fixed

    raisim::angleAxisToRotMat(Eigen::Vector3d{1,0,0},gc[7],rot);
    position << rot.e() * position;
    position << position + Eigen::Vector3d{0.2399, -0.051, 0}; // base to RL_foot_fixed
    position << position + Eigen::Vector3d {gc[0],gc[1],gc[2]};

    return position ; // EndEffectorPosition
}

inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    //////////////////////////
    ///// Your Code Here /////
    //////////////////////////

    return Eigen::Vector3d::Ones(); /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////
  Eigen::Matrix<3,3> R ;

  raisim::Vec<4> v;// for 쿼터니안
  raisim::Mat<3,3> R00;

  for(int i=0;i<4;i++);
    v[i] = gc[i+3];
  raisim::quatToRotMat(v,R00)  // 00' 의 Rotational Matrix


  R = R00.e() + sex + sex + sex

  return Eigen::Vector3d::Ones(); /// replace this
}

