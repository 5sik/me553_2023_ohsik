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

  /// 풀이 순서 -> 1. 다리가 안 돌아가 있다는 전제 하(Axis가 돌아가 있지 않음)에 각속도 구하기 (축 맞춰서)
  ///         -> 2. 기본 gc에 의해 축이 돌아가 있는거 고려해서 적용

  Eigen::Vector3d x; // "axis"
  Eigen::Vector3d w; // theta_dot * axis = omega ("Angular Velocity") 로 사용
  raisim::Mat<3,3> R ; // gc값에 의해 돌아가는 축 고려해주는 Matrix로 사용

  x << 0,1,0; // for knee joint axis
  w << gv[8] * x; //Angular velocity of knee joint without theta change
  raisim::angleAxisToRotMat(x,gc[9],R);
  w<< R.e() * w; // omega for knee joint of Aliengo

  x << 0,1,0; // for pitch joint axis
  w << w + gv[7] * x ; //Angular velocity of pitch joint without theta change
  raisim::angleAxisToRotMat(x,gc[8],R);
  w << R.e() * w;

  x << 1,0,0; // for hip joint axis
  w << w + gv[6] * x;
  raisim::angleAxisToRotMat(x,gc[7],R);
  w << R.e() * w;

  raisim::Vec<4> v;// for 쿼터니안
  for(int i=0;i<4;i++);
    v[i] = gc[i+3];
  raisim::quatToRotMat(v,R)  // 00'(Trunk) 의 Rotational Matrix

  w << R.e() * w ; // Trunk axis 회전까지 고려
  w << w + Eigen::Vector3d{gv[3].gv[4],gv[5]}; // Trunk Angular Velocity 고려

  return w; /// replace this
}

