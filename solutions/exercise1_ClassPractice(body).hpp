#ifndef ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_

#include <Eigen/Core>
#include "raisim/math.hpp"
//////////////////////////////////
////////// 더 완성해야함!! //////////
//////////////////////////////////

///////// Rotation Matrix /////////
Eigen::Matrix3d RotM(const std::string &axis, const double &angle) {
  Eigen::Matrix3d Rot;
  if (axis == "x") {
    Rot << 1, 0, 0,
        0, cos(angle), -sin(angle),
        0, sin(angle), cos(angle);
  } else if (axis == "y") {
    Rot << cos(angle), 0, sin(angle),
        0, 1, 0,
        -sin(angle), 0, cos(angle);
  } else if (axis == "z") {
    Rot << cos(angle), -sin(angle), 0,
        sin(angle), cos(angle), 0,
        0, 0, 1;
  }
  return Rot;
}

////// Quaternion To Rotation ////// quaternion값 맞나 확인 필요!!!!!!!!!!///////
Eigen::Matrix3d QtoR(const Eigen::VectorXd &gc) {
  Eigen::Matrix3d Rot;
  Eigen::Vector4d q; // quaternion vector
  for (int i = 0; i < 4; i++) q[i] = gc[i + 3];

  Rot << 2 * (pow(q[0], 2) + pow(q[1], 2)) - 1, 2 * (q[1] * q[2] - q[0] * q[3]), 2 * (q[1] * q[3] + q[0] * q[2]),
      2 * (q[1] * q[2] + q[0] * q[3]), 2 * (pow(q[0], 2) + pow(q[2], 2)) - 1, 2 * (q[2] * q[3] - q[0] * q[1]),
      2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 2 * (pow(q[0], 2) + pow(q[3], 2)) - 1;

  return Rot;
}


class Body {
public:
  struct Joint {
    /// variables
    Eigen::Vector3d jointPosition_W;
    Eigen::Vector3d jointAxis_W;

    /// definition
    Eigen::Vector3d jointAxis_P;
    Eigen::Vector3d jointPosition_P;
    Eigen::Matrix3d jointRotation_P;
    enum struct Type {
      fixed,
      floating,
      prismatic,
      revolute
    } type = Type::revolute;
  };

  Body(const Joint::Type type,
       const Eigen::Vector3d &jointAxis_P,
       const Eigen::Vector3d &jointPosition_P,
       const Eigen::Matrix3d &jointRotation_P) {
    joint_.type = type;
    Body(jointAxis_P, jointPosition_P, jointRotation_P);
  }

  Body(const Eigen::Vector3d &jointAxis_P,
       const Eigen::Vector3d &jointPosition_P,
       const Eigen::Matrix3d &jointRotation_P) {
    joint_.jointAxis_P = jointAxis_P;
    joint_.jointPosition_P = jointPosition_P;
    joint_.jointRotation_P = jointRotation_P;
  }


  void setChildren(const std::vector<Body> &children) {
    children_ = children;
  }

  Joint joint_;
  std::vector<Body> children_;
  Body *parent_ = nullptr;


};


/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition(const Eigen::VectorXd &gc) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////
  Eigen::Matrix3d TrunkRot = QtoR(gc);
  Eigen::Matrix3d HipRot = RotM("x", gc[7]);

  Body Trunk(
      {0, 0, 1},
      {gc[0], gc[1], gc[2]},
      TrunkRot);
  Body Hip(Body::Joint::Type::revolute,
           {1, 0, 0},
           {0.2399, -0.051, 0},
           HipRot);
  Body thigh(Body::Joint::Type::revolute,
             {0, 1, 0},
             {0, -0.083, 0})

  return Eigen::Vector3d::Ones();
}

#endif // ME553_2022_SOLUTIONS_EXERCISE1_STUDENTID_HPP_
