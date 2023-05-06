#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"

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


/// do not change the name of the method
class Body {
public:
  class Joint {
  public:
    /// variables
    Eigen::Vector3d jointPosition_W;
    Eigen::Vector3d jointAxis_W;
    Eigen::Matrix3d rotation_W;

    /// definition
    Eigen::Vector3d jointAxis_P;
    Eigen::Vector3d jointPosition_P;
    Eigen::Matrix3d rotation_P;
    enum struct Type {
      fixed,
      floating,
      revolute,
      prismatic
    } type = Type::revolute;
  };

  /////

  const double& angle ;

  /// basic constructor
  Body(const Joint::Type type,
       const Eigen::Vector3d& jointAxis_P,
       const Eigen::Vector3d& jointPosition_P,
       const double& angle,
       const Eigen::Matrix3d& jointRotation_P) {
    joint_.jointAxis_P = jointAxis_P;
    joint_.jointPosition_P = jointPosition_P;
  }

  void computeKinematicsDownTheTree(std::vector<double>& gc){
    if(joint_.type==Joint::Type::floating){
      joint_.jointPosition_P[0] = gc[0];
      joint_.jointPosition_P[1] = gc[1];
      joint_.jointPosition_P[2] = gc[2];

      /// convert quaternion to rotation matrix
      joint_.rotation_P = QtoR(Eigen::Vector4d {gc[3],gc[4],gc[5],gc[6]});
    }
    else if (joint_.type==Joint::Type::revolute){
      if(joint_.jointAxis_P = Eigen::Vector3d {1,0,0}){
        joint_.rotation_P << 1, 0, 0,
                              0, cos(angle), -sin(angle),
                              0, sin(angle), cos(angle);
      }
      else if(joint_.jointAxis_P = Eigen::Vector3d {0,1,0}){
        joint_.rotation_P << cos(angle), 0, sin(angle),
                              0, 1, 0,
                              -sin(angle), 0, cos(angle);
      }
      else if(joint_.jointAxis_P = Eigen::Vector3d {0,0,1}){
        joint_.rotation_P << cos(angle), -sin(angle), 0,
                        sin(angle), cos(angle), 0,
                        0, 0, 1;
      }

    }
    else if (joint_.type==Joint::Type::prismatic){

    }

    for (auto& child : children_){
      computeKinematicsDownTheTree(~~~~~ 쭉만들기);
    }
  }

  void setChildren(const std::vector<Body>& children){
    children_ = children;
  }

  Joint joint_;
  std::vector<Body> children_;
  Body* parent_ = nullptr;

};


inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
//  Body trunk(Body::Joint::Type::floating,{0,0,0},{})


  return Eigen::MatrixXd::Ones(18,18);
}