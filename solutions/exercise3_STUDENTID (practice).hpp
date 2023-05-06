#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"
////// Get Rotation Matrix by using axis and angle
explicit Eigen::Matrix3d RotM(const std::string &axis, const double &angle) {
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
///// Vector to SkewSymMat
Eigen::Matrix3d skewSymMat(const Eigen::Vector3d &vec) {
  Eigen::Matrix3d mat;
  mat << 0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0;
  return mat;
}

/// do not change the name of the method
class Body {
public:
  class Joint {
  public:
    enum struct Type {
      fixed,
      floating,
      revolute,
      prismatic
    };

    /// variables
    Eigen::Vector3d jointPosition_W;
    Eigen::Vector3d jointAxis_W;
    Eigen::Matrix3d rotation_W;

    /// definition
    Eigen::Vector3d jointAxis_P;
    Eigen::Vector3d jointPosition_P;
    Eigen::Matrix3d rotation_P;
    Type type_;
  };

  /// basic constructor
  Body(const Joint::Type type,
       const Eigen::Vector3d &jointPosition_P,
       const Eigen::Vector3d &jointAxis_P, const double &angle, const Eigen::Matrix3d &jointRotation_P,
       const Eigen::Vector3d &massPos, const double &mass, const Eigen::Matrix3d &inertia,
       Body *parent = nullptr, Body *children = nullptr) {
    joint_.type_ = type;
    joint_.jointPosition_P = jointPosition_P;
    joint_.jointAxis_P = jointAxis_P;    angle_ = angle;    joint_.rotation_P = jointRotation_P;
    massPos_P = massPos; mass_ = mass;    inertia_ = inertia;
    parent_ = parent;    children_ = children;
  }

  void setParent(Body *parent){
    parent_ = parent;
  }

  void setChildren(Body *children) {
    children_ = children;
    children_->setParent(this);
  }

  Joint joint_;
  double angle_;

  double mass_;
  Eigen::Vector3d massPos_P; // center of the mass position
  Eigen::Matrix3d inertia_;


  Body *children_ = nullptr;
  Body *parent_ = nullptr;
};

Body getCompositeBody(const Body &parent, const Body &children, double gc) { /// gc is theta of the joint
  Body superBody = parent;
  superBody.mass_ = parent.mass_ + children.mass_;

  Eigen::Matrix3d rotM;
  rotM = RotM(children.joint_.jointAxis_P(0) == 1 ? "x" : children.joint_.jointAxis_P(1) == 1 ? "y" : "z", gc);
  Eigen::Matrix3d rotChildToParent = children.joint_.rotation_P * rotM;
  Eigen::Vector3d childrenComPos_parent = children.joint_.jointPosition_P + rotChildToParent * children.massPos_P;

};
class ArticulatedSystem {
public:
  ArticulatedSystem() : bodies

  void computeKinematicsDownTheTree(std::vector<double> &gc) {
    if (joint_.type_ == Joint::Type::floating) {
      joint_.jointPosition_P[0] = gc[0];
      joint_.jointPosition_P[1] = gc[1];
      joint_.jointPosition_P[2] = gc[2];

      /// convert quaternion to rotation matrix
      joint_.rotation_P = QtoR(Eigen::Vector4d{gc[3], gc[4], gc[5], gc[6]});
    } else if (joint_.type_ == Joint::Type::revolute) {
      if (joint_.jointAxis_P == Eigen::Vector3d{1, 0, 0}) {
        joint_.rotation_P << 1, 0, 0,
            0, cos(angle_), -sin(angle_),
            0, sin(angle_), cos(angle_);
      } else if (joint_.jointAxis_P == Eigen::Vector3d{0, 1, 0}) {
        joint_.rotation_P << cos(angle_), 0, sin(angle_),
            0, 1, 0,
            -sin(angle_), 0, cos(angle_);
      } else if (joint_.jointAxis_P == Eigen::Vector3d{0, 0, 1}) {
        joint_.rotation_P << cos(angle_), -sin(angle_), 0,
            sin(angle_), cos(angle_), 0,
            0, 0, 1;
      }

    } else if (joint_.type_ == Joint::Type::prismatic) {

    }

    for (auto &child: children_) {
      computeKinematicsDownTheTree(~~~~~쭉만들기);
    }
  }

private:
  std::vector<Body> bodies;
};

inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
//  Body trunk(Body::Joint::Type::floating,{0,0,0},{})


  return Eigen::MatrixXd::Ones(18, 18);
}