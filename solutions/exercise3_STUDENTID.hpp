#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"


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

////// Quaternion To Rotation
Eigen::Matrix3d QtoR(Eigen::Vector4d q) {
  Eigen::Matrix3d Rot;

  Rot << 2 * (pow(q[0], 2) + pow(q[1], 2)) - 1, 2 * (q[1] * q[2] - q[0] * q[3]), 2 * (q[1] * q[3] + q[0] * q[2]),
      2 * (q[1] * q[2] + q[0] * q[3]), 2 * (pow(q[0], 2) + pow(q[2], 2)) - 1, 2 * (q[2] * q[3] - q[0] * q[1]),
      2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 2 * (pow(q[0], 2) + pow(q[3], 2)) - 1;

  return Rot;
}

////// skew Symmetric Matrix
Eigen::Matrix3d skewSymMat(const Eigen::Vector3d &vec) {
  Eigen::Matrix3d mat;
  mat << 0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0;
  return mat;
}

///// Make the Inertia Matrix using the half of the entries
Eigen::Matrix3d GetInertiaMatrix(double ixx, double ixy, double ixz, double iyy, double iyz, double izz) {
  Eigen::Matrix3d I;
  I << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
  return I;
}
///// Enumerate all of the joint types
enum struct JointType {
  REVOLUTE,
  PRISMATIC,
  FIXED,
  FLOATING,
} jointType = JointType::REVOLUTE;
///// call all parameters of the Body
class Body {
public:
  Body() {}

  Body(Eigen::Vector3d comPos_B, double mass, Eigen::Matrix3d inertia_B, size_t parent,
       Eigen::Vector3d jointPos_B, Eigen::Matrix3d jointRot_B, Eigen::Vector3d jointAxis_B,
       JointType jointType) {
    comPos_B_ = comPos_B;
    mass_ = mass;
    inertia_B_ = inertia_B;
    parent_ = parent;
    jointPos_B_ = jointPos_B;
    jointRot_B_ = jointRot_B;
    jointAxis_B_ = jointAxis_B;
    jointType_ = jointType;
  }

//private:
  JointType jointType_;

  // link
  Eigen::Vector3d comPos_B_; // Center of Mass
  double mass_;
  Eigen::Matrix3d inertia_B_;
  size_t parent_;  /// trunk = 0 hip = 1 thigh = 2

  // joint
  Eigen::Vector3d jointPos_B_;
  Eigen::Matrix3d jointRot_B_;
  Eigen::Vector3d jointAxis_B_;

  // variables
  Eigen::Vector3d pos_W_;
  Eigen::Matrix3d rot_W_;

}; // Body를 구성하는데 필요한 component 다 부르기 ~
///// change two objects into one integral object
Body getCompositeBody(const Body &body1, const Body &body2, double angle) {// frame of merged body is equal to parent(body1)
  Eigen::Matrix3d skew, rotMat;
  Body compositeBody = body1;
  compositeBody.mass_ = body1.mass_ + body2.mass_;

  rotMat = RotM(body2.jointAxis_B_(0) == 1 ? "x" : body2.jointAxis_B_(1) == 1 ? "y" : "z", angle);
  Eigen::Matrix3d rotParentToChild = body2.jointRot_B_ * rotMat; // ex) Roo' * Ro'1 느낌
  Eigen::Vector3d childCOMPos_wrt_parentOrigin = body2.jointPos_B_ + rotParentToChild * body2.comPos_B_;
  compositeBody.comPos_B_ = (body1.mass_ * body1.comPos_B_ + body2.mass_ * childCOMPos_wrt_parentOrigin) / compositeBody.mass_;

  skew = skewSymMat(body1.comPos_B_ - compositeBody.comPos_B_);
  compositeBody.inertia_B_ = body1.inertia_B_ - body1.mass_ * skew * skew;
  skew = skewSymMat(childCOMPos_wrt_parentOrigin - compositeBody.comPos_B_);
  compositeBody.inertia_B_ +=
      rotParentToChild * body2.inertia_B_ * rotParentToChild.transpose() - body2.mass_ * skew * skew;

  return compositeBody;
}
///// get spatial Inertia Matrix of each body
Eigen::MatrixXd getSpatialInertiaMatrix (const Body& body){
  Eigen::Matrix3d skew;
  Eigen::MatrixXd compositeMassInertia(6,6); // this matrix is 6 by 6
  compositeMassInertia.setZero();

  // pill in the each entries
  compositeMassInertia.topLeftCorner(3,3) = body.mass_ * Eigen::Matrix3d::Identity();
  skew = skewSymMat(body.rot_W_ * body.comPos_B_);
  compositeMassInertia.topRightCorner(3,3) = - body.mass_ * skew;
  compositeMassInertia.bottomLeftCorner(3,3) = - compositeMassInertia.topRightCorner(3,3);
  compositeMassInertia.bottomRightCorner(3,3) = body.rot_W_ * body.inertia_B_ * body.rot_W_.transpose()
                                                              - body.mass_ * skew * skew;

  return compositeMassInertia;
}
///// 대망의 articulated system 만들기
class ArticulatedSystem {
private:
  Eigen::VectorXd gc_;
  Eigen::VectorXd gv_;
  std::vector<Body> bodies_;
public:
  ArticulatedSystem(std::vector<Body> bodies) { bodies_ = bodies;}

  void  computeForwardKinematics(const Eigen::VectorXd& gc){
    Eigen::Matrix3d skew;
    Eigen::MatrixXd M; /// for mass matrix
    Eigen::MatrixXd compositeMassInertia; /// spatial mass inertia matrix
    M.setZero(18,18);
    Body compositeBody = bodies_.at(bodies_.size()-1);
    compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
    for (int j = bodies_.size()-1; j>0; j--){ /// from branch to root

    }

  }

};

inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
//  Body trunk(Body::Joint::Type::floating,{0,0,0},{})


  return Eigen::MatrixXd::Ones(18, 18);
}