#pragma once

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


class Joint {
private:
  /// definition
  std::string jointAxis_P;
  Eigen::Vector3d jointPosition_P;
  Eigen::Matrix3d jointRotation_P;

  double jointcoordinate_P; // Generalized Coordinate at each joint
  double jointVelocity_P; //Generalized Velocity at each joint

  Joint *parentJoint_;
  Joint *childJoint_;

public:
  Joint() {}
  /// Joint Constructor // 상대적인 조인트 position,pos,gc,gv,회전행렬,부모&자식joint 설정
  Joint(const Eigen::Vector3d &jointPosition,
        const std::string &jointAxis,
        const double &gc,
        const double &gv,
        Joint *parentJoint = nullptr,
        Joint *childJoint = nullptr) {
    jointAxis_P = jointAxis;
    jointPosition_P = jointPosition;
    jointcoordinate_P = gc;
    jointVelocity_P = gv;

    Eigen::Matrix3d Rot;
    if (jointAxis == "x") {
      Rot << 1, 0, 0,
          0, cos(gc), -sin(gc),
          0, sin(gc), cos(gc);
    } else if (jointAxis == "y") {
      Rot << cos(gc), 0, sin(gc),
          0, 1, 0,
          -sin(gc), 0, cos(gc);
    } else if (jointAxis == "z") {
      Rot << cos(gc), -sin(gc), 0,
          sin(gc), cos(gc), 0,
          0, 0, 1;
    } else {
      Rot << Eigen::Matrix3d::Identity();
    }
    jointRotation_P = Rot;
    parentJoint_ = parentJoint;
    childJoint_ = childJoint;
  }

  /// Set Childjoint
  void setChildJoint(Joint *childJoint){
    childJoint_ = childJoint;
    childJoint_ -> setParentJoint(this);
  }

  /// Set Parentjoint
  void setParentJoint(Joint *parentJoint){
    parentJoint_ = parentJoint;
  }

  /// Set Rotation matrix at each axis
  void setRotMat(const Eigen::Matrix3d& RotMat){
    jointRotation_P = RotMat;
}

  /// Convert a Quaternion to a Rotation Matrix
  Eigen::Matrix3d QtoR(const Eigen::VectorXd &gc) {
    Eigen::Matrix3d Rot;
    Eigen::Vector4d q; // quaternion vector
    for (int i = 0; i < 4; i++) q[i] = gc[i + 3];

    Rot << 2 * (pow(q[0], 2) + pow(q[1], 2)) - 1, 2 * (q[1] * q[2] - q[0] * q[3]), 2 * (q[1] * q[3] + q[0] * q[2]),
        2 * (q[1] * q[2] + q[0] * q[3]), 2 * (pow(q[0], 2) + pow(q[2], 2)) - 1, 2 * (q[2] * q[3] - q[0] * q[1]),
        2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 2 * (pow(q[0], 2) + pow(q[3], 2)) - 1;

    return Rot;
  }

  /// Find the Position from World to the point which I want to find(?)
  Eigen::Vector3d getWorldPosition (const Eigen::Vector3d& PointPosition){
    if (parentJoint_ == nullptr){
      return  jointPosition_P + jointRotation_P * PointPosition;
    }
    else {
      return parentJoint_ -> getWorldPosition()
    }
  }
  /// Find the Angular Velocity
  Eigen::Vector3d getAngularVelocity (const Eigen::Vector3d& AngularVelocity){
    if (parentJoint_ == nullptr){
      return jointRotation_P * AngularVelocity;
    }
    else {
      return parentJoint_ -> getAngularVelocity(jointRotation_P * AngularVelocity);
    }
  }

  Eigen::Vector3d wav(){
    Eigen::Vector3d ang =
  }

};


/// do not change the name of the method
inline Eigen::Vector3d getEndEffectorPosition(const Eigen::VectorXd &gc,const Eigen::VectorXd &gv) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////

  Joint Trunk({gc[0],gc[1],gc[2]},"x", 0, 0);
  Joint Hip({0.2399,-0.051,0},"x",gc[7],gv[8]);
  Joint Thigh({0,-0.083,0},"y",gc[8],gv[9]);
  Joint Calf({0,0,-0.25},"y",gc[9],gv[10]);
  Joint EndEffector({0,0,-0.25},"fixed",0,0);

  Trunk.setChildJoint(&Hip);
  Hip.setChildJoint(&Thigh);
  Thigh.setChildJoint(&Calf);
  Calf.setChildJoint(&EndEffector);



  return Eigen::Vector3d::Ones();
}