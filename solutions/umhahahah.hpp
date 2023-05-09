//
// Created by jemin on 23. 4. 20.
//
#pragma once

#ifndef ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_
#include <Eigen/Core>
#include "raisim/math.hpp"


enum class JointType : int {
  REVOLUTE = 0,
  PRISMATIC = 1,
  ROOT = 2
};

class Body {
public:
  explicit Body(Eigen::Vector3d comPos_B, double mass, Eigen::Matrix3d inertia_B, size_t parent,
                Eigen::Vector3d jointPos_B, Eigen::Matrix3d jointRot_B, Eigen::Vector3d jointAxis_B,
                JointType jointType = JointType::REVOLUTE) :
      comPos_B_(std::move(comPos_B)), mass_(mass), inertia_B_(std::move(inertia_B)), parent_(parent),
      jointPos_B_(std::move(jointPos_B)), jointRot_B_(std::move(jointRot_B)), jointAxis_B_(std::move(jointAxis_B)),
      jointType_(jointType){};

  JointType jointType_;

  // link
  Eigen::Vector3d comPos_B_;
  double mass_;
  Eigen::Matrix3d inertia_B_;
  size_t parent_;

  // joint
  Eigen::Vector3d jointPos_B_;
  Eigen::Matrix3d jointRot_B_;  // rotation matrix from parent to child frame when the generalized coordinates are zero
  Eigen::Vector3d jointAxis_B_;

  // variables
  Eigen::Vector3d pos_W_;
  Eigen::Matrix3d rot_W_;
  Eigen::Vector3d jointAxis_W_;
  Eigen::Vector3d comPos_W_;
  Eigen::Matrix3d inertia_W_;
};

class ArticulatedSystems {
public:
  explicit ArticulatedSystems(std::vector<Body> bodies) : bodies_(std::move(bodies)) {};

  void computeForwardKinematics() {
    raisim::Mat<3, 3> rotMat;

    for (int i=0; i<bodies_.size(); i++) {
      switch (bodies_[i].jointType_)
      {
        case JointType::REVOLUTE:
          bodies_[i].pos_W_ = bodies_[bodies_[i].parent_].pos_W_ + bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].jointPos_B_;
          raisim::angleAxisToRotMat(bodies_[i].jointAxis_B_, gc_[i-1], rotMat);
          bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].jointRot_B_ * rotMat.e();
          break;

        case JointType::PRISMATIC:
          bodies_[i].pos_W_ = bodies_[bodies_[i].parent_].pos_W_ + bodies_[bodies_[i].parent_].rot_W_ * (bodies_[i].jointPos_B_ + bodies_[i].jointAxis_B_ * gc_[i-1]);
          bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].jointRot_B_;
          break;

        case JointType::ROOT:
          bodies_[i].pos_W_ = bodies_[i].jointPos_B_;
          bodies_[i].rot_W_ = bodies_[i].jointRot_B_;

        default:
          break;
      }

      bodies_[i].jointAxis_W_ = bodies_[i].rot_W_ * bodies_[i].jointAxis_B_;
      bodies_[i].comPos_W_ = computeCOM_W(i);
      bodies_[i].inertia_W_ = bodies_[i].rot_W_ * bodies_[i].inertia_B_ * bodies_[i].rot_W_.transpose();
    }
  }

  void updateGC(const Eigen::VectorXd& gc) { gc_ = gc; }
  void updateGV(const Eigen::VectorXd& gv) { gv_ = gv; }
  Eigen::Vector3d computeCOM_W(int bodyId) { return bodies_[bodyId].pos_W_ + bodies_[bodyId].rot_W_ * bodies_[bodyId].comPos_B_; }

  std::vector<Eigen::MatrixXd> computeJacobian(int bodyId, const Eigen::Vector3d& goal_W) {  // jacobian in the world frame
    Eigen::MatrixXd posJacobian, angJacobian;
    posJacobian.setZero(3, gc_.size());
    angJacobian.setZero(3, gc_.size());

    int id = bodyId;
    Eigen::Vector3d posOffset;
    while (bodies_[id].jointType_ != JointType::ROOT) {
      switch (bodies_[id].jointType_)
      {
        case JointType::REVOLUTE:
          posOffset = goal_W - bodies_[id].pos_W_;
          posJacobian.block(0, id-1, 3, 1) = bodies_[id].jointAxis_W_.cross(posOffset);
          angJacobian.block(0, id-1, 3, 1) = bodies_[id].jointAxis_W_;
          break;

        case JointType::PRISMATIC:
          posJacobian.block(0, id-1, 3, 1) = bodies_[id].jointAxis_W_;
          break;

        default:
          break;
      }

      id = bodies_[id].parent_;
    }

    std::vector<Eigen::MatrixXd> jacobian;
    jacobian.push_back(posJacobian);
    jacobian.push_back(angJacobian);
    return jacobian;
  }

  Eigen::MatrixXd computeMassMatrix (const Eigen::VectorXd& gc) {
    Eigen::MatrixXd massMatrix;
    massMatrix.setZero(gc.size(), gc.size());

    std::vector<Eigen::MatrixXd> jacobian;
    Eigen::MatrixXd posJacobian, angJacobian;

    updateGC(gc);
    computeForwardKinematics();

    for (int i=0; i<bodies_.size(); i++) {
      jacobian = computeJacobian(i, bodies_[i].comPos_W_);
      posJacobian = jacobian[0]; angJacobian = jacobian[1];
      massMatrix += bodies_[i].mass_ * posJacobian.transpose() * posJacobian
                    + angJacobian.transpose() * bodies_[i].inertia_W_ * angJacobian;
    }
    return massMatrix;
  }

private:
  std::vector<Body> bodies_;
  Eigen::VectorXd gc_, gv_;
};

inline Eigen::Matrix3d getInertiaMatrix (double ixx, double ixy, double ixz, double iyy, double iyz, double izz){
  Eigen::Matrix3d I;
  I << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
  return I;
}

ArticulatedSystems setArticulatedSystems() {
  Body sliderBar(Eigen::Vector3d{0., 0., 0.}, 0., getInertiaMatrix(1., 0., 0., 1., 0., 0.), 0,
                 Eigen::Vector3d{0., 0., 5.}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(), JointType::ROOT);
  Body slider(Eigen::Vector3d{0., 0., 0.}, 2., getInertiaMatrix(2., 0., 0., 1., 0., 2.), 0,
              Eigen::Vector3d{0., 0., 0.}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0}, JointType::PRISMATIC);
  Body rod(Eigen::Vector3d{0., 0., 0.5}, 5., getInertiaMatrix(1., 0., 0., 1., 0., 1.), 1,
           Eigen::Vector3d{0., 0., 0.}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0}, JointType::REVOLUTE);

  std::vector<Body> bodies;
  bodies.push_back(sliderBar);
  bodies.push_back(slider);
  bodies.push_back(rod);

  ArticulatedSystems robot(bodies);
  return robot;
}

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  ArticulatedSystems cartpole = setArticulatedSystems();
  return cartpole.computeMassMatrix(gc);

}

#endif //ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_