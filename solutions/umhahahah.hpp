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

--------------------------------------------------------
잠시 저장
Eigen::MatrixXd getMassMatrix() {
  Eigen::MatrixXd M; // final mass matrix
  Eigen::MatrixXd compositeMassInertia, compMassInertia ; // 6x6 composite mass inertia matrix
  Eigen::MatrixXd skew;
  std::vector<Body> first_matrix;

  M.setZero(gc_.size() -1 ,gc_.size() -1 ); // 이거는 나중에 gv들어가면 gc.size()-1 을 gv.size()로 바꾸기
  //////////여기 밑에서부터 index맞춰가면서 for문 하나 더 만들어야할듯
  for (int leg=3; leg>-1; leg--) {
    Body compositeBody = bodies_[3*leg +3];
    compositeMassInertia = getSpatialInertiaMatrix(compositeBody);

    for (int j = 3*leg+3; j > 3*leg; j--) { // j>0은 이유는 body_link는 spatial로 만들어놨고 그 위로 만들려고
      for (int i = j; i >= 3*leg+1; i--) { // j는 구하려는 COM_link, i는 j기준으로 joint 한칸씩 올라올려고
        Eigen::MatrixXd ItoJMatrix(6, 6);
        ItoJMatrix.setIdentity();
        skew = skewSymMat(bodies_[j].pos_W_ - bodies_[i].pos_W_);
        ItoJMatrix.topRightCorner(3, 3) = -skew;
        M(i + 5, j + 5) = bodies_[j].joint_.S.transpose() * compositeMassInertia * ItoJMatrix * bodies_[i].joint_.S;
        M(j + 5, i + 5) = M(i + 5, j + 5); // 대칭part도 만들어줌
        skew = skewSymMat(bodies_[j].pos_W_ - bodies_[0].pos_W_);
        ItoJMatrix.topRightCorner(3, 3) = -skew;
        M.block(j + 5, 0, 1, 6) = bodies_[j].joint_.S.transpose() * compositeMassInertia * ItoJMatrix *
                                  Eigen::MatrixXd::Identity(6, 6); // Trunk의 subspace matrix는 6x6 indentity
        M.block(0, j + 5, 6, 1) = M.block(j + 5, 0, 1, 6).transpose();
      }
      if ((j==1|j==4|j==7|j==10)){
        first_matrix.push_back(compositeBody) ;
      }else{
        compositeBody = getCompositeBody(bodies_[j-1], compositeBody, gc_[j + 6]);
        compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
      }
    }
  } //////// 제일 바깥 for문은 leg를 돌릴려고 우선!!!

  Body comp;
  for (int i=0 ; i<4 ; i++){
    comp = getCompositeBody(bodies_[0],first_matrix[i],gc_[i+6]);
  }
  compMassInertia = getSpatialInertiaMatrix(comp);

  M.topLeftCorner(6,6) = compMassInertia;

  return M;
}

---------------------------------
////<main함수쪽 다리 properties>

std::vector<Body> bodies;
Body trunk(Eigen::Vector3d{0.008465, 0.004045, -0.000763}, 9.041,
           GetInertiaMatrix(0.033260231, -0.000451628, 0.000487603, 0.16117211, 4.8356e-05, 0.17460442), 0,
           Eigen::Vector3d::Zero(), Eigen::Matrix3d::Zero(), Eigen::Vector3d::Zero(),
           Body::Joint::Type::floating);
Body imu_link(Eigen::Vector3d{0, 0, 0}, 0.001,
              GetInertiaMatrix(0.0001, 0, 0, 0.000001, 0, 0.0001), 0,
              Eigen::Vector3d{0, 0, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
              Body::Joint::Type::fixed);
bodies.push_back(getCompositeBody(trunk, imu_link, 0));

for(int leg=0; leg<4; leg++) {
// leg=0:FR / leg=1:FL / leg=2:RR / leg=3:RL
double front = 1 - 2 * (leg / 2); // front=1 , Rear=-1
double right = 1 - 2 * (leg & 2); // right=1 , left=-1

Body hip(Eigen::Vector3d{-front*0.022191, -right*0.015144, -1.5e-05}, 1.993,
         GetInertiaMatrix(0.002903894, front*right*7.185e-05, -front*1.262e-06, 0.004907517, right*1.75e-06, 0.005586944), 0,
         Eigen::Vector3d{front*0.2399, -right*0.051, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
         (leg==3 ? Body::Joint::Type::prismatic : Body::Joint::Type::revolute));
bodies.push_back(hip);

Body thigh(Eigen::Vector3d{-0.005607, right*0.003877, -0.048199}, 0.639,
           GetInertiaMatrix(0.005666803, -right*3.597e-06, 0.000491446, 0.005847229, -right*1.0086e-05, 0.000369811), 3*leg+1,
           Eigen::Vector3d{0, -right*0.083, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
           ((leg==1|leg==3) ? Body::Joint::Type::prismatic : Body::Joint::Type::revolute));
bodies.push_back(thigh);

Body calf(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}, 0.207,
          GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05), 3*leg+2,
          Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
          (leg==3 ? Body::Joint::Type::prismatic : Body::Joint::Type::revolute));
Body foot(Eigen::Vector3d{0, 0, 0}, 0.06,
          GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05), 3*leg+2,
          Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
          Body::Joint::Type::fixed);
bodies.push_back(getCompositeBody(calf, foot, 0));
}

ArticulatedSystem railab(bodies);
railab.computeForwardKinematics(gc);

------------------


Body compositeBody = bodies_[bodies_.size() - 1];
compositeMassInertia = getSpatialInertiaMatrix(compositeBody);

for (int j = bodies_.size()-1 ; j>0 ; j--){ // j>0은 이유는 body_link는 spatial로 만들어놨고 그 위로 만들려고
for (int i=1; i <= j; i++) { // j는 구하려는 COM_link, i는 j기준으로 joint 한칸씩 올라올려고
Eigen::MatrixXd ItoJMatrix(6,6);
ItoJMatrix.setIdentity();
skew = skewSymMat(bodies_[j].pos_W_ - bodies_[i].pos_W_);
ItoJMatrix.topRightCorner(3,3) = -skew;
M(i+5,j+5) = bodies_[j].joint_.S.transpose() * compositeMassInertia * ItoJMatrix * bodies_[i].joint_.S;
M(j+5,i+5) = M(i+5,j+5); // 대칭part도 만들어줌

skew = skewSymMat(bodies_[j].pos_W_ - bodies_[0].pos_W_);
ItoJMatrix.topRightCorner(3,3) = -skew;
M.block(j+5,0,1,6) = bodies_[j].joint_.S.transpose() * compositeMassInertia * ItoJMatrix * Eigen::MatrixXd::Identity(6,6); // Trunk의 subspace matrix는 6x6 indentity
M.block(0,j+5,6,1) = M.block(j+5,0,1,6).transpose();
}
std::cout<< "bodies_" <<j <<"\n" << bodies_[j].pos_W_.transpose() << std::endl;
compositeBody = getCompositeBody(bodies_[j - 1], compositeBody, gc_[j + 6]); // j+6
compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
}
M.topLeftCorner(6,6) = compositeMassInertia;
return M;

---------------------


void computeForwardKinematics(const Eigen::VectorXd &gc) {
  gc_ = gc;
  Eigen::Matrix3d rotMat; // temporary variable to save the matrix

  for (int i = 0; i < bodies_.size(); i++) { //사실상 world 기준 joint pos & rot 구하는거
    switch (bodies_[i].joint_.type_) { // // bodies_[bodies_[i].parent_] : body의 index가 parent 때문에 index 하나 내려감 (parent index로 적힘)
      case (Body::Joint::Type::floating) :
        bodies_[i].pos_W_ = gc_.head(3);
        bodies_[i].rot_W_ = QtoR(gc.segment(3, 4));
        break;
      case (Body::Joint::Type::fixed) :
        bodies_[i].pos_W_ = bodies_[bodies_[i].parent_].pos_W_ + bodies_[i].joint_.jointPos_B_;
        bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].joint_.jointRot_B_;
        break;
      case (Body::Joint::Type::revolute) :
        bodies_[i].pos_W_ =
            bodies_[bodies_[i].parent_].pos_W_ + bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].joint_.jointPos_B_;
        bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].joint_.jointRot_B_
                            * RotM(bodies_[i].joint_.jointAxis_B_(0) == 1 ? "x" : bodies_[i].joint_.jointAxis_B_(1) == 1 ? "y" : "z",gc_[i + 6]);
        bodies_[i].joint_.S.tail(3) = bodies_[i].rot_W_ * bodies_[i].joint_.jointAxis_B_;
        break;
      case (Body::Joint::Type::prismatic) :
        bodies_[i].pos_W_ = bodies_[bodies_[i].parent_].pos_W_ + bodies_[bodies_[i].parent_].rot_W_ * (bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i + 6]);
        bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].joint_.jointRot_B_ * Eigen::Matrix3d::Identity(); // 마지막은 prismatic이니까 Identity matrix
        bodies_[i].joint_.S.head(3) = bodies_[i].rot_W_ * bodies_[i].joint_.jointAxis_B_;
        break;
    }
  }
}