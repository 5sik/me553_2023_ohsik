#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "raisim/math.hpp"
#include <iostream>

///////// Rotation Matrix /////////
inline Eigen::Matrix3d RotM(const std::string &axis, const double &angle) {
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
inline Eigen::Matrix3d QtoR(Eigen::Vector4d q) {
  Eigen::Matrix3d Rot;

  Rot << 2 * (pow(q[0], 2) + pow(q[1], 2)) - 1, 2 * (q[1] * q[2] - q[0] * q[3]), 2 * (q[1] * q[3] + q[0] * q[2]),
      2 * (q[1] * q[2] + q[0] * q[3]), 2 * (pow(q[0], 2) + pow(q[2], 2)) - 1, 2 * (q[2] * q[3] - q[0] * q[1]),
      2 * (q[1] * q[3] - q[0] * q[2]), 2 * (q[2] * q[3] + q[0] * q[1]), 2 * (pow(q[0], 2) + pow(q[3], 2)) - 1;

  return Rot;
}

////// skew Symmetric Matrix
inline Eigen::Matrix3d skewSymMat(const Eigen::Vector3d &vec) {
  Eigen::Matrix3d mat;
  mat << 0, -vec[2], vec[1], vec[2], 0, -vec[0], -vec[1], vec[0], 0;
  return mat;
}

///// Make the Inertia Matrix using the half of the entries
inline Eigen::Matrix3d GetInertiaMatrix(double ixx, double ixy, double ixz, double iyy, double iyz, double izz) {
  Eigen::Matrix3d I;
  I << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
  return I;
}

class Body {
public:
  class Joint {
  public:
    enum struct Type {
      fixed=0,
      floating,
      revolute,
      prismatic
    };

    /// variables
    Eigen::Vector3d jointPos_W_;
    Eigen::Vector3d jointAxis_W_;
    Eigen::Matrix3d jointRot_W_;
    Eigen::VectorXd S;

    /// definition
    Eigen::Vector3d jointAxis_B_;
    Eigen::Vector3d jointPos_B_;
    Eigen::Matrix3d jointRot_B_;
    Type type_;
  };

  // link
  Eigen::Vector3d comPos_B_;
  double mass_;
  Eigen::Matrix3d inertia_B_;
  size_t parent_;

  // joint
//  Joint::Type type_;
  Joint joint_;

  // variables
  Eigen::Vector3d pos_W_; /// nan 원인으로 얘 초기화 해줘야할듯
  Eigen::Matrix3d rot_W_;

  /// basic constructor
  Body() {}

  Body(Eigen::Vector3d comPos_B, double mass, Eigen::Matrix3d inertia_B, size_t parent,
       Eigen::Vector3d jointPos_B, Eigen::Matrix3d jointRot_B, Eigen::Vector3d jointAxis_B,
       Joint::Type type) {
    comPos_B_ = comPos_B;
    mass_ = mass;
    inertia_B_ = inertia_B;
    parent_ = parent;
    joint_.jointPos_B_ = jointPos_B;
    joint_.jointRot_B_ = jointRot_B;
    joint_.jointAxis_B_ = jointAxis_B;
    joint_.type_ = type;
    joint_.S = Eigen::VectorXd::Zero(6);
  }

};

inline Body getCompositeBody(const Body &parent, const Body &child, double gc) { //여기서 gc는 child.joint가 rotate joint면 돌려야하니깐 필요 각도
  Eigen::Matrix3d skew1, skew2, rot;
  Body compositeBody = parent;
  compositeBody.mass_ = parent.mass_ + child.mass_;

  // child.joint_.jointRot_B_ matrix는 사실 Identity, because 보통 붙어 있는 body들은 축이 돌아가 있지 않고 같은 방향을 보임 (child.axis가 시작부터 돌아있는 각을 뜻함)
  // 일반적으로는 축이 돌아가 있지를 않아서 Identity

  Eigen::Matrix3d childRot_parent_origin;
  Eigen::Vector3d childComPos_parent_origin;
//   parent joint 기준점으로 child.joint 거리에다가 회전이 반영된 child.body_com 까지 거리
  if ( child.joint_.type_ == Body::Joint::Type::prismatic) {
    childRot_parent_origin = child.joint_.jointRot_B_;
    childComPos_parent_origin = (child.joint_.jointPos_B_ + child.joint_.jointRot_B_ * child.joint_.jointAxis_B_ * gc) + child.joint_.jointRot_B_ * child.comPos_B_;
  } else {
    rot = RotM(child.joint_.jointAxis_B_(0) == 1 ? "x" : child.joint_.jointAxis_B_(1) == 1 ? "y" : "z", gc);  // rot은 child body에 붙어있는 기준 joint가 rotate joint 이면 돌려야하니깐
    childRot_parent_origin = child.joint_.jointRot_B_ * rot;
    childComPos_parent_origin = (child.joint_.jointPos_B_) + childRot_parent_origin * child.comPos_B_;
  }

  compositeBody.comPos_B_ =
      (parent.mass_ * parent.comPos_B_ + child.mass_ * childComPos_parent_origin) / compositeBody.mass_;

  skew1 = skewSymMat(parent.comPos_B_ - compositeBody.comPos_B_);
  skew2 = skewSymMat(childComPos_parent_origin - compositeBody.comPos_B_);
  compositeBody.inertia_B_ =
      parent.inertia_B_ + childRot_parent_origin * child.inertia_B_ * childRot_parent_origin.transpose()
      - parent.mass_ * skew1 * skew1 - child.mass_ * skew2 * skew2;

  return compositeBody;
}

inline Eigen::MatrixXd getSpatialInertiaMatrix(const Body &body) {
  Eigen::Matrix3d skew;
  Eigen::MatrixXd spatial_inertia_matrix(6, 6);
  spatial_inertia_matrix.setZero();

  // pill in the each entries
  spatial_inertia_matrix.topLeftCorner(3, 3) = body.mass_ * Eigen::Matrix3d::Identity();
  skew = skewSymMat(body.rot_W_ * body.comPos_B_);
  spatial_inertia_matrix.bottomLeftCorner(3, 3) = body.mass_ * skew;
  spatial_inertia_matrix.topRightCorner(3, 3) = -spatial_inertia_matrix.bottomLeftCorner(3, 3);
  spatial_inertia_matrix.bottomRightCorner(3, 3) = body.rot_W_ * body.inertia_B_ * body.rot_W_.transpose()
                                                   - body.mass_ * skew * skew;

  return spatial_inertia_matrix;
}

class ArticulatedSystem {
public:
  ArticulatedSystem(std::vector<Body> bodies) { bodies_ = bodies; }

  void computeForwardKinematics(const Eigen::VectorXd &gc) {
    gc_ = gc;
    Eigen::Matrix3d rotMat; // temporary variable to save the matrix
    bodies_[0].pos_W_ = gc_.head(3);
    bodies_[0].rot_W_ = QtoR(gc.segment(3, 4));

    for (int i = 0; i < bodies_.size(); i++) { //사실상 world 기준 joint pos & rot 구하는거
      switch (bodies_[i].joint_.type_) { // // bodies_[bodies_[i].parent_] : body의 index가 parent 때문에 index 하나 내려감 (parent index로 적힘)
//        case (Body::Joint::Type::floating) :
//          bodies_[i].pos_W_ = gc_.head(3);
//          bodies_[i].rot_W_ = QtoR(gc.segment(3, 4));
//          break;
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
//          std::cout<< "bodies_[parent]"<<bodies_[bodies_[i].parent_].pos_W_.transpose() <<std::endl;
//          std::cout<< "bodies_[joint "<<i<<"]"<<bodies_[i].joint_.jointPos_B_.transpose() <<std::endl;
//          std::cout<< "bodies_["<<i<<"]"<<bodies_[i].pos_W_.transpose() <<std::endl;
          break;
        case (Body::Joint::Type::prismatic) :
          bodies_[i].pos_W_ = bodies_[bodies_[i].parent_].pos_W_ + bodies_[bodies_[i].parent_].rot_W_ * (bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i + 6]);
          bodies_[i].rot_W_ = bodies_[bodies_[i].parent_].rot_W_ * bodies_[i].joint_.jointRot_B_ * Eigen::Matrix3d::Identity(); // 마지막은 prismatic이니까 Identity matrix
          bodies_[i].joint_.S.head(3) = bodies_[i].rot_W_ * bodies_[i].joint_.jointAxis_B_;
//          std::cout<< "bodies_["<<i<<"]"<<bodies_[i].pos_W_.transpose() <<std::endl;
          break;
      }
    }
  }
  Eigen::MatrixXd getMassMatrix() {
    Eigen::MatrixXd M; // final mass matrix
    Eigen::MatrixXd compositeMassInertia; // 6x6 composite mass inertia matrix
    Eigen::MatrixXd skew;
    Body compositeBody;
    std::vector<Body> temporary; // temperory buffer

    M.setZero(gc_.size() - 1 ,gc_.size() - 1 );
    compositeMassInertia.setZero(6,6);

    for (int leg=3; leg>-1 ; leg--){
    compositeBody = bodies_[3*leg+3];
    compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
      for (int j =3*leg+3 ; j>3*leg ; j--){ // j>0은 이유는 body_link는 spatial로 만들어놨고 그 위로 만들려고
        for (int i=3*leg+1; i<=j; i++) { // j는 구하려는 COM_link, i는 j기준으로 joint 한칸씩 올라올려고
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
        if(j==1|j==4|j==7|j==10){
          temporary.push_back(compositeBody);
          compositeBody = getCompositeBody(bodies_[0], compositeBody, gc_[j + 6]); // j+6
          compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
        } else {
          compositeBody = getCompositeBody(bodies_[j - 1], compositeBody, gc_[j + 6]); // j+6
          compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
        }
        // compositeBody = getCompositeBody((j==1|j==4|j==7|j==10) ? bodies_[0] : bodies_[j - 1], compositeBody, gc_[j + 6]); // j+6
        // compositeMassInertia = getSpatialInertiaMatrix(compositeBody);
      }
    }

    compositeBody = getCompositeBody(bodies_[0],temporary[0],gc_[16]); // RL->RR->FL->FR 순서의 힙 gc필요
    compositeBody = getCompositeBody(compositeBody,temporary[1],gc_[13]);
    compositeBody = getCompositeBody(compositeBody,temporary[2],gc_[10]);
    compositeBody = getCompositeBody(compositeBody,temporary[3],gc_[7]);
    compositeMassInertia = getSpatialInertiaMatrix(compositeBody);

    M.topLeftCorner(6,6) = compositeMassInertia;
    return M;
  }


private:
  std::vector<Body> bodies_;
  Eigen::VectorXd gc_;

};


inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

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

  ////FR
  Body FR_hip(Eigen::Vector3d{-0.022191, -0.015144, -1.5e-05}, 1.993,
              GetInertiaMatrix(0.002903894, 7.185e-05, -1.262e-06, 0.004907517, 1.75e-06, 0.005586944), 0,
              Eigen::Vector3d{0.2399, -0.051, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
              Body::Joint::Type::revolute);
  bodies.push_back(FR_hip);

  Body FR_thigh(Eigen::Vector3d{-0.005607, 0.003877, -0.048199}, 0.639,
                GetInertiaMatrix(0.005666803, -3.597e-06, 0.000491446, 0.005847229, -1.0086e-05, 0.000369811), 1,
                Eigen::Vector3d{0, -0.083, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
                Body::Joint::Type::revolute);
  bodies.push_back(FR_thigh);

  Body FR_calf(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}, 0.207,
               GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05), 2,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
               Body::Joint::Type::revolute);
  Body FR_foot(Eigen::Vector3d{0, 0, 0}, 0.06,
               GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05),2,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
               Body::Joint::Type::fixed);
  bodies.push_back(getCompositeBody(FR_calf, FR_foot, 0));

      ////FL
  Body FL_hip(Eigen::Vector3d{-0.022191, 0.015144, -1.5e-05}, 1.993,
                  GetInertiaMatrix(0.002903894, -7.185e-05, -1.262e-06, 0.004907517, -1.75e-06, 0.005586944), 0,
                  Eigen::Vector3d{0.2399, 0.051, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
                  Body::Joint::Type::revolute);
  bodies.push_back(FL_hip);

  Body FL_thigh(Eigen::Vector3d{-0.005607, -0.003877, -0.048199}, 0.639,
                GetInertiaMatrix(0.005666803, 3.597e-06, 0.000491446, 0.005847229, 1.0086e-05, 0.000369811), 4,
                Eigen::Vector3d{0, 0.083, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
                Body::Joint::Type::prismatic);
  bodies.push_back(FL_thigh);

  Body FL_calf(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}, 0.207,
               GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05), 5,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
               Body::Joint::Type::revolute);
  Body FL_foot(Eigen::Vector3d{0, 0, 0}, 0.06,
               GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05),5,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
               Body::Joint::Type::fixed);
  bodies.push_back(getCompositeBody(FL_calf, FL_foot, 0));

      ////RR
  Body RR_hip(Eigen::Vector3d{0.022191, -0.015144, -1.5e-05}, 1.993,
                  GetInertiaMatrix(0.002903894, -7.185e-05, 1.262e-06, 0.004907517, 1.75e-06, 0.005586944), 0,
                  Eigen::Vector3d{-0.2399, -0.051, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
                  Body::Joint::Type::revolute);
  bodies.push_back(RR_hip);

  Body RR_thigh(Eigen::Vector3d{-0.005607, 0.003877, -0.048199}, 0.639,
                GetInertiaMatrix(0.005666803, -3.597e-06, 0.000491446, 0.005847229, -1.0086e-05, 0.000369811), 7,
                Eigen::Vector3d{0, -0.083, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
                Body::Joint::Type::revolute);
  bodies.push_back(RR_thigh);

  Body RR_calf(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}, 0.207,
               GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05), 8,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
               Body::Joint::Type::revolute);
  Body RR_foot(Eigen::Vector3d{0, 0, 0}, 0.06,
               GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05),8,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
               Body::Joint::Type::fixed);
  bodies.push_back(getCompositeBody(RR_calf, RR_foot, 0));
      ////RL
  Body RL_hip(Eigen::Vector3d{0.022191, 0.015144, -1.5e-05}, 1.993,
                  GetInertiaMatrix(0.002903894, 7.185e-05, 1.262e-06, 0.004907517, -1.75e-06, 0.005586944), 0,
                  Eigen::Vector3d{-0.2399, 0.051, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
                  Body::Joint::Type::prismatic);
  bodies.push_back(RL_hip);

  Body RL_thigh(Eigen::Vector3d{-0.005607, -0.003877, -0.048199}, 0.639,
                GetInertiaMatrix(0.005666803, 3.597e-06, 0.000491446, 0.005847229, 1.0086e-05, 0.000369811), 10,
                Eigen::Vector3d{0, 0.083, 0}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
                Body::Joint::Type::prismatic);
  bodies.push_back(RL_thigh);

  Body RL_calf(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}, 0.207,
               GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05), 11,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
               Body::Joint::Type::prismatic);
  Body RL_foot(Eigen::Vector3d{0, 0, 0}, 0.06,
               GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05),11,
               Eigen::Vector3d{0, 0, -0.25}, Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero(),
               Body::Joint::Type::fixed);
  bodies.push_back(getCompositeBody(RL_calf, RL_foot, 0));

  ArticulatedSystem railab(bodies);
  railab.computeForwardKinematics(gc);

//  for (int i=0;i<bodies.size();i++){
//    std::cout << "bodies["<<i<<"] = "<<bodies[i].pos_W_.transpose()<<std::endl;
//  std::cout << static_cast<int>(bodies[i].joint_.type_) <<std::endl;
//  }



  return railab.getMassMatrix();
}