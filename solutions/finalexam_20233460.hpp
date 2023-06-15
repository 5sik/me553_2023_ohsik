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

///// Make the Inertia Matrix by using the half of the entries
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
    Eigen::Vector3d jointLinVel_W_;
    Eigen::Vector3d jointAngVel_W_;
    Eigen::Vector3d jointLinAcc_W_;
    Eigen::Vector3d jointAngAcc_W_;
    Eigen::VectorXd S = Eigen::VectorXd::Zero(6) ;
    Eigen::VectorXd S_dot = Eigen::VectorXd::Zero(6);
    Eigen::MatrixXd S_trunk = Eigen::MatrixXd::Zero(6,6) ;
    Eigen::MatrixXd S_dot_trunk = Eigen::MatrixXd::Zero(6,6);
    Eigen::VectorXd W = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd W_dot = Eigen::VectorXd::Zero(6);

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
  Joint joint_;

  // variables
  Eigen::Vector3d comPos_W_;
  Eigen::Matrix3d comRot_W_;
  Eigen::Matrix3d inertia_W_;
  Eigen::MatrixXd X_BP = Eigen::MatrixXd::Zero(6,6);
  Eigen::MatrixXd X_BP_dot = Eigen::MatrixXd::Zero(6,6);
  Eigen::MatrixXd articulated_M;
  Eigen::VectorXd articulated_b;

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
  skew = skewSymMat(body.joint_.jointRot_W_ * body.comPos_B_);
  spatial_inertia_matrix.bottomLeftCorner(3, 3) = body.mass_ * skew;
  spatial_inertia_matrix.topRightCorner(3, 3) = -spatial_inertia_matrix.bottomLeftCorner(3, 3);
  spatial_inertia_matrix.bottomRightCorner(3, 3) = body.joint_.jointRot_W_ * body.inertia_B_ * body.joint_.jointRot_W_.transpose()
                                                   - body.mass_ * skew * skew;

  return spatial_inertia_matrix;
}

inline Eigen::VectorXd getFictitiousForces(const Body &body) {
  Eigen::VectorXd FictitiousForces = Eigen::VectorXd::Zero(6);
  Eigen::Matrix3d omega_skew = skewSymMat(body.joint_.jointAngVel_W_);
  Eigen::Matrix3d pos_skew = skewSymMat(body.comPos_W_);
  FictitiousForces << body.mass_* omega_skew* omega_skew* body.comPos_W_, omega_skew*(body.inertia_W_-body.mass_*pos_skew*pos_skew)*body.joint_.jointAngVel_W_;

  return FictitiousForces;
}

class ArticulatedSystem {
public:
  ArticulatedSystem(std::vector<Body> bodies) { bodies_ = bodies; }

  void computeForwardKinematics(const Eigen::VectorXd &gc, const Eigen::VectorXd &gv, const Eigen::VectorXd &gf) {
    gc_ = gc;
    gv_ = gv;
    gf_ = gf;
    Eigen::Matrix3d rotMat; // temporary variable to save the matrix
    Eigen::Matrix3d parent_omega = Eigen::Matrix3d::Zero(); // Angular Velocity of parent joint
    Eigen::Matrix3d parent_alpha = Eigen::Matrix3d::Zero(); // Angular Acceleration of parent joint
    Eigen::Vector3d RelativePos_W_ = Eigen::Vector3d::Zero(); // Relative Position between Parent and Child Joint

    for (int i = 0; i < bodies_.size(); i++) { //사실상 world 기준 joint pos & rot 구하는거
      RelativePos_W_ = bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointPos_B_;
      parent_omega = skewSymMat(bodies_[bodies_[i].parent_].joint_.jointAngVel_W_);
      parent_alpha = skewSymMat(bodies_[bodies_[i].parent_].joint_.jointAngAcc_W_);
      switch (bodies_[i].joint_.type_) { // // bodies_[bodies_[i].parent_] : body의 index가 parent 때문에 index 하나 내려감 (parent index로 적힘)
        case (Body::Joint::Type::fixed) :
          bodies_[i].joint_.jointPos_W_.setZero();
          bodies_[i].joint_.jointRot_W_.setIdentity(3,3);
          bodies_[i].joint_.S_trunk = Eigen::MatrixXd::Zero(6,6); // Trunk = Base Frame
          bodies_[i].X_BP.setIdentity(6,6);
          bodies_[i].joint_.jointLinVel_W_.setZero();
          bodies_[i].joint_.jointAngVel_W_.setZero();
          bodies_[i].joint_.W << bodies_[i].joint_.jointLinVel_W_ , bodies_[i].joint_.jointAngVel_W_;
          bodies_[i].joint_.S_dot_trunk = Eigen::MatrixXd::Zero(6,6);
          bodies_[i].joint_.jointLinAcc_W_ = Eigen::Vector3d {0,0,9.81};
          bodies_[i].joint_.jointAngAcc_W_.setZero();
          bodies_[i].joint_.W_dot << bodies_[i].joint_.jointLinAcc_W_ , bodies_[i].joint_.jointAngAcc_W_;
          break;
        case (Body::Joint::Type::revolute) :
          bodies_[i].joint_.jointPos_W_ = //bodies_[bodies_[i].parent_].joint_.jointPos_W_ + RelativePos_W_;
              bodies_[bodies_[i].parent_].joint_.jointPos_W_ + bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointPos_B_;
          bodies_[i].joint_.jointRot_W_ =
              bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointRot_B_
              * RotM(bodies_[i].joint_.jointAxis_B_(0) == 1 ? "x" : bodies_[i].joint_.jointAxis_B_(1) == 1 ? "y" : "z",gc_[i-1]);
          bodies_[i].joint_.S.tail(3) = bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_;
          bodies_[i].X_BP.setIdentity(6,6);
          bodies_[i].X_BP.bottomLeftCorner(3,3) = skewSymMat(bodies_[i].joint_.jointPos_W_-bodies_[bodies_[i].parent_].joint_.jointPos_W_);
          bodies_[i].joint_.jointLinVel_W_ =
              bodies_[bodies_[i].parent_].joint_.jointLinVel_W_ + parent_omega * RelativePos_W_;
          bodies_[i].joint_.jointAngVel_W_ = bodies_[bodies_[i].parent_].joint_.jointAngVel_W_ +
                                             bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_ * gv_[i-1];
          bodies_[i].joint_.W << bodies_[i].joint_.jointLinVel_W_ , bodies_[i].joint_.jointAngVel_W_;
          bodies_[i].X_BP_dot.setZero(6,6);
          bodies_[i].X_BP_dot.bottomLeftCorner(3,3) = skewSymMat(bodies_[i].joint_.jointLinVel_W_ - bodies_[bodies_[i].parent_].joint_.jointLinVel_W_);
          bodies_[i].joint_.jointLinAcc_W_ =
              bodies_[bodies_[i].parent_].joint_.jointLinAcc_W_ + parent_alpha * RelativePos_W_ + parent_omega * parent_omega * RelativePos_W_;
          bodies_[i].joint_.jointAngAcc_W_ = bodies_[bodies_[i].parent_].joint_.jointAngAcc_W_ +
                                             parent_omega * bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[bodies_[i].parent_].joint_.jointAxis_B_ * gv_[i-1];
          bodies_[i].joint_.S_dot.tail(3) =
              parent_omega * bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_;
          break;
        case (Body::Joint::Type::prismatic) :
          bodies_[i].joint_.jointPos_W_ = bodies_[bodies_[i].parent_].joint_.jointPos_W_ +
                                          bodies_[bodies_[i].parent_].joint_.jointRot_W_ * (bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i-1]);
          bodies_[i].joint_.jointRot_W_ =
              bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointRot_B_ *
              Eigen::Matrix3d::Identity(); // 마지막은 prismatic이니까 Identity matrix
          bodies_[i].joint_.S.head(3) = bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_;
          bodies_[i].X_BP.setIdentity(6,6);
          bodies_[i].X_BP.bottomLeftCorner(3,3) = skewSymMat(bodies_[i].joint_.jointPos_W_-bodies_[bodies_[i].parent_].joint_.jointPos_W_);
          bodies_[i].joint_.jointLinVel_W_ = bodies_[bodies_[i].parent_].joint_.jointLinVel_W_
                                             + parent_omega * bodies_[bodies_[i].parent_].joint_.jointRot_W_ * (bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i-1])
                                             + bodies_[i].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_ * gv_[i-1];
          bodies_[i].joint_.jointAngVel_W_ = bodies_[bodies_[i].parent_].joint_.jointAngVel_W_;
          bodies_[i].joint_.W << bodies_[i].joint_.jointLinVel_W_ , bodies_[i].joint_.jointAngVel_W_;
          bodies_[i].X_BP_dot.setZero(6,6);
          bodies_[i].X_BP_dot.bottomLeftCorner(3,3) = skewSymMat(bodies_[i].joint_.jointLinVel_W_ - bodies_[bodies_[i].parent_].joint_.jointLinVel_W_);
          bodies_[i].joint_.jointLinAcc_W_ = bodies_[bodies_[i].parent_].joint_.jointLinAcc_W_ + parent_alpha * bodies_[bodies_[i].parent_].joint_.jointRot_W_ * (bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i-1])
                                             + parent_omega * parent_omega * bodies_[bodies_[i].parent_].joint_.jointRot_W_ *(bodies_[i].joint_.jointPos_B_ + bodies_[i].joint_.jointAxis_B_ * gc_[i-1])
                                             + 2* parent_omega * bodies_[bodies_[i].parent_].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_ * gv_[i-1];
          bodies_[i].joint_.jointAngAcc_W_ = bodies_[bodies_[i].parent_].joint_.jointAngAcc_W_;
          bodies_[i].joint_.S_dot.head(3) =
              skewSymMat(bodies_[i].joint_.jointAngVel_W_) * bodies_[i].joint_.jointRot_W_ * bodies_[i].joint_.jointAxis_B_;
          break;
      }
      bodies_[i].comPos_W_ = bodies_[i].joint_.jointRot_W_ * bodies_[i].comPos_B_;
      bodies_[i].inertia_W_ =
          bodies_[i].joint_.jointRot_W_ * bodies_[i].inertia_B_ * bodies_[i].joint_.jointRot_W_.transpose();
    }
  }


  void computeArticulatedVariable(){
    //alias for the child body notation
    Eigen::MatrixXd M_arti;    Eigen::VectorXd b_arti;
    Eigen::VectorXd W_;
    Eigen::VectorXd S_;    Eigen::VectorXd S_dot_;
    Eigen::MatrixXd X_;    Eigen::MatrixXd X_dot_;
    std::vector<Eigen::MatrixXd> M_buffer;
    std::vector<Eigen::VectorXd> b_buffer;
    Eigen::MatrixXd M_temp = Eigen::MatrixXd::Zero(6,6);
    Eigen::VectorXd b_temp = Eigen::VectorXd::Zero(6);

    for (int i=3; i>1; i--){
      bodies_[3].articulated_M = getSpatialInertiaMatrix(bodies_[2]);
      bodies_[3].articulated_b = getFictitiousForces(bodies_[2]);
      M_arti = bodies_[i].articulated_M; // articulated Mass of the child body
      b_arti = bodies_[i].articulated_b; // articulated Nonlinearities term of the child body
      W_ = bodies_[bodies_[i].parent_].joint_.W; // W_parent
      S_ = bodies_[i].joint_.S; // S_child
      S_dot_ = bodies_[i].joint_.S_dot; // S_dot_child
      X_ = bodies_[i].X_BP; // X_BP_parent2child
      X_dot_ = bodies_[i].X_BP_dot; // X_BP_dot_parent2child
      M_temp = X_ * M_arti *
               (-S_ * (S_.transpose() * M_arti * S_).inverse() * (S_.transpose() * M_arti * X_.transpose()) + X_.transpose());
      b_temp = X_ * (M_arti * (S_ * (S_.transpose() * M_arti * S_).inverse() *
                               (gf_[i-1]-S_.transpose() * M_arti * (S_dot_ * gv_[i-1] + X_dot_.transpose() * W_)
                                - S_.transpose() * b_arti) + S_dot_ * gv_[i-1] + X_dot_.transpose() * W_) + b_arti);

      bodies_[bodies_[i].parent_].articulated_M = getSpatialInertiaMatrix(bodies_[bodies_[i].parent_]) + M_temp;
      bodies_[bodies_[i].parent_].articulated_b = getFictitiousForces(bodies_[bodies_[i].parent_]) + b_temp;
      std::cout<< "bodies_["<<bodies_[i].parent_<<"_M] \n"<< bodies_[bodies_[i].parent_].articulated_M <<std::endl;
      std::cout<< "bodies_["<<bodies_[i].parent_<<"_b] \n"<< bodies_[bodies_[i].parent_].articulated_b.transpose() <<std::endl;
    }
  }

  Eigen::VectorXd computeGeneralizedAcceleration(){
    Eigen::MatrixXd M_arti;    Eigen::VectorXd b_arti;
    Eigen::VectorXd S_;    Eigen::VectorXd S_dot_;
    Eigen::MatrixXd X_;    Eigen::MatrixXd X_dot_;
    Eigen::VectorXd W_;    Eigen::VectorXd W_dot_; // W_dot_ is each joint linear & Angular Acceleration of the Articulated bodies w.r.t world frame
    Eigen::VectorXd ga_ = Eigen::VectorXd::Zero(3);; // ga is gerneralized acceleration
    Eigen::MatrixXd S_trunk = bodies_[0].joint_.S_trunk;
    Eigen::MatrixXd S_dot_trunk = bodies_[0].joint_.S_dot_trunk;

    for (int i=1; i<4; i++){
      M_arti = bodies_[i].articulated_M;
      b_arti =  bodies_[i].articulated_b;
      S_ = bodies_[i].joint_.S;
      S_dot_ = bodies_[i].joint_.S_dot;
      X_ = bodies_[i].X_BP;
      X_dot_ = bodies_[i].X_BP_dot;
      W_ = bodies_[bodies_[i].parent_].joint_.W;
      W_dot_ = bodies_[bodies_[i].parent_].joint_.W_dot;
      Eigen::VectorXd temp = Eigen::VectorXd::Ones(1); // for making double format the matrix format
      ga_[i-1] = (S_.transpose()*M_arti*S_).inverse()*
                 (temp*gf_[i-1]-S_.transpose()*M_arti*(S_dot_*gv_[i-1]+X_dot_.transpose()*W_+X_.transpose()*W_dot_)-S_.transpose()*b_arti) ;
      bodies_[i].joint_.W_dot = S_*ga_[i-1]+S_dot_*gv_[i-1]+X_dot_.transpose()*W_+X_.transpose()*W_dot_;
    }

    return ga_;
  }

private:
  std::vector<Body> bodies_;
  Eigen::VectorXd gc_;
  Eigen::VectorXd gv_;
  Eigen::VectorXd gf_;

};

/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, const Eigen::VectorXd& gf) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

  std::vector<Body> bodies;
  Body link1(Eigen::Vector3d{0., 0., 0.}, 0,
             GetInertiaMatrix(0., 0, 0, 0., 0, 0.), 0,
             Eigen::Vector3d::Zero(), Eigen::Matrix3d::Zero(), Eigen::Vector3d::Zero(),
             Body::Joint::Type::fixed);
  bodies.push_back(link1);

  Body link2(Eigen::Vector3d{0, 0, 0.2}, 1,
             GetInertiaMatrix(0.001, 0, 0, 0.001, 0, 0.001), 0,
             Eigen::Vector3d{0, 0, 0.3}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
             Body::Joint::Type::revolute);
  bodies.push_back(link2);

  Body link3(Eigen::Vector3d{0, 0, 0.2}, 1,
             GetInertiaMatrix(0.001, 0, 0, 0.001, 0, 0.001), 1,
             Eigen::Vector3d{0, 0, 0.3}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{1, 0, 0},
             Body::Joint::Type::revolute);
  bodies.push_back(link3);

  Body link4(Eigen::Vector3d{0, 0, 0.2}, 1,
             GetInertiaMatrix(0.001, 0, 0, 0.001, 0, 0.001), 2,
             Eigen::Vector3d{0, 0, 0.3}, Eigen::Matrix3d::Identity(), Eigen::Vector3d{0, 1, 0},
             Body::Joint::Type::prismatic);
  bodies.push_back(link4);

  ArticulatedSystem railab(bodies);
  railab.computeForwardKinematics(gc, gv, gf);
  railab.computeArticulatedVariable();

  return railab.computeGeneralizedAcceleration();
}