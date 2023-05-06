#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include "raisim/math.hpp"
#include <iostream>
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

////// Get Jacobian cosisting Jaco.Pos + Jaco. Ang
std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::MatrixXd>> getJacobian(const Eigen::VectorXd& gc){
  Eigen::Vector3d position;
  Eigen::Vector3d axis;
  Eigen::Matrix3d rot, temp;
  std::vector<Eigen::MatrixXd> JP;  // jacobian for linear Vel.
  std::vector<Eigen::MatrixXd> JA; // jacobian for Angular Vel.

  for (int link=0 ; link<18 ; link++){ // 1-몸통, 1-imu(mass있음), 12-힙,피치,니, 4-foot
    JP.push_back(Eigen::MatrixXd::Zero(3,18));
    JA.push_back(Eigen::MatrixXd::Zero(3,18));
  }

  // For trunk, imu link
  rot = QtoR(gc.segment(3,4));
  position << 0.008465, 0.004045, -0.000763; // distance between floating_base_joint and trunk_link_COM
  position << rot * position; // trunk com쪽이랑 quaternion 고려한 pos
  // 이제 만들어놓은 JP, JA채우기 여기서는 trunk joint frame 기준으로
  JP.at(0).middleCols(0,3) = rot.transpose(); // world frame으로 보는게 아니라서 quaternion 만큼 돌아간거 다시 되돌려 놔야함.
  temp = skewSymMat(position);
  JP.at(0).middleCols(3,3) = -rot.transpose() * temp; // -는 외적 순서를 바꿨기 때문에 ( jaco * Gv 계산 위해 )
  JA.at(0).middleCols(3,3) = rot.transpose();
  // 이제 만들어놓은 JP, JA채우기 여기서는 imu joint frame 기준으로
  JP.at(1).middleCols(0,3) = rot.transpose();  // for IMU , 여기도 trunk joint에 의해 돌아간거 회복
  JA.at(1).middleCols(3,3) = rot.transpose();  // for IMU

  for (int leg=0; leg<4; leg++){ // leg = 0:FR/ 1:FL/ 2:RR/ 3: RL
    double front = 1 - 2 * (leg / 2); // front = 1 , hind(Rear) = -1
                                      // (leg/2)는 int로써 leg = 0,1 일때 front (=1) / leg=2,3 일때 Rear(=-1)
    double right = 1 - 2 * (leg & 2); // right = 1 , legt = -1
                                      // (leg&2)는 int로써 leg = 0,2 일때 right (=1) / leg=1,3 일때 left(=-1)
    std::vector<Eigen::Vector3d> linkPos; // 순서는 foot->calf->thigh->hip
    std::vector<Eigen::Matrix3d> linkRot;


    // w.r.t the calf_joint ///RL_calf는 prismatic임 조심조심
    linkPos.push_back(Eigen::Vector3d{0,0,-0.25}); // calf_joint 에서 foot_joint까지 거리
    linkRot.push_back(Eigen::Matrix3d::Identity()); // foot joint는 따로 돌아가는게 없어서 identity
    linkPos.push_back(Eigen::Vector3d{0.002781, 6.3e-05, -0.142518}); // calf_joint에서 calf com 까지 거리 & 다리 4짝 똑같음
    linkRot.push_back(Eigen::Matrix3d::Identity()); // calf COM은 따로 돌아가는게 없어서 identity
    axis << 0,1,0;

/////이 for문 뜯어 보기 안에 들어가는 의미 분석
    for (int link=0; link <linkPos.size() ; link++){  // linkPos.size는 foot이랑 calf_link pos하나씩 들어가서 2개임
      // (RL_calf)  leg = 3 prismatic
      if(leg == 3){ /// 1-몸통, 1-imu(mass있음), 12-힙,피치,니, 4-foot
        JP.at(2+4*leg+3-link).col(6+3*leg+2)=linkRot.at(link).transpose() * axis;
        //(2+4*leg+3-link)는 2+4*leg+3 이것은 기본적으로 foot but link가 하나씩 늘수록 foot->calf->thigh->hip 순서로 감
        //(6+3*leg+2)는 leg 별로 theta_dot_calf를 뜻함 in GV에서
        JA.at(2+4*leg+3-link).col(6+3*leg+2)=Eigen::Vector3d::Zero(); // O_3x1 짜리 들어와야함 due to prismatic joint
      }
      else { ///(RL_calf가 아니면 = revolute 표현)
        JP.at(2+4*leg+3-link).col(6+3*leg+2)=linkRot.at(link).transpose() * axis.cross(linkPos.at(link));
        JA.at(2+4*leg+3-link).col(6+3*leg+2)=linkRot.at(link).transpose() * axis;
      }
    }

    // w.r.t the thigh_joint // FL & RL thigh 는 prismatic임 조심조심
    rot = RotM(axis(0) == 1 ? "x" : axis(1) == 1 ? "y" : "z",7 + 3 * leg +2); // calf_joint에서 돌아가는거 check
    position << 0,0,-0.25; // thigh_joint 에서 calf_joint 까지 거리

    for (int link=0; link<linkPos.size(); link++){ // foot이랑 calf_link 가 calf_joint에 의해 돌아간거 계산
      linkPos.at(link) << rot * linkPos.at(link) + position; // calf 돌아간거 고려해서 컨테이너에 새로 계산해서 집어넣음
      linkRot.at(link) << rot * linkRot.at(link);
    }
    linkPos.push_back(Eigen::Vector3d{-0.005607, right * 0.003877, -0.048199});
    linkRot.push_back(Eigen::Matrix3d::Identity());
    axis<< 0,1,0;

    for (int link=0; link<linkPos.size(); link++){
      // FL_thigh leg = 1 &  RL_thigh leg = 3 은 prismatic joint
      if (( leg == 1 ) || ( leg == 3 )) {
        JP.at(2 + 4 * leg + 3 - link).col(6 + 3 * leg + 1) = linkRot.at(link).transpose() * axis;
        JA.at(2 + 4 * leg + 3 - link).col(6 + 3 * leg + 1) = Eigen::Vector3d::Zero();
      }
      else {
        JP.at(2+4*leg+3-link).col(6+3*leg+1)=linkRot.at(link).transpose() * axis.cross(linkPos.at(link));
        JA.at(2+4*leg+3-link).col(6+3*leg+1)=linkRot.at(link).transpose() * axis;
      }
    }

    // w.r.t the Hip_joint // RL hip 는 prismatic임 조심조심
    rot = RotM(axis(0) == 1 ? "x" : axis(1) == 1 ? "y" : "z",7 + 3 * leg +1); // thigh_joint에서 돌아가는거 check
    position << 0,-right*0.083,0; // Hip_joint 에서 thigh_joint 까지 거리

    for (int link=0; link<linkPos.size(); link++){ // thigh_joint에 의해 돌아간거 계산
      linkPos.at(link) << rot * linkPos.at(link) + position; // thigh 돌아간거 고려해서 컨테이너에 새로 계산해서 집어넣음
      linkRot.at(link) << rot * linkRot.at(link);
    }

    linkPos.push_back(Eigen::Vector3d{-front*0.022191, -right*0.015144, -1.5e-05});
    linkRot.push_back(Eigen::Matrix3d::Identity());
    axis << 1,0,0;

    for (int link=0; link<linkPos.size(); link++){
      // RL_hip leg = 3 은 prismatic joint
      if ( leg == 3 ){
        JP.at(2+4*leg+3-link).col(6+3*leg+0)=linkRot.at(link).transpose() * axis;
        JA.at(2+4*leg+3-link).col(6+3*leg+0)=Eigen::Vector3d::Zero();
      }
      else {
        JP.at(2+4*leg+3-link).col(6+3*leg+0)=linkRot.at(link).transpose() * axis.cross(linkPos.at(link));
        JA.at(2+4*leg+3-link).col(6+3*leg+0)=linkRot.at(link).transpose() * axis;
      }
    }

    // w.r.t the Body
    rot = RotM(axis(0) == 1 ? "x" : axis(1) == 1 ? "y" : "z",7 + 3 * leg +0); // hip_joint에서 돌아가는거 check
    position << front*0.2399, -right*0.051, 0; // Hip_joint 에서 thigh_joint 까지 거리

    for(int link=0; link<linkPos.size();link++){
      linkPos.at(link) << rot * linkPos.at(link) + position; // hip 돌아간거 고려해서 컨테이너에 새로 계산해서 집어넣음
      linkRot.at(link) << rot * linkRot.at(link);
    }

    // w.r.t the world
    rot = QtoR(gc.segment(3,4));
    for(int link=0; link<linkPos.size();link++){
      linkPos.at(link) << rot * linkPos.at(link); // body, quaternion으로 인해 돌아간거 고려해서 컨테이너에 새로 계산해서 집어넣음
      linkRot.at(link) << rot * linkRot.at(link);
    }

    for(int link=0; link<linkPos.size(); link++){
      JP.at(2+4*leg+3-link).middleCols(0,3) = linkRot.at(link).transpose();
      temp = skewSymMat(linkPos.at(link));
      JP.at(2+4*leg+3-link).middleCols(3,3) = -linkRot.at(link).transpose() * temp;
      JA.at(2+4*leg+3-link).middleCols(3,3) = -linkRot.at(link).transpose();
    }
  }
 return {JP, JA};
}
///// Make the Inertia Matrix using the half of the entries
Eigen::Matrix3d GetInertiaMatrix (double ixx,double ixy,double ixz,double iyy,double iyz,double izz){
  Eigen::Matrix3d I;
  I << ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz;
  return I;
}

inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::MatrixXd M;
  M.setZero(18,18);
  auto [JPs, JAs] = getJacobian(gc);
  Eigen::MatrixXd JP, JA;
  double m; // mass
  Eigen::Matrix3d I; // Inertia Matrix

  // Trunk
  JP = JPs.at(0); JA = JAs.at(0);
  m = 9.041; I << GetInertiaMatrix(0.033260231,-0.000451628,0.000487603,0.16117211,4.8356e-05,0.17460442);
  M += JP.transpose()*m*JP + JA.transpose()*I*JA;

  // IMU
  JP = JPs.at(1); JA = JAs.at(1);
  m = 0.001; I << GetInertiaMatrix(0.0001,0,0,0.000001,0,0.0001);
  M += JP.transpose()*m*JP + JA.transpose()*I*JA;

  for(int leg=0 ; leg<4 ; leg++){
    double front = 1 - 2 * (leg / 2);
    double right = 1 - 2 * (leg % 2);

    JP = JPs.at(2+4*leg+0); JA = JAs.at(2+4*leg+0);//for Hip
    m = 1.993; I << GetInertiaMatrix(0.002903894,front*right*7.185e-05,-front*1.262e-06,0.004907517,right*1.75e-06,0.005586944);
    M += JP.transpose()*m*JP + JA.transpose()*I*JA;

    JP = JPs.at(2+4*leg+1); JA = JAs.at(2+4*leg+1);//for Thigh
    m = 0.639; I << GetInertiaMatrix(0.005666803,-right*3.597e-06,0.000491446,0.005847229,-right*1.0086e-05,0.000369811);
    M += JP.transpose()*m*JP + JA.transpose()*I*JA;

    JP = JPs.at(2+4*leg+2); JA = JAs.at(2+4*leg+2);//for calf
    m = 0.207; I << GetInertiaMatrix(0.006341369,-3e-09,-8.7951e-05,0.006355157,-1.336e-06,3.9188e-05);
    M += JP.transpose()*m*JP + JA.transpose()*I*JA;

    JP = JPs.at(2+4*leg+3); JA = JAs.at(2+4*leg+3);//for foot
    m = 0.06; I << GetInertiaMatrix(9.6e-06,0.0,0.0,9.6e-06,0.0,9.6e-06);
    M += JP.transpose()*m*JP + JA.transpose()*I*JA;

  }

  return M;
}