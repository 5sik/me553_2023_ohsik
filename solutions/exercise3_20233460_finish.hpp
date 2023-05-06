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
Eigen::Matrix3d QtoR(const Eigen::VectorXd &gc) {
  Eigen::Matrix3d Rot;
  Eigen::Vector4d q; // quaternion vector
  for (int i = 0; i < 4; i++) q[i] = gc[i + 3];

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

Eigen::Vector3d FR_HipToThigh, FR_FBaseToHip, FR_CalfToFoot, FR_ThighToCalf;
Eigen::Vector3d FL_HipToThigh, FL_FBaseToHip, FL_CalfToFoot, FL_ThighToCalf;
Eigen::Vector3d RR_HipToThigh, RR_FBaseToHip, RR_CalfToFoot, RR_ThighToCalf;
Eigen::Vector3d RL_HipToThigh, RL_FBaseToHip, RL_CalfToFoot, RL_ThighToCalf;
Eigen::Vector3d TrunkCOM;
Eigen::Vector3d ImuCOM;
Eigen::Vector3d TrunkImuCOM;
Eigen::Vector3d p1, p2, p3;
Eigen::Vector3d FR_HipCOM, FR_ThighCOM, FR_CalfCOM, FR_FootCOM, FR_CalfFootCOM;
Eigen::Vector3d FL_HipCOM, FL_ThighCOM, FL_CalfCOM, FL_FootCOM, FL_CalfFootCOM;
Eigen::Vector3d RR_HipCOM, RR_ThighCOM, RR_CalfCOM, RR_FootCOM, RR_CalfFootCOM;
Eigen::Vector3d RL_HipCOM, RL_ThighCOM, RL_CalfCOM, RL_FootCOM, RL_CalfFootCOM;
const double HipMass = 1.993, ThighMass = 0.639, CalfMass = 0.207, FootMass = 0.06, CalfFootMass = 0.207 + 0.06;
const double TrunkMass = 9.041, ImuMass = 0.001, TrunkImuMass = 9.042;


Eigen::Vector3d RootP;
Eigen::Matrix3d FR_CalfRot, FR_ThighRot, FR_HipRot;
Eigen::Matrix3d FL_CalfRot, FL_ThighRot, FL_HipRot;
Eigen::Matrix3d RR_CalfRot, RR_ThighRot, RR_HipRot;
Eigen::Matrix3d RL_CalfRot, RL_ThighRot, RL_HipRot;
Eigen::Matrix3d Roota;


Eigen::Matrix3d FR_Calf_I, FR_Thigh_I, FR_Hip_I, FR_Foot_I, FR_CalfFoot_I;
Eigen::Matrix3d FL_Calf_I, FL_Thigh_I, FL_Hip_I, FL_Foot_I, FL_CalfFoot_I;
Eigen::Matrix3d RR_Calf_I, RR_Thigh_I, RR_Hip_I, RR_Foot_I, RR_CalfFoot_I;
Eigen::Matrix3d RL_Calf_I, RL_Thigh_I, RL_Hip_I, RL_Foot_I, RL_CalfFoot_I;
Eigen::Matrix3d Trunk_I, Imu_I, TrunkImu_I;

void getCurrentState(const Eigen::VectorXd &gc) { ///get robot configuration parameter form URDF & GC
  p1 << 1, 0, 0;
  p2 << 0, 1, 0;
  p3 << 0, 1, 0;

  TrunkCOM << 0.008465, 0.004045, -0.000763;
  ImuCOM << 0, 0, 0;
  TrunkImuCOM = (TrunkMass * TrunkCOM + ImuMass * ImuCOM) / TrunkImuMass;

  FR_HipCOM << -0.022191, -0.015144, -1.5e-05;
  FR_ThighCOM << -0.005607, 0.003877, -0.048199;
  FR_CalfCOM << 0.002781, 6.3e-05, -0.142518;
  FR_FootCOM << 0, 0, -0.25;
  FR_CalfFootCOM = (CalfMass * FR_CalfCOM + FootMass * FR_FootCOM) / CalfFootMass;

  FL_HipCOM << -0.022191, 0.015144, -1.5e-05;
  FL_ThighCOM << -0.005607, -0.003877, -0.048199;
  FL_CalfCOM << 0.002781, 6.3e-05, -0.142518;
  FL_FootCOM << 0, 0, -0.25;
  FL_CalfFootCOM = (CalfMass * FL_CalfCOM + FootMass * FL_FootCOM) / CalfFootMass;

  RR_HipCOM << 0.022191, -0.015144, -1.5e-05;
  RR_ThighCOM << -0.005607, 0.003877, -0.048199;
  RR_CalfCOM << 0.002781, 6.3e-05, -0.142518;
  RR_FootCOM << 0, 0, -0.25;
  RR_CalfFootCOM = (CalfMass * RR_CalfCOM + FootMass * RR_FootCOM) / CalfFootMass;

  RL_HipCOM << 0.022191, 0.015144, -1.5e-05;
  RL_ThighCOM << -0.005607, -0.003877, -0.048199;
  RL_CalfCOM << 0.002781, 6.3e-05, -0.142518;
  RL_FootCOM << 0, 0, -0.25;
  RL_CalfFootCOM = (CalfMass * RL_CalfCOM + FootMass * RL_FootCOM) / CalfFootMass;

  ///joint rotation
  FR_HipRot = RotM("x", gc(7));
  FR_ThighRot = RotM("y", gc(8));
  FR_CalfRot = RotM("y", gc(9));
  FL_HipRot = RotM("x", gc(10));
  FL_ThighRot = RotM("y", 0); // Prismatic joint gc(11)
  FL_CalfRot = RotM("y", gc(12));
  RR_HipRot = RotM("x", gc(13));
  RR_ThighRot = RotM("y", gc(14));
  RR_CalfRot = RotM("y", gc(15));
  RL_HipRot = RotM("x", 0); // Prismatic joint gc(16)
  RL_ThighRot = RotM("y", 0); // Prismatic joint gc(17)
  RL_CalfRot = RotM("y", 0); // Prismatic joint gc(18)

  ///distance between two joints
  FR_FBaseToHip << 0.2399, -0.051, 0;
  FR_HipToThigh << 0, -0.083, 0;
  FR_ThighToCalf << 0, 0, -0.25;
  FR_CalfToFoot << 0, 0, -0.25;
  FL_FBaseToHip << 0.2399, 0.051, 0;
  FL_HipToThigh << 0, 0.083, 0;
  FL_ThighToCalf << 0, 0 +gc(11) , -0.25;
  FL_CalfToFoot << 0, 0, -0.25;
  RR_FBaseToHip << -0.2399, -0.051, 0;
  RR_HipToThigh << 0, -0.083, 0;
  RR_ThighToCalf << 0, 0, -0.25;
  RR_CalfToFoot << 0, 0, -0.25;
  RL_FBaseToHip << -0.2399, 0.051, 0;
  RL_HipToThigh << 0 +gc(16), 0.083, 0;
  RL_ThighToCalf << 0, 0+gc(17), -0.25;
  RL_CalfToFoot << 0, 0+gc(18), -0.25;
  RootP << gc(0), gc(1), gc(2);
  Roota = QtoR(gc);

  ///Inertia Tensor
  Eigen::VectorXd T, I;
  Eigen::VectorXd FR_H, FR_TH, FR_C, FR_F;
  Eigen::VectorXd FL_H, FL_TH, FL_C, FL_F;
  Eigen::VectorXd RR_H, RR_TH, RR_C, RR_F;
  Eigen::VectorXd RL_H, RL_TH, RL_C, RL_F;

  Trunk_I = GetInertiaMatrix(0.033260231, -0.000451628, 0.000487603, 0.16117211, 4.8356e-05, 0.17460442);
  Imu_I = GetInertiaMatrix(0.0001, 0, 0, 0.000001, 0, 0.0001);
  Eigen::Vector3d r1 = TrunkCOM - TrunkImuCOM;
  Eigen::Vector3d r2 = ImuCOM - TrunkImuCOM;
  TrunkImu_I =
      Trunk_I + Imu_I - TrunkMass * skewSymMat(r1) * skewSymMat(r1) - ImuMass * skewSymMat(r2) * skewSymMat(r2);

  FR_Hip_I = GetInertiaMatrix(0.002903894, 7.185e-05, -1.262e-06, 0.004907517, 1.75e-06, 0.005586944);
  FR_Thigh_I = GetInertiaMatrix(0.005666803, -3.597e-06, 0.000491446, 0.005847229, -1.0086e-05, 0.000369811);
  FR_Calf_I = GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05);
  FR_Foot_I = GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05);
  Eigen::Vector3d r3 = FR_CalfCOM - FR_CalfFootCOM;
  Eigen::Vector3d r4 = FR_FootCOM - FR_CalfFootCOM;
  FR_CalfFoot_I =
      FR_Calf_I + FR_Foot_I - CalfMass * skewSymMat(r3) * skewSymMat(r3) - FootMass * skewSymMat(r4) * skewSymMat(r4);

  FL_Hip_I = GetInertiaMatrix(0.002903894, -7.185e-05, -1.262e-06, 0.004907517, -1.75e-06, 0.005586944);
  FL_Thigh_I = GetInertiaMatrix(0.005666803, 3.597e-06, 0.000491446, 0.005847229, 1.0086e-05, 0.000369811);
  FL_Calf_I = GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05);
  FL_Foot_I = GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05);
  Eigen::Vector3d r5 = FL_CalfCOM - FL_CalfFootCOM;
  Eigen::Vector3d r6 = FL_FootCOM - FL_CalfFootCOM;
  FL_CalfFoot_I =
      FL_Calf_I + FL_Foot_I - CalfMass * skewSymMat(r5) * skewSymMat(r5) - FootMass * skewSymMat(r6) * skewSymMat(r6);

  RR_Hip_I = GetInertiaMatrix(0.002903894, -7.185e-05, 1.262e-06, 0.004907517, 1.75e-06, 0.005586944);
  RR_Thigh_I = GetInertiaMatrix(0.005666803, -3.597e-06, 0.000491446, 0.005847229, -1.0086e-05, 0.000369811);
  RR_Calf_I = GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05);
  RR_Foot_I = GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05);
  Eigen::Vector3d r7 = RR_CalfCOM - RR_CalfFootCOM;
  Eigen::Vector3d r8 = RR_FootCOM - RR_CalfFootCOM;
  RR_CalfFoot_I =
      RR_Calf_I + RR_Foot_I - CalfMass * skewSymMat(r7) * skewSymMat(r7) - FootMass * skewSymMat(r8) * skewSymMat(r8);


  RL_Hip_I = GetInertiaMatrix(0.002903894, 7.185e-05, 1.262e-06, 0.004907517, -1.75e-06, 0.005586944);
  RL_Thigh_I = GetInertiaMatrix(0.005666803, 3.597e-06, 0.000491446, 0.005847229, 1.0086e-05, 0.000369811);
  RL_Calf_I = GetInertiaMatrix(0.006341369, -3e-09, -8.7951e-05, 0.006355157, -1.336e-06, 3.9188e-05);
  RL_Foot_I = GetInertiaMatrix(1.6854e-05, 0.0, 0.0, 1.6854e-05, 0.0, 1.6854e-05);
  Eigen::Vector3d r9 = RL_CalfCOM - RL_CalfFootCOM;
  Eigen::Vector3d r10 = RL_FootCOM - RL_CalfFootCOM;
  RL_CalfFoot_I =
      RL_Calf_I + RL_Foot_I - CalfMass * skewSymMat(r9) * skewSymMat(r9) - FootMass * skewSymMat(r10) * skewSymMat(r10);

}

double Mass(int index) { // TrunkImuMass & CalfFootMass are composite bodies.
  if (index == 0) {
    return TrunkImuMass;
  } else if (index == 1 || index == 4 || index == 7 || index == 10) {
    return HipMass;
  } else if (index == 2 || index == 5 || index == 8 || index == 11) {
    return ThighMass;
  } else if (index == 3 || index == 6 || index == 9 || index == 12) {
    return CalfFootMass;
  } else {
    std::cout << "invalid mass index";
    return 0;
  }
}

/// joint Position with repsect to world frame
Eigen::Vector3d getJointAbsPos(const Eigen::VectorXd &gc, int index1, int index2) {
  //index1 : FR #1 FL #2 RR #3 RL #4
  //index2 : Hip #1 Thigh #2 Calf #3 Foot #4

  Eigen::Vector3d FBaseToHip;
  Eigen::Vector3d HipToThigh;
  Eigen::Vector3d CalfToFoot;
  Eigen::Vector3d ThighToCalf;

  Eigen::Matrix3d HipRot;
  Eigen::Matrix3d ThighRot;
  Eigen::Matrix3d CalfRot;

  if (index1 == 0) {
    return RootP;
  } else if (index1 == 1) {
    FBaseToHip = FR_FBaseToHip;
    HipToThigh = FR_HipToThigh;
    CalfToFoot = FR_CalfToFoot;
    ThighToCalf = FR_ThighToCalf;
    HipRot = FR_HipRot;
    ThighRot = FR_ThighRot;
    CalfRot = FR_CalfRot;
  } else if (index1 == 2) {
    FBaseToHip = FL_FBaseToHip;
    HipToThigh = FL_HipToThigh;
    CalfToFoot = FL_CalfToFoot;
    ThighToCalf = FL_ThighToCalf;
    HipRot = FL_HipRot;
    ThighRot = FL_ThighRot;
    CalfRot = FL_CalfRot;
  } else if (index1 == 3) {
    FBaseToHip = RR_FBaseToHip;
    HipToThigh = RR_HipToThigh;
    CalfToFoot = RR_CalfToFoot;
    ThighToCalf = RR_ThighToCalf;
    HipRot = RR_HipRot;
    ThighRot = RR_ThighRot;
    CalfRot = RR_CalfRot;
  } else if (index1 == 4) {
    FBaseToHip = RL_FBaseToHip;
    HipToThigh = RL_HipToThigh;
    CalfToFoot = RL_CalfToFoot;
    ThighToCalf = RL_ThighToCalf;
    HipRot = RL_HipRot;
    ThighRot = RL_ThighRot;
    CalfRot = RL_CalfRot;
  } else {
    return Eigen::Vector3d::Zero();
  }

  if (index2 == 1) { ///hip
    return RootP + Roota * FBaseToHip;
  } else if (index2 == 2) { ///thigh
    return RootP + Roota * (FBaseToHip + HipRot * HipToThigh);
  } else if (index2 == 3) { ///calf
    return RootP + Roota * (FBaseToHip + HipRot * (HipToThigh + ThighRot * ThighToCalf));
  } else if (index2 == 4) { ///foot
    return RootP + Roota * (FBaseToHip + HipRot * (HipToThigh + ThighRot * (ThighToCalf + CalfRot * CalfToFoot)));
  } else {
    return RootP;
//    return Eigen::Vector3d::Zero();
  }
}

/// COM Position of link with repsect to world frame
Eigen::Vector3d getCOMAbsPos(const Eigen::VectorXd &gc, int index1, int index2) {
  //index1 : Trunk #0 FR #1 FL #2 RR #3 RL #4
  //index2 : Hip #1 Thigh #2 CalfFoot #3

  Eigen::Vector3d HipCOM;
  Eigen::Vector3d ThighCOM;
  Eigen::Vector3d CalfFootCOM;

  Eigen::Matrix3d HipRot;
  Eigen::Matrix3d ThighRot;
  Eigen::Matrix3d CalfRot;

  if (index1 == 0) {
    return RootP + Roota * TrunkImuCOM;
  } else if (index1 == 1) {
    HipRot = FR_HipRot;
    ThighRot = FR_ThighRot;
    CalfRot = FR_CalfRot;
    HipCOM = FR_HipCOM;
    ThighCOM = FR_ThighCOM;
    CalfFootCOM = FR_CalfFootCOM;
  } else if (index1 == 2) {
    HipRot = FL_HipRot;
    ThighRot = FL_ThighRot;
    CalfRot = FL_CalfRot;
    HipCOM = FL_HipCOM;
    ThighCOM = FL_ThighCOM + Eigen::Vector3d {0, gc(11),0};
    CalfFootCOM = FL_CalfFootCOM;
  } else if (index1 == 3) {
    HipRot = RR_HipRot;
    ThighRot = RR_ThighRot;
    CalfRot = RR_CalfRot;
    HipCOM = RR_HipCOM;
    ThighCOM = RR_ThighCOM;
    CalfFootCOM = RR_CalfFootCOM;
  } else if (index1 == 4) {
    HipRot = RL_HipRot;
    ThighRot = RL_ThighRot;
    CalfRot = RL_CalfRot;
    HipCOM = RL_HipCOM + Eigen::Vector3d {gc(16),0,0};
    ThighCOM = RL_ThighCOM + Eigen::Vector3d {0, gc(17),0};
    CalfFootCOM = RL_CalfFootCOM + Eigen::Vector3d {0, gc(18),0};
  } else {
    return Eigen::Vector3d::Zero();
  }

  if (index2 == 1) { ///hip
    return getJointAbsPos(gc, index1, index2) + Roota * HipRot * HipCOM; // 해당 joint까지 오고 거기서 COM만큼 더감.
  } else if (index2 == 2) { ///thigh
    return getJointAbsPos(gc, index1, index2) + Roota * HipRot * ThighRot * ThighCOM; // 해당 joint까지 오고 거기서 COM만큼 더감.
  } else if (index2 == 3) { ///calf
    return getJointAbsPos(gc, index1, index2) +
           Roota * HipRot * ThighRot * CalfRot * CalfFootCOM; // 해당 joint까지 오고 거기서 COM만큼 더감.
  } else {
    return Eigen::Vector3d::Zero();
  }
}

/// Axis of each joints w.r.t world frame after all rotation
Eigen::Vector3d getRotAxis_w(int index1, int index2) {
  //index1 : FR #1 FL #2 RR #3 RL #4
  //index2 : Hip #1 Thigh #2 Calf #3 Foot #4
  //Prismatic (2,2/4,1/4,2/4,3)
  if (index1 == 1) {

    if (index2 == 1) return Roota * p1;
    else if (index2 == 2) return Roota * FR_HipRot * p2;
    else if (index2 == 3) return Roota * FR_HipRot * FR_ThighRot * p3;
    else {
      std::cout << "invalid rot index";
      return Eigen::Vector3d::Zero();
    }

  } else if (index1 == 2) {

    if (index2 == 1) return Roota * p1;
    else if (index2 == 2) return Roota * FL_HipRot * p2;
    else if (index2 == 3) return Roota * FL_HipRot * FL_ThighRot * p3;
    else {
      std::cout << "invalid rot index";
      return Eigen::Vector3d::Zero();
    }

  } else if (index1 == 3) {

    if (index2 == 1) return Roota * p1;
    else if (index2 == 2) return Roota * RR_HipRot * p2;
    else if (index2 == 3) return Roota * RR_HipRot * RR_ThighRot * p3;
    else {
      std::cout << "invalid rot index";
      return Eigen::Vector3d::Zero();
    }

  } else if (index1 == 4) {

    if (index2 == 1) return Roota * p1;
    else if (index2 == 2) return Roota * RL_HipRot * p2;
    else if (index2 == 3) return Roota * RL_HipRot * RL_ThighRot * p3;
    else {
      std::cout << "invalid rot index";
      return Eigen::Vector3d::Zero();
    }

  } else {
    std::cout << "invalid rot index";
    return Eigen::Vector3d::Zero();
  }
}

/// Get Inertia changing from the each Frame to the world Frame
Eigen::Matrix3d getInertiaTensor_w(const Eigen::VectorXd &gc, int index) {
  /// For Trunk /// I(B) = R_(B/A) * I(A) * R_(B/A)^T /// B/A 는 from Frame(A) to Frame(B)
  if (index == 0) {
    return Roota * TrunkImu_I * Roota.transpose();
  }
    /// For FR
  else if (index == 1) {
    return (Roota * FR_HipRot) * FR_Hip_I * (Roota * FR_HipRot).transpose();
  } else if (index == 2) {
    return (Roota * FR_HipRot * FR_ThighRot) * FR_Thigh_I * (Roota * FR_HipRot * FR_ThighRot).transpose();
  } else if (index == 3) {
    return (Roota * FR_HipRot * FR_ThighRot * FR_CalfRot) * FR_CalfFoot_I *
           (Roota * FR_HipRot * FR_ThighRot * FR_CalfRot).transpose();
  }/// For FL
  else if (index == 4) {
    return (Roota * FL_HipRot) * FL_Hip_I * (Roota * FL_HipRot).transpose();
  } else if (index == 5) {
    return (Roota * FL_HipRot * FL_ThighRot) * FL_Thigh_I * (Roota * FL_HipRot * FL_ThighRot).transpose();
  } else if (index == 6) {
    return (Roota * FL_HipRot * FL_ThighRot * FL_CalfRot) * FL_CalfFoot_I *
           (Roota * FL_HipRot * FL_ThighRot * FL_CalfRot).transpose();
  }/// For RR
  else if (index == 7) {
    return (Roota * RR_HipRot) * RR_Hip_I * (Roota * RR_HipRot).transpose();
  } else if (index == 8) {
    return (Roota * RR_HipRot * RR_ThighRot) * RR_Thigh_I * (Roota * RR_HipRot * RR_ThighRot).transpose();
  } else if (index == 9) {
    return (Roota * RR_HipRot * RR_ThighRot * RR_CalfRot) * RR_CalfFoot_I *
           (Roota * RR_HipRot * RR_ThighRot * RR_CalfRot).transpose();
  }/// For RL
  else if (index == 10) {
    return (Roota * RL_HipRot) * RL_Hip_I * (Roota * RL_HipRot).transpose();
  } else if (index == 11) {
    return (Roota * RL_HipRot * RL_ThighRot) * RL_Thigh_I * (Roota * RL_HipRot * RL_ThighRot).transpose();
  } else if (index == 12) {
    return (Roota * RL_HipRot * RL_ThighRot * RL_CalfRot) * RL_CalfFoot_I *
           (Roota * RL_HipRot * RL_ThighRot * RL_CalfRot).transpose();
  } else {
    std::cout << "invalid joint index";
    return Eigen::Matrix3d::Zero();
  }
}

/// Get Positional Jacobian at each Joints which are hip,thigh and calf joints
Eigen::MatrixXd getJointPosJ(const Eigen::VectorXd &gc, int index) {
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); /// For pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d r1, r2, r3;

  switch (index) {
    case 0: /// Base
      J << I3, Z3, Z3, Z3, Z3, Z3;
      break;
    case 1: /// FR_H : 1-1
      r1 = getJointAbsPos(gc,1,1)- getJointAbsPos(gc,0,0); // 몸통에서 FR_H까지 거리
      J << I3,-skewSymMat(r1),Z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      r1 = getJointAbsPos(gc,1,2)- getJointAbsPos(gc,0,0);
      r2 = getJointAbsPos(gc,1,2)- getJointAbsPos(gc,1,1);
      J<< I3,-1*skewSymMat(r1), -1*skewSymMat(r2)*getRotAxis_w(1,1),z3,z3,Z3,Z3,Z3; // P_vector is written w.r.t. the world frame.
    case 3: ///FR_C 1-3
      r1 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 1, 1);
      r3 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 1, 2);
      J << I3, -1*skewSymMat(r1), -1*skewSymMat(r2)*getRotAxis_w(1,1), -1*skewSymMat(r3)* getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      r1 = getJointAbsPos(gc, 2, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skewSymMat(r1),Z3,Z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      r1 = getJointAbsPos(gc, 2, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 2, 2) - getJointAbsPos(gc, 2, 1);
      J << I3, -1*skewSymMat(r1),Z3, -1*skewSymMat(r2)*getRotAxis_w(2,1),z3,z3,Z3,Z3;
      break;
    case 6: ///FL_C 2-3
      r1 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 2, 1);
      r3 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 2, 2);
      J << I3, -1*skewSymMat(r1),Z3, -1*skewSymMat(r2)*getRotAxis_w(2,1), getRotAxis_w(2,2),z3,Z3,Z3; // thigh_prismatic 적용
      break;
    case 7: ///RR_H 3-1
      r1 = getJointAbsPos(gc, 3, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skewSymMat(r1),Z3,Z3,Z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      r1 = getJointAbsPos(gc, 3, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 3, 2) - getJointAbsPos(gc, 3, 1);
      J << I3, -1*skewSymMat(r1),Z3,Z3, -1*skewSymMat(r2)*getRotAxis_w(3,1),z3,z3,Z3;
      break;
    case 9: ///RR_C 3-3
      r1 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 3, 1);
      r3 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 3, 2);
      J << I3, -1*skewSymMat(r1),Z3,Z3, -1*skewSymMat(r2)*getRotAxis_w(3,1), -1*skewSymMat(r3)* getRotAxis_w(3,2),z3,Z3;
      break;
    case 10: ///RL_H 4-1
      r1 = getJointAbsPos(gc, 4, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skewSymMat(r1),Z3,Z3,Z3,Z3;
      break;
    case 11: ///RL_TH 4-2
      r1 = getJointAbsPos(gc, 4, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 4, 2) - getJointAbsPos(gc, 4, 1); // prismatic이라 상대적 거리는 안씀
      J << I3, -1*skewSymMat(r1),Z3,Z3,Z3, getRotAxis_w(4,1),z3,z3;
      break;
    case 12: ///RL_C 4-3
      r1 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 4, 1);
      r3 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 4, 2);
      J << I3, -1*skewSymMat(r1),Z3,Z3,Z3, getRotAxis_w(4,1), getRotAxis_w(4,2),z3;
      break;
    default:
      std::cout<<"invalid joint index";
      break;
  }
  return J;
}

/// Get Positional Jacobian at COM Of each Links which are hip,thigh and calf Links
Eigen::MatrixXd getLinkCOMPosJ(const Eigen::VectorXd &gc, int index){
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d r1, r2, r3, r4;

  switch(index) {
    case 0: ///TrunkImu
      r1 = getCOMAbsPos(gc, 0, 0) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skewSymMat(r1),Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      r1 = getCOMAbsPos(gc, 1, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc,1,1) - getJointAbsPos(gc,1,1);
      J << I3, -1 * skewSymMat(r1),-1*skewSymMat(r2)*getRotAxis_w(1,1), z3,z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      r1 = getCOMAbsPos(gc, 1, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 1);
      r3 = getCOMAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 2);
      J << I3, -1*skewSymMat(r1), -1*skewSymMat(r2)*getRotAxis_w(1,1),-1*skewSymMat(r3)* getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C+F 1-3
      r1 = getCOMAbsPos(gc, 1, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 1);
      r3 = getCOMAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 2);
      r4 = getCOMAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 3);
      J << I3, -1*skewSymMat(r1), -1*skewSymMat(r2)*getRotAxis_w(1,1), -1*skewSymMat(r3)* getRotAxis_w(1,2),-1*skewSymMat(r4)* getRotAxis_w(1,3),Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      r1 = getCOMAbsPos(gc, 2, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc,2,1) - getJointAbsPos(gc,2,1);
      J << I3, -1 * skewSymMat(r1), Z3, -1*skewSymMat(r2)*getRotAxis_w(2,1), z3,z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      r1 = getCOMAbsPos(gc, 2, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 1);
      r3 = getCOMAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 2);
      J << I3, -1*skewSymMat(r1), Z3, -1*skewSymMat(r2)*getRotAxis_w(2,1),getRotAxis_w(2,2),z3,Z3,Z3;
      break;
    case 6: ///FL_C+F 2-3
      r1 = getCOMAbsPos(gc, 2, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 1);
      r3 = getCOMAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 2);
      r4 = getCOMAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 3);
      J << I3, -1*skewSymMat(r1), Z3, -1*skewSymMat(r2)*getRotAxis_w(2,1), getRotAxis_w(2,2),-1*skewSymMat(r4)* getRotAxis_w(2,3),Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      r1 = getCOMAbsPos(gc, 3, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc,3,1) - getJointAbsPos(gc,3,1);
      J << I3, -1 * skewSymMat(r1), Z3, Z3,-1*skewSymMat(r2)*getRotAxis_w(3,1), z3,z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      r1 = getCOMAbsPos(gc, 3, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 1);
      r3 = getCOMAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 2);
      J << I3, -1*skewSymMat(r1), Z3, Z3, -1*skewSymMat(r2)*getRotAxis_w(3,1),-1*skewSymMat(r3)* getRotAxis_w(3,2),z3,Z3;
      break;
    case 9: ///RR_C+F 3-3
      r1 = getCOMAbsPos(gc, 3, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 1);
      r3 = getCOMAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 2);
      r4 = getCOMAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 3);
      J << I3, -1*skewSymMat(r1), Z3, Z3, -1*skewSymMat(r2)*getRotAxis_w(3,1), -1*skewSymMat(r3)* getRotAxis_w(3,2),-1*skewSymMat(r4)* getRotAxis_w(3,3),Z3;
      break;
    case 10: ///RL_H 4-1
      r1 = getCOMAbsPos(gc, 4, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc,4,1) - getJointAbsPos(gc,4,1);
      J << I3, -1 * skewSymMat(r1), Z3, Z3,Z3, getRotAxis_w(4,1), z3,z3;
      break;
    case 11: ///RL_TH 4-2
      r1 = getCOMAbsPos(gc, 4, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 1);
      r3 = getCOMAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 2);
      J << I3, -1*skewSymMat(r1), Z3, Z3,Z3, getRotAxis_w(4,1),getRotAxis_w(4,2),z3;
      break;
    case 12: ///RL_C+F 4-3
      r1 = getCOMAbsPos(gc, 4, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getCOMAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 1);
      r3 = getCOMAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 2);
      r4 = getCOMAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 3);
      J << I3, -1*skewSymMat(r1), Z3, Z3,Z3, getRotAxis_w(4,1), getRotAxis_w(4,2),getRotAxis_w(4,3);
      break;
    default:
      std::cout<<"invalid body index";
      break;
  }
  return J;
}

/// Get Angular Jacobian at each Joints which are hip,thigh and calf joints
Eigen::MatrixXd getJointAngJ(const Eigen::VectorXd &gc, int index) {
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  switch (index) {
    case 0: ///Base
      J << Z3, I3, Z3, Z3, Z3, Z3;
      break;
    case 1: ///FR_H 1-1
      J << Z3, I3, getRotAxis_w(1, 1), z3, z3, Z3, Z3, Z3;
      break;
    case 2: ///FR_TH 1-2
      J << Z3, I3, getRotAxis_w(1, 1), getRotAxis_w(1, 2), z3, Z3, Z3, Z3;
      break;
    case 3: ///FR_C 1-3
      J << Z3, I3, getRotAxis_w(1, 1), getRotAxis_w(1, 2), getRotAxis_w(1, 3), Z3, Z3, Z3;
      break;
    case 4: ///FL_H 2-1
      J << Z3, I3, Z3, getRotAxis_w(2, 1), z3, z3, Z3, Z3;
      break;
    case 5: ///FL_TH 2-2
      J << Z3, I3, Z3, getRotAxis_w(2, 1), z3, z3, Z3, Z3;
      break;
    case 6: ///FL_C 2-3
      J << Z3, I3, Z3, getRotAxis_w(2, 1), z3, getRotAxis_w(2, 3), Z3, Z3;
      break;
    case 7: ///RR_H 3-1
      J << Z3, I3, Z3, Z3, getRotAxis_w(3, 1), z3, z3, Z3;
      break;
    case 8: ///RR_TH 3-2
      J << Z3, I3, Z3, Z3, getRotAxis_w(3, 1), getRotAxis_w(3, 2), z3, Z3;
      break;
    case 9: ///RR_C 3-3
      J << Z3, I3, Z3, Z3, getRotAxis_w(3, 1), getRotAxis_w(3, 2), getRotAxis_w(3, 3), Z3;
      break;
    case 10: ///RL_H 4-1
      J << Z3, I3, Z3, Z3, Z3, Z3;
      break;
    case 11: ///RL_TH 4-2
      J << Z3, I3, Z3, Z3, Z3, Z3;
      break;
    case 12: ///RL_C 4-3
      J << Z3, I3, Z3, Z3, Z3, Z3;
      break;
    default:
      std::cout << "invalid joint index";
      break;
  }
  return J;
}

/// Get the Velocity of each joint, Vector3d = {X,Y,Z}
Eigen::Vector3d getFrameVel(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::Vector3d V;

  V = getJointPosJ(gc,index) * gv;
  return V;
}

/// Get the Velocity of COM of the Link, Vector3d = {X,Y,Z}
Eigen::Vector3d getCOMVel(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::Vector3d VC;

  VC = getLinkCOMPosJ(gc,index) * gv;
  return VC;
}

/// Get the Angular Velocity of each joint, Vector3d = {W_x, W_y, W_z}
Eigen::Vector3d getFrameAngVel(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::Vector3d W;

  W = getJointAngJ(gc,index) * gv;
  return W;
}


inline Eigen::MatrixXd getMassMatrix(const Eigen::VectorXd &gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

  getCurrentState(gc);
  Eigen::MatrixXd M(18,18);
  M.setZero();
  for (int i=0; i<13; i++){

    M += ((getLinkCOMPosJ(gc,i).transpose()*Mass(i)*getLinkCOMPosJ(gc,i)) + (getJointAngJ(gc,i).transpose()*getInertiaTensor_w(gc,i)*getJointAngJ(gc,i)));
  }

  return M;
}