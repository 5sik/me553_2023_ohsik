#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
Eigen::Vector3d CalfToFoot, ThighToCalf;
Eigen::Vector3d FR_HipToThigh, FR_FBaseToHip;
Eigen::Vector3d FL_HipToThigh, FL_FBaseToHip;
Eigen::Vector3d RR_HipToThigh, RR_FBaseToHip;
Eigen::Vector3d RL_HipToThigh, RL_FBaseToHip;
Eigen::Vector3d TrunkCOM;
Eigen::Vector3d ImuCOM;
Eigen::Vector3d TrunkImuCOM;
Eigen::Vector3d p1, p2, p3;
Eigen::Vector3d FR_HipCOM, FR_ThighCOM, FR_CalfCOM, FR_FootCOM, FR_CalfFootCOM;
Eigen::Vector3d FL_HipCOM, FL_ThighCOM, FL_CalfCOM, FL_FootCOM, FL_CalfFootCOM;
Eigen::Vector3d RR_HipCOM, RR_ThighCOM, RR_CalfCOM, RR_FootCOM, RR_CalfFootCOM;
Eigen::Vector3d RL_HipCOM, RL_ThighCOM, RL_CalfCOM, RL_FootCOM, RL_CalfFootCOM;
const double HipMass=0.696, ThighMass=1.013, CalfMass=0.166, FootMass=0.06, CalfFootMass=0.166+0.06;
const double TrunkMass=4.713 , ImuMass= 0.001, TrunkImuMass=4.714;


Eigen::Vector3d RootP;
Eigen::Matrix3d FR_CalfRot, FR_ThighRot, FR_HipRot;
Eigen::Matrix3d FL_CalfRot, FL_ThighRot, FL_HipRot;
Eigen::Matrix3d RR_CalfRot, RR_ThighRot, RR_HipRot;
Eigen::Matrix3d RL_CalfRot, RL_ThighRot, RL_HipRot;
Eigen::Matrix3d Roota;


Eigen::Matrix3d FR_Calf_I, FR_Thigh_I, FR_Hip_I, FR_Foot_I,FR_CalfFoot_I;
Eigen::Matrix3d FL_Calf_I, FL_Thigh_I, FL_Hip_I, FL_Foot_I,FL_CalfFoot_I;
Eigen::Matrix3d RR_Calf_I, RR_Thigh_I, RR_Hip_I, RR_Foot_I,RR_CalfFoot_I;
Eigen::Matrix3d RL_Calf_I, RL_Thigh_I, RL_Hip_I, RL_Foot_I,RL_CalfFoot_I;
Eigen::Matrix3d Trunk_I, Imu_I, TrunkImu_I;



Eigen::Matrix3d rot_x(double angle){
  Eigen::Matrix3d X;
  X << 1, 0, 0,
       0, cos(angle),-sin(angle),
       0, sin(angle), cos(angle);
  return X;
}
Eigen::Matrix3d rot_y(double angle){
  Eigen::Matrix3d Y;
  Y << cos(angle), 0, sin(angle),
                   0, 1,0,
      -sin(angle), 0, cos(angle);
  return Y;
}

Eigen::Matrix3d skew(Eigen::Vector3d& V){
  Eigen::Matrix3d skew_V;
  skew_V << 0,        -V(2),         V(1),
      V(2),          0,      -V(0),
      -V(1),       V(0),          0;
  return skew_V;
}
Eigen::Matrix3d Urdf_I(Eigen::VectorXd& i){
  Eigen::Matrix3d I;
  I << i(0), i(1), i(2),
      i(1), i(3), i(4),
      i(2),i(4),i(5);
  return I;
}
void getCurrentState(const Eigen::VectorXd& gc){ ///get robot configuration parameter form URDF & GC
  ///rotation axis of HTC
  p1 << 1, 0, 0;
  p2 << 0, 1, 0;
  p3 << 0, 1, 0;

  ///com position
  TrunkCOM << 0.012731, 0.002186, 0.000515;
  ImuCOM << 0,0,0;
  TrunkImuCOM = (TrunkMass*TrunkCOM+ImuMass*ImuCOM)/TrunkImuMass;

  FR_HipCOM << -0.003311, -0.000635, 3.1e-05;
  FR_ThighCOM << -0.003237, 0.022327, -0.027326;
  FR_CalfCOM << 0.006435, 0.0, -0.107388;
  FR_FootCOM << 0, 0 ,-0.2;
  FR_CalfFootCOM = (CalfMass*FR_CalfCOM+FootMass*FR_FootCOM)/CalfFootMass;

  FL_HipCOM << -0.003311, 0.000635, 3.1e-05;
  FL_ThighCOM << -0.003237 ,-0.022327, -0.027326;
  FL_CalfCOM << 0.006435, 0.0, -0.107388;
  FL_FootCOM << 0, 0 ,-0.2;
  FL_CalfFootCOM = (CalfMass*FL_CalfCOM+FootMass*FL_FootCOM)/CalfFootMass;

  RR_HipCOM << 0.003311, -0.000635, 3.1e-05;
  RR_ThighCOM << -0.003237, 0.022327, -0.027326;
  RR_CalfCOM << 0.006435, 0.0, -0.107388;
  RR_FootCOM << 0, 0 ,-0.2;
  RR_CalfFootCOM = (CalfMass*RR_CalfCOM+FootMass*RR_FootCOM)/CalfFootMass;

  RL_HipCOM << 0.003311, 0.000635, 3.1e-05;
  RL_ThighCOM << -0.003237, -0.022327, -0.027326;
  RL_CalfCOM << 0.006435, 0.0, -0.107388;
  RL_FootCOM << 0, 0 ,-0.2;
  RL_CalfFootCOM = (CalfMass*RL_CalfCOM+FootMass*RL_FootCOM)/CalfFootMass;
  ///joint rotation
  FR_HipRot   = rot_x(gc(7));
  FR_ThighRot = rot_y(gc(8));
  FR_CalfRot  = rot_y(gc(9));
  FL_HipRot   = rot_x(gc(10));
  FL_ThighRot = rot_y(gc(11));
  FL_CalfRot  = rot_y(gc(12));
  RR_HipRot   = rot_x(gc(13));
  RR_ThighRot = rot_y(gc(14));
  RR_CalfRot  = rot_y(gc(15));
  RL_HipRot   = rot_x(gc(16));
  RL_ThighRot = rot_y(gc(17));
  RL_CalfRot  = rot_y(gc(18));
  ///link length
  FR_FBaseToHip << 0.183, -0.047,0;
  FR_HipToThigh << 0, -0.08505, 0;
  FL_FBaseToHip << 0.183, 0.047,0;
  FL_HipToThigh << 0, 0.08505, 0;
  RR_FBaseToHip << -0.183, -0.047,0;
  RR_HipToThigh << 0, -0.08505, 0;
  RL_FBaseToHip << -0.183, 0.047,0;
  RL_HipToThigh << 0, 0.08505, 0;
  ThighToCalf << 0, 0, -0.2;
  CalfToFoot << 0, 0, -0.2;
  RootP << gc(0), gc(1), gc(2);
  Eigen::Quaterniond RootA(gc(3),gc(4),gc(5),gc(6));
  Roota = RootA.normalized().toRotationMatrix();
  ///inertia tensor
  Eigen::VectorXd T(6),I(6);
  Eigen::VectorXd FR_H(6),FR_TH(6),FR_C(6),FR_F(6);
  Eigen::VectorXd FL_H(6),FL_TH(6),FL_C(6),FL_F(6);
  Eigen::VectorXd RR_H(6),RR_TH(6),RR_C(6),RR_F(6);
  Eigen::VectorXd RL_H(6),RL_TH(6),RL_C(6),RL_F(6);
  T<<0.01683993, 8.3902e-05, 0.000597679, 0.056579028, 2.5134e-05, 0.064713601;
  I<< 0.0001, 0, 0, 0.000001, 0, 0.0001;
  FR_H<<0.000469246, 9.409e-06, -3.42e-07, 0.00080749, 4.66e-07, 0.000552929;
  FR_TH<<0.005529065,-4.825e-06, 0.000343869,0.005139339,-2.2448e-05,0.001367788;
  FR_C<<0.002997972, 0.0, -0.000141163, 0.003014022, 0.0, 3.2426e-05;
  FR_F<<9.6e-06, 0.0, 0.0, 9.6e-06, 0.0, 9.6e-06;

  FL_H <<0.000469246,-9.409e-06, -3.42e-07, 0.00080749, -4.66e-07, 0.000552929;
  FL_TH <<0.005529065, 4.825e-06, 0.000343869, 0.005139339, 2.2448e-05, 0.001367788;
  FL_C << 0.002997972, 0.0, -0.000141163, 0.003014022, 0.0, 3.2426e-05;
  FL_F << 9.6e-06, 0.0, 0.0, 9.6e-06, 0.0, 9.6e-06;

  RR_H << 0.000469246, -9.409e-06, 3.42e-07, 0.00080749, 4.66e-07, 0.000552929;
  RR_TH << 0.005529065, -4.825e-06, 0.000343869, 0.005139339, -2.2448e-05, 0.001367788;
  RR_C << 0.002997972, 0.0, -0.000141163, 0.003014022, 0.0, 3.2426e-05;
  RR_F << 9.6e-06, 0.0, 0.0, 9.6e-06, 0.0, 9.6e-06;

  RL_H << 0.000469246, 9.409e-06, 3.42e-07, 0.00080749, -4.66e-07, 0.000552929;
  RL_TH << 0.005529065, 4.825e-06, 0.000343869, 0.005139339, 2.2448e-05, 0.001367788;
  RL_C << 0.002997972, 0.0, -0.000141163, 0.003014022, 0.0, 3.2426e-05;
  RL_F << 9.6e-06, 0.0, 0.0, 9.6e-06, 0.0, 9.6e-06;

  Trunk_I = Urdf_I(T);
  Imu_I = Urdf_I(I);
  Eigen::Vector3d r1=TrunkCOM-TrunkImuCOM;
  Eigen::Vector3d r2=ImuCOM-TrunkImuCOM;
  TrunkImu_I= Trunk_I+Imu_I-TrunkMass*skew(r1)*skew(r1)-ImuMass*skew(r2)*skew(r2);

  FR_Hip_I = Urdf_I(FR_H);
  FR_Thigh_I = Urdf_I(FR_TH);
  FR_Calf_I = Urdf_I(FR_C);
  FR_Foot_I = Urdf_I(FR_F);
  Eigen::Vector3d r3=FR_CalfCOM-FR_CalfFootCOM;
  Eigen::Vector3d r4=FR_FootCOM-FR_CalfFootCOM;
  FR_CalfFoot_I=FR_Calf_I+FR_Foot_I-CalfMass*skew(r3)*skew(r3)-FootMass*skew(r4)*skew(r4);

  FL_Hip_I = Urdf_I(FL_H);
  FL_Thigh_I = Urdf_I(FL_TH);
  FL_Calf_I = Urdf_I(FL_C);
  FL_Foot_I = Urdf_I(FL_F);
  Eigen::Vector3d r5=FL_CalfCOM-FL_CalfFootCOM;
  Eigen::Vector3d r6=FL_FootCOM-FL_CalfFootCOM;
  FL_CalfFoot_I=FL_Calf_I+FL_Foot_I-CalfMass*skew(r5)*skew(r5)-FootMass*skew(r6)*skew(r6);

  RR_Hip_I = Urdf_I(RR_H);
  RR_Thigh_I = Urdf_I(RR_TH);
  RR_Calf_I = Urdf_I(RR_C);
  RR_Foot_I = Urdf_I(RR_F);
  Eigen::Vector3d r7=RR_CalfCOM-RR_CalfFootCOM;
  Eigen::Vector3d r8=RR_FootCOM-RR_CalfFootCOM;
  RR_CalfFoot_I=RR_Calf_I+RR_Foot_I-CalfMass*skew(r7)*skew(r7)-FootMass*skew(r8)*skew(r8);


  RL_Hip_I = Urdf_I(RL_H);
  RL_Thigh_I = Urdf_I(RL_TH);
  RL_Calf_I = Urdf_I(RL_C);
  RL_Foot_I = Urdf_I(RL_F);
  Eigen::Vector3d r9=RL_CalfCOM-RL_CalfFootCOM;
  Eigen::Vector3d r10=RL_FootCOM-RL_CalfFootCOM;
  RL_CalfFoot_I=RL_Calf_I+RL_Foot_I-CalfMass*skew(r9)*skew(r9)-FootMass*skew(r10)*skew(r10);

}



double Mass(int index){
  if(index==0){
    return TrunkImuMass;
  }else if (index==1 || index==4 || index==7 || index==10){
    return HipMass;
  }else if (index==2 || index==5 || index==8 || index==11){
    return ThighMass;
  }else if (index==3 || index==6 || index==9 || index==12){
    return CalfFootMass;
  }else {
    std::cout<<"invalid mass index";
    return 0;
  }
}

///debugging-> OK
Eigen::Vector3d getJointAbsPos(const Eigen::VectorXd& gc, int index1, int index2){
//index1 : FR #1 FL #2 RR #3 RL #4
//index2 : Hip #1 Thigh #2 Calf #3 Foot #4

  Eigen::Vector3d FBaseToHip;
  Eigen::Vector3d HipToThigh;
  Eigen::Matrix3d HipRot;
  Eigen::Matrix3d ThighRot;
  Eigen::Matrix3d CalfRot;

  if(index1==0){
    return RootP;
  }else if (index1==1){
    FBaseToHip=FR_FBaseToHip;
    HipToThigh=FR_HipToThigh;
    HipRot=FR_HipRot;
    ThighRot=FR_ThighRot;
    CalfRot=FR_CalfRot;
  }else if (index1==2){
    FBaseToHip=FL_FBaseToHip;
    HipToThigh=FL_HipToThigh;
    HipRot=FL_HipRot;
    ThighRot=FL_ThighRot;
    CalfRot=FL_CalfRot;
  }else if (index1==3){
    FBaseToHip=RR_FBaseToHip;
    HipToThigh=RR_HipToThigh;
    HipRot=RR_HipRot;
    ThighRot=RR_ThighRot;
    CalfRot=RR_CalfRot;
  }else if (index1==4){
    FBaseToHip=RL_FBaseToHip;
    HipToThigh=RL_HipToThigh;
    HipRot=RL_HipRot;
    ThighRot=RL_ThighRot;
    CalfRot=RL_CalfRot;
  }else{
    return Eigen::Vector3d::Zero();
  }

  if(index2==1){ ///hip
    return RootP+Roota*FBaseToHip;
  }else if(index2==2){ ///thigh
    return RootP+Roota*(FBaseToHip+HipRot*HipToThigh);
  }else if(index2==3){ ///calf
    return RootP+Roota*(FBaseToHip+HipRot*(HipToThigh+ThighRot*ThighToCalf));
  }else if(index2==4){///foot
    return RootP+Roota*(FBaseToHip+HipRot*(HipToThigh+ThighRot*(ThighToCalf+CalfRot*CalfToFoot)));
  }else{
    return Eigen::Vector3d::Zero();
  }

}

Eigen::Vector3d getComAbsPos(const Eigen::VectorXd& gc, int index1, int index2){
//index1 : FR #1 FL #2 RR #3 RL #4 Trunk #0
//index2 : Hip #1 Thigh #2 CalfFoot #3

  Eigen::Vector3d HipCOM;
  Eigen::Vector3d ThighCOM;
  Eigen::Vector3d CalfFootCOM;

  Eigen::Matrix3d HipRot;
  Eigen::Matrix3d ThighRot;
  Eigen::Matrix3d CalfRot;

  if(index1==0){
    return RootP+Roota*TrunkImuCOM;
  }else if (index1==1){
    HipRot=FR_HipRot;
    ThighRot=FR_ThighRot;
    CalfRot=FR_CalfRot;
    HipCOM=FR_HipCOM;
    ThighCOM=FR_ThighCOM;
    CalfFootCOM=FR_CalfFootCOM;
  }else if (index1==2){
    HipRot=FL_HipRot;
    ThighRot=FL_ThighRot;
    CalfRot=FL_CalfRot;
    HipCOM=FL_HipCOM;
    ThighCOM=FL_ThighCOM;
    CalfFootCOM=FL_CalfFootCOM;
  }else if (index1==3){
    HipRot=RR_HipRot;
    ThighRot=RR_ThighRot;
    CalfRot=RR_CalfRot;
    HipCOM=RR_HipCOM;
    ThighCOM=RR_ThighCOM;
    CalfFootCOM=RR_CalfFootCOM;
  }else if (index1==4){
    HipRot=RL_HipRot;
    ThighRot=RL_ThighRot;
    CalfRot=RL_CalfRot;
    HipCOM=RL_HipCOM;
    ThighCOM=RL_ThighCOM;
    CalfFootCOM=RL_CalfFootCOM;
  }else{
    return Eigen::Vector3d::Zero();
  }

  if(index2==1){ ///hip
    return getJointAbsPos(gc,index1,index2)+Roota*HipRot*HipCOM;
  }else if(index2==2){ ///thigh
    return getJointAbsPos(gc,index1,index2)+Roota*HipRot*ThighRot*ThighCOM;
  }else if(index2==3){ ///calf
    return getJointAbsPos(gc,index1,index2)+Roota*HipRot*ThighRot*CalfRot*CalfFootCOM;
  }else{
    return Eigen::Vector3d::Zero();
  }

}

Eigen::Vector3d getRotAxis_w (int index1, int index2){
if (index1==1){

  if(index2==1) return Roota*p1;
  else if(index2==2) return Roota*FR_HipRot*p2;
  else if(index2==3) return Roota*FR_HipRot*FR_ThighRot*p3;
  else {std::cout<<"invalid rot index";
    return Eigen::Vector3d::Zero();}

}else if (index1==2){

  if(index2==1) return Roota*p1;
  else if(index2==2) return Roota*FL_HipRot*p2;
  else if(index2==3) return Roota*FL_HipRot*FL_ThighRot*p3;
  else {std::cout<<"invalid rot index";
    return Eigen::Vector3d::Zero();}

}else if (index1==3){

  if(index2==1) return Roota*p1;
  else if(index2==2) return Roota*RR_HipRot*p2;
  else if(index2==3) return Roota*RR_HipRot*RR_ThighRot*p3;
  else {std::cout<<"invalid rot index";
    return Eigen::Vector3d::Zero();}

}else if (index1==4){

  if(index2==1) return Roota*p1;
  else if(index2==2) return Roota*RL_HipRot*p2;
  else if(index2==3) return Roota*RL_HipRot*RL_ThighRot*p3;
  else {std::cout<<"invalid rot index";
    return Eigen::Vector3d::Zero();}

}else{ std::cout<<"invalid rot index";
  return Eigen::Vector3d::Zero();}
}

Eigen::Matrix3d getInertiaTensor_w (const Eigen::VectorXd& gc, int index){
  if (index==0){
    return Roota*TrunkImu_I*Roota.transpose();
  } else if (index==1){
    return (Roota*FR_HipRot)*FR_Hip_I*(Roota*FR_HipRot).transpose();
  } else if (index==2){
    return (Roota*FR_HipRot*FR_ThighRot)*FR_Thigh_I*(Roota*FR_HipRot*FR_ThighRot).transpose();
  } else if (index==3){
    return (Roota*FR_HipRot*FR_ThighRot*FR_CalfRot)*FR_CalfFoot_I*(Roota*FR_HipRot*FR_ThighRot*FR_CalfRot).transpose();
  }///FR
  else if (index==4){
    return (Roota*FL_HipRot)*FL_Hip_I*(Roota*FL_HipRot).transpose();
  }else if (index==5){
    return (Roota*FL_HipRot*FL_ThighRot)*FL_Thigh_I*(Roota*FL_HipRot*FL_ThighRot).transpose();
  }else if (index==6){
    return (Roota*FL_HipRot*FL_ThighRot*FL_CalfRot)*FL_CalfFoot_I*(Roota*FL_HipRot*FL_ThighRot*FL_CalfRot).transpose();
  }///FL
  else if (index==7){
    return (Roota*RR_HipRot)*RR_Hip_I*(Roota*RR_HipRot).transpose();
  }else if (index==8){
    return (Roota*RR_HipRot*RR_ThighRot)*RR_Thigh_I*(Roota*RR_HipRot*RR_ThighRot).transpose();
  }else if (index==9){
    return (Roota*RR_HipRot*RR_ThighRot*RR_CalfRot)*RR_CalfFoot_I*(Roota*RR_HipRot*RR_ThighRot*RR_CalfRot).transpose();
  }///RR
  else if (index==10){
    return (Roota*RL_HipRot)*RL_Hip_I*(Roota*RL_HipRot).transpose();
  }else if (index==11){
    return (Roota*RL_HipRot*RL_ThighRot)*RL_Thigh_I*(Roota*RL_HipRot*RL_ThighRot).transpose();
  }else if (index==12){
    return (Roota*RL_HipRot*RL_ThighRot*RL_CalfRot)*RL_CalfFoot_I*(Roota*RL_HipRot*RL_ThighRot*RL_CalfRot).transpose();
  }///RL
  else{
    std::cout<<"invalid joint index";
    return Eigen::Matrix3d::Zero();
  }
}

Eigen::MatrixXd getJointPosJ (const Eigen::VectorXd& gc, int index) {
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d r1, r2, r3;
  switch (index) {
    case 0: ///Base
      J << I3,Z3,Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      r1 = getJointAbsPos(gc, 1, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skew(r1),Z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      r1 = getJointAbsPos(gc, 1, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 1, 2) - getJointAbsPos(gc, 1, 1);
      J << I3, -1*skew(r1), -1*skew(r2)*getRotAxis_w(1,1),z3,z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C 1-3
      r1 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 1, 1);
      r3 = getJointAbsPos(gc, 1, 3) - getJointAbsPos(gc, 1, 2);
      J << I3, -1*skew(r1), -1*skew(r2)*getRotAxis_w(1,1), -1*skew(r3)* getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      r1 = getJointAbsPos(gc, 2, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skew(r1),Z3,Z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      r1 = getJointAbsPos(gc, 2, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 2, 2) - getJointAbsPos(gc, 2, 1);
      J << I3, -1*skew(r1),Z3, -1*skew(r2)*getRotAxis_w(2,1),z3,z3,Z3,Z3;
      break;
    case 6: ///FL_C 2-3
      r1 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 2, 1);
      r3 = getJointAbsPos(gc, 2, 3) - getJointAbsPos(gc, 2, 2);
      J << I3, -1*skew(r1),Z3, -1*skew(r2)*getRotAxis_w(2,1), -1*skew(r3)* getRotAxis_w(2,2),z3,Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      r1 = getJointAbsPos(gc, 3, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skew(r1),Z3,Z3,Z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      r1 = getJointAbsPos(gc, 3, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 3, 2) - getJointAbsPos(gc, 3, 1);
      J << I3, -1*skew(r1),Z3,Z3, -1*skew(r2)*getRotAxis_w(3,1),z3,z3,Z3;
      break;
    case 9: ///RR_C 3-3
      r1 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 3, 1);
      r3 = getJointAbsPos(gc, 3, 3) - getJointAbsPos(gc, 3, 2);
      J << I3, -1*skew(r1),Z3,Z3, -1*skew(r2)*getRotAxis_w(3,1), -1*skew(r3)* getRotAxis_w(3,2),z3,Z3;
      break;
    case 10: ///RL_H 4-1
      r1 = getJointAbsPos(gc, 4, 1) - getJointAbsPos(gc, 0, 0);
      J << I3, -1 * skew(r1),Z3,Z3,Z3,Z3;
      break;
    case 11: ///RL_TH 4-2
      r1 = getJointAbsPos(gc, 4, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 4, 2) - getJointAbsPos(gc, 4, 1);
      J << I3, -1*skew(r1),Z3,Z3,Z3, -1*skew(r2)*getRotAxis_w(4,1),z3,z3;
      break;
    case 12: ///RL_C 4-3
      r1 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 4, 1);
      r3 = getJointAbsPos(gc, 4, 3) - getJointAbsPos(gc, 4, 2);
      J << I3, -1*skew(r1),Z3,Z3,Z3, -1*skew(r2)*getRotAxis_w(4,1), -1*skew(r3)* getRotAxis_w(4,2),z3;
      break;
    default:
      std::cout<<"invalid joint index";
      break;
  }
  return J;
}

Eigen::MatrixXd getLinkCOMPosJ (const Eigen::VectorXd& gc, int index){
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d r1, r2, r3, r4;
  switch (index) {
    case 0: ///TrunkImu
      r1 = getComAbsPos(gc, 0, 0) - getJointAbsPos(gc, 0, 0);
      J << I3,
      -1 * skew(r1),Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      r1 = getComAbsPos(gc, 1, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,1,1) - getJointAbsPos(gc,1,1);
      J << I3, -1 * skew(r1),-1*skew(r2)*getRotAxis_w(1,1), z3,z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      r1 = getComAbsPos(gc, 1, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 1);
      r3 = getComAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 2);
      J << I3, -1*skew(r1), -1*skew(r2)*getRotAxis_w(1,1),-1*skew(r3)* getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C+F 1-3
      r1 = getComAbsPos(gc, 1, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 1);
      r3 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 2);
      r4 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 3);
      J << I3, -1*skew(r1), -1*skew(r2)*getRotAxis_w(1,1), -1*skew(r3)* getRotAxis_w(1,2),-1*skew(r4)* getRotAxis_w(1,3),Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      r1 = getComAbsPos(gc, 2, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,2,1) - getJointAbsPos(gc,2,1);
      J << I3, -1 * skew(r1), Z3, -1*skew(r2)*getRotAxis_w(2,1), z3,z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      r1 = getComAbsPos(gc, 2, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 1);
      r3 = getComAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 2);
      J << I3, -1*skew(r1), Z3, -1*skew(r2)*getRotAxis_w(2,1),-1*skew(r3)* getRotAxis_w(2,2),z3,Z3,Z3;
      break;
    case 6: ///FL_C+F 2-3
      r1 = getComAbsPos(gc, 2, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 1);
      r3 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 2);
      r4 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 3);
      J << I3, -1*skew(r1), Z3, -1*skew(r2)*getRotAxis_w(2,1), -1*skew(r3)* getRotAxis_w(2,2),-1*skew(r4)* getRotAxis_w(2,3),Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      r1 = getComAbsPos(gc, 3, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,3,1) - getJointAbsPos(gc,3,1);
      J << I3, -1 * skew(r1), Z3, Z3,-1*skew(r2)*getRotAxis_w(3,1), z3,z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      r1 = getComAbsPos(gc, 3, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 1);
      r3 = getComAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 2);
      J << I3, -1*skew(r1), Z3, Z3, -1*skew(r2)*getRotAxis_w(3,1),-1*skew(r3)* getRotAxis_w(3,2),z3,Z3;
      break;
    case 9: ///RR_C+F 3-3
      r1 = getComAbsPos(gc, 3, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 1);
      r3 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 2);
      r4 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 3);
      J << I3, -1*skew(r1), Z3, Z3, -1*skew(r2)*getRotAxis_w(3,1), -1*skew(r3)* getRotAxis_w(3,2),-1*skew(r4)* getRotAxis_w(3,3),Z3;
      break;
    case 10: ///RL_H 4-1
      r1 = getComAbsPos(gc, 4, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,4,1) - getJointAbsPos(gc,4,1);
      J << I3, -1 * skew(r1), Z3, Z3,Z3, -1*skew(r2)*getRotAxis_w(4,1), z3,z3;
      break;
    case 11: ///RL_TH 4-2
      r1 = getComAbsPos(gc, 4, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 1);
      r3 = getComAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 2);
      J << I3, -1*skew(r1), Z3, Z3,Z3, -1*skew(r2)*getRotAxis_w(4,1),-1*skew(r3)* getRotAxis_w(4,2),z3;
      break;
    case 12: ///RL_C+F 4-3
      r1 = getComAbsPos(gc, 4, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 1);
      r3 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 2);
      r4 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 3);
      J << I3, -1*skew(r1), Z3, Z3,Z3, -1*skew(r2)*getRotAxis_w(4,1), -1*skew(r3)* getRotAxis_w(4,2),-1*skew(r4)* getRotAxis_w(4,3);
      break;
    default:
      std::cout<<"invalid body index";
      break;
  }
  return J;
}

Eigen::MatrixXd getJointAngJ (const Eigen::VectorXd& gc, int index){
  Eigen::MatrixXd J(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  switch (index) {
    case 0: ///Base
      J << Z3,I3,Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      J << Z3,I3,getRotAxis_w(1,1),z3,z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      J << Z3,I3, getRotAxis_w(1,1),getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C 1-3
      J << Z3,I3, getRotAxis_w(1,1),getRotAxis_w(1,2),getRotAxis_w(1,3),Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      J << Z3,I3,Z3,getRotAxis_w(2,1),z3,z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      J << Z3,I3,Z3, getRotAxis_w(2,1),getRotAxis_w(2,2),z3,Z3,Z3;
      break;
    case 6: ///FL_C 2-3
      J << Z3,I3,Z3, getRotAxis_w(2,1),getRotAxis_w(2,2),getRotAxis_w(2,3),Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      J << Z3,I3,Z3,Z3,getRotAxis_w(3,1),z3,z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      J << Z3,I3,Z3,Z3, getRotAxis_w(3,1),getRotAxis_w(3,2),z3,Z3;
      break;
    case 9: ///RR_C 3-3
      J << Z3,I3,Z3,Z3, getRotAxis_w(3,1),getRotAxis_w(3,2),getRotAxis_w(3,3),Z3;
      break;
    case 10: ///RL_H 4-1
      J << Z3,I3,Z3,Z3,Z3,getRotAxis_w(4,1),z3,z3;
      break;
    case 11: ///RL_TH 4-2
      J << Z3,I3,Z3,Z3,Z3, getRotAxis_w(4,1),getRotAxis_w(4,2),z3;
      break;
    case 12: ///RL_C 4-3
      J << Z3,I3,Z3,Z3,Z3, getRotAxis_w(4,1),getRotAxis_w(4,2),getRotAxis_w(4,3);
      break;
    default:
      std::cout<<"invalid joint index";
      break;
  }
  return J;
}

Eigen::Vector3d getFrameVel ( const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index) {
  Eigen::Vector3d V;

  V = getJointPosJ(gc,index)*gv;
  return V;
}

Eigen::Vector3d getComVel(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::Vector3d VC;

  VC = getLinkCOMPosJ(gc,index)*gv;
  return VC;
}

Eigen::Vector3d getFrameAngVel(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::Vector3d W;

  W = getJointAngJ(gc,index)*gv;
  return W;
}


Eigen::MatrixXd getLinkCOMPosJDot (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){
  Eigen::MatrixXd JDot(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d r1, r2, r3, r4;
  Eigen::Vector3d v1, v2, v3, v4;
  Eigen::Vector3d w1, w2, w3;
  switch (index) {
    case 0: ///TrunkImu
      r1 = getComAbsPos(gc, 0, 0) - getJointAbsPos(gc, 0, 0);
      v1 = getComVel(gc,gv,0)-getFrameVel(gc,gv,0);
      JDot << Z3,
      -1 * skew(v1),Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      r1 = getComAbsPos(gc, 1, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,1,1) - getJointAbsPos(gc,1,1);
      v1 =  getLinkCOMPosJ(gc,1)*gv-getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,1)-getFrameVel(gc,gv,1);
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,
      -1 * skew(v1),
      -1*skew(v2)*getRotAxis_w(1,1) -1*skew(r2)*(skew(w1)*getRotAxis_w(1,1)),
      z3,z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      r1 = getComAbsPos(gc, 1, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 1);
      r3 = getComAbsPos(gc, 1, 2)- getJointAbsPos(gc, 1, 2);
      v1 = getComVel(gc,gv,2)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,2)- getFrameVel(gc,gv,1);
      v3 = getComVel(gc,gv,2)- getFrameVel(gc,gv,2);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,1);
      JDot << Z3,
      -1 * skew(v1),
      -1*skew(v2)* getRotAxis_w(1,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(1,1)),
      -1*skew(v3)* getRotAxis_w(1,2)-1*skew(r3)*(skew(w2)*getRotAxis_w(1,2)),
      z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C+F 1-3
      r1 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 1);
      r3 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 2);
      r4 = getComAbsPos(gc, 1, 3)- getJointAbsPos(gc, 1, 3);
      v1 = getComVel(gc,gv,3)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,3)- getFrameVel(gc,gv,1);
      v3 = getComVel(gc,gv,3)- getFrameVel(gc,gv,2);
      v4 = getComVel(gc,gv,3)- getFrameVel(gc,gv,3);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,1);
      w3 = getFrameAngVel(gc,gv,2);
      JDot << Z3,
      -1 * skew(v1),
      -1*skew(v2)* getRotAxis_w(1,1)-1*skew(r2)*(skew(w1)* getRotAxis_w(1,1)),
      -1*skew(v3)* getRotAxis_w(1,2)-1*skew(r3)*(skew(w2)* getRotAxis_w(1,2)),
      -1*skew(v4)* getRotAxis_w(1,3)-1*skew(r4)*(skew(w3)* getRotAxis_w(1,3)),
      Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      r1 = getComAbsPos(gc,2,1) - getJointAbsPos(gc,0, 0);
      r2 = getComAbsPos(gc,2,1) - getJointAbsPos(gc,2, 1);
      v1 = getComVel(gc,gv,4)-getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,4)-getFrameVel(gc,gv,4);
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,
      -1 * skew(v1),Z3,
      -1*skew(v2)*getRotAxis_w(2,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(2,1)),
      z3,z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      r1 = getComAbsPos(gc, 2, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 1);
      r3 = getComAbsPos(gc, 2, 2)- getJointAbsPos(gc, 2, 2);
      v1 = getComVel(gc,gv,5)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,5)- getFrameVel(gc,gv,4);
      v3 = getComVel(gc,gv,5)- getFrameVel(gc,gv,5);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,4);
      JDot << Z3,
      -1 * skew(v1),Z3,
      -1*skew(v2)* getRotAxis_w(2,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(2,1)),
      -1*skew(v3)* getRotAxis_w(2,2)-1*skew(r3)*(skew(w2)*getRotAxis_w(2,2)),
      z3,Z3,Z3;
      break;
    case 6: ///FL_C+F 2-3
      r1 = getComAbsPos(gc, 2, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 1);
      r3 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 2);
      r4 = getComAbsPos(gc, 2, 3)- getJointAbsPos(gc, 2, 3);
      v1 = getComVel(gc,gv,6)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,6)- getFrameVel(gc,gv,4);
      v3 = getComVel(gc,gv,6)- getFrameVel(gc,gv,5);
      v4 = getComVel(gc,gv,6)- getFrameVel(gc,gv,6);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,4);
      w3 = getFrameAngVel(gc,gv,5);
      JDot << Z3,
      -1 * skew(v1),Z3,
      -1*skew(v2)* getRotAxis_w(2,1)-1*skew(r2)*(skew(w1)* getRotAxis_w(2,1)),
      -1*skew(v3)* getRotAxis_w(2,2)-1*skew(r3)*(skew(w2)* getRotAxis_w(2,2)),
      -1*skew(v4)* getRotAxis_w(2,3)-1*skew(r4)*(skew(w3)* getRotAxis_w(2,3)),
      Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      r1 = getComAbsPos(gc, 3, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,3,1) - getJointAbsPos(gc,3,1);
      v1 = getComVel(gc,gv,7)-getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,7)-getFrameVel(gc,gv,7);
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,
          -1*skew(v2)*getRotAxis_w(3,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(3,1)),
          z3,z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      r1 = getComAbsPos(gc, 3, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 1);
      r3 = getComAbsPos(gc, 3, 2)- getJointAbsPos(gc, 3, 2);
      v1 = getComVel(gc,gv,8)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,8)- getFrameVel(gc,gv,7);
      v3 = getComVel(gc,gv,8)- getFrameVel(gc,gv,8);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,7);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,
          -1*skew(v2)* getRotAxis_w(3,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(3,1)),
          -1*skew(v3)* getRotAxis_w(3,2)-1*skew(r3)*(skew(w2)*getRotAxis_w(3,2)),
          z3,Z3;
      break;
    case 9: ///RR_C+F 3-3
      r1 = getComAbsPos(gc, 3, 3) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 1);
      r3 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 2);
      r4 = getComAbsPos(gc, 3, 3)- getJointAbsPos(gc, 3, 3);
      v1 = getComVel(gc,gv,9)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,9)- getFrameVel(gc,gv,7);
      v3 = getComVel(gc,gv,9)- getFrameVel(gc,gv,8);
      v4 = getComVel(gc,gv,9)- getFrameVel(gc,gv,9);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,7);
      w3 = getFrameAngVel(gc,gv,8);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,
          -1*skew(v2)* getRotAxis_w(3,1)-1*skew(r2)*(skew(w1)* getRotAxis_w(3,1)),
          -1*skew(v3)* getRotAxis_w(3,2)-1*skew(r3)*(skew(w2)* getRotAxis_w(3,2)),
          -1*skew(v4)* getRotAxis_w(3,3)-1*skew(r4)*(skew(w3)* getRotAxis_w(3,3)),
          Z3;
      break;
    case 10: ///RL_H 4-1
      r1 = getComAbsPos(gc, 4, 1) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc,4,1) - getJointAbsPos(gc,4,1);
      v1 = getComVel(gc,gv,10)-getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,10)-getFrameVel(gc,gv,10);
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,Z3,
          -1* skew(v2)*getRotAxis_w(4,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(4,1)),
          z3,z3;
      break;
    case 11: ///RL_TH 4-2
      r1 = getComAbsPos(gc, 4, 2) - getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 1);
      r3 = getComAbsPos(gc, 4, 2)- getJointAbsPos(gc, 4, 2);
      v1 = getComVel(gc,gv,11)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,11)- getFrameVel(gc,gv,10);
      v3 = getComVel(gc,gv,11)- getFrameVel(gc,gv,11);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,10);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,Z3,
          -1*skew(v2)* getRotAxis_w(4,1)-1*skew(r2)*(skew(w1)*getRotAxis_w(4,1)),
          -1*skew(v3)* getRotAxis_w(4,2)-1*skew(r3)*(skew(w2)*getRotAxis_w(4,2)),
          z3;
      break;
    case 12: ///RL_C+F 4-3
      r1 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 0, 0);
      r2 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 1);
      r3 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 2);
      r4 = getComAbsPos(gc, 4, 3)- getJointAbsPos(gc, 4, 3);
      v1 = getComVel(gc,gv,12)- getFrameVel(gc,gv,0);
      v2 = getComVel(gc,gv,12)- getFrameVel(gc,gv,10);
      v3 = getComVel(gc,gv,12)- getFrameVel(gc,gv,11);
      v4 = getComVel(gc,gv,12)- getFrameVel(gc,gv,12);
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,10);
      w3 = getFrameAngVel(gc,gv,11);
      JDot << Z3,
          -1 * skew(v1),Z3,Z3,Z3,
          -1*skew(v2)* getRotAxis_w(4,1)-1*skew(r2)*(skew(w1)* getRotAxis_w(4,1)),
          -1*skew(v3)* getRotAxis_w(4,2)-1*skew(r3)*(skew(w2)* getRotAxis_w(4,2)),
          -1*skew(v4)* getRotAxis_w(4,3)-1*skew(r4)*(skew(w3)* getRotAxis_w(4,3));
      break;
    default:
      std::cout<<"invalid body index";
      break;
  }
  return JDot;
}


Eigen::MatrixXd getLinkCOMAngJDot (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv, int index){

  Eigen::MatrixXd JDot(3, 18);
  Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d Z3 = Eigen::Matrix3d::Zero(); ///zero matrix for pedding
  Eigen::Vector3d z3 = Eigen::Vector3d::Zero();
  Eigen::Vector3d w1,w2,w3;
  switch (index) {
    case 0: ///Base
      JDot<< Z3,Z3,Z3,Z3,Z3,Z3;
      break;
    case 1: ///FR_H 1-1
      w1 = getFrameAngVel(gc,gv,0);
      JDot<< Z3,Z3,skew(w1)*getRotAxis_w(1,1),z3,z3,Z3,Z3,Z3;
      break;
    case 2: ///FR_TH 1-2
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,1);
      JDot << Z3,Z3, skew(w1)*getRotAxis_w(1,1),skew(w2)*getRotAxis_w(1,2),z3,Z3,Z3,Z3;
      break;
    case 3: ///FR_C 1-3
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,1);
      w3 = getFrameAngVel(gc,gv,2);
      JDot << Z3,Z3, skew(w1)*getRotAxis_w(1,1),skew(w2)*getRotAxis_w(1,2),skew(w3)*getRotAxis_w(1,3),Z3,Z3,Z3;
      break;
    case 4: ///FL_H 2-1
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,Z3,Z3,skew(w1)*getRotAxis_w(2,1),z3,z3,Z3,Z3;
      break;
    case 5: ///FL_TH 2-2
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,4);
      JDot << Z3,Z3,Z3,skew(w1)*getRotAxis_w(2,1),skew(w2)*getRotAxis_w(2,2),z3,Z3,Z3;
      break;
    case 6: ///FL_C 2-3
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,4);
      w3 = getFrameAngVel(gc,gv,5);
      JDot << Z3,Z3,Z3, skew(w1)*getRotAxis_w(2,1),skew(w2)*getRotAxis_w(2,2),skew(w3)*getRotAxis_w(2,3),Z3,Z3;
      break;
    case 7: ///RR_H 3-1
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,Z3,Z3,Z3,skew(w1)*getRotAxis_w(3,1),z3,z3,Z3;
      break;
    case 8: ///RR_TH 3-2
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,7);
      JDot << Z3,Z3,Z3,Z3,skew(w1)*getRotAxis_w(3,1),skew(w2)*getRotAxis_w(3,2),z3,Z3;
      break;
    case 9: ///RR_C 3-3
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,7);
      w3 = getFrameAngVel(gc,gv,8);
      JDot << Z3,Z3,Z3,Z3, skew(w1)*getRotAxis_w(3,1),skew(w2)*getRotAxis_w(3,2),skew(w3)*getRotAxis_w(3,3),Z3;
      break;
    case 10: ///RL_H 4-1
      w1 = getFrameAngVel(gc,gv,0);
      JDot << Z3,Z3,Z3,Z3,Z3,skew(w1)*getRotAxis_w(4,1),z3,z3;
      break;
    case 11: ///RL_TH 4-2
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,10);
      JDot << Z3,Z3,Z3,Z3,Z3, skew(w1)*getRotAxis_w(4,1),skew(w2)*getRotAxis_w(4,2),z3;
      break;
    case 12: ///RL_C 4-3
      w1 = getFrameAngVel(gc,gv,0);
      w2 = getFrameAngVel(gc,gv,10);
      w3 = getFrameAngVel(gc,gv,11);
      JDot << Z3,Z3,Z3,Z3,Z3,skew(w1)*getRotAxis_w(4,1),skew(w2)*getRotAxis_w(4,2),skew(w3)*getRotAxis_w(4,3);
      break;
    default:
      std::cout<<"invalid joint index";
      break;
  }
  return JDot;
}

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  getCurrentState(gc);
  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::MatrixXd M(18, 18);
  M.setZero();
  for (int i=0; i < 13; i++){
    M=M+((getLinkCOMPosJ(gc,i).transpose()*Mass(i)*getLinkCOMPosJ(gc,i)) + (getJointAngJ(gc,i).transpose()*getInertiaTensor_w(gc,i)*getJointAngJ(gc,i)));
  }
  return M;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearities (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  getCurrentState(gc);
  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  Eigen::VectorXd C(18);
  C.setZero(18);
  Eigen::Vector3d G{0,0,-9.81};
  Eigen::Vector3d w;
  for (int i=0; i < 13; i++){
    w= getJointAngJ(gc,i)*gv;
//    std::cout<<"i: "<<i<<std::endl;
    C+= ((getLinkCOMPosJ(gc,i).transpose()*Mass(i)*getLinkCOMPosJDot(gc,gv,i))*gv).head(18);
    C+= ((getJointAngJ(gc,i).transpose()*getInertiaTensor_w(gc,i)*getLinkCOMAngJDot(gc,gv,i))*gv).head(18);
    C+= (getJointAngJ(gc,i).transpose()*(w.cross((getInertiaTensor_w(gc,i)*w)))).head(18);
    C+= -1*getLinkCOMPosJ(gc,i).transpose()*Mass(i)*G;
//      std::cout << "plus term: " << ((getLinkCOMPosJ(gc, i).transpose() * Mass(i) * getLinkCOMPosJDot(gc, gv, i)) * gv).head(18)
//          + ((getJointAngJ(gc, i).transpose() * getInertiaTensor_w(gc, i) * getLinkCOMAngJDot(gc, gv, i)) * gv).head(18)
//          + (getJointAngJ(gc, i).transpose() * (w.cross((getInertiaTensor_w(gc, i) * w)))).head(18) << std::endl;
    //    C+= -1*getLinkCOMPosJ(gc,i).transpose()*G;
  }
  return C;
}

