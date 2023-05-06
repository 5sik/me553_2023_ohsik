#pragma once

#include<Eigen/Core>
#include<Eigen/Dense>
#include <iostream>

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

struct Joint{
  Eigen::Matrix3d jointRot;
  Eigen::Vector3d jointPos_B;
  Eigen::Vector3d jointPos_W;
  Eigen::Matrix3d jointRot_W;
  Eigen::Vector3d jointAxis_B;
  Eigen::Vector3d jointAxis_W;
  Eigen::VectorXd S;
};

struct Body{
  Eigen::Vector3d comPos_B;
  Eigen::Vector3d comPos_W;
  Eigen::Matrix3d I_B;
  Eigen::Matrix3d I_W;
  double mass;

};

struct Comp{
  Eigen::Vector3d c_com_W;
  double c_mass;
  Eigen::Matrix3d c_I_W;
  Eigen::MatrixXd M_spatial;
};

std::vector<Joint> joint;
std::vector<Body> body;
std::vector<Comp> comp;

void initialize_joint(const Eigen::VectorXd& gc){
  joint.resize(4);
  joint[0].jointPos_W << gc(0), gc(1), gc(2);
  Eigen::Quaterniond quat(gc(3),gc(4),gc(5),gc(6));
  joint[0].jointRot_W=quat.normalized().toRotationMatrix();
  joint[0].jointRot=quat.normalized().toRotationMatrix();
  joint[1].jointRot= rot_x(gc(7));
  joint[2].jointRot= rot_y(gc(8));
  joint[3].jointRot= rot_y(gc(9));
  joint[1].jointPos_B << 0.183, -0.047,0;
  joint[2].jointPos_B << 0, -0.08505, 0;
  joint[3].jointPos_B << 0, 0, -0.2;
  joint[1].jointAxis_B << 1, 0, 0;
  joint[2].jointAxis_B << 0, 1, 0;
  joint[3].jointAxis_B << 0, 1, 0;
}

void forward_joint(const Eigen::VectorXd& gc){
  ///joint orientation
  for(int i=1;i<4;i++){
    joint[i].jointRot_W=joint[i-1].jointRot_W*joint[i].jointRot;
  }
  ///joint pos
  for(int i=1;i<4;i++){
    joint[i].jointPos_W=joint[i-1].jointPos_W+joint[i-1].jointRot_W*joint[i].jointPos_B;
  }
  ///joint axis_w
  for(int i=1;i<4;i++){
    joint[i].jointAxis_W=joint[i-1].jointRot_W*joint[i].jointAxis_B;
  }
  ///S matrix
    for(int i=1; i<4; i++){
      joint[i].S.setZero(6);
      joint[i].S.tail(3)=joint[i].jointAxis_W;
//      std::cout<<i<<std::endl<<joint[i].S<<std::endl;
    }
}
void initialize_body(){
  body.resize(4);
  Eigen::Vector3d TrunkCOM,ImuCOM,TrunkImuCOM,  FR_CalfCOM, FR_FootCOM, FR_CalfFootCOM;
  Eigen::Matrix3d FR_Calf_I, FR_Thigh_I, FR_Hip_I, FR_Foot_I,FR_CalfFoot_I;
  Eigen::Matrix3d Trunk_I, Imu_I, TrunkImu_I;
  const double TrunkMass=4.713 , ImuMass= 0.001, TrunkImuMass=4.714,HipMass=0.696, ThighMass=1.013, CalfMass=0.166, FootMass=0.06, CalfFootMass=0.166+0.06;
  TrunkCOM << 0.012731, 0.002186, 0.000515;
  ImuCOM << 0,0,0;
  FR_CalfCOM << 0.006435, 0.0, -0.107388;
  FR_FootCOM << 0, 0 ,-0.2;
  TrunkImuCOM = (TrunkMass*TrunkCOM+ImuMass*ImuCOM)/TrunkImuMass;
  FR_CalfFootCOM = (CalfMass*FR_CalfCOM+FootMass*FR_FootCOM)/CalfFootMass;
  body[0].mass=TrunkImuMass;
  body[1].mass=HipMass;
  body[2].mass=ThighMass;
  body[3].mass=CalfFootMass;
  body[0].comPos_B= TrunkImuCOM;
  body[1].comPos_B<< -0.003311, -0.000635, 3.1e-05;
  body[2].comPos_B<< -0.003237, 0.022327, -0.027326;
  body[3].comPos_B= FR_CalfFootCOM;

  Eigen::VectorXd T(6),I(6);
  Eigen::VectorXd FR_H(6),FR_TH(6),FR_C(6),FR_F(6);
  T<<0.01683993, 8.3902e-05, 0.000597679, 0.056579028, 2.5134e-05, 0.064713601;
  I<< 0.0001, 0, 0, 0.000001, 0, 0.0001;
  FR_H<<0.000469246, 9.409e-06, -3.42e-07, 0.00080749, 4.66e-07, 0.000552929;
  FR_TH<<0.005529065,-4.825e-06, 0.000343869,0.005139339,-2.2448e-05,0.001367788;
  FR_C<<0.002997972, 0.0, -0.000141163, 0.003014022, 0.0, 3.2426e-05;
  FR_F<<9.6e-06, 0.0, 0.0, 9.6e-06, 0.0, 9.6e-06;
  Trunk_I = Urdf_I(T);
  Imu_I = Urdf_I(I);
  FR_Calf_I = Urdf_I(FR_C);
  FR_Foot_I = Urdf_I(FR_F);
  Eigen::Vector3d r1=TrunkCOM-TrunkImuCOM;
  Eigen::Vector3d r2=ImuCOM-TrunkImuCOM;
  TrunkImu_I= Trunk_I+Imu_I-TrunkMass*skew(r1)*skew(r1)-ImuMass*skew(r2)*skew(r2);
  Eigen::Vector3d r3=FR_CalfCOM-FR_CalfFootCOM;
  Eigen::Vector3d r4=FR_FootCOM-FR_CalfFootCOM;
  FR_CalfFoot_I=FR_Calf_I+FR_Foot_I-CalfMass*skew(r3)*skew(r3)-FootMass*skew(r4)*skew(r4);
  body[0].I_B = TrunkImu_I;
  body[1].I_B = Urdf_I(FR_H);
  body[2].I_B = Urdf_I(FR_TH);
  body[3].I_B = FR_CalfFoot_I;
}

void forward_body(){
  ///body com pos
  for (int i=0; i<4; i++){
    body[i].comPos_W=joint[i].jointPos_W+joint[i].jointRot_W*body[i].comPos_B;
  }
  ///Inertia to world frame
  for(int i=0; i<4; i++){
    body[i].I_W= joint[i].jointRot_W * body[i].I_B * joint[i].jointRot_W.transpose();
  }
}

void composite(){
  comp.resize(4);
  Eigen::Vector3d temp_com, rbc;
  Eigen::Matrix3d temp_I, I3;
  I3.setIdentity();
  double temp_mass=0;
  temp_I.setZero();
  temp_com.setZero();
  for(int i = 0 ; i<4; i++){
    comp[3-i].c_mass=temp_mass+body[3-i].mass;
    comp[3-i].c_com_W=(temp_com*temp_mass+body[3-i].comPos_W*body[3-i].mass)/comp[3-i].c_mass;
    rbc= comp[3-i].c_com_W-joint[3-i].jointPos_W;
    Eigen::Vector3d r1=temp_com-comp[3-i].c_com_W;
    Eigen::Vector3d r2=body[3-i].comPos_W-comp[3-i].c_com_W;
    comp[3-i].c_I_W= temp_I+body[3-i].I_W-temp_mass*skew(r1)*skew(r1)-body[3-i].mass*skew(r2)*skew(r2);
    comp[3-i].M_spatial.setZero(6,6);
    comp[3-i].M_spatial.block(0,0,3,3)=comp[3-i].c_mass*I3;
    comp[3-i].M_spatial.block(0,3,3,3)=-comp[3-i].c_mass*skew(rbc);
    comp[3-i].M_spatial.block(3,0,3,3)=comp[3-i].c_mass*skew(rbc);
    comp[3-i].M_spatial.block(3,3,3,3)=comp[3-i].c_I_W-comp[3-i].c_mass*skew(rbc)*skew(rbc);
    temp_com=comp[3-i].c_com_W;
    temp_mass=comp[3-i].c_mass;
    temp_I=comp[3-i].c_I_W;
  }
}

/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {
  Eigen::MatrixXd M(9,9);
  M.setZero();
  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  initialize_joint(gc);
  forward_joint(gc);
  initialize_body();
  forward_body();
  composite();
  M.block(0,0,6,6)=comp[0].M_spatial;
  for(int i=1; i<4; i++){
    Eigen::MatrixXd acc(6,6);
    Eigen::Vector3d rij;
    rij=joint[0].jointPos_W-joint[i].jointPos_W;
    acc.setIdentity();
    acc.block(0,3,3,3)=skew(rij);
    M.block(5+i,0,1,6)=joint[i].S.transpose()*comp[i].M_spatial*acc;
    M.block(0,5+i,6,1)=M.block(5+i,0,1,6).transpose();
  }
  for(int i=1; i<4; i++){
    for(int j=i; j<4;j++){
      Eigen::MatrixXd acc(6,6);
      Eigen::Vector3d rij;
      rij=joint[j].jointPos_W-joint[i].jointPos_W;
      acc.setIdentity();
      acc.block(0,3,3,3)=-skew(rij);
      M.block(5+j,5+i,1,1)=joint[j].S.transpose()*comp[j].M_spatial*acc*joint[i].S;
      M.block(5+i,5+j,1,1)=joint[j].S.transpose()*comp[j].M_spatial*acc*joint[i].S;
    }
  }
  return M;
}