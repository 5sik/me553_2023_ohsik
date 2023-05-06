//
// Created by jemin on 5/30/22.
//

#ifndef ME553_2022_SOLUTIONS_EXERCISE7_HPP_
#define ME553_2022_SOLUTIONS_EXERCISE7_HPP_

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

Eigen::Matrix3d rot_z(double angle){
  Eigen::Matrix3d Z;
  Z << cos(angle), -sin(angle),0,
       sin(angle), cos(angle),0,
       0, 0, 1;
  return Z;
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
  ///for child
  Eigen::Vector3d vel_joint;
  Eigen::Vector3d avel_joint;

  Eigen::VectorXd S;
  Eigen::VectorXd S_dot;

  Eigen::Vector3d acc_joint;
  Eigen::Vector3d aacc_joint;
  Eigen::VectorXd acc;
  double ga;
  Eigen::VectorXd wrench;
};

struct Body{
  Eigen::Vector3d comPos_B;
  Eigen::Vector3d comPos_W;
  Eigen::Matrix3d I_B;
  Eigen::Matrix3d I_W;
  double mass;
  Eigen::MatrixXd M_spatial;
  Eigen::VectorXd bias;
  Eigen::MatrixXd M_a;
  Eigen::VectorXd b_a;
};
std::vector<Joint> joint;
std::vector<Body> body;
Eigen::Matrix3d rpy_to_rot(double r, double p, double y){

  return rot_z(y)*rot_y(p)*rot_x(r);
}
void initialize_joint(const Eigen::VectorXd& gc){
  joint.resize(6);
  joint[0].jointRot_W= rpy_to_rot(0, 3.14159265359, 0)*rot_z(gc(0));
  joint[0].jointRot=   rpy_to_rot(0, 3.14159265359, 0)*rot_z(gc(0));
  joint[1].jointRot= rpy_to_rot(-1.57079632679, 0, 3.14159265359)*rot_z(gc(1));
  joint[2].jointRot= rpy_to_rot(0, 3.14159265359, 0)*rot_z(gc(2));
  joint[3].jointRot= rpy_to_rot(-1.57079632679, 0, 3.14159265359)*rot_z(gc(3));
  joint[4].jointRot= rpy_to_rot(1.57079632679, 0, 3.14159265359)*rot_z(gc(4));
  joint[5].jointRot= rpy_to_rot(-1.57079632679, 0, 3.14159265359)*rot_z(gc(5));
  joint[0].jointPos_W << 0, 0, 0.15675;
  joint[1].jointPos_B << 0, 0.0016, -0.11875;
  joint[2].jointPos_B << 0, -0.410, 0;
  joint[3].jointPos_B << 0, 0.2073, -0.0114;
  joint[4].jointPos_B << 0, 0, -0.10375;
  joint[5].jointPos_B << 0, 0.10375, 0;
  for (int i=0; i<6; i++){
    joint[i].jointAxis_B << 0,0,1;
  }
}

void forward_joint(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv){
  ///joint orientation
  for(int i=1;i<6;i++){
    joint[i].jointRot_W=joint[i-1].jointRot_W*joint[i].jointRot;
  }
  ///joint pos
  for(int i=1;i<6;i++){
    joint[i].jointPos_W=joint[i-1].jointPos_W+joint[i-1].jointRot_W*joint[i].jointPos_B;
  }
  ///joint axis_w
  for(int i=0;i<6;i++) {
     joint[i].jointAxis_W=joint[i].jointRot_W.col(2);
  }
  ///S matrix
  for(int i=0; i<6; i++) {
    joint[i].S.setZero(6);
    joint[i].S.tail(3) = joint[i].jointAxis_W;
  }
  ///joint vel
  joint[0].vel_joint << 0,0,0;
  joint[0].avel_joint << (joint[0].S*gv(0)).tail(3);
  for(int i=1; i<6; i++){
    joint[i].vel_joint = joint[i-1].vel_joint+ skew(joint[i-1].avel_joint)*(joint[i-1].jointRot_W*joint[i].jointPos_B);
    joint[i].avel_joint = joint[i-1].avel_joint+ (joint[i].S*gv(i)).tail(3);
  }
  ///S_dot
  for(int i=0; i<6; i++){
    joint[i].S_dot.setZero(6);
    joint[i].S_dot.tail(3) = skew(joint[i].avel_joint)*joint[i].jointAxis_W;
  }
}

void initialize_body(){
  body.resize(6);
  Eigen::Vector3d endCOM, link6COM,end6COM;
  link6COM << 0, 0, -0.06; endCOM << 0, 0, -0.1600;
  end6COM= (link6COM*1.327+endCOM*0.01)/(1.327+0.01);
  body[0].mass=0.7477;
  body[1].mass=0.99;
  body[2].mass=0.6763;
  body[3].mass=0.463;
  body[4].mass=0.463;
  body[5].mass=1.327+0.01;

  body[0].comPos_B << 0, -0.002, -0.0605;
  body[1].comPos_B<< 0, -0.2065, -0.01;
  body[2].comPos_B<< 0, 0.081, -0.0086;
  body[3].comPos_B<< 0, 0.0028848942, -0.0541932613;
  body[4].comPos_B<<  0, 0.0497208855, -0.0028562765;
  body[5].comPos_B= end6COM;

  Eigen::VectorXd SIX(6),END(6);
  Eigen::Matrix3d SIX_I, END_I, SIXEND_I;
  Eigen::VectorXd ONE(6),TWO(6),THREE(6),FOUR(6), FIVE(6);
  ONE<<0.00152031725204, 0,  0, 0.00152031725204, 0, 0.00059816;
  TWO<< 0.010502207991,0, 0, 0.000792, 0 ,0.010502207991;
  THREE<<0.00142022431908, 0, 0, 0.000304335, 0, 0.00142022431908;
  FOUR<<0.0004321316048, 0, 0, 0.0004321316048, 0, 9.26e-05;
  FIVE<<0.0004321316048, 0, 0, 9.26e-05, 0, 0.0004321316048;
  SIX<<0.0004403232387, 0, 0, 0.0004403232387, 0 ,0.0007416;
  END<<0.01, 0 ,0 ,0.01, 0 ,0.01;
  SIX_I = Urdf_I(SIX);
  END_I = rpy_to_rot(3.14159265359, 0, 0)*Urdf_I(END)* rpy_to_rot(3.14159265359, 0, 0).transpose();
  Eigen::Vector3d r1=link6COM-end6COM;
  Eigen::Vector3d r2=endCOM-end6COM;
  SIXEND_I= SIX_I+END_I-1.327*skew(r1)*skew(r1)-0.01*skew(r2)*skew(r2);
  body[0].I_B = Urdf_I(ONE);
  body[1].I_B = Urdf_I(TWO);
  body[2].I_B = Urdf_I(THREE);
  body[3].I_B = Urdf_I(FOUR);
  body[4].I_B = Urdf_I(FIVE);
  body[5].I_B = SIXEND_I;
}

void forward_body(){
  Eigen::Matrix3d I3;
  I3.setIdentity();
  Eigen::Vector3d rac;
  ///body com pos
  for (int i=0; i<6; i++){
    body[i].comPos_W=joint[i].jointPos_W+joint[i].jointRot_W*body[i].comPos_B;
  }
  ///Inertia to world frame
  for(int i=0; i<6; i++){
    body[i].I_W= joint[i].jointRot_W * body[i].I_B * joint[i].jointRot_W.transpose();
  }
  ///Spatial, bias
  for (int i=0; i<6; i++){
    rac= joint[i].jointRot_W*body[i].comPos_B;
    body[i].M_spatial.setZero(6,6);
    body[i].M_spatial.block(0,0,3,3) = body[i].mass*I3;
    body[i].M_spatial.block(0,3,3,3) = -body[i].mass*skew(rac);
    body[i].M_spatial.block(3,0,3,3) =  body[i].mass*skew(rac);
    body[i].M_spatial.block(3,3,3,3) = body[i].I_W-body[i].mass*skew(rac)*skew(rac);
    body[i].bias.setZero(6);
    body[i].bias.head(3)=body[i].mass*skew(joint[i].avel_joint)*skew(joint[i].avel_joint)*rac;
    body[i].bias.tail(3)= skew(joint[i].avel_joint)*(body[i].I_W-body[i].mass*skew(rac)*skew(rac))*joint[i].avel_joint;
  }
}

void getArticulatedBody(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv){
  Eigen::Matrix3d I3;
  I3.setIdentity();
  body[5].M_a=body[5].M_spatial;
  body[5].b_a=body[5].bias;
  for (int i=1; i<6; i++){
    Eigen::MatrixXd X, X_Dot;
    X.setIdentity(6,6);
    X_Dot.setZero(6,6);
    Eigen::Vector3d rpb, vpb;
    rpb= joint[6-i].jointPos_W-joint[5-i].jointPos_W;
    vpb= joint[6-i].vel_joint-joint[5-i].vel_joint;
    Eigen::VectorXd W_p(6);
    W_p.setZero();
    X.block(3,0,3,3)=skew(rpb);
    X_Dot.block(3,0,3,3)=skew(vpb);
    W_p.head(3)=joint[5-i].vel_joint;
    W_p.tail(3)=joint[5-i].avel_joint;
    Eigen::VectorXd S(6), bba(6), S_Dot(6);
    Eigen::MatrixXd Mba(6,6);
    double tau=0;
    Mba.setZero(); S.setZero(); bba.setZero();S_Dot.setZero();
    Mba=body[6-i].M_a;
    S=joint[6-i].S;
    S_Dot=joint[6-i].S_dot;
    bba=body[6-i].b_a;
//    body[5-i].M_a=body[5-i].M_spatial+X*Mba*((S*(-1*S.transpose()*Mba*S).inverse()*(S.transpose()*Mba*X.transpose()))+X.transpose());
    body[5-i].M_a=body[5-i].M_spatial+X*Mba*X.transpose()+X*Mba*S*(S.transpose()*Mba*S).inverse()*(-S.transpose()*Mba*X.transpose());
   body[5-i].b_a=body[5-i].bias+X*(Mba*(S*(S.transpose()*Mba*S).inverse()*(tau-S.transpose()*Mba*(S_Dot*gv(6-i)+X_Dot.transpose()*W_p)-S.transpose()*bba)+S_Dot*gv(6-i)+X_Dot.transpose()*W_p)+bba);
  }
}

void getAcc(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  Eigen::Matrix3d I3;
  I3.setIdentity();
  for (int i=0; i < 6; i++) {
  Eigen::VectorXd W_p(6);
  Eigen::VectorXd W_dot(6);
  W_p.setZero();
  Eigen::MatrixXd X, X_Dot;
  X.setIdentity(6, 6);
  X_Dot.setZero(6, 6);
  Eigen::Vector3d rpb, vpb;
  Eigen::VectorXd S(6), bba(6), S_Dot(6);
  Eigen::MatrixXd Mba(6, 6);
  double tau = 0;
    if (i==0){
      W_dot << 0, 0, 9.81, 0, 0, 0;
    }else {
      W_dot = joint[i - 1].acc;
      rpb= joint[i].jointPos_W-joint[i-1].jointPos_W;
      vpb= joint[i].vel_joint-joint[i-1].vel_joint;
      X.block(3,0,3,3)=skew(rpb);
      X_Dot.block(3,0,3,3)=skew(vpb);
      W_p.head(3)=joint[i-1].vel_joint;
      W_p.tail(3)=joint[i-1].avel_joint;
    }
    Mba = body[i].M_a;
    S = joint[i].S;
    S_Dot = joint[i].S_dot;
    bba = body[i].b_a;
    joint[i].ga = (1 / (S.transpose() * Mba * S))
        * (tau - S.transpose() * Mba * (S_Dot * gv(i)+ X_Dot.transpose()*W_p+ X.transpose()*W_dot) - S.transpose() * bba);
    joint[i].acc = S*joint[i].ga+S_Dot*gv(i)+X_Dot.transpose()*W_p+X.transpose()*W_dot;
  }
}
/// do not change the name of the method
inline Eigen::VectorXd computeGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
  initialize_joint(gc);
  forward_joint(gc,gv);
  initialize_body();
  forward_body();
  getArticulatedBody(gc,gv);
  getAcc(gc,gv);
  Eigen::VectorXd ga(6);
//  std::cout<<"GA"<<joint[0].ga<<std::endl;
  for (int i=0; i<6; i++){
    ga(i)=joint[i].ga;
  }
//  std::cout<<joint[1].jointAxis_W<<std::endl;
  return ga;
}

#endif //ME553_2022_SOLUTIONS_EXERCISE7_HPP_
