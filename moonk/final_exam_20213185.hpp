#pragma once
Eigen::Matrix3d quat_to_rot(Eigen::Vector4d quat){
  double q0=quat[0];
  double q1=quat[1];
  double q2=quat[2];
  double q3=quat[3];
  Eigen::Matrix3d R;
  R << 2*(pow(q0,2)+pow(q1,2))-1, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2),
          2*(q1*q2+q0*q3), 2*(pow(q0,2)+pow(q2,2))-1, 2*(q2*q3-q0*q1),
          2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 2*(pow(q0,2)+pow(q3,2))-1;

  return R;
}
static Eigen::Matrix3d hat(Eigen::Vector3d v){
  Eigen::Matrix3d vv;
  vv << 0, -v(2), v(1),
          v(2), 0, -v(0),
          -v(1), v(0),0;
  return vv;
}
struct Joint{
    Eigen::Vector3d rpy;
    Eigen::Vector3d xyz, xyz_axis;
    Eigen::Vector3d axis, pos, v,w;
    Eigen::VectorXd a_parent,a_child;
    Eigen::VectorXd vp,ap;
    Eigen::MatrixXd S,X,Xd,Sd;
    Eigen::MatrixXd M_spatial;
    Eigen::VectorXd F_fictitious, Force_sum;
    Eigen::Vector3d ga;
    double type; //0 if revolute, 1 if prismatic

};
struct Body{
    double mass;
    Eigen::Vector3d xyz; //COM의 위치
    Eigen::Matrix3d inertia;
    Eigen::Matrix3d inertia_w;
    Eigen::Matrix3d R;
    Eigen::VectorXd bp;
    Eigen::Vector3d rac_w;
    Eigen::MatrixXd Mp;
    Eigen::Vector3d pos_com_w;
    Eigen::Vector3d w;
};
struct Compos_body{
    double mass;
    Eigen::Vector3d xyz_com; //COM의 위치
    Eigen::Matrix3d I_w;
    Eigen::MatrixXd M_compos;
};
struct MpAba{
    Eigen::Matrix<double,6,6> M;
};
struct bpAba{
    Eigen::VectorXd b;
};
std::vector<Joint> joint;
std::vector<Body> body;
std::vector<Compos_body> Compos_body;
Eigen::VectorXd gc(12),gv(9),gf(9);
Eigen::MatrixXd angles(4,5);
Eigen::MatrixXd w_sphere(3,5), forces_sphere(3,5);

std::vector<MpAba> MpAba;
std::vector<bpAba> bpAba;
Eigen::Matrix3d Rotx(double angle) {
  Eigen::Matrix3d R_x;

  R_x << 1, 0, 0,
          0, cos(angle), -sin(angle),
          0, sin(angle), cos(angle);

  return R_x;
}
Eigen::Matrix3d Roty(double angle) {
  Eigen::Matrix3d R_y;

  R_y << cos(angle), 0, sin(angle),
          0, 1, 0,
          - sin(angle), 0, cos(angle);

  return R_y;
}
Eigen::Matrix3d Rotz(double angle) {
  Eigen::Matrix3d R_z;
  R_z.setZero();

  R_z << cos(angle), -sin(angle), 0,
          sin(angle), cos(angle), 0,
          0, 0, 1;

  return R_z;
}
Eigen::Matrix3d rpy_to_RotationM(Eigen::Vector3d rpy) {
  double x=rpy(0);
  double y=rpy(1);
  double z=rpy(2);
  Eigen::Matrix3d R_roll;
  Eigen::Matrix3d R_pitch;
  Eigen::Matrix3d R_yaw;
  R_roll << 1, 0, 0,
          0, cos(x), -sin(x),
          0, sin(x), cos(x);
  R_pitch << cos(y), 0, sin(y),
          0, 1, 0,
          - sin(y), 0, cos(y);
  R_yaw << cos(z), -sin(z), 0,
          sin(z), cos(z), 0,
          0, 0, 1;

  return R_yaw * R_pitch * R_roll;
}
Eigen::Matrix3d getRotItoI_(int i){
  double angle;
  Eigen::Matrix3d R,I3;
  I3.setIdentity();
  R.setIdentity();
  if(i==1) R=R;
  else R= quat_to_rot(angles.col(i));
  return R;
}
Eigen::Matrix3d getRotI_toII(int i){
  Eigen::Matrix3d R;
  R.setIdentity();
  R= rpy_to_RotationM(joint[i+1].rpy);
  return R;
}
Eigen::Matrix3d getRotI(int i){
  Eigen::Matrix3d R;
  R.setIdentity();
  for (int j=1;j<i;j++){
    R=R*getRotItoI_(j)* getRotI_toII(j);
  }
  return R;
}
Eigen::Matrix3d getRotI_(int i){
  Eigen::Matrix3d R;
  R.setIdentity();
  R= getRotI(i) * getRotItoI_(i);
  return R;
}
Eigen::Vector3d getPosI_toII(int i){
  return joint[i+1].xyz;
}
Eigen::Vector3d getPosI_toIcom(int i){
  return body[i].xyz;
}
Eigen::Vector3d getPosI(int i){
  Eigen::Vector3d Vec;
  Vec=joint[1].xyz;
  for (int j = 1; j < i; j++) {
    Vec += getRotI_(j) * getPosI_toII(j);
  }
  return Vec;
}
Eigen::Vector3d getPosI_(int i){

  return getPosI(i);
}
Eigen::MatrixXd getJp(int i){
  Eigen::MatrixXd J_p(3,6);
  J_p.setZero();
  for (int j=1 ;j<i+1 and j<7 ;j++){
    J_p.block(0,j-1,3,1)
            =-hat(getPosI(i)- getPosI(j))*getRotI(j).col(2);
  }
  return J_p;
}
Eigen::MatrixXd getJa(int i){
  Eigen::MatrixXd J_a(3,6);
  J_a.setZero();
  for (int j=1 ;j<i+1 and j<7 ;j++){
    J_a.block(0,j-1,3,1)=getRotI(j).col(2);
  }
  return J_a;
}
Eigen::Vector3d getPosIcom(int i){
  Eigen::Vector3d Vec;
  Vec= getPosI(i)+ getRotI_(i)*getPosI_toIcom(i);
  return Vec;
}
Eigen::Vector3d quat_to_rpy(Eigen::VectorXd quat){
  Eigen::Vector3d rpy;
  double w=quat(0);
  double x=quat(1);
  double y=quat(2);
  double z=quat(3);
  double roll,pitch,yaw;

  double sinr_cosp = 2 * (w * x + y * z);
  double cosr_cosp = 1 - 2 * (x * x + y * y);
  roll = std::atan2(sinr_cosp, cosr_cosp);

  // pitch (y-axis rotation)
  double sinp = 2 * (w * y - z * x);
  if (std::abs(sinp) >= 1)
    pitch = std::copysign(M_PI / 2, sinp); // use 90 degrees if out of range
  else
    pitch = std::asin(sinp);

  // yaw (z-axis rotation)
  double siny_cosp = 2 * (w * z + x * y);
  double cosy_cosp = 1 - 2 * (y * y + z * z);
  yaw = std::atan2(siny_cosp, cosy_cosp);

  rpy<<roll,pitch,yaw;

  return rpy;
}
void initialize_joint(){
  joint.resize(5);
  // joint type 지정 0 for revolute, 1 for prismatic
  ///TODO
  angles.setZero();
  angles.block(0,2,4,1)<<gc(0),gc(1),gc(2),gc(3);
  angles.block(0,3,4,1)<<gc(4),gc(5),gc(6),gc(7);
  angles.block(0,4,4,1)<<gc(8),gc(9),gc(10),gc(11);
  w_sphere.setZero();
  w_sphere.block(0,2,3,1)<<gv(0),gv(1),gv(2);
  w_sphere.block(0,3,3,1)<<gv(3),gv(4),gv(5);
  w_sphere.block(0,4,3,1)<<gv(6),gv(7),gv(8);

  forces_sphere.setZero();
  forces_sphere.block(0,2,3,1)<<gf(0),gf(1),gf(2);
  forces_sphere.block(0,3,3,1)<<gf(3),gf(4),gf(5);
  forces_sphere.block(0,4,3,1)<<gf(6),gf(7),gf(8);
  joint[1].rpy = {0, 0, 0};
  joint[1].xyz = {0,0,10};

  joint[2].rpy = {0,0,0};
  joint[2].xyz = {0,0,-0.24};
  joint[3].rpy = {0, 0, 0};
  joint[3].xyz = {0,0,-0.24};
  joint[4].rpy = {0, 0, 0};
  joint[4].xyz = {0,0,-0.24};


  for (int i=1;i<=4;i++){
    joint[i].X.setIdentity(6,6);
    joint[i].Xd.setZero(6,6);
    joint[i].pos= getPosI(i);
    joint[i].v= getJp(i)*gv.tail(6);
    joint[i].axis= getRotI_(i)*joint[i].xyz_axis;
    if (i == 1) joint[1].w = {0, 0, 0};
    else joint[i].w = joint[i - 1].w + getRotI_(i)*w_sphere.col(i);
    joint[i].S.setZero(6,3);
    joint[i].Sd.setZero(6,3);
    joint[i].S.block(3, 0, 3, 3) = getRotI_(i);
    joint[i].Sd.block(3, 0, 3, 3) = hat(joint[i].w)*getRotI_(i);
    if (i>0) {
      joint[i].X.block(3, 0, 3, 3) = hat(getPosI(i) - getPosI(i-1));
//      joint[i].Xd.block(3, 0, 3, 3) = hat(joint[i].v - joint[i-1].v);
      joint[i].Xd.block(3, 0, 3, 3) = hat(hat(joint[i-1].w)* (getPosI(i)- getPosI(i-1)));
    }
    joint[i].vp.setZero(6);
    joint[i].ap.setZero(6);
    joint[i].vp.block(0,0,3,1)=joint[i].v;
    joint[i].vp.block(3,0,3,1)=joint[i].w;
  }

  joint[1].vp<<0,0,0,0,0,0;
  joint[1].ap<<0,0,9.81,0,0,0;

  for (int i=1;i<=4;i++){
    joint[i].a_child.setZero(6);
    joint[i].a_parent.setZero(6);
    joint[i].M_spatial.setZero(6,6);
    joint[i].F_fictitious.setZero(6,1);
  }

}
void initialize_body() {
  body.resize(5);

  body[1].mass = 1;
  body[1].xyz = {0, 0, -0.12};
  body[1].inertia << 0.002, 0, 0, 0, 0.002, 0, 0, 0, 0.002;

  body[2].mass = 1;
  body[2].xyz = {0, 0, -0.12};
  body[2].inertia << 0.002, 0, 0, 0, 0.002, 0, 0, 0, 0.002;

  body[3].mass = 1;
  body[3].xyz = {0, 0, -0.12};
  body[3].inertia << 0.002, 0, 0, 0, 0.002, 0, 0, 0, 0.002;

  body[4].mass = 1;
  body[4].xyz = {0, 0, -0.12};
  body[4].inertia << 0.002, 0, 0, 0, 0.002, 0, 0, 0, 0.002;


  Eigen::Vector3d w_spherical;
  for (int i = 1; i <= 4; i++) {
    body[i].R = getRotI_(i);
    body[i].inertia_w = getRotI_(i) * body[i].inertia * getRotI_(i).transpose();
    body[i].pos_com_w = getPosI_(i) + getRotI_(i) * body[i].xyz;
    body[i].Mp.setZero(6,6);
    body[i].bp.setZero(6);
    body[i].rac_w= getRotI_(i)*getPosI_toIcom(i);
    body[i].Mp.block(0,0,3,3)=body[i].mass*Eigen::Matrix3d::Identity();
    body[i].Mp.block(0,3,3,3)=-body[i].mass*hat(body[i].rac_w);
    body[i].Mp.block(3,0,3,3)=body[i].mass*hat(body[i].rac_w);
    body[i].Mp.block(3,3,3,3)=body[i].inertia_w-body[i].mass*hat(body[i].rac_w)*hat(body[i].rac_w);
    body[i].bp.segment(0,3)=body[i].mass*hat(joint[i].w)*hat(joint[i].w)*body[i].rac_w;
    body[i].bp.segment(3,3)=hat(joint[i].w)*(body[i].inertia_w-body[i].mass*hat(body[i].rac_w)*hat(body[i].rac_w))*joint[i].w;
    if (i == 1) body[1].w = {0, 0, 0};
    else body[i].w = body[i - 1].w + w_sphere.col(i);
  }
}
void getSpatial(){
  Eigen::Vector3d rac;
  Eigen::Matrix3d I3;
  I3.setIdentity();
  for (int i=2;i<=4;i++){
    rac=body[i].R*body[i].xyz;
    joint[i].M_spatial.block(0,0,3,3)=body[i].mass* I3;
    joint[i].M_spatial.block(0,3,3,3)=-body[i].mass* hat(rac);
    joint[i].M_spatial.block(3,0,3,3)=body[i].mass* hat(rac);
    joint[i].M_spatial.block(3,3,3,3)=body[i].inertia_w-body[i].mass* hat(rac)* hat(rac);
  }
}
void getFictitious(){
  Eigen::Vector3d rac;
  for (int i=2;i<=4;i++){
    rac= getPosIcom(i)- getPosI_(i);
    joint[i].F_fictitious.block(0,0,3,1)=body[i].mass*hat(joint[i].w)*hat(joint[i].w)*rac;
    joint[i].F_fictitious.block(3,0,3,1)=hat(joint[i].w)*(body[i].inertia_w-body[i].mass* hat(rac)* hat(rac))*joint[i].w;
  }
}
void initializeCompos_body(){
  Compos_body.resize(5);
  Eigen::Vector3d rac;
  Compos_body[4].mass=body[4].mass;
  Compos_body[4].xyz_com=body[4].pos_com_w;
  Compos_body[4].I_w=body[4].inertia_w;

  for (int i=3;i>=1;i--){
    Compos_body[i].mass=body[i].mass+Compos_body[i+1].mass;
    Compos_body[i].xyz_com=(body[i].mass*body[i].pos_com_w+Compos_body[i+1].mass*Compos_body[i+1].xyz_com)/Compos_body[i].mass;
    Compos_body[i].I_w= body[i].inertia_w +Compos_body[i+1].I_w
                        -body[i].mass*hat(body[i].pos_com_w-Compos_body[i].xyz_com)*hat(body[i].pos_com_w-Compos_body[i].xyz_com)
                        -Compos_body[i+1].mass*hat(Compos_body[i+1].xyz_com-Compos_body[i].xyz_com)*hat(Compos_body[i+1].xyz_com-Compos_body[i].xyz_com);
  }
  for (int i=4;i>=1;i--){
    Compos_body[i].M_compos.setZero(6,6);
    rac=Compos_body[i].xyz_com-joint[i].pos;
    Compos_body[i].M_compos.block(0,0,3,3)=Compos_body[i].mass* Eigen::Matrix3d::Identity();
    Compos_body[i].M_compos.block(0,3,3,3)=-Compos_body[i].mass* hat(rac);
    Compos_body[i].M_compos.block(3,0,3,3)=Compos_body[i].mass* hat(rac);
    Compos_body[i].M_compos.block(3,3,3,3)=Compos_body[i].I_w-Compos_body[i].mass* hat(rac)* hat(rac);
  }

}
Eigen::MatrixXd CRBA(){
  Eigen::MatrixXd M(9,9);
  Eigen::MatrixXd M_type2(6,6);
  Eigen::MatrixXd M_compos_j(6,6);
  Eigen::MatrixXd M_type1(3,6);
  Eigen::Matrix3d I3;
  Eigen::MatrixXd X(6,6);
  I3.setIdentity();
  Eigen::Vector3d rac;
  Eigen::VectorXd base_acc(6);
  Eigen::Vector3d rbj,rij;

  //type 3: joint와 joint의 관계
  for (int i=2;i<=4;i++){
    for (int j=i;j<=4;j++){
      X.setIdentity();
      rij=joint[j].pos-joint[i].pos;
      X.block(3,0,3,3)=hat(rij);
      M.block(3*(j-2),3*(i-2),3,3)=(joint[j].S.transpose()*Compos_body[j].M_compos*X.transpose()*joint[i].S);
    }
  }
  for (int i=3;i<=8;i++){
    for (int j=0;j<=2;j++){
      M(j,i)=M(i,j);
    }
  }
  for (int i=6;i<=8;i++){
    for (int j=3;j<=5;j++){
      M(j,i)=M(i,j);
    }
  }
//  std::cout<<"M: "<<M<<std::endl;
  return M;
}
Eigen::VectorXd RNE(){
  Eigen::VectorXd b(9);
  Eigen::Vector3d alpha, rij;
  Eigen::VectorXd cori;
  cori.setZero(6);
  joint[1].a_child<<0,0,9.81,0,0,0;
  // 아래로 가면서 a와 alpha들을 쫙 구한다.
  for (int i=1;i<4;i++) {
    alpha = joint[i].a_child.tail(3);
    rij = getPosI(i+1) - getPosI(i);

    // 같은 바디에서 이동

    joint[i + 1].a_parent.head(3) = joint[i].a_child.head(3) + hat(alpha) * rij + hat(joint[i].w) * hat(joint[i].w) * rij;// + hat(body[i].w)*cori.head(3);
    joint[i + 1].a_parent.tail(3) = joint[i].a_child.tail(3);

    // parent 에서 child로 이동
    joint[i + 1].a_child = joint[i + 1].a_parent + joint[i + 1].Sd * w_sphere.col(i+1);

//    std::cout<<"Sd of joint "<<i+1<<": "<< joint[i+1].Sd<<std::endl;
//    std::cout<<"w of link "<<i+1<<": "<< body[i+1].w<<std::endl;
//    std::cout<<"acc of parent "<<i+1<<" :"<<joint[i+1].a_parent<<std::endl;
//    std::cout<<"acc of child "<<i+1<<" :"<<joint[i+1].a_child<<std::endl;
  }
//  joint[4].a_child.head(3)=joint[3].a_child.head(3)+hat(joint[3].a_child.tail(3))*(getPosI(4)- getPosI(3))

  getSpatial();
  getFictitious();
  joint[4].Force_sum=joint[4].M_spatial*joint[4].a_child+joint[4].F_fictitious;

  for (int i=3;i>=1;i--) {
    Eigen::MatrixXd Transport(6, 6);
    Transport.setIdentity(6, 6);
    Transport.block(3, 0, 3, 3) = hat(getPosI_(i+1) - getPosI_(i));
    joint[i].Force_sum =
            joint[i].M_spatial * joint[i].a_child + joint[i].F_fictitious + Transport * (joint[i + 1].Force_sum);
  }
  b.segment(6,3)=(joint[4].S.transpose()*joint[4].Force_sum);
  b.segment(3,3)=(joint[3].S.transpose()*joint[3].Force_sum);
  b.segment(0,3)=(joint[2].S.transpose()*joint[2].Force_sum);

  return b;
}
void Step2(){
  MpAba.resize(5);
  bpAba.resize(5);
  Eigen::MatrixXd S(6,6),Sd(6,3);
  Eigen::VectorXd bBA(6);
  Eigen::MatrixXd MBA(6,6), XBP(6,6),XBPd(6,6);
  Eigen::Vector3d uB;
  for(int i=4;i>0;i--){
    MpAba[i].M.setZero(6,6);
    bpAba[i].b.setZero(6);
    Eigen::VectorXd v(6);
    v.block(0,0,3,1)=joint[i].v;
    v.block(3,0,3,1)=joint[i].w;

    if(i==4){
      MpAba[i].M=body[i].Mp;
      bpAba[i].b=body[i].bp;
    }
    else{
      S=joint[i+1].S;
      Sd=joint[i+1].Sd;
      XBP=joint[i+1].X;
      XBPd=joint[i+1].Xd;
      bBA=bpAba[i+1].b;
      MBA=MpAba[i+1].M;
      uB=w_sphere.col(i);

      MpAba[i].M=body[i].Mp+XBP*MBA*XBP.transpose()+XBP*MBA*S*(S.transpose()*MBA*S).inverse()*(-S.transpose()*MBA*XBP.transpose());
      bpAba[i].b=body[i].bp+XBP*(MBA*(Sd*uB+XBPd.transpose()*joint[i].vp)+bBA)
                 +XBP*MBA*S*(S.transpose()*MBA*S).inverse()*(forces_sphere.col(i+1)-S.transpose()*(MBA*(Sd*uB+XBPd.transpose()*joint[i].vp)+bBA));
    }

//    std::cout<<"body 4부터 "<<i<<" 까지의 articulated inertia M : "<<std::endl<<MpAba[i].M<<std::endl;
//    std::cout<<"body 4부터 "<<i<<" 까지의 articulated inertia b : "<<std::endl<<bpAba[i].b<<std::endl;
  }
}
void Step3(){
  Eigen::VectorXd b(6), Vp(6), ap(6);
  Eigen::MatrixXd M(6,6), X(6,6),Xd(6,6),S(6,6),Sd(6,6);
  Eigen::Vector3d uB;
  Eigen::Vector3d tow;
  for (int i=1;i<=3;i++){
    S=joint[i+1].S;
    Sd=joint[i+1].Sd;
    X=joint[i+1].X;
    Xd=joint[i+1].Xd;
    b=bpAba[i+1].b;
    M=MpAba[i+1].M;
    Vp=joint[i].vp;
    ap=joint[i].ap;
    uB=w_sphere.col(i+1);
    tow=forces_sphere.col(i+1);
//    std::cout<<"Vp0: "<<Vp<<std::endl;
    joint[i+1].ga=(S.transpose()*M*S).inverse()*(tow-S.transpose()*(M*(Sd*uB+Xd.transpose()*Vp+X.transpose()*ap)+b));
    joint[i+1].ap=S*joint[i+1].ga+Sd*uB+Xd.transpose()*Vp+X.transpose()*ap;
//    std::cout<<"body "<<i<<"의 M_a: "<<std::endl<<MpA[i].M<<std::endl;
//    std::cout<<"body "<<i<<"의 b_a: "<<std::endl<<bpA[i].b<<std::endl;
//    std::cout<<"joint "<<i<<" 의 rotmat: "<<std::endl<<getRotI(i)<<std::endl;
//    std::cout<<"joint "<<i<<" 의 pos : "<<std::endl<<getPosI(i)<<std::endl;
  }
}
/// do not change the name of the method
inline Eigen::MatrixXd getMassMatrixUsingCRBA (const Eigen::VectorXd& gc_) {
  Eigen::MatrixXd M(9,9);
  gc=gc_;
  initialize_joint();
  initialize_body();
  initializeCompos_body();
  M=CRBA();
  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return M;
}

/// do not change the name of the method
inline Eigen::VectorXd getNonlinearitiesUsingRNE (const Eigen::VectorXd& gc_, const Eigen::VectorXd& gv_) {
  Eigen::VectorXd b(9);
  gc=gc_;
  gv=gv_;
  initialize_joint();
  initialize_body();
  initializeCompos_body();
  b=RNE();
  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!


  return b;
}

/// do not change the name of the method
inline Eigen::MatrixXd getGaUsingABA (const Eigen::VectorXd& gc_, const Eigen::VectorXd& gv_, const Eigen::VectorXd& gf_) {

  gc<<gc_;
  gv<<gv_;
  gf<<gf_;
  initialize_joint();
  initialize_body();
  Step2();
  Step3();
  Eigen::VectorXd Ga;
  Ga.setZero(9);
  for (int i=2;i<=4;i++){
    Ga.segment(3*(i-2),3)=joint[i].ga;
  }
  std::cout<<"Ga: "<<Ga<<std::endl;
  std::cout<<"ajdflksjafksdjflkdjsafkladf"<<std::endl;
  return Ga;

  return Eigen::VectorXd::Ones(9);
}