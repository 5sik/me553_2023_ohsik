#pragma once


#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

Eigen::Vector3d getJointAbsPos(const Eigen::VectorXd& gc, int index){
//index
//hip #1
//tight #2
//calf #3
//foot #4
Eigen::Vector3d CalfToFoot;
Eigen::Vector3d ThighToCalf;
Eigen::Vector3d HipToThigh;
Eigen::Vector3d FBaseToHip;
Eigen::Vector3d RootP;
Eigen::Vector3d Result;
Eigen::Vector3d posFoot;
Eigen::Vector3d posCalf;
Eigen::Vector3d posThigh;
Eigen::Vector3d posHip;
CalfToFoot << 0, 0, -0.2;
ThighToCalf << 0, 0, -0.2;
HipToThigh << 0, 0.08505, 0;
FBaseToHip << -0.183,0.047,0;
RootP << gc(0), gc(1), gc(2);

Eigen::Matrix3d CalfRot;
CalfRot << cos(gc(18)), 0, sin(gc(18)),
      0, 1,           0,
      -sin(gc(18)), 0, cos(gc(18));
Eigen::Matrix3d ThighRot;
ThighRot << cos(gc(17)), 0, sin(gc(17)),
      0, 1,           0,
      -sin(gc(17)), 0, cos(gc(17));
Eigen::Matrix3d HipRot;
HipRot << 1,        0,         0,
      0, cos(gc(16)), -sin(gc(16)),
      0, sin(gc(16)),  cos(gc(16));
Eigen::Quaterniond RootA(gc(3),gc(4),gc(5),gc(6));
Eigen::Matrix3d Roota = RootA.normalized().toRotationMatrix();

  if(index==0){
    return RootP;
  }else if(index==1){ ///hip
    return RootP+Roota*FBaseToHip;
  }else if(index==2){ ///thigh
    return RootP+Roota*(FBaseToHip+HipRot*HipToThigh);
  }else if(index==3){ ///calf
    return RootP+Roota*(FBaseToHip+HipRot*(HipToThigh+ThighRot*ThighToCalf));
  }else if(index==4){///foot
    return RootP+Roota*(FBaseToHip+HipRot*(HipToThigh+ThighRot*(ThighToCalf+CalfRot*CalfToFoot)));
  }else{
    return Eigen::Vector3d::Ones();
  }

}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    //////////////////////////
    ///// Your Code Here /////
    //////////////////////////

  Eigen::Matrix3d CalfRot;
  CalfRot << cos(gc(18)), 0, sin(gc(18)),
      0, 1,           0,
      -sin(gc(18)), 0, cos(gc(18));
  Eigen::Matrix3d ThighRot;
  ThighRot << cos(gc(17)), 0, sin(gc(17)),
      0, 1,           0,
      -sin(gc(17)), 0, cos(gc(17));
  Eigen::Matrix3d HipRot;
  HipRot << 1,        0,         0,
      0, cos(gc(16)), -sin(gc(16)),
      0, sin(gc(16)),  cos(gc(16));
  Eigen::Quaterniond RootA(gc(3),gc(4),gc(5),gc(6));
  Eigen::Matrix3d Roota = RootA.normalized().toRotationMatrix();

  Eigen::Vector3d Vbase;
  Eigen::Vector3d r0e, r1e, r2e, r3e; ///position difference to end effector
  Eigen::Vector3d p1, p2, p3; ///rot axis vectors in joint frame
  Eigen::Vector3d P1, P2, P3; ///rot axis vectors in world frame
  Eigen::Vector3d w0, w11, w22, w33 ;///angular velocity of joints in world frame

  r0e = getJointAbsPos(gc,  4)-getJointAbsPos(gc,  0);
  r1e = getJointAbsPos(gc,  4)-getJointAbsPos(gc,  1);
  r2e = getJointAbsPos(gc,  4)-getJointAbsPos(gc,  2);
  r3e = getJointAbsPos(gc,  4)-getJointAbsPos(gc,  3);
  p1 << 1, 0, 0;
  p2 << 0, 1, 0;
  p3 << 0, 1, 0;
  P1= Roota*p1;
  P2= Roota*HipRot*p2;
  P3= Roota*HipRot*ThighRot*p3;
  Vbase << gv(0),gv(1),gv(2);
  w0 << gv(3),gv(4),gv(5);
  w11= gv(15)*P1;
  w22= gv(16)*P2;
  w33= gv(17)*P3;


  return Vbase+w0.cross(r0e)+w11.cross(r1e)+w22.cross(r2e)+w33.cross(r3e);
   /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////

  Eigen::Vector3d p1, p2, p3; ///rot axis vectors in joint frame
  Eigen::Vector3d P1, P2, P3; ///rot axis vectors in world frame
  Eigen::Vector3d w0, w11, w22, w33 ;///angular velocity of joints in world frame
  Eigen::Matrix3d CalfRot;
  CalfRot << cos(gc(18)), 0, sin(gc(18)),
      0, 1,           0,
      -sin(gc(18)), 0, cos(gc(18));
  Eigen::Matrix3d ThighRot;
  ThighRot << cos(gc(17)), 0, sin(gc(17)),
      0, 1,           0,
      -sin(gc(17)), 0, cos(gc(17));
  Eigen::Matrix3d HipRot;
  HipRot << 1,        0,         0,
      0, cos(gc(16)), -sin(gc(16)),
      0, sin(gc(16)),  cos(gc(16));
  Eigen::Quaterniond RootA(gc(3),gc(4),gc(5),gc(6));
  Eigen::Matrix3d Roota = RootA.normalized().toRotationMatrix();
  p1 << 1, 0, 0;
  p2 << 0, 1, 0;
  p3 << 0, 1, 0;
  P1= Roota*p1;
  P2= Roota*HipRot*p2;
  P3= Roota*HipRot*ThighRot*p3;
  w0 << gv(3),gv(4),gv(5);
  w11= gv(15)*P1;
  w22= gv(16)*P2;
  w33= gv(17)*P3;
  return w0+w11+w22+w33;

}