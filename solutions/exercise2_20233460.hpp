#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"

/// do not change the name of the method

Eigen::Vector3d CrossProduct(const Eigen::Vector3d & omega,const Eigen::Vector3d & position)
{
    Eigen::Vector3d velocity;
    Eigen::Matrix3d omega_skew;
    omega_skew << 0,-omega[2],omega[1],
                 omega[2],0,-omega[0],
                 -omega[1],omega[0],0;
    velocity << omega_skew * position;

    return velocity; // 아 외적 함수가 있음 Eigen에 ->>>> a X b => a.cross(b)
}

inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    //////////////////////////
    ///// Your Code Here /////
    //////////////////////////

    Eigen::Vector3d x; // "axis"
    Eigen::Vector3d p; // "position"
    Eigen::Vector3d w; // "Angular velocity"
    Eigen::Vector3d v; // "Velocity"
    raisim::Mat<3,3> R ;

    x << 0,1,0;  // "Knee"
    p << 0,0,-0.25;
    w << x * gv[8];
    v << CrossProduct(w,p);
    raisim::angleAxisToRotMat(x,gc[9],R);
    p << R.e()*p;
    v << R.e()*v;

    x << 0,1,0; // "Pitch"
    p << p + Eigen::Vector3d{0,0,-0.25};
    w << x * gv[7];
    v << v + CrossProduct(w,p);
    raisim::angleAxisToRotMat(x,gc[8],R);
    p << R.e()*p;
    v << R.e()*v;

    x << 1,0,0; // "Hip"
    p << p + Eigen::Vector3d{0,-0.083,0};
    w << x * gv[6];
    v << v + CrossProduct(w,p);
    raisim::angleAxisToRotMat(x,gc[7],R);
    p << R.e()*p;
    v << R.e()*v;

    p << p + Eigen::Vector3d{0.2399, -0.051, 0};
    w << Eigen::Vector3d{gv[3],gv[4],gv[5]};
    raisim::Vec<4> v1;// for 쿼터니안
    for(int i=0;i<4;i++)
        v1[i] = gc[i+3];
    raisim::quatToRotMat(v1,R);  // 00'(Trunk) 의 Rotational Matrix
    v << R.e()*v + CrossProduct(w,R.e()*p);
    v << v + Eigen::Vector3d{gv[0],gv[1],gv[2]};


    return v; /// replace this
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  //////////////////////////
  ///// Your Code Here /////
  //////////////////////////

  /// 풀이 순서 -> 1. 다리가 안 돌아가 있다는 전제 하(Axis가 돌아가 있지 않음)에 각속도 구하기 (축 맞춰서)
  ///         -> 2. 기본 gc에 의해 축이 돌아가 있는거 고려해서 적용

  Eigen::Vector3d x; // "axis"
  Eigen::Vector3d w; // theta_dot * axis = omega ("Angular Velocity") 로 사용
  raisim::Mat<3,3> R; // gc값에 의해 돌아가는 축 고려해주는 Matrix로 사용

  x << 0,1,0; // for knee joint axis
  w << gv[8] * x; //Angular velocity of knee joint without theta change
  raisim::angleAxisToRotMat(x,gc[9],R);
  w<< R.e() * w; // omega for knee joint of Aliengo

  x << 0,1,0; // for pitch joint axis
  w << w + gv[7] * x ; //Angular velocity of pitch joint without theta change
  raisim::angleAxisToRotMat(x,gc[8],R);
  w << R.e() * w;

  x << 1,0,0; // for hip joint axis
  w << w + gv[6] * x;
  raisim::angleAxisToRotMat(x,gc[7],R);
  w << R.e() * w;

  raisim::Vec<4> v;// for 쿼터니안
  for(int i=0;i<4;i++)
    v[i] = gc[i+3];
  raisim::quatToRotMat(v,R);  // 00'(Trunk) 의 Rotational Matrix
  w << R.e() * w ; // Trunk axis 회전까지 고려
  w << w + Eigen::Vector3d{gv[3],gv[4],gv[5]}; // Trunk Angular Velocity 고려

  return w; /// replace this
}

