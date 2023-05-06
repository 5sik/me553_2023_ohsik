#pragma once

#include<cmath>
#include<Eigen/Core>
#include <Eigen/Dense>

/// do not change the name of the method
inline Eigen::Vector2d getGeneralizedAcceleration (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
  Eigen::Matrix2d Mass;
  Eigen::Vector2d Non;
  double sin2 = std::sin(gc(1));
  double cos2 = std::cos(gc(1));
  Mass << 7, 2.5*cos2, 2.5*cos2, 2.25;
  Non << 2.5* pow(gv(1),2)*sin2, 2.5*9.81*sin2;



  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!

  return Mass.inverse()*Non;
}
