#pragma once

#include <Eigen/Core>
#include <vector>
#include "raisim/math.hpp"

void raisim::angleAxisToRotMat(const Vec<3> &a1, const double theta, Mat<3, 3> &rotMat) {}
class Joint
{
public:
  Joint(
      double x, double y, double z,
      const std::string &axis,
      const double angle,
      const double velocity,
      Joint *parentJoint = nullptr,
      Joint *childJoint = nullptr)
  {
    jointPosition_p_ = Eigen::Vector3d(x, y, z);
    jointAxis_p_ = Eigen::Matrix3d::Identity();
    velocity_p_ = velocity;
    axis_p_ = axis;

    if (axis == "x")
    {
      jointAxis_p_ << 1, 0, 0,
          0, cos(angle), -sin(angle),
          0, sin(angle), cos(angle);
    }
    else if (axis == "y")
    {
      jointAxis_p_ << cos(angle), 0, sin(angle),
          0, 1, 0,
          -sin(angle), 0, cos(angle);
    }
    else if (axis == "z")
    {
      jointAxis_p_ << cos(angle), -sin(angle), 0,
          sin(angle), cos(angle), 0,
          0, 0, 1;
    }
    rotation_p_ = jointAxis_p_;

    parentJoint_ = parentJoint;
    childJoint_ = childJoint;
  }

  void setParentJoint(Joint *parentJoint)
  {
    parentJoint_ = parentJoint;
  }

  void setChildJoint(Joint *childJoint)
  {
    childJoint_ = childJoint;
    childJoint_->setParentJoint(this);
  }

  void setJointAxis(const Eigen::Matrix3d &jointAxis)
  {
    jointAxis_p_ = jointAxis;
    rotation_p_ = jointAxis_p_;
  }

  Eigen::Vector3d f_worldPosition(Eigen::Vector3d thisPosition)
  {
    if (parentJoint_ == nullptr)
    {
      return jointPosition_p_ + rotation_p_ * thisPosition;
    }
    else
    {
      return parentJoint_->f_worldPosition(jointPosition_p_ + rotation_p_ * thisPosition);
    }
  }

  // alias
  Eigen::Vector3d wp(Eigen::Vector3d thisPosition)
  {
    return f_worldPosition(thisPosition);
  }

  Eigen::Vector3d tp()
  {
    return f_worldPosition(Eigen::Vector3d(0, 0, 0));
  }

  Eigen::Vector3d f_worldAngularVelocity(Eigen::Vector3d thisAngularVelocity)
  {
    if (parentJoint_ == nullptr)
    {
      return jointAxis_p_ * thisAngularVelocity;
    }
    else
    {
      return parentJoint_->f_worldAngularVelocity(jointAxis_p_ * thisAngularVelocity);
    }
  }
  // alias
  Eigen::Vector3d wav()
  {
    Eigen::Vector3d _ret = velocity_p_ * jointAxis_p_.col(axis_p_ == "x" ? 0 : axis_p_ == "y" ? 1
                                                                                              : 2);
    return f_worldAngularVelocity(_ret);
  }

private:
  Eigen::Matrix3d jointAxis_p_;
  Eigen::Vector3d jointPosition_p_;
  Eigen::Matrix3d rotation_p_;

  double velocity_p_;
  std::string axis_p_;

  Joint *parentJoint_;
  Joint *childJoint_;
};

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity(const Eigen::VectorXd &gc, const Eigen::VectorXd &gv)
{

  Joint r23(
      0.0, 0.0, -0.25,
      "y", gc[9], gv[8]);

  Joint r12(
      0.0, -0.083, 0.0,
      "y", gc[8], gv[7]);

  Joint r01(
      0.2399, -0.051, 0.0,
      "x", gc[7], gv[6]);

  Joint r00(
      gc[0], gc[1], gc[2], "x", 0.0, 0.0);

  r00.setChildJoint(&r01);
  r01.setChildJoint(&r12);
  r12.setChildJoint(&r23);

  raisim::Mat<3, 3> _tmp;
  raisim::Vec<4> v;
  for (int i = 0; i < 4; i++)
    v[i] = gc[i + 3];
  raisim::quatToRotMat(v, _tmp);
  r00.setJointAxis(_tmp.e());

  Eigen::Vector3d _r03, _r13, _r23, _r33;

  Eigen::Vector3d ep = Eigen::Vector3d(0.0, 0.0, -0.25);

  _r33 = r23.wp(ep) - r23.tp();
  _r23 = r23.wp(ep) - r12.tp();
  _r13 = r23.wp(ep) - r01.tp();
  _r03 = r23.wp(ep) - r00.tp();

  return gv.segment(0, 3) + Eigen::Vector3d(gv.segment(3, 3)).cross(_r03) + r01.wav().cross(_r13) + r12.wav().cross(_r23) + r23.wav().cross(_r33);
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity(const Eigen::VectorXd &gc, const Eigen::VectorXd &gv)
{

  Joint r23(
      0.0, 0.0, -0.25,
      "y", gc[9], gv[8]);

  Joint r12(
      0.0, -0.083, 0.0,
      "y", gc[8], gv[7]);

  Joint r01(
      0.2399, -0.051, 0.0,
      "x", gc[7], gv[6]);

  Joint r00(
      gc[0], gc[1], gc[2], "x", 0.0, 0.0);

  r00.setChildJoint(&r01);
  r01.setChildJoint(&r12);
  r12.setChildJoint(&r23);

  raisim::Mat<3, 3> _tmp;
  raisim::Vec<4> v;
  for (int i = 0; i < 4; i++)
    v[i] = gc[i + 3];
  raisim::quatToRotMat(v, _tmp);
  r00.setJointAxis(_tmp.e());

  return gv.segment(3, 3) + r01.wav() + r12.wav() + r23.wav();
}