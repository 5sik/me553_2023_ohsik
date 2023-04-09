#pragma once
#include <Eigen/Core>
#include "raisim/math.hpp"


/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d axis;
    raisim::Mat<3, 3> rot;

    position << 0, 0, -0.25; // RL_calf_joint to RL_foot_fixed
    axis << 0, 1, 0;
    velocity << axis.cross(position) * gv(8); // velocity in RL_calf_joint

    raisim::angleAxisToRotMat(axis, gc(9), rot);
    position << rot.e() * position;
    velocity << rot.e() * velocity;
    position << position + Eigen::Vector3d{0, 0, -0.25}; // RL_thigh_joint to RL_foot_fixed
    axis << 0, 1, 0;
    velocity << velocity + axis.cross(position) * gv(7); // velocity in RRL_thigh_joint

    raisim::angleAxisToRotMat(axis, gc(8), rot);
    position << rot.e() * position;
    velocity << rot.e() * velocity;
    position << position + Eigen::Vector3d{0, -0.083, 0}; // RL_hip_joint to RL_foot_fixed
    axis << 1, 0, 0;
    velocity << velocity + axis.cross(position) * gv(6); // velocity in RL_hip_joint

    raisim::angleAxisToRotMat(axis, gc(7), rot);
    position << rot.e() * position;
    velocity << rot.e() * velocity;
    position << position + Eigen::Vector3d{0.2399, -0.051, 0}; // base to RL_foot_fixed

    raisim::quatToRotMat(gc.segment(3,4), rot);
    position << rot.e() * position;
    velocity << rot.e() * velocity;
    axis << gv.segment(3,3);
    velocity << velocity + axis.cross(position) + gv.head(3); // velocity in world
    return velocity;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d velocity;
    Eigen::Vector3d axis;
    raisim::Mat<3, 3> rot;

    axis << 0, 1, 0;
    velocity << axis * gv(8); // velocity in RL_calf_joint

    raisim::angleAxisToRotMat(axis, gc(9), rot);
    velocity << rot.e() * velocity;
    axis << 0, 1, 0;
    velocity << velocity + axis * gv(7); // velocity in RRL_thigh_joint

    raisim::angleAxisToRotMat(axis, gc(8), rot);
    velocity << rot.e() * velocity;
    axis << 1, 0, 0;
    velocity << velocity + axis * gv(6); // velocity in RL_hip_joint

    raisim::angleAxisToRotMat(axis, gc(7), rot);
    velocity << rot.e() * velocity;
    raisim::quatToRotMat(gc.segment(3,4), rot);
    velocity << rot.e() * velocity;
    velocity << velocity + gv.segment(3, 3); // velocity in world
    return velocity;
}
