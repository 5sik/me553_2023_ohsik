#pragma once

#include <Eigen/Core>
#include "raisim/math.hpp"

Eigen::Matrix3d RotationMatrix(double x, double y, double z) {
    Eigen::Matrix3d R_x, R_y, R_z;
    R_x << 1, 0, 0,
            0, cos(x), -sin(x),
            0, sin(x), cos(x);
    R_y << cos(y), 0, sin(y),
            0, 1, 0,
            - sin(y), 0, cos(y);
    R_z << cos(z), -sin(z), 0,
            sin(z), cos(z), 0,
            0, 0, 1;

    return R_x * R_y * R_z;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Matrix3d ROB;
    Eigen::Vector3d wOB(gv[3],gv[4],gv[5]); /// baseframe angular velocity wrt worldframe
    Eigen::Vector3d wB1(gv[6],0,0); /// LF_HAA angular velocity wrt baseframe
    Eigen::Vector3d w12(0,gv[7],0); /// LF_HFE angular velocity wrt LF_HAA frame
    Eigen::Vector3d w23(0,gv[8],0); /// LF_ADAPTER_TO_FOOT angular velocity wrt LF_HFE frame
    ROB << -2*gc[5]*gc[5]-2*gc[6]*gc[6]+1, 2*gc[4]*gc[5]-2*gc[3]*gc[6],  2*gc[3]*gc[5]+2*gc[4]*gc[6], /// Quaternoin to Rotation matrix
            2*gc[3]*gc[6]+2*gc[4]*gc[5], -2*gc[4]*gc[4]-2*gc[6]*gc[6]+1, 2*gc[5]*gc[6]-2*gc[3]*gc[4],
            2*gc[4]*gc[6]-2*gc[3]*gc[5], 2*gc[3]*gc[4]+2*gc[5]*gc[6], -2*gc[4]*gc[4]-2*gc[5]*gc[5]+1;

    Eigen::Vector3d vOB(gv[0],gv[1],gv[2]); /// baseframe linear velocity wrt worldframe
    Eigen::Vector3d rBHAA(0.277, 0.116, 0.0); /// displacement from baseframe to LF_HAA frame
    Eigen::Vector3d rHAA1(0.0635,0.041,0); /// displacement from LF_HAA frame to LF_HFE frame
    Eigen::Vector3d r12(0.0, 0.109, -0.25); /// displacement from LF_HFE frame to LF_KFE frame
    Eigen::Vector3d r23(0.1, -0.02, 0.0); /// displacement from LF_KFE frame to LF_SHANK_TO_ADAPTER frame
    Eigen::Vector3d r3e(0.0, 0.0, -0.32125); /// displacement from LF_SHANK_TO_ADAPTER frame to LF_ADAPTER_TO_FOOT frame
    Eigen::Vector3d rBe = rBHAA+RotationMatrix(gc[7],0,0)*(rHAA1+RotationMatrix(0,gc[8],0)*(r12+RotationMatrix(0,gc[9],0)*(r23+r3e))); /// displacement from baseframe to LF_ADAPTER_TO_FOOT frame
    Eigen::Vector3d rHAAe = rHAA1+RotationMatrix(0,gc[8],0)*(r12+RotationMatrix(0,gc[9],0)*(r23+r3e)); /// displacement from LF_HAA frame to LF_ADAPTER_TO_FOOT frame
    Eigen::Vector3d r1e = r12+RotationMatrix(0,gc[9],0)*(r23+r3e); /// displacement from LF_HFE frame to LF_ADAPTER_TO_FOOT frame
    Eigen::Vector3d r2e = r23+r3e; /// displacement from LF_KFE frame to LF_ADAPTER_TO_FOOT frame
    Eigen::Vector3d vBe = wOB.cross(ROB*rBe); /// linear velocity contributed by baseframe rotation
    Eigen::Vector3d v1e = (ROB*wB1).cross(ROB*RotationMatrix(gc[7],0,0)*rHAAe); /// linear velocity contributed by LF_HAA frame rotation
    Eigen::Vector3d v2e = (ROB*RotationMatrix(gc[7],0,0)*w12).cross(ROB*RotationMatrix(gc[7],0,0)*RotationMatrix(0,gc[8],0)*r1e); /// linear velocity contributed by LF_HEF frame rotation
    Eigen::Vector3d v3e = (ROB*RotationMatrix(gc[7],0,0)*RotationMatrix(0,gc[8],0)*w23).cross(ROB*RotationMatrix(gc[7],0,0)*RotationMatrix(0,gc[8],0)*RotationMatrix(0,gc[9],0)*r2e); /// linear velocity contributed by LF_KEF frame rotation

    Eigen::Vector3d vOe = vOB+vBe+v1e+v2e+v3e;
    return vOe;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Matrix3d ROB;
    Eigen::Vector3d wOB(gv[3],gv[4],gv[5]); /// angular velocity contributed by baseframe rotation
    Eigen::Vector3d wB1(gv[6],0,0);
    Eigen::Vector3d w12(0,gv[7],0);
    Eigen::Vector3d w23(0,gv[8],0);
    ROB << -2*gc[5]*gc[5]-2*gc[6]*gc[6]+1, 2*gc[4]*gc[5]-2*gc[3]*gc[6],  2*gc[3]*gc[5]+2*gc[4]*gc[6],
            2*gc[3]*gc[6]+2*gc[4]*gc[5], -2*gc[4]*gc[4]-2*gc[6]*gc[6]+1, 2*gc[5]*gc[6]-2*gc[3]*gc[4],
            2*gc[4]*gc[6]-2*gc[3]*gc[5], 2*gc[3]*gc[4]+2*gc[5]*gc[6], -2*gc[4]*gc[4]-2*gc[5]*gc[5]+1;
    Eigen::Vector3d wO1 = ROB*wB1; /// angular velocity contributed by LF_HAA rotation
    Eigen::Vector3d wO2 = ROB*RotationMatrix(gc[7],0,0)*w12; /// angular velocity contributed by LF_HEF rotation
    Eigen::Vector3d wO3 = ROB*RotationMatrix(gc[7],0,0)*RotationMatrix(0,gc[8],0)*w23; /// angular velocity contributed by LF_KEF rotation
    Eigen::Vector3d wOE = wOB+wO1+wO2+wO3;


    return wOE;
}