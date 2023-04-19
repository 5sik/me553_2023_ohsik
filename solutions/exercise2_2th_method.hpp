#pragma once
#include <string>
#include <Eigen/Core>
#include "raisim/math.hpp"

Eigen::Matrix3d RotationMatrix(std::string axis, double angle) {
    Eigen::Matrix3d RotM;

    if (axis == "x") {
        RotM << 1, 0, 0,
                0, cos(angle), -sin(angle),
                0, sin(angle), cos(angle);
    }
    else if (axis == "y") {
        RotM << cos(angle), 0, sin(angle),
                0, 1, 0,
                - sin(angle), 0, cos(angle);
    }
    else if (axis == "z") {
        RotM << cos(angle), -sin(angle), 0,
                sin(angle), cos(angle), 0,
                0, 0, 1;
    }
    return RotM;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootLinearVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {
    Eigen::Vector3d r_T2H, r_H2P, r_P2K, r_K2E;
    r_T2H<< 0.2399,-0.051,0; r_H2P<< 0,-0.083,0;
    r_P2K<< 0,0,-0.25; r_K2E<< 0,0,-0.25;

    Eigen::Vector3d r_T2E, r_H2E, r_P2E;
    Eigen::Matrix3d R_HH, R_PP, R_KK;
    raisim::Mat<3,3> R_TT;
    raisim::Vec<4> v;
    for (int i=0;i<4;i++) v[i] = gc[i+3];
    raisim::quatToRotMat(v,R_TT);
    R_HH << RotationMatrix("x",gc[7]);
    R_PP << RotationMatrix("y",gc[8]);
    R_KK << RotationMatrix("y",gc[9]);

    /// Position
    r_T2E << R_TT.e() * (r_T2H + R_HH * (r_H2P + R_PP * (r_P2K + R_KK * r_K2E)) );
    r_H2E << R_TT.e() * R_HH * (r_H2P + R_PP * (r_P2K + R_KK * r_K2E));
    r_P2E << R_TT.e() * R_HH * R_PP * (r_P2K + R_KK * r_K2E);
    r_K2E << R_TT.e() * R_HH * R_PP * R_KK * r_K2E;
    /// Angular Velocity
    Eigen::Vector3d W_TT,W_HH,W_PP,W_KK;
    W_TT << gv.segment(3,3);
    W_HH << R_TT.e() * R_HH * gv[6] * Eigen::Vector3d{1,0,0};
    W_PP << R_TT.e() * R_HH * R_PP * gv[7] * Eigen::Vector3d{0,1,0};
    W_KK << R_TT.e() * R_HH * R_PP * R_KK * gv[8] * Eigen::Vector3d{0,1,0};
    /// Linear Velocity
    Eigen::Vector3d V;
    V << gv.segment(0,3) + W_TT.cross(r_T2E) + W_HH.cross(r_H2E) + W_PP.cross(r_P2E) + W_KK.cross(r_K2E);

    return V;
}

/// do not change the name of the method
inline Eigen::Vector3d getFootAngularVelocity (const Eigen::VectorXd& gc, const Eigen::VectorXd& gv) {

    Eigen::Matrix3d R_HH, R_PP, R_KK;
    raisim::Mat<3,3> R_TT;
    raisim::Vec<4> v;
    for (int i=0;i<4;i++) v[i] = gc[i+3];
    raisim::quatToRotMat(v,R_TT);
    R_HH << RotationMatrix("x",gc[7]);
    R_PP << RotationMatrix("y",gc[8]);
    R_KK << RotationMatrix("y",gc[9]);

    /// Angular Velocity
    Eigen::Vector3d W_TT,W_HH,W_PP,W_KK,W;
    W_TT << gv.segment(3,3);
    W_HH << R_TT.e() * R_HH * gv[6] * Eigen::Vector3d{1,0,0};
    W_PP << R_TT.e() * R_HH * R_PP * gv[7] * Eigen::Vector3d{0,1,0};
    W_KK << R_TT.e() * R_HH * R_PP * R_KK * gv[8] * Eigen::Vector3d{0,1,0};

    W = W_TT + W_HH + W_PP + W_KK;

    return W;
}