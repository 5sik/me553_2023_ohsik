//
// Created by jemin on 23. 4. 20.
//

#include <string>
#include <Eigen/Core>

#ifndef ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_
#define ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_

inline Eigen::MatrixXd getMassMatrix (const Eigen::VectorXd& gc) {

  /// !!!!!!!!!! NO RAISIM FUNCTIONS HERE !!!!!!!!!!!!!!!!!
    Eigen::Matrix3d  i1, i2 ;
    Eigen::Matrix<double,3,2> j1_p, j1_a, j2_p, j2_a;
    float m1 = 2;
    float m2 = 5;

    i1 << 2,0,0,
        0,1,0,
        0,0,2;
    i2 << 1,0,0,
        0,1,0,
        0,0,1;
    j1_p << 1,0,
            0,0,
            0,0;
    j1_a << 0,0,
            0,0,
            0,0;
    j2_p << 1,0.5*cos(gc[1]),
            0,0,
            0,-0.5*sin(gc[1]);
    j2_a << 0,0,
            0,1,
            0,0;
    Eigen::MatrixXd j1_p_t = j1_p.transpose();
//    std::cout<< j1_p << j1_p_t << std::endl;
    Eigen::MatrixXd j1_a_t = j1_a.transpose();
    Eigen::MatrixXd j2_p_t = j2_p.transpose();
    Eigen::MatrixXd j2_a_t = j2_a.transpose();

    Eigen::Matrix2d Mass = j1_p_t*m1*j1_p + j1_a_t*i1*j1_a + j2_p_t*m2*j2_p + j2_a_t*i2*j2_a;

    return Mass;
}

#endif //ME553_2022_SOLUTIONS_MIDTERM_STUDENTID_HPP_
