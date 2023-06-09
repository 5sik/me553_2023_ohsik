//
// Created by Jemin Hwangbo on 2022/04/08.
//

#include "raisim/RaisimServer.hpp"
#include "exercise4_20233460.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world

  // kinova
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo_modified.urdf");

  // kinova configuration
  Eigen::VectorXd gc(aliengo->getGeneralizedCoordinateDim()), gv(aliengo->getDOF());
  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  aliengo->setState(gc, gv);

  /// if you are using an old version of Raisim, you need this line
  world.integrate1();
  aliengo->getMassMatrix();
  aliengo->setComputeInverseDynamics(true);


  std::string frameName = "floating_base";
  raisim::Vec<3> lin_vel, ang_vel, lin_acc, ang_acc;
  aliengo->getFrameVelocity(frameName,lin_vel);
  aliengo->getFrameAngularVelocity(frameName,ang_vel);

//    std::cout<< "Linear Velocity : "<< lin_vel.e().transpose()<< std::endl;
//    std::cout<< "Angular Velocity : "<< ang_vel.e().transpose()<< std::endl;


//  std::cout<< ((getNonlinearities(gc, gv) - aliengo->getNonlinearities({0,0,-9.81}).e())).transpose() <<std::endl;
//  std::cout<<"my nonlinearities is  \n"<<getNonlinearities(gc, gv).transpose()<<std::endl;
  std::cout<<"nonlinearities should be \n"<< aliengo->getNonlinearities({0,0,-9.81}).e().transpose()<<std::endl;

  if((getNonlinearities(gc, gv) - aliengo->getNonlinearities({0,0,-9.81}).e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  std::cout<< "\n"<< std::endl;
  aliengo->getFrameAcceleration(frameName,lin_acc);
  std::cout<< "Linear Acceleration : "<< lin_acc.e().transpose()<< std::endl;
  return 0;
}
