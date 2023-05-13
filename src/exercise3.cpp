//
// Created by Jemin Hwangbo on 2022/04/08.
//

#include "raisim/RaisimServer.hpp"
#include "raisim/World.hpp"
#include "exercise3_20233460.hpp"
//#include "exercise3_20233460(hardcoding).hpp"
//#include "test.hpp"
//#include "test_mass.hpp"

#define _MAKE_STR(x) __MAKE_STR(x)
#define __MAKE_STR(x) #x

int main(int argc, char* argv[]) {
  auto binaryPath = raisim::Path::setFromArgv(argv[0]);

  // create raisim world
  raisim::World world; // physics world
  raisim::RaisimServer server(&world);

  // kinova
  auto aliengo = world.addArticulatedSystem(std::string(_MAKE_STR(RESOURCE_DIR)) + "/aliengo/aliengo_modified.urdf");

  // kinova configuration
  Eigen::VectorXd gc(aliengo->getGeneralizedCoordinateDim()), gv(aliengo->getDOF());
//  gc << 0, 0, 0.54, 1.0, 0.0, 0.0, 0.0, 0.03, 0.4, -0.8, -0.03, 0.4, -0.8, 0.03, -0.4, 0.8, -0.03, -0.4, 0.8; /// Jemin: I'll randomize the gc, gv when grading
  gc << 0, 0, 0.54, 1.0,0.0,0.0,0.0, 0.5, 0.1, -0.28,  0.5, 0.1, -0.28,  0.5, 0.1, -0.28,  0.5, 0.1, -0.28;
  gv << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8;
  aliengo->setState(gc, gv);


  /// if you are using an old version of Raisim, you need this line
  world.integrate1();

//  auto footFrameIndex = aliengo->getFrameIdxByLinkName("RL_thigh");
//  raisim::Vec<3> footPosition;
//  aliengo->getFramePosition(footFrameIndex,footPosition);
//
//  std::cout<< "pos : "<< footPosition.e().transpose()<< std::endl;


  std::cout<<"mass matrix which I found is \n"<< getMassMatrix(gc) <<std::endl;
//
//  std::cout<<"mass matrix should be \n"<< aliengo->getMassMatrix().e()<<std::endl;
//  std::cout<<"\n"<<"----------------------------"<<std::endl;
  std::cout<<"delta : "<<std::endl<<(getMassMatrix(gc) - aliengo->getMassMatrix().e())<<std::endl;

  if((getMassMatrix(gc) - aliengo->getMassMatrix().e()).norm() < 1e-8)
    std::cout<<"passed "<<std::endl;
  else
    std::cout<<"failed "<<std::endl;

  server.launchServer();

  while (true) {
    RS_TIMED_LOOP(1)
  }

  server.closeConnection();

  return 0;
}
