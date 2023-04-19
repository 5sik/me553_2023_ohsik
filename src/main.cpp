#include <iostream>
#include "raisim/RaisimServer.hpp"

#define MAKE_STR(x) _MAKE_STR(x)
#define _MAKE_STR(x) #x

int main() {
  std::string rsc_dir = std::string(MAKE_STR(RSC_DIR));
  std::cout<<rsc_dir<<std::endl;
  /// create raisim world
  raisim::World world;
  world.setTimeStep(0.001);

  /// create objects
  world.addGround();
  auto cartPole =
      world.addArticulatedSystem("/home/ohsik/raisim/me553_2023/resource/doubleCartPole/doubleCartPole.urdf");
  cartPole->setIntegrationScheme(raisim::ArticulatedSystem::IntegrationScheme::RUNGE_KUTTA_4);

  /// cartPole state
  Eigen::VectorXd gc(cartPole->getGeneralizedCoordinateDim()), gv(cartPole->getDOF());
  gc.setZero();
  gc[1] = 0.01;
  cartPole->setGeneralizedCoordinate(gc);

  /// launch raisim server
  raisim::RaisimServer server(&world);
  server.launchServer();
  server.focusOn(cartPole);

  /// mass matrix
  gc << 0.1, 0.2, 0.3;
  gv << 0.1, 0.2, 0.3;
  cartPole->setState(gc, gv);

  for (int i = 0; i < 200000; i++) {
    std::this_thread::sleep_for(std::chrono::microseconds(10));
    server.integrateWorldThreadSafe();
  }

  server.killServer();

  return 0;
}
