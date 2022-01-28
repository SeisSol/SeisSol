#include "Solver/time_stepping/AbstractTimeCluster.h"
namespace seissol::unit_test {
using namespace time_stepping;

class MockTimeCluster : public time_stepping::AbstractTimeCluster {
public:
  MockTimeCluster(double maxTimeStepSize,
                  double timeTolerance,
                  long timeStepRate) :
  AbstractTimeCluster(maxTimeStepSize, timeTolerance, timeStepRate) { }

  MAKE_MOCK0(start, void(void), override);
  MAKE_MOCK0(predict, void(void), override);
  MAKE_MOCK0(correct, void(void), override);
  MAKE_MOCK1(handleAdvancedPredictionTimeMessage, void(const NeighborCluster&), override);
  MAKE_MOCK1(handleAdvancedCorrectionTimeMessage, void(const NeighborCluster&), override);
  MAKE_MOCK1(printTimeoutMessage, void(std::chrono::seconds), override);
};


TEST_CASE("TimeCluster") {
  auto cluster = MockTimeCluster(1.0, 1e-10, 1);
  // TODO(Lukas) Maybe the following two functions should be bundled.
  cluster.updateSyncTime(10);
  cluster.reset();

  SUBCASE("Cluster start synced") {
    REQUIRE(cluster.synced());
    REQUIRE(cluster.getState() == ActorState::Synced);
  }

  SUBCASE("Cluster predicts after sync") {
    REQUIRE(cluster.getNextLegalAction() == ActorAction::RestartAfterSync);
    REQUIRE_CALL(cluster, start());
    auto result = cluster.act();
    REQUIRE(result.isStateChanged);
    REQUIRE(cluster.getState() == ActorState::Corrected);

    REQUIRE(cluster.getNextLegalAction() == ActorAction::Predict);
    REQUIRE_CALL(cluster, predict());
    result = cluster.act();
    REQUIRE(result.isStateChanged);
    REQUIRE(cluster.getState() == ActorState::Predicted);
  }

}
} // namespace seissol::unit_test