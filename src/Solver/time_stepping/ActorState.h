#ifndef SEISSOL_ACTORSTATE_H
#define SEISSOL_ACTORSTATE_H

namespace seissol {
namespace time_stepping {
  enum class ActorState {
    Corrected,
    PredictedLocal,
    MaybeReady,
    Synced,
    Finished
  };
}
}

#endif //SEISSOL_ACTORSTATE_H
