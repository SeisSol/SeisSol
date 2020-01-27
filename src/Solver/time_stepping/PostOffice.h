#ifndef SEISSOL_POSTOFFICE_H
#define SEISSOL_POSTOFFICE_H

#include <queue>

namespace seissol {
namespace time_stepping {



class ParcelFactory {
private:
  std::queue<Parcel> emptyParcels;
  size_t parcelSize;
  int parcelCounter;
public:
  explicit ParcelFactory(size_t size) : parcelSize(parcelSize),
                                        parcelCounter(0) {}

  Parcel make() {
    Parcel parcel;
    if (emptyParcels.empty()) {
      parcel = emptyParcels.front();
      emptyParcels.pop();
    } else {
        constexpr int maxParcels = 2; // TODO(Lukas) Bit random
        assert(parcelCounter < maxParcels);
        void *mem = allocator.allocateMemory(parcelSize, PAGESIZE_HEAP, seissol::memory::HighBandwidth);
        parcel.payload = static_cast<real *>(mem);
        ++parcelCounter;
      }
      parcel.status = ParcelStatus::Filling;
      return parcel;
    }
};

struct Parcel {
  real* payload;
};

using ParcelHandle = int;

class Buffer {
private:
    real* data;
    size_t size;
    bool inUse = false;
    memory::ManagedAllocator allocator;

public:
  explicit Buffer(size_t size)
    : size(size) {
    data = static_cast<real*>(
              allocator.allocateMemory(size, PAGESIZE_HEAP, memory::HighBandwidth)
            );
  }

  bool isInUse() const {
    return inUse;
  }

  void setInUse(bool use) {
    inUse = use;
  }
};

class LocalBufferView {
private:
    Buffer* buffer;
public:
    LocalBufferView(Buffer* buffer)
      : buffer(buffer) {}

    void send() {
      assert(!isInUse());
      // MPI_Isend ...
      buffer->setInUse(true);
    }

    bool isInUse() {
      return buffer->isInUse();
    }

    bool test() {
      // MPI_Test ...
      // If returns true
      // buffer->setInUse(false);
      return true;
    }

    bool free() {
      // if (test())
      return true;
    }
};

class NeighborBufferView {
private:
    Buffer* buffer;
public:
    NeighborBufferView(Buffer* buffer)
        : buffer(buffer) {}

    void recv() {
      // MPI_Irecv ...
    }

    bool test() {
      // ready = MPI_Test ...
      // If ready
      // return ready
      return buffer->isInUse();
    }

    bool free() {
      // if (test)..
      buffer->setInUse(false);
      return true;
    }
};

class PostOffice {
private:
  ParcelFactory parcelFactory;
  std::queue<ParcelHandle> queue;
  seissol::memory::ManagedAllocator allocator;
  size_t parcelSize;
  std::vector<Parcel> parcels;
  ParcelHandle lastHandle = 0;

public:
  explicit PostOffice(size_t parcelSize, int numberOfBoxes)
      : parcelSize(parcelSize),
        parcels(numberOfBoxes) {

    for (auto& parcel : parcels) {
      parcel.payload = static_cast<real*>(
          allocator.allocateMemory(parcelSize,
                                   PAGESIZE_HEAP,
                                   memory::HighBandwidth));
    }
  }

  ParcelHandle getNextParcelHandle() {
    lastHandle = (lastHandle+1) % parcels.size();
    return lastHandle;
  }

  real* getStoragePointer(ParcelHandle handle) {
    return parcels[handle].payload;
  }

  void markParcelAsReady(ParcelHandle handle) {
    queue.push(handle);
  }

  bool hasParcelReady() {
    return !queue.empty();
  }

  real* getParcel() {
    assert(!queue.empty());
    const auto msg = queue.front();
    queue.pop();
    return msg;
  }
};
}
}

#endif //SEISSOL_POSTOFFICE_H
