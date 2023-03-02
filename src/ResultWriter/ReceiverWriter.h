#ifndef RESULTWRITER_RECEIVERWRITER_H_
#define RESULTWRITER_RECEIVERWRITER_H_

#include <string_view>
#include <vector>


#include "async/Module.h"
#include <Eigen/Dense>
#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/Lut.hpp>
#include <Kernels/Receiver.h>
#include <Modules/Module.h>
#include <Monitoring/Stopwatch.h>

struct LocalIntegrationData;
struct GlobalData;
namespace seissol::writer {
Eigen::Vector3d parseReceiverLine(const std::string& line);
std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName);

struct ReceiverWriterInitParam {
  int timestep; //?
};

struct ReceiverWriterParam {
  double time;
};

class ReceiverWriterExecutor {
public:
    enum class BufferIds {
      // TODO(Lukas) Is this necessary? Can be enum class?
      POINTS
    };

  ReceiverWriterExecutor() : comm(MPI_COMM_NULL) {}

  void execInit(const async::ExecInfo& info,
                const ReceiverWriterInitParam& param);

  void exec(const async::ExecInfo& info,
            const ReceiverWriterParam& param);

  void finalize() {
    stopwatch.printTime("Time receiver writer backend:", comm);
  }
  private:
     MPI_Comm comm;

     Stopwatch stopwatch;

     // https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
     std::vector<Eigen::Vector3d> points;
};



    class ReceiverWriter :
        private async::Module<ReceiverWriterExecutor, ReceiverWriterInitParam, ReceiverWriterParam>,
        public seissol::Module {
    public:
      void setUp() override {
        setExecutor(executor);
        // TODO(Lukas) Pin
        // TODO(Lukas) When is this called?
      }

      void init(std::string receiverFileName,
                std::string fileNamePrefix,
                double syncPointInterval,
                double samplingInterval,
                bool computeRotation,
                const MeshReader& mesh,
                const seissol::initializers::Lut& ltsLut,
                const seissol::initializers::LTS& lts,
                const GlobalData* global);

      void addPoints(
          const MeshReader& mesh,
          const seissol::initializers::Lut& ltsLut,
          const seissol::initializers::LTS& lts,
          const GlobalData* global);

      kernels::ReceiverCluster* receiverCluster(unsigned clusterId, LayerType layer) {
        assert(layer != Ghost);
        assert(m_receiverClusters.find(layer) != m_receiverClusters.end());
        auto& clusters = m_receiverClusters[layer];
        if (clusterId < clusters.size()) {
          return &clusters[clusterId];
        }
        return nullptr;
      }

      void write(double time);
      void close();
      void tearDown() override;

      //
      // Hooks
      //
      void syncPoint(double) override;

    private:
      [[nodiscard]] std::string fileName(unsigned pointId) const;
      void writeHeader(unsigned pointId, Eigen::Vector3d const& point);

      std::vector<Eigen::Vector3d> points;
      std::string m_receiverFileName;
      std::string m_fileNamePrefix;
      double      m_samplingInterval;
      bool        m_computeRotation;
      // Map needed because LayerType enum casts weirdly to int.
      std::unordered_map<LayerType, std::vector<kernels::ReceiverCluster>> m_receiverClusters;
      Stopwatch   stopwatch;

      ReceiverWriterExecutor executor;
    };
  } // namespace seissol::writer

#endif
