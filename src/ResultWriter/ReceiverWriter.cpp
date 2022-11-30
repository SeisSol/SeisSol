#include "ReceiverWriter.h"

#include <cctype>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include <Parallel/MPI.h>
#include <Modules/Modules.h>

#include <sstream>
#include <string>
#include <fstream>
#include <regex>

namespace seissol::writer {
Eigen::Vector3d parseReceiverLine(const std::string& line) {
  std::regex rgx("\\s+");
  std::sregex_token_iterator iter(line.begin(),
                                  line.end(),
                                  rgx,
                                  -1);
  std::sregex_token_iterator end;
  Eigen::Vector3d coordinates{};
  unsigned numberOfCoordinates = 0;
  for (; iter != end; ++iter, ++numberOfCoordinates) {
    if (numberOfCoordinates >= coordinates.size()) {
      throw std::runtime_error("Too many coordinates in line " + line + ".");
    }
    coordinates[numberOfCoordinates] = std::stod(*iter);
  }
  if (numberOfCoordinates != coordinates.size()) {
    throw std::runtime_error("To few coordinates in line " + line + ".");
  }
  return coordinates;
}

std::vector<Eigen::Vector3d> parseReceiverFile(const std::string& receiverFileName) {
  std::vector<Eigen::Vector3d> points{};

  std::ifstream file{receiverFileName};
  std::string line{};
  while (std::getline(file, line)) {
    bool onlyWhiteSpace = std::all_of(line.begin(), line.end(), [](auto& c) {
      return std::isspace(c); // This lambda is needed (opt. argument)
    });
    if (!onlyWhiteSpace) points.emplace_back(parseReceiverLine(line));
  }
  return points;
}

void ReceiverWriterExecutor::execInit(const async::ExecInfo& info,
                                      const ReceiverWriterInitParam& param) {
  // TODO

  unsigned int numberOfPoints = info.bufferSize(static_cast<int>(BufferIds::POINTS)) /
                                (sizeof(Eigen::Vector3d));

  logInfo() << "Receiver writer executor received" << numberOfPoints << "points";
  const auto* pointsPtr = static_cast<const Eigen::Vector3d*>(
      info.buffer(static_cast<int>(BufferIds::POINTS)));
  // But then, is the original ptr aligned?
  // Why does this segfault?
  auto points = std::move(std::vector< Eigen::Vector3d>(pointsPtr, pointsPtr + numberOfPoints));
  for (const auto& point : points) {
    std::cout << point << "\n";
  }
  this->points = points;
}

void ReceiverWriterExecutor::exec(const async::ExecInfo& info,
                                  const seissol::writer::ReceiverWriterParam& param) {
  stopwatch.start();

  // TODO(Lukas) Write output

  stopwatch.pause();
}

std::string ReceiverWriter::fileName(unsigned pointId) const {
  std::stringstream fns;
  fns << std::setfill('0') << m_fileNamePrefix << "-receiver-" << std::setw(5) << (pointId+1);
#ifdef PARALLEL
  fns << "-" << std::setw(5) << seissol::MPI::mpi.rank();
#endif
  fns << ".dat";
  return fns.str();
}

void ReceiverWriter::writeHeader( unsigned               pointId,
                                                   Eigen::Vector3d const& point   ) {
  auto name = fileName(pointId);

  std::vector<std::string> names({"xx", "yy", "zz", "xy", "yz", "xz", "v1", "v2", "v3"});
#ifdef USE_POROELASTIC
  std::array<std::string, 4> additionalNames({"p", "v1_f", "v2_f", "v3_f"});
  names.insert(names.end() ,additionalNames.begin(), additionalNames.end());
#endif

  /// \todo Find a nicer solution that is not so hard-coded.
  struct stat fileStat;
  // Write header if file does not exist
  if (stat(name.c_str(), &fileStat) != 0) {
    std::ofstream file;
    file.open(name);
    file << "TITLE = \"Temporal Signal for receiver number " << std::setfill('0') << std::setw(5) << (pointId+1) << "\"" << std::endl;
    file << "VARIABLES = \"Time\"";
#ifdef MULTIPLE_SIMULATIONS
    for (unsigned sim = init::QAtPoint::Start[0]; sim < init::QAtPoint::Stop[0]; ++sim) {
      for (auto const& name : names) {
        file << ",\"" << name << sim << "\"";
      }
    }
#else
    for (auto const& name : names) {
      file << ",\"" << name << "\"";
    }
#endif
    file << std::endl;
    for (int d = 0; d < 3; ++d) {
      file << "# x" << (d+1) << "       " << std::scientific << std::setprecision(12) << point[d] << std::endl;
    }
    file.close();
  }
}

void ReceiverWriter::syncPoint(double time)
{
  write(time);
  if (m_receiverClusters.empty()) {
    return;
  }

  stopwatch.start();

  for (auto& [layer, clusters] : m_receiverClusters) {
    for (auto& cluster : clusters) {
      auto ncols = cluster.ncols();
      for (auto &receiver : cluster) {
        assert(receiver.output.size() % ncols == 0);
        size_t nSamples = receiver.output.size() / ncols;

        std::ofstream file;
        file.open(fileName(receiver.pointId), std::ios::app);
        file << std::scientific << std::setprecision(15);
        for (size_t i = 0; i < nSamples; ++i) {
          for (size_t q = 0; q < ncols; ++q) {
            file << "  " << receiver.output[q + i * ncols];
          }
          file << std::endl;
        }
        file.close();
        receiver.output.clear();
      }
    }
  }

  auto timeToWriteReceiver = stopwatch.stop();
  int const rank = seissol::MPI::mpi.rank();
  logInfo(rank) << "Wrote receivers in" << timeToWriteReceiver << "seconds.";
}
void ReceiverWriter::init(std::string receiverFileName, std::string fileNamePrefix,
                double syncPointInterval, double samplingInterval,
                const MeshReader& mesh,
                const seissol::initializers::Lut& ltsLut,
                const seissol::initializers::LTS& lts,
                const GlobalData* global) {
  m_receiverFileName = std::move(receiverFileName);
  m_fileNamePrefix = std::move(fileNamePrefix);
  m_samplingInterval = samplingInterval;

  addPoints(mesh, ltsLut, lts, global);

  // Add buffers
  int bufferIdPoints = addSyncBuffer(points.data(), points.size() * sizeof(points[0]));
  // Or maybe just remove/add new buffer when size doesn't work

  sendBuffer(bufferIdPoints);

  ReceiverWriterInitParam param{};
  // TODO(Lukas) Add following?
  // param.timestep = seissol::SeisSol::main.checkPointManager().header().value(m_timestepComp);
  callInit(param);

  removeBuffer(bufferIdPoints);

  // Note: Buffer size not nec. constant due to sync points
  // Hence, need to check if size changed (and potentially point)

  setSyncInterval(syncPointInterval);
  // TODO(Lukas) Sim start?
  Modules::registerHook(*this, SYNCHRONIZATION_POINT);
}

void ReceiverWriter::write(double time) {
  // ??
}
void ReceiverWriter::close() {
  wait();
  finalize();
  stopwatch.printTime("Receiver writer frontend:");
  // ?
}
void ReceiverWriter::tearDown() {
  // ?
}


void ReceiverWriter::addPoints(MeshReader const& mesh,
                                                const seissol::initializers::Lut& ltsLut,
                                                const seissol::initializers::LTS& lts,
                                                const GlobalData* global ) {
  const auto rank = seissol::MPI::mpi.rank();
  // Only parse if we have a receiver file
  if (!m_receiverFileName.empty()) {
    points = parseReceiverFile(m_receiverFileName);
    logInfo(rank) << "Record points read from" << m_receiverFileName;
    logInfo(rank) << "Number of record points =" << points.size();
  } else {
    logInfo(rank) << "No record points read.";
  }

  unsigned numberOfPoints = points.size();
  std::vector<short> contained(numberOfPoints);
  std::vector<unsigned> meshIds(numberOfPoints);
  
  // We want to plot all quantities except for the memory variables
  const int n = NUMBER_OF_QUANTITIES - 6*NUMBER_OF_RELAXATION_MECHANISMS;
  std::vector<unsigned> quantities(n);
  std::iota(quantities.begin(), quantities.end(), 0);

  logInfo(rank) << "Finding meshIds for receivers...";
  initializers::findMeshIds(points.data(), mesh, numberOfPoints, contained.data(), meshIds.data());
#ifdef USE_MPI
  logInfo(rank) << "Cleaning possible double occurring receivers for MPI...";
  initializers::cleanDoubles(contained.data(), numberOfPoints);
#endif

  logInfo(rank) << "Mapping receivers to LTS cells...";
  m_receiverClusters[Interior].clear();
  m_receiverClusters[Copy].clear();
  for (unsigned point = 0; point < numberOfPoints; ++point) {
    if (contained[point] == 1) {
      unsigned meshId = meshIds[point];
      unsigned cluster = ltsLut.cluster(meshId);
      LayerType layer = ltsLut.layer(meshId);

      auto& clusters = m_receiverClusters[layer];
      // Make sure that needed empty clusters are initialized.
      for (unsigned c = clusters.size(); c <= cluster; ++c) {
        clusters.emplace_back(global, quantities, m_samplingInterval, syncInterval());
      }

      writeHeader(point, points[point]);
      m_receiverClusters[layer][cluster].addReceiver(meshId, point, points[point], mesh, ltsLut, lts);
    }
  }
}
}
