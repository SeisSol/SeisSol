// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Trace.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <limits>
#include <memory>
#include <mpi.h>
#include <mutex>
#include <ratio>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utils/env.h>
#include <utils/logger.h>
#include <vector>

namespace seissol::monitoring {

namespace {

struct Span {
  std::uint32_t name;
  std::uint32_t track;
  double begin;
  double end;
};

/// the spans one thread has recorded
struct ThreadBuffer {
  std::vector<Span> spans;
};

struct State {
  Trace::Settings settings;
  std::chrono::steady_clock::time_point origin;
  // the window, in microseconds since the origin
  double windowBegin{0};
  double windowEnd{0};
  std::atomic<std::size_t> spans{0};
  std::atomic<bool> truncated{false};

  // the names, the tracks and the buffers of the threads
  std::mutex mutex;
  std::vector<std::string> names;
  std::unordered_map<std::string, std::uint32_t> nameIds;
  std::vector<std::string> tracks;
  std::unordered_map<const void*, std::uint32_t> trackIds;
  std::vector<std::unique_ptr<ThreadBuffer>> buffers;
  std::uint32_t threads{0};
};

State& state() {
  static State state;
  return state;
}

ThreadBuffer& threadBuffer() {
  // the buffer outlives its thread, until the trace gets written
  thread_local ThreadBuffer* buffer = nullptr;
  if (buffer == nullptr) {
    auto& state = monitoring::state();
    const std::scoped_lock lock(state.mutex);
    buffer = state.buffers.emplace_back(std::make_unique<ThreadBuffer>()).get();
  }
  return *buffer;
}

/// a string as a JSON string, with the quotes
std::string jsonString(std::string_view text) {
  std::string result = "\"";
  for (const char c : text) {
    if (c == '"' || c == '\\') {
      result += '\\';
      result += c;
    } else if (static_cast<unsigned char>(c) < 0x20) {
      result += ' ';
    } else {
      result += c;
    }
  }
  result += '"';
  return result;
}

/// writes `count` bytes at `offset`, in pieces that MPI can count
void writeAt(MPI_File file, MPI_Offset offset, const char* data, std::size_t count) {
  constexpr std::size_t Piece = 1UL << 30U;
  for (std::size_t written = 0; written < count; written += Piece) {
    const auto size = std::min(Piece, count - written);
    MPI_File_write_at(file,
                      offset + static_cast<MPI_Offset>(written),
                      data + written,
                      static_cast<int>(size),
                      MPI_BYTE,
                      MPI_STATUS_IGNORE);
  }
}

} // namespace

std::atomic<bool>& Trace::recordingFlag() {
  static std::atomic<bool> recording{false};
  return recording;
}

Trace::Settings Trace::settingsFromEnv(utils::Env& env) {
  Settings settings;
  settings.enabled = env.get<bool>("TRACE", false);
  settings.start = std::max(0.0, env.get<double>("TRACE_START", 0.0));
  settings.duration =
      std::max(0.0, env.get<double>("TRACE_DURATION", std::numeric_limits<double>::infinity()));
  settings.maxSpans = env.get<std::size_t>("TRACE_SPANS", std::size_t{settings.maxSpans});
  return settings;
}

void Trace::start(const Settings& settings, MPI_Comm comm) {
  if (!settings.enabled) {
    return;
  }
  auto& state = monitoring::state();
  state.settings = settings;
  state.spans.store(0, std::memory_order_relaxed);
  state.truncated.store(false, std::memory_order_relaxed);
  state.windowBegin = settings.start * 1e6;
  state.windowEnd = (settings.start + settings.duration) * 1e6;
  if (std::isinf(settings.duration)) {
    logInfo() << "Recording a trace of the time stepping from" << settings.start
              << "s on, of at most" << settings.maxSpans << "spans per rank.";
  } else {
    logInfo() << "Recording a trace of the time stepping from" << settings.start << "s on for"
              << settings.duration << "s, of at most" << settings.maxSpans << "spans per rank.";
  }
  MPI_Barrier(comm);
  state.origin = std::chrono::steady_clock::now();
  recordingFlag().store(true, std::memory_order_release);
}

double Trace::time(std::chrono::steady_clock::time_point point) {
  return std::chrono::duration<double, std::micro>(point - state().origin).count();
}

std::uint32_t Trace::name(std::string_view name) {
  auto& state = monitoring::state();
  const std::scoped_lock lock(state.mutex);
  const auto [entry, added] =
      state.nameIds.try_emplace(std::string(name), static_cast<std::uint32_t>(state.names.size()));
  if (added) {
    state.names.emplace_back(name);
  }
  return entry->second;
}

std::uint32_t Trace::track(const void* key, std::string_view name) {
  auto& state = monitoring::state();
  const std::scoped_lock lock(state.mutex);
  const auto [entry, added] =
      state.trackIds.try_emplace(key, static_cast<std::uint32_t>(state.tracks.size()));
  if (added) {
    state.tracks.emplace_back(name);
  }
  return entry->second;
}

std::uint32_t Trace::threadTrack() {
  thread_local std::uint32_t track = std::numeric_limits<std::uint32_t>::max();
  if (track == std::numeric_limits<std::uint32_t>::max()) {
    auto& state = monitoring::state();
    const std::scoped_lock lock(state.mutex);
    track = static_cast<std::uint32_t>(state.tracks.size());
    state.tracks.emplace_back("host thread " + std::to_string(state.threads++));
  }
  return track;
}

void Trace::span(std::uint32_t name, std::uint32_t track, double begin, double end) {
  auto& state = monitoring::state();
  if (!recording() || end < state.windowBegin || begin > state.windowEnd) {
    return;
  }
  if (state.spans.fetch_add(1, std::memory_order_relaxed) >= state.settings.maxSpans) {
    state.truncated.store(true, std::memory_order_relaxed);
    return;
  }
  threadBuffer().spans.push_back({name, track, begin, end});
}

void Trace::write(const std::string& path, MPI_Comm comm) {
  auto& state = monitoring::state();
  if (!state.settings.enabled) {
    return;
  }
  recordingFlag().store(false, std::memory_order_release);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  // all elements but the first one start with a comma; the first one comes from rank 0
  std::ostringstream text;
  text << std::fixed << std::setprecision(3);
  const auto processName = [&](int process) {
    return R"({"name":"process_name","ph":"M","pid":)" + std::to_string(process) +
           R"(,"args":{"name":"rank )" + std::to_string(process) + R"("}})";
  };
  const std::string header = R"({"traceEvents":[)"
                             "\n" +
                             processName(0);
  const std::string footer = "\n"
                             R"(],"displayTimeUnit":"ms"})"
                             "\n";
  std::size_t spans = 0;
  {
    const std::scoped_lock lock(state.mutex);
    if (rank != 0) {
      text << ",\n" << processName(rank);
    }
    text << ",\n"
         << R"({"name":"process_sort_index","ph":"M","pid":)" << rank << R"(,"args":{"sort_index":)"
         << rank << "}}";
    for (std::size_t track = 0; track < state.tracks.size(); ++track) {
      text << ",\n"
           << R"({"name":"thread_name","ph":"M","pid":)" << rank << R"(,"tid":)" << track
           << R"(,"args":{"name":)" << jsonString(state.tracks[track]) << "}}";
      text << ",\n"
           << R"({"name":"thread_sort_index","ph":"M","pid":)" << rank << R"(,"tid":)" << track
           << R"(,"args":{"sort_index":)" << track << "}}";
    }
    for (const auto& buffer : state.buffers) {
      for (const auto& span : buffer->spans) {
        // (a span that began before the window shows from its start)
        const auto begin = std::max(span.begin, state.windowBegin);
        text << ",\n"
             << R"({"name":)" << jsonString(state.names[span.name]) << R"(,"ph":"X","pid":)" << rank
             << R"(,"tid":)" << span.track << R"(,"ts":)" << begin << R"(,"dur":)"
             << std::max(0.0, span.end - begin) << "}";
      }
      spans += buffer->spans.size();
      buffer->spans.clear();
      buffer->spans.shrink_to_fit();
    }
  }
  const auto chunk = text.str();

  // the chunks of the processes follow each other, in the order of their ranks
  unsigned long long size = chunk.size();
  unsigned long long offset = 0;
  MPI_Exscan(&size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
  if (rank == 0) {
    offset = 0;
  }
  unsigned long long total = 0;
  MPI_Allreduce(&size, &total, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);

  MPI_File file{};
  MPI_File_open(comm, path.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
  MPI_File_set_size(file, 0);
  if (rank == 0) {
    writeAt(file, 0, header.data(), header.size());
    writeAt(file, static_cast<MPI_Offset>(header.size() + total), footer.data(), footer.size());
  }
  writeAt(file, static_cast<MPI_Offset>(header.size() + offset), chunk.data(), chunk.size());
  MPI_File_close(&file);

  std::array<unsigned long long, 2> counts{spans, state.truncated.load() ? 1ULL : 0ULL};
  MPI_Allreduce(MPI_IN_PLACE, counts.data(), counts.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
  logInfo() << "Wrote the trace of" << counts[0] << "spans to" << path.c_str();
  if (counts[1] > 0) {
    logWarning() << "On" << counts[1] << "ranks, the trace reached its limit of"
                 << state.settings.maxSpans
                 << "spans (SEISSOL_TRACE_SPANS) and recorded no more after that.";
  }
}

} // namespace seissol::monitoring
