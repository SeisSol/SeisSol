#ifndef TRACKER_INCLUDED
#define TRACKER_INCLUDED

#include <atomic>
#include <stdint.h>
#include <mpi.h>
#include <string>
#include <sstream>
#include <cassert>

static std::atomic<int64_t> messages_sent = 0;
static std::atomic<int64_t> messages_received = 0;
static std::atomic<int64_t> test_call = 0;
static std::atomic<int64_t> total_bytes_sent = 0;
static std::atomic<int64_t> total_bytes_received = 0;

void track_send(int send_size);

void track_recv(int recv_size);

void track_test();

std::string communication_stats();

#endif