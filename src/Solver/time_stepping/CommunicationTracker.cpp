#include "CommunicationTracker.h"

    void track_send(int send_size)
    {
        assert(send_size > 0);
        messages_sent += 1;
        total_bytes_sent += send_size;
    }

    void track_recv(int recv_size)
    {
        assert(recv_size > 0);
        messages_received += 1;
        total_bytes_received += recv_size;
    }

    void track_test()
    {
        test_call++;
    }

    std::string communication_stats()
    {
        int rank = 0;
        
        #ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        #endif

        std::stringstream ss;
        ss << "Rank " << rank << " statistics of Isend and Irecv calls\n";
        ss << "Rank " << rank << " has sent: " << messages_sent.load() << " messages\n";
        ss << "Rank " << rank << " has sent: " << total_bytes_sent.load() << " bytes\n";
        ss << "Rank " << rank << " sent message size is on average: " << static_cast<double>(total_bytes_sent.load()) / static_cast<double>(messages_sent.load()) << "\n";
        ss << "Rank " << rank << " has received: " << messages_received.load() << " messages\n";
        ss << "Rank " << rank << " has received: " << total_bytes_received.load() << " bytes\n";
        ss << "Rank " << rank << " received message size is on average: " << static_cast<double>(total_bytes_received.load()) / static_cast<double>(messages_received.load()) << "\n";
        ss << "Rank " << rank << " called test: " << test_call.load();
        return ss.str();
    }

