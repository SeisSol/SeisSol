#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    fs::path p{argv[0]};
    std::cout << fs::canonical(p) << std::endl;
    return 0;
}
