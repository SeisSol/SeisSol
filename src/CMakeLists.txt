
include(sources.cmake)

add_subdirectory(IO)
target_link_libraries(SeisSol-lib PUBLIC seissol-io)
