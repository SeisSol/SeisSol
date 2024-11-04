
include(sources.cmake)

add_subdirectory(IO)
target_link_libraries(seissol-lib PUBLIC seissol-io)
