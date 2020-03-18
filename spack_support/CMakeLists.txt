cmake_minimum_required(VERSION 3.5)
project(seissol-spack-support)

install(DIRECTORY packages 
        DESTINATION spack_support)

install(FILES repo.yaml 
        DESTINATION spack_support)

