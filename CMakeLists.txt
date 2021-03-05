cmake_minimum_required(VERSION 3.10)
project(parse_BLAST_result)
include_directories(./include)
aux_source_directory(./src DIR_SRCS)
set(CMAKE_CXX_STANDARD 11)

add_executable(parse_BLAST_result main.cpp src/object_partition.cpp include/object_partition.h src/global.cpp include/global.h)