cmake_minimum_required(VERSION 3.10)
project(parse_BLAST_result)

set(CMAKE_CXX_STANDARD 14)

add_executable(parse_BLAST_result main.cpp object_partition.cpp object_partition.h)