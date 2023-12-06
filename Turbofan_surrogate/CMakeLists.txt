cmake_minimum_required(VERSION 3.26)
project(Turbofan_surrogate)

set(CMAKE_CXX_STANDARD 14)

add_executable(Turbofan_surrogate main.cpp
        Design.h
        Off_Design.h)
