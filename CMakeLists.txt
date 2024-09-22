cmake_minimum_required(VERSION 3.16)
project(compileFromCpp)

set(CMAKE_CXX_STANDARD 14)

add_executable(compileFromCpp
        double_slit2_FFT_tri.c
        TestEFunction.c)
