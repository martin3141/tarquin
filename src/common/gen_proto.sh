#! /bin/bash

# A handy script for regenerating the tarquin_pb.cpp and hpp files
# from the tarquin.proto description, MW 09 Mar 11.

protoc tarquin.proto --cpp_out=./
sed 's/tarquin.pb.h/tarquin_pb.hpp/g' tarquin.pb.cc > tarquin_pb.cpp
rm tarquin.pb.cc
mv tarquin.pb.h tarquin_pb.hpp
