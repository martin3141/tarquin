#ifndef __protobuf__
#define __protobuf__

#include "tarquin.hpp"
#include "tarquin_pb.hpp"
#include "CFID.hpp"
#include "CBasis.hpp"

// save results to protobuf
bool save_results(std::string file_name, tarquin::Workspace results);

// load results to protobuf
bool load_results(std::string file_name, tarquin::Workspace& results);

// save basis to protobuf
bool save_basis(std::string file_name, tarquin::CBasis basis);

// load basis from protobuf
bool load_basis(std::string file_name, tarquin::CBasis& basis);

// save fid to protobuf file
bool save_fid(std::string file_name, tarquin::CFID fid);

// load fid from protobuf file
bool load_fid(std::string file_name, tarquin::CFID& fid);

// save fid to protobuf object 
bool serialise_fid(sln::fid* sln_fid, tarquin::CSignal::fid_iterator& fid);
bool serialise_fid(sln::fid* sln_fid, tarquin::CFID fid);

// save basis to protobuf object 
bool serialise_basis(sln::basis* sln_basis, tarquin::CBasis& basis);

// save options to protobuf object 
bool serialise_options(sln::options* sln_options, tarquin::Options& options);

// copy protobuf options object to CObject
void copy_options(sln::options sln_options, tarquin::Options& options);

// copy protobuf fid object to CFID
void copy_fid(sln::fid sln_fid, tarquin::CFID& fid);

// copy protobuf basis object to CBasis
void copy_basis(sln::basis basis, tarquin::CBasis& tbasis); 

#endif
