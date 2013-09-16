#include "test_util.hpp"
#include <boost/test/unit_test.hpp>

cvm::cvector read_ref_data(std::string filename)
{
	size_t file_size = tarquin::GetFileSize(filename);
	
	// make sure that we could find the file we are using as the reference data
	BOOST_CHECK_GT(file_size, 0);

	// how many complex numbers are there in the file?
	size_t num_cplx = file_size / 8;

	// now read them in
	std::ifstream fin(filename.c_str(), std::ios_base::binary);

	cvm::cvector out(num_cplx);
	
	for( size_t i = 0; i < num_cplx; ++i )
	{
		float x = 0;
		float y = 0;
		fin.read((char*)&x, 4);
		fin.read((char*)&y, 4);

		out[i+1] = tarquin::tcomplex(x, y);
	}

    //tarquin::plot(out);

	return out;
}

void write_ref_data(cvm::cvector cvec, std::string filename)
{

    //tarquin::plot(cvec);

	std::ofstream fout(filename.c_str(), std::ios_base::binary);

	if( fout.bad() )
		throw std::runtime_error("failed to open output file for writing");

	for( int n = 0; n < cvec.size(); ++n )
	{
		float x = real(cvec[n+1]);
		float y = imag(cvec[n+1]);

		fout.write((char*)&x, 4);
		fout.write((char*)&y, 4);
	}
}

void write_ref_data_eigen(Eigen::VectorXcd cvec, std::string filename)
{
	std::ofstream fout(filename.c_str(), std::ios_base::binary);

	if( fout.bad() )
		throw std::runtime_error("failed to open output file for writing");

	for( int n = 0; n < cvec.size(); ++n )
	{
		float x = real(cvec[n]);
		float y = imag(cvec[n]);

		fout.write((char*)&x, 4);
		fout.write((char*)&y, 4);
	}
}

Eigen::VectorXcd read_ref_data_eigen(std::string filename)
{
	size_t file_size = tarquin::GetFileSize(filename);
	
	// make sure that we could find the file we are using as the reference data
	BOOST_CHECK_GT(file_size, 0);

	// how many complex numbers are there in the file?
	size_t num_cplx = file_size / 8;

	// now read them in
	std::ifstream fin(filename.c_str(), std::ios_base::binary);

	Eigen::VectorXcd out(num_cplx);
	
	for( size_t i = 0; i < num_cplx; ++i )
	{
		float x = 0;
		float y = 0;
		fin.read((char*)&x, 4);
		fin.read((char*)&y, 4);

		out[i] = tarquin::tcomplex(x, y);
	}

	return out;
}
