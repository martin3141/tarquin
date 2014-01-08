#ifndef __MRI__
#define __MRI__

#include <iostream>
#include <QtGui>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SVD>

struct Triplet 
{
	double x;
	double y;
	double z;
	
	Triplet() : 
        x(0), y(0), z(0) { }

	Triplet(double x_, double y_, double z_) :
		x(x_), y(y_), z(z_) { }
};

class MRI
{
    private:

        size_t rows;
        size_t cols;
        size_t slices;

        size_t ba;
        size_t bs;

        double row_dim;
        double col_dim;
        double slice_dim;

        double ww; // Window width
        double wc; // Window center

        std::vector<Eigen::Vector3d> row_dirn;
        std::vector<Eigen::Vector3d> col_dirn;
        std::vector<Eigen::Vector3d> slice_dirn;
        std::vector<Eigen::Vector3d> pos;

        std::vector<QImage> image_data;

        std::vector<std::vector<std::vector<double> > > pixel_data;

        double max_val;
        double min_val;

        double rescale_slope;
        double rescale_int;

        std::string series_uid;

    public:

        MRI();

        ~MRI();

        bool load_from_file(std::string strFilename);
        
        bool match_uid(std::string strFilename);

        std::vector<std::vector<std::vector<double> > > get_pixel_data()
        {
            return pixel_data;
        }

         // load_from_file_list(std::vector<std::string>);

        void set_ww(double new_ww)
        {
            if ( new_ww < 1 )
                new_ww = 1;

            this->ww = new_ww;
        }

        double get_ww()
        {
            return ww;
        }
        
        double get_rows()
        {
            return rows; 
        }
        
        double get_cols()
        {
            return cols; 
        }

        void set_wc(double new_wc)
        {
            this->wc = new_wc;
        }
        
        double get_wc()
        {
            return wc;
        }

        double get_max_val()
        {
            return this->max_val;
        }

        double get_min_val()
        {
            return this->min_val;
        }

        double get_row_fov()
        {
            return rows * row_dim;
        }

        double get_col_fov()
        {
            return cols * col_dim;
        }

        void generate_slices();
        
        void reorder_slices();
        
        size_t get_closest_slice(Eigen::Vector3d pos);
        
        void generate_slice(size_t slice);

        QImage& get_slice(size_t slice);
        
        void append_slices(MRI& image);
        
        size_t get_num_slices()
        {
            return slices;
        }

        double get_row_dim()
        {
            return row_dim;
        }

        double get_col_dim()
        {
            return col_dim;
        }

        double get_slice_dim()
        {
            return slice_dim;
        }

        std::string get_uid()
        {
            return series_uid;
        }
        
        Eigen::Vector3d get_row_dirn(size_t n)
        {
            return row_dirn[n-1];
        }

        Eigen::Vector3d get_col_dirn(size_t n)
        {
            return col_dirn[n-1];
        }
        
        Eigen::Vector3d get_pos(size_t n)
        {
            return pos[n-1];
        }

        std::vector<Eigen::Vector3d> get_row_dirn()
        {
            return row_dirn;
        }

        std::vector<Eigen::Vector3d> get_col_dirn()
        {
            return col_dirn;
        }

        std::vector<Eigen::Vector3d> get_slice_dirn()
        {
            return slice_dirn;
        }
        
        std::vector<Eigen::Vector3d> get_pos()
        {
            return pos;
        }

        void print_paras()
        {
            std::cout << std::endl;
            std::cout << "MRI details" << std::endl;
            std::cout << "-----------" << std::endl;
            std::cout << "Rows      : " << this->rows << std::endl;
            std::cout << "Cols      : " << this->cols << std::endl;
            std::cout << "Slices    : " << this->slices << std::endl;
            std::cout << "Row dim   : " << this->row_dim << std::endl;
            std::cout << "Col dim   : " << this->col_dim << std::endl;
            std::cout << "Slice dim : " << this->slice_dim << std::endl;
            std::cout << "Bits all. : " << this->ba << std::endl;
            std::cout << "Bits sto. : " << this->bs << std::endl;
        }
};

#endif
