#include "MRI.hpp"
#include "fidio/CDICOMFile.hpp"
#include <stdint.h>
#include "common.hpp" 
#include <sstream>


typedef std::pair<double,size_t> double_ind_pair;
bool comparator ( const double_ind_pair& l, const double_ind_pair& r)
    { return l.first < r.first; }

MRI::MRI()
{

}

MRI::~MRI()
{

}

size_t MRI::get_closest_slice(Eigen::Vector3d pos)
{
    
    double min_dist = std::numeric_limits<double>::infinity();
    size_t closest_slice = 0;

    for ( size_t slice = 0; slice < this->slices; slice++ )
    {
        double a = this->slice_dirn[slice](0);
        double b = this->slice_dirn[slice](1);
        double c = this->slice_dirn[slice](2);
        double d = - this->pos[slice](0)*a - this->pos[slice](1)*b - this->pos[slice](2)*c;
        double dist = abs(a*pos(0)+b*pos(1)+c*pos(2)+d);
        if ( dist < min_dist )
        {
            min_dist = dist;
            closest_slice = slice;
        }
    }

    return closest_slice+1;
}

void MRI::generate_slices()
{
    this->image_data.clear();

    QVector<QRgb> greyscale;
    for(int i = 0; i < 256; i++) 
        greyscale.push_back(QColor(i,i,i).rgb());

    if ( ( this->wc == -1 ) && ( this-> ww == -1 ) )
    {
        this->wc = ((rescale_slope*this->max_val+rescale_int) - (rescale_slope*this->min_val+rescale_int))/2.0;
        this->ww = ((rescale_slope*this->max_val+rescale_int) - (rescale_slope*this->min_val+rescale_int));
    }

    double ymin = 0;
    double ymax = 255;

    for ( size_t slice = 0; slice < this->slices; slice++ )
    {
        QImage img = QImage(this->cols, this->rows, QImage::Format_Indexed8);
        img.setColorTable(greyscale);
        for ( size_t row = 0 ; row < this->rows; row++ )
            for ( size_t col = 0; col < this->cols; col++ )
            {
                double pixel; 
                pixel = this->pixel_data[slice][row][col];

                pixel = rescale_slope*pixel + rescale_int;
                if ( pixel <= (this->wc - 0.5 - (this->ww-1)/2) ) 
                    pixel = ymin;
                else if ( pixel > (this->wc - 0.5 + (this->ww-1)/2) ) 
                    pixel = ymax;
                else
                    pixel = ((pixel - (this->wc - 0.5))/(this->ww-1) + 0.5) * (ymax-ymin) + ymin;

                img.setPixel(col, row, pixel);
            }
        this->image_data.push_back(img);
    }
}

void MRI::generate_slice(size_t gen_slice)
{
    QVector<QRgb> greyscale;
    for(int i = 0; i < 256; i++) 
        greyscale.push_back(QColor(i,i,i).rgb());

    if ( ( this->wc == -1 ) && ( this-> ww == -1 ) )
    {
        this->wc = ((rescale_slope*this->max_val+rescale_int) - (rescale_slope*this->min_val+rescale_int))/2.0;
        this->ww = ((rescale_slope*this->max_val+rescale_int) - (rescale_slope*this->min_val+rescale_int));
    }

    double ymin = 0;
    double ymax = 255;

    QImage img = QImage(this->cols, this->rows, QImage::Format_Indexed8);
    img.setColorTable(greyscale);
    for ( size_t row = 0 ; row < this->rows; row++ )
        for ( size_t col = 0; col < this->cols; col++ )
        {
            double pixel; 
            pixel = this->pixel_data[gen_slice-1][row][col];

            pixel = rescale_slope*pixel + rescale_int;
            if ( pixel <= (this->wc - 0.5 - (this->ww-1)/2) ) 
                pixel = ymin;
            else if ( pixel > (this->wc - 0.5 + (this->ww-1)/2) ) 
                pixel = ymax;
            else
                pixel = ((pixel - (this->wc - 0.5))/(this->ww-1) + 0.5) * (ymax-ymin) + ymin;

            img.setPixel(col, row, pixel);
        }
    this->image_data[gen_slice-1] = img;
}


QImage& MRI::get_slice(size_t slice)
{
    // TODO checks on slice value
    return this->image_data[slice-1];
}

bool MRI::match_uid(std::string strFilename)
{
    tarquin::CDICOMFile file;
    file.Open(strFilename);
    
    long series_id_bytes = file.MoveToTag("0020","000E");
    if( -1 == series_id_bytes )
    {
        return false;
    }
    else
    {
        std::vector<char> series_id(series_id_bytes, 0);
        file.GetFileStream().read(&series_id[0], series_id_bytes);
        std::string series_id_str(series_id.begin(), series_id.end());
        if ( this->series_uid == series_id_str.substr(0,series_id_bytes) )
            return true;
        else
            return false;
    }
}


bool MRI::load_from_file(std::string strFilename)
{
    tarquin::CDICOMFile file;
    file.Open(strFilename);

    // find rows tag
    long bytes_rows = file.MoveToTag("0028", "0010"); 
    if( -1 != bytes_rows )
    {
        uint16_t img_rows = 0;
        file.GetFileStream().read((char*)&img_rows, 2);
        this->rows = static_cast<size_t>(img_rows);
    }
    else
        return false;

    // find cols tag
    long bytes_cols = file.MoveToTag("0028", "0011"); 
    if( -1 != bytes_cols )
    {
        uint16_t img_cols = 0;
        file.GetFileStream().read((char*)&img_cols, 2);
        this->cols = static_cast<size_t>(img_cols);
    }
    else
        return false;

    // find frames tag
    long frame_bytes = file.MoveToTag("0028", "0008"); 
    if( -1 != frame_bytes )
    {
        std::vector<char> frames_char(frame_bytes, 0);
        file.GetFileStream().read(&frames_char[0], frame_bytes);
        std::string frames_str(frames_char.begin(), frames_char.end());
        std::istringstream iss(frames_str, std::istringstream::in);
        int frames = 0;
        iss >> frames;
        this->slices = static_cast<size_t>(frames);
    }
    else
    {
        this->slices = 1;
        file.Close();
	    file.Open(strFilename);
    }

    //return false;

    // SliceThickness
    long thick_bytes = file.MoveToTag("0018","0050");
    if( -1 == thick_bytes )
        return false;

    std::vector<char> thick(thick_bytes, 0);
	file.GetFileStream().read(&thick[0], thick_bytes);
    std::string thick_str(thick.begin(), thick.end());
    //std::cout << "Thickness : " << thick_str.substr(0,thick_bytes) << std::endl;

    // convert slicethickness to double
    std::istringstream istr_slice(thick_str.substr(0,thick_bytes));
    double slice_thick = 0.0;
    istr_slice >> slice_thick;
    this->slice_dim = slice_thick;

    // find pixel spacing tag
    long dims_bytes = file.MoveToTag("0028","0030");
    if( -1 == dims_bytes )
        return false;

    std::vector<char> dims(dims_bytes, 0);
	file.GetFileStream().read(&dims[0], dims_bytes);
    std::string dims_str(dims.begin(), dims.end());
    //std::cout << "Pixel Spacing : " << dims_str.substr(0,dims_bytes) << std::endl;
    // convert pixel spacing to doubles
    std::string str = dims_str.substr(0,dims_bytes);
    std::vector<double> row_col;
    tarquin::str2rvec(str, row_col);
    this->row_dim = row_col[0];
    this->col_dim = row_col[1];

    // find position tags
    long pos_bytes = file.MoveToTag("0020","0032");
    while ( -1 != pos_bytes )
    {
        //std::cout << pos_bytes << std::endl << std::flush;
        std::vector<char> pos(pos_bytes, 0);
        file.GetFileStream().read(&pos[0], pos_bytes);
        std::string pos_str(pos.begin(), pos.end());
        std::string pos_str_cut = pos_str.substr(0,pos_bytes);
        std::vector<double> pos_vec;
        tarquin::str2rvec(pos_str_cut, pos_vec);
        Eigen::Array3d temp_pos;
        for ( size_t n = 0; n < 3; n++ )
        {
            temp_pos[n] = pos_vec[n];
        }

        this->pos.push_back(temp_pos);
        //std::cout << "Position : " << pos_str_cut << std::endl << std::flush;
        pos_bytes = file.MoveToTag("0020","0032",false);
    }
    
    // this may be required, probally because the reader gets past the eof
	file.Close();
	file.Open(strFilename);
    
    long ori_bytes = file.MoveToTag("0020","0037");
    while ( -1 != ori_bytes )
    {
        std::vector<char> ori(ori_bytes, 0);
        file.GetFileStream().read(&ori[0], ori_bytes);
        std::string ori_str(ori.begin(), ori.end());
        std::string ori_str_cut = ori_str.substr(0,ori_bytes);
        //std::cout << "Orientation : " << ori_str_cut << std::endl << std::flush;

        std::vector<double> ori_vec;
        tarquin::str2rvec(ori_str_cut, ori_vec);

        Eigen::Vector3d temp_row_ori;
        for ( size_t n = 0; n < 3; n++ )
            temp_row_ori[n] = ori_vec[n];

        this->row_dirn.push_back(temp_row_ori);

        Eigen::Vector3d temp_col_ori;
        for ( size_t n = 3; n < 6; n++ )
            temp_col_ori[n-3] = ori_vec[n];

        this->col_dirn.push_back(temp_col_ori);

        Eigen::Vector3d temp_slice_ori = temp_row_ori.cross(temp_col_ori);
        temp_slice_ori.normalize();

        this->slice_dirn.push_back(temp_slice_ori);

        ori_bytes = file.MoveToTag("0020","0037",false);
    }
    
    // this may be required, probally because the reader gets past the eof
	file.Close();
	file.Open(strFilename);

    // find data tag
    long bytes = file.MoveToTag("7FE0","0010"); 
    if( -1 == bytes ) 
    {
        std::cout << "Rows   : " << this->rows << std::endl;
        std::cout << "Cols   : " << this->cols << std::endl;
        std::cout << "Slices : " << this->slices << std::endl;
        std::cout << "MRI pixel data is : " << bytes << " bytes." << std::endl;
        std::cout << "MRI file did not contain the magic tag (7FE0, 0010)." << std::endl;
        return false;
    }
    
    // if it's not the right size hunt for another, Siemens sometimes have 2 data tags
    if ( bytes != 2*this->slices*this->rows*this->cols )
    {
        file.GetFileStream().seekg(bytes, std::ios_base::cur);

        bytes = file.MoveToTag("7FE0","0010",false); 

        if( -1 == bytes || bytes != 2*this->slices*this->rows*this->cols ) 
        {
            std::cout << "Rows   : " << this->rows << std::endl;
            std::cout << "Cols   : " << this->cols << std::endl;
            std::cout << "Slices : " << this->slices << std::endl;
            std::cout << "MRI pixel data is : " << bytes << " bytes." << std::endl;
            std::cout << "MRI pixel data was found but not the right size." << std::endl;
            return false;
        }
    }

    //std::cout << "MRI pixel data is : " << bytes << " bytes." << std::endl;

    // read into pixel_data
    this->pixel_data.resize(this->slices);
    for (int i = 0; i < this->slices; ++i)
    {
        pixel_data[i].resize(this->rows);
        for (int j = 0; j < this->rows; ++j)
            pixel_data[i][j].resize(this->cols);
    }
   
    this->max_val = -std::numeric_limits<double>::infinity();
    this->min_val = std::numeric_limits<double>::infinity();
    for ( size_t slice = 0; slice < this->slices; slice++ )
        for ( size_t row = 0 ; row < this->rows; row++ )
            for ( size_t col = 0; col < this->cols; col++ )
            {
                //unsigned short int pixel; 
                //file.GetFileStream().read((char*)&pixel, 2);

                unsigned short int pixel;  // 2 bytes unsigned
                file.GetFileStream().read((char*)&pixel, 2);
                
                // bit mask up to 12 bits
                pixel = pixel & ((1 << 12)-1);

                this->pixel_data[slice][row][col] = static_cast<double>(pixel);
                if ( static_cast<double>(pixel) > this->max_val )
                    this->max_val = static_cast<double>(pixel);

                if ( static_cast<double>(pixel) < this->min_val )
                    this->min_val = static_cast<double>(pixel);
            }


    // find rescale intercept tag
    long int_bytes = file.MoveToTag("0028","1052");
    if( -1 == int_bytes )
    {
        this->rescale_int = 0;
    }
    else
    {
        std::vector<char> rescale_int(int_bytes, 0);
        file.GetFileStream().read(&rescale_int[0], int_bytes);
        std::string rescale_int_str(rescale_int.begin(), rescale_int.end());
        std::string str = rescale_int_str.substr(0,int_bytes);
        std::stringstream str_stream(str);
        str_stream >> this->rescale_int;
    }
    
    // find rescale slope tag
    long slope_bytes = file.MoveToTag("0028","1053");
    if( -1 == slope_bytes )
    {
        this->rescale_slope = 1;
    }
    else
    {
        std::vector<char> rescale_slope(slope_bytes, 0);
        file.GetFileStream().read(&rescale_slope[0], slope_bytes);
        std::string rescale_slope_str(rescale_slope.begin(), rescale_slope.end());
        std::string str = rescale_slope_str.substr(0,slope_bytes);
        std::stringstream str_stream(str);
        str_stream >> this->rescale_slope;
    }
    
    file.Close();
	file.Open(strFilename);

    long series_id_bytes = file.MoveToTag("0020","000E");
    if( -1 == series_id_bytes )
    {
        std::cout << "Warning series UID not found" << std::endl;
    }
    else
    {
        std::vector<char> series_id(series_id_bytes, 0);
        file.GetFileStream().read(&series_id[0], series_id_bytes);
        std::string series_id_str(series_id.begin(), series_id.end());
        this->series_uid = series_id_str.substr(0,series_id_bytes);
    }


    /*long ww_bytes = file.MoveToTag("0028","1051");
    if( -1 == ww_bytes )
    {
        this->ww = -1;
    }
    else
    {
        std::vector<char> ww_char(ww_bytes, 0);
        file.GetFileStream().read(&ww_char[0], ww_bytes);
        std::string ww_str(ww_char.begin(), ww_char.end());
        std::string str = ww_str.substr(0,ww_bytes);
        std::stringstream ww_stream(str);
        ww_stream >> this->ww;
    }
    
    long wc_bytes = file.MoveToTag("0028","1050");
    if( -1 == wc_bytes )
    {
        this->wc = -1;
    }
    else
    {
        std::vector<char> wc_char(wc_bytes, 0);
        file.GetFileStream().read(&wc_char[0], wc_bytes);
        std::string wc_str(wc_char.begin(), wc_char.end());
        std::string str = wc_str.substr(0,wc_bytes);
        std::stringstream wc_stream(str);
        wc_stream >> this->wc;
    }
    */

    this->wc = -1;
    this->ww = -1;

    // do some final checks
    if ( this->slices != this->pos.size() )
        return false;

    if ( this->slices != this->row_dirn.size() )
        return false;

    if ( this->slices != this->col_dirn.size() )
        return false;

    this->generate_slices();

    return true;
}

void MRI::append_slices(MRI& image)
{
    // image check rows and cols match
    if ( ( image.get_rows() != this->rows ) || ( image.get_cols() != this->cols ) )
        return;

    this->slices = this->slices + image.get_num_slices();
    
    // append pixel data
    std::vector<std::vector<std::vector<double> > > temp_pixel_data = image.get_pixel_data();
    this->pixel_data.insert(this->pixel_data.end(), temp_pixel_data.begin(), temp_pixel_data.end() );

    std::vector<Eigen::Vector3d> temp_row_dirn = image.get_row_dirn();
    this->row_dirn.insert(this->row_dirn.end(), temp_row_dirn.begin(), temp_row_dirn.end());

    std::vector<Eigen::Vector3d> temp_col_dirn = image.get_col_dirn();
    this->col_dirn.insert(this->col_dirn.end(), temp_col_dirn.begin(), temp_col_dirn.end());
    
    std::vector<Eigen::Vector3d> temp_slice_dirn = image.get_slice_dirn();
    this->slice_dirn.insert(this->slice_dirn.end(), temp_slice_dirn.begin(), temp_slice_dirn.end());

    std::vector<Eigen::Vector3d> temp_pos = image.get_pos();
    this->pos.insert(this->pos.end(), temp_pos.begin(), temp_pos.end());

    // adjust min_val and max_val
    if ( image.get_max_val() > this->max_val )
        this->max_val = image.get_max_val();

    if ( image.get_min_val() < this->min_val )
        this->min_val = image.get_min_val();
    
}

void MRI::reorder_slices()
{
    Eigen::Vector3d slice_norm = this->slice_dirn[0];
    
    std::vector<double_ind_pair> dist;

    for ( size_t n = 0; n < this->slices; n++ )
    {
        double_ind_pair temp_dist;
        temp_dist = std::make_pair(slice_norm.dot(pos[n]), n);
        dist.push_back(temp_dist);
        //std::cout << slice_norm.dot(pos[n]) << std::endl;
    }
    
    // sort to get an index
    std::sort( dist.begin(), dist.end() );

    std::vector<std::vector<std::vector<double> > > new_pixel_data;

    std::vector<Eigen::Vector3d> new_row_dirn;
    std::vector<Eigen::Vector3d> new_col_dirn;
    std::vector<Eigen::Vector3d> new_slice_dirn;
    std::vector<Eigen::Vector3d> new_pos;

    for ( size_t n = 0; n < this->slices; n++ )
    {
        //std::cout << dist[n].second << std::endl;
        new_pixel_data.push_back(this->pixel_data[dist[n].second]);
        new_row_dirn.push_back(this->row_dirn[dist[n].second]);
        new_col_dirn.push_back(this->col_dirn[dist[n].second]);
        new_slice_dirn.push_back(this->slice_dirn[dist[n].second]);
        new_pos.push_back(this->pos[dist[n].second]);
    }

    this->pixel_data = new_pixel_data;
    this->row_dirn = new_row_dirn;
    this->col_dirn = new_col_dirn;
    this->slice_dirn = new_slice_dirn;
    this->pos = new_pos;
}
