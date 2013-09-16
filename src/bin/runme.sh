
#./bin/tarquin --input /media/disk/PELMI/WS/RAW --format lcm --input_w /media/disk/PELMI/W/RAW --fs 2000 --ft 63e6 --basis_xml ../basis/phd2000.xml --show_pre true 

#./bin/tarquin --input /media/disk/HASLU/WS/RAW --format lcm --input_w /media/disk/HASLU/W/RAW --fs 2000 --ft 63e6 --basis_xml ../basis/phd2000.xml --show_pre true 

#./tarquin --input /media/disk/HASLU/WS/RAW --format lcm --input_w /media/disk/HASLU/W/RAW --fs 1000 --ft 63e6 --basis_csv ../basis/phd --output_basis ../basis/phd1000.xml

#./tarquin --input /media/disk/HASLU/WS/RAW --format lcm --input_w /media/disk/HASLU/W/RAW --fs 1000 --ft 63e6 --basis_xml ../basis/phd1000.xml --start_pnt 3 --end_pnt 512 --water_eddy true --show_pre true

#./tarquin --input ../../data/LCM/WS/RAW --format lcm --input_w ../../data/LCM/W/RAW --echo 0.030 --fs 1000 --ft 63e6 --basis_xml ../../basis/LCM.xml --output_txt results.txt --output_csv results.csv --output_image plot.png

#./tarquin --input ../../data/LCM/WS/RAW --format lcm --input_w ../../data/LCM/W/RAW --echo 0.030 --fs 1000 --ft 63.3684e6 --basis_csv ../../basis/LCModel --output_basis ../../basis/LCM.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false 

#./tarquin --input ../../data/LCM/WS/RAW --format lcm --input_w ../../data/LCM/W/RAW --echo 0.030 --fs 1000 --ft 63.3684e6 --basis_xml ../../basis/LCM.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy true --show_pre true

#./tarquin --input ../../data/real_mb.dpt --format dpt --echo 0.030 --fs 1000 --ft 63.3684e6 --basis_xml ../../basis/phd1000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --max_iters 1

./tarquin --input ../../data/real_mb.dpt --format dpt --echo 0.030 --fs 1000 --ft 63.3684e6 --basis_xml ../../basis/phd1000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --ref 4.7 --max_iters 15

#./tarquin --input /media/disk/MARTINLCM2/WS/RAW --format lcm --input_w /media/disk/MARTINLCM2/W/RAW --echo 0.030 --fs 2000 --ft 63.3684e6 --basis_xml ../../basis/LCM2000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false

#./tarquin --input /media/disk/MARTINLCM2/WS/RAW --format lcm --input_w /media/disk/MARTINLCM2/W/RAW --echo 0.030 --fs 2000 --ft 63.3684e6 --basis_csv ../../basis/LCModel --output_basis ../../basis/LCM2000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false

#./tarquin --input /media/disk/MARTINLCM2/WS/RAW --format lcm --input_w /media/disk/MARTINLCM2/W/RAW --echo 0.030 --fs 2000 --ft 63.3684e6 --basis_csv ../../basis/LCModel --output_basis ../../basis/LCM2000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false

#./tarquin --input /media/disk/MARTINLCM2/WS/RAW --format lcm --input_w /media/disk/MARTINLCM2/W/RAW --echo 0.03 --fs 2000 --ft 63.3684e6 --basis_xml ../../basis/LCM2000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false

#./tarquin --input /media/disk/MARTINLCM2/WS/RAW --format lcm --input_w /media/disk/MARTINLCM2/W/RAW --echo 30e-3 --fs 2000 --ft 63.3684e6 --basis_xml ../../basis/LCM2000.xml --output_txt results.txt --output_csv results.csv --output_image plot.eps --water_eddy false
