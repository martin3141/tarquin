# TARQUIN compilation instructions

## Ubuntu 18.04 and MX Linux 18

### install packages
```
sudo apt-get install -y build-essential cmake liblapack-dev gfortran libf2c2-dev libfftw3-dev gnuplot libboost-filesystem-dev libboost-thread-dev libboost-test-dev libqt4-dev libqwt5-qt4-dev git
```

### download TARQUIN source core
```
git clone https://github.com/martin3141/tarquin.git
cd tarquin/src
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../ 
make
```

## Arch Linux and Manjaro
```
git clone https://aur.archlinux.org/tarquin.git
cd tarquin
makepkg -si
```
