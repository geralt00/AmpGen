source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_98python3 x86_64-centos7-gcc10-opt
main_dir=`pwd`
rm -rf build
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$main_dir -DCMAKE_CXX_STANDARD=14 ..
##
make -j80
make install -j80
cd $main_dir
