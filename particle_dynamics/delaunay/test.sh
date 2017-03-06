cd build
rm -rf *
rm -rf ../main/libparallel_delaunay.a ../main/exec
cmake ..
make

#cp ./libparallel_delaunay.a ../main
#cd ../main
#g++ main.cxx libparallel_delaunay.a -lCGAL -lgmp -lboost_thread -lboost_system -lpthread -ltbb -ltbbmalloc -o exec
#./exec
