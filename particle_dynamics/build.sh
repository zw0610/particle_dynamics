cd delaunay
sh test.sh

cp ./build/libparallel_delaunay.a ../

cd ..

nvcc main.cu libparallel_delaunay.a -std=c++11 -lcublas -lCGAL -lgmp -lboost_thread -lboost_system -lpthread -ltbb -ltbbmalloc -o exec
