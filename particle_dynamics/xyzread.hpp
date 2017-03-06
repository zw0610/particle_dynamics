#ifndef xyz_reader_hpp
#define xyz_reader_hpp

#include <string>
#include <vector>

#include <iostream>
#include <fstream>

#include <cstdio>

template<typename T> class node {
public:
    const T x, y, z;
    node(const T inx, const T iny, const T inz) : x(inx), y(iny), z(inz) {};
    bool operator==(const node<T>& sn2) const;
    template<typename T2> node<T2> convert(void);

    void print(void) const;

};

template<typename T>
bool node<T>::operator==(const node<T>& sn2) const {
	return ( x==sn2.x && y==sn2.y && z==sn2.z );
}

template<typename T>
template<typename T2>
node<T2> node<T>::convert() {
	node<T2> t(stof(x), stof(y), stof(z));
	return t;
}

template<typename T>
void node<T>::print(void) const {
    std::cout << x << " " << y << " " << z << std::endl;
}

namespace reader {

	class xyz_reader {
	public:
        xyz_reader() {};
        ~xyz_reader() {};

		void read(std::string filename, std::vector< node<float> > &vertex_vec) {

            // -- Start of I/O -- //
    		std::ifstream ifile;
    		ifile.open(filename.c_str(), std::fstream::in);

    		if (!ifile) {
    			std::cout << "Cannot find " << filename << "!!!!" << std::endl;
    		} else {
    			std::cout << "Start to read " << filename << std::endl;
    		}


    		std::string temp_s;
    		while(std::getline(ifile,temp_s)) {
                float x, y, z;
    			if (sscanf(temp_s.c_str(), "%*s %f  %f  %f", &x, &y, &z)>=0) {
                    //std::cout <<  "Found!" << std::endl;
                    //std::cout <<  x << std::endl;
                    //std::cout <<  y << std::endl;
                    //std::cout <<  z << std::endl;
    				//try {

    					vertex_vec.push_back( node<float>( x, y, z ) );

    				//} catch (std::invalid_argument&) { std::cout << "cannot convert to float!" << std::endl; }
    			}
    		}

    		ifile.close();
    		std::cout << "Finished reading." << std::endl;

    		{
    			std::cout << "Total number of points is " << vertex_vec.size() << std::endl;
    		}

    		// -- End of I/O -- //

        }

	};

}


#endif
