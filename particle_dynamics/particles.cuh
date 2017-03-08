#include "xyzread.hpp"
#include "delaunay_triangulation.hpp"
#include "kernel_functions.cuh"
#include "kernel_helper.cuh"

#include <cublas_v2.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>

//Struct used for compund data struct, uint2, sorting
struct compare_uint2
{
	__host__ __device__ bool operator() (uint2 a, uint2 b) { return (a.x!=b.x)?(a.x<b.x):(a.y<b.y); }
};
//Struct used for compund data struct, uint2, unique
struct equal_uint2
{
	__host__ __device__ bool operator() (uint2 a, uint2 b) { return (a.x==b.x)&&(a.y==b.y); }
};

// Using Quaternion

// -- Transformation Quaternion -- //
namespace quaternion {

    void qt_init(float * hq, float *dq) {
        hq = (float*) malloc (4*4*sizeof(float));
        for (size_t i = 0; i<4; i++) {
            for (size_t j = 0; j<4; j++) {
                if (i==j) hq[IDX2C(i,j,4)]=1.0f;
            }
        }
        cudaMalloc ((void **)&dq ,4*4*sizeof(float));
    }

    class trans_quat {
        float * host_quat;
        float * device_quat;
    public:
        trans_quat() {
            qt_init(host_quat, device_quat);
        }

        ~trans_quat() {
            free(host_quat);
            cudaFree(device_quat);
        }

        void set(float translation [] ) {
            qt_init(host_quat, device_quat);
            for (size_t i = 0; i<3; i++) {
                host_quat[IDX2C(i,3,4)] = translation[i];
            }
        }

        void set(float degree, float rotation_axis [] );

        float * dev(void) {
            cublasSetMatrix(4, 4, sizeof(float), host_quat, 4, device_quat, 4);
            return device_quat;
        }

    };

}

namespace solid {

    class particles {
        std::vector< quaternion::trans_quat > state_list;
        float * host_list;
        float * device_list;
        float * device_vel;
        float * device_acc;
        size_t num_points;
        cublasStatus_t synchronize(void);
        uint2 *dev_edges;
        size_t num_edges;
        unsigned int *head_and_tail;
        unsigned int *empty_list;


    public:
        particles(std::string filename);
        ~particles();
        cublasStatus_t apply_trans(quaternion::trans_quat & qt);
        void print_all(bool update = false);
        void delaunay(void);
		void evolve(float dt);

    };


}




namespace solid {

    particles::particles(std::string filename) {

        std::vector< node<float> > vn;

        reader::xyz_reader xyzr;
    	xyzr.read(filename, vn);

        //for (auto item : vn) { item.print();}

        host_list = (float*) malloc (4*vn.size()*sizeof(float));
        for (size_t i = 0; i<vn.size(); i++) {
            //std::cout << IDX2C(i,0,4) << " " IDX2C(i,1,4) << " " << IDX2C(i,2,4) << " " << IDX2C(i,3,4) << std::endl;
            host_list[IDX2C(0,i,4)]=vn[i].x;
            host_list[IDX2C(1,i,4)]=vn[i].y;
            host_list[IDX2C(2,i,4)]=vn[i].z;
            host_list[IDX2C(3,i,4)]=0.0f;

        }
        cudaMalloc ((void **)&device_list ,4*vn.size()*sizeof(float));
        cublasSetMatrix(4, vn.size(), sizeof(float), host_list, 4, device_list, 4);
        num_points = vn.size();

        cudaMalloc ((void **)&device_acc ,4*vn.size()*sizeof(float));
        cudaMalloc ((void **)&device_vel ,4*vn.size()*sizeof(float));
        float * host_temp = new float[4*vn.size()];
        for (size_t i = 0; i<4*vn.size(); i++) { host_temp[i] = 0.0f; }
        cudaMemcpy(device_acc, host_temp, 4*vn.size()*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(device_vel, host_temp, 4*vn.size()*sizeof(float), cudaMemcpyHostToDevice);

        cudaMalloc ((void **)&head_and_tail ,num_points * sizeof(unsigned int));

        empty_list = new unsigned int[num_points];
        for (size_t i = 0; i<num_points; i++) { empty_list[i] = 0; }

    }

    particles::~particles() {
        free(host_list);
        delete [] empty_list;
        cudaFree(device_list);
        cudaFree(head_and_tail);
        cudaFree(dev_edges);
        cudaFree(device_acc);
        cudaFree(device_vel);
    }

    cublasStatus_t particles::apply_trans(quaternion::trans_quat & qt) {
        cublasStatus_t  stat;                //  CUBLAS  functions  status
        cublasHandle_t  handle;              //  CUBLAS  context
        cublasCreate (& handle );

        float * device_temp;
        cudaMalloc((void **)&device_temp ,4*num_points*sizeof(float));
        cublasSetMatrix(4, num_points, sizeof(float), host_list, 4, device_temp, 4);

        float  al=1.0f;                      // al=1
    	float  bet =1.0f;                    //bet=1
        stat = cublasSgemm(handle,CUBLAS_OP_N,CUBLAS_OP_N, 4, num_points, 4, &al, qt.dev(), 4, device_temp, 4, &bet, device_list, 4);

        cudaFree(device_temp);
        stat = cublasGetMatrix(4, num_points, sizeof(float), device_list, 4, host_list, 4);
        stat = cublasDestroy(handle);               //  destroy  CUBLAS  context

        return stat;
    }

    cublasStatus_t particles::synchronize(void) {
        const cublasStatus_t stat = cublasGetMatrix(4, num_points, sizeof(float), device_list, 4, host_list, 4);
        return stat;
    }

    void particles::print_all(bool update) {
        if (update) {synchronize();}
        for (size_t i = 0; i<(4*num_points); i++) {
            if ( (i%4)!=3 ) {
                std::cout << host_list[i] << " ";
            } else {
                std::cout << std::endl;
            }
        }
    }
    // delaunay member function is all checked
    void particles::delaunay(void) {

        synchronize();
        std::vector<unsigned int> hos_teh_list;
        unsigned int num_valid_cells = delaunay_triangulation( host_list, num_points, hos_teh_list );

        //If you wish to print out the index of all cells
/*
        for (size_t i = 0; i<num_valid_cells; i++) {
		for (size_t j = 0; j<4; j++) {
			std::cout << hos_teh_list[i*4+j] << " ";
		}
		std::cout << std::endl;
*/

        //Dessolve Teh List into Edge List
        thrust::device_vector< unsigned int > dev_vec_teh = hos_teh_list;
        thrust::device_vector< uint2 >        dev_vec_edg(num_valid_cells*6*2);

        Dissolve<<<num_valid_cells, 6>>>(   thrust::raw_pointer_cast( dev_vec_teh.data() ) ,
                                            thrust::raw_pointer_cast( dev_vec_edg.data() ) );
/*
        std::cout << " num_valid_cells = " << num_valid_cells << std::endl;
        std::cout << " hos_teh_list.size() = " << hos_teh_list.size() << std::endl;
        std::cout << " dev_vec_teh.size() = " << dev_vec_teh.size() << std::endl;
        std::cout << " dev_vec_edg.size) = " << dev_vec_edg.size() << std::endl;
        std::cout << "delaunay" << std::endl;
        for (size_t i = 0; i<hos_teh_list.size(); i++) {
            std::cout << hos_teh_list[i] << " ";
            if ( i%4==0 ) { std::cout << std::endl; }
        }
        std::cout << std::endl;
*/
        //Sort with uint2
        compare_uint2 comp;
        thrust::sort( dev_vec_edg.begin(), dev_vec_edg.end(), comp);

        equal_uint2 uiequal;
	    auto new_end = thrust::unique(thrust::device, dev_vec_edg.begin(), dev_vec_edg.end(), uiequal);
        num_edges = std::distance( dev_vec_edg.begin(), new_end);
/*
        std::cout << "cleaned edges" << std::endl;
        for (size_t i = 0; i<num_edges; i++) {
            uint2 item = dev_vec_edg[i];
            std::cout << i << " " << item.x << " " << item.y << std::endl;
        }
*/
        //Store the data into member varialb uint2 dev_edges
        cudaFree(dev_edges);
        cudaMalloc( (void**)&dev_edges, sizeof(uint2) * num_edges);
        thrust::copy(dev_vec_edg.begin(), new_end, thrust::device_pointer_cast(dev_edges));

        //Determin the start index and end index in dev_edges for each point.
        cudaMemcpy(head_and_tail, empty_list, num_points*sizeof(unsigned int), cudaMemcpyHostToDevice);
        LookLeft<<<(num_edges>>10)+1, 1024>>>(dev_edges, head_and_tail, num_edges);
        //std::cout << "Head and Tail" << std::endl;
        //print_dev_array<unsigned int>(head_and_tail, num_points, 1);
	}

	void particles::evolve(float dt) {
        // float4 * pos, float4 * vel, float4 * acc, unsigned int * tail, size_t num_points, uint2 * edge_list, size_t num_edges, float odt
		GetAcceleration<<<num_points, 1>>>( (float4 *) device_list, (float4 *) device_vel, (float4 *) device_acc, head_and_tail, num_points, dev_edges, num_edges, 1.0f/dt );
		//( float4 * pos, float4 * vel, float4 * acc, float dt)
        UpdateVelPos<<<num_points, 1>>>( (float4 *) device_list, (float4 *) device_vel, (float4 *) device_acc, dt );
		//( float4 * pos, float4 * vel)
        CheckWall<<<num_points, 1>>>( (float4 *) device_list, (float4 *) device_vel );
	}





}
