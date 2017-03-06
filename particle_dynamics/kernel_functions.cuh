#ifndef kernel_functions_cuh
#define kernel_functions_cuh

#include "helper_math.h"

#define  IDX2C(i,j,ld) (((j)*(ld))+( i ))

#define Rebound 0.01f
#define gx      9.8f
#define gy      0.0f
#define gz      0.0f

__global__ void VertexToCell(float4 * pt_list, unsigned int * ibm_list) {
    float px = pt_list[blockIdx.x].x;
    float py = pt_list[blockIdx.x].y;
    float pz = pt_list[blockIdx.x].z;
    ibm_list[blockIdx.x] = (int)(px+py+pz);
}


__device__ void assign( uint2 &edge, const unsigned int id0, const unsigned int id1 ) {
    edge.x = id0;
    edge.y = id1;
}

__global__ void Dissolve(unsigned int * tehs, uint2 * edges) {
    const size_t idx_tehs = 4 * blockIdx.x;
    const size_t idx_edge = 6 * blockIdx.x;
    const size_t idx_local = threadIdx.x;   //BlockDim.x = 6;

    unsigned int vex_id0 = 0;
    unsigned int vex_id1 = 0;

    switch (idx_local) {
        case 0:
            vex_id0 = tehs[idx_tehs + 0];
            vex_id1 = tehs[idx_tehs + 1];
            break;
        case 1:
            vex_id0 = tehs[idx_tehs + 0];
            vex_id1 = tehs[idx_tehs + 2];
            break;
        case 2:
            vex_id0 = tehs[idx_tehs + 0];
            vex_id1 = tehs[idx_tehs + 3];
            break;
        case 3:
            vex_id0 = tehs[idx_tehs + 1];
            vex_id1 = tehs[idx_tehs + 2];
            break;
        case 4:
            vex_id0 = tehs[idx_tehs + 1];
            vex_id1 = tehs[idx_tehs + 3];
            break;
        case 5:
            vex_id0 = tehs[idx_tehs + 2];
            vex_id1 = tehs[idx_tehs + 3];
            break;
    }

    assign( edges[idx_edge + idx_local    ], vex_id0, vex_id1 );
    assign( edges[idx_edge + idx_local + 6], vex_id1, vex_id0 );

}

__global__ void LookLeft (uint2 * edge_list, unsigned int * tail, size_t num_edges) {

    const size_t idx = (blockIdx.x<<10) + threadIdx.x;

    if ( idx < num_edges /* Valid Thread */ ) {

        if ( idx == 0 ) {
            tail[0] = 0;
        } else {
            if ( edge_list[idx-1].x != edge_list[idx].x ) {
                /*This is the head index*/
                tail[edge_list[idx].x] = idx;  }
        }

    }

}


__device__ void cross_f4(float4 &u, float4 &v, float4 &out) {
    out.x = u.y*v.z - u.z*v.y;
    out.y = u.z*v.x - u.x*v.z;
    out.z = u.x*v.y - u.y*v.x;
}
__device__ void norm_f4(float4 & vec) {
    float r = 1.0f/sqrtf((vec.x*vec.x + vec.y*vec.y + vec.z*vec.z));
    vec.x = vec.x*r; vec.y = vec.y*r; vec.z = vec.z*r;
}
__device__ float dot_f4(float4 & u, float4 & v) {
    return ( u.x*v.x + u.y*v.y + u.z*v.z);
}
__device__ void scale_f4(float4 & u, float s) {
    u.x *= s;
    u.y *= s;
    u.z *= s;
}

__global__ void GetAcceleration( float4 * pos, float4 * vel, float4 * acc, unsigned int * tail, size_t num_points, uint2 * edge_list, size_t num_edges, float odt) {
    const size_t idx = threadIdx.x;
    const size_t start = tail[idx];
    const size_t stop  = (idx==(num_points-1))?num_edges:tail[idx+1];

    float4 local_acc = make_float4(0.0f, 0.0f, 0.0f, 0.0f);

    for (size_t i = start; i<stop; i++) {
        const size_t idx2 = edge_list[i].y;
        //if (idx!=edge_list[i].x) { printf("%s\n", "Error!"); }

        float4 oo = pos[idx] - pos[idx2];
        /*Here, r1 = r2 = 0.5*/
        if ( dot_f4(oo, oo)>=(0.5+0.5)*(0.5+0.5)  ) {
            continue;
        }

        norm_f4(oo);
        float4 v2 = vel[idx2] - vel[idx];
/*
        float4 np;
        cross_f4(oo, v2, np);
        norm_f4(np);

        float4 nr;
        cross_f4(oo, np, nr);
        norm_f4(nr);
*/
        float mag1 = dot_f4(v2, oo);
        //float mag2 = dot_f4(v2, nr);

        scale_f4(oo, mag1*odt);
        local_acc += oo;
    }

    acc[idx] = local_acc;
}

__global__ void UpdateVelPos( float4 * pos, float4 * vel, float4 * acc, float dt) {
    const size_t idx = blockIdx.x;

    float4 temp_acc = acc[idx] + make_float4(gx, gy, gz, 0.0f);
    scale_f4(temp_acc, dt);

    float4 temp_vel = vel[idx] + temp_acc;
    vel[idx] = temp_vel;

    scale_f4(temp_vel, dt);
    pos[idx] += temp_vel;
}

//Here We assume this is a cube with point 0,0,0 to 10,10,10
#define XMAX 10.0f
#define YMAX 10.0f
#define ZMAX 10.0f

__device__ void OneDCheck( float pos, float dis, float & vel) {

    //Still the radius is 0.5
    if ( (pos-0.5f)<0.0f && (vel<0) ) {
        vel = -vel;
    } else if ( (pos+0.5f)>dis && (vel>0) ) {
        vel = -vel;
    }

}

__global__ void CheckWall ( float4 * pos, float4 * vel) {

    const size_t idx = blockIdx.x;

    // For X Direction
    OneDCheck( pos[idx].x, XMAX, vel[idx].x );
    OneDCheck( pos[idx].y, YMAX, vel[idx].y );
    OneDCheck( pos[idx].z, ZMAX, vel[idx].z );

}

#endif
