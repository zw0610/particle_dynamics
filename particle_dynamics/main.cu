#include "particles.cuh"

int main(int argc, char const *argv[]) {

    auto pset = solid::particles("test.xyz");

    pset.print_all();

    pset.delaunay();

    pset.evolve(0.01f);

    return 0;
}
