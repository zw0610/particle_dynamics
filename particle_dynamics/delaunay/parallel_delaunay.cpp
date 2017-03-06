#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <set>

//Type Define Area

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_data_structure_3<
CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K>,
CGAL::Triangulation_cell_base_3<K>,
CGAL::Parallel_tag>                          Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds> Triangulation;

typedef K::Point_3          Point;



int delaunay_triangulation( float *coord4, const unsigned int num_vertex, std::vector<unsigned int> & teh ) {
    
    std::vector<Point> V;
    V.reserve(num_vertex);
    for (int i = 0; i != num_vertex; ++i) {
        V.push_back( Point(coord4[i*4], coord4[i*4+1], coord4[i*4+2]) );
    }
    
    // Construct the locking data-structure, using the bounding-box of the points
    //Triangulation::Lock_data_structure locking_ds( CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
    // Construct the triangulation in parallel
    Triangulation T; //(V.begin(), V.end()/*, &locking_ds*/);

	for (size_t i = 0; i<num_vertex; i++) {
		T.insert(V[i]);
	}
    
    assert(T.is_valid());

	std::cout << "Total number of vertices is " << T.number_of_vertices() << std::endl;
	std::cout << "Total number of cells is " << T.number_of_cells() << std::endl;
	//std::cout << "Total number of faces is " << T.number_of_facets() << std::endl;
	//std::cout << "Total number of edges is " << T.number_of_edges() << std::endl;

	//Triangulation::Edge_iterator it = T.edges_begin();
	unsigned int count = 0;
	for ( Triangulation::Finite_vertices_iterator it = T.finite_vertices_begin(); it != T.finite_vertices_end(); ++it) {
		it->info() = count;		
		std::cout << "idx: " << it->info() << "  coord: " << it->point().x() << " " << it->point().y() << " " << it->point().z() << std::endl;
		count++;
	}
/*
	std::cout << "Total count = " << count << std::endl;

	for ( Triangulation::Edge_iterator it = T.edges_begin() ; it != T.edges_end() ; ++it) {

		Triangulation::Edge e=(*it);

		//std::cout << it->first->vertex( (it->second+1)%3 ) << std::endl;

		int i1 = e.first->vertex( (it->second+1)%3 )->info();
        int i2 = e.first->vertex( (it->second+2)%3 )->info();

		std::cout << i1 << "  " << i2 << std::endl;  
	}
*/
	for ( Triangulation::Cell_iterator it = T.cells_begin() ; it != T.cells_end() ; ++it) {
		//std::cout << it->first->vertex( (it->second+1)%3 ) << std::endl;
		assert(it->is_valid());
		unsigned int i0 = it->vertex( 0 )->info();
        unsigned int i1 = it->vertex( 1 )->info();
		unsigned int i2 = it->vertex( 2 )->info();
		unsigned int i3 = it->vertex( 3 )->info();

		if ( i0!=i1 && i0!=i2 && i0!=i3 && i1!=i2 && i1!=i3 && i2!=i3 ) {
			teh.push_back(i0);
			teh.push_back(i1);
			teh.push_back(i2);
			teh.push_back(i3);
		}

		//std::cout << i0 << "  " << i1 << "  " << i2 << "  " << i3 << std::endl;  
	}

	//std::ofstream oFile("output",std::ios::out);
	//oFile << T;

    //std::cout << T.number_of_edges() << std::endl;
  
    //return T.number_of_;

	return (teh.size())/4;
    
}
