#pragma once
#include <vector>
#include <glm/glm.hpp>

#include <CGAL/Simple_cartesian.h>
#include <boost/functional/hash.hpp>
#include <CGAL/Delaunay_triangulation_2.h>
#include <unordered_map>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

// Define the kernel and the Delaunay triangulation types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;  // Kernel for geometry
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;  // Delaunay triangulation
typedef Delaunay::Face_handle Face_handle;  // Define Face_handle after Delaunay

// Additional types
typedef K::Point_2 Point_2;  // Point type in 2D for the triangulation


struct Vertex
{
    float x, y, z;       // Position
    float r, g, b;       // Color
    float u, v;          // Texture coordinates
    float normalx, normaly, normalz; // Normals
    unsigned int index;
	float friction; // Friction coefficient
};


class Mesh
{
public:
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> CubeMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> SphereMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> CylinderMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> ConeMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> TorusMesh(glm::vec3 color);
	std::pair<std::vector<Vertex>, std::vector<unsigned int>> PointCloud(glm::vec3 color);
	std::pair<std::vector<Vertex>, std::vector<unsigned int>> BSplineSurface(glm::vec3 color);
void MakeBiquadraticSurface(const int n_u, const int n_v, int d_u, int d_v, std::vector<float> mu, std::vector<float> mv, std::vector<glm::vec3> mc);
private:
	glm::vec3 deBoorSurface(int d_u, int d_v, std::vector<float> mu, std::vector<float> mv, std::vector<glm::vec3> mc, float u, float v, const int n_u, const int n_v);
    glm::vec3 deBoor(int k, int degree, const std::vector<float>& knots, std::vector<glm::vec3> controlPoints, float t);

    std::vector<Vertex> Readfile(const char* fileName, glm::vec3 color);


    float angleBetween(glm::vec3 v1, glm::vec3 v2);

    std::vector<Point_2> findSharedEdge(Point_2 v0_1, Point_2 v1_1, Point_2 v2_1, Point_2 v0_2, Point_2 v1_2, Point_2 v2_2);

    void swapEdge(Point_2 v0, Point_2 v1, Delaunay::Face_handle& triangle1, Delaunay::Face_handle& triangle2);

    bool isPoorQuality(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2);


    void flipEdgeIfNecessary(std::vector<Vertex>& vertices, Delaunay& triangulation, std::map<Point_2, unsigned int>& point_to_index);

    void flipTriangleEdges(Delaunay::Face_handle& triangle1, Delaunay::Face_handle& triangle2);

    void laplacianSmoothing(std::vector<Vertex>& vertices, std::vector<std::vector<unsigned int>>& vertexNeighbors, float lambda = 0.0f, int iterations = 1);

	std::vector<Vertex> m_vertices;
	std::vector<unsigned int> m_indices;
};
