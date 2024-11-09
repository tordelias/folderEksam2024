#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/IO/Polyhedron_iostream.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> BaseVertex;
typedef CGAL::Triangulation_vertex_base_2<K, CGAL::Triangulation_ds_vertex_base_2<BaseVertex>> VertexBase;
typedef CGAL::Triangulation_data_structure_2<VertexBase> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;

// Define a custom vertex with index
typedef Delaunay::Vertex_handle Vertex_handle_with_index;
typedef K::Point_2 Point_2;

// Custom hash function for std::pair<float, float>
struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        return std::hash<T1>()(p.first) ^ (std::hash<T2>()(p.second) << 1);
    }
};

// Cube mesh generation
std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::CubeMesh(glm::vec3 color) {
    std::vector<Vertex> vertices = {
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f }, // Vertex 0 (Front face)
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f }, // Vertex 1
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f, 0.0f, 0.0f, -1.0f }, // Vertex 2
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f, 0.0f, 0.0f, -1.0f }, // Vertex 3
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f, 0.0f, 0.0f,  1.0f }, // Vertex 4 (Back face)
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f, 0.0f, 0.0f,  1.0f }, // Vertex 5
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f, 0.0f, 0.0f,  1.0f }, // Vertex 6
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f, 0.0f, 0.0f,  1.0f }  // Vertex 7
    };

    std::vector<unsigned int> indices = {
        0, 1, 2, 2, 3, 0,   // Front face
        4, 5, 6, 6, 7, 4,   // Back face
        0, 1, 5, 5, 4, 0,   // Bottom face
        2, 3, 7, 7, 6, 2,   // Top face
        1, 2, 6, 6, 5, 1,   // Right face
        3, 0, 4, 4, 7, 3    // Left face
    };

    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::SphereMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;


    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::CylinderMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;


    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::ConeMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;



    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::TorusMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;



    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::PointCloud(glm::vec3 color) {
    std::vector<Vertex> vertices = Readfile("Data/PointCloud.txt", color);
    std::vector<unsigned int> indices;

    // Step 1: Convert 3D points to 2D points (using x, z)
    std::vector<Point_2> points2D;
    std::unordered_map<Point_2, unsigned int> pointIndexMap;

    // Create the 2D points and the index map (mapping from 2D points to 3D indices)
    for (unsigned int i = 0; i < vertices.size(); ++i) {
        Point_2 p2d(vertices[i].x, vertices[i].z);
        points2D.push_back(p2d);
        pointIndexMap[p2d] = i;  // Store the original index
    }

    // Step 2: Perform Delaunay Triangulation on the 2D points
    Delaunay dt;
    dt.insert(points2D.begin(), points2D.end());

    // Step 3: Generate indices based on the triangulation
    for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
        // Get the three vertices of the triangle in the 2D plane
        Point_2 p0 = f->vertex(0)->point();
        Point_2 p1 = f->vertex(1)->point();
        Point_2 p2 = f->vertex(2)->point();

        // Look up the indices of the original vertices in the pointIndexMap
        unsigned int idx0 = pointIndexMap[p0];
        unsigned int idx1 = pointIndexMap[p1];
        unsigned int idx2 = pointIndexMap[p2];

        // Add indices to the list (three vertices per face)
        indices.push_back(idx0);
        indices.push_back(idx1);
        indices.push_back(idx2);
    }

    // Debug print the generated indices
    //std::cout << "Generated " << indices.size() / 3 << " triangles:" << std::endl;
    //for (size_t i = 0; i < indices.size(); i += 3) {
    //    std::cout << "Triangle: (" << indices[i] << ", " << indices[i + 1] << ", " << indices[i + 2] << ")" << std::endl;
    //}

    return { vertices, indices };
}


// Updated Readfile function with CGAL vertex indices
std::vector<Vertex> Mesh::Readfile(const char* fileName, glm::vec3 color) {
    std::ifstream inputFile(fileName);
    std::vector<Vertex> pointCloud;

    if (inputFile.is_open()) {
        std::string line;
        std::getline(inputFile, line);  // Skip header line if there is one
        Vertex point;
        char comma;
        int skip = 0;
        float prevY = 0;

        while (inputFile >> point.x >> comma >> point.z >> comma >> point.y) {
            if (skip == 1) {
                point.x -= 608016.02;
                point.y -= 336.8007;
                point.z -= 6750620.771;

                if (point.y > prevY) {
                    point.r = 0.0f;
                    point.g = 1.0f;
                    point.b = 0.0f;
                }

            else 
                {
                point.r = 1.0f;
                point.g = 0.0f;
                point.b = 0.0f;
                 }


                // Add the point to the pointCloud
				prevY = point.y;
                pointCloud.push_back(point);
                skip = 1;
            }
            else {
                ++skip;
            }
        }
        inputFile.close();
    }
    else {
        std::cerr << "Unable to open the input file for reading." << std::endl;
    }

    return pointCloud;
}
