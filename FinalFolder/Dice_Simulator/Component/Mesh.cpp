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
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
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

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::PointCloud(glm::vec3 color)
{
    std::vector<Vertex> vertices = Readfile("Data/PointCloud.txt", color);
    std::vector<unsigned int> indices;

    int numCols = sqrt(vertices.size());
    int numRows = sqrt(vertices.size());

    // Generate indices for grid pattern (two triangles per quad)
    for (int y = 0; y < numRows - 1; ++y) {
        for (int x = 0; x < numCols - 1; ++x) {
            int topLeft = y * numCols + x;
            int topRight = topLeft + 1;
            int bottomLeft = (y + 1) * numCols + x;
            int bottomRight = bottomLeft + 1;

            // First triangle
            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            // Second triangle
            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }


    return { vertices, indices };
}

// Readfile function
std::vector<Vertex> Mesh::Readfile(const char* fileName, glm::vec3 color) {
    std::ifstream inputFile(fileName);
    std::vector<Vertex> pointCloud;
    std::unordered_map<std::pair<float, float>, float, PairHash> heightMap;

    if (inputFile.is_open()) {
        std::string line;
        std::getline(inputFile, line);  // Skip header line if there is one
        Vertex point;
        char comma;
        float prevY = 0;
		int skip = 0;

        while (inputFile >> point.x >> comma >> point.z >> comma >> point.y) {
            if (skip == 10) {
                point.x -= 608016.02;
                point.y -= 336.8007;
                point.z -= 6750620.771;

                // Set color based on height comparison
                if (point.y > prevY) {
                    point.r = 0.0f;
                    point.g = 1.0f;
                    point.b = 0.0f;
                }
                else {
                    point.r = 1.0f;
                    point.g = 0.0f;
                    point.b = 0.0f;
                }

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

    // Perform Delaunay Triangulation to filter points based on (x, z)
    std::vector<Vertex> filteredPoints;
    std::vector<Point_2> cgalPoints;

    // Convert pointCloud to CGAL 2D points for triangulation
    for (const Vertex& v : pointCloud) {
        cgalPoints.push_back(Point_2(v.x, v.z));
    }

    Delaunay dt;
    dt.insert(cgalPoints.begin(), cgalPoints.end());

    // Iterate over triangulation vertices and add boundary points to filteredPoints
    for (auto v = dt.finite_vertices_begin(); v != dt.finite_vertices_end(); ++v) {
        Vertex filteredVertex;
        filteredVertex.x = v->point().x();
        filteredVertex.z = v->point().y();

        // Look up the original y (height) value from the heightMap
        auto it = heightMap.find({ filteredVertex.x, filteredVertex.z });
        if (it != heightMap.end()) {
            filteredVertex.y = it->second;
        }
        else {
            filteredVertex.y = 0.0f;  // Fallback if y not found
        }

        // Set color for filtered points if desired
        filteredVertex.r = 1.0f;
        filteredVertex.g = 1.0f;
        filteredVertex.b = 1.0f;

        filteredPoints.push_back(filteredVertex);
    }

    std::cout << "Filtered Point Cloud Size: " << filteredPoints.size() << std::endl;
    return filteredPoints;
}
