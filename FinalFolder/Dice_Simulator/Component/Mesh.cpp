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
    // Read the point cloud from the file (with progress bar in Readfile)
    std::vector<Vertex> vertices = Readfile("Data/32-2-516-156-31.txt", color);
    std::vector<unsigned int> indices;

    // Step 1: Convert 3D points to 2D points (using x, z)
    std::vector<Point_2> points2D;
    std::unordered_map<Point_2, unsigned int> pointIndexMap;

    // Create the 2D points and the index map (mapping from 2D points to 3D indices)
    size_t totalVertices = vertices.size();
    size_t processedVertices = 0;

    // Progress bar setup
    int barWidth = 50;
    std::cout << "adding points to Delaunay..." << std::endl;

    for (unsigned int i = 0; i < vertices.size(); ++i) {
        Point_2 p2d(vertices[i].x, vertices[i].z);
        points2D.push_back(p2d);
        pointIndexMap[p2d] = i;  // Store the original index

        // Update the progress bar every 1000 points or at the last point
        processedVertices++;
        if (processedVertices % 10000 == 0 || processedVertices == totalVertices) {
            float progress = static_cast<float>(processedVertices) / totalVertices;
            int pos = barWidth * progress;
            std::cout << "\r[";  // Start overwriting the same line
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

    // After the first progress bar is complete, print a newline to clear the line
    std::cout << std::endl;

    // Step 2: Perform Delaunay Triangulation on the 2D points
    Delaunay dt;
    dt.insert(points2D.begin(), points2D.end());

    // Step 3: Generate indices based on the triangulation
    size_t totalFaces = 0;
    for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
        totalFaces++;
    }

    size_t processedFaces = 0;
    std::cout << "Triangulating faces..." << std::endl;

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

        // Update the progress bar every 100 faces or at the last face
        processedFaces++;
        if (processedFaces % 10000 == 0 || processedFaces == totalFaces) {
            float progress = static_cast<float>(processedFaces) / totalFaces;
            int pos = barWidth * progress;
            std::cout << "\r[";  // Overwrite the same line for progress bar
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

    // After the second progress bar is complete, print a newline to finish
    std::cout << std::endl;

    return { vertices, indices };
}


std::vector<Vertex> Mesh::Readfile(const char* fileName, glm::vec3 color) {
    std::ifstream inputFile(fileName);
    std::vector<Vertex> pointCloud;
    float min_x = -816.02, max_x = 783.98;
    float min_z = -620.771, max_z = 579.229;
    int processedLines = 0;
    int totalLines = 2531030;  // Total number of lines in the file

    if (inputFile.is_open()) {
        std::string line;
        std::getline(inputFile, line);  // Skip header line if there is one
        Vertex point;

        point.r = color.x;
        point.g = color.y;
        point.b = color.z;

        // Progress bar setup
        int barWidth = 50;
        std::cout << "Loading points..." << std::endl;

        while (std::getline(inputFile, line))
        {
            if (sscanf_s(line.c_str(), "%f %f %f", &point.x, &point.z, &point.y) == 3)
            {
                point.x -= 608016.02;
                point.y -= 336.8007;
                point.z -= 6750620.771;

                point.u = (point.x - min_x) / (max_x - min_x);  // Normalize x coordinate
                point.v = (point.z - min_z) / (max_z - min_z);  // Normalize z coordinate

                pointCloud.push_back(point);
            }

            processedLines++;
            if (processedLines % 100000 == 0 || processedLines == totalLines) {
                float progress = static_cast<float>(processedLines) / totalLines;
                int pos = barWidth * progress;
                std::cout << "\r[";  // Start overwriting the same line
                for (int i = 0; i < barWidth; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << int(progress * 100.0f) << "%";
                std::flush(std::cout);
            }
        }

        // Ensure the progress bar is at 100% when done
        float progress = 100;
        int pos = barWidth * progress;
        std::cout << "\r[";  // Start overwriting the same line
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress) << "%";
        std::flush(std::cout);

        // After the progress bar is complete, print a newline
        std::cout << std::endl;
        inputFile.close();
    }
    else {
        std::cerr << "Unable to open the input file for reading." << std::endl;
    }

    return pointCloud;
}



