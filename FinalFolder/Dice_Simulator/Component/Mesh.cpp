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
        // Front face (z = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  0.0f, -1.0f, 0 }, // Bottom-left
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  0.0f, -1.0f, 1 }, // Bottom-right
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  0.0f, -1.0f, 2 }, // Top-right
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  0.0f, -1.0f, 3 }, // Top-left

        // Back face (z = 0.5f)
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  0.0f,  1.0f, 4 }, // Bottom-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  0.0f,  1.0f, 5 }, // Bottom-left
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  0.0f,  1.0f, 6 }, // Top-left
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  0.0f,  1.0f, 7 }, // Top-right

        // Left face (x = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f, -1.0f,  0.0f,  0.0f, 8 }, // Bottom-right
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f, -1.0f,  0.0f,  0.0f, 9 }, // Top-right
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f, -1.0f,  0.0f,  0.0f, 10 }, // Top-left
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f, -1.0f,  0.0f,  0.0f, 11 }, // Bottom-left

        // Right face (x = 0.5f)
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  1.0f,  0.0f,  0.0f, 12 }, // Bottom-left
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  1.0f,  0.0f,  0.0f, 13 }, // Top-left
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  1.0f,  0.0f,  0.0f, 14 }, // Top-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  1.0f,  0.0f,  0.0f, 15 }, // Bottom-right

        // Top face (y = 0.5f)
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  1.0f,  0.0f, 16 }, // Top-left
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  1.0f,  0.0f, 17 }, // Top-right
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  1.0f,  0.0f, 18 }, // Bottom-right
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  1.0f,  0.0f, 19 }, // Bottom-left

        // Bottom face (y = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f, -1.0f,  0.0f, 20 }, // Bottom-left
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f, -1.0f,  0.0f, 21 }, // Bottom-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f, -1.0f,  0.0f, 22 }, // Top-right
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f, -1.0f,  0.0f, 23 }  // Top-left
    };

    std::vector<unsigned int> indices = {
        0, 1, 2, 2, 3, 0,   // Front
        4, 5, 6, 6, 7, 4,   // Back
        8, 9, 10, 10, 11, 8, // Left
        12, 13, 14, 14, 15, 12, // Right
        16, 17, 18, 18, 19, 16, // Top
        20, 21, 22, 22, 23, 20  // Bottom
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
            std::cout << "\r[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

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
            std::cout << "\r[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

    std::cout << std::endl;

    // Step 4: Compute normals for each vertex
    for (auto& vertex : vertices) {
        vertex.normalx = 0;
        vertex.normaly = 0;
        vertex.normalz = 0;
    }

    for (size_t i = 0; i < indices.size(); i += 3) {
        unsigned int idx0 = indices[i];
        unsigned int idx1 = indices[i + 1];
        unsigned int idx2 = indices[i + 2];

        Vertex& v0 = vertices[idx0];
        Vertex& v1 = vertices[idx1];
        Vertex& v2 = vertices[idx2];

        glm::vec3 p0(v0.x, v0.y, v0.z);
        glm::vec3 p1(v1.x, v1.y, v1.z);
        glm::vec3 p2(v2.x, v2.y, v2.z);

        glm::vec3 edge1 = p1 - p0;
        glm::vec3 edge2 = p2 - p0;

        glm::vec3 faceNormal = glm::normalize(glm::cross(edge1, edge2));

        v0.normalx += faceNormal.x;
        v0.normaly += faceNormal.y;
        v0.normalz += faceNormal.z;

        v1.normalx += faceNormal.x;
        v1.normaly += faceNormal.y;
        v1.normalz += faceNormal.z;

        v2.normalx += faceNormal.x;
        v2.normaly += faceNormal.y;
        v2.normalz += faceNormal.z;
    }

    for (auto& vertex : vertices) {
        glm::vec3 normal(vertex.normalx, vertex.normaly, vertex.normalz);
        normal = glm::normalize(normal);
        vertex.normalx = normal.x;
        vertex.normaly = normal.y;
        vertex.normalz = normal.z;
    }

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



