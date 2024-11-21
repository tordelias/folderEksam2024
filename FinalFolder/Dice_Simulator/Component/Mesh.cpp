#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <omp.h>

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

  //Compute normals for each vertex
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

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::BSplineSurface(glm::vec3 color) {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    // Control points for the surface (example control grid)
    std::vector<Vertex> controlPoints = Readfile("Data/32-2-516-156-31.txt", color);

    int numU = 1600;
    int numV = 1200;
    int degree = 2;

    // Check if there are enough control points for the grid
    if (controlPoints.size() < (numU + 1) * (numV + 1)) {
        std::cerr << "Error: Not enough control points!" << std::endl;
        return { vertices, indices }; // Return empty result
    }

    // Corrected knot vectors
    std::vector<float> uKnots;
    std::vector<float> vKnots;

    // U Knot Vector
    for (int i = 0; i < (numU + degree + 1); i++) {
        if (i <= degree) {
            uKnots.push_back(0.0f);
        }
        else if (i >= numU) {
            uKnots.push_back(1.0f);  // Normalize the end knot to 1
        }
        else {
            uKnots.push_back((float)(i - degree) / (float)(numU - degree));  // Interior knots
        }
    }

    // V Knot Vector
    for (int i = 0; i < (numV + degree + 1); i++) {
        if (i <= degree) {
            vKnots.push_back(0.0f);
        }
        else if (i >= numV) {
            vKnots.push_back(1.0f);  // Normalize the end knot to 1
        }
        else {
            vKnots.push_back((float)(i - degree) / (float)(numV - degree));  // Interior knots
        }
    }

    // Resample control points into a grid-like structure
    std::vector<Vertex> gridControlPoints((numU + 1) * (numV + 1));
    int idx = 0;

    // Resample or bin the control points into the grid
    for (int i = 0; i <= numU; ++i) {
        for (int j = 0; j <= numV; ++j) {
            // Find a representative control point from the point cloud for each grid cell
            // Example method: pick the closest point to the grid position
            glm::vec3 gridPos(i / float(numU), j / float(numV), 0.0f);  // For simplicity, treat this as a 2D grid

            float minDist = std::numeric_limits<float>::max();
            int closestPointIdx = -1;

            // Loop over all points in the point cloud to find the closest one
            for (int k = 0; k < controlPoints.size(); ++k) {
                glm::vec3 controlPos(controlPoints[k].x, controlPoints[k].y, controlPoints[k].z);
                float dist = glm::distance(controlPos, gridPos);
                if (dist < minDist) {
                    minDist = dist;
                    closestPointIdx = k;
                }
            }

            // Assign the closest point to the grid control point
            gridControlPoints[idx++] = controlPoints[closestPointIdx];
        }
    }

    // Generate surface vertices using resampled control points
    for (int i = 0; i <= numU; ++i) {
        for (int j = 0; j <= numV; ++j) {
            float u = (float)i / numU;
            float v = (float)j / numV;

            glm::vec3 point(0.0f);

            // Use the resampled gridControlPoints to compute the B-spline surface
            for (int m = 0; m <= numU; ++m) {
                for (int n = 0; n <= numV; ++n) {
                    float basisU = BSplineBasis(m, degree, u, uKnots);
                    float basisV = BSplineBasis(n, degree, v, vKnots);

                    // Access the resampled control point grid
                    point += glm::vec3(gridControlPoints[m * (numV + 1) + n].x,
                        gridControlPoints[m * (numV + 1) + n].y,
                        gridControlPoints[m * (numV + 1) + n].z) * basisU * basisV;
                }
            }

            Vertex vertex;
            vertex.x = point.x;
            vertex.y = point.y;
            vertex.z = point.z;
            vertex.r = color.r;
            vertex.g = color.g;
            vertex.b = color.b;
            vertex.u = u;
            vertex.v = v;
            vertex.normalx = 0.0f;
            vertex.normaly = 0.0f;
            vertex.normalz = 1.0f;
            vertex.index = vertices.size();

            vertices.push_back(vertex);

            // Create indices for triangles
            if (i < numU && j < numV) {
                unsigned int idx1 = i * (numV + 1) + j;
                unsigned int idx2 = idx1 + 1;
                unsigned int idx3 = (i + 1) * (numV + 1) + j + 1;
                unsigned int idx4 = (i + 1) * (numV + 1) + j;

                indices.push_back(idx1);
                indices.push_back(idx2);
                indices.push_back(idx3); // Triangle 1

                indices.push_back(idx1);
                indices.push_back(idx3);
                indices.push_back(idx4); // Triangle 2
            }
        }
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

float Mesh::BSplineBasis(int i, int p, float t, const std::vector<float>& knots) {
    // Out-of-range check
    if (i < 0 || i >= knots.size() - 1)
    {
		std::cerr << "Index out of range!" << std::endl;
        return 0.0f;
    }

    // Base case: zero-degree basis function
    if (p == 0) {
        // Handle special case for the last knot span
        if (i == knots.size() - 2) {
            return (t >= knots[i] && t <= knots[i + 1]) ? 1.0f : 0.0f;
        }
        return (t >= knots[i] && t < knots[i + 1]) ? 1.0f : 0.0f;
    }

    // Recursive case: compute alpha and beta
    float alpha = 0.0f;
    if (knots[i + p] != knots[i]) {
        alpha = (t - knots[i]) / (knots[i + p] - knots[i]) * BSplineBasis(i, p - 1, t, knots);
    }

    float beta = 0.0f;
    if (knots[i + p + 1] != knots[i + 1]) {
        beta = (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * BSplineBasis(i + 1, p - 1, t, knots);
    }

    return alpha + beta;
}





std::vector<glm::vec3> Mesh::BarycentricCoordinates(std::vector<unsigned int> indices, std::vector<Vertex> vertices)
{
    std::vector<glm::vec3> result;
	float u, v, w;
    if (vertices.empty()) {
		std::cerr << "Vertices list is empty!" << std::endl;
		return std::vector<glm::vec3>();
    }
    for (int i = 0; i < indices.size(); i += 3) 
    {
        int index0 = indices[i];
        int index1 = indices[i + 1];
        int index2 = indices[i + 2];
		glm::vec3 v0(vertices[index0].x, vertices[index0].y, vertices[index0].z);
		glm::vec3 v1(vertices[index1].x, vertices[index1].y, vertices[index1].z);
		glm::vec3 v2(vertices[index2].x, vertices[index2].y, vertices[index2].z);

        glm::vec3 cpoint = (v0 + v1 + v2) / 3.0f; // get center of triangle 

        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0v2 = v2 - v0;
        glm::vec3 v0p = cpoint - v0;

        // Computing dot products
        double dot00 = glm::dot(v0v1, v0v1);
        double dot01 = glm::dot(v0v1, v0v2);
        double dot02 = glm::dot(v0v1, v0p);
        double dot11 = glm::dot(v0v2, v0v2);
        double dot12 = glm::dot(v0v2, v0p);

        // Computing barycentric coordinates
        double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        double v = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double w = (dot00 * dot12 - dot01 * dot02) * invDenom;
        double u = 1 - v - w;
		result.push_back(glm::vec3(u, v, w));
    }



    return result;
}




