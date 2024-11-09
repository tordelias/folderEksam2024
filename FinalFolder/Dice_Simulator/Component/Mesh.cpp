#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <string>


// Cube mesh generation
std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::CubeMesh(glm::vec3 color)
{
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



std::vector<Vertex> Mesh::Readfile(const char* fileName, glm::vec3 color)
{
    std::ifstream inputFile(fileName);
    std::vector<Vertex> pointCloud;
    if (inputFile.is_open()) {

        std::string line;
        std::getline(inputFile, line);
        Vertex vertex;
        char comma; // to capture the commas in the file
        Vertex point;
		int skip = 0;
        float prevY = 0;

        while (inputFile >> point.x >> comma >> point.z >> comma >> point.y) {
                // Set the color and add the point to the cloud
            if (skip == 100)
            {
                point.x -= 608016.02;
                point.y -= 336.8007;
                point.z -= 6750620.771;

                if (point.y > prevY)
                {
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
                prevY = point.y;
                pointCloud.push_back(point);
				skip = 1;
            }
            else
                ++skip;
        }

        inputFile.close();

    }
    else {
        std::cerr << "Unable to open the input file for reading." << std::endl;
    }
    std::cout << "point Cloud " << pointCloud.size();
    return pointCloud;
}

