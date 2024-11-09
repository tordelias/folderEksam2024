#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>

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
