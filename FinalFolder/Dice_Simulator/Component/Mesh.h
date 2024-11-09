#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Vertex
{
    float x, y, z;       // Position
    float r, g, b;       // Color
    float u, v;          // Texture coordinates
    float normalx, normaly, normalz; // Normals
};

class Mesh
{
public:
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> CubeMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> SphereMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> CylinderMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> ConeMesh(glm::vec3 color);
    std::pair<std::vector<Vertex>, std::vector<unsigned int>> TorusMesh(glm::vec3 color);
};
