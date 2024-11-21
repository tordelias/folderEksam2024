#pragma once
#include <vector>
#include <glm/glm.hpp>

struct Vertex
{
    float x, y, z;       // Position
    float r, g, b;       // Color
    float u, v;          // Texture coordinates
    float normalx, normaly, normalz; // Normals
    unsigned int index;
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
private:
    std::vector<Vertex> Readfile(const char* fileName, glm::vec3 color);
	float BSplineBasis(int i, int p, float t, const std::vector<float>& knots);
    std::vector<glm::vec3> BarycentricCoordinates(std::vector<unsigned int> indices, std::vector<Vertex> vertices);
};
