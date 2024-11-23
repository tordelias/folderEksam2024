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
	float BSplineBasis(int i, int p, float t, const std::vector<float>& knots);
    std::vector<glm::vec3> BarycentricCoordinates(std::vector<unsigned int> indices, std::vector<Vertex> vertices);

	std::vector<Vertex> m_vertices;
	std::vector<unsigned int> m_indices;
};
