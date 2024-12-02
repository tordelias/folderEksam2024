#pragma once
#include <glm/glm.hpp>
#include "Mesh.h"

class Component
{
public:
    virtual ~Component() = default;
};

struct TransformComponent : public Component
{
    glm::vec3 position;
    glm::vec3 rotation;
	glm::vec3 velocity;
    glm::vec3 angularVelocity; 
    glm::vec3 scale;

    TransformComponent(const glm::vec3& pos = glm::vec3(0.0f),
        const glm::vec3& rot = glm::vec3(0.0f),
        const glm::vec3& scl = glm::vec3(1.0f))
        : position(pos), rotation(rot), scale(scl), velocity(glm::vec3(0, 0, 0)), angularVelocity(glm::vec3(0, 0, 0)) {}
};

struct MeshComponent : public Component {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    const char* TexturePath;

    MeshComponent(const char* figure = "", const glm::vec3& color = glm::vec3(1,1,1), const char* texture = "");
};
class Entity;
struct SplineComponent : public Component
{
public:

	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<glm::vec3> controllpoints;
    std::shared_ptr<Entity> owner; 
    std::shared_ptr<Entity> spline;
    glm::vec3 startpos; 
	SplineComponent(std::shared_ptr<Entity> ownerptr);
    void CalculateBSpline();
};
