#pragma once
#include <glm/glm.hpp>
#include "Mesh.h"
#include "../System/ParticleSystem.h"

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

struct ParticleComponent : public Component
{
public:
	glm::vec3 position;
	glm::vec3 velocity;
	glm::vec3 acceleration;
	glm::vec3 Size;
	ParticleComponent(const glm::vec3& pos = glm::vec3(0.0f),
		const glm::vec3& vel = glm::vec3(0.f),
		const glm::vec3& acc = glm::vec3(0.f, 9.81f, 0.f),
		const glm::vec3& size = glm::vec3(1.f))
		: position(pos), velocity(vel), acceleration(acc), Size(size)
	{
		particleSystem = std::make_shared<ParticleSystem>(position, acceleration, Size, 500);
	}
	std::shared_ptr<ParticleSystem> particleSystem;
};