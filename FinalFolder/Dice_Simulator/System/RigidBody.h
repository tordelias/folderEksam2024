#pragma once
#include <glm/glm.hpp>
#include <memory>
#include <vector>

class Entity;

class RigidBody
{
public:
    RigidBody();
    ~RigidBody();

    void applyForce(glm::vec3 force, std::shared_ptr<Entity> entity);
    void applyForce(glm::vec3 force, float deltaTime, std::shared_ptr<Entity> entity);
    void applyAngularForce(glm::vec3 force, glm::vec3 pointOfImpact, std::shared_ptr<Entity> entity);
    void applyGravity(std::shared_ptr<Entity> entity, float deltaTime);
    void applyRandomForce(std::vector<std::shared_ptr<Entity>> entities);
    void Update(std::vector<std::shared_ptr<Entity>> entities, float deltaTime); // Added deltaTime parameter
    void normalForceGround(std::shared_ptr<Entity> entity, float deltaTime);
    glm::vec3 CalculateNormalForce(std::shared_ptr<Entity> entity, double height, float dt);

private:
	void BarycentricCoordinates(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, float dt);
    glm::vec3 CalculateGravity(glm::vec3 v0, glm::vec3 v1);
    float gravity = 9.81f;
    float mass = 1.0f;
    const float bounceThreshold = 2.f;       // Minimum velocity to trigger bounce
    const float stopBounceThreshold = 0.2f;   // Threshold below which we stop bouncing
	const float dampingFactor = 0.4f;		// Damping factor for bounce
	const float angularDampingFactor = 0.7f; // Damping factor for angular velocity
};
