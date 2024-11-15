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
    glm::vec3 CalculateNormalForce(std::shared_ptr<Entity> entity);

private:
    float gravity = 9.81f;
    float mass = 1.0f;
    const float bounceThreshold = 0.5f;       // Minimum velocity to trigger bounce
    const float stopBounceThreshold = 0.3f;   // Threshold below which we stop bouncing
    const float dampingFactor = 0.5f;
};
