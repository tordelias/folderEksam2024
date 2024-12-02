#pragma once
#include <glm/glm.hpp>
#include <memory>
#include <vector>

class Entity;
class Grid;

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
    void Update(std::vector<std::shared_ptr<Entity>> entities, std::shared_ptr<Grid> grid, float deltaTime); // Added deltaTime parameter
	void AddIndicesToCell(std::shared_ptr<Grid> Grid, std::shared_ptr<Entity> ground);

private:
	void UpdateCollision(std::shared_ptr<Grid> grid, float dt);
	void CheckCollision(std::shared_ptr<Entity>& object, std::vector<std::shared_ptr<Entity>>& objectToCheck, int startingIndex, float dt);
    void  SphereCollison(std::shared_ptr<Entity>& objA, std::shared_ptr<Entity>& objB, float DeltaTime);
	void ObjectCollisionResponse(std::shared_ptr<Entity>& objA, std::shared_ptr<Entity>& objB);

	void BarycentricCoordinates(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, std::shared_ptr<Grid> grid, float dt);
    glm::vec3 CalculateGravity(float inclineAngle, glm::vec3 slopeVector, glm::vec3 normal, float frictionCoefficient);


    float gravity = 9.81f;
    float mass = 1.0f;
    const float bounceThreshold = 2.f;       // Minimum velocity to trigger bounce
    const float stopBounceThreshold = 0.2f;   // Threshold below which we stop bouncing
	const float dampingFactor = 0.4f;		// Damping factor for bounce
	const float angularDampingFactor = 0.7f; // Damping factor for angular velocity
};
