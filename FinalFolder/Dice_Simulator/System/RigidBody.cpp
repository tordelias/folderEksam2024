#include "RigidBody.h"
#include "../Component/Component.h"
#include "../Entity.h"
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "../System/Grid.h"

RigidBody::RigidBody()
{
}

RigidBody::~RigidBody()
{
}

void RigidBody::applyForce(glm::vec3 force, std::shared_ptr<Entity> entity)
{
    if (auto transform = entity->GetComponent<TransformComponent>()) {
        transform->velocity += force; // Apply linear force
    }
}

void RigidBody::applyForce(glm::vec3 force, float deltaTime, std::shared_ptr<Entity> entity)
{
    if (auto transform = entity->GetComponent<TransformComponent>()) {
        transform->velocity += force * deltaTime; // Apply force scaled by deltaTime
    }
}

void RigidBody::applyAngularForce(glm::vec3 force, glm::vec3 pointOfImpact, std::shared_ptr<Entity> entity)
{
    if (auto transform = entity->GetComponent<TransformComponent>()) {
        glm::vec3 r = pointOfImpact - transform->position;
        glm::vec3 torque = glm::cross(r, force);
        float momentOfInertia = 1.0f; // Assuming unit moment of inertia
        glm::vec3 angularAcceleration = torque / momentOfInertia;
        transform->angularVelocity += angularAcceleration;
    }
}

void RigidBody::applyGravity(std::shared_ptr<Entity> entity, float deltaTime)
{
    applyForce(glm::vec3(0, -mass * gravity, 0), entity);
}

void RigidBody::applyRandomForce(std::vector<std::shared_ptr<Entity>> entities)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(10.0f, 30.0f);

    for (auto& entity : entities) {
        glm::vec3 randomForce(0, dis(gen), 0);
        glm::vec3 randomAngularForce(dis(gen), dis(gen), dis(gen));

        applyForce(randomForce, entity);

        if (auto transform = entity->GetComponent<TransformComponent>()) {
            glm::vec3 pointOfImpact = transform->position + glm::vec3(0, 0, 1);
            applyAngularForce(randomAngularForce, pointOfImpact, entity);
        }
    }
}

// Calculate the gravity force for the entity
void RigidBody::Update(std::vector<std::shared_ptr<Entity>> entities, std::shared_ptr<Grid> grid, float deltaTime)
{
	UpdateCollision(grid, deltaTime);
    for (auto& entity : entities) {
        if (entity->GetEntityID() == 0) continue;

        BarycentricCoordinates(entity, entities[0], grid, deltaTime);

        if (auto transform = entity->GetComponent<TransformComponent>()) {
            transform->position += transform->velocity * deltaTime;
            transform->rotation += transform->angularVelocity * deltaTime;
        }
    }
}

// Add the indices of the ground to the grid cells
void RigidBody::AddIndicesToCell(std::shared_ptr<Grid> grid, std::shared_ptr<Entity> ground) {
    if (auto groundMesh = ground->GetComponent<MeshComponent>()) {
        for (int i = 0; i < groundMesh->indices.size(); i += 3) {
            int index0 = groundMesh->indices[i];
            int index1 = groundMesh->indices[i + 1];
            int index2 = groundMesh->indices[i + 2];

            // Apply the transformation to the vertices in world space
            glm::mat4 transformation = glm::translate(glm::mat4(1.0f), ground->GetComponent<TransformComponent>()->position) *
                glm::scale(glm::mat4(1.0f), ground->GetComponent<TransformComponent>()->scale);

            glm::vec3 v0 = glm::vec3(transformation * glm::vec4(groundMesh->vertices[index0].x, groundMesh->vertices[index0].y, groundMesh->vertices[index0].z, 1.f));
            glm::vec3 v1 = glm::vec3(transformation * glm::vec4(groundMesh->vertices[index1].x, groundMesh->vertices[index1].y, groundMesh->vertices[index1].z, 1.f));
            glm::vec3 v2 = glm::vec3(transformation * glm::vec4(groundMesh->vertices[index2].x, groundMesh->vertices[index2].y, groundMesh->vertices[index2].z, 1.f));

            // Calculate the bounds of the triangle
            glm::vec3 minBounds = glm::min(glm::min(v0, v1), v2);
            glm::vec3 maxBounds = glm::max(glm::max(v0, v1), v2);

            // Determine grid cell bounds for the triangle without offsetting twice
            int minX = static_cast<int>(std::floor((minBounds.x) / grid->m_cellSize));
            int minY = static_cast<int>(std::floor((minBounds.z) / grid->m_cellSize));
            int maxX = static_cast<int>(std::floor((maxBounds.x) / grid->m_cellSize));
            int maxY = static_cast<int>(std::floor((maxBounds.z) / grid->m_cellSize));

            // Add the triangle indices to all intersecting grid cells
            for (int x = minX; x <= maxX; ++x) {
                for (int y = minY; y <= maxY; ++y) {
                    // Make sure to check if the cell is valid
                    Cell* cell = grid->getCell(x, y);
                    if (cell) {
                        // Add the indices to the cell's groundIndices
                        cell->groundIndices.push_back(index0);
                        cell->groundIndices.push_back(index1);
                        cell->groundIndices.push_back(index2);
                    }
                }
            }
        }
    }
}


// Update the collision detection for the rigid body
void RigidBody::UpdateCollision(std::shared_ptr<Grid> grid, float dt)
{
    for (int i = 0; i < grid->m_cells.size(); ++i)
    {
        int x = i % grid->m_numXCells;
        int y = i / grid->m_numXCells;

        Cell& cell = grid->m_cells[i];

        for (int j = 0; j < cell.entities.size(); ++j)
        {
            std::shared_ptr<Entity> object = cell.entities[j];

            // Check for nullptr on the object
            if (!object)
                continue;

            CheckCollision(object, cell.entities, j + 1, dt);

            // Check for neighbors safely
            if (x > 0)
            {
                Cell* leftCell = grid->getCell(x - 1, y);
                if (leftCell) {
                    CheckCollision(object, leftCell->entities, 0, dt);
                }

                if (y > 0)
                {
                    Cell* bottomLeftCell = grid->getCell(x - 1, y - 1);
                    if (bottomLeftCell) {
                        CheckCollision(object, bottomLeftCell->entities, 0, dt);
                    }
                }

                if (y < grid->m_numYCells - 1)
                {
                    Cell* topLeftCell = grid->getCell(x - 1, y + 1);
                    if (topLeftCell) {
                        CheckCollision(object, topLeftCell->entities, 0, dt);
                    }
                }
            }

            if (y > 0)
            {
                Cell* bottomCell = grid->getCell(x, y - 1);
                if (bottomCell) {
                    CheckCollision(object, bottomCell->entities, 0, dt);
                }
            }
        }
    }
}

// Check for collision between objects
void RigidBody::CheckCollision(std::shared_ptr<Entity>& object, std::vector<std::shared_ptr<Entity>>& objectToCheck, int startingIndex, float dt)
{
    for (int i = startingIndex; i < objectToCheck.size(); ++i)
    {
        if (objectToCheck[i]) {  // Ensure the object is not nullptr
            SphereCollison(object, objectToCheck[i], dt);
        }
        else {
            std::cerr << "Error: objectToCheck[" << i << "] is nullptr!" << std::endl;
        }
    }
}
// Check for collision between two spheres
void RigidBody::SphereCollison(std::shared_ptr<Entity>& objA, std::shared_ptr<Entity>& objB, float DeltaTime)
{
    // Ensure both entities have TransformComponents before accessing
    auto transformA = objA->GetComponent<TransformComponent>();
    auto transformB = objB->GetComponent<TransformComponent>();

    if (!transformA) {
        std::cerr << "Error: objA does not have TransformComponent!" << std::endl;
        return;
    }
    if (!transformB) {
        std::cerr << "Error: objB does not have TransformComponent!" << std::endl;
        return;
    }

    glm::vec3 posA = transformA->position;
    glm::vec3 posB = transformB->position;
    float distance_centers = glm::length(posA - posB);
	// Check if the distance between the centers of the two spheres is less than the sum of their radius
    if (distance_centers <= (transformA->scale.x + transformB->scale.x)) {
        float minimuntranslation = transformA->scale.x + transformB->scale.x - distance_centers;
        auto dirvec = glm::normalize(posA - posB);
        transformA->position = (transformA->position + dirvec * minimuntranslation);
        ObjectCollisionResponse(objA, objB);
    }
}
// Calculate the response to the collision
void RigidBody::ObjectCollisionResponse(std::shared_ptr<Entity>& objA, std::shared_ptr<Entity>& objB)
{
    auto transformA = objA->GetComponent<TransformComponent>();
    auto transformB = objB->GetComponent<TransformComponent>();

    if (!transformA) {
        std::cerr << "Error: objA does not have TransformComponent!" << std::endl;
        return;
    }
    if (!transformB) {
        std::cerr << "Error: objB does not have TransformComponent!" << std::endl;
        return;
    }

    float massA = mass;  // Use the mass of the RigidBody
    float massB = mass;

    glm::vec3 posA = transformA->position;
    glm::vec3 posB = transformB->position;
    glm::vec3 velocityA = transformA->velocity;
    glm::vec3 velocityB = transformB->velocity;

	// Calculate the normal vector
    glm::vec3 normal = glm::normalize(posB - posA);
	// Calculate the relative velocity
    glm::vec3 relativeVelocity = velocityA - velocityB;
	// Calculate the relative velocity along the normal
    float velocityAlongNormal = glm::dot(relativeVelocity, normal);

	float restitution = 1.00f;  // 1 == Perfectly elastic collision, 0 == Perfectly inelastic collision
    float impulse = (-(1 + restitution) * velocityAlongNormal) / (1 / massA + 1 / massB);

	// Apply the impulse to the objects
    glm::vec3 impulseVector = impulse * normal;

	// Update the velocities
    glm::vec3 newVelocityA = velocityA + (impulseVector / massA);
    glm::vec3 newVelocityB = velocityB - (impulseVector / massB);

    transformA->velocity = newVelocityA;
    transformB->velocity = newVelocityB;
}

// Calculate the gravity force
void RigidBody::BarycentricCoordinates(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, std::shared_ptr<Grid> grid, float dt)
{
    auto transform = entity->GetComponent<TransformComponent>();
    auto planeTransform = planeEntity->GetComponent<TransformComponent>();
    auto planeMesh = planeEntity->GetComponent<MeshComponent>();

    if (!transform) {
        std::cerr << "Error: entity does not have TransformComponent!" << std::endl;
        return;
    }
    if (!planeTransform) {
        std::cerr << "Error: planeEntity does not have TransformComponent!" << std::endl;
        return;
    }
    if (!planeMesh) {
        std::cerr << "Error: planeEntity does not have MeshComponent!" << std::endl;
        return;
    }

    auto gridCell = grid->getCell(transform->position);
    if (!gridCell) {
        std::cerr << "Error: Grid cell not found for position: " << transform->position.x << ", " << transform->position.y << ", " << transform->position.z << std::endl;
        applyGravity(entity, dt);
        return;
    }

    auto indices = gridCell->groundIndices;
    if (indices.empty()) {
        std::cerr << "Warning: No groundIndices found in grid cell!" << std::endl;
        applyGravity(entity, dt);
        return;
    }

    const float epsilon = 1e-4f;
    glm::vec3 point = transform->position;
    glm::vec3 ballSize = transform->scale;
    std::vector<Vertex>& planeVertices = planeMesh->vertices;
    float groundThreshold = ballSize.y;

    if (planeVertices.empty()) {
        std::cerr << "Error: Plane mesh vertices are empty!" << std::endl;
        return;
    }

    for (int i = 0; i < indices.size(); i += 3) {
        if (i + 2 >= indices.size()) {
            std::cerr << "Error: Index out of bounds in groundIndices: " << i << ", " << i + 1 << ", " << i + 2 << std::endl;
            continue;
        }

        int index0 = indices[i];
        int index1 = indices[i + 1];
        int index2 = indices[i + 2];

        if (index0 >= planeVertices.size() || index1 >= planeVertices.size() || index2 >= planeVertices.size()) {
            std::cerr << "Error: Index out of bounds in planeVertices: " << index0 << ", " << index1 << ", " << index2 << std::endl;
            continue;
        }

        glm::mat4 transformation = glm::translate(glm::mat4(1.0f), planeTransform->position) *
            glm::scale(glm::mat4(1.0f), planeTransform->scale);

        glm::vec3 v0 = glm::vec3(transformation * glm::vec4(planeVertices[index0].x, planeVertices[index0].y, planeVertices[index0].z, 1.0f));
        glm::vec3 v1 = glm::vec3(transformation * glm::vec4(planeVertices[index1].x, planeVertices[index1].y, planeVertices[index1].z, 1.0f));
        glm::vec3 v2 = glm::vec3(transformation * glm::vec4(planeVertices[index2].x, planeVertices[index2].y, planeVertices[index2].z, 1.0f));

        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0v2 = v2 - v0;
        glm::vec3 v0p = point - v0;

        double dot00 = glm::dot(v0v1, v0v1);
        double dot01 = glm::dot(v0v1, v0v2);
        double dot02 = glm::dot(v0v1, v0p);
        double dot11 = glm::dot(v0v2, v0v2);
        double dot12 = glm::dot(v0v2, v0p);

        double invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
        double v = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double w = (dot00 * dot12 - dot01 * dot02) * invDenom;
        double u = 1.0f - v - w;

        if (u < -epsilon || v < -epsilon || w < -epsilon) {
            continue;
        }

        if (glm::length(glm::cross(v0v1, v0v2)) < epsilon) {
            continue;
        }

        if (u >= 0 && v >= 0 && w >= 0) {
            float height = v0.y * u + v1.y * v + v2.y * w;

            glm::vec3 currentVelocity = transform->velocity;
            if (transform->position.y < height + groundThreshold) {
                // Stopping downward motion and applying corrective force if sinking
                if (currentVelocity.y < 0) {
                    transform->velocity.y *= -0.1f;
                }

                transform->velocity = currentVelocity;

                // Correcting position if sinking into the ground
                if (transform->position.y < height + groundThreshold) {
                    transform->position.y = height + groundThreshold;
					currentVelocity.y = 0.0f;
                }

                glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
                if (glm::length(normal) < epsilon) continue; // Degenerate triangle check


                glm::vec3 slopeVector = glm::normalize(glm::vec3(normal.x, 0.0f, normal.z)); // Slope direction

				// Calculate the average friction coefficient of the triangle
                float friction0 = planeVertices[index0].friction;
                float friction1 = planeVertices[index1].friction;
                float friction2 = planeVertices[index2].friction;
                float frictionCoefficient = (friction0 + friction1 + friction2) / 3;

				//calculates the acceleration due to gravity along the slope
                glm::vec3 gravityAlongSlope = CalculateGravity(0.f, slopeVector, normal, frictionCoefficient);

                // Apply both gravity and friction if there is velocity
                if (glm::length(currentVelocity) > 0.0f) {
                    glm::vec3 velocityDirection = glm::normalize(currentVelocity);

                    // Calculate friction force: opposite direction to velocity
                    glm::vec3 frictionForce = -frictionCoefficient * glm::length(normal) * velocityDirection;

                    // Cap the friction force to not exceed the current velocity's magnitude
                    if (glm::length(frictionForce) > glm::length(currentVelocity)) 
                    {
                        frictionForce = -currentVelocity; 
                    }
                    // Update velocity considering both gravity and friction
                    glm::vec3 newVelocity = currentVelocity + gravityAlongSlope + frictionForce;

                    // Stop the object if its velocity becomes very small
                    if (glm::length(newVelocity) < 0.001f) {
                        newVelocity = glm::vec3(0.0f);
                    }
                    transform->velocity = newVelocity;
                }
                else {
                    // If the velocity is zero, only gravity acts
                    transform->velocity += gravityAlongSlope * dt;
                }

                return;
            }

            // If the object isn't on the ground, apply gravity normally
            applyGravity(entity, dt);



            return;
        }
    }
    applyGravity(entity, dt);
}





glm::vec3 RigidBody::CalculateGravity(float inclineAngle, glm::vec3 slopeVector, glm::vec3 normal, float frictionCoefficient)
{
    slopeVector = glm::normalize(slopeVector);
    glm::vec3 gravityForce(0.0f, -gravity, 0.0f);
	//project the forces on the normal and slope vectors
    float normalForceMagnitude = glm::dot(gravityForce, normal); 
    glm::vec3 normalForce = normal * normalForceMagnitude;
    //Gx = G - N 
    glm::vec3 gravityParallel = gravityForce - normalForce;
	//project the forces on the slope vector
    glm::vec3 gravityAlongSlope = glm::dot(gravityParallel, slopeVector) * slopeVector;


    glm::vec3 frictionForce = -frictionCoefficient * glm::normalize(gravityAlongSlope) * glm::length(normalForce);
    if (glm::length(frictionForce) > glm::length(gravityAlongSlope)) {
        frictionForce = glm::vec3(0.0f); 
        gravityAlongSlope = glm::vec3(0.0f);
    }

    gravityAlongSlope += frictionForce;

    // Applying the force along the slope
    return gravityAlongSlope;
}


