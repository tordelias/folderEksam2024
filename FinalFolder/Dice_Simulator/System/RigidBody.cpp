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
    applyForce(glm::vec3(0, -mass * gravity, 0) * deltaTime, entity);
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

void RigidBody::Update(std::vector<std::shared_ptr<Entity>> entities, std::shared_ptr<Grid> grid, float deltaTime)
{
    for (auto& entity : entities) {
        if (entity->GetEntityID() == 0) continue;

        BarycentricCoordinates(entity, entities[0], grid, deltaTime);

        if (auto transform = entity->GetComponent<TransformComponent>()) {
            transform->position += transform->velocity * deltaTime;
            transform->rotation += transform->angularVelocity * deltaTime;
        }
    }
}

glm::vec3 RigidBody::CalculateNormalForce(std::shared_ptr<Entity> entity, double height, float dt) {
    glm::vec3 normalForce = glm::vec3(0, 0, 0);
    float gravitationalForce = mass * gravity;

    if (auto transform = entity->GetComponent<TransformComponent>()) {
        // Check if the object is at or below ground level (y <= 0)
        if (transform->position.y <= height) {
            // Clamp the object's position to prevent it from going below ground level
            transform->position.y = height + 0.5;

            float impactVelocity = abs(transform->velocity.y);

            if (impactVelocity > bounceThreshold) {
                // Exponential sinusoidal bounce
                float impactForce = mass * impactVelocity * dampingFactor; // Reduce force using damping factor
                normalForce.y = gravitationalForce + impactForce * sin(impactVelocity);

                // Apply torque for rotation at point of impact
                glm::vec3 pointOfImpact = transform->position + glm::vec3(0, 1, 0);
                applyAngularForce(glm::vec3(0, impactForce, 0), pointOfImpact, entity);

                // Reduce vertical velocity for the next bounce using damping
                transform->velocity.y *= -dampingFactor;
            }
            else if (impactVelocity > stopBounceThreshold) {
                // Apply only gravity to counter minor bounces without further sinusoidal effect
                normalForce.y = gravitationalForce;
            }
            else {
                // When velocity is minimal, stop bouncing entirely by setting velocity to zero
                normalForce.y = gravitationalForce;
                transform->velocity.y = 0.0f; // Stop small residual vertical motion
            }

            // Apply angular damping to reduce rotation gradually
            float angularDampingFactor = 0.95f;
            transform->angularVelocity *= angularDampingFactor;

            // Check if angular velocity is effectively zero
            if (glm::length(transform->angularVelocity) < 0.01f) {
                transform->angularVelocity = glm::vec3(0.0f);
            }
        }
    }

    return normalForce;
}

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







void RigidBody::BarycentricCoordinates(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, std::shared_ptr<Grid> grid, float dt) {
    auto transform = entity->GetComponent<TransformComponent>();
    auto planeTransform = planeEntity->GetComponent<TransformComponent>();
    auto planeMesh = planeEntity->GetComponent<MeshComponent>();
    auto indices = grid->getCell(transform->position)->groundIndices;

    if (!transform || !planeTransform || !planeMesh) return;

    if (indices.empty()) {
        //applyGravity(entity, dt);
        return;
    }

    const float epsilon = 1e-4f;

    glm::vec3 point = transform->position;
    glm::vec3 ballSize = transform->scale;
    std::vector<Vertex>& planeVertices = planeMesh->vertices;
    float groundThreshold = ballSize.y;

    if (planeVertices.empty()) return;

    for (int i = 0; i < indices.size(); i += 3) {
        int index0 = indices[i];
        int index1 = indices[i + 1];
        int index2 = indices[i + 2];

        if (index0 >= planeVertices.size() || index1 >= planeVertices.size() || index2 >= planeVertices.size()) {
            std::cerr << "Index out of bounds!" << std::endl;
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

        if (u < -epsilon || v < -epsilon || w < -epsilon)
            continue;

        // Skip degenerate triangles
        if (glm::length(glm::cross(v0v1, v0v2)) < epsilon)
            continue;

        // If the point is inside the triangle (u, v, w > 0)
        if (u >= 0 && v >= 0 && w >= 0) {
            double height = v0.y * u + v1.y * v + v2.y * w;

            glm::vec3 currentVelocity = transform->velocity;
			glm::vec3 normal = glm::normalize(glm::normalize(glm::cross(v0v1, v0v2)));

            if (transform->position.y < height + groundThreshold) {
                //glm::vec3 velocityNormal = glm::dot(currentVelocity, normal) * normal;
                //glm::vec3 velocityTangent = currentVelocity - velocityNormal;

                //// Reflect velocityNormal if object is falling towards the surface
                //if (glm::dot(currentVelocity, normal) < 0) {
                //    velocityNormal = -velocityNormal * 0.5f; // Slight energy loss
                //}
                //float friction = 0.1f; // Adjust for realistic sliding
                //velocityTangent *= (1 - friction);
                //transform->velocity = velocityNormal + velocityTangent;
                transform->position.y = height + groundThreshold;

                float inclineAngle = std::acos(normal.y);
                glm::vec3 slopeVector = glm::normalize(glm::vec3(normal.x, 0.0f, normal.z)); // Slope direction
                float speedAdjustment = glm::dot(currentVelocity, slopeVector);
                if (currentVelocity.y > 0) { // Ball is moving upward
                    currentVelocity.y -= speedAdjustment * sin(inclineAngle);

                    // Ensuring ball doesn't go through the floor
                    if (transform->position.y < height + groundThreshold) {
                        transform->position.y = height + groundThreshold;
                        currentVelocity.y = 0; // Stopping upward motion
                    }
                }
                else if (currentVelocity.y < 0) { // Ball is moving downward
                    currentVelocity.y += speedAdjustment * sin(inclineAngle);

                    // Ensuring ball doesn't go through the floor
                    if (transform->position.y < height + groundThreshold)
                    {
                        transform->position.y = height + groundThreshold;

                        currentVelocity.y = 0; // Stopping downward motion
                    }
                }

                if (glm::abs(normal.y) < 1.0f) 
                {
                    glm::vec3 slopeVector = glm::normalize(glm::vec3(normal.x, 0.0f, normal.z));
                    glm::vec3 gravityAlongSlope = CalculateGravity(0.f, slopeVector, normal);
                    transform->velocity += gravityAlongSlope;
                    //applyForce(gravityAlongSlope, entity);
                }

                return;
            }

            applyGravity(entity, dt);
            return;
        }
    }

    applyGravity(entity, dt);
}



glm::vec3 RigidBody::CalculateGravity(float inclineAngle, glm::vec3 slopeVector, glm::vec3 normal)
{
    slopeVector = glm::normalize(slopeVector);


    glm::vec3 gravityForce(0.0f, -gravity, 0.0f);

    // Calculating normal force (perpendicular to the slope)
    float normalForceMagnitude = glm::dot(gravityForce, normal); // Gravity along the normal
    glm::vec3 normalForce = normal * normalForceMagnitude;

    // Calculating gravitational force acting parallel to the slope (slope vector)
    glm::vec3 gravityParallel = gravityForce - normalForce; // Parallel force along the slope

    // Projecting this parallel gravity onto the slope's horizontal direction (slopeVector)
    glm::vec3 gravityAlongSlope = glm::dot(gravityParallel, normal) * normal;
	//float angle1 = acos(normal.z / glm::length(normal));
	//float angle2 = atan(normal.x / normal.z);
 //   float ax = gravity * sin(angle1) * sin(angle2) * cos(angle1); 
	//float az = gravity * sin(angle1) * cos(angle2) * cos(angle1);
 //   float ay = gravity * ((cos(angle1) * cos(angle1)) - 1);
	//glm::vec3 gravityAlongSlope = glm::vec3(ax, ay, az);

    // Applying the force along the slope
    return gravityAlongSlope;

}
