#include "RigidBody.h"
#include "../Component/Component.h"
#include "../Entity.h"
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

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

void RigidBody::Update(std::vector<std::shared_ptr<Entity>> entities, float deltaTime)
{
    for (auto& entity : entities) {
        if (entity->GetEntityID() == 0) continue;

        BarycentricCoordinates(entity, entities[0], deltaTime);

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


void RigidBody::BarycentricCoordinates(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, float dt) {
    auto transform = entity->GetComponent<TransformComponent>();
    auto planeTransform = planeEntity->GetComponent<TransformComponent>();
    auto planeMesh = planeEntity->GetComponent<MeshComponent>();
    if (!transform || !planeTransform || !planeMesh) return;

    glm::vec3 point = transform->position;
    glm::vec3 ballSize = transform->scale;
    std::vector<Vertex>& planeVertices = planeMesh->vertices;
    float groundThreshold = 0.1f;  // Tolerance for "ground" contact

    if (planeVertices.empty()) return;

    bool collided = false;  // Flag to track if a collision has occurred

    for (int i = 0; i < planeMesh->indices.size(); i += 3) {
        int index0 = planeMesh->indices[i];
        int index1 = planeMesh->indices[i + 1];
        int index2 = planeMesh->indices[i + 2];

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

        // Check if the barycentric coordinates are valid
        if (u < -0.001f || v < -0.001f || w < -0.001f) {
            continue;  // Skip invalid triangles
        }

        // Normalize the barycentric coordinates to ensure they sum to 1
        float sum = u + v + w;
        u /= sum;
        v /= sum;
        w /= sum;

        // Calculate the height at the point using barycentric interpolation
        float height = v0.y * u + v1.y * v + v2.y * w;

        // The point is considered on the ground if it is within the threshold
        if (u >= -0.001f && v >= -0.001f && w >= -0.001f) {
            glm::vec3 currentVelocity = transform->velocity;
            collided = true;  // Mark as collided
            float targetHeight = height + groundThreshold;

            // Only adjust the position if it's actually above the ground by more than the threshold
            if (transform->position.y > targetHeight) {
                transform->position.y = targetHeight;  // Snap it to the ground
            }

            // Prevent the object from "bouncing" back up
            if (transform->position.y < targetHeight) {
                transform->position.y = targetHeight;
            }

            // Apply friction if collision occurred
            if (collided) {
                glm::vec3 normal = glm::normalize(glm::cross(v0v1, v0v2));
                if (glm::length(normal) == 0.0f) continue;  // Skip degenerate triangles

                glm::vec3 velocityAlongSurface = currentVelocity - glm::dot(currentVelocity, normal) * normal;
                glm::vec3 slopeVector = glm::normalize(glm::vec3(normal.x, 0.0f, normal.z)); // Slope direction

                // Weighted friction coefficients
                float friction0 = planeVertices[index0].friction;
                float friction1 = planeVertices[index1].friction;
                float friction2 = planeVertices[index2].friction;
                float frictionCoefficient = u * friction0 + v * friction1 + w * friction2;

                // Apply friction force to horizontal velocity
                glm::vec3 frictionForce = -velocityAlongSurface * frictionCoefficient;
                transform->velocity += frictionForce * dt;

                // Apply gravity along the slope if the object is on an incline
                if (glm::length(slopeVector) > 0.00000001f) {
                    glm::vec3 gravityAlongSlope = CalculateGravity(v0v1, v0v2);
                    applyForce(gravityAlongSlope, dt, entity);
                }

                return;
            }

            applyGravity(entity, dt);
            return;  // Object is now on the ground, stop further gravity application
        }
    }

    // If no collision happened, continue applying gravity
    applyGravity(entity, dt);
}











glm::vec3 RigidBody::CalculateGravity(glm::vec3 v0, glm::vec3 v1)
{
    glm::vec3 normal = glm::normalize(glm::cross(v0, v1));
    glm::vec3 gravityForce(0.0f, -gravity, 0.0f);

    float normalForceMagnitude = glm::dot(gravityForce, normal);
    glm::vec3 normalForce = normal * normalForceMagnitude;

    glm::vec3 gravityParallel = gravityForce - normalForce;
    return gravityParallel;
}
