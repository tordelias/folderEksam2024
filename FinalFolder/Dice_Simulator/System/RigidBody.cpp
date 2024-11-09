#include "RigidBody.h"
#include "../Component/Component.h"
#include <iostream>
#include <vector>
#include <memory>
#include <random>
#include "../Entity.h"

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
        transform->velocity += (force * deltaTime); // Apply force scaled by deltaTime
    }
}

void RigidBody::applyAngularForce(glm::vec3 force, glm::vec3 pointOfImpact, std::shared_ptr<Entity> entity)
{
    if (auto transform = entity->GetComponent<TransformComponent>()) {
        // Calculate the torque (r x F)
        glm::vec3 r = pointOfImpact - transform->position;
        glm::vec3 torque = glm::cross(r, force);

        // Assuming unit moment of inertia for simplicity
        float momentOfInertia = 1.0f; // Change this to match your object’s properties
        glm::vec3 angularAcceleration = torque / momentOfInertia; // realistic cases

        // Update angular velocity
        transform->angularVelocity += angularAcceleration;
    }
}

void RigidBody::applyGravity(std::shared_ptr<Entity> entity, float deltaTime)
{
    applyForce(glm::vec3(0, -mass * gravity * deltaTime, 0), entity);
}

void RigidBody::applyRandomForce(std::vector<std::shared_ptr<Entity>> entities)
{
    for (auto& entity : entities)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        //std::uniform_real_distribution<float> dis(50.0f, 100.0f); // Adjusted to ensure significant forces
        std::uniform_real_distribution<float> dis1(10.0f, 30.0f); // Adjusted to ensure significant forces

        glm::vec3 randomForce = glm::vec3(0, dis1(gen), 0); // Random upward force
        glm::vec3 randomAngularForce = glm::vec3(dis1(gen), dis1(gen), dis1(gen)); // Random angular force

        // Apply forces
        applyForce(randomForce, entity);

        // Calculate a point of impact slightly away from the center for torque effect
        auto transform = entity->GetComponent<TransformComponent>();
        glm::vec3 pointOfImpact = transform->position + glm::vec3(0, 0, 1); // Change this to a reasonable offset

        applyAngularForce(randomAngularForce, pointOfImpact, entity);
    }
}

void RigidBody::Update(std::vector<std::shared_ptr<Entity>> entities, float deltaTime)
{
    for (auto& entity : entities) {
        applyGravity(entity, deltaTime);
        normalForceGround(entity, deltaTime);

        auto transform = entity->GetComponent<TransformComponent>();

        if (transform) {
            // Update position and rotation based on velocity and deltaTime
            transform->position += transform->velocity * deltaTime; // Scale by deltaTime

            // Update rotation based on angular velocity
            transform->rotation += transform->angularVelocity * deltaTime; // Scale by deltaTime
        }
    }
}

void RigidBody::normalForceGround(std::shared_ptr<Entity> entity, float deltaTime)
{
    applyForce(CalculateNormalForce(entity), deltaTime, entity);
}


glm::vec3 RigidBody::CalculateNormalForce(std::shared_ptr<Entity> entity) {
    glm::vec3 normalForce = glm::vec3(0, 0, 0);
    float gravitationalForce = mass * gravity;

    // Define the bounce parameters
    const float bounceThreshold = 0.5f;  // Minimum velocity to trigger bounce
    const float stopBounceThreshold = 0.1f; // Threshold below which we stop bouncing
    const float dampingFactor = 0.8f;  // Exponential decay factor for each bounce (0 < damping < 1)

    if (auto transform = entity->GetComponent<TransformComponent>()) {
        // Check if the object is at or below ground level (y <= 0)
        if (transform->position.y <= 0) {
            // Clamp the object's position to prevent it from going below ground level
            transform->position.y = 0.0f;

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





