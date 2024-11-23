#include "RigidBody.h"
#include "../Component/Component.h"
#include "../Entity.h"
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <glm/glm.hpp>

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
        if (entity->GetEntityID() <= 1) continue;

        applyGravity(entity, deltaTime);
        BarycentricCoordinates(entity, entities[1], deltaTime);

        if (auto transform = entity->GetComponent<TransformComponent>()) {
            transform->position += transform->velocity * deltaTime;
            transform->rotation += transform->angularVelocity * deltaTime;
        }
    }
}

glm::vec3 RigidBody::CalculateNormalForce(std::shared_ptr<Entity> entity, std::shared_ptr<Entity> planeEntity, glm::vec3 contactNormal, float dt)
{
    glm::vec3 normalForce(0.0f);
    float gravitationalForce = mass * gravity;

    if (auto transform = entity->GetComponent<TransformComponent>()) {
        if (transform->position.y <= 0) { // Adjust this to use contact height if necessary
            float impactVelocity = std::abs(transform->velocity.y);

            if (impactVelocity > bounceThreshold) {
                // Calculate impact force
                float impactForce = mass * impactVelocity * dampingFactor;
                normalForce = gravitationalForce * contactNormal +
                    impactForce * glm::normalize(contactNormal);

                glm::vec3 pointOfImpact = transform->position + glm::vec3(0, 1, 0); // Offset for center of mass
                applyAngularForce(glm::vec3(0, impactForce, 0), pointOfImpact, entity);

                // Invert and dampen vertical velocity
                transform->velocity.y *= -dampingFactor;
            }
            else if (impactVelocity > stopBounceThreshold) {
                // Minimal bounce, just counteract gravity
                normalForce = gravitationalForce * contactNormal;
            }
            else {
                // Stop all motion along the normal direction
                normalForce = gravitationalForce * contactNormal;
                transform->velocity.y = 0.0f;
            }
             
            // Apply angular damping
            transform->angularVelocity *= angularDampingFactor;
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
    if (!transform || !planeTransform || !planeMesh) {
        return; // Ensuring all components exist
    }

    glm::vec3 point = transform->position;
    glm::vec3 ballSize = transform->scale; // Assuming scale gives the size of the entity
    std::vector<Vertex>& planeVertices = planeMesh->vertices;
    float groundThreshold = ballSize.y * 0.5; // or an even smaller value


    if (planeVertices.empty()) {
        return; // Return early if no vertices are present in the plane
    }

    for (int i = 0; i < planeMesh->indices.size(); i += 3) {
        int index0 = planeMesh->indices[i];
        int index1 = planeMesh->indices[i + 1];
        int index2 = planeMesh->indices[i + 2];

        // Transform the plane's vertices into world space
        glm::vec3 v0 = glm::vec3(
            (planeVertices[index0].x * planeTransform->scale.x) + planeTransform->position.x,
            (planeVertices[index0].y * planeTransform->scale.y) + planeTransform->position.y,
            (planeVertices[index0].z * planeTransform->scale.z) + planeTransform->position.z);

        glm::vec3 v1 = glm::vec3(
            (planeVertices[index1].x * planeTransform->scale.x) + planeTransform->position.x,
            (planeVertices[index1].y * planeTransform->scale.y) + planeTransform->position.y,
            (planeVertices[index1].z * planeTransform->scale.z) + planeTransform->position.z);

        glm::vec3 v2 = glm::vec3(
            (planeVertices[index2].x * planeTransform->scale.x) + planeTransform->position.x,
            (planeVertices[index2].y * planeTransform->scale.y) + planeTransform->position.y,
            (planeVertices[index2].z * planeTransform->scale.z) + planeTransform->position.z);

        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0v2 = v2 - v0;
        glm::vec3 v0p = point - v0;

        // Compute dot products for barycentric coordinates
        double dot00 = glm::dot(v0v1, v0v1);
        double dot01 = glm::dot(v0v1, v0v2);
        double dot02 = glm::dot(v0v1, v0p);
        double dot11 = glm::dot(v0v2, v0v2);
        double dot12 = glm::dot(v0v2, v0p);

        // Compute barycentric coordinates
        double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        double v = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double w = (dot00 * dot12 - dot01 * dot02) * invDenom;
        double u = 1 - v - w;

        // Check if the point is within the triangle
        if (u >= 0 && v >= 0 && w >= 0) {
            float height = v0.y * u + v1.y * v + v2.y * w;

            if (transform->position.y < height + groundThreshold) {
                if (transform->position.y < height - groundThreshold) {
                    transform->position.y = height + groundThreshold;
                }

                // Prevent the entity from sinking into the plane
                applyForce(glm::vec3(0, gravity, 0), dt, entity);

                // Calculate the normal of the triangle
                glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));
                if (glm::length(normal) < 1e-6f) return; // Skip degenerate triangles

                // Calculate the normal force
                glm::vec3 normalForce = CalculateNormalForce(entity, planeEntity, normal, dt);

                // Adjust the normal force for steep slopes
                glm::vec3 gravityForce(0.0f, -gravity, 0.0f);
                glm::vec3 normalComponent = glm::dot(gravityForce, normal) * normal;
                glm::vec3 gravityAlongSlope = gravityForce - normalComponent;

                // Scale gravity along the slope to enhance uphill/downhill behavior
                const float gravitySlopeScale = 1.2f; // Experiment with this value
                gravityAlongSlope *= gravitySlopeScale;

                // Apply forces
                transform->velocity += (normalForce / mass) * dt; // Apply normal force
                transform->velocity += gravityAlongSlope * dt;    // Apply gravity along the slope

                // Apply friction
                float avgFriction = (planeVertices[index0].friction + planeVertices[index1].friction + planeVertices[index2].friction) / 3.0f;
                glm::vec3 horizontalVelocity = glm::vec3(transform->velocity.x, 0.0f, transform->velocity.z);
                float velocityMagnitude = glm::length(horizontalVelocity);

                if (velocityMagnitude > 0.01f && avgFriction > 0.0f) {
                    glm::vec3 frictionForce = -glm::normalize(horizontalVelocity) * (avgFriction * glm::length(normalForce));

                    // Scale friction to prevent over-damping
                    const float frictionScale = 0.8f; // Experiment with this value
                    frictionForce *= frictionScale;

                    transform->velocity += frictionForce * dt;

                    // Stop very slow movement
                    if (glm::length(transform->velocity) < 0.01f) {
                        transform->velocity = glm::vec3(0.0f);
                    }
                }
            }

            return;
        }

    }
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
