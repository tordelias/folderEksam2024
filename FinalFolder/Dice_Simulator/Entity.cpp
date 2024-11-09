#pragma once
#include "Entity.h"
#define GLM_ENABLE_EXPERIMENTAL
#include "Component/Component.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include "Resources/Shaders/EBO.h"
#include "Resources/Shaders/VBO.h"
#include "Resources/Shaders/VAO.h";
#include <iostream>


Entity::Entity() : vao(std::make_shared<VAO>()), vbo(std::make_shared<VBO>()), ebo(std::make_shared<EBO>()), EntityID(0), components()
{
}

Entity::~Entity()
{
}

void Entity::Destroy()
{
}

bool Entity::IsActive() const
{
    return false;
}

void Entity::ListAllComponents() const
{
}

void Entity::initalize()
{
    MeshComponent* meshComponent = GetComponent<MeshComponent>();
    if (!meshComponent) {
        std::cerr << "Entity does not have a MeshComponent!" << std::endl;
        return;
    }

    // Binding the VAO
    vao->Bind();
    vbo->Bind();

    glBufferData(GL_ARRAY_BUFFER, meshComponent->vertices.size() * sizeof(Vertex), meshComponent->vertices.data(), GL_STATIC_DRAW);

    vao->LinkAttrib(*vbo, 0, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, x)); // Position
    vao->LinkAttrib(*vbo, 1, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, r)); // Color
    vao->LinkAttrib(*vbo, 2, 2, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, u)); // TexCoords

    ebo->Bind();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, meshComponent->indices.size() * sizeof(unsigned int), meshComponent->indices.data(), GL_STATIC_DRAW);

    // Unbinding VAO, VBO, EBO
    vao->Unbind();
    vbo->Unbind();
    ebo->Unbind();

    std::cout << "Entity initialized" << std::endl;
}

void Entity::render()
{
    auto transform = GetComponent<TransformComponent>();
    auto meshComponent = GetComponent<MeshComponent>();

    glm::mat4 model = glm::mat4(1.0f);
    glm::quat quaternion = glm::quat(transform->rotation);
    glm::mat4 rotationMatrix = glm::toMat4(quaternion);

    model = glm::translate(model, transform->position);
    model = glm::scale(model, transform->scale);
    model *= rotationMatrix;

    vao->Bind();
    vbo->Bind();
    ebo->Bind();

    glDrawElements(GL_TRIANGLES, meshComponent->indices.size(), GL_UNSIGNED_INT, 0);

    vao->Unbind();
    vbo->Unbind();
    ebo->Unbind();
}
