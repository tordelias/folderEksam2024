#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include "EntityManager.h"
#include "../Resources/Shaders/shaderClass.h"
#include "../Component/Component.h"
#include "../Resources/Shaders/EBO.h"
#include "../Resources/Shaders/VBO.h"
#include "../Resources/Shaders/VAO.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>      
#include <glm/gtx/quaternion.hpp>     
#include "../Resources/Texture/Texture.h"
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <memory>
#include "../System/RigidBody.h"
#include "../System/Grid.h"



EntityManager::EntityManager(std::shared_ptr<Shader> shaderprogram) : EntityCount(0), shader(shaderprogram), rigidbody(std::make_shared<RigidBody>()), 
grid(std::make_shared<Grid>(10000, 10000, 50))
{
}

EntityManager::~EntityManager()
{
}

void EntityManager::Update()
{
	for (auto& entity : entities)
	{
        Cell* newCell = grid->getCell(entity->GetComponent<TransformComponent>()->position);
		if (newCell == nullptr)
		{
			continue;
		}
        if (newCell != entity->ownerCell)
        {
            grid->RemoveBallFromCell(entity, entity->ownerCell);
            grid->AddBaLL(entity, newCell);
            entity->ownerCell = newCell;
        }
	}
}


void EntityManager::Render(glm::mat4 viewproj, float dt) {
    Update();
    rigidbody->Update(entities, grid, dt);
    int textureCount = 0;

    for (auto& entity : entities) {
        auto meshComponent = entity->GetComponent<MeshComponent>();
        if (meshComponent == nullptr)
            continue;

        auto transform = entity->GetComponent<TransformComponent>();
        if (transform == nullptr) continue;

        glm::mat4 model = glm::mat4(1.0f);
        glm::quat quaternion = glm::quat(transform->rotation);
        glm::mat4 rotationMatrix = glm::toMat4(quaternion);

        model = glm::translate(model, transform->position);
        model = glm::scale(model, transform->scale);
        model *= rotationMatrix;

        glUniformMatrix4fv(glGetUniformLocation(shader->ID, "camMatrix"), 1, GL_FALSE, glm::value_ptr(viewproj * model));
        bool hasTexture = std::string(meshComponent->TexturePath) != "";
        glUniform1i(glGetUniformLocation(shader->ID, "useTexture"), hasTexture);

        if (hasTexture) {
            glBindTexture(GL_TEXTURE_2D, textures[textureCount]->texture);
            textureCount++;
        }
        else {
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        entity->vao->Bind();
        entity->vbo->Bind();
        entity->ebo->Bind();

        glDrawElements(GL_TRIANGLES, meshComponent->indices.size(), GL_UNSIGNED_INT, 0);

        entity->vao->Unbind();
        entity->vbo->Unbind();
        entity->ebo->Unbind();
    }

    // Render splines and reset state
    RenderSpline(viewproj, dt);
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);
    glUseProgram(shader->ID); // Ensure the correct shader program is still in use
}


void EntityManager::RenderSpline(glm::mat4 viewproj, float dt)
{
    for (auto spline : splines)
    {
        auto splineComponent = spline->GetComponent<SplineComponent>();
        if (splineComponent == nullptr)
            continue;
		splineComponent->CalculateBSpline();
        auto transform = spline->GetComponent<TransformComponent>();
        if (transform == nullptr) continue;
        // Model matrix calculations
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, glm::vec3(0.f));
        model = glm::scale(model, transform->scale);
        glUniformMatrix4fv(glGetUniformLocation(shader->ID, "camMatrix"), 1, GL_FALSE, glm::value_ptr(viewproj * model));
            glBindTexture(GL_TEXTURE_2D, 0);

        // Bind VAO, VBO, and EBO for drawing
        splineComponent->spline->vao->Bind();
        splineComponent->spline->vbo->Bind();
        splineComponent->spline->ebo->Bind();
        glDrawElements(GL_LINES, splineComponent->indices.size(), GL_UNSIGNED_INT, 0);
        // Unbind VAO, VBO, and EBO after drawing
        splineComponent->spline->vao->Unbind();
        splineComponent->spline->vbo->Unbind();
        //std::cout << "Drawing spline with indices size: " << splineComponent->indices.size() << std::endl;

        splineComponent->spline->ebo->Unbind();
    }

    glUseProgram(shader->ID); // Rebind the shader program for further rendering.
    glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray(0);

}




void EntityManager::ClearData()
{
}

bool EntityManager::HasNoEntities()
{
    return false;
}

void EntityManager::AddEntity(std::shared_ptr<Entity>& entity)
{
    if (entity->GetComponent<MeshComponent>() != nullptr)
    {
        initalizeMesh(entity);
    }
    entity->SetEntityID(EntityCount);
    ++EntityCount;
	entities.push_back(entity); //use smart pointers spawnsystem breakes here!
    if (entity->GetEntityID() == 0)
    {
		rigidbody->AddIndicesToCell(grid, entity);
    }
    else
	{
		grid->AddBaLL(entity);
	}
	if (entity->HasComponent<SplineComponent>())
	{
		splines.push_back(entity);
	}
}

std::vector<std::shared_ptr<Entity>> EntityManager::GetEntities() const
{
    return entities;
}

void EntityManager::RemoveLastEntity()
{
	if (EntityCount > 0)
	{
		if (entities.back()->GetComponent<MeshComponent>() != nullptr && entities.back()->GetComponent<MeshComponent>()->TexturePath != "")
		{
            textures.pop_back();
		}
		auto entity = entities.back();
        entities.pop_back();
		entity->Destroy();
        --EntityCount;
	}
	else
	{
		std::cerr << "No entities to remove!" << std::endl;
	}
}

void EntityManager::lauchDice()
{
	rigidbody->applyRandomForce(entities);
}

void EntityManager::initalizeMesh(std::shared_ptr<Entity>& entity)
{
    MeshComponent* meshComponent = entity->GetComponent<MeshComponent>();
    if (!meshComponent) {
        std::cerr << "Entity does not have a MeshComponent!" << std::endl;
        return;
    }
	initalizeTexture(entity);

    // Binding the VAO
    entity->vao->Bind();
    entity->vbo->Bind();

    glBufferData(GL_ARRAY_BUFFER, meshComponent->vertices.size() * sizeof(Vertex), meshComponent->vertices.data(), GL_STATIC_DRAW);

    entity->vao->LinkAttrib(*entity->vbo, 0, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, x)); // Position
    entity->vao->LinkAttrib(*entity->vbo, 1, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, r)); // Color
    entity->vao->LinkAttrib(*entity->vbo, 2, 2, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, u)); // TexCoords
    entity->vao->LinkAttrib(*entity->vbo, 3, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, normalx)); // TexCoords

    entity->ebo->Bind();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, meshComponent->indices.size() * sizeof(unsigned int), meshComponent->indices.data(), GL_STATIC_DRAW);

    // Unbinding VAO, VBO, EBO
    entity->vao->Unbind();
    entity->vbo->Unbind();
    entity->ebo->Unbind();

    if (entity->HasComponent<SplineComponent>()) {
        auto splineComp = entity->GetComponent<SplineComponent>();

        // Bind Spline VAO and VBO
        splineComp->spline->vao->Bind();
        splineComp->spline->vbo->Bind();

        // Upload vertex data
        glBufferData(GL_ARRAY_BUFFER, splineComp->vertices.size() * sizeof(Vertex), splineComp->vertices.data(), GL_STATIC_DRAW);

        // Link attributes
        splineComp->spline->vao->LinkAttrib(*splineComp->spline->vbo, 0, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, x)); // Position
        splineComp->spline->vao->LinkAttrib(*splineComp->spline->vbo, 1, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, r)); // Color
        splineComp->spline->vao->LinkAttrib(*splineComp->spline->vbo, 2, 2, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, u)); // TexCoords
        splineComp->spline->vao->LinkAttrib(*splineComp->spline->vbo, 3, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, normalx)); // Normals

        // Bind and upload indices
        splineComp->spline->ebo->Bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, splineComp->indices.size() * sizeof(unsigned int), splineComp->indices.data(), GL_STATIC_DRAW);

        // Unbind VAO/VBO/EBO for safety
        splineComp->spline->vao->Unbind();
        splineComp->spline->vbo->Unbind();
        splineComp->spline->ebo->Unbind();
    }

}

void EntityManager::initalizeTexture(std::shared_ptr<Entity>& entity)
{
	MeshComponent* meshComponent = entity->GetComponent<MeshComponent>();
	if (!meshComponent) {
		std::cerr << "Entity does not have a MeshComponent!" << std::endl;
		return;
	}

	//std::shared_ptr<Texture> texture(meshComponent->TexturePath, shader);
    std::shared_ptr<Texture> texture = std::make_shared<Texture>(meshComponent->TexturePath, shader);
	textures.push_back(texture);
}

