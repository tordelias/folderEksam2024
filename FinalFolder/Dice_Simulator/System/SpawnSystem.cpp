#include "SpawnSystem.h"
#include <iostream>
#include <string>
#include "../Entity.h"
#include "../Component/Component.h"
#include <glm/glm.hpp>
#include <GLFW/glfw3.h>
#include <memory>
#include <random>


SpawnSystem::SpawnSystem(std::shared_ptr<EntityManager> entityManager) : manager(entityManager)
{
}

SpawnSystem::~SpawnSystem()
{
}

void SpawnSystem::input(GLFWwindow* window)
{
	bool isQPressed = glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS;
	if (isQPressed && !spacePressedLastFrame)
	{
		 SpawnEntity(offset, 0, -5);
	}
	spacePressedLastFrame = isQPressed;

	bool isFPressed = glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS;
	if (isFPressed && !FPressedLastFrame)
	{
		deletelastEntity();

	}
	FPressedLastFrame = isFPressed;

	if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS && !blaunchDiece)
	{
		manager->lauchDice();
		blaunchDiece = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_L) == GLFW_RELEASE)
	{
		blaunchDiece = false;
	}
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS && !bPPressed)
	{
		if (bFillTriangles)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			bFillTriangles = false;
		}
		else
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			bFillTriangles = true;
		}
		bPPressed = true;
	}
	else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_RELEASE)
	{
		bPPressed = false;
	}
}

void SpawnSystem::SpawnEntity()
{
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f));
	cube->AddComponent<MeshComponent>("Cube", glm::vec3(1.0f, 1.0f, 1.0f), "");
	manager->AddEntity(cube);
	offset += offsetAmount;

}

#include <random> // Include for random number generation

void SpawnSystem::SpawnEntity(int x, int y, int z)
{
	// Set up random number generation for color
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> dis(0.0, 1.0);

	// Generate random color values
	float r = dis(gen);
	float g = dis(gen);
	float b = dis(gen);

	// Create the entity and components
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(x, y, z), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f));

	// Set random color for the MeshComponent
	cube->AddComponent<MeshComponent>("Cube", glm::vec3(r, g, b), "");

	// Add entity to the manager
	manager->AddEntity(cube);
	offset += offsetAmount;
}


void SpawnSystem::SpawnEntity(int x, int y, int z, const char* texturePath)
{
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(x, y, z), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f));
	cube->AddComponent<MeshComponent>("Cube", glm::vec3(1.0f, 1.0f, 1.0f), texturePath);
	manager->AddEntity(cube);
	offset += offsetAmount;
}

void SpawnSystem::SpawnEntity(int x, int y, int z, const char* texturePath, const char* meshType)
{
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(x, y, z), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(1.0f));
	cube->AddComponent<MeshComponent>(meshType, glm::vec3(1.0f, 1.0f, 1.0f), texturePath);
	manager->AddEntity(cube);
	offset += offsetAmount;
}

void SpawnSystem::SpawnEntity(int x, int y, int z, const char* texturePath, const char* meshType, float scale)
{
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(x, y, z), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(scale));
	cube->AddComponent<MeshComponent>(meshType, glm::vec3(1, 1, 1), texturePath);
	manager->AddEntity(cube);
	offset += offsetAmount;
}

void SpawnSystem::SpawnEntity(int x, int y, int z, const char* texturePath, const char* meshType, float scale, glm::vec3 rotation)
{
	std::shared_ptr<Entity> cube = std::make_shared<Entity>();
	cube->AddComponent<TransformComponent>(glm::vec3(x, y, z), rotation, glm::vec3(scale));
	cube->AddComponent<MeshComponent>(meshType, glm::vec3(1,1,1), texturePath);
	manager->AddEntity(cube);
	offset += offsetAmount;
}

void SpawnSystem::deletelastEntity()
{
	manager->RemoveLastEntity();
	offset -= offsetAmount;
}
