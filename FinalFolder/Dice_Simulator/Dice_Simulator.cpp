#pragma once
#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "Resources/Camera/Camera.h"
#include "Resources/Shaders/shaderClass.h"
#include <memory>
#include <string>
#include "Resources/Window.h"

//includes
#include "Manager/EntityManager.h"
#include "Entity.h" 
#include "Component/Component.h"
#include "System/SpawnSystem.h"

void processInput(Window window);

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 800;
bool bPPressed = false;

int main()
{
	// ---------------------------------------------------------------------------------------------------------------------------
	//                                                        Window & Camera & Shader
	// ---------------------------------------------------------------------------------------------------------------------------
	Window window;

    std::shared_ptr<Camera> camera = std::make_shared<Camera>(SCR_WIDTH, SCR_HEIGHT, glm::vec3(0.0f, 2.0f, 0.0f));
	window.CreateWindow(SCR_WIDTH, SCR_HEIGHT, "Dice Simulator", camera);

    std::shared_ptr<Shader> shaderProgram = std::make_shared<Shader>("Resources/Shaders/default.vert", "Resources/Shaders/default.frag");
	glm::vec3 lightPos(0,0,0);
	glm::mat4 lightModel = glm::mat4(1.0f);
	lightModel = glm::translate(lightModel, lightPos);

    shaderProgram->Activate();

	glUniformMatrix4fv(glGetUniformLocation(shaderProgram->ID, "model"), 1, GL_FALSE, glm::value_ptr(lightModel));
	glUniform3fv(glGetUniformLocation(shaderProgram->ID, "lightColor"), 1, glm::value_ptr(glm::vec3(1.0f, 1.0f, 1.0f)));
	glUniform3fv(glGetUniformLocation(shaderProgram->ID, "lightPos"), 1, glm::value_ptr(lightPos));

    // ---------------------------------------------------------------------------------------------------------------------------
	//                                                        Initialize bellow
    // ---------------------------------------------------------------------------------------------------------------------------
	std::shared_ptr<EntityManager> manager = std::make_shared<EntityManager>(shaderProgram);
	std::shared_ptr<SpawnSystem> spawnSystem = std::make_shared<SpawnSystem>(manager);

	std::shared_ptr<Entity> entity = std::make_shared<Entity>();
	entity->AddComponent<TransformComponent>(glm::vec3(0,0,-5), glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(10.f));
	entity->AddComponent<MeshComponent>("Torus", glm::vec3(0.f, 1.f, 0.f), "");
	manager->AddEntity(entity);
	glPointSize(1.0f);


    // ---------------------------------------------------------------------- -----------------------------------------------------
    //                                                        other
    // ---------------------------------------------------------------------------------------------------------------------------

    glEnable(GL_DEPTH_TEST);
    double lastTime = glfwGetTime();  // Initialize lastTime to current time
    float deltaTime;
    double currentTime;

    // ---------------------------------------------------------------------------------------------------------------------------
    //                                                        Main Loop
    // ---------------------------------------------------------------------------------------------------------------------------
    while (!glfwWindowShouldClose(window.GetWindow()))
    {
         currentTime = glfwGetTime();   // Get the current time
        deltaTime = static_cast<float>(currentTime - lastTime); // Calculate deltaTime
        lastTime = currentTime;

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        processInput(window);
        camera->Inputs(window.GetWindow());
        glm::mat4 viewproj = camera->Matrix(45.0f, 0.1f, 30000.0f, *shaderProgram, "camMatrix");        //Set render distance and FOV
    // ---------------------------------------------------------------------------------------------------------------------------
       
		spawnSystem->input(window.GetWindow(), camera);
		manager->Render(viewproj, deltaTime);



    // ---------------------------------------------------------------------------------------------------------------------------
        glfwSwapBuffers(window.GetWindow());
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void processInput(Window window)
{
    if (glfwGetKey(window.GetWindow(), GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window.GetWindow(), true);
}
