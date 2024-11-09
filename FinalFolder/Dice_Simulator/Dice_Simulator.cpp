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

void processInput(Window window, std::shared_ptr<Camera> camera);

const unsigned int SCR_WIDTH = 1920;
const unsigned int SCR_HEIGHT = 1080;
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
    shaderProgram->Activate();

    // ---------------------------------------------------------------------------------------------------------------------------
	//                                                        Initialize bellow
    // ---------------------------------------------------------------------------------------------------------------------------
	std::shared_ptr<EntityManager> manager = std::make_shared<EntityManager>(shaderProgram);
	std::shared_ptr<SpawnSystem> spawnSystem = std::make_shared<SpawnSystem>(manager);

  //  spawnSystem->SpawnEntity(0, 0, 0, "Resources/Texture/Textures/skybox.jpg", "Cube", 1000.0f);
 //   spawnSystem->SpawnEntity(0, 0, -10, "Resources/Texture/Textures/beako.png");
	//spawnSystem->SpawnEntity(2, 0, -10, "Resources/Texture/Textures/beako.png");
	//spawnSystem->SpawnEntity(-2, 0, -10, "");


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

        processInput(window, camera);
        camera->Inputs(window.GetWindow());
        glm::mat4 viewproj = camera->Matrix(45.0f, 0.1f, 1000.0f, *shaderProgram, "camMatrix");        //Set render distance and FOV
    // ---------------------------------------------------------------------------------------------------------------------------
       
		spawnSystem->input(window.GetWindow());
		manager->Render(viewproj, deltaTime);




    // ---------------------------------------------------------------------------------------------------------------------------
        glfwSwapBuffers(window.GetWindow());
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

void processInput(Window window, std::shared_ptr<Camera> camera)
{
    if (glfwGetKey(window.GetWindow(), GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window.GetWindow(), true);
	if (glfwGetKey(window.GetWindow(), GLFW_KEY_P) == GLFW_PRESS && !bPPressed)
	{
		bPPressed = true;
		window.ResizeWindow(400, 800);
	} 
    else if (glfwGetKey(window.GetWindow(), GLFW_KEY_P) == GLFW_RELEASE)
	{
		bPPressed = false;
    }
}
