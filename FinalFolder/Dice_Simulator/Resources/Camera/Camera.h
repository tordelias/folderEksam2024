#define GLM_ENABLE_EXPERIMENTAL
#ifndef CAMERA_CLASS_H
#define CAMERA_CLASS_H

#include<glad/glad.h>
#include<GLFW/glfw3.h>
#include<glm/glm.hpp>
#include<glm/gtc/matrix_transform.hpp>
#include<glm/gtc/type_ptr.hpp>
#include<glm/gtx/rotate_vector.hpp>
#include<glm/gtx/vector_angle.hpp>
#include "../Shaders/shaderClass.h"

//Code is from:Compulsory-2-3DProg/OpenGLSession0/LSM.h
// Github Link: https://github.com/HansPluss/Compulsory-2-3DProg.git 

class Camera
{
public:
	// Stores the main vectors of the camera
	glm::vec3 Position;
	glm::vec3 Orientation = glm::vec3(0.f, -0.2f, -1.0f);
	glm::vec3 Up = glm::vec3(0.0f, 1.0f, 0.0f);

	// Prevents the camera from jumping around when first clicking left click
	bool firstClick = true;

	// Stores the width and height of the window
	int width;
	int height;

	// Adjust the speed of the camera and it's sensitivity when looking around
	float speed = 0.1f;
	float sensitivity = 100.0f;

	// Camera constructor to set up initial values
	Camera(int width, int height, glm::vec3 position);


	// Updates and exports the camera matrix to the Vertex Shader
	glm::mat4 Matrix(float FOVdeg, float nearPlane, float farPlane, Shader& shader, const char* uniform);
	// Handles camera inputs
	void Inputs(GLFWwindow* window);

	void UpdateWindow(unsigned int width, unsigned int height)
	{
		this->width = width;
		this->height = height;
	}

};
#endif
