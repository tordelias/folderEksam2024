#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include  <string>
#include <memory>
#include <glm/glm.hpp>

class Camera;
class Window
{
public:
	Window() : window(nullptr)
	{
		glfwInit();
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
		glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	};

	void CreateWindow(int width, int height, const char* title, std::shared_ptr<Camera>& mainCam);
	void ResizeWindow(int width, int height);
	glm::vec2 GetWindowSize();

	static void framebuffer_size_callback(GLFWwindow* window, int width, int height);  // Dynamic window size

	GLFWwindow* GetWindow() { return window; }
private:
	GLFWwindow* window;
	std::shared_ptr<Camera> camera;
};

