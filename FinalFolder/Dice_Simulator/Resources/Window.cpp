#include "Window.h"
#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "Camera/Camera.h"

void Window::CreateWindow(int width, int height, const char* title, std::shared_ptr<Camera>& mainCam)
{
    // Directly set the camera, without condition
    camera = mainCam;

    // Create the GLFW window
    window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return;
    }

    glfwMakeContextCurrent(window);

    // Set the user pointer to this Window instance
    glfwSetWindowUserPointer(window, this);

    // Register the framebuffer size callback
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return;
    }

    // Initialize the viewport to the window size
    glViewport(0, 0, width, height);
}

void Window::ResizeWindow(int width, int height)
{
    if (!window)
    {
        std::cerr << "Window not created" << std::endl;
        return;
    }
    std::cout << "Resizing window to " << width << "x" << height << std::endl;

    // Update the camera dimensions
    camera->UpdateWindow(width, height);

    // Resize the GLFW window and update the OpenGL viewport
    glfwSetWindowSize(window, width, height);
    glViewport(0, 0, width, height);
}

glm::vec2 Window::GetWindowSize()
{
    if (window)
    {
        int width, height;
        glfwGetWindowSize(window, &width, &height);
        std::cout << "Window size is " << width << "x" << height << std::endl;
        return glm::vec2(width, height);
    }
    std::cerr << "Window not created" << std::endl;
    return glm::vec2();
}

void Window::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // Set the OpenGL viewport
    glViewport(0, 0, width, height);

    // Retrieve the Window instance from the user pointer
    Window* windowInstance = static_cast<Window*>(glfwGetWindowUserPointer(window));

    // Ensure that windowInstance and camera are valid
    if (windowInstance && windowInstance->camera) {
        // Directly update the camera with the new framebuffer size
        windowInstance->camera->UpdateWindow(width, height);
    }
}

double Window::GetMouseXPos()
{
    if (window)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        return xpos;
    }
    std::cerr << "Window not created" << std::endl;
    return 0.0;
}

// Get the mouse Y position
double Window::GetMouseYPos()
{
    if (window)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        return ypos;
    }
    std::cerr << "Window not created" << std::endl;
    return 0.0;
}

// Get the combined mouse position as glm::vec2
glm::vec2 Window::GetMousePos()
{
    if (window)
    {
        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);
        return glm::vec2(xpos, ypos);
    }
    std::cerr << "Window not created" << std::endl;
    return glm::vec2(0.0, 0.0);
}