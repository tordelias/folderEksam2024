#version 330 core
//Code is from:Compulsory-2-3DProg/OpenGLSession0/Resources/Shaders/default.vert
// Github Link: https://github.com/HansPluss/Compulsory-2-3DProg.git 
// Positions/Coordinates
layout (location = 0) in vec3 aPos;
// Colors
layout (location = 1) in vec3 aColor;
// Texture Coords
layout (location = 2) in vec2 aTexCoords; 

// Outputs the color for the Fragment Shader
out vec3 color;
out vec2 TexCoords; 

// Controls the scale of the vertices
uniform float scale;
uniform mat4 camMatrix;

void main() {
    // Outputs the positions/coordinates of all vertices
    gl_Position = camMatrix * vec4(aPos, 1.0);
    // Assigns the colors from the Vertex Data to "color"
    color = aColor;
    TexCoords = aTexCoords;
}
