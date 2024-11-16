#version 330 core

layout (location = 0) in vec3 aPos;       // Position
layout (location = 1) in vec3 aColor;     // Color
layout (location = 2) in vec2 aTexCoords; // Texture coordinates
layout (location = 3) in vec3 aNormal;    // Normal vector

out vec3 FragPos;    // Position in world space
out vec3 Normal;     // Normal vector
out vec3 color;      // Vertex color
out vec2 TexCoords;  // Texture coordinates

uniform mat4 model;      // Model matrix
uniform mat4 camMatrix;  // Combined View-Projection matrix

void main() {
    // Transform the vertex position into world space
    FragPos = vec3(model * vec4(aPos, 1.0));

    // Transform the normal to world space
    Normal = mat3(transpose(inverse(model))) * aNormal;

    // Pass through other attributes
    color = aColor;
    TexCoords = aTexCoords;

    // Compute final vertex position
    gl_Position = camMatrix * vec4(FragPos, 1.0);
}
