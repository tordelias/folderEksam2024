#version 330 core
// Outputs colors in RGBA
out vec4 FragColor;

// Inputs from the Vertex Shader
in vec3 color;
in vec2 TexCoords;

// Uniforms
uniform sampler2D ourTexture;
uniform bool useTexture; // New uniform to toggle texture usage

struct Material {
    vec3 diffuse;
};

void main() {
    Material mat;

    // Check if a texture should be used
    if (useTexture) {
        mat.diffuse = vec3(texture(ourTexture, TexCoords)); // Use texture color
    } else {
        mat.diffuse = color; // Use vertex color directly
    }

    vec3 finalColor = mat.diffuse; // Final color for output

    FragColor = vec4(finalColor, 1.0f); // Output final color as RGBA
}
