#version 330 core

// Outputs to the screen
out vec4 FragColor;

// Inputs from Vertex Shader
in vec3 FragPos;    // Position in world space
in vec3 Normal;     // Normal vector
in vec3 color;      // Vertex color
in vec2 TexCoords;  // Texture coordinates

// Uniforms
uniform sampler2D ourTexture;
uniform bool useTexture;
uniform vec3 lightPos;    // Light position
uniform vec3 lightColor;  // Light color
uniform vec3 viewPos;     // Camera position

void main() {
    // Ambient light
    vec3 ambient = 0.9 * lightColor;

    // Diffuse light
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Specular light
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = spec * lightColor;

    // Combine lighting effects
    vec3 lighting = ambient + diffuse + specular;

    // Apply texture or vertex color
    vec3 baseColor = useTexture ? texture(ourTexture, TexCoords).rgb : color;
    vec3 finalColor = lighting * baseColor;

    FragColor = vec4(finalColor, 1.0);
}
