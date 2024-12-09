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
    // Debugging normals by outputting them directly
    // Uncomment the next line to visualize normals as RGB colors:
    // FragColor = vec4(normalize(Normal), 1.0);

    // Ambient light (increase a bit to brighten the scene)
    vec3 ambient = 0.9 * lightColor; // Raised ambient intensity to 0.5

    // Diffuse light
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Blinn-Phong specular light (softened specular)
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 halfDir = normalize(viewDir + lightDir); // Halfway vector for Blinn-Phong
    float spec = pow(max(dot(norm, halfDir), 0.0), 16.0); // Reduced shininess exponent
    vec3 specular = 0.2 * spec * lightColor; // Lower specular intensity

    // Combine lighting effects (more ambient to brighten the world)
    vec3 lighting = ambient + diffuse + specular;

    // Clamp the lighting values to avoid over-bright spots
    lighting = clamp(lighting, 0.0, 1.0);

    // Apply texture or vertex color
    vec3 baseColor = useTexture ? texture(ourTexture, TexCoords).rgb : color;

    // Final color (apply lighting and base color)
    vec3 finalColor = (lighting * baseColor);

    // Output the final color without gamma correction
    FragColor = vec4(finalColor, 1.0);
}
