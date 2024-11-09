#ifndef MODEL_H
#define MODEL_H

#include <string>
#include <vector>
#include <memory>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <glm/glm.hpp>
#include <glad/glad.h>
#include "Resources/Shaders/shaderClass.h"

struct Vertex {
    glm::vec3 Position;
    glm::vec3 Normal;
    glm::vec2 TexCoords;
};

struct TextureData {
    unsigned int id;
    std::string type;
    std::string path;
};

class Model {
public:
    Model(const std::string& path);

    void loadModel(const std::string& path);

    void draw(const std::shared_ptr<Shader>& shader);

    std::vector<unsigned int> loadMaterialTextures(aiMaterial* mat, aiTextureType type, std::string typeName);

private:
    void processNode(aiNode* node, const aiScene* scene);

    void processMesh(aiMesh* mesh, const aiScene* scene);

    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    std::string directory;
    std::vector<unsigned int> textures;
};

#endif // MODEL_H
