#pragma once
#include <glm/glm.hpp>
#include "vector"
#include "memory"
#include "../Component/Mesh.h"
class VAO;
class VBO;
class EBO;
class Shader;
class ParticleSystem
{
public:

	ParticleSystem(glm::vec3 pos, glm::vec3 acc,glm::vec3 scale, int maxparticles);
	~ParticleSystem();
	void emit();
	void update(glm::mat4 viewproj, std::shared_ptr<Shader> shader, glm::vec3 Opos, float dt);
	void draw(glm::mat4 viewproj, std::shared_ptr<Shader> shader);
	void reset(glm::vec3& Ppos, glm::vec3& Pvel, float& life);
	void updateacceleration(glm::vec3 acc);
private: 
	std::vector<glm::vec3> position;
	std::vector<glm::vec3> velocity;
	std::vector<glm::vec3> acceleration;
	std::vector<float> lifetime;
	std::vector<std::shared_ptr<VAO>> Vao;
	std::vector<std::shared_ptr<VBO>> Vbo;
	std::vector<std::shared_ptr<EBO>> Ebo;

	glm::vec3 pos; 
	glm::vec3 acc; 
	int maxParticles = 500;
	std::vector<unsigned int> indices;
	std::vector<Vertex> vertices;
	glm::vec3 color = glm::vec3(1, 1, 1);
	glm::vec3 scale = glm::vec3(1.f);
	glm::vec3 size = glm::vec3(0.03f);
	int maxLife = 8;
	int minLife = 2;
};

