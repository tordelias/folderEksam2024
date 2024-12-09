#include "ParticleSystem.h"
#define GLM_ENABLE_EXPERIMENTAL
#include <random>
#include "../Resources/Shaders/shaderClass.h"
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include "../Resources/Shaders/EBO.h"
#include "../Resources/Shaders/VBO.h"
#include "../Resources/Shaders/VAO.h"



ParticleSystem::ParticleSystem(glm::vec3 pos, glm::vec3 acc, glm::vec3 radius, glm::vec3 particleSize,glm::ivec2 life, glm::vec3 Color, int maxparticles, int figure = 0) : pos(pos), acc(acc), maxParticles(maxparticles), scale(radius), size(particleSize), color(Color)
{
	Mesh mesh;
	if (life.x < life.y) {
		minLife = life.x;
		maxLife = life.y;
	}
	else
	{
		std::cerr << "Error: minLife is greater than maxLife" << std::endl;
	}
	if(figure == 0)
	{
		auto [cubeVertices, cubeIndices] = mesh.SphereMesh(color);
		vertices = cubeVertices;
		indices = cubeIndices;
	}
	else
	{
		auto [cubeVertices, cubeIndices] = mesh.CubeMesh(color);
		vertices = cubeVertices;
		indices = cubeIndices;
	}
	for (int i = 0; i < maxParticles; ++i)
	{
		emit();
	}
}

ParticleSystem::~ParticleSystem()
{
}

void ParticleSystem::emit()
{
	 std::random_device rd;
	 std::mt19937 gen(rd());
	 std::uniform_real_distribution<> disX(-scale.x, scale.x);
	 std::uniform_real_distribution<> disY(-scale.y, scale.y);
	 std::uniform_real_distribution<> disZ(-scale.z, scale.z);
	 std::uniform_real_distribution<> dis(minLife, maxLife);

	std::shared_ptr<VAO> vao = std::make_shared<VAO>();
	std::shared_ptr<VBO> vbo = std::make_shared<VBO>();
	std::shared_ptr<EBO> ebo = std::make_shared<EBO>();

	position.push_back(pos + glm::vec3(0,scale.y * 12.f,0) + glm::vec3(disX(gen), disY(gen), disZ(gen)));
	velocity.push_back(glm::vec3(0, 0, 0));
	acceleration.push_back(acc * glm::vec3(0, dis(gen), 0));
	lifetime.push_back(dis(gen));

	vao->Bind();
	vbo->Bind();

	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

	vao->LinkAttrib(*vbo, 0, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, x)); // Position
	vao->LinkAttrib(*vbo, 1, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, r)); // Color
	ebo->Bind();
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	// Unbinding VAO, VBO, EBO
	vao->Unbind();
	vbo->Unbind();
	ebo->Unbind();

	Vao.push_back(vao);
	Vbo.push_back(vbo);
	Ebo.push_back(ebo);  // Add EBO to the vector as well
}

void ParticleSystem::update(glm::mat4 viewproj, std::shared_ptr<Shader> shader,glm::vec3 Opos, float dt)
{
	for (int i = 0; i < maxParticles; i++)
	{
		velocity[i] += acceleration[i] * dt;
		position[i] += velocity[i] * dt;
		lifetime[i] -= dt;
		pos = Opos;
		if (lifetime[i] <= 0)
		{
			reset(position[i], velocity[i], lifetime[i]);
		}
	}
	draw(viewproj, shader);
}

void ParticleSystem::draw(glm::mat4 viewproj, std::shared_ptr<Shader> shader)
{
	shader->Activate();
	for (int i = 0; i < maxParticles; i++)
	{
		glm::mat4 model = glm::mat4(1.0f);
		model = glm::translate(model, position[i]);
		model = glm::scale(model, size);
		glUniformMatrix4fv(glGetUniformLocation(shader->ID, "camMatrix"), 1, GL_FALSE, glm::value_ptr(viewproj * model));
		Vao[i]->Bind();
		Vbo[i]->Bind();
		Ebo[i]->Bind();

		glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);

		Vao[i]->Unbind();
		Vbo[i]->Unbind();
		Ebo[i]->Unbind();
	}
}

void ParticleSystem::reset(glm::vec3& Ppos, glm::vec3& Pvel, float& life)
{
	 std::random_device rd;
	 std::mt19937 gen(rd());
	 std::uniform_real_distribution<> disX(-scale.x, scale.x);
	 std::uniform_real_distribution<> disY(-scale.y, scale.y);
	 std::uniform_real_distribution<> disZ(-scale.z, scale.z);
	 std::uniform_real_distribution<> dis(minLife, maxLife);

	Ppos = pos + glm::vec3(0,scale.y * 10.f,0) + glm::vec3(disX(gen), disY(gen), disZ(gen));
	Pvel = glm::vec3(0, 0, 0);
	 life = dis(gen);
}

void ParticleSystem::updateacceleration(glm::vec3 acc)
{
	for (int i = 0; i < maxParticles; i++)
	{
		acceleration[i] = acc;
	}
}

