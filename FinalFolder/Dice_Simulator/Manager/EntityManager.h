#pragma once
#include "../Entity.h"
#include <glm/glm.hpp>
class Shader;
class Texture; 
class TransformComponent;
class RigidBody;
class Grid; 
class EntityManager
{
public:
	EntityManager(std::shared_ptr<Shader> shader);
	~EntityManager();
	void Update();
	void Render(glm::mat4 viewproj, float dt);
	void RenderSpline(glm::mat4 viewproj, float dt);
	void ClearData();
	bool HasNoEntities();
	void AddEntity(std::shared_ptr<Entity>& entity);
	std::vector<std::shared_ptr<Entity>> GetEntities() const;
	unsigned int GetEntityCount() { return EntityCount; };
	void RemoveLastEntity();

	void lauchDice();

private:
	void initalizeMesh(std::shared_ptr<Entity>& entity);
	void initalizeTexture(std::shared_ptr<Entity>& entity);
	std::vector<std::shared_ptr<Entity>> entities;
	std::vector<std::shared_ptr<Entity>> splines;
	int EntityCount;
	std::vector<std::shared_ptr<Texture>> textures;
	std::shared_ptr<Shader> shader;	
	std::shared_ptr<RigidBody> rigidbody;
	std::shared_ptr<Grid> grid;

};

