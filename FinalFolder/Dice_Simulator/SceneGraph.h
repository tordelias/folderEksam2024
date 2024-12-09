#pragma once
#include <iostream>
#include <vector>
#include <memory>
#include "Entity.h"
#include "Component/Component.h"

//class SceneGraph {
//public:
//    void AddEntity(std::shared_ptr<Entity> entity) {
//        entities.push_back(entity);
//    }
//    void RemoveEntity(std::shared_ptr<Entity> entity) {
//        entities.erase(std::remove(entities.begin(), entities.end(), entity), entities.end());
//
//        // Remove entity from its parent (if it has one)
//        if (auto parent = entity->GetParent()) {
//            parent->RemoveChild(entity);
//        }
//    }
//
//    // Updates all entities in the scene graph
//    void Update() {
//        for (auto& entity : entities) {
//            auto parentTransform = entity->GetComponent<TransformComponent>();
//            entity->Update(parentTransform);
//        }
//    }
//
//    // Renders all entities in the scene graph
//    void Render() const {
//        for (const auto& entity : entities) {
//            entity->Render();
//        }
//    }
//
//private:
//    std::vector<std::shared_ptr<Entity>> entities;
//};
