#pragma once
#include "Component.h"
#include <memory>
#include <iostream>
#include <vector>

MeshComponent::MeshComponent(const char* figure, const glm::vec3& color, const char* texture)
    : TexturePath(texture)
{
    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();

    if (figure == "Cube") {
        auto [cubeVertices, cubeIndices] = mesh->CubeMesh(color);
        vertices = cubeVertices;
        indices = cubeIndices;
    }
    else if (figure == "Sphere") {
        auto [sphereVertices, sphereIndices] = mesh->SphereMesh(color);
        vertices = sphereVertices;
        indices = sphereIndices;
    }
    else if (figure == "Cylinder") {
        auto [cylinderVertices, cylinderIndices] = mesh->CylinderMesh(color);
        vertices = cylinderVertices;
        indices = cylinderIndices;
    }
    else if (figure == "Cone") {
        auto [coneVertices, coneIndices] = mesh->ConeMesh(color);
        vertices = coneVertices;
        indices = coneIndices;
    }
    else if (figure == "Torus") {
        auto [torusVertices, torusIndices] = mesh->TorusMesh(color);
        vertices = torusVertices;
        indices = torusIndices;
	}
	else if (figure == "PointCloud") 
    {
		auto [pointCloudVertices, pointCloudIndices] = mesh->PointCloud(color);
		vertices = pointCloudVertices;
		indices = pointCloudIndices;
	}
    else {
        std::cerr << "Invalid figure name" << std::endl;
    }
}
