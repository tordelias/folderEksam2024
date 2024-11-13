#pragma once
#include "Component.h"
#include <memory>
#include <iostream>
#include <vector>
#include <string>

MeshComponent::MeshComponent(const char* figure, const glm::vec3& color, const char* texture)
    : TexturePath(texture)
{
    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();

    // Convert figure to std::string for reliable comparison
    std::string figureStr(figure);

    if (figureStr == "Cube") {
        auto [cubeVertices, cubeIndices] = mesh->CubeMesh(color);
        vertices = cubeVertices;
        indices = cubeIndices;
    }
    else if (figureStr == "Sphere") {
        auto [sphereVertices, sphereIndices] = mesh->SphereMesh(color);
        vertices = sphereVertices;
        indices = sphereIndices;
    }
    else if (figureStr == "Cylinder") {
        auto [cylinderVertices, cylinderIndices] = mesh->CylinderMesh(color);
        vertices = cylinderVertices;
        indices = cylinderIndices;
    }
    else if (figureStr == "Cone") {
        auto [coneVertices, coneIndices] = mesh->ConeMesh(color);
        vertices = coneVertices;
        indices = coneIndices;
    }
    else if (figureStr == "Torus") {
        auto [torusVertices, torusIndices] = mesh->TorusMesh(color);
        vertices = torusVertices;
        indices = torusIndices;
    }
    else if (figureStr == "PointCloud") {
        auto [pointCloudVertices, pointCloudIndices] = mesh->PointCloud(color);
        vertices = pointCloudVertices;
        indices = pointCloudIndices;
    }
    else {
        std::cerr << "Invalid figure name" << std::endl;
    }
}
