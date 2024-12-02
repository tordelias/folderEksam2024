#pragma once
#include "Component.h"
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include "../Entity.h"
#include "../Resources/Shaders/VAO.h"
#include "../Resources/Shaders/VBO.h"
#include "../Resources/Shaders/EBO.h"


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
	else if (figureStr == "BSplineSurface") {
		auto [bsplineVertices, bsplineIndices] = mesh->BSplineSurface(color);
		vertices = bsplineVertices;
		indices = bsplineIndices;
	}
	else if (figureStr == "testing")
	{
		auto [twoHillVertices, twoHillIndices] = mesh->TwoHillFlatMiddle(color);
		vertices = twoHillVertices;
		indices = twoHillIndices;
	}
    else {
        std::cerr << "Invalid figure name" << std::endl;
    }
}

SplineComponent::SplineComponent(std::shared_ptr<Entity> ownerptr)
{
    owner = ownerptr;
	startpos = owner->GetComponent<TransformComponent>()->position;
	spline = std::make_shared<Entity>();
}

void SplineComponent::CalculateBSpline()
{
    std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
    if (glm::length(owner->GetComponent<TransformComponent>()->velocity) > 0.f)
    {
        auto [bsplineVertices, bsplineIndices] = mesh->Spline(glm::vec3(1.0f, 1.0f, 1.0f), owner->GetComponent<TransformComponent>()->position, controllpoints);
        vertices = bsplineVertices;
        indices = bsplineIndices;

        spline->vao->Bind();
        spline->vbo->Bind();
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(), GL_STATIC_DRAW);

        spline->vao->LinkAttrib(*spline->vbo, 0, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, x)); // Position
        spline->vao->LinkAttrib(*spline->vbo, 1, 3, GL_FLOAT, sizeof(Vertex), (void*)offsetof(Vertex, r)); // Color

        spline->ebo->Bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

        spline->vao->Unbind();
    }

}



