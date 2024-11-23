#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <omp.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/IO/Polyhedron_iostream.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> BaseVertex;
typedef CGAL::Triangulation_vertex_base_2<K, CGAL::Triangulation_ds_vertex_base_2<BaseVertex>> VertexBase;
typedef CGAL::Triangulation_data_structure_2<VertexBase> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;

// Define a custom vertex with index
typedef Delaunay::Vertex_handle Vertex_handle_with_index;
typedef K::Point_2 Point_2;


// Cube mesh generation
std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::CubeMesh(glm::vec3 color) {
    std::vector<Vertex> vertices = {
        // Front face (z = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  0.0f, -1.0f, 0 }, // Bottom-left
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  0.0f, -1.0f, 1 }, // Bottom-right
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  0.0f, -1.0f, 2 }, // Top-right
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  0.0f, -1.0f, 3 }, // Top-left

        // Back face (z = 0.5f)
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  0.0f,  1.0f, 4 }, // Bottom-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  0.0f,  1.0f, 5 }, // Bottom-left
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  0.0f,  1.0f, 6 }, // Top-left
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  0.0f,  1.0f, 7 }, // Top-right

        // Left face (x = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f, -1.0f,  0.0f,  0.0f, 8 }, // Bottom-right
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f, -1.0f,  0.0f,  0.0f, 9 }, // Top-right
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f, -1.0f,  0.0f,  0.0f, 10 }, // Top-left
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f, -1.0f,  0.0f,  0.0f, 11 }, // Bottom-left

        // Right face (x = 0.5f)
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  1.0f,  0.0f,  0.0f, 12 }, // Bottom-left
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  1.0f,  0.0f,  0.0f, 13 }, // Top-left
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  1.0f,  0.0f,  0.0f, 14 }, // Top-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  1.0f,  0.0f,  0.0f, 15 }, // Bottom-right

        // Top face (y = 0.5f)
        { -0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f,  1.0f,  0.0f, 16 }, // Top-left
        {  0.5f,  0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f,  1.0f,  0.0f, 17 }, // Top-right
        {  0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f,  1.0f,  0.0f, 18 }, // Bottom-right
        { -0.5f,  0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f,  1.0f,  0.0f, 19 }, // Bottom-left

        // Bottom face (y = -0.5f)
        { -0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 0.0f, 0.0f,  0.0f, -1.0f,  0.0f, 20 }, // Bottom-left
        {  0.5f, -0.5f, -0.5f, color.r, color.g, color.b, 1.0f, 0.0f,  0.0f, -1.0f,  0.0f, 21 }, // Bottom-right
        {  0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 1.0f, 1.0f,  0.0f, -1.0f,  0.0f, 22 }, // Top-right
        { -0.5f, -0.5f,  0.5f, color.r, color.g, color.b, 0.0f, 1.0f,  0.0f, -1.0f,  0.0f, 23 }  // Top-left
    };

    std::vector<unsigned int> indices = {
        0, 1, 2, 2, 3, 0,   // Front
        4, 5, 6, 6, 7, 4,   // Back
        8, 9, 10, 10, 11, 8, // Left
        12, 13, 14, 14, 15, 12, // Right
        16, 17, 18, 18, 19, 16, // Top
        20, 21, 22, 22, 23, 20  // Bottom
    };

    return { vertices, indices };
}


std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::SphereMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    int radius = 1.0f;
    int circumferenceTile = 18;
    int layerTile = 18;
    int slices = (int)(circumferenceTile + 0.5f);
    if (slices < 4) {
        slices = 4;
    }

    int half_slices = slices / 2;
    int layerCount = (int)(layerTile + 0.5f);
    if (layerCount < 2)
    {
        layerCount = 2;
    }
    float pi = 3.1415f;
    for (int layerIndex = 0; layerIndex <= layerCount; layerIndex++)
    {
        float v = (1.0 - (float)layerIndex / layerCount);
        float heightOffset = std::sin((1.0 - 2.0 * layerIndex / layerCount) * pi / 2.0);
        float cosUp = sqrt(1.0 - heightOffset * heightOffset);
        float z = heightOffset;
        for (int i = 0; i <= half_slices; i++)
        {
            float u = (float)i / (float)half_slices;
            float angle = 2 * pi * u; // pi * 2 to get full sphere
            float x = std::cos(angle) * cosUp;
            float y = std::sin(angle) * cosUp;
            Vertex V1 = Vertex{ x * radius, y * radius, z * radius, x, y, z, u, v };
            vertices.push_back(V1);
        }

    }
    for (int layer = 0; layer < layerCount; layer++)
    {
        for (int i = 0; i < half_slices; i++)
        {
            // Index for the current layer and the next layer
            int currentRow = layer * (half_slices + 1) * 2;
            int nextRow = (layer + 1) * (half_slices + 1) * 2;

            // Creating two triangles (quad) between each pair of vertices in adjacent layers
            indices.push_back(currentRow + i);        // 1st triangle: curRow, nextRow, nextRow+1
            indices.push_back(nextRow + i);
            indices.push_back(nextRow + i + 1);

            indices.push_back(currentRow + i);        // 2nd triangle: curRow, nextRow+1, curRow+1
            indices.push_back(nextRow + i + 1);
            indices.push_back(currentRow + i + 1);
        }
    }


    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::CylinderMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;


    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::ConeMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;



    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::TorusMesh(glm::vec3 color)
{
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    // Define the size of the terrain grid
    int width = 10;   // Number of vertices along the X axis
    int depth = 10;   // Number of vertices along the Z axis
    float spacing = 1.0f;  // Space between vertices

    // Center the grid around (0, 0, 0)
    float halfWidth = width / 2.0f;
    float halfDepth = depth / 2.0f;
	int count = 0;

    // Generate vertices with height variation along the Y axis
    for (int z = 0; z < depth; ++z) {
        for (int x = 0; x < width; ++x) {
            Vertex vertex;

            // Position (X and Z are used for the grid, Y is the height)
            vertex.x = (x - halfWidth) * spacing;  // Shift the x position to center
            vertex.z = (z - halfDepth) * spacing;  // Shift the z position to center
            vertex.y = std::sin(x * 0.5f) * std::cos(z * 0.5f) * 2.0f;  // Height variation using sine and cosine

            // Color
            vertex.r = color.r;
            vertex.g = color.g;
            vertex.b = color.b;

            // Texture coordinates (simple mapping)
            vertex.u = (float)x / (width - 1);
            vertex.v = (float)z / (depth - 1);

            // Normals (simplified, pointing upwards)
            vertex.normalx = 0.0f;
            vertex.normaly = 1.0f;
            vertex.normalz = 0.0f;

            // Index (for use in rendering)
            vertex.index = vertices.size();

            // Friction (random value, for example)
			if (count % 6 == 0)
            {
                vertex.friction = 0.9f;
                vertex.r = 1.f;
				vertex.g = 0.f;
				vertex.b = 0.f;
            }
			else
                vertex.friction = 0.2f;
			count++;
            // Add vertex to the list
            vertices.push_back(vertex);
        }
    }

    // Generate indices for a grid of triangles
    for (int z = 0; z < depth - 1; ++z) {
        for (int x = 0; x < width - 1; ++x) {
            // Create two triangles for each grid square
            unsigned int topLeft = z * width + x;
            unsigned int topRight = z * width + (x + 1);
            unsigned int bottomLeft = (z + 1) * width + x;
            unsigned int bottomRight = (z + 1) * width + (x + 1);

            // First triangle (top-left, top-right, bottom-left)
            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            // Second triangle (top-right, bottom-left, bottom-right)
            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    // Return the vertices and indices as a pair
    return { vertices, indices };
}

std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::PointCloud(glm::vec3 color) {
    // Read the point cloud from the file (with progress bar in Readfile)
    std::vector<Vertex> vertices = Readfile("Data/32-2-516-156-31.txt", color);
    std::vector<unsigned int> indices;

    // Step 1: Convert 3D points to 2D points (using x, z)
    std::vector<Point_2> points2D;
    std::unordered_map<Point_2, unsigned int> pointIndexMap;

    size_t totalVertices = vertices.size();
    size_t processedVertices = 0;

    // Progress bar setup
    int barWidth = 50;
    std::cout << "adding points to Delaunay..." << std::endl;

    for (unsigned int i = 0; i < vertices.size(); ++i) {
        Point_2 p2d(vertices[i].x, vertices[i].z);
        points2D.push_back(p2d);
        pointIndexMap[p2d] = i;  // Store the original index

        // Update the progress bar every 1000 points or at the last point
        processedVertices++;
        if (processedVertices % 10000 == 0 || processedVertices == totalVertices) {
            float progress = static_cast<float>(processedVertices) / totalVertices;
            int pos = barWidth * progress;
            std::cout << "\r[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

    std::cout << std::endl;

    // Step 2: Perform Delaunay Triangulation on the 2D points
    Delaunay dt;
    dt.insert(points2D.begin(), points2D.end());

    // Step 3: Generate indices based on the triangulation
    size_t totalFaces = 0;
    for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
        totalFaces++;
    }

    size_t processedFaces = 0;
    std::cout << "Triangulating faces..." << std::endl;

    for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
        // Get the three vertices of the triangle in the 2D plane
        Point_2 p0 = f->vertex(0)->point();
        Point_2 p1 = f->vertex(1)->point();
        Point_2 p2 = f->vertex(2)->point();

        // Look up the indices of the original vertices in the pointIndexMap
        unsigned int idx0 = pointIndexMap[p0];
        unsigned int idx1 = pointIndexMap[p1];
        unsigned int idx2 = pointIndexMap[p2];

        // Add indices to the list (three vertices per face)
        indices.push_back(idx0);
        indices.push_back(idx1);
        indices.push_back(idx2);

        // Update the progress bar every 100 faces or at the last face
        processedFaces++;
        if (processedFaces % 10000 == 0 || processedFaces == totalFaces) {
            float progress = static_cast<float>(processedFaces) / totalFaces;
            int pos = barWidth * progress;
            std::cout << "\r[";
            for (int j = 0; j < barWidth; ++j) {
                if (j < pos) std::cout << "=";
                else if (j == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0f) << "%";
            std::flush(std::cout);
        }
    }

    std::cout << std::endl;

  //Compute normals for each vertex
    for (auto& vertex : vertices) {
        vertex.normalx = 0;
        vertex.normaly = 0;
        vertex.normalz = 0;
    }

    for (size_t i = 0; i < indices.size(); i += 3) {
        unsigned int idx0 = indices[i];
        unsigned int idx1 = indices[i + 1];
        unsigned int idx2 = indices[i + 2];

        Vertex& v0 = vertices[idx0];
        Vertex& v1 = vertices[idx1];
        Vertex& v2 = vertices[idx2];

        glm::vec3 p0(v0.x, v0.y, v0.z);
        glm::vec3 p1(v1.x, v1.y, v1.z);
        glm::vec3 p2(v2.x, v2.y, v2.z);

        glm::vec3 edge1 = p1 - p0;
        glm::vec3 edge2 = p2 - p0;

        glm::vec3 faceNormal = glm::normalize(glm::cross(edge1, edge2));

        v0.normalx += faceNormal.x;
        v0.normaly += faceNormal.y;
        v0.normalz += faceNormal.z;

        v1.normalx += faceNormal.x;
        v1.normaly += faceNormal.y;
        v1.normalz += faceNormal.z;

        v2.normalx += faceNormal.x;
        v2.normaly += faceNormal.y;
        v2.normalz += faceNormal.z;
    }

    for (auto& vertex : vertices) {
        glm::vec3 normal(vertex.normalx, vertex.normaly, vertex.normalz);
        normal = glm::normalize(normal);
        vertex.normalx = normal.x;
        vertex.normaly = normal.y;
        vertex.normalz = normal.z;
    }

    return { vertices, indices };
}









std::pair<std::vector<Vertex>, std::vector<unsigned int>> Mesh::BSplineSurface(glm::vec3 color) {
    std::vector<Vertex> vertices = Readfile("Data/32-2-516-156-31.txt", color);
    std::vector<unsigned int> indices;

    // Control points for the surface (example control grid)
    std::vector<glm::vec3> mc; // = Readfile("Data/32-2-516-156-31.txt", color);

    int pointcloudSize = 2531030; 
    pointcloudSize /= 1000; 

   const int numU =sqrt(pointcloudSize);
   const int numV = sqrt(pointcloudSize);
   const int degree = 2;


    // Corrected knot vectors
    std::vector<float> uKnots;
    std::vector<float> vKnots;

    // U Knot Vector
    for (int i = 0; i < (numU + degree + 1); i++) {
        if (i <= degree) {
            uKnots.push_back(0.0f);
        }
        else if (i >= numU) {
            uKnots.push_back(numU - degree);  // Normalize the end knot to 1
        }
        else {
            uKnots.push_back((float)(i - degree));  // Interior knots
        }
    }

    // V Knot Vector
    for (int i = 0; i < (numV + degree + 1); i++) {
        if (i <= degree) {
            vKnots.push_back(0.0f);
        }
        else if (i >= numV) {
            vKnots.push_back(numV - degree);  // Normalize the end knot to 1
        }
        else {
            vKnots.push_back(i-degree);  // Interior knots
        }
    }

	for (auto v : vertices)
	{
		mc.push_back(glm::vec3(v.x, v.y, v.z));
	}

	std::cout << "uKnots: ";
	for (auto u : uKnots)
	{
		std::cout << u << " ";
	}
	std::cout << std::endl;
	std::cout << "vKnots: ";
	for (auto v : vKnots)
	{
		std::cout << v << " ";
	}
	std::cout << std::endl;


    std::vector<std::vector<glm::vec3>> c(numU, std::vector<glm::vec3>(numV));
	for (int i = 0; i < numU; i++)
	{
		for (int j = 0; j < numV; j++)
		{
			c[i][j] = mc[j * numU + i];
		}
	}

	MakeBiquadraticSurface(numU, numV, degree, degree, uKnots, vKnots, mc);


    return { m_vertices, m_indices };
}

void Mesh::MakeBiquadraticSurface(const int n_u, const int n_v, int d_u, int d_v, std::vector<float> mu, std::vector<float> mv, std::vector<glm::vec3> mc)
{
    float h = 0.1f; // Spacing
    int nu = static_cast<int>((mu[n_u] - mu[d_u]) / h);  // Calculate the number of steps in u
    int nv = static_cast<int>((mv[n_v] - mv[d_v]) / h);  // Calculate the number of steps in v

    // Iterate through v and u to generate surface points
    for (int i = 0; i < nv; ++i)
    {
        for (int j = 0; j < nu; ++j)
        {
            float u = j * h;
            float v = i * h;

            // Find the corresponding knot intervals for u and v
            //int my_u = FindKnotInterval(mu, d_u, n_u, u);
            //int my_v = FindKnotInterval(mv, d_v, n_v, v);

            // Calculate the basis function coefficients for the current u and v
            //auto koeff_par = B2(u, v, my_u, my_v);

            // Evaluate the biquadratic surface at the current u and v
            glm::vec3 surfacePoint = deBoorSurface(d_u, d_v, mu, mv, mc, u, v, n_u, n_v);

            Vertex vertex;

            // Assign the position values from the evaluated surface point
            vertex.x = surfacePoint.x;
            vertex.y = surfacePoint.y;
            vertex.z = surfacePoint.z;

            // Assign color and texture coordinates
            vertex.r = 1.0f;
            vertex.g = 1.0f;
            vertex.b = 1.0f;

            vertex.u = static_cast<float>(j) / (nu - 1); // Column index normalized
            vertex.v = static_cast<float>(i) / (nv - 1); // Row index normalized

            vertex.normalx = 0.0f;
            vertex.normaly = 1.0f;
            vertex.normalz = 0.0f;

            // Push the computed surface point into the vertices array
            m_vertices.push_back(vertex);
        }
    }

    // Generate indices for the triangle mesh
    for (int i = 0; i < nv - 1; ++i) {
        for (int j = 0; j < nu - 1; ++j) {
            int idx1 = i * nu + j;
            int idx2 = idx1 + 1;
            int idx3 = idx1 + nu;
            int idx4 = idx3 + 1;

            // First triangle (idx1, idx2, idx3)
            m_indices.push_back(idx1);
            m_indices.push_back(idx2);
            m_indices.push_back(idx3);

            // Second triangle (idx2, idx4, idx3)
            m_indices.push_back(idx2);
            m_indices.push_back(idx4);
            m_indices.push_back(idx3);
        }
    }


    // Generate indices for the triangle mesh
    for (int i = 0; i < nv - 1; ++i) {
        for (int j = 0; j < nu - 1; ++j) {
            int idx1 = i * nu + j;
            int idx2 = idx1 + 1;
            int idx3 = idx1 + nu;
            int idx4 = idx3 + 1;

            // First triangle (idx1, idx2, idx3)
            m_indices.push_back(idx1);
            m_indices.push_back(idx2);
            m_indices.push_back(idx3);

            // Second triangle (idx2, idx4, idx3)
            m_indices.push_back(idx2);
            m_indices.push_back(idx4);
            m_indices.push_back(idx3);
        }
    }
}



glm::vec3 Mesh::deBoorSurface(int d_u, int d_v, std::vector<float> mu, std::vector<float> mv, std::vector<glm::vec3> ControlsPoints, float u, float v, const int KnotsUsize, const int KnotsVsize)
{
    // Step 1: Apply de Boor along the u-direction
    std::vector<glm::vec3> tempPoints;
    for (int i = 0; i < KnotsVsize; ++i) {
        // Extract the i-th row of control points (nu points per row)
        std::vector<glm::vec3> row;
        for (int j = 0; j < KnotsUsize; ++j) {
            row.push_back(ControlsPoints[i * KnotsUsize + j]);
        }

        // Apply de Boor along the u-direction for this row
        glm::vec3 t = deBoor(d_u, d_u, mu, row, u);
        tempPoints.push_back(t); // Store the result for v-direction interpolation
    }

    // Now apply de Boor along the v-direction on the result from the u-direction
    glm::vec3 f = deBoor(d_v, d_v, mv, tempPoints, v);
    return f; // Interpolate along v with the new points
}


glm::vec3 Mesh::deBoor(int k, int degree, const std::vector<float>& knots, std::vector<glm::vec3> controlPoints, float t)
{
    // Find the knot span index
    int span = -1;
    for (int i = degree; i < knots.size() - 1; ++i) {
        if (t >= knots[i] && t < knots[i + 1]) {
            span = i;
            break;
        }
    }

    if (span == -1) {
        std::cout << "Could not find knot span" << std::endl;
        return glm::vec3(0.0f); // Return a default value if not found
    }

    // Initialize d as a copy of control points for the relevant knot span
    std::vector<glm::vec3> d(degree + 1);
    for (int i = 0; i <= degree; ++i) {
        d[i] = controlPoints[span - degree + i];
    }

    // Perform de Boor recursion
    for (int r = 1; r <= degree; ++r) {
        for (int j = degree; j >= r; --j) {
            float alpha = (t - knots[span - degree + j]) / (knots[span + 1 + j - r] - knots[span - degree + j]);
            d[j] = (1.0f - alpha) * d[j - 1] + alpha * d[j];
        }
    }

    // The evaluated point is now stored in d[degree]
    return d[degree];
}

std::vector<Vertex> Mesh::Readfile(const char* fileName, glm::vec3 color) {
    std::ifstream inputFile(fileName);
    std::vector<Vertex> pointCloud;
    float min_x = -816.02, max_x = 783.98;
    float min_z = -620.771, max_z = 579.229;
    int processedLines = 0;
    int totalLines = 2531030;  // Total number of lines in the file
	int numPointsIncreacedFriction = 0;

    if (inputFile.is_open()) {
        std::string line;
        std::getline(inputFile, line);  // Skip header line if there is one
        Vertex point;


        // Progress bar setup
        int barWidth = 50;
        std::cout << "Loading points..." << std::endl;

        while (std::getline(inputFile, line))
        {
            if (sscanf_s(line.c_str(), "%f %f %f", &point.x, &point.z, &point.y) == 3 && processedLines % 200 == 1)
            {
                point.x -= 608016.02;
                point.y -= 336.8007;
                point.z -= 6750620.771;

                point.u = (point.x - min_x) / (max_x - min_x);  // Normalize x coordinate
                point.v = (point.z - min_z) / (max_z - min_z);  // Normalize z coordinate

				if (numPointsIncreacedFriction < 50)
				{
					point.friction = 0.9f;
					point.r = 1.0f;
					point.g = 0.0f;
					point.b = 0.0f;
					numPointsIncreacedFriction++;
				}
				else
				{
                    point.r = color.x;
                    point.g = color.y;
                    point.b = color.z;
					point.friction = 0.1f;
					numPointsIncreacedFriction++;
				}
				if (numPointsIncreacedFriction == 300)
				{
					numPointsIncreacedFriction = 0;
				}

                pointCloud.push_back(point);
            }

            processedLines++;
            if (processedLines % 100000 == 0 || processedLines == totalLines) {
                float progress = static_cast<float>(processedLines) / totalLines;
                int pos = barWidth * progress;
                std::cout << "\r[";  // Start overwriting the same line
                for (int i = 0; i < barWidth; ++i) {
                    if (i < pos) std::cout << "=";
                    else if (i == pos) std::cout << ">";
                    else std::cout << " ";
                }
                std::cout << "] " << int(progress * 100.0f) << "%";
                std::flush(std::cout);
            }
        }

        // Ensure the progress bar is at 100% when done
        float progress = 100;
        int pos = barWidth * progress;
        std::cout << "\r[";  // Start overwriting the same line
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress) << "%";
        std::flush(std::cout);

        // After the progress bar is complete, print a newline
        std::cout << std::endl;
        inputFile.close();
    }
    else {
        std::cerr << "Unable to open the input file for reading." << std::endl;
    }

    return pointCloud;
}

float Mesh::BSplineBasis(int i, int p, float t, const std::vector<float>& knots) {
    // Out-of-range check
    if (i < 0 || i >= knots.size() - 1)
    {
		std::cerr << "Index out of range!" << std::endl;
        return 0.0f;
    }

    // Base case: zero-degree basis function
    if (p == 0) {
        // Handle special case for the last knot span
        if (i == knots.size() - 2) {
            return (t >= knots[i] && t <= knots[i + 1]) ? 1.0f : 0.0f;
        }
        return (t >= knots[i] && t < knots[i + 1]) ? 1.0f : 0.0f;
    }

    // Recursive case: compute alpha and beta
    float alpha = 0.0f;
    if (knots[i + p] != knots[i]) {
        alpha = (t - knots[i]) / (knots[i + p] - knots[i]) * BSplineBasis(i, p - 1, t, knots);
    }

    float beta = 0.0f;
    if (knots[i + p + 1] != knots[i + 1]) {
        beta = (knots[i + p + 1] - t) / (knots[i + p + 1] - knots[i + 1]) * BSplineBasis(i + 1, p - 1, t, knots);
    }

    return alpha + beta;
}





std::vector<glm::vec3> Mesh::BarycentricCoordinates(std::vector<unsigned int> indices, std::vector<Vertex> vertices)
{
    std::vector<glm::vec3> result;
	float u, v, w;
    if (vertices.empty()) {
		std::cerr << "Vertices list is empty!" << std::endl;
		return std::vector<glm::vec3>();
    }
    for (int i = 0; i < indices.size(); i += 3) 
    {
        int index0 = indices[i];
        int index1 = indices[i + 1];
        int index2 = indices[i + 2];
		glm::vec3 v0(vertices[index0].x, vertices[index0].y, vertices[index0].z);
		glm::vec3 v1(vertices[index1].x, vertices[index1].y, vertices[index1].z);
		glm::vec3 v2(vertices[index2].x, vertices[index2].y, vertices[index2].z);

        glm::vec3 cpoint = (v0 + v1 + v2) / 3.0f; // get center of triangle 

        glm::vec3 v0v1 = v1 - v0;
        glm::vec3 v0v2 = v2 - v0;
        glm::vec3 v0p = cpoint - v0;

        // Computing dot products
        double dot00 = glm::dot(v0v1, v0v1);
        double dot01 = glm::dot(v0v1, v0v2);
        double dot02 = glm::dot(v0v1, v0p);
        double dot11 = glm::dot(v0v2, v0v2);
        double dot12 = glm::dot(v0v2, v0p);

        // Computing barycentric coordinates
        double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
        double v = (dot11 * dot02 - dot01 * dot12) * invDenom;
        double w = (dot00 * dot12 - dot01 * dot02) * invDenom;
        double u = 1 - v - w;
		result.push_back(glm::vec3(u, v, w));
    }



    return result;
}




