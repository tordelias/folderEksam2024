#include "Mesh.h"
#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include <chrono>
#include <thread>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/convex_hull_2.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Euclidean_distance.h>
//#include <CGAL/IO/Polyhedron_iostream.h>
//
//
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Triangulation_vertex_base_2<K> BaseVertex;
//typedef CGAL::Triangulation_vertex_base_2<K, CGAL::Triangulation_ds_vertex_base_2<BaseVertex>> VertexBase;
//typedef CGAL::Triangulation_data_structure_2<VertexBase> Tds;
//typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
//
//// Define a custom vertex with index
//typedef Delaunay::Vertex_handle Vertex_handle_with_index;
//typedef K::Point_2 Point_2;

void showLoadingBar(size_t current, size_t total, const std::string& taskName = "", int barWidth = 50) {
    float progress = static_cast<float>(current) / total;
    int pos = barWidth * progress;
    std::cout << "\r" << taskName << " [";
    for (int j = 0; j < barWidth; ++j) {
        if (j < pos) std::cout << "=";
        else if (j == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0f) << "%";
    std::flush(std::cout);
}
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
            Vertex V1 = Vertex{ x * radius, y * radius, z * radius, color.r, color.g, color.b, u, v };
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
                vertex.friction = 0.0f;
                vertex.r = 1.f;
				vertex.g = 0.f;
				vertex.b = 0.f;
            }
			else
                vertex.friction = 0.0f;
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
    // Read the point cloud with full Vertex data
    std::vector<Vertex> vertices = Readfile("Data/SmallData_centered.txt", color);
    std::vector<unsigned int> indices;
    std::vector<glm::vec3> vertexNormals(vertices.size(), glm::vec3(0.0f));

    // Create a mapping from Point_2 to vertex indices using std::map for compatibility with CGAL
    std::map<Point_2, unsigned int> point_to_index;
    std::vector<Point_2> points;
    points.reserve(vertices.size());

    for (unsigned int i = 0; i < vertices.size(); ++i) {
        Point_2 pt(vertices[i].x, vertices[i].z); // Using x and z as 2D coordinates
        points.push_back(pt);
        point_to_index[pt] = i;
    }

    // Perform Delaunay triangulation
    Delaunay delaunay;
    delaunay.insert(points.begin(), points.end());

    // Reserve space for indices upfront to minimize reallocations
    indices.reserve(delaunay.number_of_faces() * 3);

    // Compute normals in parallel
#pragma omp parallel
    {
        std::vector<glm::vec3> localNormals(vertices.size(), glm::vec3(0.0f));
        std::vector<unsigned int> localIndices;

#pragma omp for nowait
        for (auto face = delaunay.finite_faces_begin(); face != delaunay.finite_faces_end(); ++face) {
            // Retrieve vertex indices
            unsigned int i0 = point_to_index[face->vertex(0)->point()];
            unsigned int i1 = point_to_index[face->vertex(1)->point()];
            unsigned int i2 = point_to_index[face->vertex(2)->point()];

            // Store triangle indices
            localIndices.insert(localIndices.end(), { i0, i1, i2 });

            // Calculate the normal vector for the triangle
            const glm::vec3 p0(vertices[i0].x, vertices[i0].y, vertices[i0].z);
            const glm::vec3 p1(vertices[i1].x, vertices[i1].y, vertices[i1].z);
            const glm::vec3 p2(vertices[i2].x, vertices[i2].y, vertices[i2].z);

            glm::vec3 edge1 = p1 - p0;
            glm::vec3 edge2 = p2 - p0;
            glm::vec3 normal = glm::cross(edge1, edge2);

            // Handle degenerate triangles
            if (glm::dot(normal, normal) < 0.f) {
                // Skip degenerate triangles (too small or collinear)
                continue;
            }


            // Accumulate normals
            localNormals[i0] += normal;
            localNormals[i1] += normal;
            localNormals[i2] += normal;
        }

        // Merge local results into the global data structures
#pragma omp critical
        {
            indices.insert(indices.end(), localIndices.begin(), localIndices.end());
            for (unsigned int i = 0; i < vertices.size(); ++i) {
                vertexNormals[i] += localNormals[i];
            }
        }
    }

    // Normalize accumulated normals and assign them to vertices
#pragma omp parallel for
    for (int i = 0; i < static_cast<int>(vertices.size()); ++i) {
        glm::vec3& normal = vertexNormals[i];
        if (glm::dot(normal, normal) > 1e-12f) {
            normal = glm::normalize(normal);
        }
        else {
            normal = glm::vec3(0.0f, 1.0f, 0.0f); // Fallback for degenerate cases
        }

        vertices[i].normalx = normal.x;
        vertices[i].normaly = normal.y;
        vertices[i].normalz = normal.z;
    }
    flipEdgeIfNecessary(vertices, delaunay, point_to_index);
    return std::make_pair(std::move(vertices), std::move(indices));
}


// flipEdgeIfNecessary(vertices, delaunay, point_to_index);

void Mesh::laplacianSmoothing(std::vector<Vertex>& vertices, std::vector<std::vector<unsigned int>>& vertexNeighbors, float lambda, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<glm::vec3> newPositions(vertices.size(), glm::vec3(0.0f));

        // Parallelize Laplacian smoothing
#pragma omp parallel for
        for (unsigned int i = 0; i < vertices.size(); ++i) {
            glm::vec3 smoothedPosition(0.0f);
            unsigned int num_neighbors = vertexNeighbors[i].size();

            for (unsigned int neighbor : vertexNeighbors[i]) {
                smoothedPosition += glm::vec3(vertices[neighbor].x, vertices[neighbor].y, vertices[neighbor].z);
            }

            if (num_neighbors > 0) smoothedPosition /= num_neighbors;
            glm::vec3 originalPosition(vertices[i].x, vertices[i].y, vertices[i].z);
            newPositions[i] = originalPosition + lambda * (smoothedPosition - originalPosition);
        }

        // Update vertex positions
        for (unsigned int i = 0; i < vertices.size(); ++i) {
            vertices[i].x = newPositions[i].x;
            vertices[i].y = newPositions[i].y;
            vertices[i].z = newPositions[i].z;
        }
    }
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
        Vertex point;

        point.r = color.x;
        point.g = color.y;
        point.b = color.z;


        if (std::getline(inputFile, line))
                totalLines = std::stoi(line);

		std::cout << "Total number of Points: " << totalLines << std::endl;
		std::cout << "Processing Points" << std::endl;

        while (std::getline(inputFile, line))
        {
            if (sscanf_s(line.c_str(), "%f %f %f", &point.x, &point.z, &point.y) == 3 )
            {



                point.u = (point.x - min_x) / (max_x - min_x);  // Normalize x coordinate
                point.v = (point.z - min_z) / (max_z - min_z);  // Normalize z coordinate

				//if (numPointsIncreacedFriction < 500)
				//{
				//	point.friction = 0.9f;
				//	point.r = 1.0f;
				//	point.g = 0.0f;
				//	point.b = 0.0f;
				//	numPointsIncreacedFriction++;
				//}
				//else
				//{
    //                point.r = color.x;
    //                point.g = color.y;
    //                point.b = color.z;
				//	point.friction = 0.f;
				//	numPointsIncreacedFriction++;
				//}
				//if (numPointsIncreacedFriction  >= totalLines % 100)
				//{
				//	numPointsIncreacedFriction = 0;
				//}

                pointCloud.push_back(point);
            }
            processedLines++;
            if (processedLines % 100000 == 0 || processedLines == totalLines) 
                showLoadingBar(processedLines, totalLines);

        }
        showLoadingBar(totalLines, totalLines, "");
        std::cout << std::endl; 
        inputFile.close();
    }
    else {
        std::cerr << "Unable to open the input file for reading." << std::endl;
    }
    return pointCloud;
}


void Mesh::flipEdgeIfNecessary(std::vector<Vertex>& vertices, Delaunay& triangulation, std::map<Point_2, unsigned int>& point_to_index) {
    // Loop over all faces in the triangulation
    for (auto face = triangulation.finite_faces_begin(); face != triangulation.finite_faces_end(); ++face) {
        // Dereference the iterator to get the Face_handle
        Delaunay::Face_handle faceHandle = face;

        // Get the vertices of the current triangle
        Point_2 v0 = faceHandle->vertex(0)->point();
        Point_2 v1 = faceHandle->vertex(1)->point();
        Point_2 v2 = faceHandle->vertex(2)->point();

        // Convert points to vertex positions (assuming point_to_index maps Point_2 -> index of the vertex)
        glm::vec3 p0(vertices[point_to_index[v0]].x, vertices[point_to_index[v0]].y, vertices[point_to_index[v0]].z);
        glm::vec3 p1(vertices[point_to_index[v1]].x, vertices[point_to_index[v1]].y, vertices[point_to_index[v1]].z);
        glm::vec3 p2(vertices[point_to_index[v2]].x, vertices[point_to_index[v2]].y, vertices[point_to_index[v2]].z);

        // If the triangle has poor quality, we need to flip its edges
        if (isPoorQuality(p0, p1, p2)) {
            // Check the shared edges and find the neighbor
            for (unsigned int i = 0; i < 3; ++i) {
                unsigned int idx1 = i;
                unsigned int idx2 = (i + 1) % 3;

                // Get the neighboring face for this edge (i)
                auto neighbor = faceHandle->neighbor(i);  // `neighbor(i)` returns a Face_handle

                // Check if the neighboring face is valid and not infinite (boundary)
                if (neighbor != triangulation.infinite_face()) {
                    // Get the vertices of the neighboring face
                    Point_2 v3 = neighbor->vertex(0)->point();
                    Point_2 v4 = neighbor->vertex(1)->point();
                    Point_2 v5 = neighbor->vertex(2)->point();

                    // Check if the edge is shared (one of the vertices should be the same)
                    if ((v0 == v3 || v0 == v4 || v0 == v5) && (v1 == v3 || v1 == v4 || v1 == v5)) {
                        // Perform edge flip: Replace one triangle with its flipped version
                        flipTriangleEdges(faceHandle, neighbor); // Pass Face_handle directly (no need to dereference)


                        break;
                    }
                }
            }
        }
    }
}


void Mesh::flipTriangleEdges(Delaunay::Face_handle& triangle1, Delaunay::Face_handle& triangle2) {
    // Flip the shared edge between the two triangles
    Point_2 v0_1 = triangle1->vertex(0)->point();
    Point_2 v1_1 = triangle1->vertex(1)->point();
    Point_2 v2_1 = triangle1->vertex(2)->point();

    Point_2 v0_2 = triangle2->vertex(0)->point();
    Point_2 v1_2 = triangle2->vertex(1)->point();
    Point_2 v2_2 = triangle2->vertex(2)->point();

    // Find the shared edge between the two triangles
    std::vector<Point_2> sharedEdge = findSharedEdge(v0_1, v1_1, v2_1, v0_2, v1_2, v2_2);

    if (sharedEdge.size() == 2) {
        Point_2 sharedV1 = sharedEdge[0];
        Point_2 sharedV2 = sharedEdge[1];

        // Flip the edge (swap the diagonals)
        if (sharedV1 == v0_1 && sharedV2 == v1_1) {
            // Swap the diagonal between the two triangles
            swapEdge(v0_2, v2_2, triangle1, triangle2);
        }
        else if (sharedV1 == v1_1 && sharedV2 == v2_1) {
            swapEdge(v0_2, v1_2, triangle1, triangle2);
        }
        else if (sharedV1 == v2_1 && sharedV2 == v0_1) {
            swapEdge(v1_2, v2_2, triangle1, triangle2);
        }
    }
}

std::vector<Point_2> Mesh::findSharedEdge(Point_2 v0_1, Point_2 v1_1, Point_2 v2_1, Point_2 v0_2, Point_2 v1_2, Point_2 v2_2) {
    std::vector<Point_2> sharedEdge;
    if ((v0_1 == v0_2 || v0_1 == v1_2 || v0_1 == v2_2) &&
        (v1_1 == v0_2 || v1_1 == v1_2 || v1_1 == v2_2)) {
        sharedEdge.push_back(v0_1);
        sharedEdge.push_back(v1_1);
    }
    if ((v1_1 == v0_2 || v1_1 == v1_2 || v1_1 == v2_2) &&
        (v2_1 == v0_2 || v2_1 == v1_2 || v2_1 == v2_2)) {
        sharedEdge.push_back(v1_1);
        sharedEdge.push_back(v2_1);
    }
    return sharedEdge;
}

void Mesh::swapEdge(Point_2 v0, Point_2 v1, Delaunay::Face_handle& triangle1, Delaunay::Face_handle& triangle2) {
    // Perform edge swap by updating the Delaunay triangulation
    // Swap the diagonals between the two triangles
    triangle1->set_vertex(0, triangle2->vertex(1));
    triangle2->set_vertex(1, triangle1->vertex(2));
    triangle2->set_vertex(2, triangle1->vertex(0));

    // Update vertex positions and recompute face normals as necessary
}


bool Mesh::isPoorQuality(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2) {
    glm::vec3 edge1 = v1 - v0;
    glm::vec3 edge2 = v2 - v0;
    glm::vec3 edge3 = v2 - v1;

    // Check for large angles (indicating a very long or thin triangle)
    float angle1 = angleBetween(edge1, edge2);
    float angle2 = angleBetween(edge2, edge3);
    float angle3 = angleBetween(edge3, edge1);

    return (angle1 < 0.1f || angle2 < 0.1f || angle3 < 0.1f); // Threshold for "bad" triangles
}

float Mesh::angleBetween(glm::vec3 v1, glm::vec3 v2) {
    return acos(glm::dot(glm::normalize(v1), glm::normalize(v2)));
}









