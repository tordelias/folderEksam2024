#include "Grid.h"
#include "iostream"
#include "../Component/Component.h"
#include "../Entity.h"

Grid::Grid(int width, int height, int cellSize)
    : m_width(width), m_height(height), m_cellSize(cellSize)
{
    // Adjust to cover negative and positive space
    m_numXCells = ceil(static_cast<float>(m_width) / cellSize);
    m_numYCells = ceil(static_cast<float>(m_height) / cellSize);

    // Double the size to account for negative and positive axes
    m_cells.resize((m_numYCells * 2) * (m_numXCells * 2));
}

Grid::~Grid()
{
}

void Grid::AddBaLL(std::shared_ptr<Entity> entity)
{
    Cell* cell = getCell(entity->GetComponent<TransformComponent>()->position);
    if (cell)
    {
        cell->entities.push_back(entity);
        entity->ownerCell = cell;
        entity->cellvectorindex = cell->entities.size() - 1;
    }
}

void Grid::AddBaLL(std::shared_ptr<Entity> entity, Cell* cell)
{
    if (cell)
    {
        cell->entities.push_back(entity);
        entity->ownerCell = cell;
        entity->cellvectorindex = cell->entities.size() - 1;
    }
}

Cell* Grid::getCell(int x, int y)
{
    // Offset x and y by the center of the grid
    int offsetX = x + (m_numXCells);
    int offsetY = y + (m_numYCells);

    if (offsetX < 0 || offsetX >= m_numXCells * 2 || offsetY < 0 || offsetY >= m_numYCells * 2)
    {
        std::cerr << "Grid index out of bounds!  -GetCell" << std::endl;
        return nullptr;
    }

    return &m_cells[offsetY * (2 * m_numXCells) + offsetX];
}

Cell* Grid::getCell(const glm::vec3& pos)
{
    // Calculate cell indices and account for negative coordinates
    int cellX = static_cast<int>(std::floor(pos.x / m_cellSize));
    int cellY = static_cast<int>(std::floor(pos.z / m_cellSize));

    return getCell(cellX, cellY);
}

void Grid::RemoveBallFromCell(std::shared_ptr<Entity> entity, Cell* cell)
{
    if (cell)
    {
        auto& balls = cell->entities;
        balls.erase(std::remove(balls.begin(), balls.end(), entity), balls.end());
    }
}
