#pragma once
#include <vector>
#include <glm/glm.hpp>
#include <memory>

class Entity;

struct Cell {
public: 
    std::vector<std::shared_ptr<Entity>> entities;
	std::vector<unsigned int> groundIndices;
};

class Grid
{
    //friend class Collision;
public:
    Grid(int width, int height, int cellSize);
    ~Grid();
    //|-----------------------------------------------------------------------------|
    //|                                Public Functions                             |
    //|-----------------------------------------------------------------------------|
    void AddBaLL(std::shared_ptr<Entity> entity);
    void AddBaLL(std::shared_ptr<Entity>, Cell* entity);
    void RemoveBallFromCell(std::shared_ptr<Entity> entity, Cell* cell);
    //|-----------------------------------------------------------------------------|
    //|                                Getters                                      |
    //|-----------------------------------------------------------------------------|

    Cell* getCell(int x, int y);
    Cell* getCell(const glm::vec3& pos);

    int m_numXCells;
    int m_numYCells;
    int m_cellSize;

private:
    //|-----------------------------------------------------------------------------|
    //|                                Private variables                            |
    //|-----------------------------------------------------------------------------|
    std::vector<Cell> m_cells;
    int m_width;
    int m_height;
};
