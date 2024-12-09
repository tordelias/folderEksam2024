#pragma once
#include <glm/glm.hpp>
#include "Mesh.h"
#include "../System/ParticleSystem.h"
#include <string>

class Component
{
public:
    virtual ~Component() = default;
};

struct TransformComponent : public Component
{
    glm::vec3 position;
    glm::vec3 rotation;
	glm::vec3 velocity;
    glm::vec3 angularVelocity; 
    glm::vec3 scale;

    TransformComponent(const glm::vec3& pos = glm::vec3(0.0f),
        const glm::vec3& rot = glm::vec3(0.0f),
        const glm::vec3& scl = glm::vec3(1.0f))
        : position(pos), rotation(rot), scale(scl), velocity(glm::vec3(0, 0, 0)), angularVelocity(glm::vec3(0, 0, 0)) {}
};

struct MeshComponent : public Component {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;
    const char* TexturePath;

    MeshComponent(const char* figure = "", const glm::vec3& color = glm::vec3(1,1,1), const char* texture = "");
};
class Entity;
struct SplineComponent : public Component
{
public:

	std::vector<Vertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<glm::vec3> controllpoints;
    std::shared_ptr<Entity> owner; 
    std::shared_ptr<Entity> spline;
    glm::vec3 startpos; 
	SplineComponent(std::shared_ptr<Entity> ownerptr);
    void CalculateBSpline();
};

struct ParticleComponent : public Component
{
public:
	glm::vec3 position;
	glm::vec3 velocity;
	glm::vec3 acceleration;
	glm::vec3 Size;
	ParticleComponent(const glm::vec3& pos = glm::vec3(0.0f),
		const glm::vec3& vel = glm::vec3(0.f),
		const glm::vec3& acc = glm::vec3(0.f, 9.81f, 0.f),
		const glm::vec3& size = glm::vec3(1.f))
		: position(pos), velocity(vel), acceleration(acc), Size(size)
	{
		particleSystem = std::make_shared<ParticleSystem>(position, acceleration, Size, glm::vec3(0.1f), glm::ivec2(0, 2), glm::vec3(0.5, 0.5, 0.5), 50, 1);
	}
	std::shared_ptr<ParticleSystem> particleSystem;
};

extern "C" 
{
#include "../lua54/include/lua.h"
#include "../lua54/include/lauxlib.h"
#include "../lua54/include/lualib.h"
}

#ifdef _WIN32
#pragma comment(lib, "lua54/lua54.lib")
#endif

#include <iostream>>
#include <filesystem>
#include <memory>

class LuaComponent : public Component {
public:
    lua_State* L;                  // Lua state
    std::shared_ptr<ParticleSystem> particleSystem; // ParticleSystem instance
    std::string luaFilePath;       // Path to the Lua script file
	std::filesystem::file_time_type lastModifiedTime;  // Last modified time of the Lua file

public:
    // Constructor initializes the Lua state, registers functions, and handles Lua file
    LuaComponent(std::shared_ptr<ParticleSystem> system, const std::string& luaFile)
        : L(luaL_newstate()), particleSystem(system), luaFilePath(luaFile)
    {
        luaL_openlibs(L);  // Open Lua libraries
        lastModifiedTime = {};  // Initialize last modified time to zero

        if (luaL_dofile(L, luaFilePath.c_str()) != LUA_OK) {
            std::cerr << "Error loading Lua file: " << lua_tostring(L, -1) << std::endl;
            lua_pop(L, 1);  // Pop the error message from the Lua stack
        }
        else {
            RegisterFunctions();  // Register functions after Lua script is loaded
        }
    }

    ~LuaComponent() {
        if (L) {
            lua_close(L);  // Close the Lua state when done
        }
    }

    void RegisterFunctions() {
        // Push 'this' as light userdata
        lua_pushlightuserdata(L, this);

        // Push the C++ function with upvalue
        lua_pushcclosure(L, LuaComponent::updateAccelerationWrapper, 1);

        // Set as a global Lua function
        lua_setglobal(L, "updateAcceleration");

        // Debug: Check if registered properly
        lua_getglobal(L, "updateAcceleration");
        if (lua_isfunction(L, -1)) {
            //std::cout << "updateAcceleration function registered successfully in Lua state." << std::endl;
        }
        else {
            //std::cerr << "Error: updateAcceleration not found in Lua state!" << std::endl;
        }
        lua_pop(L, 1); // Pop the checked value from the stack
    }



    static int updateAccelerationWrapper(lua_State* L) {
        // Retrieve 'this' pointer from upvalue
        LuaComponent* luaComponent = static_cast<LuaComponent*>(lua_touserdata(L, lua_upvalueindex(1)));

        if (!luaComponent) {
            std::cerr << "Error: luaComponent is nullptr!" << std::endl;
            return 0; // Early exit
        }

        // Ensure we have 3 arguments (x, y, z)
        if (lua_gettop(L) != 3) {
            luaL_error(L, "Expected 3 arguments (x, y, z) for updateAcceleration.");
            return 0;
        }

        // Get the new acceleration values from Lua arguments
        glm::vec3 newAcc(
            luaL_checknumber(L, 1),  // x
            luaL_checknumber(L, 2),  // y
            luaL_checknumber(L, 3)   // z
        );

        // Call the instance method
        luaComponent->updateAcceleration(newAcc);

        return 0; // No return values
    }


    // Update acceleration
    void updateAcceleration(const glm::vec3& newAcc) {
        if (particleSystem) {
            particleSystem->updateacceleration(newAcc);  // Call the ParticleSystem's update method
            std::cout << "Acceleration updated to: (" << newAcc.x << ", " << newAcc.y << ", " << newAcc.z << ")\n";
        }
    }

    // Function to check if the Lua file has been modified and reload it
    void CheckAndReloadLuaFile() {
        try {
            // Get the last write time using std::filesystem
            auto currentModifiedTime = std::filesystem::last_write_time(luaFilePath);

            // Compare the current file's last write time with the saved one
            if (currentModifiedTime != lastModifiedTime) {
                // File has changed, reload it
                std::cout << "Lua file has been modified. Reloading..." << std::endl;
                ReloadLuaFile();
                lastModifiedTime = currentModifiedTime;  // Update the last modified time
            }
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error checking Lua file: " << e.what() << std::endl;
        }
    }


    void ReloadLuaFile() {
        // Clear Lua state
        lua_close(L);
        L = luaL_newstate();
        luaL_openlibs(L);

        // Register functions and reload the Lua script
        RegisterFunctions();

        if (luaL_dofile(L, luaFilePath.c_str()) != LUA_OK) {
            std::cerr << "Error loading Lua file: " << lua_tostring(L, -1) << std::endl;
            lua_pop(L, 1);
        }
        else {
            std::cout << "Lua file loaded successfully: " << luaFilePath << std::endl;
        }
    }


};

