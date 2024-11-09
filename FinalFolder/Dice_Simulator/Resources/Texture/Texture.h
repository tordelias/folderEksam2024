#pragma once
#include "stb/stb_image.h"
#include"../Shaders/shaderClass.h"
#include <string>
class Texture
{
	int ID; 
public:
	unsigned int texture;

	Texture(const char* texture1, const std::shared_ptr<Shader>& shaderProgram);

	
};

