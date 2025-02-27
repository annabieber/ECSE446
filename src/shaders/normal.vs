/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 normalMat;

layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

out vec3 vNormal;

void main()
{
    // TODO(A1): Implement this
	//transform triangles to match camera view 
	gl_Position = projection * view * model * vec4(position, 1.0);
	//return the normal to the fragment shader
	vNormal = normal;
}
