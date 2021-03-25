/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core


in vec3 vNormal;

out vec3 color;

void main()
{
    // TODO(A1): Implement this
	//the value of the color is the absolute value of the normal calculate in the vertex shader
	color = abs(vNormal);
}
