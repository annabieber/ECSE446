/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core
#define PI 3.14159265358979323846f

uniform vec3 cameraPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 albedo;

in vec3 vPos;
in vec3 vNormal;

out vec3 color;

void main()
{
	color = vec3(0.f);
    // TODO(A2): Implement this
	//need to find cosWi
	vec3 wi = lightPos - vPos;
	vec3 wo = normalize(cameraPos - vPos);
	float cosWi = dot(vNormal, normalize(wi));
	float cosWo = dot(wo, vNormal);

	//need to check that cosWi and cosWo bigger than zero
	if (cosWi > 0) {
		//multiply by Li 
		vec3 light = lightIntensity / dot(wi, wi);
		color = light * (albedo / PI) * cosWi;
	}
	
}
