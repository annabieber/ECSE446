/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#version 330 core

uniform vec3 camPos;
uniform vec3 lightPos;
uniform vec3 lightIntensity;
uniform vec3 rho_d;
uniform vec3 rho_s;
uniform float exponent;

in vec3 vPos;
in vec3 vNormal;

out vec3 color;

#define INV_PI 0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f

void main() {
	
	//set the color to black (or absence of light)
	color = vec3(0.f);
	//wi and wo are calculated with the light and camera positions and the intersection point 
	vec3 wi = lightPos - vPos;
	vec3 wo = camPos - vPos;
	//obtain cosWi and cosWo by doing the dot product of the normal and the wi/wo
	float cosWi = dot(vNormal, normalize(wi));
	float cosWo = dot(vNormal, normalize(wo)); 
	 
	 //check that both are positif and calulate the color
	if(cosWi > 0 && cosWo > 0){
		//obtain the specular direction (wr) by reflecting wi and the normal, then get alpha by doing the dot prodct of wr and wo
		vec3 specularDirection = reflect(normalize(wi), normalize(vNormal));
		float alpha = dot(specularDirection, normalize(wo));
		float maxValue = max(0.f,pow(alpha, exponent));
		//calculate Li 
		vec3 light = lightIntensity / dot(wi, wi);
		//calculate the color and multiply by Li and cosWi
		color = light * (rho_d * INV_PI + rho_s * (exponent + 2) * INV_TWOPI * maxValue) * cosWi;
	}

	

}