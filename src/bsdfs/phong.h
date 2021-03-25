/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
	struct PhongBSDF : BSDF {

	std::unique_ptr<Texture < v3f>> specularReflectance;
	std::unique_ptr<Texture < v3f>> diffuseReflectance;
	std::unique_ptr<Texture < float>> exponent;
	float specularSamplingWeight;
	float scale;

	PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
		const tinyobj::material_t& mat = scene.materials[matID];

		if (mat.specular_texname.empty())
			specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
		else
			specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

		if (mat.diffuse_texname.empty())
			diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
		else
			diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

		exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

		//get scale value to ensure energy conservation
		v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
		float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
		scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

		float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
		float sAvg = getLuminance(specularReflectance->getAverage() * scale);
		specularSamplingWeight = sAvg / (dAvg + sAvg);

		components.push_back(EGlossyReflection);
		components.push_back(EDiffuseReflection);

		combinedType = 0;
		for (unsigned int component : components)
			combinedType |= component;
	}

	inline float getExponent(const SurfaceInteraction& i) const override {
		return exponent->eval(worldData, i);
	}

	inline v3f reflect(const v3f& d) const {
		return v3f(-d.x, -d.y, d.z);
	}

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);

		// TODO(A2): Implement this

		//obtain cosWo, cosWi, rhoD, rhoS and n (exponent)
		float cosWo = Frame::cosTheta(glm::normalize(i.wo));
		float cosWi = Frame::cosTheta(glm::normalize(i.wi));
		v3f rhoD = diffuseReflectance->eval(worldData, i) * PhongBSDF::scale; //scale to ensure energy conservation
		v3f rhoS = specularReflectance->eval(worldData, i) * PhongBSDF::scale; //scale to ensure energy conservation
		float n = exponent->eval(worldData, i);


		//check that both are positif and calulate the color
		if (cosWi > 0 && cosWo > 0) {
			//alpha can also be obtained by getting the specular direction with i.wo and doing the dot production with i.wi
			v3f specularDirection = (PhongBSDF::reflect(glm::normalize(i.wi)));
			float alpha = glm::dot(specularDirection, glm::normalize(i.wo));
			float max = glm::max(0.f, pow(alpha, n));
			//calculate the value with the phong BRDF
			val = ((rhoD * INV_PI) + (rhoS * (n + 2) * INV_TWOPI * max)) * cosWi;
		}
		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;

		// TODO(A3): Implement this
		float n = exponent->eval(worldData, i);
		//get i.wi in local coordinates
		Frame f = Frame(reflect(i.wo));
		v3f wi = f.toLocal(i.wi);
		//compute the probability density of the mixture
		pdf = specularSamplingWeight * Warp::squareToPhongLobePdf(wi, n) + (1-specularSamplingWeight) * Warp::squareToCosineHemispherePdf(i.wi);
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, Sampler& sampler, float* pdf) const override {
		v3f val(0.f);

		// TODO(A3): Implement this
		float randomNumber = sampler.next();

		//pick a random number, if less than the specular weight do phong sampling
		//else do diffuse sampling
		//this allows for randomness in the sampling methods
		if (randomNumber < specularSamplingWeight) {
			//get sample direction and sample pdf
			float n = exponent->eval(worldData, i);
			v3f sampleDirection = Warp::squareToPhongLobe(sampler.next2D(), n);

			//create frame with reflected wo direction and set wi in world coordinates
			Frame f = Frame(reflect(i.wo));
			i.wi = f.toWorld(sampleDirection);

		}
		else {
			// get sample direction and sample pdf
			v3f sampleDirection = Warp::squareToCosineHemisphere(sampler.next2D());
			//assign the sampling diretion to the surface 
			i.wi = sampleDirection;
		}

		*pdf = PhongBSDF::pdf(i);
		//Evaluate the surface interaction (phong or diffuse)
		if (pdf != 0) {
			val = eval(i);
			return val;
		}
		else {
			return val;
		}
	}

	std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END