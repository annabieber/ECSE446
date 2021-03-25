/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Ambient occlusion integrator
 */
	struct AOIntegrator : Integrator {

	// Use this in your switch statement to select the sampling type 
	ESamplingType m_samplingStrategy;

	explicit AOIntegrator(const Scene& scene) : Integrator(scene) {
		m_samplingStrategy = scene.config.integratorSettings.ao.sampling_type;
	}

	v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		/*
		Use the m_sampling_type variable to set wi and the corresponding pdf
		appropriately for sphere, hemisphere, or cosine sampling.

		You can use a switch statement or an if/else block.

		The m_sampling_type variable is an enum. The different values of the enum
		can be accessed through:
		ESamplingType::ESpherical
		ESamplingType::EHemispherical
		ESamplingType::ECosineHemispherical
		*/

		// TODO(A3): Implement this

		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		if (intersection) {
			v3f sampleDirection;
			float sampledPdf;

			//get sampleDirection and pdf depending on the sampling strategy
			if (m_samplingStrategy == ESamplingType::ESpherical) {
				sampleDirection = glm::normalize(Warp::squareToUniformSphere(sampler.next2D()));
				sampledPdf = Warp::squareToUniformSpherePdf();
			}
			else if (m_samplingStrategy == ESamplingType::EHemispherical) {
				sampleDirection = glm::normalize(Warp::squareToUniformHemisphere(sampler.next2D()));
				sampledPdf = Warp::squareToUniformHemispherePdf(sampleDirection);
			}
			else if (m_samplingStrategy == ESamplingType::ECosineHemispherical) {
				sampleDirection = glm::normalize(Warp::squareToCosineHemisphere(sampler.next2D()));
				sampledPdf = Warp::squareToCosineHemispherePdf(sampleDirection);
			}

			//set wi in world coordinates and calculate cosTheta
			tempSurface.wi = tempSurface.frameNs.toWorld(sampleDirection);
			float cosTheta = glm::dot(tempSurface.frameNs.n, tempSurface.wi);
			float albedo = 1.0;

			//get shadow ray with max_t being half of the scene's bounding sphere radius
			Ray shadowRay = Ray(tempSurface.p, tempSurface.wi, Epsilon, scene.aabb.getBSphere().radius * 0.5);
			bool shadowIntersection = scene.bvh->intersect(shadowRay, tempSurface);

			//if there are no interesection with the shadow ray and the surface then we should calculate Li
			if (!shadowIntersection) {
				Li = v3f(albedo * INV_PI * (glm::max(0.f, cosTheta) / sampledPdf));
			}
		}

		return Li;
	}
};

TR_NAMESPACE_END