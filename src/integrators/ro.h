/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once
#include <random>
#include "bsdfs/phong.h"

TR_NAMESPACE_BEGIN

/**
 * Reflection occlusion integrator
 */
struct ROIntegrator : Integrator {

    float m_exponent;

    explicit ROIntegrator(const Scene& scene) : Integrator(scene) {
        m_exponent = scene.config.integratorSettings.ro.exponent;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }


    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);

        // TODO(A3): Implement this

		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		if (intersection) {
			//get sampledirection and pdf from warp functions
			v3f sampleDirection = Warp::squareToPhongLobe(sampler.next2D(), m_exponent);
			float sampledPdf = Warp::squareToPhongLobePdf(sampleDirection, m_exponent);
		
			//cosine lobe is originally aligned with the normal
			//rotate it to align with the sampling direction
			glm::quat rot(v3f(0.0f, 0.0f, 1.0f), reflect(tempSurface.wo));
			tempSurface.wi = glm::toMat4(rot) * v4f(sampleDirection, 1.f);
			tempSurface.wi = normalize(tempSurface.wi);

			//set wi in world coordinates and calculate cosTheta
			tempSurface.wi = tempSurface.frameNs.toWorld(tempSurface.wi);
			float cosTheta = glm::dot(tempSurface.frameNs.n, tempSurface.wi);

			//get shadow ray with max_t being the scene's bounding sphere radius
			Ray shadowRay = Ray(tempSurface.p, tempSurface.wi, Epsilon, scene.aabb.getBSphere().radius);
			bool shadowIntersection = scene.bvh->intersect(shadowRay, tempSurface);

			//if there are no interesection with the shadow ray and the surface then we should calculate Li
			if (!shadowIntersection) {
				float max = glm::max(0.f, pow(sampleDirection.z, m_exponent));
				Li = v3f((m_exponent + 2) * INV_TWOPI * glm::max(0.f, cosTheta) * max / sampledPdf);
			}
		}

        return Li;
    }
};

TR_NAMESPACE_END