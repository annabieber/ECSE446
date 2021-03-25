/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Simple direct illumination integrator.
 */
struct SimpleIntegrator : Integrator {
    explicit SimpleIntegrator(const Scene& scene) : Integrator(scene) { }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        v3f Li(0.f);

        // TODO(A2): Implement this
		SurfaceInteraction temp;
		bool intersection = scene.bvh->intersect(ray, temp);

		if (intersection) { 
			//retrieve the light position p and its intensity I
			v3f position = scene.getFirstLightPosition();
			v3f intensity = scene.getFirstLightIntensity();

			//transform the world - space incoming direction i.wi to the surface's local coordinates
			temp.wi = temp.frameNs.toLocal(position - temp.p);

			//retrieve the intersected surface's material using getBSDF(*)
			v3f surfaceMaterial = getBSDF(temp)->eval(temp);

			//create shadow ray from intersection point on the surface 
			Ray shadowRay = Ray(temp.p, glm::normalize(position - temp.p), Epsilon, glm::length(temp.p - position));
			bool shadowIntersection = scene.bvh->intersect(shadowRay, temp);

			//if there are no interesection with the shadow ray and the surface then we should calculate Li
			if (!shadowIntersection) {
				Li = intensity / glm::pow(glm::length(temp.p - position), 2) * surfaceMaterial;
			}
		}
        return Li;
    }
};

TR_NAMESPACE_END