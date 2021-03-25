/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
struct PathTracerIntegrator : Integrator {
    explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
        m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
        m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
        m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
        m_rrProb = scene.config.integratorSettings.pt.rrProb;
    }


    v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
        v3f Li(0.f);

        // TODO(A5): Implement this
		//temp value for implicit path tracing
		v3f Limplicit(1.f);
		//iterate until we hit max depth i.e. max number of bounces
		for (int i = 0; i <= m_maxDepth ; i++) {
			//check if we hit an emitter and if we're on the correct side of the object in the scene - from tutorial
			v3f emission = getEmission(hit);
			if (emission != v3f(0.f) && dot(v3f(0.f, 0.f, 1.f), hit.wo) > 0.f) {
				//return the accumlated BRDF and take into account the emission
				Limplicit *= emission;
				return Li = Limplicit;
			}
			//get the bsdf from the sampler in phong
			float pdf;
			v3f bsdf = getBSDF(hit)->sample(hit, sampler, &pdf);
			//when adding cosine forshortening factor I obtain weird results for my image 
			//float cosTheta = glm::dot(hit.frameNs.n, hit.wi);
			//accumulate
			Limplicit *= bsdf  / pdf;
			
			//create a new ray as the direction was updated in the sample function
			Ray ray2 = Ray(hit.p, glm::normalize(hit.frameNs.toWorld(hit.wi)));
			//intersect ray and surface
			bool intersection = scene.bvh->intersect(ray2, hit);
			//if there's no intersection we break and return black
			if (!intersection) {
				break;
			}
		}
		return Li;
    }

	//helper method that performs directIllumnation same as renderArea in direct.h minus the check for emission and the loop as 
	//this is being called at every recurrence in indirectillumination
	v3f directIllumination(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {		
		//temp value for direct illumination
		v3f Ldirect(0.f);

		//PDF to be returned
		float emitterPdf;

		//Gets emitter and data
		size_t id = selectEmitter(sampler.next(), emitterPdf);
		const Emitter& emitter = getEmitterByID(id);
		v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
		float emitterRadius = scene.getShapeRadius(emitter.shapeID);

		//Samples the light
		v3f pos;
		v3f ne;

		//get the sample direction and pdf 
		v3f wiW;
		float pdf;

		//sample position - "returns" the normal, position and pdf
		sampleEmitterPosition(sampler, emitter, ne, pos, pdf);

		//need to get wiw - look at sample function from direct.h 
		//Position in world 
		wiW = glm::normalize(pos - hit.p);
		//wi must be in local coordinates
		hit.wi = hit.frameNs.toLocal(wiW);

		//Creates a ray with the intersection point and the sampled direction
		Ray shadowRay = Ray(hit.p, wiW);
		SurfaceInteraction shadowSurface;
		bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);
		
		//checks if there is no object intersection
		if (shadowIntersection) {
			//get the emission of the shadow ray 
			v3f shadowEmission = getEmission(shadowSurface);

			//checks if hits a surface emitter
			if (shadowEmission != v3f(0)) {
				//calculate the jacobian determinant of the transformation
				float geometryTerm = glm::abs(glm::dot(ne, wiW)) / powf((glm::distance(pos, hit.p)), 2);
				//calculate the visibility term 
				float visibilityTerm = glm::dot(emitterCenter - hit.p, ne);
				float visibility = 1;
				if (visibilityTerm <= 0) {
					visibility = 0;
				}
				//Computes the value for the monte carlo estimator 
				Ldirect = shadowEmission * geometryTerm * (getBSDF(hit)->eval(hit)) / (pdf * emitterPdf);
			}
		}
		return Ldirect;
	}

	//helper method for indirect illumnation - done recursively 
	//similar to what is done in the implicit path tracing
	v3f indirectIllumination(int iterations, const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		//temp value for indirect illumination
		v3f Lindirect(0.f);
		v3f temp(0.f);

		//check if we hit an emitter and if we're on the correct side of the object in the scene 
		v3f emission = getEmission(hit);
		if (emission != v3f(0.f) && dot(v3f(0.f, 0.f, 1.f), hit.wo) > 0.f) {
			return emission;
		}

		//Continue recursion if we haven't hit the max number of bounces
		if (iterations < m_maxDepth) {
			//get the bsdf from the sampler in phong
			float pdf;
			v3f bsdf = getBSDF(hit)->sample(hit, sampler, &pdf);
			//when adding cosine forshortening factor I obtain weird results for my image 
			//float cosTheta = glm::dot(hit.frameNs.n, hit.wi);
			temp = bsdf / pdf;
			//create a new ray as the direction was updated in the sample function
			Ray ray2 = Ray(hit.p, glm::normalize(hit.frameNs.toWorld(hit.wi)));
			bool intersection = scene.bvh->intersect(ray2, hit);
			//if there's no intersection we exit the recursion and return black
			if (!intersection) {
				temp = v3f(0.f);
				return temp;
			}
		}
		//when we hit the max number of bounces stop recursion
		else {
			temp = v3f(0.f);
			return temp;
		}

		// Russian Roulette - from PBR textbook
		if (iterations > m_rrDepth) {
			float p = m_rrProb;
			if (sampler.next() > p) {
				temp = v3f(0.f);
				return temp;
			}
			temp /= 1 - p;
		}

		//do recursion 
		Lindirect = temp * (indirectIllumination(iterations++, ray, sampler, hit) + directIllumination(ray, sampler, hit));
		return Lindirect;
	}

	v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
		v3f Li(0.f);

		// TODO(A5): Implement this

		//DIRECT ILLUMINATION
		v3f Ldirect = directIllumination(ray, sampler, hit);

		//INDIRECT ILLUMINATION
		v3f Lindirect = indirectIllumination(1, ray, sampler, hit);

		Li = Lindirect + Ldirect;
		return Li;

	}

    v3f render(const Ray& ray, Sampler& sampler) const override {
        Ray r = ray;
        SurfaceInteraction hit;

        if (scene.bvh->intersect(r, hit)) {
            if (m_isExplicit)
                return this->renderExplicit(ray, sampler, hit);
            else
                return this->renderImplicit(ray, sampler, hit);
        }
        return v3f(0.0);
    }

    int m_maxDepth;     // Maximum number of bounces
    int m_rrDepth;      // When to start Russian roulette
    float m_rrProb;     // Russian roulette probability
    bool m_isExplicit;  // Implicit or explicit
};

TR_NAMESPACE_END
