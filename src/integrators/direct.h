/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        // TODO(A3): Implement this
		//sample direction with cosine to hemisphere warp function
		v3f sampleDirection = Warp::squareToCosineHemisphere(sample);
		//create frame and set wi in world coordinates
		Frame f = Frame(n);
		wiW = f.toWorld(sampleDirection);
		//calculate the pdf in the sampled direction
		pdf = Warp::squareToCosineHemispherePdf(sampleDirection);
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO(A3): Implement this
		v3f sampleDirection = Warp::squareToUniformSphere(sample);

		//Position in world
		pos = emitterCenter + sampleDirection * emitterRadius;
		wiW = glm::normalize(pos - pShading);

		//Normal of the emitter must be in local
		ne = normalize(emitterRadius * sampleDirection);

		//get pdf 
		pdf = Warp::squareToUniformSpherePdf();
		pdf = pdf / pow(emitterRadius, 2); //pdf is going to return 1/4PISquared so we want to divide by the emitter radius squared
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO(A3): Implement this
		//find max cos theta by using the emitter and shading
		float thetaMax = asinf(emitterRadius / glm::distance(pShading, emitterCenter));
		float cosThetaMax = cosf(thetaMax);

		//get sample direction and pdf
		v3f sampleDirection = Warp::squareToUniformCone(sample, cosThetaMax);
		pdf = Warp::squareToUniformConePdf(cosThetaMax);

		//wiW is aligned to the z axis
		//Rotates the cone to align it with the emitter
		glm::quat rot(v3f(0.0f, 0.0f, 1.0f), emitterCenter - pShading);
		wiW = glm::toMat4(rot) * v4f(sampleDirection, 1.f);
		wiW = normalize(wiW);
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene 
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if hits directly an emitter
			if (emission != v3f(0)) {
				Lr = emission;
				return Lr;
			}

			//iterate over the number of emitter samples
			for (int i = 0; i < m_emitterSamples; i++) {
				//PDF to be returned
				float emitterPdf;

				//Gets emitter and data
				size_t id = selectEmitter(sampler.next(), emitterPdf);
				const Emitter& emitter = getEmitterByID(id);
				v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
				float emitterRadius = scene.getShapeRadius(emitter.shapeID);

				//Samples the light
				v3f pShading = tempSurface.p;
				v3f pos;
				v3f ne;

				//get the sample direction and pdf 
				v3f wiW;
				float pdf;
				sampleSphereByArea(sampler.next2D(), pShading, emitterCenter, emitterRadius, pos, ne, wiW, pdf);
				//wi must be in local coordinates
				tempSurface.wi = tempSurface.frameNs.toLocal(wiW);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, wiW);
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//Checks if there is no object intersection
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);

					//Checks if hits a surface emitter
					if (shadowEmission != v3f(0)) {
						//calculate the jacobian determinant of the transformation
						float geometryTerm = glm::abs(glm::dot(ne, wiW)) / powf((glm::distance(pos, tempSurface.p)), 2);

						//calculate the visibility term 
						float visibilityTerm = glm::dot(emitterCenter - tempSurface.p, ne);
						float visibility = 1;
						if (visibilityTerm <= 0) {
							visibility = 0;
						}
						//Computes the value for the monte carlo estimator 
						Lr += visibility * shadowEmission * geometryTerm * (getBSDF(tempSurface)->eval(tempSurface)) / (pdf * emitterPdf);
					}
				}
			}
			Lr = Lr / m_emitterSamples;
		}
        return Lr;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if it hits an emitter
			if (emission != v3f(0)) {
				//return the emission of the emitter
				Lr = emission;
				return Lr;
			}

			//iterate over the number of bsdf samples
			for (int i = 0; i < m_bsdfSamples; i++) {

				//get the sample direction and pdf 
				float pdf;
				v3f wiW;
				sampleSphereByCosineHemisphere(sampler.next2D(), tempSurface.frameNs.n, p3f(0.f), v3f(0.f), 0.f, wiW, pdf);

				//set wi in local coordinates
				tempSurface.wi = tempSurface.frameNs.toLocal(wiW);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, wiW);
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//check if it hits a surface emitter
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);
					if (shadowEmission != v3f(0)) {
						//Computes the value for the monte carlo estimator
						Lr += shadowEmission * (getBSDF(tempSurface)->eval(tempSurface)) / pdf;
					}
				}
			}
			Lr = Lr / m_bsdfSamples;
		}
        return Lr;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if it hits an emitter
			if (emission != v3f(0)) {
				//return the emission of the emitter
				Lr = emission;
				return Lr;
			}

			//iterate over the number of bsdf samples
			for (int i = 0; i < m_bsdfSamples; i++) {

				//get the bsdf from the sampler in phong
				float pdf;
				v3f bsdf = getBSDF(tempSurface)->sample(tempSurface, sampler, &pdf);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, tempSurface.frameNs.toWorld(tempSurface.wi));
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//Checks if there is no object intersection
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);
					
					//check if it hits a surface emitter
					if (shadowEmission != v3f(0)) {
						//calculate the monte carlo estimator
						Lr += shadowEmission * (bsdf) / pdf;
					}
				}
			}
			Lr = Lr / m_bsdfSamples;
		}
        return Lr;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

        // TODO(A3): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if it hits an emitter
			if (emission != v3f(0)) {
				//return the emission of the emitter
				Lr = emission;
				return Lr;
			}

			//iterate over the number of emitter samples
			for (int i = 0; i < m_emitterSamples; i++) {

				//PDF to be returned
				float emitterPdf;

				//Gets emitter and data
				size_t id = selectEmitter(sampler.next(), emitterPdf);
				const Emitter& emitter = getEmitterByID(id);
				v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
				float emitterRadius = scene.getShapeRadius(emitter.shapeID);

				//Samples the light
				v3f pShading = tempSurface.p;

				//get the sample direction and pdf 
				v3f wiW;
				float pdf;
				sampleSphereBySolidAngle(sampler.next2D(), pShading, emitterCenter, emitterRadius, wiW, pdf);
				//wi must be in local coordinates
				tempSurface.wi = tempSurface.frameNs.toLocal(wiW);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, wiW);
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//Checks if there is no object intersection
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);

					//Checks if hits a surface emitter
					if (shadowEmission != v3f(0)) {
						//Computes the value for the monte carlo estimator 
						Lr += shadowEmission * (getBSDF(tempSurface)->eval(tempSurface)) / (pdf * emitterPdf);
					}
				}
			}
			Lr = Lr / m_emitterSamples;
		}
        return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {

        v3f Lr(0.f);

        // TODO(A4): Implement this

		//temporary Lr values for BSDF and emitter respectively
		v3f LrBSDF;
		v3f LrSolidAngle;

		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene (to display emitter in scene)
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if it hits an emitter
			if (emission != v3f(0)) {
				//return the emission of the emitter
				Lr = emission;
				return Lr;
			}

			//Taken from the BSDf renderer
			//iterate over the number of bsdf samples
			for (int i = 0; i < m_bsdfSamples; i++) {

				//get the bsdf from the sampler in phong
				float pdf;
				v3f bsdf = getBSDF(tempSurface)->sample(tempSurface, sampler, &pdf);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, tempSurface.frameNs.toWorld(tempSurface.wi));
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//Checks if there is no object intersection
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);

					//check if it hits a surface emitter
					if (shadowEmission != v3f(0)) {
						//Gets emitter and data
						size_t id = getEmitterIDByShapeID(shadowSurface.shapeID);
						const Emitter& emitter = getEmitterByID(id);
						float emitterPdf = getEmitterPdf(emitter);
						v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
						float emitterRadius = scene.getShapeRadius(emitter.shapeID);

						//find max cos theta by using the emitter and shading
						float thetaMax = asinf(emitterRadius / glm::distance(tempSurface.p, emitterCenter));
						float cosThetaMax = cosf(thetaMax);

						//sample the pdf with cos theta 
						float samplePdf = Warp::squareToUniformConePdf(cosThetaMax) ; 
						
						//calculate the weight for the bsdf
						float weightBSDF = balanceHeuristic(m_bsdfSamples, pdf, m_emitterSamples, samplePdf * emitterPdf); //call the weighting function to get the weight

						LrBSDF += weightBSDF * shadowEmission * (bsdf) / pdf;
					}
				}
			}

			//Taken from the solid angle rendering - TODO take into account the weighted function from the handout
			//iterate over the number of emitter samples
			for (int j = 0; j < m_emitterSamples; j++) {

				//PDF to be returned
				float emitterPdf;

				//Gets emitter and data
				size_t id = selectEmitter(sampler.next(), emitterPdf);
				const Emitter& emitter = getEmitterByID(id);
				v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
				float emitterRadius = scene.getShapeRadius(emitter.shapeID);

				//Samples the light
				v3f pShading = tempSurface.p;

				//get the sample direction and pdf 
				v3f wiW;
				float pdf;
				sampleSphereBySolidAngle(sampler.next2D(), pShading, emitterCenter, emitterRadius, wiW, pdf);
				//wi must be in local coordinates
				tempSurface.wi = tempSurface.frameNs.toLocal(wiW);

				//Creates a ray with the intersection point and the sampled direction
				Ray shadowRay = Ray(tempSurface.p, wiW);
				SurfaceInteraction shadowSurface;
				bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

				//Checks if there is no object intersection
				if (shadowIntersection) {
					//get the emission of the shadow ray 
					v3f shadowEmission = getEmission(shadowSurface);

					//Checks if hits a surface emitter
					if (shadowEmission != v3f(0)) {

						float weight = balanceHeuristic(m_emitterSamples, pdf * emitterPdf, m_bsdfSamples, getBSDF(tempSurface)->pdf(tempSurface));
						//Computes the value for the monte carlo estimator 
						LrSolidAngle += weight * shadowEmission * (getBSDF(tempSurface)->eval(tempSurface)) / (pdf * emitterPdf);
					}
				}
			}

			if (m_emitterSamples == 0) {
				Lr = LrBSDF / m_bsdfSamples;
				return Lr;
			}

			if (m_bsdfSamples == 0) {
				Lr = LrSolidAngle / m_emitterSamples;
				return Lr; 
			}
		}
		//return Lr where each Lr is divided by it's respective sampler 
		Lr = LrBSDF / m_bsdfSamples + LrSolidAngle / m_emitterSamples;
		
        return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == ESamplingStrategy::EMIS)
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::EArea)
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ESolidAngle)
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == ESamplingStrategy::ECosineHemisphere)
            return this->renderCosineHemisphere(ray, sampler);
        else
            return this->renderBSDF(ray, sampler);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    ESamplingStrategy m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END