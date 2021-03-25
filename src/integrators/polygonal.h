/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <tiny_obj_loader.h>
#define RAY_EPS_CV 1e-5 // Use when setting min and max dist for ray in control variates code
TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator for polygonal light sources
 * Follows Arvo '94.
 */
	struct PolygonalIntegrator : Integrator {

	float m_alpha;             // Control variates "strength"
	size_t m_visSamples;       // # of samples to estimate h - alpha*g
	bool m_traceShadows;       // Trace shadows or not
	EPolygonalMethod m_method; // Method to use (Arvo, or control variates)

	std::vector<std::vector<v3f>> m_triangles; // Data structure to store triangles

	explicit PolygonalIntegrator(const Scene& scene) : Integrator(scene) {
		m_alpha = scene.config.integratorSettings.poly.alpha;
		m_visSamples = scene.config.integratorSettings.poly.visSamples;
		m_traceShadows = scene.config.integratorSettings.poly.traceShadows;
		m_method = scene.config.integratorSettings.poly.method;

		/**
		 * 1) Get # of triangles on emitter
		 * 2) Store vertices in m_triangles
		 */
		 // TODO(A4): Implement this

		size_t shapeID = scene.getFirstLight();
		auto shape = scene.worldData.shapes[shapeID];
		size_t numVertices = shape.mesh.indices.size(); //total number of vertices, divide by 3 to get number of triangles

		//iterate over all the vertices and assign 3 for one triangle, add to list of triangles
		for (int i = 0; i < numVertices; i += 3) {
			//printf("%i", i);
			std::vector<v3f> tri(3);    // Create empty triangle (3 vertices / triangle)
			tri[0] = scene.getObjectVertexPosition(shapeID, i);       // Add vertices to triangle
			tri[1] = scene.getObjectVertexPosition(shapeID, i + 1);
			tri[2] = scene.getObjectVertexPosition(shapeID, i + 2);
			m_triangles.push_back(tri); // Append triangle vector to vector
		}
	}

	/// Reflect
	inline v3f reflect(const v3f& d) const {
		return v3f(-d.x, -d.y, d.z);
	}

	/**
	 * === PHONG BONUS ONLY ===
	 * Compute the following integral:
	 *    T(a, b, n, x) = \int_0^x [a \cos(\theta) + b \sin(\theta)]ˆn d\theta
	 * Uses a recurrent relation (see Snyder's note, 1996)
	 *
	 * Series function:
	 *    T_sum(a, b, n, x) = \sum_{i=0}ˆ{(n-1)/2} T(a, b, 2i+1, x)
	 * assuming n is _odd_
	 */
	float cosineSinePowerIntegralSum(float a, float b, int exp, float theta) const {
		if (exp % 2 == 0) exp += 1; // Make exponent odd
		float Tsum = 0.f;

		// Implementing this function may be useful if you attempt the bonus

		// TODO(A4): Implement this

		return Tsum;
	}

	/**
	 * Compute edge (v1--v2) contribution
	 * The exp term is only needed if you attempt the bonus, otherwise, you can ignore it
	 */
	float getEdgeContrib(const v3f& v1, const v3f& v2, const SurfaceInteraction& i, int exp = 0) const {
		float contrib = 0.f;

		// TODO(A4): Implement this
		v3f x = i.p;
		v3f firstTerm = (v1 - x) / glm::l2Norm(v1, x);
		v3f secondTerm = (v2 - x) / glm::l2Norm(v2, x);
		float subtendedAngle = glm::acos(glm::dot(firstTerm, secondTerm));
		v3f crossProduct = glm::cross(secondTerm, (v1 - x));
		crossProduct = crossProduct / glm::l2Norm(crossProduct);
		contrib = subtendedAngle * glm::dot(crossProduct, i.frameNs.n);
		return contrib;
	}


	/// Direct illumination using Arvo '94 analytic solution for polygonal lights
	v3f renderAnalytic(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);

		// TODO(A4): Implement this

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
			float temp = 0;
			for (int i = 0; i < m_triangles.size(); i++) {
				temp += getEdgeContrib(m_triangles.at(i).at(0), m_triangles.at(i).at(1), tempSurface);
				temp += getEdgeContrib(m_triangles.at(i).at(1), m_triangles.at(i).at(2), tempSurface);
				temp += getEdgeContrib(m_triangles.at(i).at(2), m_triangles.at(i).at(0), tempSurface);
			}

			//PDF to be returned
			float emitterPdf;

			//Gets emitter and data
			size_t id = selectEmitter(sampler.next(), emitterPdf);
			const Emitter& emitter = getEmitterByID(id);

			v3f e = emitter.getPower() / emitter.area;
			v3f LrTemp = e * INV_TWOPI * temp;
			SurfaceInteraction surface = tempSurface;
			surface.wi = v3f(0, 0, 1); 
			v3f albedo = getBSDF(surface)->eval(surface);  //returns albedo over pi 

			Lr = albedo * LrTemp;

		}
		return Lr;
	}

	/**
	 * Stand-alone estimator for h - alpha*g (with primary ray)
	 * Trace a primary ray, check for emitter hit, and then call `estimateVisDiff()`
	 * Used by polygonal render pass
	 */
	v3f estimateVisDiffRealTime(const Ray& ray, Sampler& sampler, const Emitter& em) {
		v3f D(0.f);

		SurfaceInteraction hit;
		if (!scene.bvh->intersect(ray, hit)) return D;

		const BSDF* bsdf = getBSDF(hit);
		if (bsdf->isEmissive()) return D;

		hit.wi = v3f(0, 0, 1); // Trick to get 1/pi * albedo without cosine term
		D = estimateVisDiff(sampler, hit, em);

		return D;
	}

	/// Stand-alone estimator for h - alpha*g (without primary ray)
	/// Use RAY_EPS_CV when setting min and max dist for shadow ray
	v3f estimateVisDiff(Sampler& sampler, SurfaceInteraction& i, const Emitter& em) const {
		v3f sum(0.f);

		// TODO(A4): Implement this
		//need to calculate h - alpha * g over all the samples
		v3f h(0.f);
		v3f g(0.f);
		for (int j = 0; j < m_visSamples; j++) {

			//get emitter pdf
			float emitterPdf = getEmitterPdf(em);
			
			v3f pe;    // Point on emitter
			v3f ne;    // Surface normal at point
			float pdf; // PDF of choosing point
			sampleEmitterPosition(sampler, em, ne, pe, pdf); // Sample mesh uniformly

			//calculate wi in world coordinates
			v3f wiW = glm::normalize(pe - i.p);
			//wi must be in local coordinates
			i.wi = i.frameNs.toLocal(wiW);

			v3f lightIntensity = scene.getFirstLightIntensity();

			Ray shadowRay = Ray(i.p, wiW, RAY_EPS_CV, glm::distance(pe, i.p) + RAY_EPS_CV); 

			SurfaceInteraction shadowSurface;
			bool shadowIntersection = scene.bvh->intersect(shadowRay, shadowSurface);

			//calculate the visibility term 
			float visibilityTerm = glm::dot(pe - i.p, ne);
			float visibility = 0;
			if (visibilityTerm <= 0) {
				visibility = 1;
			}

			//calculate the jacobian determinant of the transformation
			float geometryTerm = glm::abs(glm::dot(ne, wiW) / glm::pow((glm::distance(pe, i.p)), 2));

			//Checks if there is no object intersection
			if (shadowIntersection) {
				//get the emission of the shadow ray 
				v3f shadowEmission = getEmission(shadowSurface);

				//Checks if hits a surface emitter
				if (shadowEmission != v3f(0)) {
					//Computes the value of the shadow
					h += visibility * shadowEmission * geometryTerm * (getBSDF(i)->eval(i)) / (pdf * emitterPdf);
				}
			}
			//Computes the value of the MC Estimator
			g += visibility * lightIntensity * geometryTerm * (getBSDF(i)->eval(i)) / (pdf * emitterPdf);
		}
		sum = h - (m_alpha * g);
		return sum;
	}

	/// Control variates using Arvo '94 for direct illumination; ray trace shadows

	v3f renderControlVariates(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);

		// TODO(A4): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		v3f LrG(0.f);

		//Checks if ray intersects the scene (to display emitter in scene)
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if it hits an emitter
			if (emission != v3f(0)) {
				//return the emission of the emitter
				Lr = emission;
				return Lr;
			}
			float temp = 0;
			for (int i = 0; i < m_triangles.size(); i++) {
				temp += getEdgeContrib(m_triangles.at(i).at(0), m_triangles.at(i).at(1), tempSurface);
				temp += getEdgeContrib(m_triangles.at(i).at(1), m_triangles.at(i).at(2), tempSurface);
				temp += getEdgeContrib(m_triangles.at(i).at(2), m_triangles.at(i).at(0), tempSurface);
			}

			//PDF to be returned
			float emitterPdf;

			//Gets emitter and data
			size_t id = selectEmitter(sampler.next(), emitterPdf);
			const Emitter& emitter = getEmitterByID(id);

			v3f e = emitter.getPower() / emitter.area;
			v3f LrTemp = e * INV_TWOPI * temp;
			SurfaceInteraction surface = tempSurface;
			surface.wi = v3f(0, 0, 1); 
			v3f albedo = getBSDF(surface)->eval(surface);  //returns albedo over pi 
			v3f g = albedo * LrTemp;

			v3f sum = estimateVisDiff(sampler, tempSurface, emitter) / m_visSamples;

			Lr = m_alpha * g + sum;
			//from the discussion board onmycourses
			Lr = clampBelow(Lr, 0);
		}
		return Lr;
	}

	/// Direct illumination using surface area sampling
	v3f renderArea(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);

		// TODO(A4): Implement this
		SurfaceInteraction tempSurface;
		bool intersection = scene.bvh->intersect(ray, tempSurface);

		//Checks if ray intersects the scene (to display emitter in scene)
		if (intersection) {
			v3f emission = getEmission(tempSurface);

			//Checks if hits directly an emitter
			if (emission != v3f(0)) {
				Lr = emission;
				return Lr;
			}

			//iterate over the number of emitter samples
			for (int i = 0; i < m_triangles.size(); i++) {
				//PDF to be returned
				float emitterPdf;

				//Gets emitter and data
				size_t id = selectEmitter(sampler.next(), emitterPdf);
				const Emitter& emitter = getEmitterByID(id);
				v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
				float emitterRadius = scene.getShapeRadius(emitter.shapeID);

				v3f pe;    // Point on emitter
				v3f ne;    // Surface normal at point
				float pdf; // PDF of choosing point
				sampleEmitterPosition(sampler, emitter, ne, pe, pdf); // Sample mesh uniformly

				//calculate wi in world coordinates
				v3f wiW = glm::normalize(pe - tempSurface.p);
				//wi must be in local coordinates
				tempSurface.wi = tempSurface.frameNs.toLocal(wiW);

				//get intensity of scene
				v3f intensity = scene.getFirstLightIntensity();

				//trace shadow rays or not
				if (m_traceShadows) {
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
							float geometryTerm = glm::abs(glm::dot(ne, wiW) / glm::pow((glm::distance(pe, tempSurface.p)), 2));

							//calculate the visibility term 
							float visibilityTerm = glm::dot(pe - tempSurface.p, ne);
							float visibility = 0;
							if (visibilityTerm <= 0) {
								visibility = 1;
							}
							//Computes the value for the monte carlo estimator 
							Lr += visibility * shadowEmission * geometryTerm * (getBSDF(tempSurface)->eval(tempSurface)) / (pdf * emitterPdf);
						}
					}
				}
				else {
					//calculate the jacobian determinant of the transformation
					float geometryTerm = glm::abs(glm::dot(ne, wiW) / glm::pow((glm::distance(pe, tempSurface.p)), 2));

					//calculate the visibility term 
					float visibilityTerm = glm::dot(pe - tempSurface.p, ne);
					float visibility = 0;
					if (visibilityTerm <= 0) {
						visibility = 1;
					}
					//Computes the value for the monte carlo estimator 
					Lr += visibility * intensity * geometryTerm * (getBSDF(tempSurface)->eval(tempSurface)) / (pdf * emitterPdf);
				}
			}
			Lr = Lr / m_triangles.size();
		}
		return Lr;
	}

	/// Branch to corresponding method
	v3f render(const Ray& ray, Sampler& sampler) const override {
		switch (m_method) {
		case EPolygonalMethod::ESurfaceArea:
			return PolygonalIntegrator::renderArea(ray, sampler);
			break;
		case EPolygonalMethod::EControlVariates:
			return PolygonalIntegrator::renderControlVariates(ray, sampler);
			break;
		default:
			return PolygonalIntegrator::renderAnalytic(ray, sampler);
			break;
		}
	}

};

TR_NAMESPACE_END