/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/renderpass.h"
#include "tiny_obj_loader.h"

TR_NAMESPACE_BEGIN

/**
 * Analytic polygonal light source renderpass.
 * Based on the normal shading renderpass
 */
	struct PolygonalPass : RenderPass {
	GLuint diffuseShader{ 0 };
	GLuint emitterShader{ 0 };
	GLuint modelMatUniform{ 0 };
	GLuint viewMatUniform{ 0 };
	GLuint projectionMatUniform{ 0 };
	GLuint normalMatUniform{ 0 };
	bool firstPass = true;


	// Storing emitter vertex data for passing to shaders
	size_t nbVertices;
	size_t nbTriangles;
	std::vector<float> emitterVertexData; // Each point of each vertex is a float
	// Layout should be [ t1.v1.x t1.v1.y t1.v1.z t1.v2.x t1.v2.y t1.v2.z t1.v3.x t1.v3.y t1.v3.z t2.v1.x .... ]

	// Uniforms
	GLuint emitterVerticesUniform{ 0 };
	GLuint windowSizeUniform{ 0 }; // Needed for texture mapping of CV term
	GLuint lightIrradianceUniform{ 0 };
	v3f lightIrradiance{ 0 };
	GLuint nbTrianglesUniform{ 0 };

	// Allows for MAX_NUM_EMITTER_TRIANGLES triangles in polygonal light source
	int const MAX_NUM_EMITTER_TRIANGLES = 40;

	// Make sure MAX_NUM_EMITTER_TRIANGLES matches the value in polygonal.fs
	int const maxNbVertices = MAX_NUM_EMITTER_TRIANGLES * 3 * 3;

	// Emitter and sampler
	Emitter emitter;
	Sampler* sampler = nullptr;

	// Color the emitter properly
	GLuint lightIntensityUniform;

	// Storage for CV term
	std::vector<float> cvTermData;
	GLuint cvTermTexture;
	size_t numDataValues;

	// Used for computing CV term
	int samplesAccumulated = 0; // Used for moving average of CV term
	std::unique_ptr<PolygonalIntegrator> polygonalIntegrator;

	explicit PolygonalPass(const Scene& scene) : RenderPass(scene) {}

	/// Initialize renderpass
	bool init(const Config& config) override {
		RenderPass::init(config);

		// Create shader for diffuse surfaces
		GLuint vs = compileShader("polygonal.vs", GL_VERTEX_SHADER);
		GLuint fs = compileShader("polygonal.fs", GL_FRAGMENT_SHADER);
		diffuseShader = compileProgram(vs, fs);
		glDeleteShader(vs);
		glDeleteShader(fs);

		// Create shader for emitters
		vs = compileShader("polygonal.vs", GL_VERTEX_SHADER);
		fs = compileShader("emitter_polygonal.fs", GL_FRAGMENT_SHADER);
		emitterShader = compileProgram(vs, fs);
		glDeleteShader(vs);
		glDeleteShader(fs);

		shaders[DIFFUSE_SHADER_IDX] = diffuseShader;
		shaders[EMITTER_SHADER_IDX] = emitterShader;
		shaders[PHONG_SHADER_IDX] = -1; // Should not be any Phong surfaces, unless bonus

		// Initialize the polygonal integrator which will be used for computing CV term
		polygonalIntegrator = std::unique_ptr<PolygonalIntegrator>(new PolygonalIntegrator(scene));
		polygonalIntegrator->m_alpha = scene.config.integratorSettings.poly.alpha;
		polygonalIntegrator->m_visSamples = scene.config.integratorSettings.poly.visSamples;

		// Allocate storage for CV term
		numDataValues = scene.config.width * scene.config.height * 3;
		cvTermData.reserve(numDataValues);

		// Initialize texture and set it to zeros
		glGenTextures(1, &cvTermTexture);
		clearShadows();

		/**
		 * 1) Retrieve and assign `emitter`; create new `sampler` (seeded by your student ID)
		 * 2) Compute # of triangles on emitter and set `nbTriangles` and `nbVertices`
		 * 3) Compute and set `lightIrradiance`
		 * 4) Loop over each vertex and store xyz-components to `emitterVertexdata`
		 *      e.g. using: `emitterVertexData.push_back(vertex.x);`
		 */
		 // TODO(A4): Implement this

		//retrieve emitter and create new sampler
		emitter = scene.emitters.at(0);
		sampler = new Sampler(260678856);

		//retrieve the number of vertices - same as in part 1.2
		size_t shapeID = emitter.shapeID;
		auto shape = scene.worldData.shapes[shapeID];
		nbVertices = shape.mesh.indices.size();
		nbTriangles = nbVertices / 3;

		//find and set light irrandiance - same calculations as in offline part
		lightIrradiance = emitter.getPower() / emitter.area;
		lightIntensity = emitter.getRadiance();

		//loop over and store xyz components
		for (int i = 0; i < nbVertices; i++) {
			v3f vertex = scene.getObjectVertexPosition(shapeID, i);
			for (int j = 0; j < 3; j++) {
				emitterVertexData.push_back(vertex[j]); //vertex.x = vertex[0]
			}
		}

		// Create vertex buffers
		const auto& shapes = scene.worldData.shapes;
		objects.resize(shapes.size());
		for (size_t i = 0; i < objects.size(); i++) {
			buildVBO(i);
			buildVAO(i);
			assignShader(objects[i], shapes[i], scene.bsdfs);
		}
		return true;
	}

	/// Clean up after renderpass (deleter vertex buffers, etc.)
	void cleanUp() override {
		for (auto& object : objects) {
			glDeleteBuffers(1, &object.vbo);
			glDeleteVertexArrays(1, &object.vao);
		}
		RenderPass::cleanUp();
		delete sampler;
	}

	/// Event handler for shadows; clear on camera moved or draw on space bar pressed
	void handleEvents(SDL_Event& e) override {
		/**
		 * 1) Check if camera moved using `RenderPass::updateCamera(e)`: if so, _clear_ shadows
		 * 2) Check if space bar was pressed: if so, _add_ shadows
		 */
		 // TODO(A4): Implement this
		bool moved = RenderPass::updateCamera(e);

		//check if the camera has been moved
		bool cameraChange = RenderPass::updateCamera(e);
		//get the state of the all the keyboard keys - from SDL documentation
		const Uint8* state = SDL_GetKeyboardState(NULL);
		if (cameraChange) {
			clearShadows();
		}
		if (state[SDL_SCANCODE_SPACE]) {
			addShadows();
		}
	}

	/// Remove shadows; called when the camera moves
	void clearShadows() {
		cvTermData.clear();
		cvTermData.assign(numDataValues, 0.f);
		samplesAccumulated = 0;

		// Use this to bind the texture, and fill it with the updated data values
		glBindTexture(GL_TEXTURE_2D, cvTermTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene.config.width, scene.config.height, 0, GL_RGB, GL_FLOAT, cvTermData.data());
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	}

	/// Superpose ray-traced shadows on top of analytic diffuse shading
	void addShadows() {
		/**
		 * 1) Loop over pixels: height first _then_ width (order matters here!)
		 * 2) Trace ray through each pixel, just like A1
		 *    Note: Use _current_ camera settings (e.g. `camera.camera_look_at`)
		 * 3) Call `polygonalIntegrator->estimateVisibilityDiffHelper()` to estimate difference D = h - alpha*g / pdf
		 * 4) For each pixel, update entries of the vector `cvTermData` (e.g. `cvTermData[idx] = val`)
		 *    Note: Make sure to compute _moving average_ for each element using `samplesAccumulated`
		 *          Textures can only store positive values (and D < 0 always), so take -D
		 * 5) After loop, increment `samplesAccumulated`
		 * 6) Bind texture object and fill buffer with data from `cvTermData` (see `clearShadows()` for example)
		 */
		 // TODO(A4): Implement this
		float height = scene.config.height;
		float width = scene.config.width;
		v3f at = camera.camera_look_at;
		v3f up = camera.camera_up;
		v3f eye = camera.camera_position;

		//camera perspective
		float scaling = tan((deg2rad * scene.config.camera.fov) / 2.f);
		//camera-to-world transformation
		glm::mat4 inverseView = glm::inverse(glm::lookAt(eye, at, up));
		//aspect ratio
		float aspectRatio = (float)width / (float)height;

		//temp array of vector to store the different values of -diff
		std::vector<float> temp;

		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				//scaled width
				float px = ((x + 0.5) * 2 / width - 1) * scaling * aspectRatio;
				//scaled height
				float py = (((y + 0.5) * 2 / height) - 1) * scaling;
				v4f direction(px, py, -1.f, 0.f); //z=-1 because it's defined that way and w=1? because it's a vector representing a point
				//Convert to world view and normalize
				v4f directionWorld = inverseView * direction;
				directionWorld = normalize(directionWorld);
				//create a ray from the origin point of the camera and the direction calculated
				Ray ray = Ray(eye, directionWorld);

				//interest the surface with the ray
				SurfaceInteraction tempSurface;
				bool intersection = scene.bvh->intersect(ray, tempSurface);
				v3f diff(0.0f);
				//if there is an itersection we do the estimate
				if (intersection) {
					diff = polygonalIntegrator->estimateVisDiffRealTime(ray, *sampler, emitter);
				}
				
				//push the 3 components of diff to the temp array 
				//should be able to store -diff directly into cvTermData, however incorrect shadows were obtained
				for (int j = 0; j < 3; j++) {
					temp.push_back(-diff[j]);
				}
				
				//store the the coordinates of  -diff because textures can only take positive values and diff is negative  
				for (int i = 0; i < 3; i++) {
					//if this is the first iteration add the value of diff
					if (samplesAccumulated == 0) {
						cvTermData[y * width + x + i] = temp[y * width + x + i];
					}
					//if not the first iteration then do an average
					else {
						//formula found  on wikipedia - cumulative moving average
						cvTermData[y * width + x + i] = (temp[y * width + x + i] + samplesAccumulated * cvTermData[y * width + x + i]) / (samplesAccumulated + 1);
					}
				}
			}
		}
		//increment the samples
		samplesAccumulated++;

		//bind the texture, and fill it with the updated data values - similar to clearShadows()
		glBindTexture(GL_TEXTURE_2D, cvTermTexture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene.config.width, scene.config.height, 0, GL_RGB, GL_FLOAT, cvTermData.data());
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	}

	/// Main rendering routine
	void render() override {
		// Used for automated testing on TAs end
		// Make shadows appear without space bar press; used for grading
		if (scene.config.test && firstPass) {
			addShadows();
			firstPass = false;
		}

		// Standard real-time render
		glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
		glClearColor(0.f, 0.f, 0.f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_DEPTH_TEST);

		// Update camera
		glm::mat4 model, view, projection;
		camera.Update();
		camera.GetMatricies(projection, view, model);

		// Draw objects
		for (auto& object : objects) {
			// Define shader to use
			glUseProgram(object.shaderID);
			glActiveTexture(GL_TEXTURE0); // Required!
			glBindTexture(GL_TEXTURE_2D, cvTermTexture); // Need to bind texture after setting shader

			// Update camera
			camera.Update();
			camera.GetMatricies(projection, view, model);

			// Create uniforms
			modelMatUniform = GLuint(glGetUniformLocation(object.shaderID, "model"));
			viewMatUniform = GLuint(glGetUniformLocation(object.shaderID, "view"));
			projectionMatUniform = GLuint(glGetUniformLocation(object.shaderID, "projection"));
			normalMatUniform = GLuint(glGetUniformLocation(object.shaderID, "normalMat"));

			// Pass common Uniforms
			glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
			glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
			glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));
			glUniformMatrix4fv(normalMatUniform, 1, GL_FALSE, &(normalMat[0][0]));

			/**
			 * 1) Check if object is emitter or diffuse (no Phong, unless bonus)
			 * 2) Bind correct uniforms for each shader
			*/
			// TODO(A4): Implement this
			//check which shader the object has initialized to determine which uniforms needs to be passed
			//the uniforms passed are defined in the polygonal and emitter_polygonal
			if (object.shaderIdx == DIFFUSE_SHADER_IDX) {
				GLuint albedoUniform = GLuint(glGetUniformLocation(object.shaderID, "albedo"));
				glUniform3f(albedoUniform, object.albedo.x, object.albedo.y, object.albedo.z);
				emitterVerticesUniform = GLuint(glGetUniformLocation(object.shaderID, "emitterVertices"));
				glUniform1fv(emitterVerticesUniform, emitterVertexData.size(), emitterVertexData.data());
				nbTrianglesUniform = GLuint(glGetUniformLocation(object.shaderID, "nbTriangles"));
				glUniform1i(nbTrianglesUniform, nbTriangles);
				lightIrradianceUniform = GLuint(glGetUniformLocation(object.shaderID, "lightIrradiance"));
				glUniform3f(lightIrradianceUniform, lightIrradiance.x, lightIrradiance.y, lightIrradiance.z);
				windowSizeUniform = GLuint(glGetUniformLocation(object.shaderID, "windowSize"));
				glUniform2f(windowSizeUniform, width, height);
			}
			else if (object.shaderIdx == EMITTER_SHADER_IDX) {
				lightIntensityUniform = GLuint(glGetUniformLocation(object.shaderID, "lightIntensity"));
				glUniform3f(lightIntensityUniform, lightIntensity.x, lightIntensity.y, lightIntensity.z);
			}

			// Bind VAO, draw and bind to zero
			glBindVertexArray(object.vao);
			glDrawArrays(GL_TRIANGLES, 0, object.nVerts);
			glBindVertexArray(0);
		}
		RenderPass::render();
	}
};

TR_NAMESPACE_END