/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/core.h>
#include "core/renderpass.h"
#include "tiny_obj_loader.h"
#include "integrators/path.h"

TR_NAMESPACE_BEGIN

/**
 * Global Illumination baking renderpass.
 */
struct GIPass : RenderPass {
    GLuint shader{0};

    GLuint modelMatUniform{0};
    GLuint viewMatUniform{0};
    GLuint projectionMatUniform{0};

    int m_samplePerVertex;

    std::unique_ptr<PathTracerIntegrator> m_ptIntegrator;

    explicit GIPass(const Scene& scene) : RenderPass(scene) {
        m_ptIntegrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        m_ptIntegrator->m_maxDepth = scene.config.integratorSettings.gi.maxDepth;
        m_ptIntegrator->m_rrProb = scene.config.integratorSettings.gi.rrProb;
        m_ptIntegrator->m_rrDepth = scene.config.integratorSettings.gi.rrDepth;
        m_samplePerVertex = scene.config.integratorSettings.gi.samplesByVertex;
    }

    virtual void buildVBO(size_t objectIdx) override {
        GLObject& obj = objects[objectIdx];

        // TODO(A5): Implement this
		
		//same as renderpass.cpp
		obj.nVerts = scene.getObjectNbVertices(objectIdx);
		obj.vertices.resize(obj.nVerts * N_ATTR_PER_VERT);
		//sampler for function
		Sampler sampler = Sampler(260678856);

		int t = 0; 
		//need a temp variable to iterate over all the vertices other than i
		for (int i = 0; i < obj.nVerts; i ++) {
			
			//1. set an arbitrary ray direction (ωo) and shift the shading point position by ϵ along the normal to avoid self-intersections during integration.
			v3f wo = v3f(1.f, 1.f, 1.f);
			Ray ray = Ray(v3f(0.f), wo); 
			//shift shading position by E along normal
			v3f normal = scene.getObjectVertexNormal(objectIdx, i);
			v3f e = Epsilon * normal;
			v3f position = scene.getObjectVertexPosition(objectIdx, i);


			//2. Create a SurfaceInteraction and manually populate its attributes. 
			SurfaceInteraction tempSurface;
			//tempSurface.t = ;
			//tempSurface.u = ;
			//tempSurface.v = ;
			tempSurface.p = position + e;
			tempSurface.wo = wo;
			tempSurface.shapeID = objectIdx; 
			tempSurface.primID = scene.getPrimitiveID(i);
			tempSurface.frameNg = Frame(normal); 
			tempSurface.frameNs = Frame(normal);
			tempSurface.matID = scene.getMaterialID(objectIdx, tempSurface.primID);

			//3. Call your PathTracerIntegrator->renderExplicit() with the incoming “ray” and the SurfaceInteraction references.
			v3f color(0.f);
			//loop over samples per vertex
			for (int j = 0; j < m_samplePerVertex; j++) {
				color += m_ptIntegrator->renderExplicit(ray, sampler, tempSurface);
			}
			//get average
			color /= m_samplePerVertex;

			//Vertices need to contain the raw data of the vertex attributes stored one after the other
			//i.e. vertices[0 to 5] = x,y,z,R,G,B (vertex 0 data)
			obj.vertices[t + 0] = position.x; 
			obj.vertices[t + 1] = position.y;
			obj.vertices[t + 2] = position.z;
			obj.vertices[t + 3] = color.r;
			obj.vertices[t + 4] = color.g;
			obj.vertices[t + 5] = color.b;

			//N_ATTR_PER_VERT = 6
			t += N_ATTR_PER_VERT;
		}

        // VBO
        glGenVertexArrays(1, &obj.vao);
        glBindVertexArray(obj.vao);

        glGenBuffers(1, &obj.vbo);
        glBindBuffer(GL_ARRAY_BUFFER, obj.vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(GLfloat) * obj.nVerts * N_ATTR_PER_VERT,
                     (GLvoid*) (&obj.vertices[0]),
                     GL_STATIC_DRAW);
    }

    bool init(const Config& config) override {
        RenderPass::init(config);

        // Create shader
        GLuint vs = compileShader("gi.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("gi.fs", GL_FRAGMENT_SHADER);
        shader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create uniforms
        modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
        viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
        projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));

        // Create vertex buffers
        objects.resize(scene.worldData.shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
        }

        return true;
    }

    void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // TODO(A5): Implement this
		//Same as normal.h 
		
		// Define shader to use
		glUseProgram(shader);

		// Update camera
		glm::mat4 model, view, projection;
		camera.Update();
		camera.GetMatricies(projection, view, model);

		// Pass uniforms
		glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
		glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
		glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));
		//glUniformMatrix4fv(normalMatUniform, 1, GL_FALSE, &(normalMat[0][0]));

		// Draw
		for (auto& object : objects) {
			/**
			 * 1) Bind vertex array of current object.
			 * 2) Draw its triangles.
			 * 3) Bind vertex array to 0.
			 */
			 // TODO(A1): Implement this

			glBindBuffer(GL_ARRAY_BUFFER, object.vbo);
			glBindVertexArray(object.vao);
			glDrawArrays(GL_TRIANGLES, 0, object.nVerts);
			glBindVertexArray(0);
		}

        RenderPass::render();
    }

};

TR_NAMESPACE_END
