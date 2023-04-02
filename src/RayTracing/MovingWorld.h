#ifndef _MOVINGWORLD_H_
#define _MOVINGWORLD_H_

#include "bunnykiller.h"

#include "RayTracing/World.h"

#include <vector>

#ifdef _USE_EMBREE_
#if defined(__GNUC__) || defined(__GNUG__)
#include <x86intrin.h>
#elif defined(_MSC_VER)
#include <intrin.h>
#endif

#define EMBREE_STATIC_LIB
#include "External/Embree/include/embree2/rtcore.h"
#include "External/Embree/include/embree2/rtcore_ray.h"
#endif // _USE_EMBREE_

template<unsigned D, class Radiance>
class MovingWorld : public World<D, Radiance>
{
protected:
    struct Object
    {
        RTCScene g_scene;
        VectorN<D> velocity;
        Real start_time;
        Real end_time;
    };

    std::vector<Object> objects;
    
    Real light_speed = 1.;

public:
    MovingWorld() :
	{}

	//-----------------------------------------------------------
	// Geometry
#ifdef _USE_EMBREE_
	void add_triangle_mesh(const std::string &name_file, Material3D *mat, const VectorN<D> &velocity,
            Real start_time = 0, Real end_time = std::numeric_limits<Real>::infinity())
	{
		std::ifstream fin(name_file.c_str());
		if (!fin.is_open()) {
			throw std::runtime_error("Error: failed to open model file: \"" + name_file + "\"");
		}

		/* Parsing aux variables */
		std::string line;
		unsigned line_num;

		/* Triangle data */
		struct obj_Triangle {
			unsigned v1, v2, v3;
			int vt1, vt2, vt3;
			int vn1, vn2, vn3;

			obj_Triangle(unsigned v1, unsigned v2, unsigned v3) :
				v1(v1), v2(v2), v3(v3),
				vt1(-1), vt2(-1), vt3(-1),
				vn1(-1), vn2(-1), vn3(-1)
			{}

			void set_vt(unsigned _vt1, unsigned _vt2, unsigned _vt3)
			{
				vt1 = _vt1; vt2 = _vt2; vt3 = _vt3;
			}

			void set_vn(unsigned _vn1, unsigned _vn2, unsigned _vn3)
			{
				vn1 = _vn1; vn2 = _vn2; vn3 = _vn3;
			}
		};
		

		struct rtc_Vertex
		{
			Vector3f pos;
			Vector3f normal;
			Vector2f uv;

			rtc_Vertex(float x, float y, float z) :
				pos(x, y, z),
				normal(NAN, NAN, NAN),
				uv(NAN, NAN)
			{}
		};

		std::vector<obj_Triangle> obj_triangles;
		std::vector<rtc_Vertex> vertices;
		std::vector<Vector3f> normals;
		std::vector<Vector2f> uvs;

		/* OBJ indices either start at 1 or at the end of the current sequence.
		   Convert them to C++ indices during parsing. Also check bounds and
		   abort with badly formatted files */
		auto obj_index  = [&](int index, size_t last_index) -> unsigned {
			int cpp_index = (index > 0) ? (index - 1) : (last_index + index);
			if (cpp_index > (int)(last_index - 1) || cpp_index < 0) {
				throw std::out_of_range("Error: reading open model file \"" + name_file +
				                        "\": out of range index " + std::to_string(index) +
					                " at line " + std::to_string(line_num));
			}
			return (unsigned) cpp_index;
		};

		auto obj_index_v = [&](int index) -> unsigned {
			return obj_index(index, vertices.size());
		};

		auto obj_index_vt = [&](int index) -> unsigned {
			return obj_index(index, uvs.size());
		};

		auto obj_index_vn = [&](int index) -> unsigned {
			return obj_index(index, normals.size());
		};

		/* Parse termination check */
		auto check_trail = [&](char& trail) {
			if (trail != '\0') {
				throw std::runtime_error("Error: reading open model file \"" +
						         name_file + "\" at line " +
						         std::to_string(line_num));
			}
		};

		/* Parse line by line. We need to check all the possible OBJ variations
		   in vertex and face definitions. We do not check the validity of any
		   other data in the file. */
		for (line_num = 1; std::getline(fin, line); line_num++) {
			/* Last parsed item buffer */
			char trail = '\0';

			/* Skip comments */
			if (sscanf(line.c_str(), "#%c", &trail) == 1) {
				continue;
			}

			/* Store vertex if found */
			float x, y, z;
			if (sscanf(line.c_str(), "v %f %f %f %*f %c", &x, &y, &z, &trail) >= 4) {
				check_trail(trail);
				vertices.push_back(rtc_Vertex(x, y, z));
				continue;
			}

			if (sscanf(line.c_str(), "v %f %f %f %c", &x, &y, &z, &trail) >= 3) {
				check_trail(trail);
				vertices.push_back(rtc_Vertex(x, y, z));
				continue;
			}

			/* Store texture coordinate if found */
			float u, v;
			if (sscanf(line.c_str(), "vt %f %f %*f %c", &u, &v, &trail) >= 3) {
				check_trail(trail);
				/* Arbitrary UV convention which matches Blender */
				uvs.push_back(Vector2f(u, 1. - v));
				continue;
			}

			if (sscanf(line.c_str(), "vt %f %f %c", &u, &v, &trail) >= 2) {
				check_trail(trail);
				/* Arbitrary UV convention which matches Blender */
				uvs.push_back(Vector2f(u, 1. - v));
				continue;
			}

			/* Store vertex normal if found */
			if (sscanf(line.c_str(), "vn %f %f %f %c", &x, &y, &z, &trail) >= 3) {
				check_trail(trail);
				normals.push_back(Vector3f(x, y, z));
				continue;
			}

			/* Ignore other vertex data */
			if (sscanf(line.c_str(), "vp %*f %c", &trail) >= 1) {
				check_trail(trail);
				continue;
			}

			if (sscanf(line.c_str(), "vp %*f %*f %c", &trail) >= 2) {
				check_trail(trail);
				continue;
			}

			if (sscanf(line.c_str(), "vp %*f %*f %*f %c", &trail) >= 3) {
				check_trail(trail);
				continue;
			}

			/* Reject malformed vertex data */
			if (sscanf(line.c_str(), "v%c", &trail) == 1) {
				throw std::runtime_error("Error: reading open model file \"" +
							 name_file + "\" at line " +
							 std::to_string(line_num));
			}

			/* Triangle with just geometric information */
			int v1, v2, v3;
			if (sscanf(line.c_str(), "f %d %d %d %c", &v1, &v2, &v3, &trail) >= 3) {
				check_trail(trail);
				obj_triangles.push_back(
					obj_Triangle(obj_index_v(v1), obj_index_v(v2), obj_index_v(v3)));
				continue;
			}				

			/* Triangle with texture indices */
			int vt1, vt2, vt3;
			if (sscanf(line.c_str(), "f %d/%d %d/%d %d/%d %c",
					&v1, &vt1, &v2, &vt2, &v3, &vt3, &trail) >= 6) {
				check_trail(trail);
				obj_triangles.push_back(
					obj_Triangle(obj_index_v(v1), obj_index_v(v2), obj_index_v(v3)));
				obj_triangles.back().set_vt(
					obj_index_vt(vt1), obj_index_vt(vt2), obj_index_vt(vt3));
				continue;
			}

			/* Triangle with normal indices */
			int vn1, vn2, vn3;
			if (sscanf(line.c_str(), "f%d//%d %d//%d %d//%d %c",
					&v1, &vn1, &v2, &vn2, &v3, &vn3, &trail) >= 6) {
				check_trail(trail);
				obj_triangles.push_back(
					obj_Triangle(obj_index_v(v1), obj_index_v(v2), obj_index_v(v3)));
				obj_triangles.back().set_vn(
					obj_index_vn(vn1), obj_index_vn(vn2), obj_index_vn(vn3));
				continue;
			}

			/* Triangle with texture and normal indices */
			if (sscanf(line.c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %c",
					&v1, &vt1, &vn1, &v2, &vt2, &vn2, &v3, &vt3, &vn3, &trail) >= 9) {
				check_trail(trail);
				obj_triangles.push_back(
					obj_Triangle(obj_index_v(v1), obj_index_v(v2), obj_index_v(v3)));
				obj_triangles.back().set_vt(
					obj_index_vt(vt1), obj_index_vt(vt2), obj_index_vt(vt3));
				obj_triangles.back().set_vn(
					obj_index_vn(vn1), obj_index_vn(vn2), obj_index_vn(vn3));
				continue;
			}

			/* Reject not-triangular meshes */
			if (sscanf(line.c_str(), "f%c", &trail) == 1) {
				throw std::runtime_error("Error: reading open model file \"" +
							 name_file + "\": invalid face at line " +
							 std::to_string(line_num));
			}
		}

		/* Triangle buffer */
		struct rtc_Triangle
		{
			int v[3];
		};

		/* Allocate triangl buffer */
		size_t n_triangles = obj_triangles.size();

		RTCORE_ALIGN(16)
		rtc_Triangle* rtc_triangles = (rtc_Triangle*) calloc(n_triangles, sizeof(rtc_Triangle));

		/* Copy all the triangles into the Embree buffer and correct if necessary */
		for (size_t i = 0; i < n_triangles; i++) {
			/* Extract vertex indices */
			rtc_Triangle& rtc_tri = rtc_triangles[i];
			obj_Triangle& obj_tri = obj_triangles[i];

			rtc_tri.v[0] = obj_tri.v1;
			rtc_tri.v[1] = obj_tri.v2;
			rtc_tri.v[2] = obj_tri.v3;
			
			/* Geometric normal by Embree's convention */
			Vector3f p0 = vertices[rtc_tri.v[0]].pos;
			Vector3f p1 = vertices[rtc_tri.v[1]].pos;
			Vector3f p2 = vertices[rtc_tri.v[2]].pos;

			Vector3f ng = cross(p2 - p0, p1 - p0).normalized();

			/* Assign texture coordinates */
			std::array<Vector2f, 3> vt = {-1., -1., -1.};

			/* Check if face has texture coordinates. In case is has, assign them,
			   otherwise leave them as negative coordinates */
			if (obj_tri.vt1 >= 0) {
				vt = {uvs[obj_tri.vt1], uvs[obj_tri.vt2], uvs[obj_tri.vt3]};
			}

			/* Embree only allows one set of texture coords per vertex, so when two
			   triangles doesn't shared the same texture coords we need to duplicate
			   the vertex */
			for (size_t j = 0; j < 3; j++) {
				rtc_Vertex& vj = vertices[rtc_tri.v[j]];

				if (vj.uv.has_nans()) {
					/* Vertex not assigned yet */
					vj.uv = vt[j];
				} else if (vj.uv == vt[j]) {
					/* Assigned and shared */
				} else {
					/* Unshared texture coordinate in shared vertex, create new vertex */
					rtc_Vertex v = vertices[rtc_tri.v[j]];
					v.uv = vt[j];
					vertices.push_back(v);

					rtc_tri.v[j] = (int)(vertices.size() - 1);
				}
			}

			/* Assign shading normals */
			std::array<Vector3f, 3> vn = {ng, ng, ng};

			/* Check if face has shading normals. In case is has, assign them,
			   otherwise leave them as the geometric normal */
			if (obj_tri.vn1 >= 0) {
				vn = {normals[obj_tri.vn1].normalized(),
					  normals[obj_tri.vn2].normalized(),
				      normals[obj_tri.vn3].normalized()};
			}

			/* Embree only allows one set of shading normals per vertex, so when two
			   triangles doesn't shared the same shading normals we need to duplicate
			   the vertex */
			for (size_t j = 0; j < 3; j++) {
				rtc_Vertex& vj = vertices[rtc_tri.v[j]];

				if (vj.normal.has_nans()) {
					/* Vertex not assigned yet */
					vj.normal = vn[j];
				} else if (vj.normal == vn[j]) {
					/* Assigned and shared */
				} else {
					/* Unshared shading normal in shared vertex, create new vertex */
					rtc_Vertex v = vertices[rtc_tri.v[j]];
					v.normal = vn[j];
					vertices.push_back(v);

					rtc_tri.v[j] = (int)(vertices.size() - 1);
				}
			}
		}

		/* Allocate vertex data buffer */
		size_t n_vertices = vertices.size();
		
		RTCORE_ALIGN(16)
		rtc_Vertex* rtc_vertices = (rtc_Vertex*) calloc(n_vertices, sizeof(rtc_Vertex));

		/* Copy vertex data */
		std::copy(vertices.begin(), vertices.end(), rtc_vertices);

        /* Create new object, store velocity */
        Object new_obj;
        new_obj.g_scene = rtcDeviceNewScene(World::g_device, RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT1 | RTC_INTERPOLATE);
        new_obj.velocity = velocity;
        new_obj.start_time = start_time;
        new_obj.end_time = end_time;

		/* Create a Embree triangulated mesh */
		unsigned rtc_mesh = rtcNewTriangleMesh(new_obj.g_scene, RTC_GEOMETRY_STATIC, n_triangles, n_vertices);

		/* Set mesh indices buffer */
		rtcSetBuffer2(new_obj.g_scene, rtc_mesh, RTC_INDEX_BUFFER, rtc_triangles,
			0, sizeof(rtc_Triangle), n_triangles);
		/* Set mesh vertices buffer */
		rtcSetBuffer2(new_obj.g_scene, rtc_mesh, RTC_VERTEX_BUFFER, rtc_vertices,
			0, sizeof(rtc_Vertex), n_vertices);
		/* Set mesh vertex normals buffer */
		rtcSetBuffer2(new_obj.g_scene, rtc_mesh, RTC_USER_VERTEX_BUFFER0, rtc_vertices,
			sizeof(rtc_Vertex::pos), sizeof(rtc_Vertex), n_vertices);
		/* Set mesh vertex texture coordinates buffer */
		rtcSetBuffer2(new_obj.g_scene, rtc_mesh, RTC_USER_VERTEX_BUFFER1, rtc_vertices,
			sizeof(rtc_Vertex::pos) + sizeof(rtc_Vertex::normal), sizeof(rtc_Vertex), n_vertices);

		/* Set mesh user data */
		rtcSetUserData(new_obj.g_scene, rtc_mesh, (void *)mat);
	}

	inline void freeze()
	{
		/* commit changes to scene */
		rtcCommit(World::g_scene);
        for(auto &obj:objects)
            rtcCommit(obj.g_scene);

		for (unsigned int i = 0; i < light_source_list.size(); ++i)
			light_source_list[i]->setup();
	}

	inline AABB<D> get_bounding_volume() const
	{
		RTCBounds rtc_bounds;
		rtcGetBounds(World::g_scene, rtc_bounds);
        AABB<D> res(VectorN<D>(rtc_bounds.lower_x,
		                          rtc_bounds.lower_y,
		                          rtc_bounds.lower_z),
		                          VectorN<D>(rtc_bounds.upper_x,
		                                     rtc_bounds.upper_y,
		                                     rtc_bounds.upper_z));
        for(auto &obj : objects)
        {
            rtcGetBounds(obj.g_scene, rtc_bounds);
            res.add(AABB<D>(VectorN<D>(rtc_bounds.lower_x,
		                          rtc_bounds.lower_y,
		                          rtc_bounds.lower_z),
		                          VectorN<D>(rtc_bounds.upper_x,
		                                     rtc_bounds.upper_y,
		                                     rtc_bounds.upper_z)));
        }
	}

	// Return the object that first intersects `ray'
	inline bool first_intersection(Ray<D>& ray, Intersection<D> &it,
			const Real &max_length = std::numeric_limits<Real>::infinity()) const
	{
		/* Initialize ray */
		RTCRay rtc_ray;
		// Ray origin
		rtc_ray.org[0] = (float)ray.get_origin()[0];
		rtc_ray.org[1] = (float)ray.get_origin()[1];
		rtc_ray.org[2] = (float)ray.get_origin()[2];
		rtc_ray.dir[0] = (float)ray.get_direction()[0];
		rtc_ray.dir[1] = (float)ray.get_direction()[1];
		rtc_ray.dir[2] = (float)ray.get_direction()[2];
		rtc_ray.tnear = (float)_SIGMA_VISIBILITY_;
		rtc_ray.tfar = (float)max_length;
		rtc_ray.time = 0.f; // Time for motion blur
		rtc_ray.mask = 0XFFFFFFFF;
		rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
		rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

        RTCRay static_rtc_ray = rtc_ray;
		rtcIntersect(World::g_scene, rtc_ray);
        for(auto obj : objects)
        {
            RTCRay t_rtc_ray = static_rtc_ray;
            VectorN<D> origin = ray.get_origin();
            Real moving_time = ray.get_start_time() - obj.start_time;
            origin = origin + obj.velocity * moving_time;
            t_rtc_ray.org[0] = (float)origin[0];
            t_rtc_ray.org[1] = (float)origin[1];
            t_rtc_ray.org[2] = (float)origin[2];
            VectorN<D> dir = ray.get_direction();
            dir = (dir * light_speed / ray.get_ior() - obj.velocity).normalized();
            t_rtc_ray.dir[0] = (float)dir[0];
            t_rtc_ray.dir[1] = (float)dir[1];
            t_rtc_ray.dir[2] = (float)dir[2];
            rtcIntersect(obj.g_scene, t_rtc_ray);
            Real t_time = t_rtc_ray.tfar / light_speed * ray.get_ior() + ray.get_start_time();
            if(obj.start_time <= t_time && t_time <= obj.end_time && t_rtc_ray.tfar < rtc_ray.tfar)
                rtc_ray = t_rtc_ray;
        }

		if (rtc_ray.geomID != RTC_INVALID_GEOMETRY_ID){
			Vector2f uv_aux;
			rtcInterpolate2(World::g_scene, rtc_ray.geomID, rtc_ray.primID,
				rtc_ray.u, rtc_ray.v, RTC_USER_VERTEX_BUFFER1,
				(float*)&uv_aux, nullptr, nullptr, nullptr, nullptr, nullptr, 2);

			Vector2 uv(uv_aux[0], uv_aux[1]);

			Vector3f normal_aux;
			rtcInterpolate2(World::g_scene, rtc_ray.geomID, rtc_ray.primID,
				rtc_ray.u, rtc_ray.v, RTC_USER_VERTEX_BUFFER0,
				(float*)&normal_aux, nullptr, nullptr, nullptr, nullptr, nullptr, 3);

			VectorN<D> normal(normal_aux[0], normal_aux[1], normal_aux[2]);

			// Store distance in the ray
			ray.set_parameter(rtc_ray.tfar);

			// Store intersection information
			const Material<D>* mat = (const Material<D>*)rtcGetUserData(g_scene, rtc_ray.geomID);

			it.set(ray, mat, normal, uv);
			it.set_coordinate_system();

			return true;
		}
		
		return false;
	}

	inline bool intersects(Ray<D> &ray, const Real &max_length = std::numeric_limits<Real>::infinity(),
			const Real epsilon = _SIGMA_VISIBILITY_) const
	{
		/* initialize ray */
		RTCRay rtc_ray;
		// Ray origin
		rtc_ray.org[0] = (float)ray.get_origin()[0];
		rtc_ray.org[1] = (float)ray.get_origin()[1];
		rtc_ray.org[2] = (float)ray.get_origin()[2];
		rtc_ray.dir[0] = (float)ray.get_direction()[0];
		rtc_ray.dir[1] = (float)ray.get_direction()[1];
		rtc_ray.dir[2] = (float)ray.get_direction()[2];
		rtc_ray.tnear = (float)epsilon;
		rtc_ray.tfar = (float)(max_length - epsilon);
		rtc_ray.time = 0.f; // Time for motion blur
		rtc_ray.mask = 0XFFFFFFFF;
		rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
		rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

		RTCRay static_rtc_ray = rtc_ray;
		rtcIntersect(World::g_scene, rtc_ray);
        for(auto obj : objects)
        {
            RTCRay t_rtc_ray = static_rtc_ray;
            VectorN<D> origin = ray.get_origin();
            Real moving_time = ray.get_start_time() - obj.start_time;
            origin = origin + obj.velocity * moving_time;
            t_rtc_ray.org[0] = (float)origin[0];
            t_rtc_ray.org[1] = (float)origin[1];
            t_rtc_ray.org[2] = (float)origin[2];
            VectorN<D> dir = ray.get_direction();
            dir = (dir * light_speed / ray.get_ior() - obj.velocity).normalized();
            t_rtc_ray.dir[0] = (float)dir[0];
            t_rtc_ray.dir[1] = (float)dir[1];
            t_rtc_ray.dir[2] = (float)dir[2];
            rtcIntersect(obj.g_scene, t_rtc_ray);
            Real t_time = t_rtc_ray.tfar / light_speed * ray.get_ior() + ray.get_start_time();
            if(obj.start_time <= t_time && t_time <= obj.end_time && t_rtc_ray.tfar < rtc_ray.tfar)
                rtc_ray = t_rtc_ray;
        }

		return (rtc_ray.geomID != RTC_INVALID_GEOMETRY_ID);
	}

	inline bool is_visible(const VectorN<D> &v1, const VectorN<D> &v2,
			const Real epsilon = _SIGMA_VISIBILITY_) const
	{
		/* initialize ray */
		RTCRay rtc_ray;
		// Ray origin
		rtc_ray.org[0] = (float)v1[0];
		rtc_ray.org[1] = (float)v1[1];
		rtc_ray.org[2] = (float)v1[2];
		rtc_ray.dir[0] = (float)(v2 - v1)[0];
		rtc_ray.dir[1] = (float)(v2 - v1)[1];
		rtc_ray.dir[2] = (float)(v2 - v1)[2];
		rtc_ray.tnear = (float)epsilon;
		rtc_ray.tfar = (float)((v2 - v1).length() - epsilon);
		rtc_ray.time = 0.f; // Time for motion blur
		rtc_ray.mask = 0XFFFFFFFF;
		rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
		rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

		RTCRay static_rtc_ray = rtc_ray;
		rtcOccluded(World::g_scene, rtc_ray);
        for(auto obj : objects)
        {
            RTCRay t_rtc_ray = static_rtc_ray;
            VectorN<D> origin = ray.get_origin();
            Real moving_time = ray.get_start_time() - obj.start_time;
            origin = origin + obj.velocity * moving_time;
            t_rtc_ray.org[0] = (float)origin[0];
            t_rtc_ray.org[1] = (float)origin[1];
            t_rtc_ray.org[2] = (float)origin[2];
            VectorN<D> dir = ray.get_direction();
            dir = (dir * light_speed / ray.get_ior() - obj.velocity).normalized();
            t_rtc_ray.dir[0] = (float)dir[0];
            t_rtc_ray.dir[1] = (float)dir[1];
            t_rtc_ray.dir[2] = (float)dir[2];
            rtcOccluded(obj.g_scene, t_rtc_ray);
            Real t_time = t_rtc_ray.tfar / light_speed * ray.get_ior() + ray.get_start_time();
            if(obj.start_time <= t_time && t_time <= obj.end_time && t_rtc_ray.tfar < rtc_ray.tfar)
                rtc_ray = t_rtc_ray;
        }

		return (rtc_ray.geomID != RTC_INVALID_GEOMETRY_ID);
	}
#endif // _USE_EMBREE_

    Real time_of_flight(Real distance) const
	{
		return distance / light_speed; /* 0.000299792458f;*/
	}
};

#endif // _MOVINGWORLD_H_
