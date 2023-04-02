#ifndef _TRANSIENT_RENDERER_H_
#define _TRANSIENT_RENDERER_H_

#include "bunnykiller.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "RenderEngine.h"
#include "LightSource/LightSample.h"
#include "Film/StreakCameraFilm.h"
#include "Utils/Timer.h"

/** Class implementing transient rendering engine. 
	It contains the rendering main loop, and supports 
	non-progressive time-resolved rendering. 
	Note that it needs StreakCameraFilm to work properly. */
template<unsigned D, class Radiance, class RadianceAttenuation>
class MovingTransientRenderer : public RenderEngine<D, Radiance, RadianceAttenuation>
{
protected:
	using RETR = RenderEngine<D, Radiance, RadianceAttenuation>;
	using FilmR = typename RETR::FilmR;
	using RadianceSampleR = typename RETR::RadianceSampleR;
	using RadianceSampleRecordR = typename RETR::RadianceSampleRecordR;
	using RadianceSampleRecordVectorR = typename RETR::RadianceSampleRecordVectorR;
	using IntegratorR = typename RETR::IntegratorR;
	using WorldR = typename RETR::WorldR;
protected:
	bool time_sampling;
public:
	bool sensor_mode;
public:
	MovingTransientRenderer() :
		RETR(nullptr), time_sampling(false), sensor_mode(false)
	{}

	MovingTransientRenderer(FILE *_f_log, bool _time_sampling = false) :
		RETR(_f_log), time_sampling(_time_sampling), sensor_mode(false)
	{}

	virtual ~MovingTransientRenderer() {}

	virtual void render(const char *name, WorldR& world, IntegratorR *integrator,
						const Camera<D>& camera, FilmR *film, Sampler *sampler) const;

	void set_sensor_mode (bool _sensor_mode)
	{
		sensor_mode = _sensor_mode;
	}
}; // MovingTransientRenderer

template<unsigned D, class Radiance, class RadianceAttenuation>
void MovingTransientRenderer<D,Radiance,RadianceAttenuation>::render(const char *name, WorldR& world,
		IntegratorR *integrator, const Camera<D>& camera, FilmR *film, Sampler *sampler) const
{
	if (RETR::f_log) {
		fprintf(RETR::f_log, "RENDER\n");
		fprintf(RETR::f_log, "======\n");
	}
	// ----------------------------------------------------------------------
	// If no integrator has been defined, use the default...	
	if (!integrator) {
		throw std::runtime_error("Error: Need to define an integrator for rendering\n");
	}

	// This rendered can only render into streak films. Thus, it is needed to
	// test it.
	//if( static_cast<StreakCameraFilm*>(film) == 0 )
	//	throw("Transient Renderer needs an StreakCameraFilm");

	// ----------------------------------------------------------------------
	// If needed, preprocess integrator ...
	printf("Preprocessing integrator...\r");
	Timer timer;
	timer.start();

	integrator->preprocess();
	{
		Real secs = timer.get_secs();
		unsigned hours = (unsigned)(secs)/3600;
		secs -= hours*3600;
		unsigned minutes = (unsigned)(secs)/60;
		secs -= minutes*60;

		printf("Preprocessed integrator:  [%d:%d:%f]\n", hours, minutes, secs);

		if (RETR::f_log) {
			fprintf(RETR::f_log, "Preprocessed integrator: %d:%d:%f\n", hours, minutes, secs);
		}
	}

	// ----------------------------------------------------------------------
	// Go rendering!
	printf("Rendering...\r");
	timer.start();

	// Loop over all samples in the image
	Sample film_sample;
	size_t nb_samples = sampler->get_nb_samples(),
		   curr_sample = 0;
	double tstamp = -1.;

	const PolarizationFrame<D> camera_frame = camera.get_frame();
	RadianceSampleRecordVectorR samples_rec;

	// TODO: 增加枚举每帧时间功能
	while (sampler->get_next_sample(film_sample)) {
		// Update timer
		double secs = timer.get_secs();
		if ((secs - tstamp) > .5) {
			tstamp = secs;
			unsigned hours = (unsigned)(secs)/3600;
			secs -= hours*3600;
			unsigned minutes = (unsigned)(secs)/60;
			secs -= minutes*60;

			double perc = double(curr_sample)/double(nb_samples)*100;
			printf("Rendering ... %f%%:   [%d:%d:%f]             \r", perc, hours, minutes, secs);
			fflush(stdout);
		}
		curr_sample++;

		VectorN<D-1> ic;
		Vector2 meh = film->window_coords2image_coords(film_sample.position);
		
		if (D == 2) {
			ic[0] = meh[0];
		} else {
			ic[0] = meh[0];
			ic[1] = meh[1];
		}
		
		// Get the camera sample ray...
		Ray<D> r = camera.get_ray(ic);

		// ...trace the ray, compute time-resolved samples, and store them...
        world.Li(r, integrator, samples_rec, film->get_time_length(), film->get_time_resolution());

		// ...discard NaNs and infinities and align
		for (RadianceSampleR& sample : samples_rec.samples) {
			if (!sample.radiance.is_finite())
				sample.radiance = Radiance(0.);
			align_to_frame(camera_frame, sample.radiance);
		}

		// ...and finally store the samples in the film
		film->add_samples(film_sample, samples_rec);

		// No point on shrinking
		samples_rec.samples.reserve(samples_rec.samples.capacity());
		samples_rec.samples.clear();
	}

	//----------------------------------------------------------------------
	// End rendering, stop timer
	{
		double secs = timer.get_secs();
		unsigned hours = (unsigned)(secs)/3600;
		secs -= hours*3600;
		unsigned minutes = (unsigned)(secs)/60;
		secs -= minutes*60;

		printf("Rendering ... [DONE]:   [%d:%d:%f]\n", hours, minutes, secs);
		fflush(stdout);
 
		if (RETR::f_log) {
			fprintf(RETR::f_log, "Total time: %d:%d:%f\n", hours, minutes, secs);
			fflush(RETR::f_log);
		}
	}
	
	// Save the image
	film->write(name);
}

#endif //_TRANSIENT_RENDERER_H_
