#ifndef _TIMESLICE_FILM_H_
#define _TIMESLICE_FILM_H_

#include "bunnykiller.h"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "Color/Spectrum.h"
#include "Film/Film.h"
#include "Filter/UberFilter.h"
#include "Image/Image.h"
#include "Image/ImageIO.h"

/**
 * @brief The film only store one frame at one time.
 * 
 * @tparam D 
 * @tparam Radiance 
 */
template<unsigned D, class Radiance>
class TimesliceFilm : public Film<D, Radiance>
{
protected:
	using FTR = Film<D, Radiance>;
	using RadianceSampleR = typename FTR::RadianceSampleR;
	using RadianceSampleRecordR = typename FTR::RadianceSampleRecordR;
	using RadianceSampleRecordVectorR = typename FTR::RadianceSampleRecordVectorR;
protected:
	size_t m_frame_num;
	Real m_exposure_time;
	Real m_offset;

	bool camera_unwarp;

	Filter* time_filter;
protected:
	/* Buffers containing the streak images */
	std::array<Imaging::Image<Real>, Radiance::components> m_slice;

	/* Buffers for optional variance estimation */
	Imaging::Image<Real> m_slice_variance;
	Imaging::Image<Real> m_slice_nsamples;

	size_t m_current_slice_id;
protected:
	void advance_timeline();
	void average_and_save_slice(const unsigned int y, const char* name = nullptr);

public:
	TimesliceFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
		bool camera_unwarp = false, Filter*_filter = nullptr, Filter* _time_filter = nullptr,
		FilmComponents comp = FilmComponents::RADIANCE) :
			FTR(w, h, _filter, comp),
			m_frame_num(t),
			m_exposure_time(exposure),
			m_offset(0.),
			camera_unwarp(camera_unwarp),
			time_filter(_time_filter)
	{
		size_t available_slices = ((FTR::m_filter) ? FTR::m_filter->get_size() : 1);

		for (size_t i = 0; i < available_slices; i++) {

			m_slice = std::array<Imaging::Image<Real>, Radiance::components>();
			m_slice.fill(Imaging::Image<Real>(FTR::width, FTR::height, Radiance::spectral_samples));

			if (FTR::m_components && FilmComponents::VARIANCE) {
				m_slice_variance = Imaging::Image<Real>(FTR::width, FTR::height, Radiance::spectral_samples);
				m_slice_nsamples = Imaging::Image<Real>(FTR::width, FTR::height, Radiance::spectral_samples);
			}
		}
		m_current_slice_id = (unsigned int) (std::floor(m_offset / m_exposure_time));
	}

	virtual ~TimesliceFilm()
	{
	}

	virtual unsigned int get_time_resolution() const override
	{
		return m_frame_num;
	}

	virtual void add_sample(const Sample& sample, const RadianceSampleRecordR& rec) override
	{
		add_samples(sample,
				RadianceSampleRecordVectorR( { rec.sample }, rec.distance, rec.pos, rec.normal));
	}

	virtual void add_samples(const Sample& sample, const RadianceSampleRecordVectorR& rec) override
	{
		FTR::averaged = false;

		/* First the pixel coverage is calculated */
		std::vector<Vector2> film_samples;
		FTR::get_film_samples(sample.position, film_samples);

		/* We iterate over all pixels */
		std::vector<Real> temporal_samples;

		unsigned int t = (unsigned int) (std::floor(sample.current_time / m_exposure_time));
		/* Check if new timeline is needed */
		while (t > m_current_slice_id) {
			average_and_save_slice(m_current_slice_id, FTR::m_name.c_str());
			advance_timeline();
		}
		for (const Vector2& s : film_samples) {
			unsigned int x = (unsigned int) (std::floor(s[0]));
			unsigned int y = (unsigned int) (std::floor(s[1]));

			/* Apply spatial filter */
			Real w = FTR::m_filter->evaluate(s - sample.position) * sample.weight;
			FTR::m_weight.add(x, y, &w);

			/* Iterate over samples */
			for (const RadianceSampleR& r : rec.samples) {
				/* Calculate sample's time coordinate */
				Real cam_time = camera_unwarp ? (r.time - rec.distance) : r.time;
				Real pos_sensor = (cam_time - m_offset)/ m_exposure_time;

				/* If out of the range of the sensor, discard */
				if (!(pos_sensor < m_frame_num && pos_sensor >= 0.0))
					continue;
				/* Finally, just draw the radiances on the streak film, multiplied by the
				 * spatio-temporal filtering weight!. Note that the temporal weight doesn't act as
				 * a filtering weight, but as a PSF. Thus, is additive.
				 */
				add_temporal_sample(r, w, x, y, sample.current_time / m_exposure_time, pos_sensor);
				/* And store it into the accumulated image */
				FTR::draw_pixel(x, y, r.radiance * w, -1.);
			}

			/* Store components if requested */
			FTR::add_components(x, y, rec, w);
		}
	}

	void add_temporal_sample(const RadianceSampleR& r, Real w, unsigned int x, unsigned int y,
		Real tr, Real ts, bool filtered = true);

	virtual Real get_time_length() const override
	{
		return Real(m_frame_num) * m_exposure_time + m_offset;
	}

	virtual Real get_exposure_time() const override
	{
		return m_exposure_time;
	}

	void set_offset(Real offset)
	{
		m_offset = offset;
	}

	virtual Real get_time_offset() const override
	{
		return m_offset;
	}

	void set_streak(unsigned int y)
	{
		ptrdiff_t delta = m_last_available_slice - m_first_available_slice;
		m_first_available_slice = y;
		m_last_available_slice = y + delta;
	}

	virtual void average(FilmComponents comp = FilmComponents::ALL) override
	{
		average_and_save_slice(m_current_slice_id);
		FTR::average(comp);
	}

	virtual void write(const char* name) override
	{
		name = (name) ? name : FTR::m_name.c_str();
		average_and_save_slice(m_current_slice_id, name);
		FTR::write(name);
	}

public:
	Imaging::Image<Real>& streak_image(unsigned int y = 0, unsigned int i = 0)
	{
		size_t index = (size_t)std::floor(y) - m_first_available_slice;

		return m_slice[index][i];
	}
}; // TimesliceFilm

template<unsigned D, class Radiance>
void TimesliceFilm<D, Radiance>::advance_timeline()
{
	for (size_t i = 0; i < Radiance::components; i++)
		m_slice[i].clean();
	m_slice_variance.clean();
	m_slice_nsamples.clean();
		
	m_current_slice_id++;
	FTR::m_weight.clean();
}

template<unsigned D, class Radiance>
void TimesliceFilm<D, Radiance>::add_temporal_sample(const RadianceSampleR& r, Real w,
	unsigned int x, unsigned int y, Real tr, Real ts, bool filtered)
{
	Radiance rad = r.radiance;

	if (w > 0.0) {
		rad *=  w;
	}

	if (filtered) {
		rad *= time_filter->evaluate(ts - tr);
	};

	m_slice[0].add(x, y, &rad[0], Radiance::spectral_samples);

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance */
		for (size_t c = 0; c < Radiance::spectral_samples; c++) {
			Real curr_val = m_slice[0](x, y, 1);
			Real& nsamples = m_slice_nsamples(x, y, 1);
			Real& var = m_slice_variance(x, y, 1);

			/* mean = current_value/nsamples */
			Real mean = (nsamples) ? curr_val/nsamples : Real(0.0);
			nsamples++;

			Real delta_rad = rad[c] - mean;
			mean = mean + delta_rad / nsamples;

			var += delta_rad * (rad[c] - mean);
		}
	}
}

template<>
void TimesliceFilm<3, PolarizedLight<3>>::add_temporal_sample(const RadianceSampleR& r, Real w,
		unsigned int x, unsigned int y, Real tr, Real ts, bool filtered)
{
	Real wt = 1.0;
	if (filtered) {
		wt = time_filter->evaluate(ts - tr);
	};

	for (size_t c = 0; c < PolarizedLight<3>::components; c++) {
		Spectrum rad = (w > 0.0) ? r.radiance[c] * w * wt : r.radiance[c] * wt;

		m_slice[c].add(x, y, &rad[0], PolarizedLight<3>::spectral_samples);
	}

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance, over visible light */
		Spectrum rad = (w > 0.0) ? r.radiance[0] * w * wt : r.radiance[0] * wt;

		for (size_t c = 0; c < PolarizedLight<3>::spectral_samples; c++) {
			Real curr_val = m_slice[0](x, y, 1);
			Real& nsamples = m_slice_nsamples(x, y, 1);
			Real& var = m_slice_variance(x, y, 1);

			/* mean = current_value/nsamples */
			Real mean = (nsamples) ? curr_val/nsamples : Real(0.0);
			nsamples++;

			Real delta_rad = rad[c] - mean;
			mean = mean + delta_rad / nsamples;

			var += delta_rad * (rad[c] - mean);
		}
	}
}

template<unsigned D, class Radiance>
void TimesliceFilm<D, Radiance>::average_and_save_slice(const unsigned int t, const char* name)
{
	std::array<Imaging::Image<Real>, Radiance::components>& out_image = m_slice;
	Imaging::Image<Real>* out_image_var =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_slice_variance : nullptr;
	Imaging::Image<Real>* out_image_ns =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_slice_nsamples : nullptr;

	if (!FTR::averaged) {
		/* Apply spatial filter */
		for (size_t c = 0; c < Radiance::components; c++) {
			out_image[c].weight(FTR::m_weight);
		
			if (FTR::m_components && FilmComponents::VARIANCE) {
				out_image_var->weight(FTR::m_weight);
				out_image_ns->weight(FTR::m_weight);
			}
		}
	}

	/* Buffer to store the final name */
	char nimg[2048];

	/* Store streak */
	if (name) {
		for (size_t c = 0; c < Radiance::components; c++) {
			if (Radiance::components > 1) {
				sprintf(nimg, "%s_s%d_%04d.%s", name, unsigned(c), t, FTR::m_extension.c_str());
			} else {
				sprintf(nimg, "%s_%04d.%s", name, t, FTR::m_extension.c_str());
			}

			/* Workaround for hdr images less than 8 pixels wide (matlab fails reading them),
			 * store them transposed
			 */
			if (out_image[c].width() < 8) {
				out_image[c].transpose();
			}

			Imaging::save(out_image[c], nimg);
		}

		if (out_image_var) {
			sprintf(nimg, "%s_%04d_variance.%s", name, t, FTR::m_extension.c_str());
			if (out_image_var->width() < 8) {
				out_image_var->transpose();
			}
			Imaging::save(*out_image_var, nimg);
		}

		if (out_image_ns) {
			sprintf(nimg, "%s_%04d_nsamples.%s", name, t, FTR::m_extension.c_str());
			if (out_image_ns->width() < 8) {
				out_image_ns->transpose();
			}
			Imaging::save(*out_image_ns, nimg);
		}
	}
}

#endif // _TIMESLICE_FILM_H_
