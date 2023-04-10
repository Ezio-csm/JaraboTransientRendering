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
	std::vector<std::array<Imaging::Image<Real>, Radiance::components>> m_available_slices;
	size_t m_first_available_slice;
	size_t m_last_available_slice;

	/* Buffers for optional variance estimation */
	std::vector<Imaging::Image<Real>> m_available_slices_variance;
	std::vector<Imaging::Image<Real>> m_available_slices_nsamples;

	/* Weight buffer for every slice*/
	std::vector<Imaging::Image<Real>> m_weight_slices;
protected:
	void advance_timeline();

	void average_and_save_slice(const unsigned int y, const char* name = nullptr);

	void get_film_samples(const Real x, const unsigned int d, std::vector<Real>& samples) const
	{
		if (d < 2) {
			FTR::get_film_samples(x, d, samples);
		} else {
			samples.clear();

			int size_kernel = FTR::m_filter->get_size();
			int half_size = size_kernel / 2;

			int t = (int) (std::floor(x));

			unsigned int init_t = std::max<int>(0, t - half_size);
			unsigned int end_t = std::min<int>(t + half_size, int(m_frame_num - 1));

			for (size_t i = init_t; i <= end_t; i++) {
				samples.push_back(Real(i) + .5);
			}
		}
	}

public:
	TimesliceFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
		bool camera_unwarp = false, Filter*_filter = nullptr,
		FilmComponents comp = FilmComponents::RADIANCE) :
			FTR(w, h, _filter, comp),
			m_frame_num(t),
			m_exposure_time(exposure),
			m_offset(0.),
			camera_unwarp(camera_unwarp),
			time_filter(FTR::m_filter->get_subfilter(2)),
			m_first_available_slice(0),
			m_last_available_slice(0)
	{
		size_t available_slices = ((time_filter) ? time_filter->get_size() : 1);

		for (size_t i = 0; i < available_slices; i++) {

			m_available_slices.emplace_back(std::array<Imaging::Image<Real>,
					Radiance::components>());
			m_available_slices.back().fill(Imaging::Image<Real>(FTR::width, FTR::height, Radiance::spectral_samples));

			if (FTR::m_components && FilmComponents::VARIANCE) {
				m_available_slices_variance.emplace_back(Imaging::Image<Real>(
						FTR::width, FTR::height, Radiance::spectral_samples));
				m_available_slices_nsamples.emplace_back(Imaging::Image<Real>(
						FTR::width, FTR::height, Radiance::spectral_samples));
			}
		}

		m_last_available_slice = m_first_available_slice + available_slices - 1;
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
		std::vector<Real> temporal_samples;
		FTR::get_film_samples(sample.position, film_samples);
		get_film_samples(sample.current_time, 2, temporal_samples);

		for (const Real& sample_t : temporal_samples) {
			unsigned int t = (unsigned int) (std::floor(sample_t));
			while (t > m_last_available_slice) {
				average_and_save_slice(m_first_available_slice, FTR::m_name.c_str());
				advance_timeline();
			}
			size_t index = t - m_first_available_slice;

			for (const Vector2& s : film_samples) {
				unsigned int x = (unsigned int) (std::floor(s[0]));
				unsigned int y = (unsigned int) (std::floor(s[1]));

				/* Apply spatial filter */
				Real w = FTR::m_filter->evaluate(s - sample.position) * sample.weight;
				m_weight_slice[index].add(x, y, &w);

				/* Iterate over samples */
				for (const RadianceSampleR& r : rec.samples) {
					/* Calculate sample's time coordinate */
					Real cam_time = camera_unwarp ? (r.time - rec.distance) : r.time;
					Real pos_sensor = (cam_time - m_offset) / m_exposure_time;

					/* If out of the range of the sensor, discard */
					if (!(pos_sensor < m_frame_num && pos_sensor >= 0.0))
						continue;

					add_temporal_sample(r, w, x, y, sample_t, pos_sensor);

					/* And store it into the accumulated image */
					FTR::draw_pixel(x, y, r.radiance * w, -1.);
				}

				/* Store components if requested */
				FTR::add_components(x, y, rec, w);
			}
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
		for (size_t i = m_first_available_slice; i <= m_last_available_slice; i++) {
			average_and_save_slice(i);
		}

		FTR::average(comp);
	}

	virtual void write(const char* name) override
	{
		name = (name) ? name : FTR::m_name.c_str();

		for (size_t i = m_first_available_slice; i <= m_last_available_slice; i++) {
			average_and_save_slice(i, name);
		}

		FTR::write(name);
	}

public:
	Imaging::Image<Real>& streak_image(unsigned int y = 0, unsigned int i = 0)
	{
		size_t index = (size_t)std::floor(y) - m_first_available_slice;

		return m_available_slices[index][i];
	}
}; // TimesliceFilm

template<unsigned D, class Radiance>
void TimesliceFilm<D, Radiance>::advance_timeline()
{
	/* Clean the streak and move it to the back */
	for (size_t i = 0; i < Radiance::components; i++) {
		(*m_available_slices.begin())[i].clean();
	}
	(*m_weight_slices.begin()).clean();

	std::rotate(m_available_slices.begin(),
		m_available_slices.begin()++, m_available_slices.end());

	/* If present, update the variance estimation too */
	if (FTR::m_components && FilmComponents::VARIANCE) {
		(*m_available_slices_variance.begin()).clean();
		std::rotate(m_available_slices_variance.begin(),
				m_available_slices_variance.begin()++, m_available_slices_variance.end());

		(*m_available_slices_nsamples.begin()).clean();
		std::rotate(m_available_slices_nsamples.begin(),
				m_available_slices_nsamples.begin()++, m_available_slices_nsamples.end());
	}

	m_first_available_slice = std::min<size_t>(m_first_available_slice + 1, m_frame_num - 1);
	m_last_available_slice = std::min<size_t>(m_last_available_slice + 1, m_frame_num - 1);
}

template<unsigned D, class Radiance>
void TimesliceFilm<D, Radiance>::add_temporal_sample(const RadianceSampleR& r, Real w,
	unsigned int x, unsigned int y, Real tr, Real ts, bool filtered)
{
	size_t index = (size_t)std::floor(tr) - m_first_available_slice;

	Radiance rad = r.radiance;

	if (w > 0.0) {
		rad *=  w;
	}

	if (filtered) {
		rad *= time_filter->evaluate(ts - tr);
	};

	m_available_slices[index][0].add(x, y, &rad[0], Radiance::spectral_samples);

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance */
		for (size_t c = 0; c < Radiance::spectral_samples; c++) {
			Real curr_val = m_available_slices[index][0](x, y, 1);
			Real& nsamples = m_available_slices_nsamples[index](x, y, 1);
			Real& var = m_available_slices_variance[index](x, y, 1);

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
	size_t index = (size_t)std::floor(ts) - m_first_available_slice;

	Real wt = 1.0;
	if (filtered) {
		wt = time_filter->evaluate(ts - tr);
	};

	for (size_t c = 0; c < PolarizedLight<3>::components; c++) {
		Spectrum rad = (w > 0.0) ? r.radiance[c] * w * wt : r.radiance[c] * wt;

		m_available_slices[index][c].add(x, y, &rad[0], PolarizedLight<3>::spectral_samples);
	}

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance, over visible light */
		Spectrum rad = (w > 0.0) ? r.radiance[0] * w * wt : r.radiance[0] * wt;

		for (size_t c = 0; c < PolarizedLight<3>::spectral_samples; c++) {
			Real curr_val = m_available_slices[index][0](x, y, 1);
			Real& nsamples = m_available_slices_nsamples[index](x, y, 1);
			Real& var = m_available_slices_variance[index](x, y, 1);

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
	size_t index = (size_t)std::floor(t) - m_first_available_slice;

	std::array<Imaging::Image<Real>, Radiance::components>& out_image = m_available_slices[index];

	Imaging::Image<Real>* out_image_var =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_available_slices_variance[index] : nullptr;
	Imaging::Image<Real>* out_image_ns =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_available_slices_nsamples[index] : nullptr;

	if (!FTR::averaged) {
		/* Apply spatial filter */
		for (size_t c = 0; c < Radiance::components; c++) {
			out_image[c].weight(m_weight_slices[index]);
			
			if (FTR::m_components && FilmComponents::VARIANCE) {
				out_image_var->weight(m_weight_slices[index]);
				out_image_ns->weight(m_weight_slices[index]);
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
