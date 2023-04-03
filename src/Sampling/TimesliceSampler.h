#ifndef _TIMESLICE_H_
#define _TIMESLICE_H_

#include "bunnykiller.h"

#include "Sampling/Sampler.h"



class TimesliceSampler: public Sampler
{
	int time_resolution;
	Real time_exposure;

	int current_x, current_y;
	int current_time;
    Real limit_x, limit_y;
	int current_subx, current_suby;
	bool end;
public:
	TimesliceSampler(int _size_x, int _size_y, int _sqrt_spp, int _time_resolution, Real _time_exposure) :
		Sampler(_size_x, _size_y, _sqrt_spp),
		time_resolution(_time_resolution), time_exposure(_time_exposure),
        current_x(0), current_y(0), current_time(0),
        limit_x(std::nextafter(Real(1), Real(0.))), limit_y(limit_x),
        current_subx(0), current_suby(0),
		end(false)
	{}

	bool get_next_sample(Sample &sample)
	{
		if (end)
            return false;

        Real shift_x = (current_subx + Random::StdRNG.next_real())*inv_spp;
        Real shift_y = (current_suby + Random::StdRNG.next_real())*inv_spp;

        // Despite the RNG returning values in the [0, 1) range, sometimes
        // the result ends being (current + 1) due to rounding, so we clamp
        Real pos_x = std::min<Real>(current_x + shift_x, limit_x);
        Real pos_y = std::min<Real>(current_y + shift_y, limit_y);

        sample.position = Vector2(pos_x, pos_y);
		sample.weight = 1.;
		sample.current_time = current_time * time_exposure + 0.5f * time_exposure;

		if (++current_subx == spp) {
			current_subx = 0;
			if (++current_suby == spp) {
				current_suby = 0;
				if (single_pixel) {
					if (++current_time == time_resolution) {
						end = true;
						return true;
					}
				}
				
				if (++current_x == size_x) {
					current_x = 0;
					if (single_scanline) {
						if (++current_time == time_resolution) {
							end = true;
							return true;
						}
					}

					if (++current_y == size_y) {
						current_y = 0;
						if (++current_time == time_resolution) {
							end = true;
							return true;
						}
					}
				}
                limit_x = std::nextafter(Real(current_x + 1), Real(0.));
                limit_y = std::nextafter(Real(current_y + 1), Real(0.));
			}
		}

		return true;
	}

	size_t get_nb_samples() const
	{
		return (single_pixel ?
				size_t(spp)*size_t(spp) :
				size_t(size_x)*size_t(size_y)*size_t(spp)*size_t(spp)) * size_t(time_resolution);
	}

	void restart()
	{
		end = false; 
		current_time = 0;
		if (single_pixel) {
			current_x = single_pixel_x;
			current_y = single_pixel_y;
		} else if (single_scanline) {
			current_y = single_pixel_y;
			current_x = current_subx = current_suby = 0;
		} else {
			current_x = current_y = current_subx = current_suby = 0;
		}
        limit_x = std::nextafter(Real(current_x + 1), Real(0.));
        limit_y = std::nextafter(Real(current_y + 1), Real(0.));
	}

	virtual void set_scanline(int y)
	{
		Sampler::set_scanline(y);
		restart();
	}

	virtual void set_single_pixel(int x, int y)
	{
        Sampler::set_single_pixel(x, y);
		restart();
	}
}; // TimesliceSampler

#endif // _TIMESLICE_H_
