#ifndef _LIGHTFUNCTIONS_H_
#define _LIGHTFUNCTIONS_H_

#include "bunnykiller.h"
#include "Color/Spectrum.h"

template <class Radiance>
class BaseLightFunction
{
public:
    virtual Radiance getRadianceByTime(Real time) = 0;
    virtual bool isOn(Real time) = 0;
    virtual void sample(Radiance &f, Real &t, Real &p) = 0;
};

template <class Radiance>
class BoxWindowFunc : public BaseLightFunction<Radiance>
{
public:
    Real start_time;
    Real end_time;

    BoxWindowFunc(Real _start_time = 0, Real _end_time = std::numeric_limits<Real>::infinity()) : start_time(_start_time), end_time(_end_time) {}

    virtual Radiance getRadianceByTime(Real time) override
    {
        if (start_time <= time && time <= end_time)
            return Radiance(1.);
        return Radiance(0.);
    }

    virtual bool isOn(Real time)
    {
        return start_time <= time && time <= end_time;
    }

    virtual void sample(Radiance &f, Real &t, Real &p)
    {
        f = Radiance(1.);
        Sampling.jittered(t, p, start_time, end_time);
    }
};

#endif // _LIGHTFUNCTIONS_H_