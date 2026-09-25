# FDAS_demo
A small piece of code written by Charlie Li in his summer research, aiming at using **Fourier Domain Acceleration Search (FDAS)** method to search for pulsars.
FDAS pulsar search uses information in fourier domain to do matched filtering, and it can boost the signal to noise ratio significantly.
This code is written as a *science research training* and may not reach the standard of a professional code. The writer overcame a lot of difficulties as he pushed the code.
So he writes all his *confusion* and *solution* below, in which may be helpful to those who want to learn and write a similar code. He hopes it's valuable for beginners like him.
The code contains about 3 weeks of work.

## directory
[math method detail](#the-math-algorithm)

## The math algorithm
### Fourier Transform (FT) 
Fourier Transform is the math method of doing inner product of signal and periodic base function to gain the value of a Harmonic component.

$Spectrum(\omega) = Norm_factor * \int_{-\inf}^{\inf} Signal(t^*)*exp(1j*\omega*t^*)dt^*$

**Discrete Fourier Transform (DFT)** is a discreter version of FT, and it can applied to computer program. And **Fast Fourier Transfrom** can accelerate the calculate of DFT. 

DFT is the same to the process below:
- given the $\omega$ bin.
- rotate each bin in time domain by angle of $\omega * t$ in complex plane.
- sum all the after-rotation time bin complex value in complex plane.
- given the Spectrum value of $\omega$ bin. 

If there is a **periodic signal**, all the peak value will add together at exactly the $omega$ and cause an **interfence enhancement** which is a peak in Fourier domain.
And if the $omega$ doesn't match, there would be a Interference cancellation and no peaks shows off.

In this view:
- A *time shift* will cause t->t+dt and makes a rotation on phase in complex plane. see [here](#rotation-caused-by-time-shift)
- 

**If we do the FT and just extract the value around the harmonics, noise which is not periodic will be filtered, and S/N will be boosted.**






## Appendix

### Rotation caused by time shift
Let's assume some signal become a rect-like function in fourier domain. We simplify the problem by consider only 3 bin,

$(\omega-d\omega)*t, (\omega)*t, (\omega+d\omega)*t$

If we shift the time to t+dt,

$(\omega-d\omega)*(t+dt), (\omega)*(t+dt), (\omega+d\omega)*(t+dt)$

There is not only a common $\delta \phi = \omega * dt$, but also a small rotation phase $-d\omega * dt, 0, d\omega * dt$. It's a second order term.
