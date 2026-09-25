# FDAS_demo
A small piece of code written by Charlie Li in his summer research, aiming at using Fourier Domain Acceleration Search (FDAS) method to search for pulsars.
**FDAS** pulsar search uses information in fourier domain to do matched filtering, and it can boost the signal to noise ratio significantly.
This code is written as a *science research training* and may not reach the standard of a professional code. The writer overcame a lot of difficulties as he pushed the code.
So he writes all his *confusion* and *solution* below, in which may be helpful to those who want to learn and write a similar code.
The code contains about 3 weeks of work.

## **The math algorithm**
### **Fourier Transform (FT)** 
Fourier Transform is the math method of doing inner product of signal and periodic base function to gain the value of a Harmonic component.

$$
Spectrum(\omega) = Norm_factor * \int_{-\inf}^{\inf} Signal(t^*)*exp(1j*\omega*t^*)dt^*
$$

**Discrete Fourier Transform (DFT)** is a discreter version of FT, and it can applied to computer program. And **Fast Fourier Transfrom** can accelerate the calculate of DFT. 

DFT is the same to the process below:
- given the $\omega$ bin.
- rotate each bin in time domain by angle of $\omega * t$ in complex plane.
- sum all the after-rotation time bin in complex plane.
- given the Spectrum value of $\Omega$ bin. 
