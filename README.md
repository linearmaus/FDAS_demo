# FDAS_demo
A small piece of code written by Charlie Li in his summer research, aiming at using **Fourier Domain Acceleration Search (FDAS)** method to search for pulsars. FDAS pulsar search uses information in fourier domain to do matched filtering, and it can boost the signal to noise ratio significantly. This code is written as a *science research training* and may not reach the standard of a professional code. The writer overcame a lot of difficulties as he pushed the code. So he writes all his *confusion* and *solution* below, in which may be helpful to those who want to learn and write a similar code. He hopes it's valuable for beginners like him.
*The code contains about 3 weeks of work.*

## Directory
[math method detail](#the-math-algorithm)

## General View
Basically the method is doing such thing:
- doing FFT to the time domain signal, get the spectrum
- locate the fourier bin of the given frequency (f0) and (maybe) frequency derivative (f1)
- calculate the template of funcion base, and do inner product with the spectrum to get the harmonics
- using a inverse-fourier-transform to reconstruct the profile in time domain
- calculate the S/N, which should be boosted
- finding the best S/N to determine the best f0 and f1.

## The math algorithm
### Fourier Transform (FT) 
Fourier Transform is the math method of doing inner product of signal and periodic base function to gain the value of a Harmonic component.

$$
Spectrum(\omega) = Normfactor * \int_{-\inf}^{\inf} Signal(t^*)*e^{i\omega t^*}dt^*
$$

We write this as $Spectrum(\omega) = \hat F(Signal(t))$ where f hat is the Fourier transform. 

**Discrete Fourier Transform (DFT)** is a discreter version of FT, and it can applied to computer program. And **Fast Fourier Transfrom (FFT)** can accelerate the calculation of DFT. 

**DFT is the same to the process below**:
- given the $\omega$ bin.
- rotate each bin in time domain by angle of $\omega * t$ in complex plane.
- sum all the after-rotation time bin complex value in complex plane.
- given the Spectrum value of $\omega$ bin. 

If there is a **periodic signal**, all the peak value will add together at exactly the $omega$ and cause an **interfence enhancement** which is a peak in Fourier domain.
And if the $\omega$ doesn't match, there would be a Interference cancellation and no peaks shows off.

In this view:
- A *time shift* will cause t->t+dt and makes a rotation on phase in complex plane. see [here](#rotation-caused-by-time-shift).
- Assume a perodic delta-function array, this view insures it have value in all harmonics and all value is the same. (Rotation dosen't miss any peak, and everyone is added in the sum).
- DFT may cause "leak" because the frequency bin may not exactly be at the right frequency. But the information does not miss, we can reconstruct the right value using the values leaked around the center bin. see [here](#matched-flitering).

**If we do the FT and just extract the value around the harmonics, noise which is not periodic will be filtered, and S/N will be boosted. This also equivalents to a matched filtering in the Fourier domain.**

### Template
We want to reconstruct the profile, to which we must **break down the signal into base functions**:
$Signal(t) = \Sigma_n A_nf_n(t)$
Fourier Transform is **linear** which means $\hat F (af_1(t)+bf_2(t)) = a \hat F(f_1(t)) + b \hat F(f_2(t))$ 
So:
$$
\hat F(Signal(t)) = \hat F(\Sigma_n A_n f_n) = \Sigma_n A_n \hat F(f_n(t))
$$
for $f_n = rect(x)$ -> $\hat F(f_n) = sinc(x)$ which performs better for discrete signal, for which it can be seen as the sum up of all the bins- each is a small rect function- although it's because of the "leak" of fourier bin. More details in [reference paper](https://doi.org/10.1086/342285) at page 7. Things are different for chirp signal, which the frequency changes linearly.
$$
\Phi(t) = f_0 t + \frac{1}{2}f_1t^2 + \phi
$$
and $sinc(x)$ is not a good template anymore. Instead, we uses 
$$
f(t) = cos(f_0t+\frac{1}{2}f_1t^2 + \phi)
$$
as the function base. So we need to know what the $\hat F(f_{f0,f1}(t))$ is. In [reference paper](https://doi.org/10.1086/342285) page 15, it is written in the form of fourier bin index r ($r_o$ is the center bin, $\dot r$ is the derivative of r) and time bin index u: 
$$
n(u) = a * cos[2\pi(r_o u + \frac{\dot r} 2 u^2)+\phi] = \frac{a}{2}[e^{2\pi i(r_ou+\frac{\dot r}{2}u^2)}e^{i\phi}+e^{-2\pi i(r_ou+\frac{\dot r}{2}u^2)}e^{-i\phi}],~~~(34b)
$$
After some carefully math calculation we can get $\hat F(n(u))$ ：
$$
A_{r_c'} = \frac{aN}{2}e^{i\phi}\int_0^1e^{i\pi (\dot ru^2+2q_ru)}du,~~~(35)
$$
But there's a little mistake in the $formula (36)$ and $(38), (39)$ which should be
$$
\int_0^1e^{i\pi(\dot r u^2 + 2q_r u)}du = \frac1 {\sqrt{2\dot r}}e^{-i\pi \frac{q_r^2}{\dot r}}([C(Y_r)-C(Z_r)]+i[S(Y_r)-S(Z_r)]),~~~(36)
$$
So the correct template formula should be
$$
A_{r_c,\dot r} = \Sigma_{k = [r]-m/2}^{k = [r]+m/2} A_k\frac 1 {\sqrt{2\dot r}}e^{i\pi \frac{q_k^2}{\dot r}}([C(Y_r)-C(Z_r)]-i[S(Y_r)-S(Z_r)]),~~~(39)
$$
**There is a addition multiple of i in the [reference paper](https://doi.org/10.1086/342285)**.

## Code implementation

### signal generator

It's a pity that the writer didn't run his code on a real pulsar signal. Instead he writes a fake signal generator to generate phi and calculate the intensity of signal.

### FDAS minimum code
the minimum code of FDAS reconstruct profile.
process include:
- generate

'''python

'''


## Appendix

### Rotation caused by time shift
Let's assume some signal become a rect-like function in fourier domain. We simplify the problem by consider only 3 bin,
$(\omega-d\omega)*t, (\omega*t), (\omega+d\omega)*t$
If we shift the time to t+dt,
$(\omega-d\omega)*(t+dt), (\omega)*(t+dt), (\omega+d\omega)*(t+dt)$
There is not only a common $\delta \phi = \omega * dt$, but also a small rotation phase $(-d\omega * dt), 0, (d\omega * dt)$. It's a second order term.

### Matched filtering
If we have a time series $f(t)$ which is normalized, and we have a template $F(t)$ which is also normalized (They are all real), the convolution

$$
C(t^*) = \int_{-inf}^{inf}f(t)F(t-t^*)dt
$$

Only reaches **1** when the "shape" of template matches the signal and time delay t* = 0. If we do this in Fourier domain, the convolution will have a phase but its **magnitude** still reaches **1**. In this code we did not consider the convolution- we just reconstruct the profile using given frequency and f_dot. However, the code is same to the convolution when we swipes through different $t^*$. 
It's called matched filtering because it only gives the best result when the signal and template matches.
It's same when there is a "leak" or something caused by the discrete FT. It still produce a freq-series signal, and we need to do the inner product with Template (see [here](#template)). This would allow us to **collect** the values in each Fourier bin and produce the real value of the center frequency.
