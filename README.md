# Fourier-Domain Reconstruction of Chirp Signals via Matched Filtering
A small piece of code written by Charlie Li in his summer research, as a small part of code in **Fourier Domain Acceleration Search (FDAS)** method to search for pulsars. It uses information in fourier domain to do matched filtering, and it can boost the signal to noise ratio significantly. This code is written as a *science research training* and may not reach the standard of a professional code. The writer overcame a lot of difficulties as he pushed forward. So he writes all his *confusion* and *solution* below, in which may be helpful to those who want to learn and write a similar code. He hopes it's valuable for beginners like him. *The code contains about 3 weeks of work, and it was done under the guidence of Dr. Ewan Barr, group leader of Electronics, backend development in max-planck institute for radio astronomy.*

## Directory
[General View](#general-view)
[APN reference](#apn-reference)
[Math method](#the-math-algorithm)
[Code implementation](#code-implementation)
[Appendix](appendix)

## General View
Basically the method is doing such thing:
- do FFT to the time domain signal, get the spectrum using `np.fft.rfft`
- locate the fourier bin of the given frequency (f0) and (maybe) frequency derivative (f1)
- calculate the template of funcion base, and do inner product with the spectrum to get the harmonics
- using a inverse-fourier-transform to reconstruct the profile in time domain using `np.fft.irfft`
- calculate the S/N, which should be boosted
- finding the best S/N to determine the best f0 and f1.

## APN reference

**Mostly used functions are below.**

### `_chirp_recon(spectrum, df, f0, f1, n_harm, ker_half_len, phase_bin_num, *, debugarray = False, center=False)`
**using `from reconstruct_core import _chirp_recon` to import.**
**function to reconstruct the profile using given frequency and its first order derivative.** `spectrum` is the spectrum we get using rfft. `df` is the step length of spectrum in **Hz**. `f0` is the frequency at the begining of the signal in **Hz**. `f1` is the derivative of frequency in **Hz/s**.`n_harm` is the number of harmonics we want to calculate.`ker_half_len` is the half length of the kernel, in bins. The real length of the kernel is $2*ker\_half\_len+1$. `phase_bin_num` is the size of reconstructed profile, and it will be set as the maximum of the parameter sent in and the twice of the harmonic numbers. `debugarray` is chosible, set to `True` will make the function return the value of harmonics instead of the profile. `center` is chosible, set to `True` will automatically set the signal to the center of the window, which will break the phase consistency between different parts- **The reference point of the profile phase is at the center of the signal.**
**returns `profile` in length of the biggest of `2*n_harms` and `phase_bin_num`.**

### `_mad_snr_estimator(profile)`
**using`from myfuncs import _mad_snr_estimator` to import.**
**function to estimate the profile S/N**. `profile` is the 1-D array of profile.
**returns an 1D array of normalized profile at the same size.**

### `_signal_generate_jerk(tsamp, duration, f0, f1, f2, duty_cycle, begin_phase, mode)`
**using`from myfuncs import _signal_generate_jerk` to import.**
**function to generate a fake jerk signal.**`tsamp` is the duration of the sample time in **s**. `duration` is the duration of the signal, in **s**.`f0` is the begining of frequency in **Hz**. `f1` is the first-order-derivative of frequency in **Hz/s**. `f2` is the second-order derivative of frequency in **$Hz/s^2$**. `duty_cycle` is the signal duty cycle. `begin_phase` is the begin phase of signal. `mode` can be `mode=0`or`mode=1`.`mode=0` generate a signal of $sin(ft)+0.5sin(2ft)+0.25sin(3ft)$ as a test signal. `mode=1` generate a signal of gaussian with a amplitude S/N of 1.
## The math algorithm
### Fourier Transform (FT) 
Fourier Transform is the math method of doing inner product of signal and periodic base function to gain the value of a Harmonic component.

$$
Spectrum(\omega) = Normfactor * \int_{-\infty}^{\infty} Signal(t^{\*})*e^{i\omega t^{\*}}dt^{\*}
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

### The whole process
- generate the signal in time domain- '_generate_signal'
- doing FFT using 'np.fft.rfft()' to get spectrum
- *extract harmonics in spectrum using different template*
- **(correct the phase in spectrum and template)**
- **(fixing ailasing of high-order-harmonics)**
- reconstruct profile using spectrum

### Modules
#### median S/N estimator

If we want to estimate the S/N of a signal, one method is using the **median**. The method is:
- calculate the difference of the original signal and its median $\delta x$.
- if the noise corresponds to gaussian/normal distribution, absolute value of residual $|\delta x|$ should be in a distribution.
- using median of $|\delta x|$ to calculate the standard deviation $\sigma$ of the noise. In normal distribution $\sigma = 1.4826*median(|\delta x|)$
- normalize the whole profile by $\sigma$ to get the S/N.

**This method is at the limit of small S/N. It works when the signal is not too significant to change the median of the absolute residual.** It will return a value for 'pure' signal without noise, but the value makes no sence and has no meaning. 

```python
import numpy as np

def _mad_snr_estimator(profile):
    
    median = np.median(profile)
    prof_zero = profile - median # remove the baseline
    mad = np.median(np.abs(prof_zero))
    sigma = 1.4826 * mad # gaussian
    if sigma < 1e-6: # protect branch
        return profile, 0.0
    prof_norm = prof_zero / sigma
    peak_snr = np.max(prof_norm)
  
    return prof_norm, peak_snr
```

#### signal generator

It's a pity that the writer didn't run his code on a real pulsar signal. Instead he writes a fake signal generator to generate phi and calculate the intensity of signal.
Process include:
- genterate the time sequence with equal-length sampling time dt.
- calculate the true phase using given $f0, f1, f2$ in formula $\phi(t) = f_0t+\frac 1 2 f_1t^2+\frac 1 6 f_2t^3$
- using $\phi$ to generate the signal.

```python
import numpy as np

def _signal_generate_chirp(tsamp, duration, f0_true, f1_dot, duty_cycle, begin_phase, mode):

    t = np.arange(0, duration, tsamp)
    
    signal_phase = __generate_phase_f_fdot_fd2(t, f0_true, f1_dot, 0, begin_phase)

    if mode == 0:
        # this branch is used to generate the debug fake signal to correct the phase. f(x) = sinx + 0.5sin2x + 0.25sin3x. 
        pulse = __sin_wave(signal_phase) + 0.5*__sin_wave2(signal_phase) +0.25*__sin_wave3(signal_phase)
        signal_withnoise = pulse
    elif mode == 1:
        pulse = __gaussian_wave(signal_phase, 0, duty_cycle)
        signal_withnoise = __noise(t) + pulse
        
    return signal_withnoise
```

#### FDAS search code
the minimum code of FDAS reconstruct profile.
process include:
- generate fake signal and get the spectrum
- extract **harmonic values** using f0 candidate
- reconstruct the profile
- calculate the profile S/N
- **the real f0 should corresponds to the f0 that generate the biggest S/N**

the center of the code is the extract value and reconstruct profile part.

The main loop is below:
```python
...
for idx in range(search_steps):
    prof_raw = _chirp_recon(spectrum, df, f0_candidate=trial_freqs[idx], 0, n_harm, ker_half_len, phase_bin_num)
    prof_norm, snr = _mad_snr_estimator(prof_raw)
    if snr > best_snr:
        best_snr = snr
        best_freq = trial_freqs[idx]
        best_profile = prof_norm
    profile_map[idx,:] = prof_norm
    snr_map[idx] = snr
...
```
**If we want to search $\dot f$ we simply add one more loop to search the different axis.**

### Extract the harmnonic values

for each given parameter, we would know where the center frequency of the signal is- it corresponds to the **center** of the template.
there would be differernt shapes of template. Such as sinc function, and the template we calculate for chirp signal.
when f1 is small enough, the chirp template should be degenerated to the sinc template. so we use sinc as a default when it's not possible for calculating a chirp template.
```python
def __conj_template_chirp(q_k,i_bin_delta):
    
    # using cos(f*t+0.5*f_dot*t**2) in fourier domain to extract signal
    # the funcion calculate the conj of the template
    
    # varibles:
    # q_k: r_k - r_c, diff from the center bin (sub_bin)
    # i_bin_delta: bin shift in the fourier domain, '*i' for the i-th harm
    
    abs_delta = np.abs(i_bin_delta)
    sgn = np.sign(i_bin_delta)
    norm_factor = (2*abs_delta) ** -0.5
    Y_k = 2 * norm_factor * (q_k - 0.5*abs_delta)
    Z_k = 2 * norm_factor * (q_k + 0.5*abs_delta)
    S_Y, C_Y = fresnel(Y_k)
    S_Z, C_Z = fresnel(Z_k)
    integral_factor = (C_Z -C_Y) - 1j*sgn*(S_Z-S_Y)
    phase_factor = np.exp(sgn*1j*np.pi*(q_k**2)/abs_delta)
    template = norm_factor*phase_factor*integral_factor
    
    return template

def __conj_template_sinc(q_k):
    
    # protection, for f1_shift < 1e-12
    # generate the sinc(x) template- also conj
    
    return np.sinc(q_k) )
```
And next process is to calculate the spectrum part which is used to do inner product with template. there will meet with [ailase](#ailse-correction). The code is in that part, too.
Remember to use `np.round()` to calculate the index of centerbin, to better deal with the sub-bin offsets.

### Template phase correction
**If the spectrum part and the template are not in single bin (they have a length), must take the phase correction into consideration. Especially for chirp signal.**
**Complex in math understanding but not in code.**
if we want to do the inner product of the spectrum and the template, we must make sure the "shape" and the "time" is all correct. In appendix [here](#rotation-caused-by-time-shift), we can see the time shift will lead to a rotation of all the values of spectrum in complex plane.
**Since the center of the template is at the middle of chirp signal- it's also the phase reference point- if we want to do the matched filtering, we also need to rotate the signal to the center of the time-domain, to make the reference point of the signal at the middle. -That leads to a rotation in Fourier bins.**

The main process is as below:
- get the center bin index
- get the values of the spectrum around the center bin
- calculate the offsets of the bins around it
- calculate the template
- correct the phase by multiple a rotation
- do inner product to extract the value of harmonics.

```python
...
    for n in range(n_harm):
        ...
        if np.abs(n_bin_delta) < 1e-12:
            # protect branch
            template = __conj_template_sinc(sub_bin_q_k)
            arrow = 1
        else:
            template = __conj_template_chirp(sub_bin_q_k,n_bin_delta)
        # correct the phase
        phase_correction = (-1.0)**k_indices
        # it is the same as the code below:
        # phase_correction = np.exp(arrow*1j*np.pi*k_indices)
        
        spectrum_parts *= phase_correction
        # get the harmonics, harmonics_array[0] set to 0 + 0j (No DC)
        harmonics_array[n+1] = np.dot(spectrum_parts, template)
        ...
...
```

### Ailase correction

Undersampling will lead to ailasing in spectrum. That is because the biggest frequency which is determined by the sample dt $f_{max} = \frac {1} {\delta t}$ is not big enough to collect all the harmonics. That will lead to the high frequency spectrum to filp back into the low frequency spectrum. But the information is not really lost. We can extract those information by simply mirror the values of the Fourier bins.
For the math demonstration of the aliasing, click [here](#ailase-in-math). Solve ailasing is easy in math, but a little bit longer in code.
the method we use is below:
- for a percular frequency, calculate the index of the center bin
- get the length of the spectrum
- calculate the start and end index of the spectrum part
- calculate if there is a aliase in the spectrum parts. If so, slice, filp, and reconstruct the spectrum part.
- do inner product.
  
```python
for n in range(n_harm):
        spectrum_parts = np.zeros(2 * ker_half_len + 1, dtype=complex)
        begin_idx = center_bin[n] - ker_half_len
        end_idx = center_bin[n] + ker_half_len + 1
        # slice and correct the edge
        if begin_idx < 0:
            fold_length = abs(begin_idx)
            spectrum_parts[:fold_length] = np.conj(spectrum[1:fold_length + 1][::-1])
            spectrum_parts[fold_length:] = spectrum[:end_idx]
        elif end_idx > length:
            fold_length = end_idx - length
            spectrum_parts[:length - begin_idx] = spectrum[begin_idx:]
            spectrum_parts[-fold_length:] = np.conj(spectrum[-fold_length - 1:-1][::-1])
        else:
            spectrum_parts[:] = spectrum[begin_idx:end_idx]
        # for higher order aliasing
        k_indices = center_bin[n] + offsets
        if is_conj[n]:
            spectrum_parts = np.conj(spectrum_parts[::-1])
            sub_bin_q_k = (N - K_full[n]) - k_indices
        else:
            sub_bin_q_k = K_full[n] - k_indices
```

## Appendix

### Rotation caused by time shift
Let's assume some signal become a rect-like function in fourier domain. We simplify the problem by consider only 3 bin,

$(\omega-d\omega) \ast t$, $(\omega \ast t)$, $(\omega+d\omega) \ast t$

If we shift the time to t+dt,

$(\omega-d\omega) \ast (t+dt)$, $(\omega) \ast (t+dt)$, $(\omega+d\omega) \ast (t+dt)$

There is not only a common $\delta \phi = \omega * dt$, but also a small rotation phase $(-d\omega \ast dt)$, 0, $(d\omega \ast dt)$. It's a second order term.

### Matched filtering
If we have a time series $f(t)$ which is normalized, and we have a template $F(t)$ which is also normalized (They are all real), the convolution

$$
C(t^{\*}) = \int_{-\infty}^{\infty}f(t)F(t-t^*)dt
$$

Only reaches **1** when the "shape" of template matches the signal and time delay t* = 0. If we do this in Fourier domain, the convolution will have a phase but its **magnitude** still reaches **1**. In this code we did not consider the convolution- we just reconstruct the profile using given frequency and f_dot. However, the code is same to the convolution when we swipes through different $t^*$. 
It's called matched filtering because it only gives the best result when the signal and template matches.
It's same when there is a "leak" or something caused by the discrete FT. It still produce a freq-series signal, and we need to do the inner product with Template (see [here](#template)). This would allow us to **collect** the values in each Fourier bin and produce the real value of the center frequency.

### Ailase in math

When we see the cars on the road, there's a illusion that all the wheels of the car seems to be rotate **backwards**. This is caused by undersampling.
input the signal and do `np.fft.rfft` will retrun a spectrum range from $0->f_s/2$. The spectrum actually reaches $f_s$, but because of the fft it's the same for $0->f_s/2$ and $f_S/2->f_s$ by a simple mirror. when the frequency is higher and phase changes much more than $2 \pi$ in the relate two time bin ($f>f_s$), it seems equal to a signal which is at lower freqency $f-f_S$. 
If $f = f_s/2 + \delta$, because of the symmetry- it is equal to $f_s/2 - delta$. If $f = f_s + \delta$, it equals to $\delta$- which can bee seen as a signal which is filped twice.
In general, we could see ailase as signal filping around $0$ and $f_s/2$.
