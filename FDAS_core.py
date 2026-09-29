# packed function for reconstruct the profile
import numpy as np
from scipy.special import fresnel

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
    
    return np.sinc(q_k) #* np.exp(-1j * np.pi * q_k)

def _chirp_recon(spectrum, df, f0, f1, n_harm, ker_half_len, phase_bin_num, *, debugarray = False, center=False):

    # reconstruct the profile using given f and f shift, extract signal in spectrum
    
    # variables:
    # spectrum: signal in fourier domain
    # df: step length of spectrum (/s), also as 1/(signal duration)
    # f0: start frequency when the signal begins (Hz)
    # f1: derivative of the frequency, f_dot (Hz/s)
    # n_harm: harmonic numbers need to calculate
    # ker_half_len: half length of the template window (bin)
    
    phase_bin_num = max(n_harm * 2, phase_bin_num)# length of the profile
    harmonics_array = np.zeros(n_harm + 1, dtype=complex)
    freq_c = f0 + 0.5 * f1 / df # center freqency 
    bin_delta = f1 / (df**2) # bin_shift, also equals to f_dot*time**2
    offsets = np.arange(-ker_half_len, ker_half_len + 1)
    arrow = np.sign(f1)

    # to fix aliasing, we need variables below
    length = len(spectrum)
    N = 2 * (length - 1)
    m_orders = np.arange(1, n_harm + 1)
    exact_bins = (freq_c * m_orders) / df

    K_full = exact_bins % N
    is_conj = K_full > N / 2
    center_bins = np.where(is_conj, N - K_full, K_full)
    center_bin = np.round(center_bins).astype(int)

    # the debug variable
    # diary = np.zeros((n_harm, 2, 2*ker_half_len+1), dtype=complex)
    
    # the main loop to extract the harmonics
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
        # calculate the template
        n_bin_delta = (n+1) * bin_delta
        if np.abs(n_bin_delta) < 1e-12:
            template = __conj_template_sinc(sub_bin_q_k)
            # arrow = 1
        else:
            template = __conj_template_chirp(sub_bin_q_k,n_bin_delta)
        # correct the phase
        phase_correction = np.exp(arrow*1j*np.pi*k_indices)
        # phase_correction = (-1.0)**k_indices
        
        spectrum_parts *= phase_correction
        # get the harmonics, harmonics_array[0] set to 0 + 0j (No DC)
        harmonics_array[n+1] = np.dot(spectrum_parts, template)
        # diary[n,0,:] = spectrum_parts
        # diary[n,1,:] = template

    # (chosiable) roughly correct the total phase
    if center:
        theta_1 = np.angle(harmonics_array[1])
        countfactor = np.arange(len(harmonics_array))
        phase_rotator =(-1.0) ** countfactor * np.exp(-1j * countfactor * theta_1)
        harmonics_array *= phase_rotator

    if debugarray:
        return harmonics_array
        
    # reconstruct the profile
    profile_time = np.fft.irfft(harmonics_array, phase_bin_num)

    return profile_time