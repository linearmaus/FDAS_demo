import numpy as np

def __generate_phase_f_fdot_fd2(time, freq0, f_dot, fd2, psi_0):
    phase = (freq0 * time + 0.5 * f_dot * (time ** 2) + 1.0/6.0 * fd2 * (time ** 3) + psi_0) % 1.0
    return phase

def __gaussian_wave(phase, psi, duty):
    center_phase = psi + duty / 2.0
    phase_diff = np.abs(phase - center_phase)
    sigma = duty / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    phase_diff = np.minimum(phase_diff, 1.0 - phase_diff) 
    signal = np.exp(-0.5 * (phase_diff / sigma) ** 2)
    return signal

def __sin_wave(phase):
    return np.sin(2*np.pi*phase)

def __sin_wave2(phase):
    return np.sin(4*np.pi*phase)

def __sin_wave3(phase):
    return np.sin(6*np.pi*phase)

def __noise(time):
    np.random.seed(42)
    noise = np.random.normal(0.0, 1.0, len(time))
    return noise


def _signal_generate_chirp(tsamp, duration, f0_true, f1_dot, duty_cycle, begin_phase, mode):

    t = np.arange(0, duration, tsamp)
    
    signal_phase = __generate_phase_f_fdot_fd2(t, f0_true, f1_dot, 0, begin_phase)

    if mode == 0:
        pulse = __sin_wave(signal_phase) + 0.5*__sin_wave2(signal_phase) +0.25*__sin_wave3(signal_phase)
        signal_withnoise = pulse
    elif mode == 1:
        pulse = __gaussian_wave(signal_phase, 0, duty_cycle)
        signal_withnoise = __noise(t) + pulse
        
    return signal_withnoise

def _signal_generate_jerk(tsamp, duration, f0, f1, f2, duty_cycle, begin_phase, mode):
    t = np.arange(0, duration, tsamp)

    signal_phase = __generate_phase_f_fdot_fd2(t, f0, f1, f2, begin_phase)

    if mode == 0:
        pulse = __sin_wave(signal_phase) + 0.5*__sin_wave2(signal_phase) +0.25*__sin_wave3(signal_phase)
        signal_withnoise = pulse
    elif mode == 1:
        pulse = __gaussian_wave(signal_phase, 0, duty_cycle)
        signal_withnoise = __noise(t) + pulse
        
    return signal_withnoise

def _mad_snr_estimator(profile):
    
    median = np.median(profile)
    prof_zero = profile - median

    mad = np.median(np.abs(prof_zero))
    sigma = 1.4826 * mad # gaussian
    if sigma < 1e-6:
        return profile, 0.0

    prof_norm = prof_zero / sigma
    peak_snr = np.max(prof_norm)
    return prof_norm , peak_snr
