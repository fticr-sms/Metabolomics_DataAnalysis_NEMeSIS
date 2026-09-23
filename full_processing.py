import numpy as np
import pandas as pd 
import scipy as sp
from scipy.fft import fft
from scipy.signal import windows as _w
import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import pyopenms as oms

# Reading Bruker .d files

def get_data(path_dir1):
    """
    Takes the path to a Bruker raw data .d file and extracts the FID and processing parameters
    """

    method = [f for f in os.listdir(path_dir1) if f.endswith('.m')]
    if len(method) != 1:
        raise Exception("Path to method file is ambiguous")
    else:
        path_dir2 = os.path.join(path_dir1, method[0])
        m_file = [f for f in os.listdir(path_dir2) if f.endswith('.method')]
        if len(m_file) != 1:
            raise Exception("Path to method file is ambiguous")
        else:
            path_to_method = os.path.join(path_dir2, m_file[0])
            #print("Path to method:", path_to_method)
            tree = ET.parse(path_to_method)
            root = tree.getroot()
            #paramlist = root.find('paramlist')
            #print([t.attrib['name'] for t in list(paramlist)])
            ML1 = np.float64(root.find("paramlist/param[@name='ML1']/value").text)
            ML2 = np.float64(root.find("paramlist/param[@name='ML2']/value").text)
            ML3 = np.float64(root.find("paramlist/param[@name='ML3']/value").text)
            SW_h = np.float64(root.find("paramlist/param[@name='SW_h']/value").text)

            # Reading FID data
            fid = np.fromfile(os.path.join(path_dir1, 'fid'), dtype=np.int32)
            #print(fid)
            #print('FID length:', len(fid))

            data = {'FID': fid, 'ML1': ML1, 'ML2': ML2, 'ML3': ML3, 'SW_h': SW_h}

            return data

# Apodization functions

_ALIASES = {
"boxcar": "none", "rectangular": "none",
"hanning": "hann", "sine-squared": "hann", "sine-bell-squared": "hann",
"kilgour": "hann", "asymmetric-hann": "hann",
"sine-bell": "sine", "half-sine": "sine", "full-sine": "sine",
"blackmanharris": "blackman-harris",
"triangular": "bartlett",
    }

WINDOWS = (
"none",             # aliases: boxcar, rectangular
"hann",             # aliases: hanning, sine-squared, kilgour (tunable hann)
"hamming",
"blackman",
"blackman-harris",
"kaiser",           # tunable via beta
"sine",             # aliases: sine-bell, half-sine. at F = 0 it becomes cosine
"bartlett",
"welch",
"tukey",            # tunable via tukey_alpha
"gaussian",         # tunable via std_frac
)


def apodize(fid: np.ndarray, function: str = "sine", F: float = 0.5, beta: float = np.pi, tukey_alpha: float = 0.1, std_frac: float = 0.15):
    """Apply one of several possible apodization functions to the transient.
       - fid: The transient data to be apodized.
       - function: The apodization function to apply. 'none', 'hann', 'hamming', 
                   'blackman', 'blackman-harris', 'kaiser', 'sine', 'bartlett', 
                   'welch', 'tukey', or 'gaussian'. Aliases are also accepted.
       - F: The position of the window maximum as a fraction of transient length.
            The name is taken from kilgour and van orden, 2015. but it applies to all functions. 
            F = 0.5 is symmetric, F = 0.0 is falling half only.
       - beta: The shape parameter for the kaiser apodization function.
               0 = boxcar, ~4 ~ sine bell, ~6 ~ Hann-width but 12 dB lower sidelobes, ~8.6 ~ Blackman, 14 = very low sidelobes. 
               Values under 4 usually preferred for untargeted metabolomics.
       - tukey_alpha: Fraction of the window that is tapered in tukey apodization (0 = boxcar, 1 = Hann).
       - std_frac: Parameter for gaussian apodization, the standard deviation (sigma) as a fraction of the window length. 
    """
    # Determine the length of the input transient and normalize the function name
    n = len(fid)
    f = function.strip().lower().replace("_", "-")
    f = _ALIASES.get(f, f)


    # See if parameters are valid, and raise an error if not.
    if f not in WINDOWS:
        raise ValueError(
            f"Unknown apodization function {function!r}. "
            f"Valid options: {', '.join(WINDOWS)} "
            f"(plus aliases {', '.join(sorted(_ALIASES))})"
        )
    if not 0.0 <= F <= 0.5:
        raise ValueError("The value of F must be between 0.0 and 0.5")

    # Base shapes: each returns a symmetric window of length m, indices 0..m-1
    def _bases(beta: float, tukey_alpha: float, std_frac: float):
        return {
            "none": lambda m: np.ones(m),
            "hann": np.hanning,
            "hamming": np.hamming,
            "blackman": np.blackman,
            "blackman-harris": lambda m: _w.blackmanharris(m, sym=True),
            "kaiser": lambda m: np.kaiser(m, beta),
            "sine": lambda m: np.sin(np.pi * np.arange(m) / (m - 1)) if m > 1 else np.ones(1),
            "bartlett": np.bartlett,
            "welch": lambda m: 1.0 - ((np.arange(m) - (m - 1) / 2) / ((m - 1) / 2)) ** 2
                            if m > 1 else np.ones(1),
            "tukey": lambda m: _w.tukey(m, alpha=tukey_alpha, sym=True),
            "gaussian": lambda m: _w.gaussian(m, std=std_frac * m, sym=True),
        }
    
    # Place the maximum of base at index round(F*(n-1)), keeping both ends zero
    def _splice(base, n: int, F: float) -> np.ndarray:
        if F == 0.5:
            return base(n)                      
        j = int(round(F * (n - 1)))             # index of the maximum
        fall = base(2 * (n - 1 - j) + 1)[n - 1 - j:]
        if j == 0:
            return fall
        return np.concatenate([base(2 * j + 1)[:j], fall])

    window = _splice(_bases(beta, tukey_alpha, std_frac)[f], n, F)

    return fid * window

# Intensities and frequencies calculation

def fid_to_ints(fid_apo: np.ndarray, zp_factor: int=15):
    """Compute the FFT of the (apodized) transient, apply zero-padding and extract relevant intensity information.
       - fid: The transient data.
       - zp_factor: The factor by which to zero-pad the transient before FFT.
       Padding acts as a natural interpolation that makes peak picking possible in a complex spectrum.
       The initial product of the fft is a symmetrical spectrum, with identical peaks to both sides 
       representing the same chemical entities. The peak intensities are the absolute values of each peak.
    """
    fid_apo_fft = fft(fid_apo, len(fid_apo) * zp_factor)
    ints = np.abs(fid_apo_fft[:int(len(fid_apo_fft)/2)])

    ints = ints[::-1] # Inverting the array to be consistent with masses and frequencies

    return ints

def calculate_frequencies(ints: np.ndarray, sw_h: float):
    """
    Calculate frequency axis
    
    ints: array of intensities (output of fid_to_ints). 
    sw_h: spectral width in Hz
    """
    n_points = len(ints)
    # Frequency spacing
    df = sw_h / n_points
    
    # Frequency array 
    frequencies = np.arange(n_points) * df

    frequencies = frequencies[::-1] # Inverting the array so that higher frequencies (smaller m/z) appear first
    
    return frequencies


# Mass calibration (spectra may still require re-calibration with one or more known peaks)

def mz_from_frequency(frequencies, ML1: float, ML2: float, ML3: float = 0.0):
    """
    Apply the Bruker FT-ICR calibration equation

        m/z = ML1 / (f + ML2) + ML3 / (f + ML2)**2

    ML3 is zero on most files, in which case the second term vanishes.
    """
    x = np.asarray(frequencies, dtype=np.float64) + ML2
    x = np.where(x == 0, 1e-10, x)
    if ML3 != 0.0:
        return ML1 / x + ML3 / x ** 2
    return ML1 / x
 
 
def ppm_error(mz_obs, mz_ref):
    """
    Calculate mass error in ppm of the observed mass:
        error_ppm = (mz_obs - mz_ref) / mz_obs * 1e6
    """
    mz_obs = np.asarray(mz_obs, dtype=np.float64)
    return (mz_obs - np.asarray(mz_ref, dtype=np.float64)) / mz_obs * 1e6


def frequency_from_mz(mz, ML1: float, ML2: float, ML3: float = 0.0):
    """
    Inverse of the calibration equation.
    """
    mz = np.asarray(mz, dtype=np.float64)
    if ML3 == 0.0:
        x = ML1 / mz
    else:
        x = (ML1 + np.sqrt(ML1 ** 2 + 4.0 * mz * ML3)) / (2.0 * mz)
    return x - ML2


def get_raw_mzs(spectrum: dict):
    """
    Derive the profile m/z axis from the stored frequencies and the calibration
    coefficients in the spectrum dictionary.
    """
    return mz_from_frequency(spectrum['freqs'], spectrum['ML1'],
                             spectrum['ML2'], spectrum['ML3'])


# Wrapper function (takes initial .d file as input and returns mass spectrum as a dictionary) 

def raw_spectrum(path_dir: str, apodization_function: str = "sine", F: float = 0.5, beta: float = np.pi, tukey_alpha: float = 0.1, std_frac: float = 0.15, zp_factor: int = 15, mz_range: tuple = None):
    """
    Build spectrum as a dictionary from a Bruker .d file.

    - mz_range: (min_mz, max_mz) of the acquisition/analysis. 
      Later filtering is also possible with filter_spectrum()

    See other functions for the remaining parameters.
    """

    data = get_data(path_dir)

    fid = data['FID']
    ML1 = data['ML1']
    ML2 = data['ML2']
    ML3 = data['ML3']
    sw_h = data['SW_h']

    # Apodization
    fid_apo = apodize(fid, function=apodization_function, F=F, beta=beta, tukey_alpha=tukey_alpha, std_frac=std_frac)

    # Calculate intensities
    ints = fid_to_ints(fid_apo, zp_factor=zp_factor)

    # Calculate frequencies
    freqs = calculate_frequencies(ints, sw_h)

    # Convert m/z bounds to frequency bounds and determine freq to keep
    f_bounds = None
    if mz_range is not None:
        mz_low, mz_high = float(mz_range[0]), float(mz_range[1])
        if not mz_low < mz_high:
            raise ValueError("mz_range must be (min_mz, max_mz) with min_mz < max_mz")
        f_at_low_mz = frequency_from_mz(mz_low, ML1, ML2, ML3)
        f_at_high_mz = frequency_from_mz(mz_high, ML1, ML2, ML3)
        keep = (freqs >= f_at_high_mz) & (freqs <= f_at_low_mz)
        f_bounds = (float(f_at_high_mz), float(f_at_low_mz))

    if not keep.any():
        raise ValueError(f"mz_range {mz_range} leaves no points; the acquisition "
                         f"covers m/z {mz_from_frequency(sw_h, ML1, ML2, ML3):.3f} upwards")

    # Apply the bounds to both freqs and ints                     
    freqs, ints = freqs[keep], ints[keep]

    # Build spectrum dictionary (intensity, frequencies, apodization function and parameters for reference)
    raw_spectrum = {'raw_ints': ints,
                'freqs': freqs,
                'ML1': ML1,
                'ML2': ML2,
                'ML3': ML3,
                'apodization_function': apodization_function,
                'F': F,
                'beta': beta,
                'tukey_alpha': tukey_alpha,
                'std_frac': std_frac,
                'zp_factor': zp_factor,
                'mz_range': mz_range,       # as supplied, None if not trimmed
                'f_range': f_bounds,        # the same bounds in Hz, None if not trimmed
                'SW_h': sw_h}

    return raw_spectrum


### PyOpenMS centroiding

  
def pyopenms_centroiding(raw_spectrum: dict, snr: float = 3.0, mre: int = 5, spd: float = 5.0,
                         spacing_difference_gap: float = 4.0, missing: int = 1,
                         report_FWHM: bool = True, report_FWHM_unit: str = 'relative',
                         sn_max_intensity: int = -1, sn_auto_max_stdev_factor: float = 3.0,
                         sn_auto_max_percentile: int = 95, sn_auto_mode: int = 0,
                         sn_win_len: float = 200.0, sn_bin_count: int = 30,
                         sn_noise_for_empty_window: float = 1e20,
                         sn_write_log_messages: bool = True):
    """
    Centroid the spectrum using PyOpenMS.
 
    raw_spectrum: dictionary from the previous processing steps.
 
    PeakPickerHiRes parameters:
      snr: Signal-to-noise threshold for peak picking. 0 disables the filter.
      mre: 'SignalToNoise:min_required_elements', minimum number of raw points
           a noise-estimation window must hold. A window with fewer is assigned
           sn_noise_for_empty_window, which discards every peak inside it.
      spd: 'spacing_difference', maximum allowed spacing to the next raw point,
           in multiples of the minimal spacing at the apex, for the peak to be
           extended in that direction.
      spacing_difference_gap: spacing, in multiples of the minimal spacing, above
           which the flank is treated as a gap and the peak is stopped.
      missing: number of consecutive missing points tolerated on a flank.
      report_FWHM: attach the peak widths to the result ('cent_fwhm').
      report_FWHM_unit: 'relative' (ppm) or 'absolute' (Da).
      sn_*: the remaining SignalToNoise sub-parameters.
 
    Returns the spectrum dictionary with the centroid arrays and the peak
    picking parameters added. 
    """
    if report_FWHM_unit not in ('relative', 'absolute'):
        raise ValueError("report_FWHM_unit must be 'relative' or 'absolute'")
 
    # Create a PyOpenMS spectrum object.
    raw_mzs = get_raw_mzs(raw_spectrum) # Calculate masses from freqs
    oms_spectrum = oms.MSSpectrum() 
    oms_spectrum.setMSLevel(1) # This software is for MS1 data only, so this is hardcoded
    oms_spectrum.set_peaks((np.ascontiguousarray(raw_mzs, dtype=np.float64),
                            np.ascontiguousarray(raw_spectrum['raw_ints'], dtype=np.float64)))
    if not oms_spectrum.isSorted():
        raise ValueError("The m/z axis handed to PeakPickerHiRes is not ascending")
 
    # Create a peak picker object
    peak_picker = oms.PeakPickerHiRes()
    picked_spectrum = oms.MSSpectrum()
 
    # Set parameters for the peak picker
    params = peak_picker.getParameters()
    params.setValue("signal_to_noise", snr)
    params.setValue("spacing_difference", spd)
    params.setValue("spacing_difference_gap", spacing_difference_gap)
    params.setValue("missing", missing)
    params.setValue("ms_levels", [1])
    params.setValue("report_FWHM", "true" if report_FWHM else "false")
    params.setValue("report_FWHM_unit", report_FWHM_unit)
    params.setValue("SignalToNoise:max_intensity", sn_max_intensity)
    params.setValue("SignalToNoise:auto_max_stdev_factor", sn_auto_max_stdev_factor)
    params.setValue("SignalToNoise:auto_max_percentile", sn_auto_max_percentile)
    params.setValue("SignalToNoise:auto_mode", sn_auto_mode)
    params.setValue("SignalToNoise:win_len", sn_win_len)
    params.setValue("SignalToNoise:bin_count", sn_bin_count)
    params.setValue("SignalToNoise:min_required_elements", mre)
    params.setValue("SignalToNoise:noise_for_empty_window", sn_noise_for_empty_window)
    params.setValue("SignalToNoise:write_log_messages", "true" if sn_write_log_messages else "false")
    peak_picker.setParameters(params)
 
    # Perform peak picking
    peak_picker.pick(oms_spectrum, picked_spectrum)
 
    # Extract m/z and intensity from the picked spectrum
    centroids_mzs, centroids_ints = picked_spectrum.get_peaks()
 
    # Extract the peak widths. PeakPickerHiRes returns them in a float data
    # array named 'FWHM_ppm' when the unit is relative and 'FWHM' when absolute
    centroids_fwhm = None
    if report_FWHM:
        fwhm_name = 'FWHM_ppm' if report_FWHM_unit == 'relative' else 'FWHM'
        for array in picked_spectrum.getFloatDataArrays():
            name = array.getName()
            if isinstance(name, bytes):
                name = name.decode()
            if name == fwhm_name:
                centroids_fwhm = np.asarray(array.get_data(), dtype=np.float64)
                break
        if centroids_fwhm is None:
            raise RuntimeError(
                f"report_FWHM is on but PeakPickerHiRes returned no {fwhm_name!r} "
                f"array (found: {[a.getName() for a in picked_spectrum.getFloatDataArrays()]})")
        if centroids_fwhm.size != centroids_mzs.size:
            raise RuntimeError(
                f"the {fwhm_name!r} array holds {centroids_fwhm.size} values but "
                f"{centroids_mzs.size} centroids were picked")
 
    # Calculate frequencies for the centroided m/z values
    centroids_freqs = frequency_from_mz(centroids_mzs, raw_spectrum['ML1'], raw_spectrum['ML2'], raw_spectrum['ML3'])
 
    # Build centroided spectrum dictionary. Processing parameters are saved for future reference
    spectrum = dict(raw_spectrum)
    spectrum.update({
        'cent_mzs': centroids_mzs,
        'cent_ints': centroids_ints,
        'cent_freqs': centroids_freqs,
        'cent_fwhm': centroids_fwhm,
        'snr': snr,
        'mre': mre,
        'spd': spd,
        'spacing_difference_gap': spacing_difference_gap,
        'missing': missing,
        'report_FWHM': report_FWHM,
        'report_FWHM_unit': report_FWHM_unit,
        'sn_max_intensity': sn_max_intensity,
        'sn_auto_max_stdev_factor': sn_auto_max_stdev_factor,
        'sn_auto_max_percentile': sn_auto_max_percentile,
        'sn_auto_mode': sn_auto_mode,
        'sn_win_len': sn_win_len,
        'sn_bin_count': sn_bin_count,
        'sn_noise_for_empty_window': sn_noise_for_empty_window,
        'sn_write_log_messages': sn_write_log_messages
    })
 
    return spectrum


def get_cent_mzs(spectrum: dict):
    """
    Derive the centroid m/z axis from the stored frequencies and the calibration
    coefficients currently held in the spectrum dictionary.

    The m/z axis is deliberately not stored: it is a function of 'cent_freqs' and
    (ML1, ML2, ML3), so recomputing it here guarantees it always matches the
    coefficients in the dictionary, including after a re-calibration.
    """
    return mz_from_frequency(spectrum['cent_freqs'], spectrum['ML1'],
                             spectrum['ML2'], spectrum['ML3'])


_CAL_MODES = {'ML1': 1, 'ML1+ML2': 2}


def match_calibrants(spectrum: dict, mz_ref, tol_ppm: float = 3.0,
                     min_intensity: float = None, min_rel_intensity: float = 0.0):
    """
    List every centroid lying within `tol_ppm` of each reference mass.

    - spectrum : dictionary produced by pyopenms_centroiding.
    - mz_ref   : theoretical m/z of the calibrants.
    - tol_ppm  : matching tolerance. 
    - min_intensity, min_rel_intensity : optional intensity thresholds, the
                 latter as a fraction of the base peak.

    Returns a dict:
      'mz_ref'       : array of the reference masses, in the order given
      'candidates'   : list, one entry per reference mass, each a dict of arrays
                       ('index', 'mz', 'frequency', 'intensity', 'fwhm',
                       'error_ppm') sorted by |error_ppm|
      'n_candidates' : number of candidates found for each reference mass
      'unmatched'    : reference masses with no candidate at all
      'tol_ppm'      : the tolerance used
      'fwhm_unit'    : 'relative' (ppm) or 'absolute' (Da), or None if the
                       spectrum carries no peak widths
    """

    cent_freqs = spectrum['cent_freqs']
    cent_ints = spectrum['cent_ints']

    fwhm = spectrum.get('cent_fwhm')
    fwhm = fwhm if fwhm is not None else None
    fwhm_unit = spectrum.get('report_FWHM_unit') if fwhm is not None else None

    mz_axis = mz_from_frequency(cent_freqs, spectrum['ML1'], spectrum['ML2'],
                                spectrum['ML3'])

    thr = 0.0
    if min_rel_intensity:
        thr = float(min_rel_intensity) * float(cent_ints.max())
    if min_intensity is not None:
        thr = max(thr, float(min_intensity))

    candidates, n_cand, unmatched = [], [], []

    for mz_t in mz_ref:
        err = ppm_error(mz_axis, mz_t)
        hit = np.flatnonzero((np.abs(err) < tol_ppm) & (cent_ints > thr)) # see which peaks are matches
        hit = hit[np.argsort(np.abs(err[hit]))] # smallest errors first
        candidates.append({
            'index': hit,
            'mz': mz_axis[hit],
            'frequency': cent_freqs[hit],
            'intensity': cent_ints[hit],
            'fwhm': fwhm[hit] if fwhm is not None else np.full(hit.size, np.nan),
            'error_ppm': err[hit],
        })
        n_cand.append(hit.size)
        if hit.size == 0:
            unmatched.append(float(mz_t))

    return {
        'mz_ref': mz_ref,
        'candidates': candidates,
        'n_candidates': np.asarray(n_cand, dtype=int),
        'unmatched': unmatched,
        'tol_ppm': float(tol_ppm),
        'fwhm_unit': fwhm_unit,
    }


def calibrants_table(matched: dict):
    """
    Flatten the output of match_calibrants into a DataFrame, one row per
    candidate, for inspection before choosing. 'position' is the number to pass
    in `selection` to override the default (position 0, the lowest error).
    """
    rows = []
    for i, (mz_t, c) in enumerate(zip(matched['mz_ref'], matched['candidates'])):
        if c['index'].size == 0:
            rows.append({'calibrant': i, 'mz_ref': mz_t, 'position': np.nan,
                         'centroid': np.nan, 'mz_obs': np.nan, 'error_ppm': np.nan,
                         'intensity': np.nan, 'fwhm': np.nan, 'default': False})
            continue
        for pos in range(c['index'].size):
            rows.append({'calibrant': i, 'mz_ref': mz_t, 'position': pos,
                         'centroid': c['index'][pos], 'mz_obs': c['mz'][pos],
                         'error_ppm': c['error_ppm'][pos],
                         'intensity': c['intensity'][pos], 'fwhm': c['fwhm'][pos],
                         'default': pos == 0})
    return pd.DataFrame(rows)


def select_calibrants(matched: dict, selection=None):
    """
    Reduce the candidate lists to one centroid per reference mass.

    - selection : None to take the lowest-error candidate for every reference
                  mass, or a dict {calibrant_index: candidate_position}, or a
                  sequence of positions with one entry per reference mass.
                  Positions refer to the ordering from match_calibrants, so 0 is
                  the closest peak, 1 the next closest, and so on.

    Reference masses with no candidate are dropped and listed in 'unmatched'.
    """
    n_ref = matched['mz_ref'].size
    if selection is None:
        sel = {}
    elif isinstance(selection, dict):
        sel = {int(k): int(v) for k, v in selection.items()}
    else:
        seq = list(selection)
        if len(seq) != n_ref:
            raise ValueError(f"selection has {len(seq)} entries but there are "
                             f"{n_ref} reference masses")
        sel = {i: int(p) for i, p in enumerate(seq)}
    for i in sel:
        if not 0 <= i < n_ref:
            raise ValueError(f"selection refers to calibrant {i}, which does not exist")

    mz_ref, f_obs, mz_obs, err, amp, width, idx, pos_used = [], [], [], [], [], [], [], []
    unmatched = []

    for i, (mz_t, c) in enumerate(zip(matched['mz_ref'], matched['candidates'])):
        if c['index'].size == 0:
            unmatched.append(float(mz_t))
            continue
        p = sel.get(i, 0)
        if not 0 <= p < c['index'].size:
            raise ValueError(
                f"calibrant {i} (m/z {mz_t:.6f}) has {c['index'].size} candidate(s) "
                f"within {matched['tol_ppm']:g} ppm, so position {p} does not exist")
        mz_ref.append(float(mz_t)); f_obs.append(float(c['frequency'][p]))
        mz_obs.append(float(c['mz'][p])); err.append(float(c['error_ppm'][p]))
        amp.append(float(c['intensity'][p])); width.append(float(c['fwhm'][p]))
        idx.append(int(c['index'][p])); pos_used.append(p)

    return {
        'mz_ref': np.asarray(mz_ref), 'f_obs': np.asarray(f_obs),
        'mz_obs': np.asarray(mz_obs), 'error_ppm': np.asarray(err),
        'intensity': np.asarray(amp), 'fwhm': np.asarray(width),
        'index': np.asarray(idx, dtype=int),
        'position': np.asarray(pos_used, dtype=int),
        'unmatched': unmatched,
    }


def fit_calibration(f_obs, mz_ref, ML1: float, ML2: float, ML3: float = 0.0,
                    mode: str = 'auto'):
    """
    Re-calculate the calibration coefficients so that the calibrant masses
    match the expected ones by the equation:

        m/z = ML1 / (f + ML2) + ML3 / (f + ML2)**2

    - mode : 'ML1'      -> scale factor only; 1 calibrant is enough
             'ML1+ML2'  -> scale and frequency offset; needs >= 2 calibrants
             'auto'     -> 'ML1' for a single calibrant, 'ML1+ML2' otherwise
    - ML3 is never fitted. The function returns it unchanged.
    """

    # normalizing and validating the inputs
    f_obs = np.atleast_1d(np.asarray(f_obs, dtype=np.float64)) 
    mz_ref = np.atleast_1d(np.asarray(mz_ref, dtype=np.float64)) # forces 1d arrays insteasd of floats, useful for report
    if f_obs.size != mz_ref.size:
        raise ValueError("f_obs and mz_ref must have the same length")
    n = f_obs.size
    if n == 0:
        raise ValueError("No calibrants were provided")
    if not np.all(np.isfinite(f_obs)) or not np.all(np.isfinite(mz_ref)):
        raise ValueError("non-finite value among the calibrant frequencies or masses")
    if np.any(mz_ref <= 0):
        raise ValueError("reference masses must be positive")

    # determining which coefficients o fit
    if mode == 'auto':
        mode = 'ML1' if n == 1 else 'ML1+ML2'
    if mode not in _CAL_MODES:
        raise ValueError(f"Unknown mode {mode!r}. Valid options: "
                         f"{', '.join(_CAL_MODES)} (ML3 is never fitted)")
    n_par = _CAL_MODES[mode]
    if n < n_par:
        raise ValueError(f"mode {mode!r} needs at least {n_par} calibrants, got {n}")

    # calculating pre-recalibration errors
    err_before = ppm_error(mz_from_frequency(f_obs, ML1, ML2, ML3), mz_ref)

    # decomposing the equation and taking ML3 out of the way
    x_in = f_obs + ML2
    mz_lin = mz_ref - (ML3 / x_in**2 if ML3 else 0.0)

    # solving
    if mode == 'ML1':
        # ML2 is frozen, so assuming only one calibrant:
        #   ML1 = mz_ref * (f_obs + ML2)
        # If there is more than one calibrant, the optimal
        # solution will be the mean. 
        ML1_new = float(np.mean(mz_lin * x_in))
        ML2_new = float(ML2)
    else:
        # Multiplying mz_lin = ML1/(f + ML2) through by (f + ML2):
        #     mz_lin*f = ML1 - mz_lin*ML2
        #     mz_lin*f = ML1*1 + ML2*(-mz_lin) 
        # In matrix form, the equations system can be written as A*x = y
        # where x is the column of unknowns (ML1 and ML2) and A is the
        # coefficients they are multiplied for (1 for ML1 and -mz_lin for
        # ML2). y is the column of known terms (mz_lin*f)

        if np.unique(mz_ref).size < 2:
            raise ValueError("ML1+ML2 needs at least two distinct reference masses")
        
        A = np.column_stack([np.ones(n), -mz_lin])
        y = mz_lin * f_obs
        sol, *_ = np.linalg.lstsq(A, y, rcond=None)

        ML1_new, ML2_new = float(sol[0]), float(sol[1])

    # calculating post-recalibration errors
    err_after = ppm_error(mz_from_frequency(f_obs, ML1_new, ML2_new, ML3), mz_ref)

    # assembling dictionary
    return {
        'ML1': ML1_new,
        'ML2': ML2_new,
        'ML3': float(ML3),                     # carried through, never fitted
        'mode': mode,
        'n_ref': int(n),
        'mz_ref': mz_ref,
        'f_obs': f_obs,
        'error_ppm_before': err_before,
        'error_ppm_after': err_after,
        'rms_ppm_before': float(np.sqrt(np.mean(err_before ** 2))),
        'rms_ppm_after': float(np.sqrt(np.mean(err_after ** 2))),
        'max_abs_ppm_after': float(np.max(np.abs(err_after))),
    }


def apply_calibration(spectrum: dict, ML1: float, ML2: float, ML3: float = 0.0):
    """
    Return a copy of `spectrum` with the calibration coefficients replaced and
    'cent_mzs' recomputed from 'cent_freqs'. The frequency axes are the
    measurement and are left untouched, and the profile m/z axis is derived on
    demand by get_raw_mzs(), so nothing else needs updating.

    The coefficients the file was acquired with are kept under 'ML1_acq',
    'ML2_acq' and 'ML3_acq', written only the first time, so a spectrum can be
    re-calibrated repeatedly without losing them.
    """
    new = dict(spectrum)
    for key in ('ML1', 'ML2', 'ML3'):
        new.setdefault(key + '_acq', spectrum[key])
    new['ML1'], new['ML2'], new['ML3'] = float(ML1), float(ML2), float(ML3)
    if 'cent_freqs' in spectrum:
        new['cent_mzs'] = mz_from_frequency(spectrum['cent_freqs'], ML1, ML2, ML3)
    return new


def recalibrate(spectrum: dict, mz_ref, tol_ppm: float = 1.0, mode: str = 'auto',
                selection=None, min_intensity: float = None,
                min_rel_intensity: float = 0.0):
    """
    Re-calibrate a centroided spectrum against a set of exactly known masses.

    Candidates are listed by match_calibrants, one is chosen per reference mass
    by select_calibrants (lowest absolute error unless `selection` says
    otherwise), ML1 and optionally ML2 are re-fitted in the frequency domain,
    and the centroid m/z axis is recomputed. The profile axis follows
    automatically, since get_raw_mzs() derives it from 'freqs' and the
    coefficients. ML3 is passed through unchanged.

    - tol_ppm   : matching tolerance, 1 ppm by default.
    - selection : None, a dict {calibrant_index: candidate_position}, or a
                  sequence of positions, one per reference mass.
    - mode      : passed to fit_calibration.

    Every calibrant selected is used, with equal weight. Nothing is rejected or
    down-weighted; if a calibrant is unreliable, drop it from `mz_ref` or point
    `selection` at a different candidate.
    """
    matched = match_calibrants(spectrum, mz_ref, tol_ppm=tol_ppm,
                               min_intensity=min_intensity,
                               min_rel_intensity=min_rel_intensity)
    chosen = select_calibrants(matched, selection=selection)
    if chosen['mz_ref'].size == 0:
        raise ValueError(
            f"No centroid lies within {tol_ppm:g} ppm of any reference mass. "
            "Either the calibrants are absent from this spectrum or the "
            "calibration is further out than the tolerance allows.")

    ML1_0 = float(spectrum['ML1'])
    ML2_0 = float(spectrum['ML2'])
    ML3_0 = float(spectrum['ML3'])
    fit = fit_calibration(chosen['f_obs'], chosen['mz_ref'], ML1_0, ML2_0, ML3_0,
                          mode=mode)

    new = apply_calibration(spectrum, fit['ML1'], fit['ML2'], fit['ML3'])
    new['recalibration'] = {
        'mode': fit['mode'],
        'n_ref': fit['n_ref'],
        'ML1_before': ML1_0, 'ML2_before': ML2_0, 'ML3_before': ML3_0,
        'mz_ref': chosen['mz_ref'],
        'f_obs': chosen['f_obs'],
        'mz_obs_before': chosen['mz_obs'],
        'mz_obs_after': mz_from_frequency(chosen['f_obs'], fit['ML1'],
                                          fit['ML2'], fit['ML3']),
        'intensity': chosen['intensity'],
        'fwhm': chosen['fwhm'],
        'fwhm_unit': matched['fwhm_unit'],
        'centroid_index': chosen['index'],
        'position_used': chosen['position'],
        'n_candidates': matched['n_candidates'],
        'candidates': matched['candidates'],
        'error_ppm_before': fit['error_ppm_before'],
        'error_ppm_after': fit['error_ppm_after'],
        'rms_ppm_before': fit['rms_ppm_before'],
        'rms_ppm_after': fit['rms_ppm_after'],
        'max_abs_ppm_after': fit['max_abs_ppm_after'],
        'unmatched': chosen['unmatched'],
        'tol_ppm': float(tol_ppm),
        'selection': selection,
    }
    return new


def calibration_report(spectrum: dict):
    """
    Text summary of a re-calibrated spectrum, listing the calibrants used and
    every other candidate that fell inside the tolerance, with the position
    number needed to select it instead.
    """
    r = spectrum['recalibration'] if 'recalibration' in spectrum else spectrum
    unit = r.get('fwhm_unit')
    w_label = 'FWHM(ppm)' if unit == 'relative' else ('FWHM(Da)' if unit else 'FWHM')
    lines = [
        f"mode: {r['mode']}   calibrants used: {r['n_ref']}   "
        f"tolerance: {r['tol_ppm']:g} ppm",
        f"ML1: {r['ML1_before']:.12g} -> {spectrum.get('ML1', float('nan')):.12g}",
        f"ML2: {r['ML2_before']:.12g} -> {spectrum.get('ML2', float('nan')):.12g}",
        f"ML3: {r['ML3_before']:.12g} (not fitted)",
        f"RMS error: {r['rms_ppm_before']:.4f} -> {r['rms_ppm_after']:.4f} ppm"
        f"   (max |err| {r['max_abs_ppm_after']:.4f} ppm)",
    ]
    if r.get('unmatched'):
        lines.append("no candidate within tolerance: "
                     + ", ".join(f"{m:.6f}" for m in r['unmatched']))
    lines.append("")
    lines.append(f"{'m/z ref':>13} {'m/z before':>13} {'m/z after':>13} "
                 f"{'before(ppm)':>12} {'after(ppm)':>11} {'pos':>4}")
    for mzr, mb, ma, eb, ea, p in zip(r['mz_ref'], r['mz_obs_before'], r['mz_obs_after'],
                                      r['error_ppm_before'], r['error_ppm_after'],
                                      r['position_used']):
        lines.append(f"{mzr:13.6f} {mb:13.6f} {ma:13.6f} {eb:12.4f} {ea:11.4f} {p:4d}")

    extra = [(i, c) for i, c in enumerate(r.get('candidates', [])) if c['index'].size > 1]
    if extra:
        lines.append("")
        lines.append("More than one centroid inside the tolerance; pass "
                     "selection={calibrant: position} to choose another:")
        lines.append(f"  {'calibrant':>9} {'pos':>4} {'m/z':>13} {'error(ppm)':>11} "
                     f"{'intensity':>12} {w_label:>10}")
        for i, c in extra:
            for p in range(c['index'].size):
                lines.append(f"  {i:9d} {p:4d} {c['mz'][p]:13.6f} {c['error_ppm'][p]:11.4f} "
                             f"{c['intensity'][p]:12.4g} {c['fwhm'][p]:10.4g}")
    return "\n".join(lines)


### Spectrum filtering

def filter_spectrum(spectrum: dict, mz_range: tuple = None):
    """
    Filter the spectrum by m/z range.

    spectrum: dictionary from previous processing steps.
    mz_range: tuple (min_mz, max_mz) to filter the spectrum. If None, no filtering is applied.

    The profile and centroid arrays are each cut with a single mask built on
    their frequency axis, so the arrays of a given set always keep the same
    length. Every other entry of the dictionary is carried over unchanged.
    """
    if mz_range is None:
        return spectrum

    min_mz, max_mz = mz_range
    min_mz_freq = frequency_from_mz(min_mz, spectrum['ML1'], spectrum['ML2'], spectrum['ML3'])
    max_mz_freq = frequency_from_mz(max_mz, spectrum['ML1'], spectrum['ML2'], spectrum['ML3'])

    # Note: frequencies decrease as m/z increases
    freq_mask = (spectrum['freqs'] >= max_mz_freq) & (spectrum['freqs'] <= min_mz_freq)
    cent_freq_mask = (spectrum['cent_freqs'] >= max_mz_freq) & (spectrum['cent_freqs'] <= min_mz_freq)

    filtered_spectrum = dict(spectrum)
    filtered_spectrum['freqs'] = spectrum['freqs'][freq_mask]
    filtered_spectrum['raw_ints'] = spectrum['raw_ints'][freq_mask]
    filtered_spectrum['cent_mzs'] = spectrum['cent_mzs'][cent_freq_mask]
    filtered_spectrum['cent_ints'] = spectrum['cent_ints'][cent_freq_mask]
    filtered_spectrum['cent_freqs'] = spectrum['cent_freqs'][cent_freq_mask]
    if spectrum.get('cent_fwhm') is not None:
        filtered_spectrum['cent_fwhm'] = spectrum['cent_fwhm'][cent_freq_mask]

    return filtered_spectrum
