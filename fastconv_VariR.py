import numpy as np
from scipy.interpolate import interp1d
from PyAstronomy.pyasl import broadGaussFast

def fastconv_VariR(model_wave: np.ndarray, model_flux: np.ndarray, func_coeffs: tuple, quiet: bool = True) -> np.ndarray:
    """
    Convolve a spectrum with a wavelength-dependent instrumental resolution R(λ), such as JWST instruments.
    
    This function assumes a linear relationship for the resolving power: 
    R(λ) = a * λ + b. 
    It maps the wavelength grid to a uniform u-space to utilize fast Fourier transform 
    Gaussian broadening.
    
    Parameters
    ----------
    model_wave : np.ndarray
        Array of input wavelengths. Must be strictly increasing.
    model_flux : np.ndarray
        Array of input fluxes corresponding to model_wave.
    func_coeffs : tuple or list
        Coefficients (a, b) for the resolution function R(λ) = a * λ + b.
    quiet : bool, optional
        If True, suppresses print statements. Default is True.
        
    Returns
    -------
    convolved_flux : np.ndarray
        The broadened flux array evaluated on the original wavelength grid.
    """
    if not quiet:
        print('Start convolution...')

    a, b = func_coeffs
    
    # Transform to uniform u-space where FWHM is constant (=1)
    u_space = a * model_wave + b * np.log(model_wave)
    function_org = interp1d(u_space, model_flux)
    
    # Create a uniformly spaced grid with 3x the original resolution
    new_u = np.linspace(u_space[0], u_space[-1], len(u_space) * 3)
    new_y = function_org(new_u)
    
    # In u-space, FWHM is defined as 1. FWHM = 2.355 * sigma
    sigma_u_space = 1.0 / 2.355
    new_y_LR = broadGaussFast(new_u, new_y, sigma_u_space,
                              edgeHandling="firstlast", maxsig=5.0)
    
    function_new = interp1d(new_u, new_y_LR, kind='linear', bounds_error=False, fill_value="extrapolate")
    convolved_flux = function_new(u_space) # Back to the original wave
    
    if not quiet:
        print('End convolution...')
        
    return convolved_flux
