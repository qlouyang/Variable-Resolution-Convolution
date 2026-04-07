# Fast Variable Resolution Convolution for Spectra

This repository provides a highly efficient Python function to convolve stellar or exoplanetary spectra with a wavelength-dependent instrumental profile (resolving power $R(\lambda)$).

## Method

Standard convolution with a variable resolution requires computing a unique Gaussian kernel for each wavelength bin, which is computationally expensive. This code utilizes an elegant variable transformation. Assuming the resolving power follows a linear relation $R(\lambda) = a\lambda + b$, we map the wavelength space into a uniform $u$-space where the Full Width at Half Maximum (FWHM) becomes strictly 1:

$$du = \frac{R(\lambda)}{\lambda} d\lambda = \left(a + \frac{b}{\lambda}\right) d\lambda$$
$$u(\lambda) = a\lambda + b\ln(\lambda)$$

By interpolating the spectrum onto this uniform $u$-grid, we can apply an ultrafast FFT-based Gaussian broadening (`PyAstronomy.pyasl.broadGaussFast`), before interpolating back to the original wavelength grid.

## Dependencies

- `numpy`
- `scipy`
- `PyAstronomy`

## Usage

```python
import numpy as np
from convolve_spectrum_vari_R import convolve_spectrum_vari

# Example: Broaden a spectrum where R(λ) = 2.5 * λ + 1500
wavelengths = np.linspace(10000, 20000, 5000) # Angstroms
flux = np.random.normal(1.0, 0.1, 5000)       # Normalized flux
coeffs = (2.5, 1500)                          # a, b for R(λ) = a*λ + b

convolved_flux = convolve_spectrum_vari(wavelengths, flux, coeffs, quiet=False)

```
## Citation
The details can be found at **this link** (after the paper is accepted)
If you use this code in your research, please cite our paper: **this link**
