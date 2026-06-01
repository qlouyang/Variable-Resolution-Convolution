# Fast Variable Resolution Convolution for Spectra

This repository provides a highly efficient Python function to convolve stellar or exoplanetary spectra with a wavelength-dependent instrumental resolving power $$R(\lambda)$$.

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
from fastconv_VariR import fastconv_VariR

# Example: Broaden a spectrum where R(λ) = 500 * λ + 50
wavelengths = np.linspace(0.6, 3.0, 2000) # microns
flux = np.random.normal(1.0, 0.1, 5000)       # Normalized flux
coeffs = (500, 50)                          # a, b for R(λ) = a*λ + b

convolved_flux = fastconv_VariR(wavelengths, flux, coeffs, quiet=False)

```
## Citation
The details can be found at [https://arxiv.org/abs/2605.30871](https://arxiv.org/abs/2605.30871)
If you use this code in your research, please cite our paper.
