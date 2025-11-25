# QEDFeynmanDiagrams.jl

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://QEDjl-project.github.io/QEDFeynmanDiagrams.jl/dev/)
[![DOI](https://rodare.hzdr.de/badge/DOI/10.14278/rodare.3713.svg)](https://doi.org/10.14278/rodare.3713)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

Generator for QED Feynman diagrams and [`ComputableDAGs.jl`](https://github.com/ComputableDAGs/ComputableDAGs.jl) to compute scattering processes' matrix elements.

## Testing

This project can be tested using julia's `Pkg.test()`. Additionally, the behavior can be modified via some environment variables:
- `TEST_<GPU> = 1`: Enables GPU tests for a GPU vendor, available are `CUDA`, `AMDGPU`, `ONEAPI`, and `METAL`. Disabled by default. Make sure the relevant libraries are installed on the executing machine.
- `TEST_CPU = 0`: Disables CPU tests. Useful when you only want to run GPU tests, for example in the CI. Enabled by default.
- `LARGE_TESTS = 1`: Whether to run large tests. Some tests can be slow and can be enabled or disabled with this variable. This is set to `CI` by default, which enables them locally and disables them in the CI (i.e., when the `CI` environment variable is set).

## Acknowledgements and Funding

This work was partly funded by the Center for Advanced Systems Understanding (CASUS) that is financed by Germany’s Federal Ministry of Research, Technology and Space (BMFTR) and by the Saxon Ministry for Science, Culture and Tourism (SMWK) with tax funds on the basis of the budget approved by the Saxon State Parliament.

I'd also like to thank Michael Bussmann for funding the project, and Simeon Ehrig for help with the design and maintenance of the package.
