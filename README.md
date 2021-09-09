# 2-D isopycnal model for multi-oxygen-isotopologues of Pacific ocean

1. The Isopycnal_surface_gsw (folder) utilize the gsw-tool kit (https://github.com/TEOS-10/GSW-Matlab) to calculate the O<sub>2</sub> sat in a specific density in the modelled region in Pacific ocean.

2. Python version implements this 2-D isopycnal advection diffusion reaction modelt to simulate the steady state solution. However this model does not support GPU related acceleration. Advanced Numba or other scientific JIT compiler might be need to increase the running efficiency.

3. Grid search implements the proposed gird search function on mapping O<sub>2</sub> measured in isopycnal_sirface_gsw into our model regime.

4. Isopycnal_model file implements the multi-isotopologue enabled 2D-advection diffusion reaction isopycnal model (MATLAB, gpu with CUDA is a must).