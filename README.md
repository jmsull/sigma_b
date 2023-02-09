# sigma_b
*Very* simple code to compute the pesky normalization factor in super-sample covariance w/ trapz rule for a simulation box

$\sigma_{b}^{2} = \int \frac{d^{3}\mathbf{k}}{(2\pi)^3} |W(\mathbf{k})|^{2} P_{L}(k)$

where 

$W^{2}(\mathbf{k}) = j_{0}^{2}(k_{x}L/2)j_{0}^{2}(k_{y}L/2)j_{0}^{2}(k_{z}L/2)$

which is appropriately normalized to $1/V = L^{-3}$.

# Notes: 

Run this code with ``julia sigmab.jl`` 

The integration range and number of points appear to be decent at a $<3$% level  (using window error as a guide) but are actually much better than this for estimating the window integral for a typical cosmology - that integral peaks at very low k, whereas the window integral doesn't converge quickly b/c at high k we've entered the wiggle zone. If we you want a more accurate estimate you can increase kmax or N (at your own runtime risk).

Default ouput on my machine is:
``Is the window close to 1? true
Power spectrum variance in a cubic box of side length 625.0 Mpc/h is 6.756502209257932e-5.``

