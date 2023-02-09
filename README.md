# sigma_b
Compute pesky normalization factor in super-sample covariance w/ trapz rule for a simulation box

$\sigma_{b}^{2} = \frac{1}{V^2} \int \frac{d^{3}\mathbf{k}}{(2\pi)^3} |W(\mathbf{k})|^{2} P_{L}(k)$

where 

$V=\int d^{3}\mathbf{x} W(\mathbf{x})$

and

$W(\mathbf{k}) = j_{0}(k_{x})j_{0}(k_{y})j_{0}(k_{z})$
