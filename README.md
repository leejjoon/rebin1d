# Simple 1-d rebin module

- It rebins 1-d spectrum into new wavelength grid.
- Compared to interpolation,
  - fluxes are conserved.
  - handling missing values can be easier
- The points are interpolated and then integrated with new bins.
- In practice, we interpolate on the accumulated values of points.
- Two Interpolation methos is supported : nearest neighbor and linear.
- For the nearest neighbor interpolation, optional drizzling can be used.
