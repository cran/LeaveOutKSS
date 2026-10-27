# LeaveOutKSS 0.1.1

- Fixed the exact leverage branch of `leverages()` and `leverages_parallel()`.
  The `Bii_fe` and `Bii_pe` bias terms were divided by `NT - 1` instead of
  being sums of squared deviations, and `Bii_cov` used the wrong vectors.
  Results obtained with `type_algorithm = "exact"` were therefore effectively
  uncorrected for the variance components and mis-signed for the covariance.
  The default `type_algorithm = "JLA"` path was not affected. Reported by
  santiagohermo (GitHub).
- Added a regression test comparing the exact leverage branch against a dense
  ground truth.
- Capped the iteration count passed to `sanic::solve_cg()` at the dimension of
  the normal equations. 'sanic' 0.1.0 rejects larger values, which made every
  solver call fail on models with fewer than 300 (leverages) or 1000 (fixed
  effect estimation) coefficients.

# LeaveOutKSS 0.1.0

- Includes sparse-matrix graph adjustments for larger datasets.
- Improves graph construction and `lincom` warning handling.
