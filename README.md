# IterativeClosestPoint (variants)

## Introduction

Small experimental project for performing experiments with
registration of surfaces using different approaches than the
traditional Iterative Closest Point (ICP), which is used by
`vtkIterativeClosestPoint` from the VTK repository. Many variants
exists - as a start only only alternative is implemented. Instead of
minimizing point distance, the point-to-plane distance is minimized,
when registering the point $p$ to $q$.

$\|(R(\theta)p + t - q)\cdot q\|,$

where $R$ and $t$ is the rotation and translation, which are
optimized.

This is a small snapshot of larger project based on the VTK repository
and meant primarily for educational purposes.

## Implementation

**vtkImplicitPolyDataDistance2** A thread-safe version of the original
`vtkImplicitPolyDataDistance`, where it is possible to return both the
closest point and the gradient. It is pretty slow and should be
optimized using meta-programming at some point.

**vtkPolyDataCorrespondenceFilter** A custom filter, which can be used
for finding correspondences between two `vtkPolyData`. It is
implemented using `vtkSMPTools` to speed-up the performance.

**vtkICP** Registration algorithm inheriting from `vtkLinearTransform`
similarly to `vtkIterativeClosestPoint`. Right now, a point-to-point
metric, a point-to-plane metric as well as a symmetrized version of
the point-to-plane metric by Szymon Rusinkiewicz
[symmetrized](https://gfx.cs.princeton.edu/pubs/Rusinkiewicz_2019_ASO/symm_icp.pdf)
can be chosen. The first algorithm is optimized using the standard
closed form SVD and the latter ones are optimized using simple
steepest descent. A Gauss-Newton version will follow at some point
supporting various priors. For both algorithm a distance threshold can
be set as well as the maximum number of landmarks.

## Build and run

It is implemented following the directory and build structure of
VTK. No attention has been given to performance (yet).

## Example

There is a small Python example in
[RegistrationDemo](./ICP/Testing/Python/RegistrationDemo.py) demonstrating how
registration using a point-to-plane metric converges in much fewer
iterations than the conventional point-to-point metric.
