#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

if sys.version_info < (3,5):
  raise Exception("Unsupported python version")
else:
  # TODO: Do something similar to VTK
  import importlib
  # import vtkmodules.all
  all_m = importlib.import_module('icpmodules.all')

  # import icpmodules
  icpmodules_m = importlib.import_module('icpmodules')

  # make icpmodules.all act as the icpmodules package to support importing
  # other modules from icpmodules package via `vtk`.
  all_m.__path__ = icpmodules_m.__path__
  all_m.__version__ = icpmodules_m.__version__

  # replace old `icp` module with the `all` package.
  sys.modules[__name__] = all_m
