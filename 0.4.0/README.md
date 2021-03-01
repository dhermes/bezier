# `0.4.0`

## CircleCI

Build `490`, Logs gone

## AppVeyor

Build `1.0.225.master`

### Build Matrix

| Job ID             | `TARGET_ARCH` | `CONDA_PY` | `PY_CONDITION`       | `CONDA_INSTALL_LOCN`  | `TEST_TYPE`  |
|--------------------|---------------|------------|----------------------|-----------------------|--------------|
| `edb0vmxs4kg42ok8` | `x86`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda`       | `base`       |
| `6v92c8528xc1fa9s` | `x86`         | `35`       | `python=3.5`         | `C:\\Miniconda35`     | `base`       |
| `3bnvrvp8d65htua3` | `x64`         | `35`       | `python=3.5`         | `C:\\Miniconda35-x64` | `base`       |
| `i1fst2yl4oeo24vc` | `x86`         | `36`       | `python=3.6`         | `C:\\Miniconda36`     | `base`       |
| `afh95yg2rgc2nw6q` | `x64`         | `36`       | `python=3.6`         | `C:\\Miniconda36-x64` | `base`       |
| `neq7aj86yvsan0wp` | `x64`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda-x64`   | `cover`      |
| `0bujris4fx80wd1e` | `x64`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda-x64`   | `functional` |
| `qvb5r6eappti0jqi` | `x64`         | `35`       | `python=3.5`         | `C:\\Miniconda35-x64` | `doctest`    |
