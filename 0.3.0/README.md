# `0.3.0`

## CircleCI

Build `321`, Logs gone

## AppVeyor

Build `1.0.54.master`

### Build Matrix

| Job ID             | `TARGET_ARCH` | `CONDA_PY` | `PY_CONDITION`       | `CONDA_INSTALL_LOCN`  | `TEST_TYPE`  |
|--------------------|---------------|------------|----------------------|-----------------------|--------------|
| `g14ut3aavk81f1u5` | `x86`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda`       | `base`       |
| `9cb5ay6x1kwddkeh` | `x86`         | `34`       | `python=3.4`         | `C:\\Miniconda3`      | `base`       |
| `1cwo8mk0u3w8s6ub` | `x64`         | `34`       | `python=3.4`         | `C:\\Miniconda3-x64`  | `base`       |
| `2l4m3smre0xpoarm` | `x86`         | `35`       | `python=3.5`         | `C:\\Miniconda35`     | `base`       |
| `7yu9in1n1l7nqd6d` | `x64`         | `35`       | `python=3.5`         | `C:\\Miniconda35-x64` | `base`       |
| `lm74avs2bkrc979x` | `x64`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda-x64`   | `cover`      |
| `qtm245qf9eqy74c9` | `x64`         | `27`       | `python >=2.7.13,<3` | `C:\\Miniconda-x64`   | `functional` |
| `35cjrlpbwicop266` | `x64`         | `35`       | `python=3.5`         | `C:\\Miniconda35-x64` | `doctest`    |
