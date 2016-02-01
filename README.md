# Non-linear Schrodinger Equation Solver(NLS)

## Description

Package to solve non-linear Schrodinger equation written with Fortran 95 and Python 2.

### Features

- core written with native Fortran
- wrapped with Python interaface

## Usage

### Fortran executable

1. Change directory to `src/`
2. Run `make`
3. Change directory to `test/`
4. Run `make` and then run `./test_nls`
5. Change directory to `bin/`
6. Run `./solve` in order to start calculation
7. Visualize solution with `python2 src/visualize.py`

### Python module

1. Change directory to `src/`
2. Run `make glue` in order to build python lib in `bin/` directory that incapsulate native code
3. Use python to import native module in a way


```python
    from nls.native import nls
    print(nls.__doc__)
```

4. Or use python to import python wrap in a way

```python
    from nls import Problem
    model = Problem().model()
    model.solve()
    model.visualize()
    model.show()
    model.report()
```

## Testing

## Requirements

1. `make`
2. `gfortran`
3. `numpy` with `f2py` extension
4. `matplotlib`

## Credits 

(c) Daniel Bershatsky <dainel.bershatsky@skolkovotech.ru>, 2015


