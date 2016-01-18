# Non-linear Schrodinger Equation Solver(NLS)

## Description

Package to solve non-linear Schrodinger equation written with Fortran 95 and Python 2.

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
2. Run `make glue` in order to build python lib in `bin/` directory
3. Use python to import the module in a way

    from nls import nls
    print(nls.__doc__)

## Testing

## Requirements

1. `make`
2. `gfortran`
3. `numpy` with `f2py` extension
4. `jupythor`

## Credits 

(c) Daniel Bershatsky <dainel.bershatsky@skolkovotech.ru>, 2015


