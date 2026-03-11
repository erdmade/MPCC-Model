# MPCCmodel

MPCCmodel is a single column climate model that runs C code for atmospheric and radiative transfer calculations. Atmospheric layer number, trace gas concentrations, and clouds can be user-defined.

## Dependencies

- **GCC**  
- **NetCDF C library** 

### Installing NetCDF dependency (Ubuntu/Debian)

```bash
sudo apt update
sudo apt install libnetcdf-dev 
```

## Compilation

To build the program, run:

```bash
make
```