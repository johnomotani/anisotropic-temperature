To build (replacing the path with the path to a built copy of BOUT++):
```
$ cmake . -B build -DCMAKE_PREFIX_PATH=/path/to/built/BOUT++
$ cmake --build build
```

To run:
```
$ cd build
$ ./anisotropic-temperature
```

To visualise the results
```
$ ./animate-results
```

Version info
------------
These examples were tested with:
BOUT++ 9445ace7f50b905c9791bbd914a11a56e92979d9
xBOUT 12c0088a27d9554b01a4f0ac4e2658b65329b9d0
