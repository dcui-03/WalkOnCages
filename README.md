# WalkOnCages

This repo contains my from-scratch implementation of Stochastic Computation of Barycentric Coordinates [de Goes and Desbrun, ACM SIGGRAPH 2024] in C++. Please see Notes below for more information.

**Build Instructions**

You only need to clone Polyscope into a folder called deps; Eigen is fetched automatically by CMake. This has been tested on Windows (MSVC), Linux, and macOS (Apple Silicon):

```
git clone https://github.com/dcui-03/WalkOnCages
cd ./WalkOnCages
mkdir deps && cd ./deps
git clone --recurse-submodules https://github.com/nmwsharp/polyscope.git
cd ..
```

*macOS / Linux*

Presuming you would like to use OpenMP for multi-threading:
```
cmake -S . -B build   -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++   -DOpenMP_ROOT=$(brew --prefix libomp)
cd ./build && make -j4
```
Otherwise,
```
mkdir build && cd ./build
cmake ..
make -j4
```

To run on Linux, use the following command from root if using a mesh-based cage:
```
./build/stochastic_bc <path to query mesh> <path to cage mesh>
```
Or use the following alternative for a curvenet-based cage:
```
./build/stochastic_bc <path to query mesh>
```

*Windows (MSVC)*
```
cmake -S . -B build
cmake --build build --config Release
```
To run on MSVC, use the following command from root:
```
.\build\Release\stochastic_bc.exe <path to query mesh> <path to cage mesh>
```

**Notes**

This code is for a personal project and learning purposes, and is therefore not entirely optimized. As a result, it is also not suitable as a replacement for commercial software: the visualization uses Polyscope [Sharp, 2019] to test basic features and functionality on toy examples. If you would like to contribute, collaborate, or provide feedback, please contact me directly.

*Please note that this is an unofficial implementation of the methods of Stochastic Computation of Barycentric Coordinates [de Goes and Desbrun, 2024], and I am not affiliated with Pixar Animation Studios nor the Walt Disney Company. Certain underlying methods may be subject to third-party patent rights. The MIT License applies to the source code in this repository and does not grant any rights under third-party patents.*
