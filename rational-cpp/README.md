# rational-cpp

[![Build](https://github.com/zfergus/rational-cpp/actions/workflows/continuous.yml/badge.svg)](https://github.com/zfergus/rational-cpp/actions/workflows/continuous.yml)

C++ Wrapper of GMP's Rational Type

## Usage

### Adding to CMake Project

The easiest way to add rational-cpp to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) (requires CMake ≥ 3.14).

For example,

```cmake
include(FetchContent)
FetchContent_Declare(
    rational_cpp
    GIT_REPOSITORY https://github.com/zfergus/rational-cpp.git
    GIT_TAG ${RATIONAL_CPP_GIT_TAG}
)
FetchContent_MakeAvailable(rational_cpp)
```

where `RATIONAL_CPP_GIT_TAG` is set to the version of the library you want to use. This will download and add rational-cpp to CMake. Rational-cpp can then be linked against using

```cmake
target_link_libraries(${PROJECT_NAME} PUBLIC rational::rational)
```

where `PROJECT_NAME` is the name of your library/binary.

### Dependencies

The only dependency is [GMP](https://gmplib.org/), which must be installed at a system level.

### Example: Reading Rational Numbers from Strings

```cpp
#include <rational/rational.hpp>

rational::Rational r("1", "2");
assert(r == 0.5);
```

### Example: Using Rational with `Eigen::Matrix`/`Eigen::Array`

To use Rational as the data type with Eigen's Matrix/Array structures, you need to use the `Eigen::DontAlign` template argument. For example:

```cpp
typedef Eigen::Matrix<Rational, 3, 1, Eigen::ColMajor | Eigen::DontAlign> Vector3r;
typedef Eigen::Array<Rational, 3, 1, Eigen::ColMajor | Eigen::DontAlign> Array3r;
```

## Acknowledgements

This is a modified version of `rational-gmp.cpp` from [CVC3](https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/rational-gmp_8cpp-source.html) by Sergey Berezin.
The licensed for which can be found [here](https://cs.nyu.edu/acsys/cvc3/releases/1.5/doc/LICENSE.html).

## License

MIT License © 2023, Zachary Ferguson (See <a href="https://github.com/zfergus/rational-cpp/blob/main/LICENSE"><code>LICENSE</code></a> for details).
