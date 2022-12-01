# TriangulationProject


I made a cmake project to better understand triangulation methods. The code is adapted from https://imkaywu.github.io/blog/.
I don't have the standard folder structure because of the small nature of the project. Maybe I'll reorganize in the future.
It wasn't necessary to create an object for the functions, I just did it to do it.


Requirements:
cmake above 3.0.0 and Eigen 3.4

You can change the cmake version and Eigen version in CMakeLists.txt and try to build if you have different versions.

Use the standard cmake commands to use the repository with the given test file in Ubuntu:

from the TriangulationProject folder
mkdir build && cd build  
cmake ..  
cmake --build . or make -j#  

./TriangulationProject to run