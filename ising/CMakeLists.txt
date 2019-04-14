add_executable(observables observables.cpp Ising2D.cpp Ising2D.h)
target_include_directories(observables PRIVATE ../pcg_random/include)
add_executable(lattice lattice.cpp Ising2D.cpp Ising2D.h)
target_include_directories(lattice PRIVATE ../pcg_random/include)
