add_executable(correlations correlations.cpp Heisenberg1D.cpp)
target_include_directories(correlations PRIVATE ../pcg_random/include)