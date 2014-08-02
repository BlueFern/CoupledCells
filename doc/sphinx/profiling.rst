Performance/Time Profiling
==========================

This is work in progress...

There are two levels of time profiling of different granularity.

 * Coarser level of profiling is provided in the function :apilink:`rksuite_solver_CT(double tnow, double tfinal, double interval, double *y, double* yp, int total, double TOL, double* thres, int file_write_per_unit_time, int line_number, checkpoint_handle *check, char* path, IO_domain_info* my_IO_domain_info)`, which calls :apilink:`checkpoint_coarse_time_profiling_data(grid_parms grid, time_stamps* t_stamp, IO_domain_info* my_IO_domain_info)`.
 * Finer level of profiling is provided by the function :apilink:`compute_with_time_profiling(time_stamps* t_stamp, grid_parms grid, celltype1** smc, celltype2** ec, conductance cpl_cef, double t, double* y, double* f)`.

.. todo:: Simple search in the source returns numerous other cases of ``MPI_Wtime()`` calls, which looks like time profiling is not limited to the functions listed above. Perhaps we need to define a CMake variable which can be used in the source whether time profiling is to be used at runtime.

.. todo:: :apilink:`compute_with_time_profiling(time_stamps* t_stamp, grid_parms grid, celltype1** smc, celltype2** ec, conductance cpl_cef, double t, double* y, double* f)` function duplicates code from :apilink:`compute(grid_parms grid, celltype1** smc, celltype2** ec, conductance cpl_cef, double t, double* y, double* f)`. Perhaps the functions can be merged and time profiling enabling can be controlled with CMake variables at configuration time.
