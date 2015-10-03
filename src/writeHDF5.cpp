#include "writeHDF5.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

/*
 * There are APIs to determine if datasets and groups are left open.
 * H5Fget_obj_count will get the number of open objects in the file, and
 * H5Fget_obj_ids will return a list of the open object identifiers.
 */

#define _2D 2

// TODO: This function can be made general enough to write any buffer to a file with a given name.
void write_HDF5_JPLC(grid_parms* grid, double *jplc_buffer, char *path)
{
	// printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/jplc_%d.h5", path, grid->branch_tag);

	printf("[%d] Writing JPLC file: %s\n", grid->universal_rank, filename);

	hid_t file_id;
	hid_t space_id;
	hid_t dset_id;

	herr_t status;

#if 0
	hid_t fapl_id;
	hid_t plist_id;

	MPI_Info info;
	// This code is only good for a single writer when we have a single trunk, not a bifurcation.
	MPI_Group world_group, writer_group;
	MPI_Comm writer_comm;
	int ranks[1] = {0};

	MPI_Info_create(&info);
	MPI_Info_set(info, "IBM_largeblock_io", "true");
	MPI_Info_set(info, "stripping_unit", "4194304");

#if 0
	// http://lists.hdfgroup.org/pipermail/hdf-forum_lists.hdfgroup.org/2014-November/008210.html
	// http://www-01.ibm.com/support/knowledgecenter/SSFK3V_1.3.0/com.ibm.cluster.pe.v1r3.pe500.doc/am107_ifopen.htm?lang=en
	MPI_INFO_SET(info,"H5F_ACS_CORE_WRITE_TRACKING_PAGE_SIZE_DEF","524288",error)
	MPI_INFO_SET(info,"ind_rd_buffer_size","41943040", error)
	MPI_INFO_SET(info,"ind_wr_buffer_size","5242880", error)
	MPI_INFO_SET(info,"romio_ds_read","disable", error)
	MPI_INFO_SET(info,"romio_ds_write","disable", error)
	MPI_INFO_SET(info,"romio_cb_write","enable", error)
	MPI_INFO_SET(info,"cb_buffer_size","4194304", error)
#endif

	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	MPI_Group_excl(world_group, 1, ranks, &writer_group);
	MPI_Comm_create(MPI_COMM_WORLD, writer_group, &writer_comm);

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, writer_comm, info);

	// Remember this at the end of this function:
	MPI_Info_free(&info);
	MPI_Group_free(&world_group);
	MPI_Group_free(&writer_group);
	MPI_Comm_free(&writer_comm);

	// Create a HDF5 file.
	// hid_t H5Fcreate( const char *name, unsigned flags, hid_t fcpl_id, hid_t fapl_id )
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
#endif

	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t dims[_2D] = {grid->num_ec_axially * grid->m, grid->num_ec_circumferentially * grid->n};

	// Create dataspace.
	// hid_t H5Screate_simple(int rank, const hsize_t * dims, const hsize_t * maxdims )
	space_id = H5Screate_simple(_2D, dims, dims);

	// Create dataset.
	// hid_t H5Dcreate( hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t dcpl_id )
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	dset_id = H5Dcreate(file_id, "/jplc", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	dset_id = H5Dcreate(file_id, "/jplc", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif

	// Write dataset.
	// herr_t H5Dwrite(hid_t dataset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t xfer_plist_id, const void * buf )
	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, jplc_buffer);
	CHECK(status, FAIL, __FUNCTION__);

	// Close everything.
	status = H5Dclose(dset_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Sclose(space_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Fclose(file_id); CHECK(status, FAIL, __FUNCTION__);

	// printf("[%d] <<<<<< Leaving %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);
}

void write_EC_data_HDF5(grid_parms* grid, ec_data_buffer *ec_buffer, int write_count, char* path)
{
	// printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/ec_data_t_%d_b_%d.h5", path, write_count, grid->branch_tag);

	printf("[%d] Writing EC file: %s\n", grid->universal_rank, filename);

	hid_t file_id;

	hid_t space_id;
	hid_t ec_Ca_id;
	hid_t ec_cpl_Ca_id;
	hid_t ec_IP3_id;
	hid_t ec_cpl_IP3_id;
	hid_t ec_Vm_id;
	hid_t ec_cpl_Vm_id;
	hid_t ec_SR_id;

	herr_t status;

	// Create a HDF5 file.
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t dims[_2D] = {grid->num_ec_axially * grid->m, grid->num_ec_circumferentially * grid->n};

	// Create dataspace.
	space_id = H5Screate_simple(_2D, dims, dims);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_Ca_id = H5Dcreate(file_id, "/EC_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_Ca_id = H5Dcreate(file_id, "/EC_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_Ca);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_cpl_Ca_id = H5Dcreate(file_id, "/EC_Ca_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_cpl_Ca_id = H5Dcreate(file_id, "/EC_Ca_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_cpl_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_Ca);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_IP3_id = H5Dcreate(file_id, "/EC_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_IP3_id = H5Dcreate(file_id, "/EC_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_IP3);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_cpl_IP3_id = H5Dcreate(file_id, "/EC_IP3_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_cpl_IP3_id = H5Dcreate(file_id, "/EC_IP3_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_cpl_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_IP3);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_Vm_id = H5Dcreate(file_id, "/EC_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_Vm_id = H5Dcreate(file_id, "/EC_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_Vm);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_cpl_Vm_id = H5Dcreate(file_id, "/EC_Vm_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_cpl_Vm_id = H5Dcreate(file_id, "/EC_Vm_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_cpl_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_Vm);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	ec_SR_id = H5Dcreate(file_id, "/EC_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	ec_SR_id = H5Dcreate(file_id, "/EC_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(ec_SR_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_SR);
	CHECK(status, FAIL, __FUNCTION__);

	// Close everything.
	status = H5Dclose(ec_Ca_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_SR_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_Vm_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_IP3_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_cpl_Ca_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_cpl_Vm_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(ec_cpl_IP3_id); CHECK(status, FAIL, __FUNCTION__);

	status = H5Sclose(space_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Fclose(file_id); CHECK(status, FAIL, __FUNCTION__);
}

void write_SMC_data_HDF5(grid_parms* grid, smc_data_buffer *smc_buffer, int write_count, char* path)
{
	// printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/smc_data_t_%d_b_%d.h5", path, write_count, grid->branch_tag);

	printf("[%d] Writing SMC file: %s\n", grid->universal_rank, filename);

	hid_t file_id;

	hid_t space_id;
	hid_t smc_Ca_id;
	hid_t smc_cpl_Ca_id;
	hid_t smc_IP3_id;
	hid_t smc_cpl_IP3_id;
	hid_t smc_Vm_id;
	hid_t smc_cpl_Vm_id;
	hid_t smc_SR_id;
	hid_t smc_W_id;

	herr_t status;

	// Create a HDF5 file.
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t dims[_2D] = {grid->num_smc_axially * grid->m, grid->num_smc_circumferentially * grid->n};

	// Create dataspace.
	space_id = H5Screate_simple(_2D, dims, dims);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_Ca_id = H5Dcreate(file_id, "/SMC_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_Ca_id = H5Dcreate(file_id, "/SMC_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_Ca);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_cpl_Ca_id = H5Dcreate(file_id, "/SMC_Ca_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_cpl_Ca_id = H5Dcreate(file_id, "/SMC_Ca_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_cpl_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_Ca);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_IP3_id = H5Dcreate(file_id, "/SMC_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_IP3_id = H5Dcreate(file_id, "/SMC_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_IP3);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_cpl_IP3_id = H5Dcreate(file_id, "/SMC_IP3_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_cpl_IP3_id = H5Dcreate(file_id, "/SMC_IP3_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_cpl_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_IP3);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_Vm_id = H5Dcreate(file_id, "/SMC_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_Vm_id = H5Dcreate(file_id, "/SMC_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_Vm);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_cpl_Vm_id = H5Dcreate(file_id, "/SMC_Vm_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_cpl_Vm_id = H5Dcreate(file_id, "/SMC_Vm_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_cpl_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_Vm);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_SR_id = H5Dcreate(file_id, "/SMC_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_SR_id = H5Dcreate(file_id, "/SMC_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_SR_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_SR);
	CHECK(status, FAIL, __FUNCTION__);

	// Create dataset.
#if H5_VERS_MAJOR == 1 && H5_VERS_MINOR > 6
	smc_W_id = H5Dcreate(file_id, "/SMC_w", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#else
	smc_W_id = H5Dcreate(file_id, "/SMC_w", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT);
#endif
	// Write dataset.
	status = H5Dwrite(smc_W_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_W);
	CHECK(status, FAIL, __FUNCTION__);

	// Close everything.
	status = H5Dclose(smc_Ca_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_SR_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_W_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_Vm_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_IP3_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_cpl_Ca_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_cpl_Vm_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Dclose(smc_cpl_IP3_id); CHECK(status, FAIL, __FUNCTION__);

	status = H5Sclose(space_id); CHECK(status, FAIL, __FUNCTION__);
	status = H5Fclose(file_id); CHECK(status, FAIL, __FUNCTION__);
}
