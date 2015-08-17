#include "writeHDF5.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define SRC_LOC __FILE__ ":" TOSTRING(__LINE__)

#define _2D 2

// TODO: This function can general enough to write any buffer to a file with a given name.
void write_HDF5_JPLC(grid_parms* grid, double *jplc_buffer, char *path)
{
	// printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/jplc_%d.h5", path, grid->branch_tag);

	printf("Writing JPLC file: %s\n", filename);

	hid_t file_id;
	hid_t space_id;
	hid_t dset_id;
	herr_t status;

	// Create a HDF5 file.
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	hsize_t dims[_2D] = {grid->num_ec_axially * grid->m, grid->num_ec_circumferentially * grid->n};

	// Create dataspace.
	space_id = H5Screate_simple(_2D, dims, NULL);

	// Create dataset.
	dset_id = H5Dcreate(file_id, "/jplc", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	// Write dataset.
	status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, jplc_buffer);

	// Close everything.
	status = H5Dclose(dset_id);
	status = H5Sclose(space_id);
	status = H5Fclose(file_id);
}

void write_EC_data_HDF5(grid_parms* grid, ec_data_buffer *ec_buffer, int write_count, char* path)
{
	printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/ec_data_t_%d_b_%d.h5", path, write_count, grid->branch_tag);

	printf("Writing EC data file: %s\n", filename);

	hid_t file_id;

	hid_t space_id;
	hid_t space_id_1;
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
	space_id = H5Screate_simple(_2D, dims, NULL);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	// Create dataset.
	ec_Ca_id = H5Dcreate(file_id, "/ec_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_Ca);

	// Create dataset.
	ec_cpl_Ca_id = H5Dcreate(file_id, "/ec_cpl_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_cpl_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_Ca);

	// Create dataset.
	ec_IP3_id = H5Dcreate(file_id, "/ec_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_IP3);

	// Create dataset.
	ec_cpl_IP3_id = H5Dcreate(file_id, "/ec_cpl_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_cpl_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_IP3);

	// Create dataset.
	ec_Vm_id = H5Dcreate(file_id, "/ec_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_Vm);

	// Create dataset.
	ec_cpl_Vm_id = H5Dcreate(file_id, "/ec_cpl_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_cpl_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_cpl_Vm);

	// Create dataset.
	ec_SR_id = H5Dcreate(file_id, "/ec_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(ec_SR_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ec_buffer->_ec_SR);

	// Close everything.
	status = H5Dclose(ec_Ca_id);
	status = H5Dclose(ec_SR_id);
	status = H5Dclose(ec_Vm_id);
	status = H5Dclose(ec_IP3_id);
	status = H5Dclose(ec_cpl_Ca_id);
	status = H5Dclose(ec_cpl_Vm_id);
	status = H5Dclose(ec_cpl_IP3_id);

	status = H5Sclose(space_id);

	status = H5Fclose(file_id);
}

void write_SMC_data_HDF5(grid_parms* grid, smc_data_buffer *smc_buffer, int write_count, char* path)
{
	printf("[%d] >>>>>> Entering %s:%s\n", grid->universal_rank, __FILE__, __FUNCTION__);

	char filename[256];
	int err = sprintf(filename, "%s/smc_data_t_%d_b_%d.h5", path, write_count, grid->branch_tag);

	printf("Writing EC data file: %s\n", filename);

	hid_t file_id;

	hid_t space_id;
	hid_t space_id_1;
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
	space_id = H5Screate_simple(_2D, dims, NULL);

	// printf("[%d] Checkpoint %s:%d\n", grid->universal_rank, __FILE__, __LINE__);

	// Create dataset.
	smc_Ca_id = H5Dcreate(file_id, "/SMC_Ca", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_Ca);

	// Create dataset.
	smc_cpl_Ca_id = H5Dcreate(file_id, "/SMC_Ca_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_cpl_Ca_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_Ca);

	// Create dataset.
	smc_IP3_id = H5Dcreate(file_id, "/SMC_IP3", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_IP3);

	// Create dataset.
	smc_cpl_IP3_id = H5Dcreate(file_id, "/SMC_IP3_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_cpl_IP3_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_IP3);

	// Create dataset.
	smc_Vm_id = H5Dcreate(file_id, "/SMC_Vm", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_Vm);

	// Create dataset.
	smc_cpl_Vm_id = H5Dcreate(file_id, "/SMC_Vm_coupling", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_cpl_Vm_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_cpl_Vm);

	// Create dataset.
	smc_SR_id = H5Dcreate(file_id, "/SMC_SR", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_SR_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_SR);

	// Create dataset.
	smc_SR_id = H5Dcreate(file_id, "/SMC_w", H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	// Write dataset.
	status = H5Dwrite(smc_W_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, smc_buffer->_smc_W);

	// Close everything.
	status = H5Dclose(smc_Ca_id);
	status = H5Dclose(smc_SR_id);
	status = H5Dclose(smc_W_id);
	status = H5Dclose(smc_Vm_id);
	status = H5Dclose(smc_IP3_id);
	status = H5Dclose(smc_cpl_Ca_id);
	status = H5Dclose(smc_cpl_Vm_id);
	status = H5Dclose(smc_cpl_IP3_id);

	status = H5Sclose(space_id);

	status = H5Fclose(file_id);
}
