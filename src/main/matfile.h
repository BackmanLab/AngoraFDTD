/* AUTORIGHTS
Copyright (C) 2006-2012  Ilker R. Capoglu

    This file is part of the Angora package.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MATFILE_H
#define MATFILE_H

//template declaration for the function that reads a material region from a file

#include "headers.h"

//base Angora exception class
#include "angora_excp.h"

#include "material/Cmat_types.h"

#include "shape/Cshape.h"

#include <fstream>

//#include "matlab.hpp"
#include "extras.h" /* Aya sub functions */

/* Aya */
#include <mpi.h>     /* Aya MPI and MPI-IO live here */
#include <stdio.h>   /* Aya all IO stuff lives here */
#include <inttypes.h>
#define MASTER_RANK 0
#define TRUE 1
#define FALSE 0
#define BOOLEAN int /* aya */
#define MBYTE 1048576
#define DEBUG

#ifndef ANGORA_MAX_NEWMAT
#define ANGORA_MAX_NEWMAT 1000
#endif

extern double dx,dt;

extern int rank;

extern int iback,ifront;
extern int jleft,jright;
extern int klower,kupper;
extern int nodes_x, nodes_y, nodes_z;
//extern int NCELLS_X,NCELLS_Y,NCELLS_Z,NPML;

//extern int num_of_distinct_eps_x,num_of_distinct_eps_y,num_of_distinct_eps_z;
//extern int num_of_distinct_mu_x,num_of_distinct_mu_y,num_of_distinct_mu_z;
//extern int num_of_distinct_cond_e_x,num_of_distinct_cond_e_y,num_of_distinct_cond_e_z;
//extern int num_of_distinct_cond_h_x,num_of_distinct_cond_h_y,num_of_distinct_cond_h_z;

extern Array<update_coeff_type,3> Ca_X,Cb_X,Cc_X,Ca_Y,Cb_Y,Cc_Y,Ca_Z,Cb_Z,Cc_Z;
extern Array<update_coeff_type,3> Da_X,Db_X,Da_Y,Db_Y,Da_Z,Db_Z;

extern Array<eps_x_type,3> eps_x_indices;
extern Array<eps_y_type,3> eps_y_indices;
extern Array<eps_z_type,3> eps_z_indices;
extern Array<mu_x_type,3> mu_x_indices;
extern Array<mu_y_type,3> mu_y_indices;
extern Array<mu_z_type,3> mu_z_indices;
extern Array<cond_e_x_type,3> cond_e_x_indices;
extern Array<cond_e_y_type,3> cond_e_y_indices;
extern Array<cond_e_z_type,3> cond_e_z_indices;
extern Array<cond_h_x_type,3> cond_h_x_indices;
extern Array<cond_h_y_type,3> cond_h_y_indices;
extern Array<cond_h_z_type,3> cond_h_z_indices;

extern Array<float,1> eps_x,eps_y,eps_z;
extern Array<float,1> mu_x,mu_y,mu_z;
extern Array<float,1> cond_e_x,cond_e_y,cond_e_z;
extern Array<float,1> cond_h_x,cond_h_y,cond_h_z;

extern Array<bool,3> dispersion_exists_at_Ex_position,dispersion_exists_at_Ey_position,dispersion_exists_at_Ez_position;
extern Array<update_coeff_type,4> alpha_X,xi_X,gamma_X,alpha_Y,xi_Y,gamma_Y,alpha_Z,xi_Z,gamma_Z;
extern Array<omega_p_x_type,4> omega_p_x_indices;
extern Array<omega_p_y_type,4> omega_p_y_indices;
extern Array<omega_p_z_type,4> omega_p_z_indices;
extern Array<tau_p_x_type,4> tau_p_x_indices;
extern Array<tau_p_y_type,4> tau_p_y_indices;
extern Array<tau_p_z_type,4> tau_p_z_indices;
extern Array<Omega_p_x_type,4> Omega_p_x_indices;
extern Array<Omega_p_y_type,4> Omega_p_y_indices;
extern Array<Omega_p_z_type,4> Omega_p_z_indices;
extern Array<float,1> omega_p_x,omega_p_y,omega_p_z,
                      tau_p_x,tau_p_y,tau_p_z,
                      Omega_p_x,Omega_p_y,Omega_p_z;


template<typename MatType> //data type for the material property
void PlaceMaterialRegionFromFile(const string& MaterialFileName, const int& xPos, const int& yPos, const int& zPos, const string& anchor, const string& constitutive_param_type, const int& max_number_of_new_materials = ANGORA_MAX_NEWMAT, const_Cshape_shared_ptr shape_mask = const_Cshape_shared_ptr(new Cuniverse()))
{//Reads rectangular-prism-shaped dielectric region from file and places into grid
// (xPos,yPos,zPos) are the x-y-z coordinates of the anchor of the region (measured in cells from the back-left-lower corner of the grid)


	if ((constitutive_param_type!="rel_permittivity")&&(constitutive_param_type!="rel_permeability")&&(constitutive_param_type!="electric_conductivity")&&(constitutive_param_type!="magnetic_conductivity")&&(constitutive_param_type!="lorentz_delta_epsilon"))
	{
#ifdef __GNUG__
//GNU C++ compiler is being used, use the nice predefined variables for the function name
//		InvalidNumericArgumentException<double> exc(__PRETTY_FUNCTION__,1);
		string func_name = __FUNCTION__;
#else
		string func_name = "";
#endif
		throw AngoraInvalidArgumentExceptionWithType<string>(func_name,constitutive_param_type,
			"(valid arguments are \"rel_permittivity\", \"rel_permeability\", \"electric_conductivity\", \"magnetic_conductivity\", and \"lorentz_delta_epsilon\")");
	}

  clock_t start;
  double duration;
  char *filename = const_cast<char*>(MaterialFileName.c_str());

  if(rank==0){
    cout << "Starting to read file: " << endl;
    cout << "\t" << filename << endl;
    start = clock();
  }

	int xExtentInCells,yExtentInCells,zExtentInCells;	//x, y and z extents of the

  // INSERT clean 3d MPI read code here

  int filename_length;
  int header_size = 3 * sizeof(int);
  int total_size[header_size];
  int NDIMS = 3;

  int my_rank, pool_size, order, file_name_length,
  array_of_distribs[NDIMS],
  array_of_dargs[NDIMS], count, read_buffer_size;  // array_of_gsizes[NDIMS], , array_of_psizes[NDIMS],

  double *read_buffer;
  BOOLEAN i_am_the_master = FALSE, input_error = FALSE,
  file_open_error = FALSE, file_write_error = FALSE, verbose = FALSE,
  my_read_error = FALSE, read_error = FALSE;
  char *file_name = NULL, message[BUFSIZ];

  /* MPI variables */

  MPI_Offset file_size;
  MPI_File fh;
  MPI_Status status;
  MPI_Datatype file_type;
  MPI_Aint file_type_extent;
  int file_type_size;
  int error_string_length;
  char error_string[BUFSIZ];
  extern int errno;


  int array_of_psizes[] = {nodes_x, nodes_y, nodes_z};//{0,0,0};
  MPI_Comm cartcomm;
  int periods[] = {0, 0,0};  /// should be 1??
  int pool_coords[NDIMS];
  const int reorder = 1; /// should be 0??

  // see: http://www.mcs.anl.gov/research/projects/mpi/usingmpi2/examples/moreio/subarray_c.htm
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  MPI_Cart_create(MPI_COMM_WORLD, NDIMS, array_of_psizes, periods, reorder, &cartcomm);
  MPI_Comm_rank(cartcomm, &my_rank);//MPI_Comm_rank(cartcomm, &my_rank);
  MPI_Cart_coords(cartcomm, my_rank, NDIMS, pool_coords);


  if (my_rank == MASTER_RANK) i_am_the_master = TRUE;
  //if (i_am_the_master) {
  //filename = const_cast<char*>(MaterialFileName.c_str());
  //}

  file_open_error = MPI_File_open(MPI_COMM_WORLD, filename,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  if (file_open_error != MPI_SUCCESS) {
    cout << "Error opening material input file " << MaterialFileName << "." << endl << endl;
    exit(-1);
  }

//MPI_Barrier(MPI_COMM_WORLD); // this might not be necessary??

  MPI_File_seek(fh, 0, MPI_SEEK_SET);
  MPI_File_read(fh, total_size, header_size, MPI_INT, &status);
  // TODO: try to cast the header values so you only have to do this once....
  // something liek this: MPI_Bcast(&total_size, 3, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

  xExtentInCells = total_size[0];
  yExtentInCells = total_size[1];
  zExtentInCells = total_size[2];

  //now, read through the file again, determine the material index for each point, and update the material indices in the main grid
  int material_offset;	//offset of the current material index beginning from the material index that was saved before the creation of new materials
  int material_index_eps,material_index_cond_e,material_index_mu,material_index_cond_h,material_index_Omega_p;		//material indices at a given point

  //calculate the coordinates of the back-left-lower corner of the region
  int xCornerPos = xPos;
  int yCornerPos = yPos;
  int zCornerPos = zPos;
  if ((anchor!="center")&&(anchor!="BLU")&&(anchor!="BLL")&&(anchor!="BRU")&&(anchor!="BRL")
                &&(anchor!="FLU")&&(anchor!="FLL")&&(anchor!="FRU")&&(anchor!="FRL"))
  {
    if (rank==0)
    {
      cout << "Invalid anchor point \"" << anchor << "\" for material input file " << MaterialFileName << " in node " << rank << endl << endl;
      exit(-1);
    }
  }

  if (anchor=="center")
  {
    xCornerPos = (int)floor(xPos-xExtentInCells/2.0);	//shift reference to the center
    yCornerPos = (int)floor(yPos-yExtentInCells/2.0);	//shift reference to the center
    zCornerPos = (int)floor(zPos-zExtentInCells/2.0);	//shift reference to the center
  }
  else if (anchor=="BLL")	//back-left-lower
  {
    //reference point is the back-left-lower corner by default
  }
  else if (anchor=="BLU")	//back-left-upper
  {
    zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
  }
  else if (anchor=="BRL")	//back-right-lower
  {
    yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
  }
  else if (anchor=="BRU")	//back-right-upper
  {
    yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
    zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
  }
  else if (anchor=="FLL")	//front-left-lower
  {
    xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
  }
  else if (anchor=="FLU")	//front-left-upper
  {
    xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
    zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
  }
  else if (anchor=="FRL")	//front-right-lower
  {
    xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
    yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
  }
  else if (anchor=="FRU")	//front-right-upper
  {
    xCornerPos = xPos - xExtentInCells;	//shift reference to the front corner
    yCornerPos = yPos - yExtentInCells;	//shift reference to the right corner
    zCornerPos = zPos - zExtentInCells;	//shift reference to the upper corner
  }

  //indices of the cell at the back-left-lower corner
  //these are simply equal to the coordinates of the back-left-lower corner plus one
  int xCornerCell = xCornerPos+1;
  int yCornerCell = yCornerPos+1;
  int zCornerCell = zCornerPos+1;


    /* Prepare for calling MPI_Type_create_darray */
  int array_of_gsizes[] = {xExtentInCells, yExtentInCells, zExtentInCells};

  for (int ii = 0; ii < NDIMS; ii++) {  // 3 dimensions
    array_of_distribs[ii] = MPI_DISTRIBUTE_BLOCK;
    array_of_dargs[ii]    = MPI_DISTRIBUTE_DFLT_DARG;
  }


  int iback3[3]  = {max(iback,xCornerCell), max(iback,xCornerCell), max(iback,xCornerCell)};
  int ifront3[3] = {min(ifront,xCornerCell + xExtentInCells-1),min(ifront+1,xCornerCell+xExtentInCells-1), min(ifront+1,xCornerCell+xExtentInCells-1)};
  int jleft3[3]  = {max(yCornerCell, jleft), max(yCornerCell, jleft), max(yCornerCell, jleft)};
  int jright3[3] = {min(yCornerCell + yExtentInCells-1, jright+1), min(yCornerCell + yExtentInCells-1, jright), min(yCornerCell + yExtentInCells-1, jright+1)};
  int klower3[3] = {max(zCornerCell, klower), max(zCornerCell, klower), max(zCornerCell, klower)};
  int kupper3[3] = {min(zCornerCell + zExtentInCells-1,kupper+1), min(zCornerCell + zExtentInCells-1,kupper+1), min(zCornerCell + zExtentInCells-1,kupper)};


  // Each parallel sub block overlaps its neigh by one. We check here if we can expand the size
  int subSize[] = {total_size[0] / array_of_psizes[0],
                   total_size[1] / array_of_psizes[1],
                   total_size[2] / array_of_psizes[2]};


   // always the same. dependent only on pool
   int start_indices[] = {min(max(iback-xCornerCell,0), xExtentInCells-1),
                          min(max(jleft-yCornerCell,0), yExtentInCells-1),
                          min(max(klower-zCornerCell,0),zExtentInCells-1)};

  // the read out is always one pixel bigger, and sometimes more (at edges) than size when divied evenly between pools
  int subSizeExtended[] = {min(max(maxSize(iback3, ifront3, NDIMS),3), total_size[0] - start_indices[0]),
                           min(max(maxSize(jleft3, jright3, NDIMS),3), total_size[1] - start_indices[1]),
                           min(max(maxSize(klower3, kupper3, NDIMS),3), total_size[2] - start_indices[2])};


// global indices pertain to i,j,k position in global coordinate system
  int start_indices_global[] = {min(iback3[0], min(iback3[1], iback3[2])) - 1,
                                min(jleft3[0], min(jleft3[1], jleft3[2])) - 1,
                                min(klower3[0], min(klower3[1], klower3[2])) - 1};

  MPI_Type_create_subarray(NDIMS, array_of_gsizes, subSizeExtended, start_indices,
                       MPI_ORDER_FORTRAN, MPI_DOUBLE, &file_type);
  MPI_Type_commit(&file_type);
  MPI_Type_size(file_type, &file_type_size);
  read_buffer_size = ( file_type_size ) / sizeof(double); // was int
  read_buffer = (double*) malloc(read_buffer_size * sizeof(double));


  MPI_File_set_view(fh, header_size, MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, read_buffer, read_buffer_size, MPI_DOUBLE, &status);


  // read in the global min and max
  double maxV = find_max(read_buffer, read_buffer_size); //read_buffer[0];
  double minV = find_min(read_buffer, read_buffer_size);//read_buffer[0];
  double globalMax;
  double globalMin;
  int recvcounts[] = {1};

  MPI_Reduce(&minV, &globalMin, 1, MPI_DOUBLE, MPI_MIN, 0,  MPI_COMM_WORLD);
  MPI_Reduce(&maxV, &globalMax, 1, MPI_DOUBLE, MPI_MAX, 0,  MPI_COMM_WORLD);

  MPI_Bcast(&globalMin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&globalMax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double param_lower_limit; //lower limit of the constitutive parameter
	if ((constitutive_param_type=="rel_permittivity")||(constitutive_param_type=="rel_permeability"))
	{
		param_lower_limit = 1;
	}
	else if ((constitutive_param_type=="electric_conductivity")||(constitutive_param_type=="magnetic_conductivity")||(constitutive_param_type=="lorentz_delta_epsilon"))
	{
		param_lower_limit = 0;
	}
	else
	{
		throw AngoraDeveloperException("Error in PlaceMaterialRegionFromFile: unknown constitutive parameter type");
	}

// TODO: GET and share all min/max to other processors - for now we will hard code it.
	MatType max_param = globalMax; // 0;	//maximum constitutive parameter value -- TODO: HARDCODE THE TRUE VALUES
	MatType min_param = globalMin;// 1e20;	//minimum constitutive parameter value -- TODO: HARDCODE THE TRUE VALUES

  //maximum number of different material types that can be extracted from the region
//	int max_num_of_materials = 1000; 	//pretty random, may have to find a more efficient way in the future
	//minimum difference in constitutive parameter between different materials
	MatType param_step = (max_param-min_param)/(max_number_of_new_materials-1);
	//add the new materials to the material list
	//before increasing the number of materials, save the current maximum material indices
	int material_index_saved_eps_x = eps_x.size()-1;
	int material_index_saved_eps_y = eps_y.size()-1;
	int material_index_saved_eps_z = eps_z.size()-1;
	int material_index_saved_mu_x = mu_x.size()-1;
	int material_index_saved_mu_y = mu_y.size()-1;
	int material_index_saved_mu_z = mu_z.size()-1;
	int material_index_saved_cond_e_x = cond_e_x.size()-1;
	int material_index_saved_cond_e_y = cond_e_y.size()-1;
	int material_index_saved_cond_e_z = cond_e_z.size()-1;
	int material_index_saved_cond_h_x = cond_h_x.size()-1;
	int material_index_saved_cond_h_y = cond_h_y.size()-1;
	int material_index_saved_cond_h_z = cond_h_z.size()-1;
	int material_index_saved_Omega_p_x = Omega_p_x.size()-1;
	int material_index_saved_Omega_p_y = Omega_p_y.size()-1;
	int material_index_saved_Omega_p_z = Omega_p_z.size()-1;


	//dummy material index
	Cmat NewMaterial;

	if ((constitutive_param_type=="rel_permittivity"))
	{//the relative permittivity of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_eps(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="rel_permeability"))
	{//the relative permeability of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_mu(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="electric_conductivity"))
	{//the electric conductivity (in S/m) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_cond_e(min_param+(i-1)*param_step);
		}
	}
	if ((constitutive_param_type=="magnetic_conductivity"))
	{//the magnetic conductivity (in Ohm/m) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_cond_h(min_param+(i-1)*param_step);
		}
	}
    if ((constitutive_param_type=="lorentz_delta_epsilon"))
	{//the magnetic conductivity (in Ohm/m) of the new material is min_param+(i-1)*param_step
		for (int i=1; i<=max_number_of_new_materials; i++)
		{
			NewMaterial.set_Omega_p(min_param+(i-1)*param_step);
		}
	}
  // So we have different limits based on if we are assigning the X, Y or Z Field components
  // we also have limits based on the x,y,z local position in the grid
  // eg. iback3 - is the x starting position for [X field, Y field and Z field]
  // at each we point, we need to check
/*
ofstream myfile1;
if(rank<20){
  printf("OPENING FILE...\n");
  myfile1.open("Debugging_MPI_tessst.txt",ios::app);
}*/

  int *temp;
  int i, j, k;

  MatType param_temp;	//constitutive parameter value that has been read
  for(int index = 0; index < read_buffer_size; index++ ){
    temp = ind2sub(subSizeExtended, index);
    i = start_indices_global[0] + temp[0] + 1; // 1 indexing
    j = start_indices_global[1] + temp[1] + 1; // 1 indexing
    k = start_indices_global[2] + temp[2] + 1; // 1 indexing
    param_temp = read_buffer[index];

    /*
    if(rank<20){
      myfile1 << my_rank << ": X: " << param_temp << " -> " << "[" <<
              i << ", " << j << ", " << k << "]""\n";
    }*/

    if (constitutive_param_type == "rel_permittivity"){
      if(isWithinBounds(i, iback3[0],  ifront3[0]) &&
         isWithinBounds(j, jleft3[0],  jright3[0]) &&
         isWithinBounds(k, klower3[0], kupper3[0])){
          // set x fields of rel-permittivity
/*
          if(rank==0){
            myfile1 << my_rank << ": X: " << param_temp << " -> " << "[" <<
                    i << ", " << j << ", " << k << "]""\n";
          }
*/
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_eps = material_index_saved_eps_x + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material at the center of the cell
          eps_x_indices(i,j,k) = material_index_eps;
          Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
          Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
      }
      if(isWithinBounds(i, iback3[1],  ifront3[1]) &&
         isWithinBounds(j, jleft3[1],  jright3[1]) &&
         isWithinBounds(k, klower3[1], kupper3[1])){
          // set y fields of rel-permittivity
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_eps = material_index_saved_eps_y + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material on the lower side of the cell
          eps_y_indices(i,j,k) = material_index_eps;
          Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
          Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));

      }
      if(isWithinBounds(i, iback3[2],  ifront3[2]) &&
         isWithinBounds(j, jleft3[2],  jright3[2]) &&
         isWithinBounds(k, klower3[2], kupper3[2])){
         // set z fields of rel-permittivity
         material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
         //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
         material_index_eps = material_index_saved_eps_z + (material_offset + 1);	//this is the absolute index in the current material list
                                           // +1 because of the range of material_offset above
         //place material on the left side of the cell
         eps_z_indices(i,j,k) = material_index_eps;
         Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
         Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
      }
    }
    else if (constitutive_param_type == "electric_conductivity"){
      if(isWithinBounds(i, iback3[0],  ifront3[0]) &&
         isWithinBounds(j, jleft3[0],  jright3[0]) &&
         isWithinBounds(k, klower3[0], kupper3[0])){
          // set x fields of electric_conductivity
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_cond_e = material_index_saved_cond_e_x + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material at the center of the cell
          cond_e_x_indices(i,j,k) = material_index_cond_e;
          Ca_X(i,j,k)=(1-dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
          Cb_X(i,j,k)=dt/eps_x(eps_x_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_x(cond_e_x_indices(i,j,k))/(2.0*eps_x(eps_x_indices(i,j,k))*epsilon_0));
      }
      if(isWithinBounds(i, iback3[1],  ifront3[1]) &&
         isWithinBounds(j, jleft3[1],  jright3[1]) &&
         isWithinBounds(k, klower3[1], kupper3[1])){
          // set y fields of electric_conductivity
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_cond_e = material_index_saved_cond_e_y + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material on the lower side of the cell
          cond_e_y_indices(i,j,k) = material_index_cond_e;
          Ca_Y(i,j,k)=(1-dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
          Cb_Y(i,j,k)=dt/eps_y(eps_y_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_y(cond_e_y_indices(i,j,k))/(2.0*eps_y(eps_y_indices(i,j,k))*epsilon_0));
      }
      if(isWithinBounds(i, iback3[2],  ifront3[2]) &&
         isWithinBounds(j, jleft3[2],  jright3[2]) &&
         isWithinBounds(k, klower3[2], kupper3[2])){
         // set z fields of electric_conductivity

         material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
         //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
         material_index_cond_e = material_index_saved_cond_e_z + (material_offset + 1);	//this is the absolute index in the current material list
                                           // +1 because of the range of material_offset above
         //place material on the left side of the cell
         cond_e_z_indices(i,j,k) = material_index_cond_e;
         Ca_Z(i,j,k)=(1-dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0))/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
         Cb_Z(i,j,k)=dt/eps_z(eps_z_indices(i,j,k))/epsilon_0/dx/(1+dt*cond_e_z(cond_e_z_indices(i,j,k))/(2.0*eps_z(eps_z_indices(i,j,k))*epsilon_0));
      }
    }
    else if (constitutive_param_type == "rel_permeability"){
      if(isWithinBounds(i, iback3[0],  ifront3[0]) &&
         isWithinBounds(j, jleft3[0],  jright3[0]) &&
         isWithinBounds(k, klower3[0], kupper3[0])){
          // set x fields of rel_permeability
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_mu = material_index_saved_mu_x + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material at the center of the cell
          mu_x_indices(i,j,k) = material_index_mu;
          Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
          Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
      }
      if(isWithinBounds(i, iback3[1],  ifront3[1]) &&
         isWithinBounds(j, jleft3[1],  jright3[1]) &&
         isWithinBounds(k, klower3[1], kupper3[1])){
          // set y fields of rel_permeability
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_mu = material_index_saved_mu_y + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material on the lower side of the cell
          mu_y_indices(i,j,k) = material_index_mu;
          Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
          Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
      }
      if(isWithinBounds(i, iback3[2],  ifront3[2]) &&
         isWithinBounds(j, jleft3[2],  jright3[2]) &&
         isWithinBounds(k, klower3[2], kupper3[2])){
        // set z fields of rel_permeability
         material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
         //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
         material_index_mu = material_index_saved_mu_z + (material_offset + 1);	//this is the absolute index in the current material list
                                           // +1 because of the range of material_offset above
         //place material on the left side of the cell
         mu_z_indices(i,j,k) = material_index_mu;
         Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
         Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
      }
    }
    else if (constitutive_param_type == "magnetic_conductivity"){
      if(isWithinBounds(i, iback3[0],  ifront3[0]) &&
         isWithinBounds(j, jleft3[0],  jright3[0]) &&
         isWithinBounds(k, klower3[0], kupper3[0])){
          // set x fields of magnetic_conductivity
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_cond_h = material_index_saved_cond_h_x + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material at the center of the cell
          cond_h_x_indices(i,j,k) = material_index_cond_h;
          Da_X(i,j,k)=(1-dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0))/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
          Db_X(i,j,k)=dt/mu_x(mu_x_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_x(cond_h_x_indices(i,j,k))/(2.0*mu_x(mu_x_indices(i,j,k))*mu_0));
      }
      if(isWithinBounds(i, iback3[1],  ifront3[1]) &&
         isWithinBounds(j, jleft3[1],  jright3[1]) &&
         isWithinBounds(k, klower3[1], kupper3[1])){
          // set y fields of magnetic_conductivity
          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_cond_h = material_index_saved_cond_h_y + (material_offset + 1);	//this is the absolute index in the current material list
                                            // +1 because of the range of material_offset above
          //place material on the lower side of the cell
          cond_h_y_indices(i,j,k) = material_index_cond_h;
          Da_Y(i,j,k)=(1-dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0))/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
          Db_Y(i,j,k)=dt/mu_y(mu_y_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_y(cond_h_y_indices(i,j,k))/(2.0*mu_y(mu_y_indices(i,j,k))*mu_0));
      }
      if(isWithinBounds(i, iback3[2],  ifront3[2]) &&
         isWithinBounds(j, jleft3[2],  jright3[2]) &&
         isWithinBounds(k, klower3[2], kupper3[2])){
             // set z fields of magnetic_conductivity

         material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
         //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
         material_index_cond_h = material_index_saved_cond_h_z + (material_offset + 1);	//this is the absolute index in the current material list
                                           // +1 because of the range of material_offset above
         //place material on the left side of the cell
         cond_h_z_indices(i,j,k) = material_index_cond_h;
         Da_Z(i,j,k)=(1-dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0))/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
         Db_Z(i,j,k)=dt/mu_z(mu_z_indices(i,j,k))/mu_0/dx/(1+dt*cond_h_z(cond_h_z_indices(i,j,k))/(2.0*mu_z(mu_z_indices(i,j,k))*mu_0));
      }
    }
    else if (constitutive_param_type == "lorentz_delta_epsilon"){
      if(isWithinBounds(i, iback3[0],  ifront3[0]) &&
         isWithinBounds(j, jleft3[0],  jright3[0]) &&
         isWithinBounds(k, klower3[0], kupper3[0])){
          // set x fields of lorentz_delta_epsilon

          material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
          //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
          material_index_Omega_p = material_index_saved_Omega_p_x + (material_offset + 1);	//this is the absolute index in the current material list
                          // +1 because of the range of material_offset above
          //place material at the center of the cell
          Omega_p_x_indices(i,j,k,0) = material_index_Omega_p;
          alpha_X(i,j,k,0) = (2-pow2(omega_p_x(omega_p_x_indices(i,j,k,0))*dt))/
                             (1+1/tau_p_x(tau_p_x_indices(i,j,k,0))*dt);
          xi_X(i,j,k,0) = (1/tau_p_x(tau_p_x_indices(i,j,k,0))*dt-1)/
                          (1/tau_p_x(tau_p_x_indices(i,j,k,0))*dt+1);
          gamma_X(i,j,k,0) = (epsilon_0*pow2(Omega_p_x(Omega_p_x_indices(i,j,k,0))*dt))/
                             (1/tau_p_x(tau_p_x_indices(i,j,k,0))*dt+1)/(2*dt);

          float gamma_p_sum(0.0);
          gamma_p_sum = 0.5*gamma_X(i,j,k,0)*(2*dt);

          Ca_X(i,j,k) = (2*epsilon_0*eps_x(eps_x_indices(i,j,k))-cond_e_x(cond_e_x_indices(i,j,k))*dt)/
              (2*epsilon_0*eps_x(eps_x_indices(i,j,k))+gamma_p_sum+cond_e_x(cond_e_x_indices(i,j,k))*dt);
          Cb_X(i,j,k) = 2*dt/dx/
              (2*epsilon_0*eps_x(eps_x_indices(i,j,k))+gamma_p_sum+cond_e_x(cond_e_x_indices(i,j,k))*dt);
          Cc_X(i,j,k) = gamma_p_sum/
              (2*epsilon_0*eps_x(eps_x_indices(i,j,k))+gamma_p_sum+cond_e_x(cond_e_x_indices(i,j,k))*dt);
      }
      if(isWithinBounds(i, iback3[1],  ifront3[1]) &&
         isWithinBounds(j, jleft3[1],  jright3[1]) &&
         isWithinBounds(k, klower3[1], kupper3[1])){
          // set y fields of lorentz_delta_epsilon

          Omega_p_y_indices(i,j,k,0) = material_index_Omega_p;
          alpha_Y(i,j,k,0) = (2-pow2(omega_p_y(omega_p_y_indices(i,j,k,0))*dt))/
                             (1+1/tau_p_y(tau_p_y_indices(i,j,k,0))*dt);
          xi_Y(i,j,k,0) = (1/tau_p_y(tau_p_y_indices(i,j,k,0))*dt-1)/
                          (1/tau_p_y(tau_p_y_indices(i,j,k,0))*dt+1);
          gamma_Y(i,j,k,0) = (epsilon_0*pow2(Omega_p_y(Omega_p_y_indices(i,j,k,0))*dt))/
                             (1/tau_p_y(tau_p_y_indices(i,j,k,0))*dt+1)/(2*dt);

          float gamma_p_sum(0.0);
          gamma_p_sum += 0.5*gamma_Y(i,j,k,0)*(2*dt);
          Ca_Y(i,j,k) = (2*epsilon_0*eps_y(eps_y_indices(i,j,k))-cond_e_y(cond_e_y_indices(i,j,k))*dt)/
              (2*epsilon_0*eps_y(eps_y_indices(i,j,k))+gamma_p_sum+cond_e_y(cond_e_y_indices(i,j,k))*dt);
          Cb_Y(i,j,k) = 2*dt/dx/
              (2*epsilon_0*eps_y(eps_y_indices(i,j,k))+gamma_p_sum+cond_e_y(cond_e_y_indices(i,j,k))*dt);
          Cc_Y(i,j,k) = gamma_p_sum/
              (2*epsilon_0*eps_y(eps_y_indices(i,j,k))+gamma_p_sum+cond_e_y(cond_e_y_indices(i,j,k))*dt);
      }
      if(isWithinBounds(i, iback3[2],  ifront3[2]) &&
         isWithinBounds(j, jleft3[2],  jright3[2]) &&
         isWithinBounds(k, klower3[2], kupper3[2])){
         // set z fields of lorentz_delta_epsilon

         material_offset = (int)((param_temp-min_param)/(max_param-min_param)*(max_number_of_new_materials-1));
         //material_offset is between 0 and (max_number_of_new_materials-1)  (both included)
         material_index_Omega_p = material_index_saved_Omega_p_z + (material_offset + 1);	//this is the absolute index in the current material list
                         // +1 because of the range of material_offset above
         //place material on the left side of the cell
         Omega_p_z_indices(i,j,k,0) = material_index_Omega_p;

         alpha_Z(i,j,k,0) = (2-pow2(omega_p_z(omega_p_z_indices(i,j,k,0))*dt))/
                            (1+1/tau_p_z(tau_p_z_indices(i,j,k,0))*dt);
         xi_Z(i,j,k,0) = (1/tau_p_z(tau_p_z_indices(i,j,k,0))*dt-1)/
                         (1/tau_p_z(tau_p_z_indices(i,j,k,0))*dt+1);
         gamma_Z(i,j,k,0) = (epsilon_0*pow2(Omega_p_z(Omega_p_z_indices(i,j,k,0))*dt))/
                            (1/tau_p_z(tau_p_z_indices(i,j,k,0))*dt+1)/(2*dt);
         float gamma_p_sum(0.0);
         gamma_p_sum += 0.5*gamma_Z(i,j,k,0)*(2*dt);
         Ca_Z(i,j,k) = (2*epsilon_0*eps_z(eps_z_indices(i,j,k))-cond_e_z(cond_e_z_indices(i,j,k))*dt)/
             (2*epsilon_0*eps_z(eps_z_indices(i,j,k))+gamma_p_sum+cond_e_z(cond_e_z_indices(i,j,k))*dt);
         Cb_Z(i,j,k) = 2*dt/dx/
             (2*epsilon_0*eps_z(eps_z_indices(i,j,k))+gamma_p_sum+cond_e_z(cond_e_z_indices(i,j,k))*dt);
         Cc_Z(i,j,k) = gamma_p_sum/
             (2*epsilon_0*eps_z(eps_z_indices(i,j,k))+gamma_p_sum+cond_e_z(cond_e_z_indices(i,j,k))*dt);
      }
    }
  }


  MPI_Barrier(MPI_COMM_WORLD);
  MPI_File_close(&fh);
  if(rank==0){
    duration = (clock() - start ) / (double) CLOCKS_PER_SEC;
    cout << "Finished reading file in " << duration << " seconds." << endl;
  }

/*
if(rank<20){
  printf("CLOSING FILE...\n");
  myfile1.close();
}*/
//  exit(-1);
}
#endif // MATFILE_H
