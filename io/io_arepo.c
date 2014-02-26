/*
 * Arepo I/O for Rockstar
 * Dylan Nelson (dnelson@cfa.harvard.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "hdf5.h" /* HDF5 required */
#include "io_arepo.h"
#include "io_util.h"
#include "../universal_constants.h"
#include "../check_syscalls.h"
#include "../config_vars.h"
#include "../config.h"
#include "../particle.h"

#define AREPO_NTYPES 6

void arepo_read_ids(char *filename, struct particle *p, int64_t num_p, int64_t nDM)
{
  hid_t HDF_FileID, HDF_GroupID, HDF_DatasetID, HDF_DataspaceID, HDF_Type;
  
  uint64_t *buffer;
  buffer = (uint64_t*) malloc(nDM * sizeof(uint64_t));
  
  HDF_Type = H5T_NATIVE_ULLONG; // 64

  // open file and group
  HDF_FileID  = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, "PartType1" );

  // open dataspace
  HDF_DatasetID   = H5Dopen(HDF_GroupID, "ParticleIDs");
  HDF_DataspaceID = H5Dget_space(HDF_DatasetID);
  H5Sselect_all(HDF_DataspaceID);
  
  // read and close
  H5Dread( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer );
  
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
  
  // write into p and free buffer
  for (int i=0; i < nDM; i++)
    p[num_p + i].id = buffer[i];
    
  free(buffer);
}

void arepo_read_posvel(char *filename, char *fieldName, struct particle *p, int64_t num_p, int64_t nDM)
{
  hid_t HDF_FileID, HDF_GroupID, HDF_DatasetID, HDF_DataspaceID, HDF_Type;
  
  int i;
  int offset = 0;
	
  if( strcmp(fieldName,"Velocities") == 0 )
    offset = 3;
    
  float *buffer;
  buffer = (float*) malloc(3 * nDM * sizeof(float));
  
  HDF_Type = H5T_NATIVE_FLOAT; //32

  HDF_FileID  = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID = H5Gopen( HDF_FileID, "PartType1" );

  // open dataspace and select all
  HDF_DatasetID   = H5Dopen(HDF_GroupID, fieldName);
  HDF_DataspaceID = H5Dget_space(HDF_DatasetID);
  H5Sselect_all(HDF_DataspaceID);
  
  // read and close
  H5Dread( HDF_DatasetID, HDF_Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer );
  
  H5Dclose( HDF_DatasetID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
  
  // write into p and free buffer
  for (i=0; i < nDM; i++)
  {
    p[num_p + i].pos[0+offset] = buffer[3*i+0];
    p[num_p + i].pos[1+offset] = buffer[3*i+1];
    p[num_p + i].pos[2+offset] = buffer[3*i+2];
  }
  
  free(buffer);
}

float HDFReadHeader_float(char *filename, char *objName)
{
  hid_t HDF_Type, HDF_FileID, HDF_GroupID, HDF_AttributeID, HDF_DataspaceID;
  HDF_Type = H5T_NATIVE_FLOAT;
  
  // open attribute
  HDF_FileID      = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID     = H5Gopen( HDF_FileID, "Header" );
  HDF_AttributeID = H5Aopen_name( HDF_GroupID, objName );
  
  // get space associated with attribute
  HDF_DataspaceID = H5Aget_space( HDF_AttributeID );
  
  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  hsize_t dimsize[ndims];
  
  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );
  
  // read data
  float data = 0.0;
  H5Aread( HDF_AttributeID, HDF_Type, &data );
  
  // close and return
  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
  
  return data;
}

void HDFReadHeader_array(char *filename, char *objName, int *npart, float *massTable)
{
  hid_t HDF_Type, HDF_FileID, HDF_GroupID, HDF_AttributeID, HDF_DataspaceID;
	
  if( npart != NULL )
    HDF_Type = H5T_NATIVE_UINT;
  else
    HDF_Type = H5T_NATIVE_FLOAT;
  
  // open attribute
  HDF_FileID      = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
  HDF_GroupID     = H5Gopen( HDF_FileID, "Header" );
  HDF_AttributeID = H5Aopen_name( HDF_GroupID, objName );

  // get space associated with attribute
  HDF_DataspaceID = H5Aget_space( HDF_AttributeID );
  
  int ndims = H5Sget_simple_extent_ndims( HDF_DataspaceID );
  hsize_t dimsize[ndims];
  
  H5Sget_simple_extent_dims( HDF_DataspaceID, dimsize, NULL );
  
  // read data
  if(npart != NULL)
    H5Aread( HDF_AttributeID, HDF_Type, npart );
  if(massTable != NULL)
    H5Aread( HDF_AttributeID, HDF_Type, massTable );
		
  // close and return
  H5Aclose( HDF_AttributeID );
  H5Sclose( HDF_DataspaceID );
  H5Gclose( HDF_GroupID );
  H5Fclose( HDF_FileID );
}

void arepo_rescale_particles(struct particle *p, int64_t p_start, int64_t nelems) {
  int64_t i, j;
  uint32_t id;
  double vel_rescale = sqrt(SCALE_NOW);
  if (LIGHTCONE) vel_rescale = 1;
	
  for (i=0; i<nelems; i++)
  {
    if (AREPO_ID_BYTES == 4) {
      exit(0); // better verify this whole routine for 32bit ids
      memcpy(&id, &(p[p_start+i].id), sizeof(uint32_t));
      p[p_start+i].id = id;
    }
    for (j=0; j<3; j++) {
      p[p_start+i].pos[j]   *= AREPO_LENGTH_CONVERSION;
      p[p_start+i].pos[j+3] *= vel_rescale;
    }
  }
	
}

void load_particles_arepo(char *filename, struct particle **p, int64_t *num_p)
{
  int npart[AREPO_NTYPES];
  float massTable[AREPO_NTYPES];
	
  // read header
  Ol = HDFReadHeader_float(filename,"OmegaLambda");
  Om = HDFReadHeader_float(filename,"Omega0");
  h0 = HDFReadHeader_float(filename,"HubbleParam");
 
  SCALE_NOW = HDFReadHeader_float(filename,"Time");
  BOX_SIZE = HDFReadHeader_float(filename,"BoxSize");
  BOX_SIZE *= AREPO_LENGTH_CONVERSION;  
  
  /*
  int npart_low[AREPO_NTYPES], npart_high[AREPO_NTYPES];
  HDFReadHeader_npart(filename,"NumPart_Total",&npart_low[0]);
  HDFReadHeader_npart(filename,"NumPart_Total_HighWord",&npart_high[0]);
  int64_t total_particles = ( ((int64_t)npart_high[AREPO_DM_PARTTYPE]) << 32 ) 
                             + (int64_t)npart_low[AREPO_DM_PARTTYPE];
  */
    
  HDFReadHeader_array(filename,"NumPart_ThisFile",&npart[0],NULL);
  HDFReadHeader_array(filename,"MassTable",NULL,&massTable[0]);
    
  TOTAL_PARTICLES = npart[AREPO_DM_PARTTYPE];
  PARTICLE_MASS   = massTable[AREPO_DM_PARTTYPE] * AREPO_MASS_CONVERSION;
  AVG_PARTICLE_SPACING = cbrt(PARTICLE_MASS / (Om*CRITICAL_DENSITY));
	
  if(RESCALE_PARTICLE_MASS) {
    printf("\n\nWARNING: RESCALE_PARTICLE_MASS probably don't want this (AREPO).\n\n");
    PARTICLE_MASS = Om*CRITICAL_DENSITY * pow(BOX_SIZE, 3) / TOTAL_PARTICLES;
  }
 
  printf("AREPO: filename:       %s\n",filename);
  printf("AREPO: BoxSize:        %g\n",BOX_SIZE);
  printf("AREPO: h0:             %g\n",h0);
  printf("AREPO: Time:           %g\n",SCALE_NOW);
  printf("AREPO: TotNumPartDM:   %" PRIu64 "\n",TOTAL_PARTICLES);
  printf("AREPO: dmPartMass:     %g\n",PARTICLE_MASS);
  printf("AREPO: avgPartSpacing: %g\n\n",AVG_PARTICLE_SPACING);
  
  if( !TOTAL_PARTICLES ) {
    printf("   SKIPPING FILE, PARTICLE COUNT ZERO.\n");
    return;
  }

  // re-allocate
  *p = (struct particle *) 
    check_realloc(*p, ((*num_p)+TOTAL_PARTICLES)*sizeof(struct particle), "Allocating particles.");

  // read IDs, pos, vel
  arepo_read_ids(filename, *p, *num_p, TOTAL_PARTICLES);
  arepo_read_posvel(filename, "Coordinates", *p, *num_p, TOTAL_PARTICLES);
  arepo_read_posvel(filename, "Velocities", *p, *num_p, TOTAL_PARTICLES);
  
  arepo_rescale_particles(*p, *num_p, TOTAL_PARTICLES);
  
  *num_p += TOTAL_PARTICLES;
}
