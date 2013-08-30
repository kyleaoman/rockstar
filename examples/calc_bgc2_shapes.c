#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include "../config_vars.h"
#include "../config.h"
#include "../rockstar.h"
#include "../io/meta_io.h"
#include "../io/io_bgc2.h"
#include "../check_syscalls.h"
#include "../io/bgc2.h"
#include "load_bgc2.h"
#include "../groupies.c"
#include "../halo.h"
#include "../potential.h"

float calc_distance2(float *p1, float *p2, float box_size) {
  double ds=0, dx;
  for (int64_t i=0; i<3; i++) {
    dx = fabs(p1[i]-p2[i]);
    if (dx > box_size/2.0) dx = box_size - dx;
    ds += dx*dx;
  }
  return ds;
}

int main(int argc, char **argv) {
  int64_t i,j,k, max_p=0;
  struct bgc2_header hdr;

  if (argc < 4) {
    printf("Usage: %s mass_definition radius_fraction file1.bgc2 ...\n", argv[0]);
    printf("E.g.: %s vir 0.5 file1.bgc2 #calculates shape within 0.5 Rvir.\n", argv[0]);
    exit(1);
  }

  float rfrac = atof(argv[2]);

  printf("#HaloID HaloParent Mass Radius Vmax X Y Z VX VY VZ Mass_within_%gxR%s A[x] A[y] A[z] b_to_a c_to_a\n", rfrac, argv[1]);
  printf("#Units: positions in Mpc/h (comoving)\n");
  printf("#Units: radii and other lengths in kpc/h (comoving)\n");
  printf("#Units: velocities in km/s (non-comoving)\n");
  printf("#HaloParent: the ID of a larger halo containing this halo, if any; -1 otherwise\n");
  printf("#Shapes calculated within %g x R%s\n", rfrac, argv[1]);
  printf("#b_to_a, c_to_a: Ratio of second and third largest shape ellipsoid axes (B and C) to largest shape ellipsoid axis (A) (dimensionless).\n#  Shapes are determined by the method in Allgood et al. (2006).\n#A[x],A[y],A[z]: Largest shape ellipsoid axis (kpc/h comoving)\n#All particles are considered (bound and unbound) so these measurements\n#  may differ from the Rockstar catalogs, especially for subhalos.\n");
  for (i=3; i<argc; i++) {
    num_groups = num_parts = 0;
    load_bgc2(argv[i], &hdr, &grps, &num_groups, &parts, &num_parts);
    if (i==3) {
      printf("#Box size: %g Mpc/h\n", hdr.box_size);
      printf("#Particle Mass: %g Msun/h\n", hdr.part_mass);
      printf("#h: %g; Om: %g; Ol: %g\n", hdr.Hubble0, hdr.Omega0, hdr.OmegaLambda);
      printf("#Redshift: %g\n", hdr.redshift);
    }
    int64_t p_start = 0;
    SCALE_NOW = 1.0/(hdr.redshift+1.0);
    h0 = hdr.Hubble0;
    Om = hdr.Omega0;
    Ol = hdr.OmegaLambda;
    PARTICLE_MASS = hdr.part_mass;
    MASS_DEFINITION = argv[1];
    calc_mass_definition();
    double dens_thresh = particle_thresh_dens[0]*(4.0*M_PI/3.0)*hdr.part_mass;
    for (j=0; j<hdr.ngroups; j++) {
      struct halo h;
      memset(&h, 0, sizeof(struct halo));
      check_realloc_var(po, sizeof(struct potential), max_p, grps[j].npart);
      for (k=0; k<grps[j].npart; k++) {
	memcpy(po[k].pos, parts[p_start+k].pos, sizeof(float)*3);
	memcpy(po[k].pos+3, parts[p_start+k].vel, sizeof(float)*3);
	po[k].pe = po[k].ke = 0;
	po[k].flags = 0;
	po[k].r2 = calc_distance2(grps[j].pos, po[k].pos, hdr.box_size);
      }
      memcpy(h.pos, grps[j].pos, sizeof(float)*3);
      memcpy(h.pos+3, grps[j].vel, sizeof(float)*3);
      qsort(po, grps[j].npart, sizeof(struct potential), dist_compare);

      h.r = 0;
      for (k=grps[j].npart-1; k>=0; k--) {
	double dens = (k+1)*PARTICLE_MASS;
	double r = sqrt(po[k].r2);
	if (dens/(r*r*r) >= dens_thresh) break;
      }
      if (k>=0) h.r = po[k].r2;
      for (k=0; k<grps[j].npart; k++)
	if (po[k].r2 > rfrac*rfrac*h.r) break;

      h.r = rfrac*sqrt(h.r)*1e3;
      h.m = (h.r) ? (k+1)*PARTICLE_MASS : 0;
      calc_shape(&h, k, 0);
      printf("%"PRId64" %"PRId64" %e %f %f %f %f %f %f %f %f %e %f %f %f %f %f\n",
	     grps[j].id, grps[j].parent_id, grps[j].mass, 1e3*grps[j].radius, grps[j].vmax, grps[j].pos[0], grps[j].pos[1], grps[j].pos[2], grps[j].vel[0], grps[j].vel[1], grps[j].vel[2], k*PARTICLE_MASS, h.A[0], h.A[1], h.A[2], h.b_to_a, h.c_to_a);
      p_start += grps[j].npart;
    }
  }
  return 0;
}
