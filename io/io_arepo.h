#ifndef _IO_AREPO_H_
#define _IO_AREPO_H_
#include <stdint.h>
#include "../particle.h"

void load_particles_arepo(char *filename, struct particle **p, int64_t *num_p);

#endif /* _IO_AREPO_H_ */
