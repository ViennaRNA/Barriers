#ifndef BARRIERS_COMPRESS_H
#define BARRIERS_COMPRESS_H

void
ini_pack_em(barrier_options opt);


char *
pack_em(const char *string);


char *
unpack_em(const char *packed);


#endif
