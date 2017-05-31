//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
//
// Abort the program when an error occurs.
//________________________________________________________
//A.Z. - 02/08 => Based on my C++ equivalent library.
//=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|
#ifndef ABORT_H
#define ABORT_H

#include <stdio.h>
#include <stdlib.h>


#define Abort( error ){ fprintf(stderr,"!!  File: %s\n!!  Line: %d\n!!  %s\n\a!!  ABORTING ...\n",__FILE__,__LINE__,error);exit(1);}

#endif /* ABORT_H */
