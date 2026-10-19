*##############################################################################
*##############  Contact elements for discontinuity surfaces  #################
*##############################################################################
// Contact elements (discontinuity surfaces)
*if( (strcmp(GenData(Problem_mode),"Full_wave")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*if(strcmp(GenData(Contact_detect),"Manual")==0)
*set elems(Prisma)
*loop elems
CE = n([*ElemsConec(1),*ElemsConec(2),*ElemsConec(3),*ElemsConec(4),*ElemsConec(5),*ElemsConec(6)]);
*end elems
*endif
*endif