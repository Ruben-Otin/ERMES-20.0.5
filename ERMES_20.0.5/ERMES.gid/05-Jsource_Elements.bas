*##############################################################################
*####################  Volumetric J source elements  ##########################
*##############################################################################
// Volumetric J source elements
*if( (strcmp(GenData(Problem_mode),"Full_wave")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*set cond FD_current_density *elems *canrepeat
*set elems(tetrahedra)
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Source))==0)
*set var MaterialNumber = MatNum
*endif
*end materials
JE(*ElemsConec(1),*ElemsConec(2),*ElemsConec(3),*ElemsConec(4),*MaterialNumber);
*end elems
*endif