*##############################################################################
*###################  Projection and integral elements  #######################
*##############################################################################
*if(strcmp(GenData(Problem_mode),"Electrostatic")==0)
*if( (strcmp(GenData(Results_mode),"GP_1")==0) || (strcmp(GenData(Results_mode),"GP_4")==0) || (strcmp(GenData(Results_mode),"GP_10")==0) )
// Electrostatic fields integration over a surface of non-smoothed fields  
*set cond Field_surface_integral *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Surface_to_integer))==0)
*set var MaterialNumber=MatNum
*endif
*end materials
PVIE(*ElemsConec(1),*ElemsConec(2),*ElemsConec(3),*ElemsConec(4),0);
*end elems

// Initiate non-smooth integration surface counter
ProblemType = Ini_Surf_Counter;

*endif
// Electrostatic fields integration over a surface
*set cond Field_surface_integral *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Surface_to_integer))==0)
*set var MaterialNumber=MatNum
*endif
*end materials	
PSIE(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Electrostatic fields integration over a volume
*set cond Field_volume_integral *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Volume_to_integer))==0)
*set var MaterialNumber=MatNum
*endif
*end materials
PVIE(*ElemsConec(1),*ElemsConec(2),*ElemsConec(3),*ElemsConec(4),*MaterialNumber);
*end elems

*else
// Projection over rectangular waveguide TE10 mode
*set cond Projection_RWTE10 *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(RWPort))==0)
*set var MaterialNumber=MatNum
*endif
*end materials	
PRWP(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Projection over coaxial TEM mode
*set cond Projection_CoaxialTEM *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(CoaxialPort))==0)
*set var MaterialNumber=MatNum
*endif
*end materials	
PCOP(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Fields integration over a surface
*set cond Field_surface_integral *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Surface_to_integer))==0)
*set var MaterialNumber=MatNum
*endif
*end materials	
PSIE(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Fields integration over a volume
*set cond Field_volume_integral *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Volume_to_integer))==0)
*set var MaterialNumber=MatNum
*endif
*end materials	
PVIE(*ElemsConec(1),*ElemsConec(2),*ElemsConec(3),*ElemsConec(4),*MaterialNumber);
*end elems

*endif