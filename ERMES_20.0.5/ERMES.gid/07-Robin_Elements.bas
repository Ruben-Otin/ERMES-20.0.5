*##############################################################################
*##############  Robin boundary condition surface elements  ###################
*##############################################################################
*if(strcmp(GenData(Problem_mode),"Electrostatic")==0)
// Electrostatic Robin boundary condition elements 
*set cond Robin_condition_Electrostatic *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if( strcmp(MatProp(0),cond(Robin_coefficients))==0 )
*set var MaterialNumber = MatNum
*endif
*end materials
GRC(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

*else
// Full-wave far-field elements
*set cond Far_field_condition_Full_wave *elems *canrepeat
*loop elems *onlyInCond
FF(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3));
*end elems

// Cold plasma far-field elements
*set cond Far_field_condition_Cold_plasma *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Plasma_wave_type))==0)
*set var MaterialNumber = MatNum
*endif
*end materials
PFF(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Full-wave Robin boundary condition elements
*set cond Robin_condition_Full_wave *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Robin_coefficients))==0)
*set var MaterialNumber = MatNum
*endif
*end materials	
GRC(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Rectangular waveguide port elements (TE10 mode)
*set cond RW_port_TE10_condition *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(RW_Port))==0)
*set var MaterialNumber = MatNum
*endif
*end materials
RWP(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

// Coaxial waveguide port elements (TEM mode)
*set cond Coaxial_TEM_condition *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Coaxial_Port))==0)
*set var MaterialNumber = MatNum
*endif
*end materials
COP(*GlobalNodes(1),*GlobalNodes(2),*GlobalNodes(3),*MaterialNumber);
*end elems

*endif 