*##############################################################################
*###########  Imported Robin boundary condition surface elements  #############
*##############################################################################
*Intformat "%i"
*set cond Robin_condition_Imported *elems *canrepeat
*loop elems *onlyInCond
*loop materials *NotUsed
*if(strcmp(MatProp(0),cond(Imported_Robin_surface))==0)
*set var RBCSurfID(int) = MatProp(IRBCS_ID,int)
*endif
*end materials	
*GlobalNodes(1) *GlobalNodes(2) *GlobalNodes(3) *RBCSurfID 
*end elems