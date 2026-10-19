*##############################################################################
*######################  Plane waves parameters  ##############################
*##############################################################################
*Realformat "%.18f"
*if( (strcmp(GenData(Problem_mode),"Full_wave")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*loop materials *NotUsed
*if(strcmp(MatProp(property),"Robin_coeffs_Plane_waves")==0)
*if(strcmp(MatProp(Multiple_waves),"On")==0)
*MatProp(Module_E_field,real) *MatProp(Phase_E_field,real) *MatProp(Polarization_X,real) *MatProp(Polarization_Y,real) *MatProp(Polarization_Z,real) *MatProp(Wave_vector_X,real) *MatProp(Wave_vector_Y,real) *MatProp(Wave_vector_Z,real)
*endif
*endif
*end materials
*endif