*##############################################################################
*##########################  Material properties ##############################
*##############################################################################
*Intformat "%i"
*Realformat "%.18f"
*#######################  IHL material properties  ############################
*loop materials 
*if(strcmp(MatProp(property),"ihl_material")==0)
// IHL material 
PROPERTIES[*MatNum].REAL_ELECTRIC_CONDUCTIVITY = *MatProp(Electric_conductivity_real,real);
PROPERTIES[*MatNum].IMAG_ELECTRIC_CONDUCTIVITY = *MatProp(Electric_conductivity_img,real);
PROPERTIES[*MatNum].REAL_ELECTRIC_PERMITTIVITY = *MatProp(Electric_permittivity_real,real);
PROPERTIES[*MatNum].IMAG_ELECTRIC_PERMITTIVITY = *MatProp(Electric_permittivity_img,real);
PROPERTIES[*MatNum].REAL_MAGNETIC_PERMEABILITY = *MatProp(Magnetic_permeability_real,real);
PROPERTIES[*MatNum].IMAG_MAGNETIC_PERMEABILITY = *MatProp(Magnetic_permeability_img,real);

*endif
*if(strcmp(MatProp(property),"cold_plasma")==0)
// Cold plasma background material
PROPERTIES[*MatNum].REAL_ELECTRIC_CONDUCTIVITY = 0.0;
PROPERTIES[*MatNum].IMAG_ELECTRIC_CONDUCTIVITY = 0.0;
PROPERTIES[*MatNum].REAL_ELECTRIC_PERMITTIVITY = 1.0;
PROPERTIES[*MatNum].IMAG_ELECTRIC_PERMITTIVITY = 0.0;
PROPERTIES[*MatNum].REAL_MAGNETIC_PERMEABILITY = 1.0;
PROPERTIES[*MatNum].IMAG_MAGNETIC_PERMEABILITY = 0.0;

*endif
*end materials
*###########  Robin conditions and J sources properties  ######################
*if(strcmp(GenData(Length_units),"custom")==0)
*set var pLengthFactor(real) = GenData(custom_L_factor,real)
*elseif(strcmp(GenData(Length_units),"m")==0)
*set var pLengthFactor(real) = 1
*elseif(strcmp(GenData(Length_units),"dm")==0)
*set var pLengthFactor(real) = 0.1
*elseif(strcmp(GenData(Length_units),"cm")==0)
*set var pLengthFactor(real) = 0.01
*elseif(strcmp(GenData(Length_units),"mm")==0)
*set var pLengthFactor(real) = 0.001
*endif
*loop materials *NotUsed
*if(strcmp(MatProp(property),"RWPort_coeffs")==0)
// Rectangular waveguide port
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(RW_Port_ID,int),0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*MatProp(Mod_RW_Ein,real),*MatProp(Phas_RW_Ein,real)];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(RW_00_X,real)),*operation(pLengthFactor*MatProp(RW_00_Y,real)),*operation(pLengthFactor*MatProp(RW_00_Z,real)),0,0,0];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*operation(pLengthFactor*MatProp(Height_X,real)),*operation(pLengthFactor*MatProp(Height_Y,real)),*operation(pLengthFactor*MatProp(Height_Z,real)),*operation(pLengthFactor*MatProp(Width_X,real)),*operation(pLengthFactor*MatProp(Width_Y,real)),*operation(pLengthFactor*MatProp(Width_Z,real))];

*endif
*if(strcmp(MatProp(property),"CoaxialPort_coeffs")==0)
// Coaxial waveguide port
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(CC_Port_ID,int),0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*MatProp(Mod_CC_Ein,real),*MatProp(Phas_CC_Ein,real)];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(CC_00_X,real)),*operation(pLengthFactor*MatProp(CC_00_Y,real)),*operation(pLengthFactor*MatProp(CC_00_Z,real)),*operation(pLengthFactor*MatProp(Inner_D,real)),*operation(pLengthFactor*MatProp(Exter_D,real)),0];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Epr_coax,real),*MatProp(Mu_coax,real),0,0,0,0];

*endif
*if(strcmp(MatProp(property),"Robin_coeffs_Plane_waves")==0)
// Plane wave Robin coefficients 
*set var MultWaves(int) = 1
*if(strcmp(MatProp(Multiple_waves),"On")==0)
*set var MultWaves(int) = 2
*endif
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(IBC_conductivity_real,real),*MatProp(IBC_conductivity_img,real)];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*MatProp(IBC_permittivity_real,real),*MatProp(IBC_permittivity_img,real)];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*MatProp(IBC_permeability_real,real),*MatProp(IBC_permeability_img,real),0,*MatProp(Module_E_field,real),*MatProp(Phase_E_field,real),*MultWaves];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Polarization_X,real),*MatProp(Polarization_Y,real),*MatProp(Polarization_Z,real),*MatProp(Wave_vector_X,real),*MatProp(Wave_vector_Y,real),*MatProp(Wave_vector_Z,real)];

*endif
*if(strcmp(MatProp(property),"Robin_coeffs_Gauss_beam")==0)
// Gaussian beam Robin coefficients 
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(GB_module_E_per,real),*MatProp(GB_phase_E_per,real)];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*operation(pLengthFactor*MatProp(GB_radius_par,real)),*operation(pLengthFactor*MatProp(GB_radius_per,real))];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(GB_00_X,real)),*operation(pLengthFactor*MatProp(GB_00_Y,real)),*operation(pLengthFactor*MatProp(GB_00_Z,real)),*MatProp(GB_module_E_par,real),*MatProp(GB_phase_E_par,real),3];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(GB_polarization_X,real),*MatProp(GB_polarization_Y,real),*MatProp(GB_polarization_Z,real),*MatProp(GB_wave_vector_X,real),*MatProp(GB_wave_vector_Y,real),*MatProp(GB_wave_vector_Z,real)];

*endif
*if(strcmp(MatProp(property),"Robin_coeffs_Quasiestatic")==0)
// Quasiestatic Robin coefficients 
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(Robin_coeff_real,real),*MatProp(Robin_coeff_img,real)];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*MatProp(Module_Flux,real),*MatProp(Phase_Flux,real)];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [0,0,0,0,0,0];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [0,0,0,0,0,0];

*endif
*if(strcmp(MatProp(property),"Robin_coeffs_Electrostatic")==0)
// Electrostatic Robin coefficients 
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(Robin_coeff,real),0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*MatProp(Flux,real),0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [0,0,0,0,0,0];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [0,0,0,0,0,0];

*endif
*if(strcmp(MatProp(property),"FDJcartesian")==0)
// Cartesian J source
PROPERTIES[*MatNum].COMPLEX_IBC                = [0,0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [0,0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [0,0,0,0,0,0];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Jx,real),*MatProp(Phx,real),*MatProp(Jy,real),*MatProp(Phy,real),*MatProp(Jz,real),*MatProp(Phz,real)];

*endif
*if(strcmp(MatProp(property),"FDJaxisymm")==0)
// Axisymmetric J source
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(Ja,real),*MatProp(Pha,real)];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [0,0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(AX_00_X,real)),*operation(pLengthFactor*MatProp(AX_00_Y,real)),*operation(pLengthFactor*MatProp(AX_00_Z,real)),*MatProp(AX_axis_X,real),*MatProp(AX_axis_Y,real),*MatProp(AX_axis_Z,real)];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [0,0,0,0,0,0];

*endif
*if(strcmp(MatProp(property),"FDJcombine")==0)
// Combined cartesian and axisymmetric J source
PROPERTIES[*MatNum].COMPLEX_IBC                = [*MatProp(Ja,real),*MatProp(Pha,real)];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [0,0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(AX_00_X,real)),*operation(pLengthFactor*MatProp(AX_00_Y,real)),*operation(pLengthFactor*MatProp(AX_00_Z,real)),*MatProp(AX_axis_X,real),*MatProp(AX_axis_Y,real),*MatProp(AX_axis_Z,real)];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Jx,real),*MatProp(Phx,real),*MatProp(Jy,real),*MatProp(Phy,real),*MatProp(Jz,real),*MatProp(Phz,real)];

*endif
*if(strcmp(MatProp(property),"FDJhelicoidal")==0)
// Helicoidal J source
PROPERTIES[*MatNum].COMPLEX_IBC                = [0,0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*operation(pLengthFactor*MatProp(HL_R_,real)),0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(HL_00_X,real)),*operation(pLengthFactor*MatProp(HL_00_Y,real)),*operation(pLengthFactor*MatProp(HL_00_Z,real)),*MatProp(HL_axis_X,real),*MatProp(HL_axis_Y,real),*MatProp(HL_axis_Z,real)];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Jpol,real),*MatProp(Php,real),*MatProp(Jtor,real),*MatProp(Pht,real),0,0];

*endif
*if(strcmp(MatProp(property),"FDJplasmam")==0)
// Plasma mode J source
PROPERTIES[*MatNum].COMPLEX_IBC                = [0,0];
PROPERTIES[*MatNum].COMPLEX_IBC_2o             = [*operation(pLengthFactor*MatProp(PM_R_,real)),0];
PROPERTIES[*MatNum].COMPLEX_NEUMANN_FLOW       = [*operation(pLengthFactor*MatProp(PM_00_X,real)),*operation(pLengthFactor*MatProp(PM_00_Y,real)),*operation(pLengthFactor*MatProp(PM_00_Z,real)),*MatProp(PM_axis_X,real),*MatProp(PM_axis_Y,real),*MatProp(PM_axis_Z,real)];
PROPERTIES[*MatNum].SINUSOIDAL_SURFACE_CURRENT = [*MatProp(Jpol,real),*MatProp(Php,real),*MatProp(Jtor,real),*MatProp(Pht,real),*MatProp(PM_M,int),1];

*endif
*end materials
*#######################  Surfaces and volumes IDs  ###########################
// Surfaces and volumes IDs
*loop materials *NotUsed
*if(strcmp(MatProp(property),"Plasma_wave_type")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [*MatProp(Plasma_Wave_Type_ID,int),0];
*endif
*if(strcmp(MatProp(property),"PBCSurface")==0)
*if(strcmp(MatProp(PBCSurface_ID),"11")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [11,*GenData(PBC_FB_pha_dif,real)];
*elseif(strcmp(MatProp(PBCSurface_ID),"12")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [12,*GenData(PBC_FB_pha_dif,real)];
*elseif(strcmp(MatProp(PBCSurface_ID),"21")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [21,*GenData(PBC_RL_pha_dif,real)];
*elseif(strcmp(MatProp(PBCSurface_ID),"22")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [22,*GenData(PBC_RL_pha_dif,real)];
*endif
*endif
*if(strcmp(MatProp(property),"SurfaceToInt")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [*MatProp(Surface_ID,int),0];
*endif
*if(strcmp(MatProp(property),"VolumeToInt")==0)
PROPERTIES[*MatNum].COMPLEX_IBC = [*MatProp(Volume_ID,int),0];
*endif
*end materials