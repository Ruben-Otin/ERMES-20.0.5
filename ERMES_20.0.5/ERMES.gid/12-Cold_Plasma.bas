*##############################################################################
*######################  Cold plasma parameters  ##############################
*##############################################################################
*Intformat "%i"
*Realformat "%.18f"
*loop materials 
*if(strcmp(MatProp(property),"cold_plasma")==0)
*##########  Material Id  #####################################################
*MatNum
*##########  Plasma geometry extrapolation ####################################
*if(strcmp(MatProp(Geometry_extrapolation),"Flat_1D")==0)
0
*elseif(strcmp(MatProp(Geometry_extrapolation),"Curved")==0)
1
*elseif(strcmp(MatProp(Geometry_extrapolation),"Full_3D")==0)
2
*elseif(strcmp(MatProp(Geometry_extrapolation),"Tensor")==0)
3
*endif
*##########  Length scales  ###################################################
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
*##########  Poloidal curvatures  #############################################
*operation(pLengthFactor*MatProp(Poloidal_ne_curvature,real))
*operation(pLengthFactor*MatProp(Poloidal_Be_curvature,real))
*##########  1D line first point  #############################################
*operation(pLengthFactor*MatProp(L1_00_X,real)) *operation(pLengthFactor*MatProp(L1_00_Y,real)) *operation(pLengthFactor*MatProp(L1_00_Z,real))
*##########  Tokamak center  ##################################################
*operation(pLengthFactor*MatProp(TK_00_X,real)) *operation(pLengthFactor*MatProp(TK_00_Y,real)) *operation(pLengthFactor*MatProp(TK_00_Z,real))
*##########  Tokamak axis  ####################################################
*MatProp(TK_axis_X,real) *MatProp(TK_axis_Y,real) *MatProp(TK_axis_Z,real)
*##########  Electron density file path  ######################################
*MatProp(ne_file_path)
*##########  External B field file path  ######################################
*MatProp(Be_file_path)
*##########  Species relative densities  ######################################
*MatProp(H1_1+_/_ne)
*MatProp(H2_1+_/_ne)
*MatProp(H3_1+_/_ne)
*MatProp(He3_1+_/_ne)
*MatProp(He3_2+_/_ne)
*MatProp(He4_1+_/_ne)
*MatProp(He4_2+_/_ne)
*MatProp(Be9_1+_/_ne)
*MatProp(Be9_2+_/_ne)
*MatProp(Be9_3+_/_ne)
*MatProp(Be9_4+_/_ne)
*MatProp(Ne20_10+_/_ne)
*MatProp(Ni60_20+_/_ne)
*MatProp(W184_25+_/_ne)
*##########  Matrix format  ###################################################
*if(strcmp(MatProp(Matrix_format),"Full_matrix")==0)
0
*elseif(strcmp(MatProp(Matrix_format),"Herm-Full")==0)
1
*elseif(strcmp(MatProp(Matrix_format),"Herm-Symm")==0)
2
*endif
*##########  Kij tolerance  ###################################################
*if(strcmp(MatProp(Kij_tolerance),"custom")==0)
*MatProp(custom_Kij_tol)
*elseif(strcmp(MatProp(Kij_tolerance),"Off")==0)
0.0
*else
*MatProp(Kij_tolerance)
*endif
*##########  E parallel tolerance  ############################################
*if(strcmp(MatProp(E_par_tolerance),"custom")==0)
*MatProp(custom_E_par_tol)
*elseif(strcmp(MatProp(E_par_tolerance),"Off")==0)
0.0
*else
*MatProp(E_par_tolerance)
*endif
*##########  Collisional damping  #############################################
*if(strcmp(MatProp(Damping_mode),"No_damping")==0)
0
*elseif(strcmp(MatProp(Damping_mode),"Background_ctc")==0)
1
*elseif(strcmp(MatProp(Damping_mode),"Background_den")==0)
2
*elseif(strcmp(MatProp(Damping_mode),"ComplexFreq_all")==0)
3
*elseif(strcmp(MatProp(Damping_mode),"ComplexFreq_ele")==0)
4
*endif
*MatProp(Damping__value)
*endif
*end materials