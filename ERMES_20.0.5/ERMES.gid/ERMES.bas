*##############################################################################
*#####################  ERMES problem definition ##############################
*##############################################################################
*Intformat "%i"
*Realformat "%.18f"
*#######################  Problem settings  ###################################
// Problem settings
*if(strcmp(GenData(Problem_mode),"Full_wave")==0)
ProblemType = Full_Wave;
*if(strcmp(GenData(Symmetry),"3D")==0)
ProblemType = Exyz_3D;
*elseif(strcmp(GenData(Symmetry),"XY")==0)
ProblemType = Exy_3D;
*elseif(strcmp(GenData(Symmetry),"ZZ")==0)
ProblemType = Ez_3D;
*elseif(strcmp(GenData(Symmetry),"AY")==0)
ProblemType = Ea_3D;
*endif
*elseif(strcmp(GenData(Problem_mode),"Cold_plasma")==0)
ProblemType = Cold_Plasma;
*elseif(strcmp(GenData(Problem_mode),"Electrostatic")==0)
ProblemType = Electrostatic;
*endif
*#######################  GiD geometric tolerance  ############################
*if(strcmp(GenData(Geo_error_tol),"0.0")==0)
ProblemType = GeoTol_0;
*elseif(strcmp(GenData(Geo_error_tol),"1e-1")==0)
ProblemType = GeoTol_1;
*elseif(strcmp(GenData(Geo_error_tol),"1e-2")==0)
ProblemType = GeoTol_1em2;
*elseif(strcmp(GenData(Geo_error_tol),"1e-3")==0)
ProblemType = GeoTol_1em3;
*elseif(strcmp(GenData(Geo_error_tol),"1e-4")==0)
ProblemType = GeoTol_1em4;
*elseif(strcmp(GenData(Geo_error_tol),"1e-5")==0)
ProblemType = GeoTol_1em5;
*elseif(strcmp(GenData(Geo_error_tol),"1e-6")==0)
ProblemType = GeoTol_1em6;
*elseif(strcmp(GenData(Geo_error_tol),"1e-7")==0)
ProblemType = GeoTol_1em7;
*elseif(strcmp(GenData(Geo_error_tol),"1e-8")==0)
ProblemType = GeoTol_1em8;
*elseif(strcmp(GenData(Geo_error_tol),"1e-9")==0)
ProblemType = GeoTol_1em9;
*elseif(strcmp(GenData(Geo_error_tol),"1e-10")==0)
ProblemType = GeoTol_1em10;
*elseif(strcmp(GenData(Geo_error_tol),"1e-11")==0)
ProblemType = GeoTol_1em11;
*elseif(strcmp(GenData(Geo_error_tol),"1e-12")==0)
ProblemType = GeoTol_1em12;
*endif
*##################  Periodic boundary condition type  ########################
*if(strcmp(GenData(PBC_symmetry),"Cyclic")==0)
ProblemType = PBC_Cyclic;
*elseif(strcmp(GenData(PBC_symmetry),"Periodic")==0)
ProblemType = PBC_Periodic;
*endif
*##########################  Solution mode  ###################################
*if(strcmp(GenData(Sweep_frequency_mode),"1")==0)
ProblemType = Release_Mode;
*elseif(strcmp(GenData(Solving_mode),"Release")==0)
ProblemType = Release_Mode;
*elseif(strcmp(GenData(Solving_mode),"Debug")==0)
ProblemType = Debug_Mode;
*elseif(strcmp(GenData(Solving_mode),"Read")==0)
ProblemType = Read_Mode;
*endif
*######################  Element type and order  ##############################
ProblemType = *GenData(_Element_type);
*####################  E-fields or AV-potentials  #############################
*if(strcmp(GenData(Potentials),"On")==0)
ProblemType = Potentials_On;
*elseif(strcmp(GenData(Potentials),"Off")==0)
ProblemType = Potentials_Off;
*endif
*###############  AV continuity on contact surfaces  ##########################
*if(strcmp(GenData(AV_continuity),"On")==0)
ProblemType = AV_Continuity_On;
*elseif(strcmp(GenData(AV_continuity),"Off")==0)
ProblemType = AV_Continuity_Off;
*endif
*####################  Visualization points  ##################################
*if(strcmp(GenData(Results_mode),"Nodes")==0)
ProblemType = Results_On_Nodes;
*elseif(strcmp(GenData(Results_mode),"GP_1")==0)
ProblemType = Results_On_1GP;
*elseif(strcmp(GenData(Results_mode),"GP_4")==0)
ProblemType = Results_On_4GP;
*endif
*#################  LL2P smoothing on results  ################################
*if(strcmp(GenData(LL2P_smooth),"On")==0)
ProblemType = LL2P_Smoothing_On;
*else
ProblemType = LL2P_Smoothing_Off;
*endif
*################  Type of surface normal averaging  ##########################
*if(strcmp(GenData(Normals_type),"Geometric_average")==0)
ProblemType = Geometric_Average_Normals;
*elseif(strcmp(GenData(Normals_type),"Area_weighted")==0)
ProblemType = Area_Weighted_Normals;
*endif
*#######################  RME stabilization  #################################
*if(strcmp(GenData(RME_stabilize),"On")==0)
ProblemType = RMED_Stab_On;
*elseif(strcmp(GenData(RME_stabilize),"Off")==0)
ProblemType = RMED_Stab_Off;
*endif
*#######################  EDG stabilization  ##################################
*if(strcmp(GenData(EDG_stabilize),"On")==0)
ProblemType = EDGE_Stab_On;
*elseif(strcmp(GenData(EDG_stabilize),"Off")==0)
ProblemType = EDGE_Stab_Off;
*endif
*#######################  LL2P stabilization  #################################
*if(strcmp(GenData(LL2P_stabilize),"On")==0)
ProblemType = LL2P_Stab_On;
*elseif(strcmp(GenData(LL2P_stabilize),"Off")==0)
ProblemType = LL2P_Stab_Off;
*endif
*#################  Iterative solver initial guess  ###########################
*if(strcmp(GenData(Initial_guess),"Read_file")==0)
ProblemType = Read_Guess_On;
*elseif(strcmp(GenData(Initial_guess),"Nil_vector")==0)
ProblemType = Read_Guess_Off;
*endif
*##################  Write solution vector on file  ###########################
*if(strcmp(GenData(Results_in_file),"No")==0)
ProblemType = Write_Solution_Off;
*elseif(strcmp(GenData(Results_in_file),"Every_step")==0)
ProblemType = Write_Solution_Steps;
*elseif(strcmp(GenData(Results_in_file),"Final_step")==0)
ProblemType = Write_Solution_Final;
*endif
*###################  Import Robin flux from files  ###########################
*if(strcmp(GenData(Import_Robin_flux),"On")==0)
ProblemType = Import_Robin_On;
*else
ProblemType = Import_Robin_Off;
*endif
*###################  Import J currents from files  ###########################
*if(strcmp(GenData(Import_J_currents),"On")==0)
ProblemType = Import_J_On;
*else
ProblemType = Import_J_Off;
*endif
*############  Import volumetric element matrices from files  #################
*if(strcmp(GenData(Import_ele_matrix),"On")==0)
ProblemType = Import_VEM_On;
*else
ProblemType = Import_VEM_Off;
*endif
*####################  Write results fields on files  #########################
*if(strcmp(GenData(_Export_fields),"On")==0)
ProblemType = Export_Fields_On;
*else
ProblemType = Export_Fields_Off;
*endif
*##########################  Output file format  ##############################
*if(strcmp(GenData(Output_format),"Ascii")==0)
ProblemType = Results_Format_Ascii;
*else
ProblemType = Results_Format_Binary;
*endif
*######################  Number of parallel procesess  ########################
ProblemType = *GenData(Processors)pr;
*############################  Problem frequency  #############################
*if(strcmp(GenData(Sweep_frequency_mode),"0")==0)
*if( (strcmp(GenData(Complex_frequency_mode),"1")==0) && (strcmp(GenData(Problem_mode),"Full_wave")==0) )
ProblemType = Load_Complex_Frequency;
*else
ProblemType = Load_Angular_Frequency;
*endif
*else
*set var dospi   (real) = 6.283185307179586476925286766559
*set var IniFreq (real) = GenData(Initial_frequency,real)
*set var EndFreq (real) = GenData(Final_frequency,real)
*set var StepFreq(real) = GenData(Step_frequency,real)
SweepFrequency = vTow(*operation(IniFreq*dospi),*operation(EndFreq*dospi),*operation(StepFreq*dospi));
*endif

// Check settings consistency
ProblemType = Check_Consistency;

// Read list of nodes
ProblemType = Read_Nodes_ID_Coord;

// Read material properties
ProblemType = Read_Material_Properties;

*if( (strcmp(GenData(Problem_mode),"Full_wave")==0) || (strcmp(GenData(Problem_mode),"Cold_plasma")==0) )
*if(strcmp(GenData(Problem_mode),"Cold_plasma")==0)
// Read plasma parameters
ProblemType = Load_Cold_Plasma_Parameters;

*endif
// Read plane waves parameters
ProblemType = Load_Plane_Waves_Parameters;

// Create high order nodes
ProblemType = Create_High_Order_Nodes;

// Make contact elements
*if(strcmp(GenData(Contact_detect),"Manual")==0)
ProblemType = Make_Contact;
*elseif(strcmp(GenData(Contact_detect),"Ignore")==0)
ProblemType = Ignr_Contact;
*endif

// Impose periodic boundary conditions 
ProblemType = Impose_PBC;

// Read PEC, PMC and TEC elements
ProblemType = Read_Dirichlet_Elements;

*if(strcmp(GenData(Problem_mode),"Cold_plasma")==0)
// Set E parallel to zero
ProblemType = Set_E_Parallel_ToZero;

*endif
// Read singular nodes list
ProblemType = Read_Nodes_Singular;

// Set un-gauged layers
ProblemType = Find_Sing_2L;
ProblemType = Find_Sing_3L;
ProblemType = Find_Sing_4L;
ProblemType = Find_Sing_5L;
ProblemType = Find_Sing_6L;

*if(strcmp(GenData(_Element_type),"LL2P_3sb")==0)
// Collapse bubbles on un-gauged layers
ProblemType = Collapse_Bubbles;
*if(strcmp(GenData(LL2P_keep_div),"On")==0)
ProblemType = Clear_Singularities;
*endif

*endif
*endif
// Read voltages on nodes 
ProblemType = Read_Nodes_Voltages;

// Initialize building 
ElementsGroup = Electromagnetic_Group;

*if(strcmp(GenData(Solving_mode),"Debug")==0)
// Print mesh
ProblemType = Print_High_Order_Mesh;

// Print "Debug" mode results
CalculateNodal(BOUNDARY_NORMALS);
PrintOnNodes(BOUNDARY_NORMALS);
CalculateNodal(CONTACT_NORMALS);
PrintOnNodes(CONTACT_NORMALS);
*if(strcmp(GenData(Problem_mode),"Cold_plasma")==0)
ProblemType = Show_All_Plasma_Parameters;
*endif

// End "Debug" mode
ProblemType = End_Debug_Mode;

*else
// Build linear system and solve
*if(strcmp(GenData(Sweep_frequency_mode),"1")==0)
ProblemType = Build;
*elseif(strcmp(GenData(Solving_mode),"Read")==0)
ProblemType = Read_Solve_Print_File;
*else
ProblemType = Build;
*endif
*endif