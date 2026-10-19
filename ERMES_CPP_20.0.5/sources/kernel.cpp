
#include "../includes/exception.h"
#include "../includes/kernel.h"
#include "../includes/process.h"
#include "../includes/gid_input.h"
#include "../includes/modeler.h"

#include "../ERMES/Utils/MatrixFileWritting.h"

namespace Kratos
{
    //**********************************************************************************************************************
    // - Node generation
    //**********************************************************************************************************************
    void Kernel::GenerateNode( IndexType NodeId, double X, double Y, double Z )
    {
        mpModeler->Generate_Node( NodeId, X, Y, Z );
    }
    
    //**********************************************************************************************************************
    //  - Creates linear solver object
    //**********************************************************************************************************************
    void Kernel::GenerateLinearSolver( const  String& Name         , 
                                       const  String& SolverType   , 
                                       int    MaxNumberOfIterations, 
                                       int    StepIterations       , 
                                       double Tolerance            )
    {
        if( mpModeler->Is_ColdPlasmaMode() )
        {
            mpModeler->Generate_LinearSolver_ColdPlasma( Name, SolverType, MaxNumberOfIterations, StepIterations, Tolerance );
        }
        else                                  
        {   
            mpModeler->Generate_LinearSolver_FullWave( Name, SolverType, MaxNumberOfIterations, StepIterations, Tolerance );
        }
    }

    //**********************************************************************************************************************
    //  - Set time on current time step
    //**********************************************************************************************************************    
    void Kernel::CreateTimeStep()
    {
        mpModel->SetCurrentTime( mTime );
    }
    
    //**********************************************************************************************************************
    //  - Creates solution step
    //**********************************************************************************************************************
    void Kernel::CreateSolutionStep()
    {
        mpModel->CreateWithCurrentSolutionStepData();
    }
    
    //**********************************************************************************************************************
    //  - Initialization and set of variables for sweep frequency mode
    //**********************************************************************************************************************
    void Kernel::SweepFrequency( double IniFreq, double EndFreq, double StepFreq )
    {
        if( mpModeler->Is_ElectrostaticMode() )
        {
            mFrequencySweep = false;

            mpModeler->Set_Frequency_Mode( mFrequencySweep );
        }
        else
        {
            mFrequencySweep = true;

            mpModeler->Set_Frequency_Mode( mFrequencySweep );
            
	        mIniFreq  = IniFreq ;
	        mEndFreq  = EndFreq ;
	        mStepFreq = StepFreq;

            mpModeler->Set_Frequency( mIniFreq );
        }
    }

    //**********************************************************************************************************************
    //  - Initialization of all boolean variables
    //**********************************************************************************************************************
    void Kernel::Ini_Boolean_Variables()
    {
        mBuild = false;   

        mFSing2L = false; 
        mFSing3L = false; 
        mFSing4L = false;
	    mFSing5L = false;
	    mFSing6L = false;

		mSetEdgeDOFs = false;
	    mPrintHOMesh = false;

	    mCreateHONodes = false; 
        mCreateContact = false;

        mPlasmaInterplte = false;
        mCollapseBubbles = false;

        mCalculate_E = false;
        mCalculate_H = false;
        mCalculate_B = false;
	    mCalculate_J = false;
	    mCalculate_JouleH = false;

        mExport_J_Electrostatic = false;

        mCalculate_E_Electrostatic = false;
        mCalculate_J_Electrostatic = false;
        mCalculate_JouleH_Electrostatic = false;
    }

    //**********************************************************************************************************************
    // - Clears all ResultsOnGPs vectors used in frequency sweep mode
    //**********************************************************************************************************************
    void Kernel::Clear_Results_OnGPs()
    {
        if( mREAL_E_OnGP.size() > 0 ) Clear_Results_OnGPs( mREAL_E_OnGP );
        if( mIMAG_E_OnGP.size() > 0 ) Clear_Results_OnGPs( mIMAG_E_OnGP );
                        
        if( mREAL_H_OnGP.size() > 0 ) Clear_Results_OnGPs( mREAL_H_OnGP );
        if( mIMAG_H_OnGP.size() > 0 ) Clear_Results_OnGPs( mIMAG_H_OnGP );
                              
        if( mREAL_J_OnGP.size() > 0 ) Clear_Results_OnGPs( mREAL_J_OnGP );
        if( mIMAG_J_OnGP.size() > 0 ) Clear_Results_OnGPs( mIMAG_J_OnGP );
    }

    //**********************************************************************************************************************
    // - Clears ResultsOnGPs vector
    //**********************************************************************************************************************
    void Kernel::Clear_Results_OnGPs( ResultsOnGPsType& ResultsOnGPs )
    {
        ResultsOnGPs.clear(); 

        Vector< Vector< Vector<double> > >().swap( ResultsOnGPs );
    }

    //**********************************************************************************************************************
    // - Assigns the corresponding ResultsOnGPsType Vector to a variable with Key == VarKey
    //**********************************************************************************************************************
    void Kernel::Assign_ResultsOnGPsVector( ResultsOnGPsType& ResultsOnGPsVector, int VarKey )
    {
        if     ( VarKey == E.getKey() ) ResultsOnGPsVector = mE_OnGP; 
        else if( VarKey == H.getKey() ) ResultsOnGPsVector = mH_OnGP; 
        else if( VarKey == J.getKey() ) ResultsOnGPsVector = mJ_OnGP; 
        else if( VarKey == B.getKey() ) ResultsOnGPsVector = mB_OnGP; 
        
        else if( VarKey == REAL_E.getKey() ) ResultsOnGPsVector = mREAL_E_OnGP; 
        else if( VarKey == IMAG_E.getKey() ) ResultsOnGPsVector = mIMAG_E_OnGP; 
        else if( VarKey ==  MOD_E.getKey() ) ResultsOnGPsVector =  mMOD_E_OnGP; 

        else if( VarKey == REAL_H.getKey() ) ResultsOnGPsVector = mREAL_H_OnGP; 
        else if( VarKey == IMAG_H.getKey() ) ResultsOnGPsVector = mIMAG_H_OnGP; 
        else if( VarKey ==  MOD_H.getKey() ) ResultsOnGPsVector =  mMOD_H_OnGP; 

        else if( VarKey == REAL_B.getKey() ) ResultsOnGPsVector = mREAL_B_OnGP; 
        else if( VarKey == IMAG_B.getKey() ) ResultsOnGPsVector = mIMAG_B_OnGP; 
        else if( VarKey ==  MOD_B.getKey() ) ResultsOnGPsVector =  mMOD_B_OnGP; 

        else if( VarKey == REAL_J.getKey() ) ResultsOnGPsVector = mREAL_J_OnGP; 
        else if( VarKey == IMAG_J.getKey() ) ResultsOnGPsVector = mIMAG_J_OnGP; 
        else if( VarKey ==  MOD_J.getKey() ) ResultsOnGPsVector =  mMOD_J_OnGP; 

        else if( VarKey == JOULE_HEATING  .getKey() ) ResultsOnGPsVector = mJOULE_HEATING_OnGP; 
        else if( VarKey == LORENTZ_FORCE  .getKey() ) ResultsOnGPsVector = mLORENTZ_FORCE_OnGP; 
        else if( VarKey == POYNTING_VECTOR.getKey() ) ResultsOnGPsVector = mPOYNTING_VECTOR_OnGP; 

        else if( VarKey == ELECTROSTATIC_ELECTRIC_FIELD .getKey() ) ResultsOnGPsVector = mElectrostatic_E_OnGP;     
        else if( VarKey == ELECTROSTATIC_CURRENT_DENSITY.getKey() ) ResultsOnGPsVector = mElectrostatic_J_OnGP;     
        else if( VarKey == ELECTROSTATIC_JOULE_HEATING  .getKey() ) ResultsOnGPsVector = mElectrostatic_JouleH_OnGP; 
    }

    //**********************************************************************************************************************
    // - Problem setting 
    //**********************************************************************************************************************
    void Kernel::SetProblemType( const String& ProblemType )
    {
	    Ini_Boolean_Variables();
        
	    if     ( ProblemType == "Build"          ) Build_Linear_System();
        else if( ProblemType == "End_Debug_Mode" ) End_Debuger_Mode   ();

        else if( ProblemType == "Read_Nodes_ID_Coord" ) Read_Nodes_ID_Coord_File();
        else if( ProblemType == "Read_Nodes_Singular" ) Read_Nodes_Singular_File();
        else if( ProblemType == "Read_Nodes_Voltages" ) Read_Nodes_Voltages_File();

        else if( ProblemType == "Read_Dirichlet_Elements" ) Read_Dirichlet_Elements_File ();
        else if( ProblemType == "Read_Projector_Elements" ) Read_Projector_Elements_File ();
        else if( ProblemType == "Set_E_Parallel_ToZero"   ) mpModeler->SetAll_Epar_ToZero();

        else if( ProblemType == "Print_High_Order_Mesh"   ) Print_High_Order_Mesh  ();
        else if( ProblemType == "Create_High_Order_Nodes" ) Create_High_Order_Nodes();

        else if( ProblemType == "Read_Solve_Print_File"  ) Read_Solve_Problem_File();
        else if( ProblemType == "Export_Static_Currents" ) Export_Static_Currents ();

        else if( ProblemType == "Impose_PBC"   ) Read_PBC_Elements_File();
        else if( ProblemType == "Make_Contact" ) Read_Contact_Elements_File( mCreateContact );
        else if( ProblemType == "Ignr_Contact" ) return;

        else if( ProblemType == "Read_Material_Properties"    ) Read_Material_Properties_File   ();
	    else if( ProblemType == "Load_Cold_Plasma_Parameters" ) Read_Cold_Plasma_Parameters_File();
        else if( ProblemType == "Load_Plane_Waves_Parameters" ) Read_Plane_Waves_Parameters_File();
        
        else if( ProblemType == "Load_Angular_Frequency" ) Read_Angular_Frequency_File();
        else if( ProblemType == "Load_Complex_Frequency" ) Read_Complex_Frequency_File();

        else if( ProblemType == "Collapse_Bubbles" ) Read_Volume_Elements_File( mCollapseBubbles );  

        else if( ProblemType == "Find_Sing_2L" ) Read_Volume_Elements_File( mFSing2L ); 
        else if( ProblemType == "Find_Sing_3L" ) Read_Volume_Elements_File( mFSing3L ); 
	    else if( ProblemType == "Find_Sing_4L" ) Read_Volume_Elements_File( mFSing4L ); 
	    else if( ProblemType == "Find_Sing_5L" ) Read_Volume_Elements_File( mFSing5L ); 
	    else if( ProblemType == "Find_Sing_6L" ) Read_Volume_Elements_File( mFSing6L );  

        else if( ProblemType == "Calculate_E" ) Calculate_E();
        else if( ProblemType == "Calculate_H" ) Calculate_H();
        else if( ProblemType == "Calculate_B" ) Calculate_B();
        else if( ProblemType == "Calculate_J" ) Calculate_J();

        else if( ProblemType == "Calculate_Static_E" ) Calculate_Static_E();
        
        else if( ProblemType == "Show_Joule_Heating"   ) Show_Joule_Heating  ();
        else if( ProblemType == "Show_Poynting_Vector" ) Show_Poynting_Vector();
        else if( ProblemType == "Show_Lorentz_Force"   ) Show_Lorentz_Force  ();

        else if( ProblemType == "Show_Static_E" ) Show_Static_E();
        else if( ProblemType == "Show_Static_J" ) Show_Static_J();

        else if( ProblemType == "Show_Static_Joule_Heating" ) Show_Static_Joule_Heating();

        else if( ProblemType == "Show_B_External"            ) Show_B_External           ();
	    else if( ProblemType == "Show_Electron_Density"      ) Show_Electron_Density     ();
	    else if( ProblemType == "Show_Permittivity_Tensor"   ) Show_Permittivity_Tensor  ();
	    else if( ProblemType == "Show_All_Plasma_Parameters" ) Show_All_Plasma_Parameters();

        else if( ProblemType == "Show_E_Parallel"       ) Show_E_Parallel       ();
        else if( ProblemType == "Show_Eparallel_t"      ) Show_E_Parallel_t     ();
        else if( ProblemType == "Show_E_Perpendicular"  ) Show_E_Perpendicular  ();
        else if( ProblemType == "Show_Eperpendicular_t" ) Show_E_Perpendicular_t();
 
	    else mpModeler->Set_Problem_Type( ProblemType );
    }

    //*********************************************************************************************
    // -  Calculates REAL_E, IMAG_E, MOD_E
    //*********************************************************************************************
    void Kernel::Calculate_E()
    {
        mpModeler->Ini_E_Derivation(); 

        if( ( mpModeler->Is_PotentialsOn    () == true  ) ||  
            ( mpModeler->Is_LL2PSmoothingOn () == true  ) || 
            ( mpModeler->Is_EdgeElementsOn  () == true  ) || 
            ( mpModeler->Is_ResultsOnNodesOn() == false )  ) 
        {
            Read_Volume_Elements_File( mCalculate_E );
 
            mpModeler->End_E_Derivation();
        }
    }

    //*********************************************************************************************
    // -  Calculates REAL_H, IMAG_H, MOD_H
    //*********************************************************************************************
    void Kernel::Calculate_H()
    {
        mpModeler->Ini_H_Derivation();

        Read_Volume_Elements_File( mCalculate_H );

        mpModeler->End_H_Derivation();
    }

    //*********************************************************************************************
    // -  Calculates REAL_B, IMAG_B, MOD_B
    //*********************************************************************************************
    void Kernel::Calculate_B()
    {
        mpModeler->Ini_B_Derivation();

        Read_Volume_Elements_File( mCalculate_B );

        mpModeler->End_B_Derivation();
    }

    //*********************************************************************************************
    // -  Calculates REAL_J, IMAG_J, MOD_J
    //*********************************************************************************************
    void Kernel::Calculate_J()
    {
        mpModeler->Ini_J_Calculation();

        Read_Source_Elements_File( mCalculate_J );  

        mpModeler->Calculate_Imported_Currents(); 

        Read_Volume_Elements_File( mCalculate_J );

        mpModeler->End_J_Calculation();
    }

    //*********************************************************************************************
    // -  Calculates the static electric field 
    //*********************************************************************************************
    void Kernel::Calculate_Static_E()
    {
        mpModeler->Ini_Electrostatic_E_Derivation();

        Read_Volume_Elements_File( mCalculate_E_Electrostatic );

        mpModeler->End_Electrostatic_E_Derivation();
    }

    //*********************************************************************************************
    // -  Prints in the results file the static electric field 
    //*********************************************************************************************
    void Kernel::Show_Static_E()
    {
        Print( ELECTROSTATIC_ELECTRIC_FIELD );
    }

    //*********************************************************************************************
    // -  Calculates and prints in the results file the static current density 
    //*********************************************************************************************
    void Kernel::Show_Static_J()
    {
        mpModeler->Ini_Electrostatic_J_Derivation();

        Read_Volume_Elements_File( mCalculate_J_Electrostatic );

        mpModeler->End_Electrostatic_J_Derivation();

        Print( ELECTROSTATIC_CURRENT_DENSITY );               
    }

    //*********************************************************************************************
    // -  Calculates and prints in the results file the static Joule heating
    //*********************************************************************************************
    void Kernel::Show_Static_Joule_Heating()
    {
        mpModeler->Ini_Electrostatic_JouleH_Derivation();

        Read_Volume_Elements_File( mCalculate_JouleH_Electrostatic );

        mpModeler->End_Electrostatic_JouleH_Derivation();

        Print( ELECTROSTATIC_JOULE_HEATING );
    }

    //*********************************************************************************************
    // - Export static currents
    //*********************************************************************************************
    void Kernel::Export_Static_Currents()
    {
        mpModeler->Ini_Exporting_Electrostatic_Current(); 

        Read_Volume_Elements_File( mExport_J_Electrostatic );

        mpModeler->End_Exporting_Electrostatic_Current();    
    }

    //*********************************************************************************************
    // - Calculates and prints Joule heating
    //*********************************************************************************************
    void Kernel::Show_Joule_Heating()
    {
        mpModeler->Ini_Joule_Heating_Calculation();
        
        Read_Volume_Elements_File( mCalculate_JouleH );  

        mpModeler->End_Joule_Heating_Calculation();

        Print( JOULE_HEATING );
    }

    //*********************************************************************************************
    // - Calculates and prints Poynting vector
    //*********************************************************************************************
    void Kernel::Show_Poynting_Vector()
    {
        Calculate( POYNTING_VECTOR );

        Print( POYNTING_VECTOR );   
    }

    //*********************************************************************************************
    // - Calculates and prints Lorentz vector
    //*********************************************************************************************
    void Kernel::Show_Lorentz_Force()
    {
        Calculate( LORENTZ_FORCE );

        Print( LORENTZ_FORCE );       
    }

    //*********************************************************************************************
    // - Prints electron density in plasma mode
    //*********************************************************************************************
    void Kernel::Show_Electron_Density()
    {
        mpModeler->Extract_ColdPlasma_ElectronDensity();

        PrintOnNodes( ELECTRON_DENSITY );
    }

    //*********************************************************************************************
    // - Prints B external in plasma mode
    //*********************************************************************************************
    void Kernel::Show_B_External()
    {
        mpModeler->Extract_ColdPlasma_Bexternal();

        PrintOnNodes( B_EXT );    
    }

    //*********************************************************************************************
    // - Prints permittivity tensor in plasma mode
    //*********************************************************************************************
    void Kernel::Show_Permittivity_Tensor()
    {
        mpModeler->Extract_ColdPlasma_PermittivityTensor();

        PrintOnNodes( PLASMA_PERMITTIVITY_S ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_D ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_P ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_R ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_L );  
    }

    //*********************************************************************************************
    // - Prints plasma parameters in plasma mode
    //*********************************************************************************************
    void Kernel::Show_All_Plasma_Parameters()
    {
        mpModeler->Extract_ColdPlasma_AllParameters();

        PrintOnNodes( ELECTRON_DENSITY      );
        PrintOnNodes( B_EXT                 );
        PrintOnNodes( PLASMA_PERMITTIVITY_S ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_D ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_P ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_R ); 
        PrintOnNodes( PLASMA_PERMITTIVITY_L ); 
    }

    //*********************************************************************************************
    // - Prints E parallel in plasma mode
    //*********************************************************************************************
    void Kernel::Show_E_Parallel()
    {
        mpModeler->Extract_ColdPlasma_Eparallel();

        PrintOnNodes( MOD_E_PARALLEL );   
    }

    //*********************************************************************************************
    // - Prints E parallel in time domain in plasma mode
    //*********************************************************************************************
    void Kernel::Show_E_Parallel_t()
    {
        CalculateNodal( E_PARALLEL_T );

        PrintOnNodes( E_PARALLEL_T );       
    }

    //*********************************************************************************************
    // - Prints E perpendicular in plasma mode
    //*********************************************************************************************
    void Kernel::Show_E_Perpendicular()
    {
        mpModeler->Extract_ColdPlasma_Eperpendicular();

        PrintOnNodes( MOD_E_PERPENDICULAR );       
    }

    //*********************************************************************************************
    // - Prints E perpendicular in time domain in plasma mode
    //********************************************************************************************* 
    void Kernel::Show_E_Perpendicular_t()
    {
        CalculateNodal( E_PERPENDICULAR_T );

        PrintOnNodes( E_PERPENDICULAR_T );       
    }

    //*********************************************************************************************
    // - Builds linear system
    //*********************************************************************************************
    void Kernel::Build_Linear_System()
    {
        std::cout << "Building linear system..." << std::endl;
   
        time_t bStart; 

        time( &bStart ); 

        if( mIniFreq == 0.0 ) 
        {
            mIniFreq += mStepFreq;
        }
        
        if( mFrequencySweep ) 
        {
            mpModeler->Set_Frequency( mIniFreq );
        }

        Read_Volume_Elements_File ( mBuild );
        Read_Source_Elements_File ( mBuild );
        Read_RobinBC_Elements_File( mBuild );

        mpModeler->Build_Imported_VolElmts();
        mpModeler->Build_Imported_RobinBCs();
        mpModeler->Build_Imported_Currents();

        time_t bEnd; 

        time( &bEnd );
        
        std::cout << "Building finished in " << difftime( bEnd, bStart )
                  << " seconds."             << std::endl << std::endl;
        
        Read_Solve_Problem_File();

        if( mFrequencySweep )
        {
            Clear_Results_OnGPs();

            mpModeler->Finish_First_Step();
            
            for( double freq = mIniFreq+mStepFreq; freq <= mEndFreq; freq += mStepFreq )
            {
                if( freq == 0.0 ) 
                {
                    continue;
                }
                
                std::cout << "Building linear system..." << std::endl;
              
                time( &bStart );
                      
                mpModeler->Set_Frequency( freq );
        
                Read_Volume_Elements_File ( mBuild );
                Read_Source_Elements_File ( mBuild );
                Read_RobinBC_Elements_File( mBuild );
        
                mpModeler->Build_Imported_VolElmts();
                mpModeler->Build_Imported_RobinBCs();
                mpModeler->Build_Imported_Currents();
        
        	    time( &bEnd );
                
        	    std::cout << "Building finished in " << difftime(bEnd, bStart)
        		          << " seconds."             << std::endl << std::endl;
                    
                Read_Solve_Problem_File();
                
                Clear_Results_OnGPs();

                mpModeler->Finish_Step();
            }
        }
    }

    //*********************************************************************************************
    // - Prints high order mesh
    //*********************************************************************************************
    void Kernel::Print_High_Order_Mesh()
    {
        mpDefaultOutput->SetFileFormat( mpModeler->Is_OutputFileFormat_ASCII() );

        if( mpModeler->Is_ElectrostaticMode() )  
        {
            mpModel->SetNumPrintableNodes( mpModel->GetNodesArraySize() );
        }

        mpDefaultOutput->PrintMeshNodes( mpModeler->Get_Element_Order() );
        
        mElementIdCounter = 0;

        Read_Volume_Elements_File( mPrintHOMesh );
        
        mpDefaultOutput->CloseHOMeshFile();

        mpDefaultOutput->PrintGaussPointsHeader();
    }

    //*********************************************************************************************
    // - Reads nodes IDs and coordinates list file
    //*********************************************************************************************
    void Kernel::Read_Nodes_ID_Coord_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-1.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads voltages on nodes list file
    //*********************************************************************************************
    void Kernel::Read_Nodes_Voltages_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-2.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads singular nodes list file
    //*********************************************************************************************
    void Kernel::Read_Nodes_Singular_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-3.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads material properties file
    //*********************************************************************************************
    void Kernel::Read_Material_Properties_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-11.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads Dirichlet elements file ( PEC, PMC and TEC boundary conditions )
    //*********************************************************************************************
    void Kernel::Read_Dirichlet_Elements_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-6.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads projection elements file ( waveguide ports and field integrals )
    //*********************************************************************************************
    void Kernel::Read_Projector_Elements_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-10.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads volume elements file 
    //*********************************************************************************************
    void Kernel::Read_Volume_Elements_File( bool& Active_Magnitude )
    {
        Active_Magnitude = true;
        
        String file_name = mBaseFileName; 

        file_name << "-4.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );

        Active_Magnitude = false;
    }

    //*********************************************************************************************
    // - Reads source elements file 
    //*********************************************************************************************
    void Kernel::Read_Source_Elements_File( bool& Active_Magnitude )
    {
        Active_Magnitude = true;
        
        String file_name = mBaseFileName; 

        file_name << "-5.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );

        Active_Magnitude = false;
    }

    //*********************************************************************************************
    // - Reads Robin boundary conditions elements file 
    //*********************************************************************************************
    void Kernel::Read_RobinBC_Elements_File( bool& Active_Magnitude )
    {
        Active_Magnitude = true;
        
        String file_name = mBaseFileName; 

        file_name << "-7.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );

        Active_Magnitude = false;
    }

    //*********************************************************************************************
    // - Reads Periodic Boundary Condition (PBC) elements file
    //*********************************************************************************************
    void Kernel::Read_PBC_Elements_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-8.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads contact elements file
    //*********************************************************************************************
    void Kernel::Read_Contact_Elements_File( bool& Active_Magnitude )
    {
        Active_Magnitude = true;
        
        String file_name = mBaseFileName; 

        file_name << "-9.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );

        Active_Magnitude = false;
    }

    //*********************************************************************************************
    // - Reads cold plasma parameters file
    //*********************************************************************************************
    void Kernel::Read_Cold_Plasma_Parameters_File()
    {
        String file_name = mBaseFileName;

        file_name << "-12.dat";
        
        mpModeler->LoadPlasmaParameters( file_name );
    }

    //*********************************************************************************************
    // - Reads plane waves parameters file
    //*********************************************************************************************
    void Kernel::Read_Plane_Waves_Parameters_File()
    {
        String file_name = mBaseFileName;

        file_name << "-13.dat";
        
        mpModeler->LoadPWavesParameters( file_name );
    }

    //*********************************************************************************************
    // - Reads angular frequency file
    //*********************************************************************************************
    void Kernel::Read_Angular_Frequency_File()
    {
        mFrequencySweep = false;

        mpModeler->Set_Frequency_Mode( mFrequencySweep );
        
        String file_name = mBaseFileName;

        file_name << "-16.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Reads complex frequency file
    //*********************************************************************************************
    void Kernel::Read_Complex_Frequency_File()
    {
        mFrequencySweep = false;
        
        String file_name = mBaseFileName;

        file_name << "-17.dat";
        
        mpModeler->LoadComplexFrequency( file_name );
    }

    //*********************************************************************************************
    // - Reads solve problem file
    //*********************************************************************************************
    void Kernel::Read_Solve_Problem_File()
    {
        String file_name = mBaseFileName; 

        file_name << "-18.dat";

        GidInput inputDataFile( file_name ); 

        inputDataFile.Read( this );
    }

    //*********************************************************************************************
    // - Ends debuger mode
    //*********************************************************************************************
    void Kernel::End_Debuger_Mode()
    {
        std::cout << std::endl;
        std::cout << "ERMES (debug mode) finished." << std::endl;
        std::cout << std::endl;
        
        mpDefaultOutput->~Output();

        exit( EXIT_SUCCESS );
    }

    //*********************************************************************************************
    // - Set problem frequency value
    //*********************************************************************************************
    void Kernel::SetFrequency( double ProblemFrequency )
    {
        mpModeler->Set_Frequency( ProblemFrequency );
    }
    
    //*********************************************************************************************
    // - Solve problem
    //*********************************************************************************************
    void Kernel::Solve( const String& SolvingStrategyName, const String& ElementsGroup )
    { 
        if( mpModeler->Is_ElectrostaticMode() ) 
        {
            mpModeler->Solve_Problem_Electrostatic(); 
        }
        else    
        {
            mpModeler->Solve_Problem_FullWave();
        }
    }
    
    //*********************************************************************************************
    // - Set DOFs 
    //*********************************************************************************************
    void Kernel::Ini_Building()
    {
        std::cout << "Done." << std::endl << std::endl;

		std::cout << "Setting DOFs and sparse matrix..." << std::endl;
		
		mpModeler->Set_Node_Elements_DOFs();

		if( mpModeler->Is_EdgeElementsOn() ) 
		{
			Read_Volume_Elements_File( mSetEdgeDOFs );
		}

		if( mpModeler->Is_Plasma_Profile_Resized() ) 
		{
			Read_Volume_Elements_File( mPlasmaInterplte ); 
		}

		mpModeler->Set_System_Matrix_Size();

		std::cout << "Done." << std::endl << std::endl;
    }

    //*********************************************************************************************
    // -  Set values
    //*********************************************************************************************
    void Kernel::SetValue( int NodeId, const Variable<double>& rVariable, const double& rValue )
    {
        if( rVariable.Name() == "T" ) mTime = rValue;

        mpModel->Value( rVariable, *(mpModel->GetNode(NodeId)) ) = rValue;
    }
    
    //*********************************************************************************************
    // -  Set values
    //*********************************************************************************************
    void Kernel::SetValue( const Variable<double>& rVariable, const double& rValue )
    {
        if( rVariable.Name() == "T" ) mTime = rValue;
        
        mpModel->Value( rVariable ) = rValue;
    }

    //*********************************************************************************************
    // -  Fix real DOFs
    //*********************************************************************************************
    void Kernel::FixDof( IndexType NodeId, const Variable<double>& rVariable, const double& Value )
    {
        mpModeler->FixDof( NodeId, rVariable, Value );
    }

    //*********************************************************************************************
    // -  Fix complex DOFs
    //*********************************************************************************************
    void Kernel::FixCDof( IndexType NodeId, const Variable<double>& rVariable, const Vector<double>& vValue )
    {
        mpModeler->FixCDof( NodeId, rVariable, vValue );
    }

    //*********************************************************************************************
    // -  Prints 1st order element
    //*********************************************************************************************
    void Kernel::PrintElement_1stOrder( int* NodesId, unsigned int PropertiesId )
    {
        int NodesIdMat[ 5 ];
    
        NodesIdMat[ 0 ] = NodesId[ 0 ];
        NodesIdMat[ 1 ] = NodesId[ 1 ];
        NodesIdMat[ 2 ] = NodesId[ 2 ];
        NodesIdMat[ 3 ] = NodesId[ 3 ];

        NodesIdMat[ 4 ] = PropertiesId;
    
        mElementIdCounter++;

        mpDefaultOutput->PrintHOElement( mElementIdCounter, NodesIdMat );
    }
        
    //*********************************************************************************************
    // -  Prints 2nd order element
    //*********************************************************************************************
    void Kernel::PrintElement_2ndOrder( int* NodesId, unsigned int PropertiesId )
    {
        Vector<int> HONodesId;
    
        mpModeler->Push_HONodes_OnVolume( NodesId, HONodesId );
    
        int NodesIdMat[ 11 ];
    
        NodesIdMat[ 0 ] = HONodesId[ 0 ];
        NodesIdMat[ 1 ] = HONodesId[ 1 ];
        NodesIdMat[ 2 ] = HONodesId[ 2 ];
        NodesIdMat[ 3 ] = HONodesId[ 3 ];
        NodesIdMat[ 4 ] = HONodesId[ 4 ];
        NodesIdMat[ 5 ] = HONodesId[ 5 ];
        NodesIdMat[ 6 ] = HONodesId[ 6 ];
        NodesIdMat[ 7 ] = HONodesId[ 7 ];
        NodesIdMat[ 8 ] = HONodesId[ 8 ];
        NodesIdMat[ 9 ] = HONodesId[ 9 ];

        NodesIdMat[ 10 ] = PropertiesId;
    
        mElementIdCounter++;

        mpDefaultOutput->PrintHOElement( mElementIdCounter, NodesIdMat );
    }

    //*********************************************************************************************
    // - Generates contact nodes pairs
    //*********************************************************************************************
    void Kernel::GenerateContactPairs( Vector<int>& CtCNodesId )
    {
        mpModeler->Generate_Contact_Pairs( CtCNodesId );
    }
    
    //*********************************************************************************************
    // - Perfect Electric Conductor (PEC)
    //*********************************************************************************************
    void Kernel::GeneratePECElement( Vector<int>& NodesId )
    {
        mpModeler->Generate_PEC_Element( NodesId );
    }
    
    //*********************************************************************************************
    // - Perfect Magnetic Conductor (PMC)
    //*********************************************************************************************
    void Kernel::GeneratePMCElement( Vector<int>& NodesId )
    {
        mpModeler->Generate_PMC_Element( NodesId );
    }
    
    //*********************************************************************************************
    // - Tranversal electric condition ( ortogonal to a PMC surface ) 
    //*********************************************************************************************
    void Kernel::GenerateTECElement( Vector<int>& NodesId )
    {
        mpModeler->Generate_TEC_Element( NodesId );
    }
   
    //*********************************************************************************************
    // - Compute Sij parameters for coaxial TEM mode
    //*********************************************************************************************
    void Kernel::GenerateProjectionCoaxialTEM( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Generate_Projection_CoaxPort_TEM( NodesId, PropertiesId );
    }
    
    //*********************************************************************************************
    // - Compute Sij parameters for rectangular waveguide TE10 mode
    //*********************************************************************************************
    void Kernel::GenerateProjectionRWPortTE10( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Generate_Projection_RectPort_TE10( NodesId, PropertiesId );
    }

    //*********************************************************************************************
    // - Generic Robin element
    //*********************************************************************************************
    void Kernel::GenerateGenericRobinElement( int* NodesId, unsigned int PropertiesId )
    {
        if( mpModeler->Is_ElectrostaticMode() ) 
        {
            mpModeler->Generate_GenericRobin_Element_Electrostatic( NodesId, PropertiesId );
        }
        else            
        {
            mpModeler->Generate_GenericRobin_Element_FullWave( NodesId, PropertiesId );
        }
    }
    
    //*********************************************************************************************
    // - Far field element
    //*********************************************************************************************
    void Kernel::GenerateFarFieldElement( int* NodesId )
    {
        mpModeler->Generate_FarField_Element_FullWave( NodesId );
    }
    
    //*********************************************************************************************
    // - Cold palsma far field element
    //*********************************************************************************************
    void Kernel::GeneratePlasmaFarFieldElement( int* NodesId, unsigned int PropertiesId )
    {  
        if( mpModeler->Is_ColdPlasmaMode() ) 
        {
            mpModeler->Generate_FarField_Element_ColdPlasma( NodesId, PropertiesId );
        }
        else
        {
            mpModeler->Generate_FarField_Element_FullWave( NodesId );
        }
    }
    
    //*********************************************************************************************
    // - Coaxial port element
    //*********************************************************************************************
    void Kernel::GenerateCoaxPortTEMElement( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Generate_CoaxPort_TEM_Element( NodesId, PropertiesId );
    }
    
    //*********************************************************************************************
    // - Rectangular waveguide element
    //*********************************************************************************************
    void Kernel::GenerateRWPortTE10Element( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Generate_RectPort_TE10_Element( NodesId, PropertiesId );
    }   

    //*********************************************************************************************
    // - J source elements
    //*********************************************************************************************
    void Kernel::GenerateSourceElement( int* NodesId, unsigned int PropertiesId )
    {
        if( mBuild ) 
        {
            mpModeler->Generate_Source_Element( NodesId, PropertiesId );
        }
        else if( mCalculate_J )   
        {
            mpModeler->Calculate_J_Imposed_OnNodes( NodesId, PropertiesId );
        }
    }

    //*********************************************************************************************
    // - Volumetric integral
    //*********************************************************************************************
    void Kernel::GenerateVolIntgElement( int* NodesId, unsigned int PropertiesId )
    {
        if( mpModeler->Is_ElectrostaticMode() ) 
        {
            mpModeler->Generate_Vol_Intgrl_Element_Electrostatic( NodesId, PropertiesId );
        }
        else  
        {
            mpModeler->Generate_Vol_Intgrl_Element_FullWave( NodesId, PropertiesId );
        }
    }
   
    //*********************************************************************************************
    // - Periodic boundary conditions
    //*********************************************************************************************
    void Kernel::GeneratePBCElement( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Generate_PBC_Element( NodesId, PropertiesId );
    }
   
    //*********************************************************************************************
    // - Surface integrals
    //*********************************************************************************************
    void Kernel::GenerateSrfIntgElement( int* NodesId, unsigned int PropertiesId )
    {
        if( mpModeler->Is_ElectrostaticMode() ) 
        {
            mpModeler->Generate_Srf_Intgrl_Element_Electrostatic( NodesId, PropertiesId );
        }
        else  
        {
            mpModeler->Generate_Srf_Intgrl_Element_FullWave( NodesId, PropertiesId );
        }
    }

    //**********************************************************************************************************************
    // - Volume elements
    //**********************************************************************************************************************
    void Kernel::GenerateVolumeElement( int* NodesId, unsigned int PropertiesId )
    { 
		if     ( mBuild         ) Build_Volume_Element    ( NodesId, PropertiesId );
        else if( mPrintHOMesh   ) Print_High_Order_Mesh   ( NodesId, PropertiesId );
        else if( mCreateHONodes ) Create_HONodes_InElement( NodesId, PropertiesId );

		else if( mFSing2L ) mpModeler->FindSing2L( NodesId );                                
        else if( mFSing3L ) mpModeler->FindSing3L( NodesId );                                
        else if( mFSing4L ) mpModeler->FindSing4L( NodesId );                               
        else if( mFSing5L ) mpModeler->FindSing5L( NodesId );                         
        else if( mFSing6L ) mpModeler->FindSing6L( NodesId );

		else if( mSetEdgeDOFs     ) mpModeler->Set_Edge_Elements_DOFs    ( NodesId );  
        else if( mCollapseBubbles ) mpModeler->Collapse_LL2P_Bubbles     ( NodesId );
        else if( mPlasmaInterplte ) mpModeler->Plasma_Profile_Interpolate( NodesId );

        else if( mExport_J_Electrostatic ) mpModeler->Export_Electrostatic_J_Element_GP( NodesId, PropertiesId );

        else if( mCalculate_E ) mpModeler->Calculate_E_Element( NodesId, PropertiesId, mREAL_E_OnGP, mIMAG_E_OnGP );
        else if( mCalculate_H ) mpModeler->Calculate_H_Element( NodesId, PropertiesId, mREAL_H_OnGP, mIMAG_H_OnGP );
        else if( mCalculate_B ) mpModeler->Calculate_B_Element( NodesId, PropertiesId, mREAL_B_OnGP, mIMAG_B_OnGP );
        else if( mCalculate_J ) mpModeler->Calculate_J_Element( NodesId, PropertiesId, mREAL_J_OnGP, mIMAG_J_OnGP );
        
        else if( mCalculate_E_Electrostatic      ) mpModeler->Calculate_E_Element_Electrostatic     ( NodesId, PropertiesId, mElectrostatic_E_OnGP      );
        else if( mCalculate_J_Electrostatic      ) mpModeler->Calculate_J_Element_Electrostatic     ( NodesId, PropertiesId, mElectrostatic_J_OnGP      );
        else if( mCalculate_JouleH               ) mpModeler->Calculate_JouleH_Element_FullWave     ( NodesId, PropertiesId, mJOULE_HEATING_OnGP        );
        else if( mCalculate_JouleH_Electrostatic ) mpModeler->Calculate_JouleH_Element_Electrostatic( NodesId, PropertiesId, mElectrostatic_JouleH_OnGP );
    }

    //**********************************************************************************************************************
    // - Prints high order mesh
    //**********************************************************************************************************************
    void Kernel::Print_High_Order_Mesh( int* NodesId, unsigned int PropertiesId )
    {
        PrintElement_1stOrder( NodesId, PropertiesId );
    }

    //**********************************************************************************************************************
    // - Creates high order nodes 
    //**********************************************************************************************************************
    void Kernel::Create_High_Order_Nodes()
    {
        mpModel->SetNumPrintableNodes( mpModel->GetNodesArraySize() );

        Read_Volume_Elements_File( mCreateHONodes );
    }

    //**********************************************************************************************************************
    // - Creates high order nodes in tetrahedral element
    //**********************************************************************************************************************
    void Kernel::Create_HONodes_InElement( int* NodesId, unsigned int PropertiesId )
    {
        mpModeler->Create_HONodes( NodesId );

        mpModeler->Contact_Properties( NodesId, PropertiesId );     
    }

    //**********************************************************************************************************************
    // - Builds linear system with volume elements
    //**********************************************************************************************************************
    void Kernel::Build_Volume_Element( int* NodesId, unsigned int PropertiesId )
    {
        if( mpModeler->Is_ElectrostaticMode() ) 
        {
            mpModeler->Generate_Volume_Element_Electrostatic( NodesId, PropertiesId );
        }
        else if( mpModeler->Is_ColdPlasmaMode() ) 
        {
            mpModeler->Generate_Volume_Element_ColdPlasma( NodesId, PropertiesId );
        }
        else 
        {
            mpModeler->Generate_Volume_Element_FullWave( NodesId, PropertiesId ); 
        }
    }

    //**********************************************************************************************************************
    // - Calculates variable on Gauss points
    //**********************************************************************************************************************
    void Kernel::CalculateOnGaussPoints( int VarKey )
    {
        // E_(t)
        if( VarKey == E.getKey() )
        {
            int NumElements    = mREAL_E_OnGP   .size();
            int NumGaussPoints = mREAL_E_OnGP[0].size();

            double time = mpModel->CurrentTime();
            double freq = mpModeler->Get_Problem_Frequency();
                
            VectorType rE( 3 );
            VectorType iE( 3 );
            VectorType Et( 3 );

            Vector2Type E_OnElement( NumGaussPoints );

            mE_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rE = mREAL_E_OnGP[ei][gp]; 
                    iE = mIMAG_E_OnGP[ei][gp];

                    Et[0] = rE[0] * cos( freq * time ) + iE[0] * sin( freq * time );                                                  
                    Et[1] = rE[1] * cos( freq * time ) + iE[1] * sin( freq * time );                                                   
                    Et[2] = rE[2] * cos( freq * time ) + iE[2] * sin( freq * time );
                
                    E_OnElement[gp] = Et;
                }

                mE_OnGP[ei] = E_OnElement ;
            }
        }
        // H_(t)
        else if( VarKey == H.getKey() )
        {
            int NumElements    = mREAL_H_OnGP   .size();
            int NumGaussPoints = mREAL_H_OnGP[0].size();
            
            double time = mpModel->CurrentTime();
            double freq = mpModeler->Get_Problem_Frequency();
            
            VectorType rH( 3 );
            VectorType iH( 3 );
            VectorType Ht( 3 );

            Vector2Type H_OnElement( NumGaussPoints );

            mH_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rH = mREAL_H_OnGP[ei][gp]; 
                    iH = mIMAG_H_OnGP[ei][gp];

                    Ht[0] = rH[0] * cos( freq * time ) + iH[0] * sin( freq * time );                                                  
                    Ht[1] = rH[1] * cos( freq * time ) + iH[1] * sin( freq * time );                                                   
                    Ht[2] = rH[2] * cos( freq * time ) + iH[2] * sin( freq * time );
                
                    H_OnElement[gp] = Ht;
                }

                mH_OnGP[ei] = H_OnElement ;
            }
        }
        // J_(t)
        else if( VarKey == J.getKey() )
        {
            int NumElements    = mREAL_J_OnGP   .size();
            int NumGaussPoints = mREAL_J_OnGP[0].size();

            double time = mpModel->CurrentTime();
            double freq = mpModeler->Get_Problem_Frequency();
            
            VectorType rJ( 3 );
            VectorType iJ( 3 );
            VectorType Jt( 3 );

            Vector2Type J_OnElement( NumGaussPoints );

            mJ_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rJ = mREAL_J_OnGP[ei][gp]; 
                    iJ = mIMAG_J_OnGP[ei][gp];

                    Jt[0] = rJ[0] * cos( freq * time ) + iJ[0] * sin( freq * time );                                                  
                    Jt[1] = rJ[1] * cos( freq * time ) + iJ[1] * sin( freq * time );                                                   
                    Jt[2] = rJ[2] * cos( freq * time ) + iJ[2] * sin( freq * time );
                
                    J_OnElement[gp] = Jt;
                }

                mJ_OnGP[ei] = J_OnElement ;
            }
        }
        // B_(t)
        else if( VarKey == B.getKey() )
        {
            int NumElements    = mREAL_B_OnGP   .size();
            int NumGaussPoints = mREAL_B_OnGP[0].size();

            double time = mpModel->CurrentTime();
            double freq = mpModeler->Get_Problem_Frequency();
            
            VectorType rB( 3 );
            VectorType iB( 3 );
            VectorType Bt( 3 );

            Vector2Type B_OnElement( NumGaussPoints );

            mB_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rB = mREAL_B_OnGP[ei][gp]; 
                    iB = mIMAG_B_OnGP[ei][gp];

                    Bt[0] = rB[0] * cos( freq * time ) + iB[0] * sin( freq * time );                                                  
                    Bt[1] = rB[1] * cos( freq * time ) + iB[1] * sin( freq * time );                                                   
                    Bt[2] = rB[2] * cos( freq * time ) + iB[2] * sin( freq * time );
                
                    B_OnElement[gp] = Bt;
                }

                mB_OnGP[ei] = B_OnElement ;
            }
        }
        // MOD_E
        else if( VarKey == MOD_E.getKey() )
        {
            int NumElements    = mREAL_E_OnGP   .size();
            int NumGaussPoints = mREAL_E_OnGP[0].size();
        
            VectorType rE( 3 );
            VectorType iE( 3 );

            VectorType ModE_OnGP( 1 );

            Vector2Type ModE_OnElement( NumGaussPoints );

            mMOD_E_OnGP.resize( NumElements );
        
            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rE = mREAL_E_OnGP[ei][gp]; 
                    iE = mIMAG_E_OnGP[ei][gp];

                    ModE_OnGP[0] = sqrt( rE[0]*rE[0] + iE[0]*iE[0] +
                                         rE[1]*rE[1] + iE[1]*iE[1] +
                                         rE[2]*rE[2] + iE[2]*iE[2] );
                
                    ModE_OnElement[gp] = ModE_OnGP;
                }

                mMOD_E_OnGP[ei] = ModE_OnElement;
            }
        }
        // MOD_H
        else if( VarKey == MOD_H.getKey() )
        {
            int NumElements    = mREAL_H_OnGP   .size();
            int NumGaussPoints = mREAL_H_OnGP[0].size();
        
            VectorType rH( 3 );
            VectorType iH( 3 );

            VectorType ModH_OnGP( 1 );

            Vector2Type ModH_OnElement( NumGaussPoints );

            mMOD_H_OnGP.resize( NumElements );            
        
            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rH = mREAL_H_OnGP[ei][gp];
                    iH = mIMAG_H_OnGP[ei][gp];

                    ModH_OnGP[0] = sqrt( rH[0]*rH[0] + iH[0]*iH[0] +
                                         rH[1]*rH[1] + iH[1]*iH[1] +
                                         rH[2]*rH[2] + iH[2]*iH[2] );
                
                    ModH_OnElement[gp] = ModH_OnGP;
                }

                mMOD_H_OnGP[ei] = ModH_OnElement;
            }
        }
        // MOD_J
        else if( VarKey == MOD_J.getKey() )
        {
            int NumElements    = mREAL_J_OnGP   .size();
            int NumGaussPoints = mREAL_J_OnGP[0].size();
        
            VectorType rJ( 3 );
            VectorType iJ( 3 );

            VectorType ModJ_OnGP( 1 );

            Vector2Type ModJ_OnElement( NumGaussPoints );

            mMOD_J_OnGP.resize( NumElements );                        
        
            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rJ = mREAL_J_OnGP[ei][gp];
                    iJ = mIMAG_J_OnGP[ei][gp];

                    ModJ_OnGP[0] = sqrt( rJ[0]*rJ[0] + iJ[0]*iJ[0] +
                                         rJ[1]*rJ[1] + iJ[1]*iJ[1] +
                                         rJ[2]*rJ[2] + iJ[2]*iJ[2] );
                
                    ModJ_OnElement[gp] = ModJ_OnGP;
                }

                mMOD_J_OnGP[ei] = ModJ_OnElement;
            }
        }
        // MOD_B
        else if( VarKey == MOD_B.getKey() )
        {
            int NumElements    = mREAL_B_OnGP   .size();
            int NumGaussPoints = mREAL_B_OnGP[0].size();
        
            VectorType rB( 3 );
            VectorType iB( 3 );

            VectorType ModB_OnGP( 1 );

            Vector2Type ModB_OnElement( NumGaussPoints );

            mMOD_B_OnGP.resize( NumElements );                          
        
            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rB = mREAL_B_OnGP[ei][gp];
                    iB = mIMAG_B_OnGP[ei][gp];

                    ModB_OnGP[0] = sqrt( rB[0]*rB[0] + iB[0]*iB[0] +
                                         rB[1]*rB[1] + iB[1]*iB[1] +
                                         rB[2]*rB[2] + iB[2]*iB[2] );
                
                    ModB_OnElement[gp] = ModB_OnGP;
                }

                mMOD_B_OnGP[ei] = ModB_OnElement;
            }
        }
        // POYNTING_VECTOR
        else if( VarKey == POYNTING_VECTOR.getKey() )
        {
            int NumElements    = mREAL_E_OnGP   .size();
            int NumGaussPoints = mREAL_E_OnGP[0].size();

            VectorType rE( 3 ), iE( 3 );
            VectorType rH( 3 ), iH( 3 );

            VectorType avgS( 3 );
        
            cVectorType cE( 3 );
            cVectorType cH( 3 );

            cVectorType cE_X_conj_H( 3 );

            Vector2Type avgS_OnElement( NumGaussPoints );

            mPOYNTING_VECTOR_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rE = mREAL_E_OnGP[ei][gp]; 
                    iE = mIMAG_E_OnGP[ei][gp];

                    cE[0] = std::complex<double>( rE[0], iE[0] );
                    cE[1] = std::complex<double>( rE[1], iE[1] );
                    cE[2] = std::complex<double>( rE[2], iE[2] );

                    rH = mREAL_H_OnGP[ei][gp]; 
                    iH = mIMAG_H_OnGP[ei][gp];

                    cH[0] = std::complex<double>( rH[0], iH[0] );
                    cH[1] = std::complex<double>( rH[1], iH[1] );
                    cH[2] = std::complex<double>( rH[2], iH[2] );

                    // E X conj(H)
                    cE_X_conj_H[0] = cE[1] * std::conj( cH[2] ) - cE[2] * std::conj( cH[1] );
                    cE_X_conj_H[1] = cE[2] * std::conj( cH[0] ) - cE[0] * std::conj( cH[2] );
                    cE_X_conj_H[2] = cE[0] * std::conj( cH[1] ) - cE[1] * std::conj( cH[0] );
                    
                    // < S > = 0.5 * real( E X conj(H) )
                    avgS[0] = 0.5 * std::real( cE_X_conj_H[0] );
                    avgS[1] = 0.5 * std::real( cE_X_conj_H[1] );
                    avgS[2] = 0.5 * std::real( cE_X_conj_H[2] );
                    
                    avgS_OnElement[gp] = avgS;
                }

                mPOYNTING_VECTOR_OnGP[ei] = avgS_OnElement;
            }            
        }
        // LORENTZ_FORCE
        else if( VarKey == LORENTZ_FORCE.getKey() )
        {
            int NumElements    = mREAL_J_OnGP   .size();
            int NumGaussPoints = mREAL_J_OnGP[0].size();

            VectorType rJ( 3 ), iJ( 3 );
            VectorType rB( 3 ), iB( 3 );

            VectorType avgF( 3 );

            cVectorType cJ( 3 );
            cVectorType cB( 3 );

            cVectorType cJ_X_conj_B( 3 );

            Vector2Type avgF_OnElement( NumGaussPoints );

            mLORENTZ_FORCE_OnGP.resize( NumElements );

            // Loop over elements
            for( int ei=0; ei<NumElements; ei++ )
            {
                // Loop over GPs inside element
                for( int gp=0; gp<NumGaussPoints; gp++ )
                {
                    rJ = mREAL_J_OnGP[ei][gp]; 
                    iJ = mIMAG_J_OnGP[ei][gp];

                    cJ[0] = std::complex<double>( rJ[0], iJ[0] );
                    cJ[1] = std::complex<double>( rJ[1], iJ[1] );
                    cJ[2] = std::complex<double>( rJ[2], iJ[2] );

                    rB = mREAL_B_OnGP[ei][gp]; 
                    iB = mIMAG_B_OnGP[ei][gp];

                    cB[0] = std::complex<double>( rB[0], iB[0] );
                    cB[1] = std::complex<double>( rB[1], iB[1] );
                    cB[2] = std::complex<double>( rB[2], iB[2] );

                    // J X conj(B)
                    cJ_X_conj_B[0] = cJ[1] * std::conj( cB[2] ) - cJ[2] * std::conj( cB[1] );
                    cJ_X_conj_B[1] = cJ[2] * std::conj( cB[0] ) - cJ[0] * std::conj( cB[2] );
                    cJ_X_conj_B[2] = cJ[0] * std::conj( cB[1] ) - cJ[1] * std::conj( cB[0] );
                    
                    // < F > = 0.5 * real( J X conj(B) )
                    avgF[0] = 0.5 * std::real( cJ_X_conj_B[0] );
                    avgF[1] = 0.5 * std::real( cJ_X_conj_B[1] );
                    avgF[2] = 0.5 * std::real( cJ_X_conj_B[2] );
                    
                    avgF_OnElement[gp] = avgF;
                }

                mLORENTZ_FORCE_OnGP[ei] = avgF_OnElement;
            }            
        }
    }

    //**********************************************************************************************************************
    // - Main
    //**********************************************************************************************************************
    void Kernel::Execute( const String& BaseFileName )
    {
        mBaseFileName = BaseFileName;

	    mpModeler->SetBaseFileName( mBaseFileName );
    
	    String results_file_name = BaseFileName;

	    mpDefaultOutput = Output::Pointer( new GidOutput( mpModel, results_file_name ) );
    
        String nodes_file_name = BaseFileName;

        nodes_file_name << ".dat";

        GidInput input( nodes_file_name );
    
        try
	    {
	        input.Read( this );
	    }
	    catch (...)
	    {
	        std::cout << std::endl;
	        std::cout << "!!!!!!!!!!!!!! External Error !!!!!!!!!!!!!!!" << std::endl;
	        std::cout << std::endl;
	        std::cout << "Please, check memory quotas or input files..." << std::endl;
	        std::cout << std::endl;
	        std::cout << "ERMES analysis finished.                     " << std::endl;
	        std::cerr << "ERMES error: exception while reading input files" << std::endl;
	        exit( EXIT_FAILURE );
	    }
    }
} 
