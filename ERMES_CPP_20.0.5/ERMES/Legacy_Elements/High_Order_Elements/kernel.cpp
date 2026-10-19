/************************************************************   
 *          
 *   Last modified by:    $Author: R. Otin $
 *   Date:                $Date: 2003/10/06 12:53:20 $
 *   Revision:            $Revision: 1.24 $
 *
 *************************************************************/

#include "../includes/exception.h"
#include "../includes/kernel.h"
#include "../includes/process.h"
#include "../includes/gid_input.h"
#include "../includes/modeler.h"

#include "../ERMES/Utils/MatrixFileWritting.h"

namespace Kratos
{
    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************
    Kernel::Kernel(): mpModel( new Model() ), mpModeler( new Modeler(mpModel) )
    {
    }
    
    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************
    Kernel::~Kernel()
    {
        KRATOS_TRACE(" Kernel::~Kernel()","Kernel destroyed");
    }
    
    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************
    void Kernel::GenerateLinearSolver( const String& Name, const String& SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
    {
        if ( mpModeler->Is_PlasmaMode() ) mpModeler->GenerateLinearSolver_Plasma( Name, SolverType, MaxNumberOfIterations, StepIterations, Tolerance );
        else                              mpModeler->GenerateLinearSolver       ( Name, SolverType, MaxNumberOfIterations, StepIterations, Tolerance );
    }

    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************    
    void Kernel::CreateTimeStep()
    {
        mpModel->SetCurrentTime(mTime);
    }
    
    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************
    void Kernel::CreateSolutionStep()
    {
        mpModel->CreateWithCurrentSolutionStepData();
    }
    
    //*******************************************************************************************************************
    //  
    //*******************************************************************************************************************
    void Kernel::SweepFrequency( double IniFreq, double EndFreq, double StepFreq )
    {
        mFrequencySweep = true;
    
        mpModeler->SetFrequencyMode( mFrequencySweep );
    
	    mIniFreq  = IniFreq;
	    mEndFreq  = EndFreq;
	    mStepFreq = StepFreq;
    }
    
    //*******************************************************************************************************************
    // - Main problem setting 
    //*******************************************************************************************************************
    void Kernel::SetProblemType( const String& ProblemType )
    {
	    mFSing2L = false;
	    mFSing3L = false;
	    mFSing4L = false;
	    mFSing5L = false;
	    mFSing6L = false;
    
	    mCreateHONodes = false; 
	    mPrintHOMesh   = false;
    
	    mCalculateH      = false;
	    mCalculateJ      = false;
	    mCalculateJouleH = false;
    
        mBuild             = false;
        mSetConnectivities = false;
    
        mCalculate_E_Static = false;
        mCalculate_J_Static = false;
    
        mExport_J_Static_GP = false;
    
        mCalculate_E_Static_GP = false;
        mCalculate_J_Static_GP = false;
    
        mCalculate_JouleH_Static_GP = false;
        mCalculate_JouleH_Static    = false;
    
	    if (ProblemType == "Build")
        {
            // Building linear system
            std::cout << "Building linear system..." << std::endl;
    
            // Initiate clock
            time_t bStart; time(&bStart); 
    
            // Set initial frequency in Frequency_Sweep mode
            if ( mFrequencySweep ) mpModeler->SetFrequency( mIniFreq );
    
            // Activate build mode
            mBuild = true;
    
            // Building volumetric elements
            String build_file_name = mBaseFileName; 
            build_file_name << "-1.dat";
    
            GidInput inputBuild(build_file_name);
            inputBuild.Read(this);
    
            // Building sources and boundary conditions
            String bound_build_file_name = mBaseFileName;
            bound_build_file_name << "-2.dat";
    
            GidInput inputBoundBuild(bound_build_file_name);
            inputBoundBuild.Read(this);
    
            // Building imported static currents
            mpModeler->Build_Imported_Static_Currents();
    
            // Finish build mode
            mBuild = false;
    
            // End clock
            time_t bEnd; time(&bEnd);
    
            std::cout << "Building finished in " << difftime(bEnd, bStart)
                      << " seconds."             << std::endl << std::endl;
    
            // Solving problem and writing results
            String solve_file_name = mBaseFileName;
            solve_file_name << "-5.dat";
    
            GidInput inputSolve(solve_file_name);
            inputSolve.Read(this);
    
            // If frequency sweep is activated then calculate next frequency
	        if ( mFrequencySweep )
	  	    {
	  	        mpModeler->FinishFirstStep();

                for( double freq = mIniFreq+mStepFreq; freq <= mEndFreq; freq += mStepFreq )
	  	  	    {
	  	  	        // Building linear system in frequency sweep mode
                    std::cout << "Building linear system..." << std::endl;
                  
                    // Initiate clock
	  	  	        time(&bStart);
                          
                    // Set frequency for this frequency step
                    mpModeler->SetFrequency( freq );
    
                    // Start building
                    mBuild = true;
            
                    // Building volumetric elements
                    GidInput inputBuild(build_file_name);
                    inputBuild.Read(this);
 
                    // Building sources and boundary conditions
	  	  		    GidInput inputBoundBuild(bound_build_file_name);
	  	  	        inputBoundBuild.Read(this);
                    
                    // Building imported static currents
                    mpModeler->Build_Imported_Static_Currents();
    
                    // Finish build mode
                    mBuild = false;
            
                    // End clock
	  	  		    time(&bEnd);
                    
	  	  		    std::cout << "Building finished in " << difftime(bEnd, bStart)
	  	  			          << " seconds."             << std::endl << std::endl;
                        
                    // Solving problem and writing results
	  	  		    GidInput inputSolve(solve_file_name);
	  	  	        inputSolve.Read(this);
                    
	  	            mpModeler->FinishStep();
	  	   	    }
	  	    }
        }
	    else if (ProblemType == "CreateHONodes")
        {    
            mCreateHONodes = true;
            
            String CHONodes_file_name = mBaseFileName;
            CHONodes_file_name << "-1.dat";
            GidInput inputDataFile(CHONodes_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "CreateQHONodes")
        { 
            String CQHONodes_file_name = mBaseFileName;
            CQHONodes_file_name << "-6.dat";
            GidInput inputDataFile(CQHONodes_file_name);
            inputDataFile.Read(this);
	    }
	    else if (ProblemType == "MakeContact")
        {    
            String contact_file_name = mBaseFileName;
            contact_file_name << "-3.dat";
            GidInput inputDataFile(contact_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "ImposePBC")
        {         
            String contact_file_name = mBaseFileName;
            contact_file_name << "-7.dat";
            GidInput inputDataFile(contact_file_name);
            inputDataFile.Read(this);
        }
        else if (ProblemType == "FindSing2L")
        {    
            mpModeler->IniEraseExtraNodes();
            
            mFSing2L = true;
                     
            String fsing_file_name = mBaseFileName;
            fsing_file_name << "-1.dat";
            GidInput inputDataFile(fsing_file_name);
            inputDataFile.Read(this);
            
            mpModeler->EndEraseExtraNodes();
        }
        else if (ProblemType == "FindSing3L")
        {    
            mFSing3L = true;
                  
            String fsing_file_name = mBaseFileName;
            fsing_file_name << "-1.dat";
            GidInput inputDataFile(fsing_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "FindSing4L")
        {    
            mFSing4L = true;
                  
            String fsing_file_name = mBaseFileName;
            fsing_file_name << "-1.dat";
            GidInput inputDataFile(fsing_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "FindSing5L")
        {    
            mFSing5L = true;
    
            String fsing_file_name = mBaseFileName;
            fsing_file_name << "-1.dat";
            GidInput inputDataFile(fsing_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "FindSing6L")
        {   
            mFSing6L = true;
                  
            String fsing_file_name = mBaseFileName;
            fsing_file_name << "-1.dat";
            GidInput inputDataFile(fsing_file_name);
            inputDataFile.Read(this);
        }
	    else if (ProblemType == "PrintHOMesh")
        {        
	        mpDefaultOutput->SetFileFormat( mpModeler->Is_OutputFileFormat_ASCII() );
	         
	        mpDefaultOutput->PrintMeshNodes(1);
    
            mPrintHOMesh = true;
    
            mElementIdCounter = 0;
    
            String printMesh_file_name = mBaseFileName;
	  	    printMesh_file_name << "-1.dat";
	  	    GidInput inputDataFile(printMesh_file_name);
	  	    inputDataFile.Read(this);
    
            mpDefaultOutput->CloseHOMeshFile();
    
            if ( mpModeler->Is_StaticMode() ) mpDefaultOutput->PrintGaussPointsHeader();
        }
        else if (ProblemType == "CalculateH")
        {    
            mCalculateH = true;
    
	  	    mpModeler->IniHDerivation();
             
	        String calH_file_name = mBaseFileName;
	  	    calH_file_name << "-1.dat";
	  	    GidInput inputDataFile(calH_file_name);
	  	    inputDataFile.Read(this);
   
            mpModeler->Average_H_In_Contact_Nodes();
            mpModeler->Apply_PEC_and_PMC_NBC_To_H();
    
            mCalculateH = false;
        }
        else if (ProblemType == "Show_E_Static_smoothed")
        {    
            mCalculate_E_Static = true;
    
            mpModeler->IniStatic_E_Derivation();
    
            String calE_file_name = mBaseFileName;
            calE_file_name << "-1.dat";
            GidInput inputDataFile(calE_file_name);
            inputDataFile.Read(this);
    
            mpModeler->EndStatic_E_Derivation();
                        
            mpDefaultOutput->Print(ELECTRIC_FIELD_SMOOTHED);
    
            mCalculate_E_Static = false;
        }
        else if (ProblemType == "Show_J_Static_smoothed")
        {    
            mCalculate_J_Static = true;
    
            mpModeler->IniStatic_J_Derivation();
    
            String calJ_file_name = mBaseFileName;
            calJ_file_name << "-1.dat";
            GidInput inputDataFile(calJ_file_name);
            inputDataFile.Read(this);
    
            mpModeler->EndStatic_J_Derivation();
    
            mpDefaultOutput->Print(CURRENT_DENSITY_SMOOTHED);
    
            mCalculate_J_Static = false;
        }
        else if (ProblemType == "Show_JouleHeating_Static_smoothed")
        {    
            mCalculate_JouleH_Static = true;
    
            mpModeler->IniStatic_JouleH_Derivation();
    
            String calJouleH_file_name = mBaseFileName;
            calJouleH_file_name << "-1.dat";
            GidInput inputDataFile(calJouleH_file_name);
            inputDataFile.Read(this);
    
            mpModeler->EndStatic_JouleH_Derivation();
    
            mpDefaultOutput->Print(JOULE_HEATING_SMOOTHED);
    
            mCalculate_JouleH_Static = false;
        }
        else if (ProblemType == "Show_E_Static_GP")
        {    
            mCalculate_E_Static_GP = true;
           
            String calE_file_name = mBaseFileName;
            calE_file_name << "-1.dat";
            GidInput inputDataFile(calE_file_name);
            inputDataFile.Read(this);
           
            mpDefaultOutput->PrintOnGaussPoints(ELECTRIC_FIELD_GP, mStatic_E_GaussPs);
    
            mCalculate_E_Static_GP = false;
        }
        else if (ProblemType == "Show_J_Static_GP")
        {    
            mCalculate_J_Static_GP = true;
    
            String calJ_file_name = mBaseFileName;
            calJ_file_name << "-1.dat";
            GidInput inputDataFile(calJ_file_name);
            inputDataFile.Read(this);
    
            mpDefaultOutput->PrintOnGaussPoints(CURRENT_DENSITY_GP, mStatic_J_GaussPs);
    
            mCalculate_J_Static_GP = false;
        }
        else if (ProblemType == "Export_Static_Currents")
        {
            mpModeler->Ini_Exporting_Static_Current();         
            
            mExport_J_Static_GP = true;
    
            String expJ_file_name = mBaseFileName;
            expJ_file_name << "-1.dat";
            GidInput inputDataFile(expJ_file_name);
            inputDataFile.Read(this);
    
            mExport_J_Static_GP = false;
            
            mpModeler->End_Exporting_Static_Current();      
        }
        else if (ProblemType == "Show_JouleHeating_Static_GP")
        {    
            mCalculate_JouleH_Static_GP = true;
    
            String calJouleH_file_name = mBaseFileName;
            calJouleH_file_name << "-1.dat";
            GidInput inputDataFile(calJouleH_file_name);
            inputDataFile.Read(this);
    
            mpDefaultOutput->PrintOnGaussPoints(JOULE_HEATING_GP, mStatic_JouleH_GaussPs);
    
            mCalculate_JouleH_Static_GP = false;
        }
	    else if (ProblemType == "CalculateJouleHeating")
        {    
            mCalculateJouleH = true;
             
	        String calJouleH_file_name = mBaseFileName;
	  	    calJouleH_file_name << "-1.dat";
	        GidInput inputDataFile(calJouleH_file_name);
	  	    inputDataFile.Read(this);
    
	  	    mpDefaultOutput->Print(JOULE_HEATING);
    
	  	    mCalculateJouleH = false;
        }
        else if (ProblemType == "CalculatePoyntingVector")
        {   
            mpModeler->CalculateNodalE3D(POYNTING_VECTOR.getKey());
            mpDefaultOutput->Print(POYNTING_VECTOR);   
        }
        else if (ProblemType == "CalculateLorentzForce")
        {   
            mpModeler->CalculateNodalE3D(LORENTZ_FORCE.getKey());
            mpDefaultOutput->Print(LORENTZ_FORCE);   
        }
	    else if (ProblemType == "CalculateJ")
        {    
            mCalculateJ = true;
            
            mpModeler->IniJCalculation();
                   
            String calJVE_file_name = mBaseFileName;
            calJVE_file_name << "-1.dat";
            GidInput inputVEFile(calJVE_file_name);
            inputVEFile.Read(this);
            
            String calJJE_file_name = mBaseFileName;
            calJJE_file_name << "-2.dat";
            GidInput inputJEFile(calJJE_file_name);
            inputJEFile.Read(this);
    
            mpModeler->Calculate_Imported_Static_Currents(); 
            
            mpModeler->EndJCalculation();
            
            mCalculateJ = false;
        }
	    else if (ProblemType == "Show_Electron_Density")
	    {
            mpModeler->Extract_ColdPlasma_ElectronDensity();
            mpDefaultOutput->Print(ELECTRON_DENSITY);
	    }
	    else if (ProblemType == "Show_B_External")
	    {
            mpModeler->Extract_ColdPlasma_Bexternal();
            mpDefaultOutput->Print(B_EXT);
	    }
	    else if (ProblemType == "Show_Permittivity_Tensor")
	    {
	        mpModeler->Extract_ColdPlasma_PermittivityTensor();
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_S); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_D); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_P); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_R); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_L); 
	    }
        else if (ProblemType == "Show_E_Parallel")
        {
            mpModeler->Extract_ColdPlasma_Eparallel();
            mpDefaultOutput->Print(MOD_E_PARALLEL);   
        }
        else if (ProblemType == "Calculate_Eparallel_t")
        {   
            mpModeler->CalculateNodalE3D(E_PARALLEL_T.getKey());
            mpDefaultOutput->Print(E_PARALLEL_T);   
        }
        else if (ProblemType == "Show_E_Perpendicular")
        {
            mpModeler->Extract_ColdPlasma_Eperpendicular();
            mpDefaultOutput->Print(MOD_E_PERPENDICULAR);   
        }
        else if (ProblemType == "Calculate_Eperpendicular_t")
        {   
            mpModeler->CalculateNodalE3D(E_PERPENDICULAR_T.getKey());
            mpDefaultOutput->Print(E_PERPENDICULAR_T);   
        }
	    else if (ProblemType == "Show_All_Plasma_Parameters")
	    {
	        mpModeler->Extract_ColdPlasma_AllParameters();
	        mpDefaultOutput->Print(ELECTRON_DENSITY);
	        mpDefaultOutput->Print(B_EXT);
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_S); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_D); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_P); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_R); 
	        mpDefaultOutput->Print(PLASMA_PERMITTIVITY_L); 
	    }
        else if (ProblemType == "Set_E_parallel_zeros")
        {
            mpModeler->Set_E_parallel_zeros();
        }
	    else if (ProblemType == "EndDebugMode")
	    {  
            std::cout << std::endl;
            std::cout << "ERMES (debug mode) finished." << std::endl;
            std::cout << std::endl;
            
            mpDefaultOutput->~Output();
            exit(EXIT_SUCCESS);
	    }
	    else if (ProblemType == "ReadSolutionFromFile")
	    {
	        String solve_file_name = mBaseFileName;
	        solve_file_name << "-5.dat";
	        GidInput inputSolve(solve_file_name);
	        inputSolve.Read(this);
	    }
	    else if (ProblemType == "LoadPlasmaParameters")
	    {
            String plasma_file_name = mBaseFileName;
            plasma_file_name << "-8.dat";
            mpModeler->LoadPlasmaParameters(plasma_file_name);
	    }
	    else
        {
	        mpModeler->SetProblemType(ProblemType);
        }
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::SetFrequency( double ProblemFrequency )
    {
        mFrequencySweep = false;

        mpModeler->SetFrequencyMode( mFrequencySweep  );
        mpModeler->SetFrequency    ( ProblemFrequency );
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::Solve( const String& SolvingStrategyName, const String& ElementsGroup )
    { 
        if ( mpModeler->Is_StaticMode() ) mpModeler->SolveStaticProblem(); 
        else                              mpModeler->SolveProblem      ();
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::IniBuilding()
    {
        mpModeler->IniBuilding();
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::Execute( const String& ProcessName, const String& ElementsGroup )
    {
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::GenerateSolvingStrategy( const String& Name, const String& StrategyType, const String& LinearSolverName )
    {
    }
    
    //*******************************************************************************************************************
    // 
    //*******************************************************************************************************************
    void Kernel::Execute( const String& BaseFileName )
    {
        mBaseFileName = BaseFileName;
	    mpModeler->SetBaseFileName(mBaseFileName);
    
	    String results_file_name = BaseFileName;
	    mpDefaultOutput = Output::Pointer(new GidOutput(mpModel, results_file_name));
    
        String nodes_file_name = BaseFileName;
        nodes_file_name << ".dat";
        GidInput input(nodes_file_name);
    
        std::cout << "Reading nodes..." << std::endl;  
         
	    try
	    {
	        input.Read(this);
	    }
	    catch (...)
	    {
	        std::cout << std::endl;
	        std::cout << "!!!!!!!!!!!! Check materials !!!!!!!!!!!!!" << std::endl;
	        std::cout << std::endl;
	        std::cout << "ERMES will not continue building.         " << std::endl;
	        std::cout << "Please, check materials and start again..." << std::endl;
	        std::cout << std::endl;
	        std::cout << "ERMES analysis finished.                  " << std::endl;
	        exit(0);
	    }
    }
} 
