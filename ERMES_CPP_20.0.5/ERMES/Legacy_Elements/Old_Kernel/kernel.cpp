/***********************************************************   
 *          
 *   Last modified by:    $Author:  $
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

  Kernel::Kernel(): mpModel(new Model()), mpModeler(new Modeler(mpModel))
  {
  }
  
  Kernel::~Kernel()
  {
       KRATOS_TRACE(" Kernel::~Kernel()","Kernel destroyed");
  }

  void Kernel::GenerateLinearSolver(const String& Name, const String& SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance)
  {
       mpModeler->GenerateLinearSolver(Name, SolverType, MaxNumberOfIterations, StepIterations, Tolerance);
  }
   
  void Kernel::CreateTimeStep()
  {
       mpModel->SetCurrentTime(mTime);
  }
  

  void Kernel::CreateSolutionStep()
  {
       mpModel->CreateWithCurrentSolutionStepData();
  }

  void Kernel::SweepFrequency(double IniFreq,double EndFreq,double StepFreq)
  {
	   mFrequencySweep = true;

	   mIniFreq  = IniFreq;
	   mEndFreq  = EndFreq;
	   mStepFreq = StepFreq;
  }

  void Kernel::SetProblemType(const String& ProblemType)
  {
       mBuild         = false;
	   mCreateHONodes = false; 
	   mPrintHOMesh   = false;

	   mFSing2L       = false;
	   mFSing3L       = false;
	   mFSing4L       = false;
	   mFSing5L       = false;
	   mFSing6L       = false;

	   mCalculateH      = false;
	   mCalculateJ      = false;
	   mCalculateJouleH = false;

	   if (ProblemType == "Build")
       { 
		    time_t bStart, bEnd;

		    if (mFrequencySweep)
			{
                mpModeler->SetFrequency(mIniFreq);

			    mBuild = true;

				std::cout << "Building linear system..." << std::endl;
				time(&bStart);
				
				String build_file_name = mBaseFileName;
				build_file_name << "-1.dat";

				GidInput inputBuild(build_file_name);
				inputBuild.Read(this);

				String bound_build_file_name = mBaseFileName;
				bound_build_file_name << "-2.dat";

				GidInput inputBoundBuild(bound_build_file_name);
				inputBoundBuild.Read(this);

				time(&bEnd);
				std::cout << "Building finished in " << difftime(bEnd, bStart) 
					      << " seconds."             << std::endl << std::endl;

				String solve_file_name = mBaseFileName;
				solve_file_name << "-5.dat";

				GidInput inputSolve(solve_file_name);
				inputSolve.Read(this);

			    mpModeler->FinishFirstStep();
                
                for(double freq = mIniFreq+mStepFreq; freq <= mEndFreq; freq += mStepFreq)
				{
					 mpModeler->SetFrequency(freq);

			         mBuild = true;

					 std::cout << "Building linear system..." << std::endl;
					 time(&bStart);

					 GidInput inputBuild(build_file_name);
					 inputBuild.Read(this);

					 GidInput inputBoundBuild(bound_build_file_name);
				     inputBoundBuild.Read(this);

					 time(&bEnd);
					 std::cout << "Building finished in " << difftime(bEnd, bStart)
						       << " seconds."             << std::endl << std::endl;

					 GidInput inputSolve(solve_file_name);
				     inputSolve.Read(this);

					 mpModeler->FinishStep();
				}
			}
			else
			{
				mBuild = true;

				std::cout << "Building linear system..." << std::endl;
				time(&bStart);

				String build_file_name = mBaseFileName;
				build_file_name << "-1.dat";
				GidInput inputBuild(build_file_name);
				inputBuild.Read(this);

				String bound_build_file_name = mBaseFileName;
				bound_build_file_name << "-2.dat";
				GidInput inputBoundBuild(bound_build_file_name);
				inputBoundBuild.Read(this);

				time(&bEnd);
				std::cout << "Building finished in " << difftime(bEnd, bStart)
					      << " seconds."             << std::endl << std::endl;

				String solve_file_name = mBaseFileName;
				solve_file_name << "-5.dat";
				GidInput inputSolve(solve_file_name);
				inputSolve.Read(this);
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
       }
	   else if (ProblemType == "CalculateH")
       {    
            mCalculateH = true;

			mpModeler->IniHDerivation();
            
		    String calH_file_name = mBaseFileName;
			calH_file_name << "-1.dat";
			GidInput inputDataFile(calH_file_name);
			inputDataFile.Read(this);

			mCalculateH = false;
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

			mpModeler->EndJCalculation();

			mCalculateJ = false;
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
	   else
       {
		    mpModeler->SetProblemType(ProblemType);
       }
  }

  void Kernel::SetFrequency(double ProblemFrequency)
  {
	  mFrequencySweep = false;
      mpModeler->SetFrequency(ProblemFrequency);
  }

  void Kernel::Solve(const String& SolvingStrategyName, const String& ElementsGroup)
  { 
      mpModeler->SolveProblem();
  }

  void Kernel::IniBuilding()
  {
      mpModeler->IniBuilding();
  }

  void Kernel::Execute(const String& ProcessName, const String& ElementsGroup)
  {
  }

  void Kernel::GenerateSolvingStrategy(const String& Name, const String& StrategyType, const String& LinearSolverName)
  {
  }

  void Kernel::Execute(const String& BaseFileName)
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
