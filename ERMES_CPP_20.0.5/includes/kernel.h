
#if !defined(KRATOS_KERNEL)
#define KRATOS_KERNEL

#include <vector>

#include "../includes/model.h"
#include "../includes/modeler.h"
#include "../includes/gid_output.h"
#include "../includes/kratos_string.h"

namespace Kratos
{
    class GidInput;
    
    class Kernel
    {
      public:
    
        typedef unsigned int IndexType;

        typedef Vector< double >                VectorType;
        typedef Vector< std::complex<double> > cVectorType;

        typedef Vector< Vector< double > >                Vector2Type;
        typedef Vector< Vector< std::complex<double> > > cVector2Type;

        typedef Vector< Vector< Vector< double > > > ResultsOnGPsType;

        Kernel() : mpModel( new Model() ), mpModeler( new Modeler( mpModel ) )
        {
        }
    
        virtual ~Kernel()
        {
            KRATOS_TRACE( " Kernel::~Kernel()", "Kernel destroyed" );
        }

        //*****************************************************************************************
        // - Set properties
        //*****************************************************************************************    
        template<class TDataType>
        void SetProperties( IndexType PropertiesId, const Variable<TDataType>& rVariable, const TDataType& Value )
        {
            mpModeler->SetProperties( PropertiesId, rVariable, Value );
        }

        //*****************************************************************************************
        // - Gets the value of a variable
        //*****************************************************************************************
        template<class TDataType>
        TDataType GetValue( const Variable<TDataType>& rVariable )
        {
            if ( rVariable.Name() == "T" ) 
            {
                return mTime;
            }
            else  
            {
                return mpModel->Value( rVariable );
            }
        }

        //*****************************************************************************************
        // - Calculates variable on nodes or Gauss points
        //*****************************************************************************************
        template<class TDataType>
        void Calculate( const Variable<TDataType>& rVariable )
        {
    	    if( mpModeler->Is_ResultsOnNodesOn() ) 
            {
                CalculateNodal( rVariable );
            }
            else 
            {
                CalculateOnGaussPoints( rVariable );
            }
        }

        //*****************************************************************************************
        // - Calculates variable on nodes
        //*****************************************************************************************
        template<class TDataType>
        void CalculateNodal( const Variable<TDataType>& rVariable )
        {
            mpModeler->Calculate_Nodal( rVariable.getKey() );
        }

        //*****************************************************************************************
        // - Calculates variable on Gauss points
        //*****************************************************************************************
        template<class TDataType>
        void CalculateOnGaussPoints( const Variable<TDataType>& rVariable )
        { 
            CalculateOnGaussPoints( rVariable.getKey() );
        }

        //*****************************************************************************************
        // - Prints variable on nodes or Gauss points
        //*****************************************************************************************
        template<class TDataType>
        void Print( const Variable<TDataType>& rVariable )
        {
            if( mpModeler->Is_ResultsOnNodesOn() ) 
            {
                PrintOnNodes( rVariable );
            }
            else   
            {
                PrintOnGaussPoints( rVariable );
            }
        }

        //*****************************************************************************************
        // - Prints variable on nodes
        //*****************************************************************************************
        template<class TDataType>
        void PrintOnNodes( const Variable<TDataType>& rVariable )
        {
            mpDefaultOutput->PrintOnNodes( rVariable );
        }

        //*****************************************************************************************
        // - Prints variable on Gauss points
        //*****************************************************************************************
        template<class TDataType>
        void PrintOnGaussPoints( const Variable<TDataType>& rVariable )
        {
            ResultsOnGPsType ResultsOnGPsVector;

            Assign_ResultsOnGPsVector( ResultsOnGPsVector, rVariable.getKey() );

            mpDefaultOutput->PrintOnGaussPoints( rVariable, ResultsOnGPsVector );
        }

        void Assign_ResultsOnGPsVector( ResultsOnGPsType& ResultsOnGPsVector, int VarKey );

        void CalculateOnGaussPoints( int VarKey );

        void Ini_Building         ();
        void End_Debuger_Mode     ();
        void Ini_Boolean_Variables();

        void Build_Linear_System    ();
        void Print_High_Order_Mesh  ();
        void Create_High_Order_Nodes();

        void Read_Nodes_ID_Coord_File();
        void Read_Nodes_Singular_File();
        void Read_Nodes_Voltages_File();

        void Read_Material_Properties_File();
        void Read_Dirichlet_Elements_File ();
        void Read_Projector_Elements_File ();

        void Read_PBC_Elements_File ();
        void Read_Solve_Problem_File();

        void Read_Cold_Plasma_Parameters_File();
        void Read_Plane_Waves_Parameters_File();

        void Read_Angular_Frequency_File();
        void Read_Complex_Frequency_File();

        void Read_Volume_Elements_File ( bool& Active_Magnitude );
        void Read_Source_Elements_File ( bool& Active_Magnitude );

        void Read_RobinBC_Elements_File( bool& Active_Magnitude );
        void Read_Contact_Elements_File( bool& Active_Magnitude );

        void Clear_Results_OnGPs();
        void Clear_Results_OnGPs( ResultsOnGPsType& ResultsOnGPs );

        void GenerateNode( IndexType NodeId, double X, double Y, double Z );

        void CreateTimeStep    ();
        void CreateSolutionStep(); 

        void Solve  ( const String& SolvingStrategyName, const String& ElementsGroup );
        void Execute( const String& BaseFileName );

        void SetProblemType( const String& ProblemType );

        void SetFrequency  ( double ProblemFrequency );
        void SweepFrequency( double IniFreq, double EndFreq, double StepFreq );

        void SetValue( int NodeId, const Variable<double>& rVariable, const double& rValue );
        void SetValue(             const Variable<double>& rVariable, const double& rValue );

        void FixDof ( IndexType NodeId, const Variable<double>& rVariable, const double&          Value );
        void FixCDof( IndexType NodeId, const Variable<double>& rVariable, const Vector<double>& vValue );

        void GenerateLinearSolver( const String& Name, const String& SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance );

        void GeneratePECElement( Vector<int>& NodesId );
        void GeneratePMCElement( Vector<int>& NodesId );
        void GenerateTECElement( Vector<int>& NodesId );
        void GeneratePBCElement( int* NodesId, unsigned int PropertiesId );

        void GenerateContactPairs( Vector<int>& CtCNodesId );

        void GenerateFarFieldElement( int* NodesId ); 

        void GenerateSourceElement( int* NodesId, unsigned int PropertiesId );
        void GenerateVolumeElement( int* NodesId, unsigned int PropertiesId );
                                          
        void GenerateVolIntgElement( int* NodesId, unsigned int PropertiesId );
        void GenerateSrfIntgElement( int* NodesId, unsigned int PropertiesId );
                                          
        void Build_Volume_Element    ( int* NodesId, unsigned int PropertiesId );
        void Print_High_Order_Mesh   ( int* NodesId, unsigned int PropertiesId );
        void Create_HONodes_InElement( int* NodesId, unsigned int PropertiesId );
                                          
        void PrintElement_1stOrder( int* NodesId, unsigned int PropertiesId );
        void PrintElement_2ndOrder( int* NodesId, unsigned int PropertiesId );

        void GenerateCoaxPortTEMElement( int* NodesId, unsigned int PropertiesId );      
        void GenerateRWPortTE10Element ( int* NodesId, unsigned int PropertiesId );
                                          
        void GenerateProjectionCoaxialTEM( int* NodesId, unsigned int PropertiesId );
        void GenerateProjectionRWPortTE10( int* NodesId, unsigned int PropertiesId );

        void GenerateGenericRobinElement  ( int* NodesId, unsigned int PropertiesId );
        void GeneratePlasmaFarFieldElement( int* NodesId, unsigned int PropertiesId ); 

        void Calculate_E();
        void Calculate_H();
        void Calculate_B();
        void Calculate_J();
        void Calculate_Static_E();

        void Show_Static_E();
        void Show_Static_J();

        void Show_Static_Joule_Heating();
        void Export_Static_Currents   ();

        void Show_Joule_Heating  ();
        void Show_Lorentz_Force  ();
        void Show_Poynting_Vector();

        void Show_B_External           ();
        void Show_Electron_Density     (); 
        void Show_Permittivity_Tensor  ();
        void Show_All_Plasma_Parameters();

        void Show_E_Parallel       ();
        void Show_E_Parallel_t     ();
        void Show_E_Perpendicular  ();
        void Show_E_Perpendicular_t();

      private:
    
        String mBaseFileName;

        Model::Pointer mpModel;

        Modeler::Pointer mpModeler;

        Output::Pointer mpDefaultOutput;

        ResultsOnGPsType mElectrostatic_E_OnGP;
        ResultsOnGPsType mElectrostatic_J_OnGP;
        ResultsOnGPsType mElectrostatic_JouleH_OnGP;

        ResultsOnGPsType mREAL_E_OnGP, mIMAG_E_OnGP, mMOD_E_OnGP;
        ResultsOnGPsType mREAL_H_OnGP, mIMAG_H_OnGP, mMOD_H_OnGP;
        ResultsOnGPsType mREAL_B_OnGP, mIMAG_B_OnGP, mMOD_B_OnGP;
        ResultsOnGPsType mREAL_J_OnGP, mIMAG_J_OnGP, mMOD_J_OnGP;

        ResultsOnGPsType mE_OnGP;
        ResultsOnGPsType mH_OnGP;
        ResultsOnGPsType mJ_OnGP;
        ResultsOnGPsType mB_OnGP;

        ResultsOnGPsType mJOULE_HEATING_OnGP;  
        ResultsOnGPsType mLORENTZ_FORCE_OnGP;
        ResultsOnGPsType mPOYNTING_VECTOR_OnGP;

        unsigned int mElementIdCounter;

        double mIniFreq ;
    	double mEndFreq ;
    	double mStepFreq;
        double mTime;

        bool mBuild;
        bool mPrintHOMesh;

        bool mCreateHONodes;
        bool mCreateContact;

        bool mCollapseBubbles;
        bool mPlasmaInterplte;

        bool mFSing2L;
        bool mFSing3L;
        bool mFSing4L;
        bool mFSing5L;
        bool mFSing6L;

		bool mSetEdgeDOFs;

        bool mCalculate_E;
        bool mCalculate_H;
        bool mCalculate_B;
        bool mCalculate_J;
        bool mCalculate_JouleH;
    
        bool mCalculate_E_Electrostatic;
        bool mCalculate_J_Electrostatic;
        bool mCalculate_JouleH_Electrostatic;
        bool mExport_J_Electrostatic;

        bool mFrequencySweep;
    }; 
} 

#endif 

