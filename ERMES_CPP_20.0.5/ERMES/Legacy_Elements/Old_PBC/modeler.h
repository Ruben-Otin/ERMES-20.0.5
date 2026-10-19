
#if !defined(KRATOS_MODELER)
#define KRATOS_MODELER

#include <time.h>
#include <math.h>
#include <complex>
#include <map>

#include "../external_libraries/boost/boost/smart_ptr.hpp"

#include "../includes/model.h"
#include "../includes/kratos_string.h"
#include "../includes/exception.h"
#include "../includes/point.h"
#include "../includes/process.h"

#include "../linear_solvers/complex_solver.h"

#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class Modeler
    {
      public:
    
        enum ProblemTypeType {E3D};

        typedef unsigned int IndexType;
    
        typedef boost::shared_ptr<Modeler> Pointer;
        typedef Vector<Node::Pointer>      NodesArrayType;

        typedef Vector< double >                VectorType;
        typedef Vector< std::complex<double> > cVectorType;

        typedef Vector< Vector< double > >                Vector2Type;
        typedef Vector< Vector< std::complex<double> > > cVector2Type;

        typedef Vector< Vector< Vector< double > > > ResultsOnGPsType;

        struct ContactPairData
	    {
	  	    short int matR;
	  	    short int matL;
	  	    unsigned int ctcNode;
	  	    Vector<double> ctcNormal;
	    };
    
        //******************************************************************************************************************
        //* - Modeler constructor
        //******************************************************************************************************************
        Modeler( Model::Pointer pModel ) : mpModel( pModel )
        {
            mNonSmooth_Surf_Counter = 0;
        }
    
        //******************************************************************************************************************
        //* - Modeler destructor 
        //******************************************************************************************************************
        virtual ~Modeler()
        {
        }
    
        //******************************************************************************************************************
        //* - Set properties
        //******************************************************************************************************************
        template<class TDataType>
	    void SetProperties( IndexType PropertiesId, const Variable<TDataType>& rVariable, const TDataType& Value )
        {
	        if ( mpModel->GetProperties( PropertiesId ).get() == 0 )
            {
	            mpModel->AddProperties( PropertiesId, Properties::Pointer( new Properties( *mpModel ) ) );
            }
    
	        typename PropertyFunction<TDataType>::Pointer constant_property( new ConstantProperty<TDataType>( Value ) );
    
	        mpModel->GetProperties( PropertiesId )->SetProperty( rVariable, constant_property );
        }
    
        //******************************************************************************************************************
        //* - Stores solution on variable ( full wave mode - complex values )
        //******************************************************************************************************************
        template<class TDataType>
        void SaveSystemSolution_FullWave( const Variable<TDataType>& rVariable )
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
	        for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
	  	        if ( ( (*itDof)->GetVariableKey() == varKey ) && ( !(*itDof)->IsFixed() ) )
                {
	  	            mpModel->Value( rVariable, *itDof ) = x_vector[ (*itDof)->EquationId() ];
                }
            }
        }
    
        //******************************************************************************************************************
        //* - Stores solution on variable ( electrostatic mode - real values )
        //******************************************************************************************************************
        template<class TDataType>
        void SaveSystemSolution_Electrostatic( const Variable<TDataType>& rVariable )
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
            for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
                if ( ( (*itDof)->GetVariableKey() == varKey ) && ( !(*itDof)->IsFixed() ) )
                {
                    mpModel->Value( rVariable, *itDof ) = std::real( x_vector[ (*itDof)->EquationId() ] );
                }
            }
        }

        //******************************************************************************************************************
        //* - Set base file name
        //******************************************************************************************************************
	    void SetBaseFileName( String baseFileName ) 
        { 
            mBaseFileName = baseFileName; 
        }

        void GenerateNode( IndexType NodeId, double X, double Y, double Z );
    
        void GeneratePECElement  ( std::vector<int>& LONodesId );
        void GeneratePSBCElement ( std::vector<int>& LONodesId );
	    void GeneratePMCElement  ( std::vector<int>& LONodesId );
        void GenerateTEPMCElement( std::vector<int>& LONodesId );
    
        Vector<double> Calculate_Geom_Averaged_Normal( std::vector< std::vector<double> >& vectorSet );
        Vector<double> Calculate_Area_Weighted_Normal( std::vector< std::vector<double> >& vectorSet );
    
        void GetGrad_forSurfInt_SmoothingOff_Electrostatic( int* NodesId );
    
        void GenerateVEQElement( std::vector<int>& QHONodesId );
    
        void GenerateVolumeElement_Electrostatic ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolumeElement_ColdPlasma    ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolumeElement_FullWave      ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateSourceElement               ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateVolIntElement               ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolIntElement_Electrostatic ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateSurfIntElement              ( int* NodesId, unsigned int PropertiesId );
        void GenerateSurfIntElement_Electrostatic( int* NodesId, unsigned int PropertiesId );
    
        void GeneratePBCElement                  ( int* NodesId, unsigned int PropertiesId );

        void Set_Volume_Element_Parameters       ( int* NodesId, Element::Pointer pElement );

        void SetPBC              ();
	    void Set_FrontBack_PBC   ();
	    void Set_RightLeft_PBC   ();

        void Set_E_parallel_zeros();
        void SetAllcEzToCero     ();
	    void SetAllcExyToCero    ();
	    void SetAllcEyzToCero    ();

        void SetDofSet ();
        void SetNormals();
    
	    void SettingPECNormals  ();
        void SettingPSBCNormals ();
	    void SettingPMCNormals  ();
	    void SettingTEPMCNormals();

        void SetProblemType  ( const String& ProblemType      );
        void SetFrequencyMode( bool          SweepingFreq     );
        void SetFrequency    ( double        ProblemFrequency );

        int    PBC_Node_To_Node( double natX, double natY );
        double PBC_Node_To_Edge( double natX, double natY, std::vector<int>& ElementNodeIds, std::vector<int>& EdgeNodeIds );

        double PBC_RotateY( Node::Pointer prNode );
        void   PBC_InvRot ( int NodeId, double& cos_A, double& sin_A );

	    void PBC_FrontBack_TKT(std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_RightLeft_TKT(std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
    
	    void PBC_FrontBack_TR (std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
	    void PBC_RightLeft_TR (std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);

        void PBC_FB_NodeToNode_cID   (int NodeId, int  local_i, std::vector<int>& HONodesId, Vector<int>& IdVector);
	    void PBC_FB_NodeToEdge_TKT   (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_FB_NodeToElement_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
    
	    void PBC_RL_NodeToNode_TKT   (int NodeId, int  local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_RL_NodeToEdge_TKT   (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_RL_NodeToElement_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
    
	    void PBC_FB_NodeToEdge_TR    (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
	    void PBC_FB_NodeToElement_TR (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
    
	    void PBC_RL_NodeToNode_TR    (int NodeId, int  local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
	    void PBC_RL_NodeToEdge_TR    (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
	    void PBC_RL_NodeToElement_TR (int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);
    
	    void PBC_IdVectorReDef       (int local_i, std::vector<Node::Pointer>& pNodes, std::vector<int>& HONodesId, Vector<int>& IdVector);
    
	    void PBC_StiffMatrixReDef    (int local_i, int oldSize, int newSize, std::vector<double>& N, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_ResVectorReDef      (int local_i, int oldSize, int newSize, std::vector<double>& N, Vector<std::complex<double> >& ResVector  );
    
	    void PBC_StiffMatrixReDef_Rot(int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Matrix<std::complex<double> >& StiffMatrix);
	    void PBC_ResVectorReDef_Rot  (int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Vector<std::complex<double> >& ResVector);
    
	    bool PBC_LeftNodeIsFrontNode (int rightNodeId, std::vector<int>& leftNodesId, std::vector<int>& FrontLeftNodesId);
    
        void PBC_NaturalBase_N_Surface( std::vector<double>& N, std::vector<double>& natCoord );
        void PBC_NaturalBase_N_Line   ( std::vector<double>& N,             double   natCoord );
        void PBC_Mapping_Parameters   ( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY );
    
        void Fix_All_DOFs_In_Node             ( Node::Pointer pNode );
        void Assign_Front_Node_To_Closest_Edge( int NodeId, std::vector< std::vector<int> >& ElementSet );
        void Assign_Right_Node_To_Closest_Edge( int NodeId, std::vector< std::vector<int> >& ElementSet );

        void UpdateFixedVOLTAGE   ();
        void UpdateAxisToCartesian();

        void UpdateCoord_Contact  ();
        void UpdateCoord_PlasmaRLP();
        void UpdateCoord_EPEC     ();
	    void UpdateCoord_EPBC     ();

        void UpdateCoord_FB_PBC(                                  );
	    void UpdateCoord_RL_PBC( std::vector<int>& FrontLeftNodes );
	    void UpdateCoord_LF_PBC( std::vector<int>& FrontLeftNodes );
    
	    void GenerateRWPortTE10Element ( int* NodesId, unsigned int PropertiesId );
	    void GenerateCoaxPortTEMElement( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void GenerateGenericRobinElement_FullWave     ( int* NodesId, unsigned int PropertiesId );
        void GenerateGenericRobinElement_Electrostatic( int* NodesId, unsigned int PropertiesId );
        void GenerateGenericRobinElement_ColdPlasma   ( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void GenerateProjectionRWPortTE10 ( int* NodesId, unsigned int PropertiesId );
	    void GenerateProjectionCoaxPortTEM( int* NodesId, unsigned int PropertiesId );
                                                                                         
        void GenerateFarFieldElement           ( int* NodesId                            );
        void GenerateFarFieldElement_ColdPlasma( int* NodesId, unsigned int PropertiesId );

        void Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix );
    
        void FixDof ( unsigned int NodeId, const Variable<double>& rVariable, const double& Value               );
        void FixCDof( unsigned int NodeId, const Variable<double>& rVariable, const std::vector<double>& vValue );
    
        void Ini_Building   ();
	    void FinishStep     ();
        void FinishFirstStep();

        void SolveLinearSystem         ();        
        void SolveProblem_FullWave     ();
        void SolveProblem_Electrostatic();
 
	    void Ortogonalization( Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto );
    
        void GenerateLinearSolver_ColdPlasma( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance );
        void GenerateLinearSolver           ( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance );
    
        void TangencialCoordinates( Vector<double>& n, Vector<double>& t, Vector<double>& b );

        void CalculateSijParameters();

        void Calculate_Joule_Heating_Element      ( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP );
        void Calculate_Joule_Heating_OnGaussPoints( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP );
	    void Calculate_Joule_Heating_OnNodes      ( int* NodesId, unsigned int PropertiesId );

        void Calculate_E_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_E_OnGP, ResultsOnGPsType& IMAG_E_OnGP );
	    void Calculate_H_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_H_OnGP, ResultsOnGPsType& IMAG_H_OnGP );
        void Calculate_B_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_B_OnGP, ResultsOnGPsType& IMAG_B_OnGP );
        void Calculate_J_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_J_OnGP, ResultsOnGPsType& IMAG_J_OnGP );

        void Split_Complex_VectorVector( Vector2Type& real_VV, Vector2Type& imag_VV, cVector2Type& complex_VV );

        void Apply_PEC_and_PMC_To_H();
        void Apply_PEC_and_PMC_To_E();
    
        void Ini_Electrostatic_E_Derivation     ();
        void Ini_Electrostatic_J_Derivation     ();
        void Ini_Electrostatic_JouleH_Derivation();
    
        void End_Electrostatic_E_Derivation     ();
        void End_Electrostatic_J_Derivation     ();
        void End_Electrostatic_JouleH_Derivation();
                                                                   
        void Calculate_E_Element_Electrostatic     ( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs );
        void Calculate_J_Element_Electrostatic     ( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs );
        void Calculate_JouleH_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs );
    
        void Ini_Exporting_Electrostatic_Current();       
        void End_Exporting_Electrostatic_Current();

        void Build_Imported_Electrostatic_Currents    ();
        void Calculate_Imported_Electrostatic_Currents();

        void Export_Electrostatic_J_Element_GP( int* NodesId, unsigned int PropertiesId );
        void Assemble_Imported_Current_Element( int* NodesId, double* Jvec, double Jpha );

	    void Ini_E_Derivation();
        void End_E_Derivation();
            
	    void Ini_H_Derivation();
        void End_H_Derivation();
                             
        void Ini_B_Derivation();
        void End_B_Derivation();

	    void Ini_J_Calculation();
	    void End_J_Calculation();

	    void Ini_Joule_Heating_Calculation();
	    void End_Joule_Heating_Calculation();

	    void FindSing2L( int* NodesId );
	    void FindSing3L( int* NodesId );
	    void FindSing4L( int* NodesId );
	    void FindSing5L( int* NodesId );
	    void FindSing6L( int* NodesId );

        void Apply_Rotation_To_Matrix( Vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix );
	    void Apply_Rotation_To_Vector( Vector<int>& HONodesId, Vector<std::complex<double> >& ResVector   );

        void TransT_K_T( Matrix< std::complex<double> >& K, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex );
        void TransT_B  ( Vector< std::complex<double> >& B, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex );

        void Get_En_R_invEp_Rt_Dn( Matrix< std::complex<double> >& R_invEp_Rt, int NodeId );
        void Get_Rotation_Matrix ( Matrix< std::complex<double> >& RotMatrix, Vector<double>& n );
    
        void ColdPlasma_TKT( std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix );
        void ColdPlasma_TR ( std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector   );
    
        void Rotation_Matrix_SDP_T_CC(Matrix<std::complex<double> >& T, unsigned int NodeId);
    
	    void PushHONodesOnVolume ( int* NodesId             , std::vector<int>& HONodesId );
	    void PushHONodesOnSurface( int* NodesId             , std::vector<int>& HONodesId );
	    void PushHONodesOnSurface( std::vector<int>& NodesId, std::vector<int>& HONodesId );
    
	    void HONodesOnEdge  ( int IdNode1, int IdNode2,              std::vector<int>& HONodesId );
	    void HONodesOnFace  ( int IdNode1, int IdNode2, int IdNode3, std::vector<int>& HONodesId );
	    void HONodesOnVolume( int* NodesId,                          std::vector<int>& HONodesId );

        void Get_Volume_Element_Nodes ( int* NodesId, std::vector<Node::Pointer>& pNodes );
        void Get_Surface_Element_Nodes( int* NodesId, std::vector<Node::Pointer>& pNodes );

        void Get_Volume_Element_Nodes ( int* NodesId, std::vector<int>& HONodesId, std::vector<Node::Pointer>& pNodes );
        void Get_Surface_Element_Nodes( int* NodesId, std::vector<int>& HONodesId, std::vector<Node::Pointer>& pNodes );
    
	    void ReorderHONodesOnFace_4th( int orderingCase, std::vector<unsigned int>& face, std::vector<int>& HONodesId );
    
	    void CreateHONodes         ( int* NodesId );
	    void CreateHONodes_OnEdges ( int* NodesId );
	    void CreateHONodes_OnFaces ( int* NodesId );
	    void CreateHONodes_OnVolume( int* NodesId );
    
	    void CreateHONodes_OnEdge_2nd( int IdNode1, int IdNode2 );
	    void CreateHONodes_OnEdge_3th( int IdNode1, int IdNode2 );
	    void CreateHONodes_OnEdge_4th( int IdNode1, int IdNode2 );
                                       
	    void CreateHONodes_OnFace    ( int IdNode1, int IdNode2, int IdNode3 );
	    void CreateHONodes_OnFace_4th( int IdNode1, int IdNode2, int IdNode3 );
   
	    void PrintProjection              ();
        void PrintProjection_Electrostatic();
    
	    void SetContacts           ();
	    void SetAvgNormalsInContact();
    
	    void JoinContactNormalsInDirichletS();
	    void ForceDirichletNormalsInContact();
    
	    void OrtogonalizeContactNormals();

	    void ContactProperties      ( int* NodesId, unsigned int PropertiesId );
        void CalculateContactNormals( Vector<int>& NodesId );

	    void GenerateContactPairs           ( Vector<int>& CtCNodesId );
        void Ini_Contact_Smoothing_Vol      ( Vector<int>& CtCNodesId );
        void Calculate_Contact_Smoothing_Srf( Vector<int>& CtCNodesId );
        void Calculate_Contact_Smoothing_Vol( int* NodesId, unsigned int PropertiesId );

        void Set_TempContactPairs_E( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L );
        void Set_TempContactPairs_A( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L );

        void Apply_Contact_TKT_AV( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix );
        void Apply_Contact_TKT_Ef( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix );

        void End_Contact_Smooting ();
	    void DelTempContactVectors();
    
	    void Apply_Contact_Matrix_TKT( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix );
	    void Apply_Contact_Vector_TR ( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector<std::complex<double> >& ResVector   );
    
        void Apply_Contact_Matrix_TKT_ColdPlasma( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix );
        void Apply_Contact_Vector_TR_ColdPlasma ( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector<std::complex<double> >& ResVector   );

        void Get_Contact_Matrix_ColdPlasma( Matrix< std::complex<double> >& CtCMatrix, int RNodesId );
        void Get_Contact_Matrix_FullWave  ( Matrix< std::complex<double> >& CtCMatrix, int RNodesId );

        void Assemble_ElementMatrix_to_Global_Symmetric       ( Vector<int>& EleIdVector, Matrix              <double>  & EleStiffMatrix );
        void Assemble_ElementMatrix_to_Global_Symmetric       ( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix );
        void Assemble_ElementMatrix_to_Global_NonSymmetric    ( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix );
        void Assemble_ElementMatrix_to_Global_Aux_Symmetric   ( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix );
        void Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix );

        void Assemble_ResidualVector_to_Global( Vector<int>& EleIdsVector, Vector< std::complex<double> >& EleResVector );
        void Assemble_ResidualVector_to_Global( Vector<int>& EleIdsVector, Vector              <double>  & EleResVector );

        void Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix );
        void Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Vector< std::complex<double> >& EleResVector   );

	    double VolumeOfTetrahedra( std::vector<Node::Pointer>& cnodes );
	    double AreaOfTriangle    ( std::vector<Node::Pointer>& cnodes );

	    double DistanceToAxis         ( unsigned int NodeId );
	    void   AxisymmetricCoordinates( unsigned int NodeId, Vector<double>& n, Vector<double>& t, Vector<double>& b );
    
	    void Axisym_TKT( std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix );
	    void Axisym_TR ( std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector   );
    
        void Calculate_J_Induced_OnNodes( int* NodesId, unsigned int PropertiesId );
        void Calculate_J_Imposed_OnNodes( int* NodesId, unsigned int PropertiesId );

        void Calculate_J_Induced_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_eddy_OnGP );
        void Calculate_J_Imposed_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_imps_OnGP );

	    void Calculate_detJ( std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<Node::Pointer>& cnodes );
	    void Calculate_detJ( std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY,                          std::vector<Node::Pointer>& cnodes );
    
	    void NaturalDerivatives3D_2ndQ( Matrix<double>& dNda, Matrix<double>& dNdb, Matrix<double>& dNdu, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ );
    
	    void LoadPlasmaParameters( String plasma_file_name );
    
	    bool Is_OutputFileFormat_BIN  ();
	    bool Is_OutputFileFormat_ASCII();
    
        bool Is_ElectrostaticMode();
        bool Is_ColdPlasmaMode   ();
        bool Is_FullWaveMode     ();

        bool Is_QuadraticMesh();
        bool Is_SmoothingOn  ();
        bool Is_PotentialsOn ();

        int    Get_Element_Order    ();
        double Get_Problem_Frequency();

        void Get_Boundary_Normal( Vector<double>& BoundaryNormal, int NodeId );
        void Get_ContactP_Normal( Vector<double>& ContactPNormal, int NodeId );

	    void Extract_ColdPlasma_ElectronDensity   ();
	    void Extract_ColdPlasma_Bexternal         ();
	    void Extract_ColdPlasma_PermittivityTensor();
        void Extract_ColdPlasma_Eparallel         ();
        void Extract_ColdPlasma_Eperpendicular    ();
	    void Extract_ColdPlasma_AllParameters     ();
    
        void CheckConsistency();
    
        void Send_Error_Msg( String ErrorMsg, String CheckMsg, int VerticalSpaces );
    
        void Show_Global_Matrix();

        void CalculateNodal( int variableKey );

        std::complex<double> Determinant( Matrix< std::complex<double> >& M, int n );
    
      private:

        Model::Pointer mpModel;
    
        CSRMatrix<std::complex<double> > A_matrix;
        CSRMatrix<std::complex<double> > A_matrix_aux;
    
        bool mA_matrix_aux_Required;
    
	    Vector<std::complex<double> > x_vector;
	    Vector<std::complex<double> > b_vector;
    
        time_t mStart;
        time_t mEnd;
    
        ProblemTypeType mProblemType;
    
        double mProblemFrequency;
    
        std::set<Dof::Pointer, ComparePDof> mDofSet;
    
        std::map<unsigned int, std::vector<double> > mNormals;
    
        // 'E' = PEC normal | 'D' = Plasma Sheath | 'H' = PMC | 'T' = TEPMC
        std::map<unsigned int, char> mType_Of_BC_Normal;
    
        std::map<unsigned int, std::vector<std::vector<double> > > mSetPECnormals  ;
        std::map<unsigned int, std::vector<std::vector<double> > > mSetPSBCnormals ;
	    std::map<unsigned int, std::vector<std::vector<double> > > mSetPMCnormals  ;
	    std::map<unsigned int, std::vector<std::vector<double> > > mSetTEPMCnormals;
    
        ComplexSolver::Pointer mpLinearSolver;
    
	    ColdPlasma::Pointer mpColdPlasma;
    
	    String mExternalSolverPath;
	    String mBaseFileName;
    
	    bool mSolveWithExternal;
    
	    bool mNormalsGeomAvg;
	    bool mNormalsAreaWtd;
    
        Vector<int> mNeighbourElements;
    
        Modeler();
    
        const String mDefaultInputFileName;
        const String mDefaultOutputFileName;
    
        Process::Pointer mpTransientLoop;
        Process::Pointer mpAlphaLoop;
    
	    int mSystemSize;

	    std::map<double, std::complex<double> > mProjection;
	    std::map<double, std::complex<double> > mProjectionNormalized;
	    std::map<double, std::complex<double> > mNormalization;
    
	    std::map<double, Vector<std::complex<double> > > mVolumeIntegralE;
	    std::map<double, Vector<std::complex<double> > > mVolumeIntegralH;
	    std::map<double, Vector<std::complex<double> > > mVolumeIntegralJ;
    
	    std::map<double, Vector<std::complex<double> > > mSurfaceIntegralE;
        std::map<double, Vector<std::complex<double> > > mSurfaceIntegralH;
        std::map<double, Vector<std::complex<double> > > mSurfaceIntegralJ;
    
        std::map<double, Vector<double> > mVolIntegralE_Electrostatic;

        std::map<double, Vector<double> > mSurfIntegralE_Electrostatic;
   
        int mNonSmooth_Surf_Counter;
    
        Vector< Vector<double> > mSurfInt_NonSmooth_E_Electrostatic;
    
	    std::map<double, double> mVolumeValue;
        std::map<double, double> mSurfaceValue;
    
	    int mElementOrder;
    
	    std::map<unsigned int, double> mSingular;
        std::map<unsigned int, double> mNSingular2L;
        std::map<unsigned int, double> mNSingular3L;
        std::map<unsigned int, double> mNSingular4L;
	    std::map<unsigned int, double> mNSingular5L;
        std::map<unsigned int, double> mNSingular6L;
    
        std::map<unsigned int, double> mFixVoltage;

        std::map<unsigned int, bool> mContactNode;
        std::map<unsigned int, int > mContactNeighbours;
	    std::map<unsigned int, Vector<std::complex<double> > > mContactEfield;
   
	    std::map<unsigned int, std::vector<std::vector<unsigned int> > > mEdgeHONodes;
        std::map<unsigned int, std::vector<std::vector<unsigned int> > > mFaceHONodes;
	    std::map<unsigned int, std::vector<std::vector<unsigned int> > > mVolmHONodes;
    
        std::map<unsigned int, std::vector<std::vector<double> > > mSetContactNormals;
    
	    std::map<unsigned int, std::vector<double> > mContactNormals;
	    std::map<unsigned int, ContactPairData>      mContactPairs  ;
    
	    std::map<unsigned int, unsigned short int> mTempNodeProperties;
	    std::map<unsigned int, unsigned int>       mTempContactPairs  ;
    
	    std::map<unsigned int, std::vector<double> > tmp_rJ;
	    std::map<unsigned int, std::vector<double> > tmp_iJ;
    
        std::map<unsigned int, std::vector<double> > tmp_ist_rJ;
	    std::map<unsigned int, std::vector<double> > tmp_ist_iJ;
    
        std::map<unsigned int, std::vector<std::vector<int> > > mPBCFrontElements;
	    std::map<unsigned int, std::vector<std::vector<int> > > mPBCBackElements ;
    
	    std::map<unsigned int, std::vector<std::vector<int> > > mPBCRightElements;
        std::map<unsigned int, std::vector<std::vector<int> > > mPBCLeftElements ;
    
	    std::map<unsigned int, int> mPBC_NodeNodePairs_Front;
	    std::map<unsigned int, int> mPBC_NodeNodePairs_Right;
    
        std::map<unsigned int, std::vector<int> > mPBC_NodeEdgePairs_Front;
	    std::map<unsigned int, std::vector<int> > mPBC_NodeEdgePairs_Right;
    
        std::map<unsigned int, std::vector<int> > mPBC_NodeElementPairs_Front;
        std::map<unsigned int, std::vector<int> > mPBC_NodeElementPairs_Right;
    
	    std::map<unsigned int, double> mPBC_NodeEdgeCoord_Front;
	    std::map<unsigned int, double> mPBC_NodeEdgeCoord_Right;
    
	    std::map<unsigned int, std::vector<double> > mPBC_NodeElementCoord_Front;
	    std::map<unsigned int, std::vector<double> > mPBC_NodeElementCoord_Right;

        std::fstream mExportCurrentFile;

	    unsigned int mItSolverNumThreads;

	    double mPBCTolerance  ;
        double mGiDTolerance  ;
        double mLL2P_hk_factor;

	    bool mIsRightPBCTilted;

	    bool mQuadraticGeometry;
    
	    bool mSetAllEzToZero ;
	    bool mSetAllExyToZero;
	    bool mAxisymmetric   ;
    
	    bool mColdPlasmaMode   ;
        bool mElectrostaticMode;
        bool mFullWaveMode     ;

        bool mPotentials_On  ;
        bool mSmoothing_On   ;
        bool mAVContinuity_On;

        bool mFrequencySweep;
    
	    bool mReadFileInitialGuess;
    
	    bool mReleaseSolutionMode;
	    bool mDebugSolutionMode  ;
        bool mReadSolutionMode   ;
    
	    bool mWriteResultEveryStep;
	    bool mWriteResultFinalStep;
    
	    bool mOutputFileFormatIsBIN  ;
	    bool mOutputFileFormatIsASCII;
    
        bool mImportCurrents;

        int mResultsOnGPs;
    };
}

#endif

