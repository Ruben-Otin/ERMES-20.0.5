
#if !defined(KRATOS_MODELER)
#define KRATOS_MODELER

#include <time.h>
#include <math.h>
#include <complex>
#include <map>

#include <memory>

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
    
        typedef std::shared_ptr<Modeler> Pointer;
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
        void Save_System_Solution_FullWave( const Variable<TDataType>& rVariable )
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
	        for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
	  	        if( ( (*itDof)->GetVariableKey() == varKey ) && ( !(*itDof)->IsFixed() ) )
                {
	  	            mpModel->Value( rVariable, *itDof ) = x_vector[ (*itDof)->EquationId() ];
                }
            }
        }
    
        //******************************************************************************************************************
        //* - Stores solution on variable ( electrostatic mode - real values )
        //******************************************************************************************************************
        template<class TDataType>
        void Save_System_Solution_Electrostatic( const Variable<TDataType>& rVariable )
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
            for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
                if( ( (*itDof)->GetVariableKey() == varKey ) && ( !(*itDof)->IsFixed() ) )
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

        void Add_NodalFieldIntegral( Vector<double>& F_real, Vector<double>& F_imag, double Intg_Ni, Vector< std::complex<double> >& VolIntg_F );
    
        void GetGrad_forSrfIntg_SmoothingOff_Electrostatic( int* NodesId );
        
        void Generate_Node( IndexType NodeId, double X, double Y, double Z );
    
        void Generate_PEC_Element( Vector<int>& LONodesId );
	    void Generate_PMC_Element( Vector<int>& LONodesId );
        void Generate_TEC_Element( Vector<int>& LONodesId );
        void Generate_PBC_Element( int* NodesId, unsigned int PropertiesId );

        void Generate_Volume_Element_Electrostatic( int* NodesId, unsigned int PropertiesId );
        void Generate_Volume_Element_ColdPlasma   ( int* NodesId, unsigned int PropertiesId );
        void Generate_Volume_Element_FullWave     ( int* NodesId, unsigned int PropertiesId );
        void Generate_Source_Element              ( int* NodesId, unsigned int PropertiesId );
    
        void Generate_Vol_Intgrl_Element_FullWave     ( int* NodesId, unsigned int PropertiesId );
        void Generate_Vol_Intgrl_Element_Electrostatic( int* NodesId, unsigned int PropertiesId );
    
        void Generate_Srf_Intgrl_Element_FullWave     ( int* NodesId, unsigned int PropertiesId );
        void Generate_Srf_Intgrl_Element_Electrostatic( int* NodesId, unsigned int PropertiesId );

        void Set_Volume_Element_Parameters( int* NodesId, Element::Pointer pElement );

        void Clear_ElecStatic_Maps();
        void Clear_Projection_Maps();
        void Clear_VolumeIntg_Maps();
        void Clear_SurfacIntg_Maps();

        void Set_PBC          ();
	    void Set_FrontBack_PBC();
	    void Set_RightLeft_PBC();

        void SetAll_Epar_ToZero();
        void SetAll_cFzo_ToZero();
	    void SetAll_cFxy_ToZero();
	    void SetAll_cFyz_ToZero();

		void Set_Contacts          ();          
		void Set_DisctAV_Nodes     ();      
	    void Set_AvgNormals_Contact();

        void Set_Node_Elements_DOFs();
        void Set_System_Matrix_Size();

        void Set_Normals    ();
	    void Set_PEC_Normals();
	    void Set_PMC_Normals();
	    void Set_TEC_Normals();

        void Set_Problem_Type  ( const String& ProblemType      );
        void Set_Frequency_Mode( bool          SweepingFreq     );
        void Set_Frequency     ( double        ProblemFrequency );

        void Fix_All_DOFs_In_Node( Node::Pointer pNode );

        void PBC_InvRot( int NodeId, double& cos_A, double& sin_A );

	    void PBC_FrontBack_TKT( Vector<int>& HONodesId, Vector<int>& IdVector, Matrix< std::complex<double> >& StiffMatrix );
	    void PBC_RightLeft_TKT( Vector<int>& HONodesId, Vector<int>& IdVector, Matrix< std::complex<double> >& StiffMatrix );
    
	    void PBC_FrontBack_TR ( Vector<int>& HONodesId, Vector<int>& IdVector, Vector< std::complex<double> >& ResVector );
	    void PBC_RightLeft_TR ( Vector<int>& HONodesId, Vector<int>& IdVector, Vector< std::complex<double> >& ResVector );

	    bool PBC_LeftNodeIsFrontNode( int rightNodeId, Vector<int>& leftNodesId, Vector<int>& FrontLeftNodesId );

        void Update_Static_Voltage();
        void Update_FullWv_Voltage();
        void UpdateAxisToCartesian();

        void UpdateCoord_Contact  ();
        void UpdateCoord_DisctAV  ();
        void UpdateCoord_PlasmaRLP();
        void UpdateCoord_EPEC     ();
	    void UpdateCoord_EPBC     ();

        void UpdateCoord_FB_PBC(                             );
	    void UpdateCoord_RL_PBC( Vector<int>& FrontLeftNodes );
	    void UpdateCoord_LF_PBC( Vector<int>& FrontLeftNodes );
    
	    void Generate_RectPort_TE10_Element( int* NodesId, unsigned int PropertiesId );
	    void Generate_CoaxPort_TEM_Element ( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void Generate_GenericRobin_Element_FullWave     ( int* NodesId, unsigned int PropertiesId );
        void Generate_GenericRobin_Element_Electrostatic( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void Generate_Projection_RectPort_TE10( int* NodesId, unsigned int PropertiesId );
	    void Generate_Projection_CoaxPort_TEM ( int* NodesId, unsigned int PropertiesId );
                                                                                         
        void Generate_FarField_Element_FullWave  ( int* NodesId                            );
        void Generate_FarField_Element_ColdPlasma( int* NodesId, unsigned int PropertiesId );

        void Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix );
    
        void FixDof ( unsigned int NodeId, const Variable<double>& rVariable, const double& Value          );
        void FixCDof( unsigned int NodeId, const Variable<double>& rVariable, const Vector<double>& vValue );
    
	    void Finish_Step      ();
        void Finish_First_Step();

        void Create_Directory( String Directory_Name );

        void Initz_VolIntgs_Fields_OnFiles();
        void Initz_SrfIntgs_Fields_OnFiles();

        void Write_VolIntgs_Fields_OnFiles();
        void Write_SrfIntgs_Fields_OnFiles();

        void Solve_Linear_System        ();        
        void Solve_Problem_FullWave     ();
        void Solve_Problem_Electrostatic();
 
	    void Ortogonalization( Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto );

        void Generate_LinearSolver_FullWave  ( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance );
        void Generate_LinearSolver_ColdPlasma( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance );

        bool Use_Direct_Solver( const String& SolverType );
        
        void TangencialCoordinates( Vector<double>& n, Vector<double>& t, Vector<double>& b );

        void Calculate_Sij_Parameters();

        void Calculate_JouleH_Element_FullWave    ( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP );
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

        void Build_Imported_Currents();
        void Build_Imported_RobinBCs();
        void Build_Imported_VolElmts();

        void Calculate_Imported_Currents();
        void Read_Imported_Nodal_Current( std::fstream& Import_J_File, Vector< Vector< std::complex<double> > >& Import_cJ_Nodal );

        void Assemble_Imported_Current_Element( int* NodesId, Vector<double>& Element_J );
        void Assemble_Imported_Current_OnNodes( Vector< Vector< std::complex<double> > >& Import_cJ_Nodal );

        void Export_Electrostatic_J_Element_GP( int* NodesId, unsigned int PropertiesId  );

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

        void Clear_Singularities();

        void Set_Edge_Elements_DOFs    ( int* NodesId );
        void Collapse_LL2P_Bubbles     ( int* NodesId );
        void Plasma_Profile_Interpolate( int* NodesId );

        void Apply_Rotation_To_Matrix( Vector<int>& HONodesId, Matrix< std::complex<double> >& StiffMatrix );
	    void Apply_Rotation_To_Vector( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector   );

        void TransT_K_T( Matrix< std::complex<double> >& K, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex );
        void TransT_B  ( Vector< std::complex<double> >& B, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex );

        void Get_Rotation_Matrix_DisctAV( Matrix< std::complex<double> >& RotMatrix , int NodeId );
        void Get_Rotation_Matrix        ( Matrix< std::complex<double> >& RotMatrix , Vector<double>& n );
    
        void ColdPlasma_TKT( Vector<int>& HONodesId, Matrix< std::complex<double> >& StiffMatrix );
        void ColdPlasma_TR ( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector   );
    
        void Rotation_Matrix_SDP_T_CC( Matrix< std::complex<double> >& T, unsigned int NodeId );
    
	    void Push_HONodes_OnVolume( int*         NodesId, Vector<int>& HONodesId );
		void Push_HONodes_OnVolume( Vector<int>& NodesId, Vector<int>& HONodesId );

	    void Push_HONodes_OnSurface( int*         NodesId, Vector<int>& HONodesId );
	    void Push_HONodes_OnSurface( Vector<int>& NodesId, Vector<int>& HONodesId );
    
	    void HONodes_OnEdge  ( int  IdNode1, int IdNode2,              Vector<int>& HONodesId );
	    void HONodes_OnFace  ( int  IdNode1, int IdNode2, int IdNode3, Vector<int>& HONodesId );
	    void HONodes_OnVolume( int* NodesId,                           Vector<int>& HONodesId );

        void Get_Volume_Element_Nodes ( int* NodesId, Vector<Node::Pointer>& pNodes );
        void Get_Surface_Element_Nodes( int* NodesId, Vector<Node::Pointer>& pNodes );

        void Get_Volume_Element_Nodes ( int* NodesId, Vector<int>& HONodesId, Vector<Node::Pointer>& pNodes );
        void Get_Surface_Element_Nodes( int* NodesId, Vector<int>& HONodesId, Vector<Node::Pointer>& pNodes );
    
	    void Reorder_HONodes_OnFace_4th( int orderingCase, Vector<unsigned int>& face, Vector<int>& HONodesId );
    
	    void Create_HONodes        ( int* NodesId );
	    void Create_HONodes_OnEdges( int* NodesId );
	    void Create_HONodes_OnFaces( int* NodesId );
	    void Create_HONodes_OnVolum( int* NodesId );
    
	    void Create_HONodes_OnEdge_2nd( int IdNode1, int IdNode2 );
	    void Create_HONodes_OnEdge_3rd( int IdNode1, int IdNode2 );
	    void Create_HONodes_OnEdge_4th( int IdNode1, int IdNode2 );
                                       
	    void Create_HONodes_OnFace_3rd( int IdNode1, int IdNode2, int IdNode3 );
	    void Create_HONodes_OnFace_4th( int IdNode1, int IdNode2, int IdNode3 );
   
	    void Print_Projection              ();
        void Print_Projection_Electrostatic();
    
	    void JoinContactNormalsInDirichletS();
	    void ForceDirichletNormalsInContact();
	    void OrtogonalizeContactNormals    ();

	    void Contact_Properties       ( int* NodesId, unsigned int PropertiesId );
	    void Generate_Contact_Pairs   ( Vector<int>& CtCNodesId );
        void Calculate_Contact_Normals( Vector<int>& NodesId    );

        void Set_TempContactPairs_E ( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L );
        void Set_TempContactPairs_A ( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L );

        void Apply_DisctAV_TKT_AV( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix );
        void Apply_Contact_TKT_AV( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix );
        void Apply_Contact_TKT_Ef( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix );

	    void DelTempContactVectors();
    
	    void Apply_Contact_Matrix_TKT( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix );
	    void Apply_Contact_Vector_TR ( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector< std::complex<double> >& ResVector   );
    
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

	    void AxisymmetricCoordinates( unsigned int NodeId, Vector<double>& n, Vector<double>& t, Vector<double>& b );
    
	    void Apply_Axisymmetric_TKT( Vector<int>& HONodesId, Matrix< std::complex<double> >& StiffMatrix );
	    void Apply_Axisymmetric_TR ( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector   );
    
        void Calculate_J_Induced_OnNodes( int* NodesId, unsigned int PropertiesId );
        void Calculate_J_Imposed_OnNodes( int* NodesId, unsigned int PropertiesId );

        void Calculate_J_Induced_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_eddy_OnGP );
        void Calculate_J_Imposed_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_imps_OnGP );

	    void LoadPlasmaParameters( String plasma_file_name );
        void LoadPWavesParameters( String pwaves_file_name );
        void LoadComplexFrequency( String cpxfrq_file_name );

        void Get_Boundary_Normal( Vector<double>& BoundaryNormal, int NodeId );
        void Get_ContactP_Normal( Vector<double>& ContactPNormal, int NodeId );

	    void Extract_ColdPlasma_ElectronDensity   ();
	    void Extract_ColdPlasma_Bexternal         ();
	    void Extract_ColdPlasma_PermittivityTensor();
        void Extract_ColdPlasma_Eparallel         ();
        void Extract_ColdPlasma_Eperpendicular    ();
	    void Extract_ColdPlasma_AllParameters     ();
    
        void Check_Consistency();
        void Print_Set_Up_Info();
    
        void Send_Error_Msg( String ErrorMsg, String CheckMsg, int VerticalSpaces );
    
        void Show_Global_Matrix();

        void Calculate_Nodal( int variableKey );

	    bool Is_OutputFileFormat_BIN  ();
	    bool Is_OutputFileFormat_ASCII();

        bool Is_Plasma_Profile_Resized();
    
        bool Is_ElectrostaticMode();
        bool Is_ColdPlasmaMode   ();
        bool Is_FullWaveMode     ();

        bool Is_PotentialsOn    ();
        bool Is_ResultsOnNodesOn();
		bool Is_EdgeElementsOn  ();
        bool Is_LL2PSmoothingOn ();

        bool Is_DisctAVNode( int NodeId );

        int Get_Element_Order();
        int Get_Number_Of_GPs();

        double Get_Problem_Frequency();

	    double VolumeOfTetrahedra( Vector<Node::Pointer>& cnodes );
	    double AreaOfTriangle    ( Vector<Node::Pointer>& cnodes );
	    double DistanceToAxis    ( unsigned int NodeId );

        std::complex<double> Determinant( Matrix< std::complex<double> >& M, int n );

        Vector<double> Calculate_Geom_Averaged_Normal( Vector< Vector<double> >& vectorSet );
        Vector<double> Calculate_Area_Weighted_Normal( Vector< Vector<double> >& vectorSet );
    
      private:

        Modeler();

        time_t mStart;
        time_t mEnd;

        Model::Pointer mpModel;

        ProblemTypeType mProblemType;  

	    Vector< std::complex<double> > x_vector;
	    Vector< std::complex<double> > b_vector;

        CSRMatrix< std::complex<double> > A_matrix;
        CSRMatrix< std::complex<double> > A_matrix_aux;
    
        std::set<Dof::Pointer, ComparePDof> mDofSet;

	    ColdPlasma::Pointer mpColdPlasma;

        ComplexSolver::Pointer mpLinearSolver;

	    String mExternalSolverPath;
	    String mBaseFileName;

        const String mDefaultInputFileName;
        const String mDefaultOutputFileName;
    
        Vector<int> mNeighbourElements;
        Vector<int> mDisctAVNodes     ;

        Vector< Vector<int> > mVolIntg_Nodes;
        Vector< Vector<int> > mSrfIntg_Nodes;

        Vector< Vector<double> > mPlaneWavesParameters;
        Vector< Vector<double> > mSurfInt_NonSmooth_E_Electrostatic;

        std::map<int, Vector< Vector<int> > > mVolIntg_Elements;
        std::map<int, Vector< Vector<int> > > mSrfIntg_Elements;

	    std::map<double, double> mVolume_m3;
        std::map<double, double> mSurfce_m2;

        std::map<double, std::complex<double> > mInputPorts;
	    std::map<double, std::complex<double> > mProjection;
        std::map<double, std::complex<double> > mNormalization;
	    std::map<double, std::complex<double> > mProjectionNormalized;
	    
	    std::map<double, Vector< std::complex<double> > > mVolIntg_E;
	    std::map<double, Vector< std::complex<double> > > mVolIntg_H;
	    std::map<double, Vector< std::complex<double> > > mVolIntg_B;
	    std::map<double, Vector< std::complex<double> > > mVolIntg_J;
        std::map<double, Vector< std::complex<double> > > mVolIntg_F;
                                 
	    std::map<double, Vector< std::complex<double> > > mSrfIntg_E;
        std::map<double, Vector< std::complex<double> > > mSrfIntg_H;
        std::map<double, Vector< std::complex<double> > > mSrfIntg_B;
        std::map<double, Vector< std::complex<double> > > mSrfIntg_J;
        std::map<double, Vector< std::complex<double> > > mSrfIntg_S;
    
        std::map<double, Vector<double> > mVolIntg_E_Electrostatic;
        std::map<double, Vector<double> > mSrfIntg_E_Electrostatic;

	    std::map<unsigned int, double> mSingular   ;
        std::map<unsigned int, double> mNSingular2L;
        std::map<unsigned int, double> mNSingular3L;
        std::map<unsigned int, double> mNSingular4L;
	    std::map<unsigned int, double> mNSingular5L;
        std::map<unsigned int, double> mNSingular6L;

        std::map<unsigned int, double               > mFix_Static_Voltage;
        std::map<unsigned int, std::complex<double> > mFix_FullWv_Current;
        std::map<unsigned int, std::complex<double> > mFix_FullWv_Voltage;
        
	    std::map<unsigned int, Vector< Vector<unsigned int> > > mEdgeHONodes;
        std::map<unsigned int, Vector< Vector<unsigned int> > > mFaceHONodes;
	    std::map<unsigned int, Vector< Vector<unsigned int> > > mVolmHONodes;
    
        std::map<unsigned int, Vector< Vector<double> > > mSetContactNormals;
    
	    std::map<unsigned int, Vector<double> > mContactNormals;
	    std::map<unsigned int, ContactPairData> mContactPairs  ;
    
	    std::map<unsigned int, unsigned short int> mTempNodeProperties;
	    std::map<unsigned int, unsigned int>       mTempContactPairs  ;
    
	    std::map<unsigned int, Vector<double> > tmp_rJ;
	    std::map<unsigned int, Vector<double> > tmp_iJ;
    
        std::map<unsigned int, Vector<double> > tmp_ist_rJ;
	    std::map<unsigned int, Vector<double> > tmp_ist_iJ;

        std::map<unsigned int, Vector<double> > mNormals   ;
        std::map<unsigned int, Vector<double> > mCtCNormals;

        std::map<unsigned int, Vector< Vector<int> > > mPBCFrontElements;
	    std::map<unsigned int, Vector< Vector<int> > > mPBCBackElements ;
    
	    std::map<unsigned int, Vector< Vector<int> > > mPBCRightElements;
        std::map<unsigned int, Vector< Vector<int> > > mPBCLeftElements ;

        std::map<unsigned int, Vector< Vector<double> > > mSetPECnormals;
	    std::map<unsigned int, Vector< Vector<double> > > mSetPMCnormals;
	    std::map<unsigned int, Vector< Vector<double> > > mSetTECnormals;

        // 'E' = PEC | 'H' = PMC | 'T' = TEC
        std::map<unsigned int, char> mType_Of_BC_Normal;
    
	    std::map<unsigned int, int> mPBC_NodeNodePairs_Front;
	    std::map<unsigned int, int> mPBC_NodeNodePairs_Right;

        std::fstream mExportCurrentFile;

        std::complex<double> mComplexFrequency;

        std::complex<double> mExp_jPha_FB;
        std::complex<double> mExp_jPha_RL;

	    unsigned int mItSolverNumThreads;

        double mProblemFrequency;
        double mGeoTolerance    ;
        double mLL2P_hk_factor  ;

        bool mPotentials    ;
        bool mAVContinuity  ;
        bool mLL2PSmoothing ;
        bool mResultsOnNodes;

        bool mComplxFreqMode;
		bool mEdgeElementsOn;
        bool mEDGE_stabilize;
        bool mRMED_stabilize;

        bool mFrequencySweep;
        bool mImportRobinBCs;
        bool mImportCurrents;
        bool mImportEleMatrx;
        bool mExportFields  ;

	    bool mNormalsGeomAvg;
	    bool mNormalsAreaWtd;
        bool mIsRightPBCTilted;

	    bool mSetAllFzoToZero;
	    bool mSetAllFxyToZero;
	    bool mAxisymmetric   ;

	    bool mColdPlasmaMode   ;
        bool mElectrostaticMode;
        bool mFullWaveMode     ;

	    bool mSolveWithExternal;
        bool mA_matrix_aux_Required;

	    bool mReleaseSolutionMode;
	    bool mDebugSolutionMode  ;
        bool mReadSolutionMode   ;

	    bool mWriteSolutionEveryStep;
	    bool mWriteSolutionFinalStep;
	    bool mReadFileInitialGuess;
    
	    bool mOutputFileFormatIsBIN  ;
	    bool mOutputFileFormatIsASCII;

	    int mSystemSize  ;
        int mResultsOnGPs;
        int mElementOrder;
        int mNonSmooth_Surf_Counter;
    };
}

#endif

