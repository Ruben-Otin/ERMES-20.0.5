
#if !defined(KRATOS_MODELER)
#define  KRATOS_MODELER

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
    
        typedef boost::shared_ptr<Modeler> Pointer;
    
        typedef Vector<Node::Pointer> NodesArrayType;
    
        typedef unsigned int IndexType;
    
        //*****************************************************************************
        //* - Modeler constructor
        //*****************************************************************************
        Modeler( Model::Pointer pModel ) : mpModel( pModel )
        {
            mMaxDistToSing = 0.00;
	        mMinDistToSing = 0.00;
    
	        mMaxDistNode1 = 0; mMaxDistNode2 = 0;
	        mMinDistNode1 = 0; mMinDistNode2 = 0;
    
            mNonSmooth_Surf_Id = 0;
        }
    
        //*****************************************************************************
        //* - Modeler destructor 
        //*****************************************************************************
        virtual ~Modeler()
        {
        }
    
        //*****************************************************************************
        //* - Generate nodes
        //*****************************************************************************
        void GenerateNode(IndexType NodeId, double X, double Y, double Z)
        {
            Node::Pointer new_node(new Node(NodeId,X,Y,Z));
    
            mpModel->AddNode( new_node );
    
            if ( mStaticMode == true )
            {
                new_node->AddVDof( mDofSet );
            }
            else
            {
                new_node->AddDofs( mDofSet );
            }
        }
    
        //*****************************************************************************
        //* - Set properties
        //*****************************************************************************
        template<class TDataType>
	    void SetProperties( IndexType PropertiesId, const Variable<TDataType>& rVariable, const TDataType& Value )
        {
	        if (mpModel->GetProperties(PropertiesId).get() == 0)
            {
	            mpModel->AddProperties(PropertiesId, Properties::Pointer(new Properties(*mpModel)));
            }
    
	        typename PropertyFunction<TDataType>::Pointer constant_property(new ConstantProperty<TDataType>(Value));
    
	        mpModel->GetProperties(PropertiesId)->SetProperty(rVariable, constant_property);
        }
    
        //*****************************************************************************
        //* - Stores solution on variable (complex values)
        //*****************************************************************************
        template<class TDataType>
        void SaveSystemSolution(const Variable<TDataType>& rVariable)
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
	        for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
	  	        if ( ((*itDof)->GetVariableKey() == varKey) && (!(*itDof)->IsFixed()) )
                {
	  	            mpModel->Value(rVariable, *itDof) = x_vector[(*itDof)->EquationId()];
                }
            }
        }
    
        //*****************************************************************************
        //* - Stores solution on variable (Static problem mode)
        //*****************************************************************************
        template<class TDataType>
        void SaveSystemSolution_Static(const Variable<TDataType>& rVariable)
        {
            std::set<Dof::Pointer,ComparePDof>::iterator itDof;
    
            int varKey = rVariable.getKey();
    
            for( itDof = mDofSet.begin(); itDof != mDofSet.end(); ++itDof )
            {
                if ( ( (*itDof)->GetVariableKey() == varKey ) && ( !(*itDof)->IsFixed() ) )
                {
                    mpModel->Value(rVariable, *itDof) = std::real( x_vector[(*itDof)->EquationId()] );
                }
            }
        }
    
        void UpdateFixedVOLTAGE();
    
       //*****************************************************************************
       //* - Calculate derived nodal magnitudes
       //*****************************************************************************
        template<class TDataType>
        void CalculateNodal(const Variable<TDataType>& rVariable)
        {
            CalculateNodalE3D(rVariable.getKey());
        }
    
        void CalculateNodalE3D   ( int variableKey );
    
        void GeneratePECElement  ( std::vector<int>& LONodesId );
        void GeneratePSBCElement ( std::vector<int>& LONodesId );
	    void GeneratePMCElement  ( std::vector<int>& LONodesId );
        void GenerateTEPMCElement( std::vector<int>& LONodesId );
    
        Vector<double> Calculate_Geom_Averaged_Normal( std::vector< std::vector<double> >& vectorSet );
        Vector<double> Calculate_Area_Weighted_Normal( std::vector< std::vector<double> >& vectorSet );
    
        void CalculateDistToSing             ( int* NodesId, std::vector<Node::Pointer>& nodes );
        void GetGrad_SurfInt_NonSmooth_Static( int* NodesId                                    );
    
        void GenerateVEQElement           ( std::vector<int>& QHONodesId );
    
        void GenerateVolumeElement_Static ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolumeElement_Plasma ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolumeElement        ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateSourceElement        ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateVolIntElement        ( int* NodesId, unsigned int PropertiesId );
        void GenerateVolIntElement_Static ( int* NodesId, unsigned int PropertiesId );
    
        void GenerateSurfIntElement       ( int* NodesId, unsigned int PropertiesId );
        void GenerateSurfIntElement_Static( int* NodesId, unsigned int PropertiesId );
    
        void GeneratePBCElement           ( int* NodesId, unsigned int PropertiesId );

        void Set_Elements_Vector          ( int* NodesId, unsigned int PropertiesId );
        void Set_Nodal_Connectivities     (                                         );
        void Clear_Nodal_Connectivities   (                                         );
        void Set_Bubble_Derivatives       (                                         );
        void Clear_Bubble_Derivatives     (                                         );
    
        void SetPBC           ();
	    void Set_FrontBack_PBC();
	    void Set_RightLeft_PBC();
    
        void PBC_Mapping_Parameters       ( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY );
        void PBC_Mapping_Parameters_Linear( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY );
        void PBC_Mapping_Parameters_Quadrt( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY );
    
        int    PBC_Node_To_Node( double natX, double natY, std::vector<int>& ElementNodeIds );
        double PBC_Node_To_Edge( double natX, double natY, std::vector<int>& ElementNodeIds, std::vector<int>& EdgeNodeIds );
    
	    double               Determinant( Matrix<       double         >& M, int n );
        std::complex<double> Determinant( Matrix< std::complex<double> >& M, int n );
    
  	    double PBC_RotateY( Node::Pointer prNode );
    
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
    
	    void PBC_InvRot( int NodeId, double& cos_A, double& sin_A );
    
	    void UpdateCoordEPBC   ();
        void UpdateCoord_FB_PBC();
    
	    void UpdateCoord_RL_PBC(std::vector<int>& FrontLeftNodes);
	    void UpdateCoord_LF_PBC(std::vector<int>& FrontLeftNodes);
    
	    void GenerateRWPortTE10Element ( int* NodesId, unsigned int PropertiesId );
	    void GenerateCoaxPortTEMElement( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void GenerateGenericRobinElement       ( int* NodesId, unsigned int PropertiesId );
        void GenerateGenericRobinElement_Static( int* NodesId, unsigned int PropertiesId );
        void GenerateGenericRobinElement_Plasma( int* NodesId, unsigned int PropertiesId );
                                                                                         
	    void GenerateProjectionRWPortTE10 ( int* NodesId, unsigned int PropertiesId );
	    void GenerateProjectionCoaxPortTEM( int* NodesId, unsigned int PropertiesId );
                                                                                         
        void GenerateFarFieldElement       ( int* NodesId                            );
        void GenerateFarFieldElement_Plasma( int* NodesId, unsigned int PropertiesId );
    
        void Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix );
    
        void FixDof  (unsigned int NodeId, const Variable<double>& rVariable, const double& Value);
        void FixCDof (unsigned int NodeId, const Variable<double>& rVariable, const std::vector<double>& vValue);
	    void FreeCDof(unsigned int NodeId, const Variable<double>& rVariable);
    
        void SolveProblem      ();
        void SolveStaticProblem();
    
        void SetDofSet ();
        void SetNormals();
    
	    void SettingPECNormals  ();
        void SettingPSBCNormals ();
	    void SettingPMCNormals  ();
	    void SettingTEPMCNormals();
    
	    void Ortogonalization(Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto);
    
        void GenerateLinearSolver_Plasma(String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance);
        void GenerateLinearSolver       (String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance);
    
        void SetProblemType  (const String& ProblemType     );
        void SetFrequencyMode(bool          SweepingFreq    );
        void SetFrequency    (double        ProblemFrequency);
    
        void SolveLinearSystem();
    
        void IniBuilding();
    
	    void UpdateCoordEPEC();
    
        void TangencialCoordinates(std::vector<double>& n, std::vector<double>& t, std::vector<double>& b);
    
	    void FinishStep     ();
        void FinishFirstStep();
    
        void CalculateSijParameters();
    
	    void CalculateJouleHElement       (int* NodesId, unsigned int PropertiesId);
        void CalculateJouleHElement_Plasma(int* NodesId, unsigned int PropertiesId);
    
	    void CalculateHElement(int* NodesId, unsigned int PropertiesId);
    
        void Average_H_In_Contact_Nodes();
        void Apply_PEC_and_PMC_NBC_To_H();
    
        void IniStatic_E_Derivation     ();
        void IniStatic_J_Derivation     ();
        void IniStatic_JouleH_Derivation();
    
        void Calculate_Static_E_Element     (int* NodesId, unsigned int PropertiesId);
        void Calculate_Static_J_Element     (int* NodesId, unsigned int PropertiesId);
        void Calculate_Static_JouleH_Element(int* NodesId, unsigned int PropertiesId);
    
        void EndStatic_E_Derivation     ();
        void EndStatic_J_Derivation     ();
        void EndStatic_JouleH_Derivation();
    
        void Calculate_Static_E_Element_GP     (Vector< Vector<double> >& Static_E_GaussPs     , int* NodesId, unsigned int PropertiesId);
        void Calculate_Static_J_Element_GP     (Vector< Vector<double> >& Static_J_GaussPs     , int* NodesId, unsigned int PropertiesId);
        void Calculate_Static_JouleH_Element_GP(Vector<     double     >& Static_JouleH_GaussPs, int* NodesId, unsigned int PropertiesId);
    
        void Ini_Exporting_Static_Current();       
        void End_Exporting_Static_Current();
        void Export_Static_J_Element_GP  (int* NodesId, unsigned int PropertiesId);
    
        void Build_Imported_Static_Currents    ();
        void Calculate_Imported_Static_Currents();
        void Assemble_Imported_Current_Element (int* NodesId, double* Jvec, double Jpha);
    
	    void IniHDerivation ();
	    void IniJCalculation();
	    void EndJCalculation();
    
	    void FindSing2L(int* NodesId);
	    void FindSing3L(int* NodesId);
	    void FindSing4L(int* NodesId);
	    void FindSing5L(int* NodesId);
	    void FindSing6L(int* NodesId);
    
	    int ElementOrder();
    
	    void HO_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix);
	    void HO_TR (std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector  );
    
        void Calculate_En_R_invEp_Rt_Dn( Matrix< std::complex<double> >& R_invEp_Rt, int NodeId );
    
        void Plasma_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix);
        void Plasma_TR (std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector  );
    
        void Rotation_Matrix_SDP_T_CC(Matrix<std::complex<double> >& T, unsigned int NodeId);
    
	    void PushHONodesOnVolume (int* NodesId, std::vector<int>& HONodesId);
    
	    void PushHONodesOnSurface(int* NodesId             , std::vector<int>& HONodesId);
	    void PushHONodesOnSurface(std::vector<int>& NodesId, std::vector<int>& HONodesId);
    
	    void HONodesOnEdge  (int IdNode1, int IdNode2,              std::vector<int>& HONodesId);
	    void HONodesOnFace  (int IdNode1, int IdNode2, int IdNode3, std::vector<int>& HONodesId);
	    void HONodesOnVolume(int* NodesId,                          std::vector<int>& HONodesId);
    
	    void ReorderHONodesOnFace_4th(int orderingCase, std::vector<unsigned int>& face, std::vector<int>& HONodesId);
    
	    void CreateHONodes         (int* NodesId);
	    void CreateHONodes_OnEdges (int* NodesId);
	    void CreateHONodes_OnFaces (int* NodesId);
	    void CreateHONodes_OnVolume(int* NodesId);
    
	    void CreateHONodes_OnEdge_2nd(int IdNode1, int IdNode2);
	    void CreateHONodes_OnEdge_3th(int IdNode1, int IdNode2);
	    void CreateHONodes_OnEdge_4th(int IdNode1, int IdNode2);
    
	    void CreateHONodes_OnFace_3th(int IdNode1, int IdNode2, int IdNode3);
	    void CreateHONodes_OnFace_4th(int IdNode1, int IdNode2, int IdNode3);
	    
	    void SetBaseFileName(String baseFileName) { mBaseFileName = baseFileName; }
    
	    void PrintProjection       ();
        void PrintProjection_Static();
    
        struct ContactPairData
	    {
	  	    short int matR;
	  	    short int matL;
	  	    unsigned int ctcNode;
	  	    Vector<double> ctcNormal;
	    };
    
	    void ContactProperties(int* NodesId, unsigned int PropertiesId);
	    void DelTempContactVectors();
    
	    void GenerateContactPairs   (Vector<int>& CtCNodesId);
        void CalculateContactNormals(Vector<int>& NodesId   );
    
	    void SetContactDiscontinuity();
	    void SetAvgNormalsInContact ();
    
	    void JoinContactNormalsInDirichletS();
	    void ForceDirichletNormalsInContact();
    
	    void OrtogonalizeContactNormals();
    
	    void ContactHO_TKT( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix );
	    void ContactHO_TR ( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector   );
    
        void ContactHO_TKT_Plasma( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix );
        void ContactHO_TR_Plasma ( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector   );
        void Get_Contact_Matrix_Plasma( Matrix< std::complex<double> >& ContactMatrix, int RNodesId );
    
	    void UpdateContactCoord();
    
	    double VolumeOfTetrahedra(std::vector<Node::Pointer>& cnodes);
	    double AreaOfTriangle    (std::vector<Node::Pointer>& cnodes);
    
	    void SetAllcEzToCero();
	    void SetAllcExyToCero();
	    void SetAllcEyzToCero();
    
	    double DistanceToAxis         (unsigned int NodeId);
	    void   AxisymmetricCoordinates(unsigned int NodeId,std::vector<double>& n, std::vector<double>& t, std::vector<double>& b);
    
	    void Axisym_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix);
	    void Axisym_TR (std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector);
    
	    void UpdateAxisToCartesian();
    
	    void CalculateJVEElement       (int* NodesId, unsigned int PropertiesId);
        void CalculateJVEElement_Plasma(int* NodesId, unsigned int PropertiesId);
	    void CalculateJJEElement       (int* NodesId, unsigned int PropertiesId);
    
	    void Calculate_detJ(std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<Node::Pointer>& cnodes);
	    void Calculate_detJ(std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY, std::vector<Node::Pointer>& cnodes);
    
	    void NaturalDerivatives3D_2ndQ(Matrix<double>& dNda, Matrix<double>& dNdb, Matrix<double>& dNdu, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ);
    
	    void LoadPlasmaParameters(String plasma_file_name);
    
	    void IniEraseExtraNodes();
	    void EndEraseExtraNodes();
    
	    bool Is_OutputFileFormat_BIN();
	    bool Is_OutputFileFormat_ASCII();
    
        bool Is_StaticMode();
        bool Is_PlasmaMode();
    
	    void Extract_ColdPlasma_ElectronDensity();
	    void Extract_ColdPlasma_Bexternal();
	    void Extract_ColdPlasma_PermittivityTensor();
        void Extract_ColdPlasma_Eparallel();
        void Extract_ColdPlasma_Eperpendicular();
	    void Extract_ColdPlasma_AllParameters();
    
        void Set_E_parallel_zeros();
    
        void UpdateCoordPlasmaRLP();
    
        void CheckConsistency();
    
        void Send_Error_Msg( String ErrorMsg, String CheckMsg, int VerticalSpaces );
    
        void Show_Global_Matrix();
    
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
    
        Vector<int> mNeigbourElements;
    
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
    
        std::map<double, Vector<double> > mVolIntegralE_Static;
        std::map<double, Vector<double> > mVolIntegralE_Static_GP;
        std::map<double, Vector<double> > mVolIntegralJ_Static;
    
        std::map<double, Vector<double> > mSurfIntegralE_Static;
        std::map<double, Vector<double> > mSurfIntegralE_Static_GP;
        std::map<double, Vector<double> > mSurfIntegralJ_Static;
    
        int mNonSmooth_Surf_Id;
    
        Vector< Vector<double> > mSurfInt_NonSmooth_E_Static;
    
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

        std::map<unsigned int, Properties::Pointer> mPropertiesMap;

        Vector< Vector<unsigned int> > mNodesConnectivity;
        Vector< Vector<unsigned int> > mElementsVector;

        Vector<        double  > mbVolumes;
        Vector< Vector<double> > mb_dNk_dX;
        Vector< Vector<double> > mb_dNk_dY;
        Vector< Vector<double> > mb_dNk_dZ; 

        std::fstream mExportCurrentFile;
    
	    double mPBCTolerance;
        double mGiDTolerance;
    
	    bool mIsRightPBCTilted;
    
	    double mMaxDistToSing;
	    double mMinDistToSing;
    
	    int mMaxDistNode1, mMaxDistNode2;
	    int mMinDistNode1, mMinDistNode2;
    
	    bool mQuadraticGeometry;
    
	    bool mSetAllEzToZero;
	    bool mSetAllExyToZero;
	    bool mAxisymmetric;
    
	    bool mPlasmaMode;
        bool mStaticMode;
    
        bool mFrequencySweep;
    
	    bool mReadFileInitialGuess;
    
	    bool mReleaseSolutionMode;
	    bool mDebugSolutionMode;
        bool mReadSolutionMode;
    
	    bool mWriteResultEveryStep;
	    bool mWriteResultFinalStep;
    
	    bool mOutputFileFormatIsBIN;
	    bool mOutputFileFormatIsASCII;
    
        bool mImportCurrents;
    
	    std::vector<bool> mUsefulNodes;
    
	    unsigned int mItSolverNumThreads;
    };
}

#endif

