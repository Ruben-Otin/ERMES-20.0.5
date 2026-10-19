
#if !defined(KRATOS_MODELER )
#define  KRATOS_MODELER

#include <time.h>
#include <math.h>
#include <complex>
#include <map>

#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

#include "../includes/model.h"
#include "../includes/kratos_string.h"
#include "../includes/exception.h"
#include "../includes/point.h"
#include "../includes/process.h"
#include "../linear_solvers/complex_solver.h"

namespace Kratos
{

  class Modeler
  {
    public:

      typedef enum ProblemTypeType {E3D};
     
      typedef boost::shared_ptr<Modeler> Pointer;
      
      typedef Vector<Node::Pointer> NodesArrayType;
      
      typedef unsigned int IndexType;


      Modeler::Modeler(Model::Pointer pModel) : mpModel(pModel)
      {
          czero = std::complex<double>(0.00,0.00);

		  mMaxDistToSing = 0.00;
	      mMinDistToSing = 0.00;

	      mMaxDistNode1 = 0; mMaxDistNode2 = 0;
	      mMinDistNode1 = 0; mMinDistNode2 = 0;
      }

      virtual ~Modeler()
      {}
     
      //*****************************************************************************
      //* - Generate nodes
      //*****************************************************************************
      void GenerateNode(IndexType NodeId, double X, double Y, double Z)
      {
          Node::Pointer new_node(new Node(NodeId,X,Y,Z));
          mpModel->AddNode(new_node);
          new_node->AddDofs(mDofSet);
      }

      //*****************************************************************************
      //* - Set properties
      //*****************************************************************************
      template<class TDataType> 
	  void SetProperties(IndexType PropertiesId, 
			             const Variable<TDataType>& rVariable, 
			             const TDataType& Value)
      {
	      if (mpModel->GetProperties(PropertiesId).get() == 0)
          {
	          mpModel->AddProperties(PropertiesId, Properties::Pointer(new Properties(*mpModel)));
          }
        
	      PropertyFunction<TDataType>::Pointer constant_property(new ConstantProperty<TDataType>(Value));
        
	      mpModel->GetProperties(PropertiesId)->SetProperty(rVariable, constant_property);
      }

      //*****************************************************************************
      //* - Like former update
      //****************************************************************************
      template<class TDataType> 
      void SaveSystemSolution(const Variable<TDataType>& rVariable)
      {
          std::set<Dof::Pointer,ComparePDof>::iterator itDof;

          int variableKey = rVariable.getKey();
	
	      for (itDof=mDofSet.begin();itDof != mDofSet.end(); ++itDof)
          {
		      if ( ((*itDof)->GetVariableKey() == variableKey) && 
                   (!(*itDof)->IsFixed()) )
              {
		          mpModel->Value(rVariable, *itDof) = x_vector[(*itDof)->EquationId()]; 
              }
          } 
      }

     //*****************************************************************************
     //* - Calculate derived nodal magnitudes 
     //*****************************************************************************
      template<class TDataType>
      void CalculateNodal(const Variable<TDataType>& rVariable)
      {
          CalculateNodalE3D(rVariable.getKey());
      }

      void CalculateNodalE3D(int variableKey);

      void GeneratePECElement(std::vector<int>& LONodesId);
	  void GeneratePMCElement(std::vector<int>& LONodesId);

	  void GenerateQPECElement(std::vector<int>& LONodesId);

	  void GenerateVEQElement(std::vector<int>& QHONodesId);

	  void GenerateSPECElement (std::vector<int>& LONodesId);
	  void GenerateTEPMCElement(std::vector<int>& LONodesId);	  
	  
	  void CalculateDistToSing(int* NodesId, std::vector<Node::Pointer> nodes);

      void GenerateVolumeElement(int* NodesId, unsigned int PropertiesId);
      void GenerateSourceElement(int* NodesId, unsigned int PropertiesId);

	  void GenerateVolIntElement (int* NodesId, unsigned int PropertiesId);
	  void GenerateSurfIntElement(int* NodesId, unsigned int PropertiesId);

	  void GeneratePBCElement (int* NodesId, unsigned int PropertiesId);
	  void SetPBC             ();

	  void Set_FrontBack_PBC();
	  void Set_RightLeft_PBC();

	  bool FrontNodeIsBackNode( Node::Pointer pfNode, unsigned int backNodeId );
	  bool RightNodeIsLeftNode( Node::Pointer prNode, unsigned int leftNodeId );
	  
	  bool NearPointIsInExtEdge( unsigned int& centralNodeId, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );

	  void nCoordToXY( double nx, double ny, double& X, double& Y, std::vector<Node::Pointer>& pNodes );
	  void nCoordToZY( double nx, double ny, double& Z, double& Y, std::vector<Node::Pointer>& pNodes );

	  void FillElemSetNearNode( unsigned int NodeId, std::vector<std::vector<unsigned int> >& elemSet, std::vector<std::vector<unsigned int> >& elemSetNearNode );

	  void PBC_FB_NearestPointInElementSet( Node::Pointer pfNode, std::vector<std::vector<unsigned int> >& elemSet, std::vector<std::vector<double> >& minNatCoordSet, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );
	  void PBC_RL_NearestPointInElementSet( Node::Pointer prNode, std::vector<std::vector<unsigned int> >& elemSet, std::vector<std::vector<double> >& minNatCoordSet, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );

	  void PBC_FB_SetPointInElement( Node::Pointer pfNode, std::vector<unsigned int>& bNodesIds, std::vector<double>& natCoord);
	  void PBC_RL_SetPointInElement( Node::Pointer prNode, std::vector<unsigned int>& lNodesIds, std::vector<double>& natCoord);

	  void PBC_FB_SetPointInEdge( Node::Pointer pfNode, std::vector<unsigned int>& bNodesIds, double nx, double ny);
	  void PBC_RL_SetPointInEdge( Node::Pointer prNode, std::vector<unsigned int>& lNodesIds, double nx, double ny);

	  void PBC_FB_SearchMinPointInElementToTol( Node::Pointer pfNode, double prec, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );
	  void PBC_RL_SearchMinPointInElementToTol( Node::Pointer prNode, double prec, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );

	  void PBC_FB_SearchMinPointInEdges( Node::Pointer pfNode, double prec, std::vector<std::vector<unsigned int> >& elemSet, std::vector<std::vector<double> >& minNatCoordSet, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );
	  void PBC_RL_SearchMinPointInEdges( Node::Pointer prNode, double prec, std::vector<std::vector<unsigned int> >& elemSet, std::vector<std::vector<double> >& minNatCoordSet, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord );

	  unsigned int PBC_FB_SearchNearestNode( Node::Pointer pfNode, std::vector<std::vector<unsigned int> >& elemSetBack );
	  unsigned int PBC_RL_SearchNearestNode( Node::Pointer prNode, std::vector<std::vector<unsigned int> >& elemSetLeft );

      double PBC_FB_MinDistPointInEle( Node::Pointer pfNode, std::vector<unsigned int>& NodesIds, double nx_ini, double nx_end, double ny_ini, double ny_end, double step, std::vector<double>& natCoord );
	  double PBC_RL_MinDistPointInEle( Node::Pointer prNode, std::vector<unsigned int>& NodesIds, double nx_ini, double nx_end, double ny_ini, double ny_end, double step, std::vector<double>& natCoord );

 	  double PBC_RotateY(Node::Pointer prNode);

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
	  void PBC_ResVectorReDef      (int local_i, int oldSize, int newSize, std::vector<double>& N, Vector<std::complex<double> >& ResVector);

	  void PBC_StiffMatrixReDef_Rot(int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Matrix<std::complex<double> >& StiffMatrix);
	  void PBC_ResVectorReDef_Rot  (int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Vector<std::complex<double> >& ResVector);

	  void PBC_InvRot(unsigned int NodeId, double& cos_A, double& sin_A );	

	  bool PBC_LeftNodeIsFrontNode(unsigned int rightNodeId, std::vector<unsigned int>& leftNodesId, std::vector<unsigned int>& FrontLeftNodesId);

	  void UpdateCoordEPBC   ();
	  void UpdateCoord_RL_PBC(std::vector<unsigned int>& FrontLeftNodes);
	  void UpdateCoord_FB_PBC();
	  void UpdateCoord_LF_PBC(std::vector<unsigned int>& FrontLeftNodes);

	  void GenerateRWPortTE10Element  (int* NodesId, unsigned int PropertiesId);
	  void GenerateCoaxPortTEMElement (int* NodesId, unsigned int PropertiesId);
	  void GenerateGenericRobinElement(int* NodesId, unsigned int PropertiesId);

      void GenerateFarFieldElement(int* NodesId);

	  void GenerateProjectionRWPortTE10 (int*  NodesId, unsigned int PropertiesId);
	  void GenerateProjectionCoaxPortTEM(int*  NodesId, unsigned int PropertiesId);

      void GenerateElement (const String& ElementName, std::vector<int>& NodesId, unsigned int PropertiesId)
      {};

      void GenerateSElement(const String& ElementName, std::vector<int>& NodesId, unsigned int PropertiesId)
      {};

      void FixDof  (unsigned int NodeId, const Variable<double>& rVariable, const double& Value);
      void FixCDof (unsigned int NodeId, const Variable<double>& rVariable, const std::vector<double>& vValue);
	  void FreeCDof(unsigned int NodeId, const Variable<double>& rVariable);
 
      void SolveProblem();

      void SetDofSet();

      void SetNormals();

	  void SettingPECNormals();
	  void SettingPMCNormals();

	  void SettingSPECNormals();
	  void SettingTEPMCNormals();

	  void Ortogonalization(Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto);
     
      void GenerateLinearSolver(String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance);
     
      void SetProblemType(const String& ProblemType);

      void SetFrequency(double ProblemFrequency);

      void SolveLinearSystem();

      void IniBuilding();

	  void UpdateCoordEPEC();

      void TangencialCoordinates(std::vector<double>& n, std::vector<double>& t, std::vector<double>& b);

	  void FinishStep();
      void FinishFirstStep();

      void CalculateSijParameters();

	  void CalculateJouleHElement(int* NodesId, unsigned int PropertiesId);
	  void CalculateHElement     (int* NodesId, unsigned int PropertiesId);

	  void IniHDerivation();
	  void IniJCalculation();
	  void EndJCalculation();

	  void FindSing2L(int* NodesId);
	  void FindSing3L(int* NodesId);
	  void FindSing4L(int* NodesId);
	  void FindSing5L(int* NodesId);
	  void FindSing6L(int* NodesId);

	  int ElementOrder();

	  void HO_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix);
	  void HO_TR (std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector);

	  void PushHONodesOnVolume (int* NodesId, std::vector<int>& HONodesId);

	  void PushHONodesOnSurface(int* NodesId             , std::vector<int>& HONodesId);
	  void PushHONodesOnSurface(std::vector<int>& NodesId, std::vector<int>& HONodesId);

	  void HONodesOnEdge (int IdNode1, int IdNode2,              std::vector<int>& HONodesId);
	  void HONodesOnFace (int IdNode1, int IdNode2, int IdNode3, std::vector<int>& HONodesId);
	  void HONodesOnVolum(int* NodesId,                          std::vector<int>& HONodesId);

	  void ReorderHONodesOnFace_4th(int orderingCase, std::vector<unsigned int>& face, std::vector<int>& HONodesId);

	  void CreateHONodes        (int* NodesId);
	  void CreateHONodes_OnEdges(int* NodesId);
	  void CreateHONodes_OnFaces(int* NodesId);
	  void CreateHONodes_OnVolum(int* NodesId);

	  void CreateHONodes_OnEdge_2nd(int IdNode1, int IdNode2);

	  void CreateHONodes_OnEdge_3th(int IdNode1, int IdNode2);
	  void CreateHONodes_OnFace_3th(int IdNode1, int IdNode2, int IdNode3);

	  void CreateHONodes_OnEdge_4th (int IdNode1, int IdNode2);
	  void CreateHONodes_OnFace_4th (int IdNode1, int IdNode2, int IdNode3);
	  void CreateHONodes_OnVolum_4th(int* IdNodes);

	  void SetBaseFileName(String baseFileName) { mBaseFileName = baseFileName; }

	  void PrintProjection();

      struct ContactPairData
	  {  
		  short int matR;
		  short int matL;
		  unsigned int ctcNode;
		  Vector<double> ctcNormal;
	  };

	  void ContactProperties(int* NodesId, unsigned int PropertiesId);
	  void DelTempContactVectors();

	  void GenerateContactPairs(Vector<int>& CtCNodesId);
	  void SetContactDiscontinuity();

	  void CalculateContactNormals(Vector<int>& NodesId);
	  void GenerateCorrectCTCNormals(Vector<int>& LONodesId);
	  void SetAvgNormalsInContact();

	  void CalculateContactQNormals  (std::vector<int>& LONodesId);
	  void GenerateCorrectCTCQNormals(std::vector<int>& LONodesId);

	  void AddContactPECNormals();
	  void ForcePECNormalsInContact();

	  void OrtogonalizeContactNormals();

	  void ContactHO_TKT(std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix);
	  void ContactHO_TR (std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector);

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

	  void CalculateJVEElement(int* NodesId, unsigned int PropertiesId);
	  void CalculateJJEElement(int* NodesId, unsigned int PropertiesId);

	  void SurfDerivativesInNodes_2ndQ(Matrix<double>& dNda, Matrix<double>& dNdb);

	  void Calculate_detJ(std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<Node::Pointer>& cnodes);
	  void Calculate_detJ(std::vector<double>& detJ, std::vector<double>& cX, std::vector<double>& cY, std::vector<Node::Pointer>& cnodes);

	  void NaturalDerivatives3D_2ndQ(Matrix<double>& dNda, Matrix<double>& dNdb, Matrix<double>& dNdu, std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ);

	  void IniEraseExtraNodes();
	  void EndEraseExtraNodes();

	  bool Is_OutputFileFormat_BIN();
	  bool Is_OutputFileFormat_ASCII();

    private:

      std::complex<double> czero; 
          
      Model::Pointer mpModel;
      
      CSRMatrix<std::complex<double> > A_matrix;
	  Vector<std::complex<double> > x_vector;
	  Vector<std::complex<double> > b_vector;

      time_t mStart;
      time_t mEnd;

      ProblemTypeType mProblemType; 

      double mProblemFrequency;

      std::set<Dof::Pointer, ComparePDof> mDofSet;

      std::map<unsigned int, std::vector<double> > mNormals;

      std::map<unsigned int, std::vector<std::vector<double> > > mSetPECnormals;
	  std::map<unsigned int, std::vector<std::vector<double> > > mSetPMCnormals;

	  std::map<unsigned int, std::vector<std::vector<double> > > mSetSPECnormals;
	  std::map<unsigned int, std::vector<std::vector<double> > > mSetTEPMCnormals;

      ComplexSolver::Pointer mpLinearSolver;
	  
	  String mExternalSolverPath;
	  String mBaseFileName;

	  bool mSolveWithExternal;
	  bool mSolveWithBiCG;
	  bool mSolveWithQMR;
	  bool mSolveWithCG;

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

	  std::map<double, double> mVolumeValue;
      std::map<double, double> mSurfaceValue;  

	  int mElementOrder;

	  std::map<unsigned int, double> mSingular;
      std::map<unsigned int, double> mNSingular2L;
      std::map<unsigned int, double> mNSingular3L;
      std::map<unsigned int, double> mNSingular4L;
	  std::map<unsigned int, double> mNSingular5L;
      std::map<unsigned int, double> mNSingular6L;

	  std::map<unsigned int, std::vector<std::vector<unsigned int> > > mEdgeHONodes;
      std::map<unsigned int, std::vector<std::vector<unsigned int> > > mFaceHONodes;
	  std::map<unsigned int, std::vector<std::vector<unsigned int> > > mVolmHONodes;

	  std::map<unsigned int, unsigned short int> mTempNodeProperties;
	  std::map<unsigned int, std::vector<std::vector<double> > > mSetContactNormals;
	  std::map<unsigned int, std::vector<std::vector<double> > > mSetCorrectedCtCNormals;
	  std::map<unsigned int, std::vector<double> > mContactNormals;
	  std::map<unsigned int, unsigned int> mTempContactPairs;
	   
	  std::map<unsigned int, ContactPairData> mContactPairs;

	  std::map<unsigned int, std::vector<double> > tmp_rJ;
	  std::map<unsigned int, std::vector<double> > tmp_iJ;

	  std::map<unsigned int, std::vector<std::vector<unsigned int> > > mPBCFrontElements;
	  std::map<unsigned int, std::vector<std::vector<unsigned int> > > mPBCBackElements;

	  std::map<unsigned int, std::vector<std::vector<unsigned int> > > mPBCRightElements;
      std::map<unsigned int, std::vector<std::vector<unsigned int> > > mPBCLeftElements;

	  std::map<unsigned int, unsigned int> mPBC_NodeNodePairs_Front;
	  std::map<unsigned int, unsigned int> mPBC_NodeNodePairs_Right;

      std::map<unsigned int, std::vector<unsigned int> > mPBC_NodeEdgePairs_Front;	  
	  std::map<unsigned int, std::vector<unsigned int> > mPBC_NodeEdgePairs_Right;

      std::map<unsigned int, std::vector<unsigned int> > mPBC_NodeElementPairs_Front;
      std::map<unsigned int, std::vector<unsigned int> > mPBC_NodeElementPairs_Right;

	  std::map<unsigned int, double> mPBC_NodeEdgeCoord_Front;
	  std::map<unsigned int, double> mPBC_NodeEdgeCoord_Right;

	  std::map<unsigned int, std::vector<double> > mPBC_NodeElementCoord_Front;
	  std::map<unsigned int, std::vector<double> > mPBC_NodeElementCoord_Right;

	  double mPBCTolerance;

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

	  bool mReadFileInitialGuess;

	  bool mReleaseSolutionMode;
	  bool mDebugSolutionMode;
      bool mReadSolutionMode;

	  bool mWriteResultEveryStep;
	  bool mWriteResultFinalStep;

	  bool mOutputFileFormatIsBIN;
	  bool mOutputFileFormatIsASCII;

	  std::vector<bool> mUsefulNodes;

	  unsigned int mItSolverNumThreads;

    }; 
  
}  

#endif 

