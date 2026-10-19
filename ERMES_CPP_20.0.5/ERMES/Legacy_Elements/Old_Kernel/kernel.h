
#if !defined(KRATOS_KERNEL )
#define  KRATOS_KERNEL

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
  
      Kernel();
     
      virtual ~Kernel();
          
      void GenerateNode(IndexType NodeId, double X, double Y, double Z)
      {
          mpModeler->GenerateNode(NodeId, X, Y, Z);   
      }
      
      template<class TDataType> 
      void SetProperties(IndexType PropertiesId, 
                         const Variable<TDataType>& rVariable, 
                         const TDataType& Value)
      {
          mpModeler->SetProperties(PropertiesId, rVariable, Value);
      }
      
      ////////////////////////////////////////////////////////////////////////////
      //
      ///////////////////////////////////////////////////////////////////////////
      void GenerateElement(const String& ElementsGroup,const String& ElementName, 
                           IndexType ElementId, Vector<int>& NodesId, 
                           IndexType PropertiesId, const String& ElementData)
      {
          if (mBuild) 
              mpModeler->GenerateElement (ElementName, NodesId, PropertiesId);
          else 
              mpModeler->GenerateSElement(ElementName, NodesId, PropertiesId);
      }

      void GenerateSElement(const String& ElementName, Vector<int>& NodesId,IndexType PropertiesId)
      {
          if (mBuild) 
              mpModeler->GenerateElement(ElementName, NodesId, PropertiesId);
          else 
              mpModeler->GenerateSElement(ElementName, NodesId, PropertiesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      //
      ///////////////////////////////////////////////////////////////////////////
      void GenerateContactPairs(Vector<int>& CtCNodesId)
      {
           mpModeler->GenerateContactPairs(CtCNodesId);
      }
	  void GenerateCorrectCTCNormals(Vector<int>& NodesId)
      {
           mpModeler->GenerateCorrectCTCNormals(NodesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      //
      ///////////////////////////////////////////////////////////////////////////
      void GeneratePECElement(Vector<int>& NodesId)
      {
           mpModeler->GeneratePECElement(NodesId);
      }
	  void GeneratePMCElement(Vector<int>& NodesId)
      {
           mpModeler->GeneratePMCElement(NodesId);
      }
	  void GenerateVEQElement(Vector<int>& NodesId)
      {
           mpModeler->GenerateVEQElement(NodesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      //
      ///////////////////////////////////////////////////////////////////////////
      void GenerateSPECElement(Vector<int>& NodesId)
      {
           mpModeler->GenerateSPECElement(NodesId);
      }
	  void GenerateTEPMCElement(Vector<int>& NodesId)
      {
           mpModeler->GenerateTEPMCElement(NodesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      //
      ///////////////////////////////////////////////////////////////////////////
      void GenerateVolumeElement(int* NodesId, unsigned int PropertiesId)
      {   
           if (mBuild)
		   {
			    mpModeler->GenerateVolumeElement(NodesId, PropertiesId);
		   }
		   else if (mCreateHONodes)
		   {
			    mpModeler->CreateHONodes(NodesId);
				mpModeler->ContactProperties(NodesId, PropertiesId);
		   }
		   else if (mPrintHOMesh)
		   {
			    int Order = mpModeler->ElementOrder();

                if      ( Order == 1 ) PrintHOElementInMesh_1stToLinear(NodesId, PropertiesId);
				else if ( Order == 2 ) PrintHOElementInMesh_2ndToLinear(NodesId, PropertiesId);
				else if ( Order == 3 ) PrintHOElementInMesh_3thToLinear(NodesId, PropertiesId);
				else if ( Order == 4 ) PrintHOElementInMesh_4thToLinear(NodesId, PropertiesId);
		   }
		   else if (mFSing2L)
		   {
			    mpModeler->FindSing2L(NodesId);
		   }
		   else if (mFSing3L)
		   {
			    mpModeler->FindSing3L(NodesId);
		   }
		   else if (mFSing4L)
		   {
			    mpModeler->FindSing4L(NodesId);
		   }
		   else if (mFSing5L)
		   {
			    mpModeler->FindSing5L(NodesId);
		   }
		   else if (mFSing6L)
		   {
			    mpModeler->FindSing6L(NodesId);
		   }
		   else if (mCalculateH)
		   {
			    mpModeler->CalculateHElement(NodesId, PropertiesId);
		   }	
		   else if (mCalculateJouleH)
		   {
			    mpModeler->CalculateJouleHElement(NodesId, PropertiesId);
		   }
		   else if (mCalculateJ)
		   {
			    mpModeler->CalculateJVEElement(NodesId, PropertiesId);
		   }	
      }

	  ////////////////////////////////////////////////////////////////////////////
      // - Special volumetric elements
      ///////////////////////////////////////////////////////////////////////////
      void GenerateSourceElement(int* NodesId, unsigned int PropertiesId)
      {   
         if (!mCalculateJ) mpModeler->GenerateSourceElement(NodesId, PropertiesId);
		 else              mpModeler->CalculateJJEElement  (NodesId, PropertiesId);
      }

	  void GenerateVolIntElement(int* NodesId, unsigned int PropertiesId)
      {   
          mpModeler->GenerateVolIntElement(NodesId, PropertiesId);
      }

	  void GeneratePBCElement(int* NodesId, unsigned int PropertiesId)
      {   
          mpModeler->GeneratePBCElement(NodesId, PropertiesId);
      }

	  void GenerateSurfIntElement(int* NodesId, unsigned int PropertiesId)
      {   
          mpModeler->GenerateSurfIntElement(NodesId, PropertiesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      // - Generate boundary elements
      ///////////////////////////////////////////////////////////////////////////
	  void GenerateBoundaryElement(int* NodesId, unsigned int PropertiesId)
      {}

	  void GenerateFarFieldElement(int* NodesId)
      {
          if (!mCalculateJ) mpModeler->GenerateFarFieldElement(NodesId);
      }

	  void GenerateGenericRobinElement(int* NodesId, unsigned int PropertiesId)
      {
		  if (!mCalculateJ) mpModeler->GenerateGenericRobinElement(NodesId, PropertiesId);
	  }

	  void GenerateCoaxPortTEMElement(int* NodesId, unsigned int PropertiesId)
      {
          if (!mCalculateJ) mpModeler->GenerateCoaxPortTEMElement(NodesId, PropertiesId);
      }

	  void GenerateRWPortTE10Element(int* NodesId, unsigned int PropertiesId)
      {
          if (!mCalculateJ) mpModeler->GenerateRWPortTE10Element(NodesId, PropertiesId);
      }

	  ////////////////////////////////////////////////////////////////////////////
      // - Projection in port planes
      ///////////////////////////////////////////////////////////////////////////
      void GenerateProjectionElement(int* NodesId, unsigned int PropertiesId)
      {}

	  void GenerateProjectionCoaxialTEM(int* NodesId, unsigned int PropertiesId)
      {
		  mpModeler->GenerateProjectionCoaxPortTEM(NodesId, PropertiesId);
      }

	  void GenerateProjectionRWPortTE10(int* NodesId, unsigned int PropertiesId)
      {
		  mpModeler->GenerateProjectionRWPortTE10(NodesId, PropertiesId);
      }
      
	  void GenerateProjectionLineElement(int* NodesId, unsigned int PropertiesId)
      {}

	  ////////////////////////////////////////////////////////////////////////////
      // -  Printing linear elements from 1st order elements
      ///////////////////////////////////////////////////////////////////////////
	  void PrintHOElementInMesh_1stToLinear(int* NodesId, unsigned int PropertiesId)
	  {
		  int NodesIdMat[5];

		  NodesIdMat[0] = NodesId[0];
          NodesIdMat[1] = NodesId[1];
          NodesIdMat[2] = NodesId[2];
		  NodesIdMat[3] = NodesId[3];
		  NodesIdMat[4] = PropertiesId;

		  //[0 1 2 3]
          mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, NodesIdMat);
	  }

	  ////////////////////////////////////////////////////////////////////////////
      // -  Printing linear elements from 2nd order elements
      ///////////////////////////////////////////////////////////////////////////
	  void Print_2ndToLinear(int* NodesId, unsigned int PropertiesId)
	  {
		  int SubElLinear[5];

		  //Material of the element(for working with layers in post-process)
		  SubElLinear[4] = PropertiesId;

		  //[0 4 6 7]
          SubElLinear[0] = NodesId[0];
          SubElLinear[1] = NodesId[4];
          SubElLinear[2] = NodesId[6];
		  SubElLinear[3] = NodesId[7];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[4 1 5 8]
          SubElLinear[0] = NodesId[4];
          SubElLinear[1] = NodesId[1];
          SubElLinear[2] = NodesId[5];
		  SubElLinear[3] = NodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[6 5 2 9]
          SubElLinear[0] = NodesId[6];
          SubElLinear[1] = NodesId[5];
          SubElLinear[2] = NodesId[2];
		  SubElLinear[3] = NodesId[9];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[7 8 9 3]
          SubElLinear[0] = NodesId[7];
          SubElLinear[1] = NodesId[8];
          SubElLinear[2] = NodesId[9];
		  SubElLinear[3] = NodesId[3];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		   //[6 4 5 8]
          SubElLinear[0] = NodesId[6];
          SubElLinear[1] = NodesId[4];
          SubElLinear[2] = NodesId[5];
		  SubElLinear[3] = NodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[6 8 5 9]
          SubElLinear[0] = NodesId[6];
          SubElLinear[1] = NodesId[8];
          SubElLinear[2] = NodesId[5];
		  SubElLinear[3] = NodesId[9];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[7 4 6 8]
          SubElLinear[0] = NodesId[7];
          SubElLinear[1] = NodesId[4];
          SubElLinear[2] = NodesId[6];
		  SubElLinear[3] = NodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[7 8 6 9]
          SubElLinear[0] = NodesId[7];
          SubElLinear[1] = NodesId[8];
          SubElLinear[2] = NodesId[6];
		  SubElLinear[3] = NodesId[9];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);
	  }

	  ////////////////////////////////////////////////////////////////////////////
      // -  Printing linear elements from 2nd order elements
      ///////////////////////////////////////////////////////////////////////////
	  void PrintHOElementInMesh_2ndToLinear(int* NodesId, unsigned int PropertiesId)
	  {
		  std::vector<int> HONodesId;

		  mpModeler->PushHONodesOnVolume(NodesId, HONodesId);

          int HOEl2nd[10];

		  //[0 1 2 3]
          HOEl2nd[0] = HONodesId[0];
          HOEl2nd[1] = HONodesId[1];
          HOEl2nd[2] = HONodesId[2];
		  HOEl2nd[3] = HONodesId[3];
          HOEl2nd[4] = HONodesId[4];
          HOEl2nd[5] = HONodesId[5];
		  HOEl2nd[6] = HONodesId[6];
          HOEl2nd[7] = HONodesId[7];
		  HOEl2nd[8] = HONodesId[8];
          HOEl2nd[9] = HONodesId[9];
		  Print_2ndToLinear(HOEl2nd, PropertiesId);
	  }

	  ////////////////////////////////////////////////////////////////////////////
      // -  Printing 3th order elements --> to linear
      ///////////////////////////////////////////////////////////////////////////
	  void PrintHOElementInMesh_3thToLinear(int* NodesId, unsigned int PropertiesId)
	  {
		  std::vector<int> HONodesId;

		  mpModeler->PushHONodesOnVolume(NodesId, HONodesId);

		  int SubElLinear[4];

		  //Material of the element(for working with layers in post-process)
		  SubElLinear[4] = PropertiesId;

		  //[0 4 9 10]
          SubElLinear[0] = HONodesId[0];
          SubElLinear[1] = HONodesId[4];
          SubElLinear[2] = HONodesId[9];
		  SubElLinear[3] = HONodesId[10];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

          //[4 10 17 9]
		  SubElLinear[0] = HONodesId[4];
          SubElLinear[1] = HONodesId[10];
          SubElLinear[2] = HONodesId[17];
		  SubElLinear[3] = HONodesId[9];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

          //[9 19 10 17]
		  SubElLinear[0] = HONodesId[9];
          SubElLinear[1] = HONodesId[19];
          SubElLinear[2] = HONodesId[10];
		  SubElLinear[3] = HONodesId[17];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[4 16 9 17]
		  SubElLinear[0] = HONodesId[4];
          SubElLinear[1] = HONodesId[16];
          SubElLinear[2] = HONodesId[9];
		  SubElLinear[3] = HONodesId[17];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[16 9 19 17]
		  SubElLinear[0] = HONodesId[16];
          SubElLinear[1] = HONodesId[9];
          SubElLinear[2] = HONodesId[19];
		  SubElLinear[3] = HONodesId[17];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[16 18 17 19]
		  SubElLinear[0] = HONodesId[16];
          SubElLinear[1] = HONodesId[18];
          SubElLinear[2] = HONodesId[17];
		  SubElLinear[3] = HONodesId[19];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[4 5 16 17]
		  SubElLinear[0] = HONodesId[4];
          SubElLinear[1] = HONodesId[5];
          SubElLinear[2] = HONodesId[16];
		  SubElLinear[3] = HONodesId[17];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[17 12 5 16]
		  SubElLinear[0] = HONodesId[17];
          SubElLinear[1] = HONodesId[12];
          SubElLinear[2] = HONodesId[5];
		  SubElLinear[3] = HONodesId[16];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[18 17 12 16]
		  SubElLinear[0] = HONodesId[18];
          SubElLinear[1] = HONodesId[17];
          SubElLinear[2] = HONodesId[12];
		  SubElLinear[3] = HONodesId[16];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[5 6 16 12]
		  SubElLinear[0] = HONodesId[5];
          SubElLinear[1] = HONodesId[6];
          SubElLinear[2] = HONodesId[16];
		  SubElLinear[3] = HONodesId[12];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[6 16 12 18]
		  SubElLinear[0] = HONodesId[6];
          SubElLinear[1] = HONodesId[16];
          SubElLinear[2] = HONodesId[12];
		  SubElLinear[3] = HONodesId[18];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[5 1 6 12]
          SubElLinear[0] = HONodesId[5];
          SubElLinear[1] = HONodesId[1];
          SubElLinear[2] = HONodesId[6];
		  SubElLinear[3] = HONodesId[12];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

          //[16 6 7 18]
          SubElLinear[0] = HONodesId[16];
          SubElLinear[1] = HONodesId[6];
          SubElLinear[2] = HONodesId[7];
		  SubElLinear[3] = HONodesId[18];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[8 7 2 14]
          SubElLinear[0] = HONodesId[8];
          SubElLinear[1] = HONodesId[7];
          SubElLinear[2] = HONodesId[2];
		  SubElLinear[3] = HONodesId[14];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[9 16 8 19]
		  SubElLinear[0] = HONodesId[9];
          SubElLinear[1] = HONodesId[16];
          SubElLinear[2] = HONodesId[8];
		  SubElLinear[3] = HONodesId[19];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[19 18 16 8]
		  SubElLinear[0] = HONodesId[19];
          SubElLinear[1] = HONodesId[18];
          SubElLinear[2] = HONodesId[16];
		  SubElLinear[3] = HONodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[19 8 14 18]
		  SubElLinear[0] = HONodesId[19];
          SubElLinear[1] = HONodesId[8];
          SubElLinear[2] = HONodesId[14];
		  SubElLinear[3] = HONodesId[18];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[7 16 18 8]
		  SubElLinear[0] = HONodesId[7];
          SubElLinear[1] = HONodesId[16];
          SubElLinear[2] = HONodesId[18];
		  SubElLinear[3] = HONodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

		  //[18 14 7 8]
		  SubElLinear[0] = HONodesId[18];
          SubElLinear[1] = HONodesId[14];
          SubElLinear[2] = HONodesId[7];
		  SubElLinear[3] = HONodesId[8];
	      mElementIdCounter++;
          mpDefaultOutput->PrintHOElement(mElementIdCounter, SubElLinear);

          int HOSubEl2nd[10];

		  //[10 12 14 3]
          HOSubEl2nd[0] = HONodesId[10];
          HOSubEl2nd[1] = HONodesId[12];
          HOSubEl2nd[2] = HONodesId[14];
		  HOSubEl2nd[3] = HONodesId[3];
          HOSubEl2nd[4] = HONodesId[17];
          HOSubEl2nd[5] = HONodesId[18];
		  HOSubEl2nd[6] = HONodesId[19];
          HOSubEl2nd[7] = HONodesId[11];
		  HOSubEl2nd[8] = HONodesId[13];
          HOSubEl2nd[9] = HONodesId[15];
		  Print_2ndToLinear(HOSubEl2nd, PropertiesId);
	  }

	  /////////////////////////////////////////////////////////////////////////////////////
      // -  Printing 4th order elements (8 2nd order elements, 64 1st order elements)
      /////////////////////////////////////////////////////////////////////////////////////
      void PrintHOElementInMesh_4thToLinear(int* NodesId, unsigned int PropertiesId)
      {
          std::vector<int> HONodesId;

		  mpModeler->PushHONodesOnVolume(NodesId, HONodesId);

          int HOSubEl2nd[10];

		  //[0 5 11 14]
          HOSubEl2nd[0] = HONodesId[0];
          HOSubEl2nd[1] = HONodesId[5];
          HOSubEl2nd[2] = HONodesId[11];
		  HOSubEl2nd[3] = HONodesId[14];
          HOSubEl2nd[4] = HONodesId[4];
          HOSubEl2nd[5] = HONodesId[22];
		  HOSubEl2nd[6] = HONodesId[12];
          HOSubEl2nd[7] = HONodesId[13];
		  HOSubEl2nd[8] = HONodesId[25];
          HOSubEl2nd[9] = HONodesId[31];
		  Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[5 1 8 17]
          HOSubEl2nd[0] = HONodesId[5];
          HOSubEl2nd[1] = HONodesId[1];
          HOSubEl2nd[2] = HONodesId[8];
		  HOSubEl2nd[3] = HONodesId[17];
          HOSubEl2nd[4] = HONodesId[6];
          HOSubEl2nd[5] = HONodesId[7];
		  HOSubEl2nd[6] = HONodesId[23];
          HOSubEl2nd[7] = HONodesId[26];
		  HOSubEl2nd[8] = HONodesId[16];
          HOSubEl2nd[9] = HONodesId[28];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[11 8 2 20]
          HOSubEl2nd[0] = HONodesId[11];
          HOSubEl2nd[1] = HONodesId[8];
          HOSubEl2nd[2] = HONodesId[2];
		  HOSubEl2nd[3] = HONodesId[20];
          HOSubEl2nd[4] = HONodesId[24];
          HOSubEl2nd[5] = HONodesId[9];
		  HOSubEl2nd[6] = HONodesId[10];
          HOSubEl2nd[7] = HONodesId[32];
		  HOSubEl2nd[8] = HONodesId[29];
          HOSubEl2nd[9] = HONodesId[19];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[14 17 20 3]
          HOSubEl2nd[0] = HONodesId[14];
          HOSubEl2nd[1] = HONodesId[17];
          HOSubEl2nd[2] = HONodesId[20];
		  HOSubEl2nd[3] = HONodesId[3];
          HOSubEl2nd[4] = HONodesId[27];
          HOSubEl2nd[5] = HONodesId[30];
		  HOSubEl2nd[6] = HONodesId[33];
          HOSubEl2nd[7] = HONodesId[15];
		  HOSubEl2nd[8] = HONodesId[18];
          HOSubEl2nd[9] = HONodesId[21];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[11 5 8 17]
          HOSubEl2nd[0] = HONodesId[11];
          HOSubEl2nd[1] = HONodesId[5];
          HOSubEl2nd[2] = HONodesId[8];
		  HOSubEl2nd[3] = HONodesId[17];
          HOSubEl2nd[4] = HONodesId[22];
          HOSubEl2nd[5] = HONodesId[23];
		  HOSubEl2nd[6] = HONodesId[24];
          HOSubEl2nd[7] = HONodesId[34];
		  HOSubEl2nd[8] = HONodesId[26];
          HOSubEl2nd[9] = HONodesId[28];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[11 17 8 20]
          HOSubEl2nd[0] = HONodesId[11];
          HOSubEl2nd[1] = HONodesId[17];
          HOSubEl2nd[2] = HONodesId[8];
		  HOSubEl2nd[3] = HONodesId[20];
          HOSubEl2nd[4] = HONodesId[34];
          HOSubEl2nd[5] = HONodesId[28];
		  HOSubEl2nd[6] = HONodesId[24];
          HOSubEl2nd[7] = HONodesId[32];
		  HOSubEl2nd[8] = HONodesId[30];
          HOSubEl2nd[9] = HONodesId[29];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[14 5 11 17]
          HOSubEl2nd[0] = HONodesId[14];
          HOSubEl2nd[1] = HONodesId[5];
          HOSubEl2nd[2] = HONodesId[11];
		  HOSubEl2nd[3] = HONodesId[17];
          HOSubEl2nd[4] = HONodesId[25];
          HOSubEl2nd[5] = HONodesId[22];
		  HOSubEl2nd[6] = HONodesId[31];
          HOSubEl2nd[7] = HONodesId[27];
		  HOSubEl2nd[8] = HONodesId[26];
          HOSubEl2nd[9] = HONodesId[34];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);

		  //[14 17 11 20]
          HOSubEl2nd[0] = HONodesId[14];
          HOSubEl2nd[1] = HONodesId[17];
          HOSubEl2nd[2] = HONodesId[11];
		  HOSubEl2nd[3] = HONodesId[20];
          HOSubEl2nd[4] = HONodesId[27];
          HOSubEl2nd[5] = HONodesId[34];
		  HOSubEl2nd[6] = HONodesId[31];
          HOSubEl2nd[7] = HONodesId[33];
		  HOSubEl2nd[8] = HONodesId[30];
          HOSubEl2nd[9] = HONodesId[32];
          Print_2ndToLinear(HOSubEl2nd, PropertiesId);
      }

      void FixDof(IndexType NodeId,const Variable<double>& rVariable,const double& Value)
      {
           mpModeler->FixDof(NodeId, rVariable, Value);    
      }

      void FixCDof(IndexType NodeId,const Variable<double>& rVariable,const Vector<double>& vValue)
      {
           mpModeler->FixCDof(NodeId, rVariable, vValue);    
      }  

	  void FreeCDof(IndexType NodeId,const Variable<double>& rVariable)
      {
           mpModeler->FreeCDof(NodeId, rVariable);    
      }  

      template<class TDataType>
      void SetSource(IndexType ElementId, 
                     SourceTarget source_target,
                     const Variable<TDataType>& rVariable, 
                     const TDataType& Value)
      {
          mpModeler->SetBoundarySource(ElementId, source_target, rVariable, Value);
      }
     
      template<class TDataType>
      void Print(const Variable<TDataType>& rVariable)
      {
          mpDefaultOutput->Print(rVariable);
      }

      template<class TDataType>
      void PrintOnGaussPoints(const Variable<TDataType>& rVariable, 
                              const String& ElementsGroup)
      {
          mpDefaultOutput->PrintOnGaussPoints(rVariable,ElementsGroup);
      }

      void PrintMesh(void)
      {
	      mpDefaultOutput->PrintMesh();
      }
     
      void CreateTimeStep();
     
      void CreateSolutionStep();

      template<class TDataType>
      void Smooth(const Variable<TDataType>& rVariable)
      {
      }

      template<class TDataType>
      void CalculateNodal(const Variable<TDataType>& rVariable)
      {
		   mpModeler->CalculateNodal(rVariable);
      }

      template<class TDataType>
      TDataType GetValue(const Variable<TDataType>& rVariable)
      {
           if(rVariable.Name() == "T") return mTime;         
           return mpModel->Value(rVariable);
      }

      void SetValue(int NodeId, 
                    const Variable<double>& rVariable, 
                    const double& rValue)
      {
           if(rVariable.Name() == "T") mTime = rValue;                  
           mpModel->Value(rVariable, *(mpModel->GetNode(NodeId))) = rValue;
      }
     
      void SetValue(const Variable<double>& rVariable, 
                    const double& rValue)
      {
           if(rVariable.Name() == "T") mTime = rValue;                  
           mpModel->Value(rVariable) = rValue;
      }
     
      void GenerateLinearSolver(const String& Name, 
                                const String& SolverType, 
                                int MaxNumberOfIterations, 
								int StepIterations,
                                double Tolerance);
     
      void GenerateSolvingStrategy(const String& Name, 
                                   const String& StrategyType, 
                                   const String& LinearSolverName);
     
      void Solve(const String& SolvingStrategyName, 
                 const String& ElementsGroup);

      void Execute(const String& ProcessName, 
                   const String& ElementsGroup);
     
      void Execute(const String& BaseFileName);

      void SetProblemType(const String& ProblemType);

      void SetFrequency(double ProblemFrequency);

      void IniBuilding();

	  void SweepFrequency(double IniFreq,
		                  double EndFreq,
						  double StepFreq);
     
    private:
          
      Model::Pointer mpModel;
      
      Modeler::Pointer mpModeler;
      
      Output::Pointer mpDefaultOutput;

	  double mTime;

      String mBaseFileName;

      bool mBuild;

	  bool mFSing2L, 
		   mFSing3L, 
		   mFSing4L,
	       mFSing5L,
	       mFSing6L;

	  bool mCalculateH;

	  bool mCalculateJ;

	  bool mCalculateJouleH;

	  bool mCreateHONodes;

	  bool mPrintHOMesh;

	  bool mFrequencySweep;

	  double mIniFreq, 
		     mEndFreq, 
			 mStepFreq;

	  unsigned int mElementIdCounter;

    }; /* Class ClassName */
  
  
}  /* namespace Kratos.*/


#endif /* KRATOS_KRATOS_KERNEL defined */

