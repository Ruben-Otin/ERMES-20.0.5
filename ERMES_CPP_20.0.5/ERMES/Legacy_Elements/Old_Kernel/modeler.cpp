
//***************************************************
//***************************************************

#include "../ERMES/VolumeElement_4th.h"
#include "../ERMES/JSource_4th.h"
#include "../ERMES/FarField_4th.h"
#include "../ERMES/RWPortTE10_4th.h"
#include "../ERMES/CoaxialPortTEM_4th.h"
#include "../ERMES/GenericRobin_4th.h"

#include "../ERMES/VolumeElement_3th.h"
#include "../ERMES/JSource_3th.h"
#include "../ERMES/FarField_3th.h"
#include "../ERMES/RWPortTE10_3th.h"
#include "../ERMES/CoaxialPortTEM_3th.h"
#include "../ERMES/GenericRobin_3th.h"

#include "../ERMES/VolumeElement_2nd.h"
#include "../ERMES/JSource_2nd.h"
#include "../ERMES/FarField_2nd.h"
#include "../ERMES/RWPortTE10_2nd.h"
#include "../ERMES/CoaxialPortTEM_2nd.h"
#include "../ERMES/GenericRobin_2nd.h"

#include "../ERMES/VolumeElement_2ndQ.h"
#include "../ERMES/JSource_2ndQ.h"
#include "../ERMES/CoaxialPortTEM_2ndQ.h"
#include "../ERMES/GenericRobin_2ndQ.h"
#include "../ERMES/FarField_2ndQ.h"

#include "../ERMES/VolumeElement_1st.h"
#include "../ERMES/JSource_1st.h"
#include "../ERMES/FarField_1st.h"
#include "../ERMES/RWPortTE10_1st.h"
#include "../ERMES/CoaxialPortTEM_1st.h"
#include "../ERMES/GenericRobin_1st.h"

#include "../ERMES/Utils/MatrixFileWritting.h"
#include "../ERMES/Utils/LagrangeBases.h"
#include "../ERMES/Utils/GaussIntegrationTables.h"

//***************************************************
//***************************************************

#include <algorithm>
#include <iostream>
#include <iomanip>
#include "../includes/modeler.h"
#include "../includes/model.h"
#include "../linear_solvers/kratos_linear_solvers.h"

//Para ejecutar achivos externos
#include <shlobj.h>
#include <shlwapi.h>

#include <stdlib.h>

namespace Kratos
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - Contact discontinuity in dielectrics surface
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //************************************************************************************************* 
    //* - Fill nodes properties vector map
    //************************************************************************************************* 
	void Modeler::ContactProperties(int* NodesId, unsigned int PropertiesId)
	{
		 for (int i=0;i<4;i++) mTempNodeProperties[NodesId[i]] = PropertiesId;	
	}

	//************************************************************************************************* 
    //* - Generate pairs of nodes in a contact surface.
    //************************************************************************************************* 
	void Modeler::GenerateContactPairs(Vector<int>& CtCNodesId)
	{   
		int i;
		
		Vector<int> sR_NodesId;
		Vector<int> sL_NodesId;

		unsigned short int R_PropertyId;
		unsigned short int L_PropertyId;

        //sR_NodesId will have the nodes with higher mod(eps) 
		//sL_NodesId will have the nodes with lower  mod(eps) 
		Properties::Pointer First3_prop = mpModel->GetProperties(mTempNodeProperties[CtCNodesId[0]]);
		Properties::Pointer Last3_prop  = mpModel->GetProperties(mTempNodeProperties[CtCNodesId[3]]);

		double freq = mProblemFrequency;
		double eo   = 8.8541878176e-12;
		
		double sigma_First3    = (*First3_prop)(IHL_ELECTRIC_CONDUCTIVITY);
		double real_eps_First3 = (*First3_prop)(REAL_ELECTRIC_PERMITTIVITY) * eo;
		double imag_eps_First3 = (*First3_prop)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		std::complex<double> cEps_First3( real_eps_First3, imag_eps_First3 + (sigma_First3/freq) );

		double mod_epc_First3 = std::abs(cEps_First3);

		double sigma_Last3     = (*Last3_prop)(IHL_ELECTRIC_CONDUCTIVITY);
		double real_eps_Last3  = (*Last3_prop)(REAL_ELECTRIC_PERMITTIVITY) * eo;
		double imag_eps_Last3  = (*Last3_prop)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		std::complex<double> cEps_Last3( real_eps_Last3, imag_eps_Last3 + (sigma_Last3/freq) );

		double mod_epc_Last3 = std::abs(cEps_Last3);

		double mod_epc_R, mod_epc_L;
 
		if (mod_epc_First3 > mod_epc_Last3)
		{
			Vector<int> OneSideNodesId(3);

			for (i=0;i<3;i++) OneSideNodesId[i] = CtCNodesId[i];
			PushHONodesOnSurface(OneSideNodesId, sR_NodesId);

			for (i=0;i<3;i++) OneSideNodesId[i] = CtCNodesId[i+3];
			PushHONodesOnSurface(OneSideNodesId, sL_NodesId);
       
			CalculateContactNormals(sR_NodesId);
			CalculateContactNormals(sL_NodesId);

			R_PropertyId = mTempNodeProperties[CtCNodesId[0]];
			L_PropertyId = mTempNodeProperties[CtCNodesId[3]];

			mod_epc_R = mod_epc_First3;
			mod_epc_L = mod_epc_Last3;
		}
		else
		{
			Vector<int> OneSideNodesId(3);

			for (i=0;i<3;i++) OneSideNodesId[i] = CtCNodesId[i];
			PushHONodesOnSurface(OneSideNodesId, sL_NodesId);

			for (i=0;i<3;i++) OneSideNodesId[i] = CtCNodesId[i+3];
			PushHONodesOnSurface(OneSideNodesId, sR_NodesId);

			CalculateContactNormals(sL_NodesId);
			CalculateContactNormals(sR_NodesId);

			L_PropertyId = mTempNodeProperties[CtCNodesId[0]];
			R_PropertyId = mTempNodeProperties[CtCNodesId[3]];

			mod_epc_L = mod_epc_First3;
			mod_epc_R = mod_epc_Last3;
		}

		Vector<int>::const_iterator R_Idit;
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for(R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit) 
		{
			int RNodeId = *R_Idit;
			int LNodeId = *L_Idit;
			
			Node::Pointer pRNode = mpModel->GetNode(RNodeId-1);
			Node::Pointer pLNode = mpModel->GetNode(LNodeId-1);
			
			if (pRNode->pDofcEx()->IsFixed())
			{   
				if (mTempContactPairs[RNodeId] != LNodeId)
				{
					int currLNodeId = mTempContactPairs[RNodeId];

					Properties::Pointer currL_prop  = mpModel->GetProperties(mTempNodeProperties[currLNodeId]);

					double sigma_currL    = (*currL_prop)(IHL_ELECTRIC_CONDUCTIVITY);
					double real_eps_currL = (*currL_prop)(REAL_ELECTRIC_PERMITTIVITY) * eo;
					double imag_eps_currL = (*currL_prop)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

					std::complex<double> cEps_currL( real_eps_currL, imag_eps_currL + (sigma_currL/freq) );

					double mod_epc_currL = std::abs(cEps_currL);

					if (mod_epc_currL > mod_epc_L)
					{
						Node::Pointer pcurrLNode = mpModel->GetNode(currLNodeId-1);
						
						pcurrLNode->pDofcEx()->FixDof();
						pcurrLNode->pDofcEy()->FreeDof();

						pLNode->pDofcEy()->FixDof();

						mTempContactPairs[RNodeId]     = LNodeId;
						mTempContactPairs[currLNodeId] = LNodeId;

						std::map<unsigned int, unsigned int>::iterator cPair_it;

						for(cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it) 
						{
							if ( cPair_it->second == currLNodeId ) mTempContactPairs[cPair_it->first] = LNodeId;
						}
					}
					else 
					{
						pLNode->pDofcEx()->FixDof();
						mTempContactPairs[LNodeId] = currLNodeId;
					}
				}
			}
			else if (pLNode->pDofcEx()->IsFixed())
			{
				if (!pRNode->pDofcEy()->IsFixed())
				{
					pRNode->pDofcEx()->FixDof();
					mTempContactPairs[RNodeId] = mTempContactPairs[LNodeId];
				}
			}
			else if (pRNode->pDofcEy()->IsFixed())
			{   
				pRNode->pDofcEx()->FixDof();
			    pRNode->pDofcEy()->FreeDof();

				pLNode->pDofcEy()->FixDof();

				mTempContactPairs[RNodeId] = LNodeId;

				std::map<unsigned int, unsigned int>::iterator cPair_it;

		        for(cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it) 
				{
					if ( cPair_it->second == RNodeId ) mTempContactPairs[cPair_it->first] = LNodeId;
				}
			}
			else
			{
				pRNode->pDofcEx()->FixDof();
				pLNode->pDofcEy()->FixDof();

				mTempContactPairs[RNodeId] = LNodeId;
			}

			mTempNodeProperties[RNodeId] = R_PropertyId;
		    mTempNodeProperties[LNodeId] = L_PropertyId;

			L_Idit++;
		}
	}

	//************************************************************************************************* 
    //* - Calculates correct contact normals in surfaces
    //*************************************************************************************************
    void Modeler::GenerateCorrectCTCNormals(Vector<int>& LONodesId)
    {  	
		if ( mQuadraticGeometry )
		{  
			GenerateCorrectCTCQNormals(LONodesId);
			return;
		}

		Node::Pointer pNode;
		std::vector<Node::Pointer> nodes;
		std::vector<int>::const_iterator id_it;

		std::vector<int> NodesId;

		PushHONodesOnSurface(LONodesId, NodesId);

		for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
		{
			 pNode = mpModel->GetNode((*id_it)-1);
			 nodes.push_back(pNode);
		}

		// Exterior normal of the element  
		Vector<double> normal(3);

		// Normal = (v2 x v1)/mod_normal
		normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
				    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
		normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
				    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
		normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
				    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());
        // Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				 mSetCorrectedCtCNormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				mSetCorrectedCtCNormals[*id_it].push_back(normal);
			}
		}
	}

	//************************************************************************************************* 
    //* - Calculates geometric averaged normal in contact surfaces
    //************************************************************************************************* 
	void Modeler::CalculateContactNormals(Vector<int>& NodesId)
	{
		if ( mQuadraticGeometry )
		{  
			CalculateContactQNormals(NodesId);
			return;
		}

		Node::Pointer pNode;
		std::vector<Node::Pointer> nodes;
		std::vector<int>::const_iterator id_it;

		for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
		{
			pNode = mpModel->GetNode((*id_it)-1);
			nodes.push_back(pNode);
		}

		// Exterior normal of the element  
		Vector<double> normal(3);

		// Normal = (v2 x v1)/mod_normal
		normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
				    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
		normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
				    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
		normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
				    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

		// Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				if (mSetCorrectedCtCNormals.find(*id_it) != mSetCorrectedCtCNormals.end()) continue;
				
				mSetContactNormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
			{
				if (mSetCorrectedCtCNormals.find(*id_it) != mSetCorrectedCtCNormals.end()) continue;

				mSetContactNormals[*id_it].push_back(normal);
			}
		}
	}

	//************************************************************************************************* 
    //* - Calculates average normal in contact nodes for quadratic meshing
    //*************************************************************************************************
    void Modeler::CalculateContactQNormals(std::vector<int>& LONodesId)
	{
		std::vector<int> NodesId;
		std::vector<int>::const_iterator id_it;

		PushHONodesOnSurface(LONodesId, NodesId);

		int numNodes = 6;
		int nN       = 0;
		int pt       = 0;

		std::vector<double> Xi(numNodes), Yi(numNodes), Zi(numNodes);

		for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
		{
			Node::Pointer pNode = mpModel->GetNode((*id_it)-1);

			Xi[nN] = pNode->X();
			Yi[nN] = pNode->Y();
			Zi[nN] = pNode->Z();

			nN++;
		}

		Matrix<double> dNda, dNdb;

		SurfDerivativesInNodes_2ndQ(dNda, dNdb);

		for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
		{
			if (mSetCorrectedCtCNormals.find(*id_it) != mSetCorrectedCtCNormals.end()) continue;

			// Exterior normal of the element  
		    Vector<double> normal(3);

			double dXda = 0.00;	double dXdb = 0.00;
			double dYda = 0.00;	double dYdb = 0.00;
			double dZda = 0.00;	double dZdb = 0.00;

			//dNi/da ,dNi/db
			for (nN=0; nN<numNodes; nN++) 
			{
				dXda += dNda[nN][pt]*Xi[nN];  dXdb += dNdb[nN][pt]*Xi[nN];
				dYda += dNda[nN][pt]*Yi[nN];  dYdb += dNdb[nN][pt]*Yi[nN];         
				dZda += dNda[nN][pt]*Zi[nN];  dZdb += dNdb[nN][pt]*Zi[nN];
			}

			normal[0] = dZda*dYdb - dYda*dZdb;
			normal[1] = dXda*dZdb - dZda*dXdb;
			normal[2] = dYda*dXdb - dXda*dYdb;

			if ( mNormalsAreaWtd )
			{
				mSetContactNormals[*id_it].push_back(normal);
			}
			else
			{
				normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));
				mSetContactNormals[*id_it].push_back(normal);
			}

			pt++;

		}//end for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
	}

	//************************************************************************************************* 
    //* - Calculates average normal in contact nodes for quadratic meshing
    //*************************************************************************************************
    void Modeler::GenerateCorrectCTCQNormals(std::vector<int>& LONodesId)
	{
		std::vector<int> NodesId;
		std::vector<int>::const_iterator id_it;

		PushHONodesOnSurface(LONodesId, NodesId);

		int numNodes = 6;
		int nN       = 0;
		int pt       = 0;

		std::vector<double> Xi(numNodes), Yi(numNodes), Zi(numNodes);

		for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
		{
			Node::Pointer pNode = mpModel->GetNode((*id_it)-1);

			Xi[nN] = pNode->X();
			Yi[nN] = pNode->Y();
			Zi[nN] = pNode->Z();

			nN++;
		}

		Matrix<double> dNda, dNdb;

		SurfDerivativesInNodes_2ndQ(dNda, dNdb);

		for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
		{
			// Exterior normal of the element  
		    Vector<double> normal(3);

			double dXda = 0.00;	double dXdb = 0.00;
			double dYda = 0.00;	double dYdb = 0.00;
			double dZda = 0.00;	double dZdb = 0.00;

			//dNi/da ,dNi/db
			for (nN=0; nN<numNodes; nN++) 
			{
				dXda += dNda[nN][pt]*Xi[nN];  dXdb += dNdb[nN][pt]*Xi[nN];
				dYda += dNda[nN][pt]*Yi[nN];  dYdb += dNdb[nN][pt]*Yi[nN];         
				dZda += dNda[nN][pt]*Zi[nN];  dZdb += dNdb[nN][pt]*Zi[nN];
			}

			normal[0] = dZda*dYdb - dYda*dZdb;
			normal[1] = dXda*dZdb - dZda*dXdb;
			normal[2] = dYda*dXdb - dXda*dYdb;

			if ( mNormalsAreaWtd )
			{
				mSetCorrectedCtCNormals[*id_it].push_back(normal);
			}
			else
			{
				normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

				mSetCorrectedCtCNormals[*id_it].push_back(normal);
			}

			pt++;

		}//end for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
	}

	//************************************************************************************************* 
    //* - Set average normals in contact surfaces
    //************************************************************************************************* 
	void Modeler::SetAvgNormalsInContact()
	{
		//Setting average normals in surface dilelectric nodes
		std::map<unsigned int, std::vector<std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetContactNormals.begin(); it_Set != mSetContactNormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                Vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

            mContactNormals[it_Set->first] = ac_normal; 
        } 

		mSetContactNormals.clear();
        mSetContactNormals.swap(std::map<unsigned int,  std::vector<std::vector<double> > >(mSetContactNormals));

		for(it_Set = mSetCorrectedCtCNormals.begin(); it_Set != mSetCorrectedCtCNormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                Vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

            mContactNormals[it_Set->first] = ac_normal; 
        } 

		mSetCorrectedCtCNormals.clear();
        mSetCorrectedCtCNormals.swap(std::map<unsigned int,  std::vector<std::vector<double> > >(mSetCorrectedCtCNormals));
	}

	//************************************************************************************************* 
    //* - Setting discontinuity in dielectric surface.
    //************************************************************************************************* 
	void Modeler::SetContactDiscontinuity()
	{
		SetAvgNormalsInContact();

		//Setting data pairs struct 
		std::map<unsigned int, unsigned int>::iterator cPair_it;

		for(cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it) 
        {
			unsigned int RNodeId = cPair_it->first;
			unsigned int LNodeId = mTempContactPairs[RNodeId];

			Node::Pointer pRNode = mpModel->GetNode(RNodeId-1);
			Node::Pointer pLNode = mpModel->GetNode(LNodeId-1);

			(*mpModel)(cEx, *pRNode) = czero;
			(*mpModel)(cEy, *pRNode) = czero;
			(*mpModel)(cEz, *pRNode) = czero;
			
			pRNode->pDofcEx()->FixDof();
			pRNode->pDofcEy()->FixDof();
			pRNode->pDofcEz()->FixDof();
				
			pLNode->pDofcEx()->FreeDof();
			pLNode->pDofcEy()->FreeDof();
			pLNode->pDofcEz()->FreeDof();

			Vector<double> AvgNormal(3,0.00);
			Vector<double> RNormal(mContactNormals[RNodeId]);
			Vector<double> LNormal(mContactNormals[LNodeId]);
            
			////////////////////////////////////////////////////////////////////////
            // - You can take any of the normals, 
			//   n or -n produces the same transformation matrix.
			////////////////////////////////////////////////////////////////////////
			AvgNormal[0] = -RNormal[0];
			AvgNormal[1] = -RNormal[1];
			AvgNormal[2] = -RNormal[2];

			//Making a contact pair
			ContactPairData& pContactPair = mContactPairs[RNodeId];
			
			pContactPair.ctcNormal = AvgNormal;
			pContactPair.ctcNode   = LNodeId;
			pContactPair.matR      = mTempNodeProperties[RNodeId] ;
			pContactPair.matL      = mTempNodeProperties[LNodeId] ;
		}
		
		DelTempContactVectors();

		/*//////////////////////////////////////////////////////////////////////////////////////
		std::map<unsigned int, ContactPairData>::iterator it;
		for(it = mContactPairs.begin(); it != mContactPairs.end(); ++it) 
		{
			std::cout<<"Node1: "<<it->first            <<"  mat1: "<<(it->second).matR<<"   "
					 <<"Node2: "<<(it->second).ctcNode <<"  mat2: "<<(it->second).matL<<"   "
					 <<"n :"<<(it->second).ctcNormal
					 <<std::endl;
		}
		/////////////////////////////////////////////////////////////////////////////////////*/
	}

	//************************************************************************************************* 
    //* - Delete temporal vectors generated along contact discontinuity making.
    //************************************************************************************************* 
	void Modeler::DelTempContactVectors()
	{
		mTempNodeProperties.clear();
        mTempNodeProperties.swap(std::map<unsigned int, unsigned short int>(mTempNodeProperties));

		mContactNormals.clear();
        mContactNormals.swap(std::map<unsigned int, std::vector<double> >(mContactNormals));

		mTempContactPairs.clear();
		mTempContactPairs.swap(std::map<unsigned int, unsigned int>(mTempContactPairs));
	}

	//************************************************************************************************* 
    //* - Adding normals to second node from the first node
    //************************************************************************************************* 
	void Modeler::AddContactPECNormals()
	{	
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		//Adding normal PEC to second node in dielectric pair
		for(itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc) 
        {
			unsigned int FirstNode  = itCtc->first;
	
			if ( mSetPECnormals.find(FirstNode) != mSetPECnormals.end() )
			{
				unsigned int SecondNode = (itCtc->second).ctcNode;
				
				std::vector<std::vector<double> > FirstVectorSet (mSetPECnormals[FirstNode]);
				std::vector<std::vector<double> >::iterator it_Fvs;

				for(it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs)
				{
					Vector<double> vectorF(*it_Fvs);
					mSetPECnormals[SecondNode].push_back(vectorF);
				}
			}

			if ( mSetSPECnormals.find(FirstNode) != mSetSPECnormals.end() ) 
			{
				unsigned int SecondNode = (itCtc->second).ctcNode;
				
				std::vector<std::vector<double> > FirstVectorSet (mSetSPECnormals[FirstNode]);
				std::vector<std::vector<double> >::iterator it_Fvs;

				for(it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs)
				{
					Vector<double> vectorF(*it_Fvs);
					mSetSPECnormals[SecondNode].push_back(vectorF);
				}
			}
		}
	}

	//************************************************************************************************* 
    //* - Force PEC normal in all contact nodes that are in a PEC surface
    //*************************************************************************************************
	void Modeler::ForcePECNormalsInContact()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		//All nodes in an intersection with the same PEC normal
		//New dielectric normal. Taking into account PEC, who dominates normal direction in singularity.
        for(itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc) 
        {
			unsigned int FirstNode  = itCtc->first;
			unsigned int SecondNode = (itCtc->second).ctcNode;

			if ( (mNormals.find(SecondNode) != mNormals.end()) )
            {
				Vector<double> PECNormal(mNormals[SecondNode]);
				Vector<double> NEWctcNormal(3,0.00);
				Vector<double> Dummy(3,0.00);

				TangencialCoordinates(PECNormal, NEWctcNormal, Dummy);

				mNormals[FirstNode] = PECNormal;

				mContactPairs[FirstNode].ctcNormal = NEWctcNormal;
			}
		}
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - Fill singularity maps
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//***********************************************************
    //* - This functions are because of an error in GiD.
	//* - If you mesh using quadrilateral contact elements, 
    //*   extra nodes which are not in any element appear.  
    //***********************************************************
	void Modeler::IniEraseExtraNodes()
	{
		if ( !mQuadraticGeometry ) return;

		unsigned int numNodes = (mpModel->GetNodesArray()).size();

		mUsefulNodes.resize(numNodes);

		for ( int i=0; i<numNodes; i++ ) mUsefulNodes[i] = false;
	}

	//***********************************************************
    //* - This functions are because of an error in GiD.
	//* - If you mesh using quadrilateral contact elements, 
    //*   extra nodes which are not in any element appear.  
    //***********************************************************
	void Modeler::EndEraseExtraNodes()
	{
		if ( !mQuadraticGeometry ) return;

		unsigned int numNodes = (mpModel->GetNodesArray()).size();

		for ( int i=0; i<numNodes; i++ ) 
		{	
			if ( !mUsefulNodes[i] )
			{
				Node::Pointer pNode = mpModel->GetNode(i);

				(*mpModel)(cEx, *pNode) = czero;
				(*mpModel)(cEy, *pNode) = czero;
				(*mpModel)(cEz, *pNode) = czero;
				
				pNode->pDofcEx()->FixDof();
				pNode->pDofcEy()->FixDof();
				pNode->pDofcEz()->FixDof();
			}
		}
			
		mUsefulNodes.clear();
        mUsefulNodes.swap(std::vector<bool>(mUsefulNodes));
	}

    //***************************************
    //* - Fill mNSingular2L map
    //***************************************
	void Modeler::FindSing2L(int* NodesId)
	{
		// Read header in function Modeler::IniEraseExtraNodes()
		if ( mQuadraticGeometry )
		{
			std::vector<int> HONodesId;
			PushHONodesOnVolume(NodesId,HONodesId);
			std::vector<int>::iterator hoit;
			for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  mUsefulNodes[(*hoit)-1] = true;
		}

		int i;
		double maxSingValue = 1.00;
		
		for(i=0; i<4; i++)
        {
            if ( mSingular.find(NodesId[i]) != mSingular.end() )
            {	
				if (mSingular[NodesId[i]] > maxSingValue) 
				{
				    maxSingValue = mSingular[NodesId[i]];
				}
            }
        }

        if ( maxSingValue > 1.00 )
        {    
            for (i=0; i<4; i++) 
            {
                if ( mSingular.find(NodesId[i]) == mSingular.end() ) 
				{
				    mNSingular2L[NodesId[i]] = maxSingValue; 
				}
            }
        }
	}

	//***************************************
    //* - Fill mNSingular3L map
    //***************************************
	void Modeler::FindSing3L(int* NodesId)
	{
		int i;
		
		for(i=0; i<4; i++)
        {
            if ( mSingular.find(NodesId[i]) != mSingular.end() ) return;
        }
		
		double maxSingValue = 2.00;
		
		for(i=0; i<4; i++)
        {
            if ( mNSingular2L.find(NodesId[i]) != mNSingular2L.end() )
            {	
				if (mNSingular2L[NodesId[i]] > maxSingValue) 
				{
				    maxSingValue = mNSingular2L[NodesId[i]];
				}
            }
        }

        if ( maxSingValue > 2.00 )
        {    
            for (i=0; i<4; i++) 
            {
                if ( mNSingular2L.find(NodesId[i]) == mNSingular2L.end() ) 
				{
				    mNSingular3L[NodesId[i]] = maxSingValue; 
				}
            }
        }       
	}

	//***************************************
    //* - Fill mNSingular4L map
    //***************************************
	void Modeler::FindSing4L(int* NodesId)
	{
		int i;
		
		for(i=0; i<4; i++)
        {
            if ( mSingular.find(NodesId[i])    != mSingular.end()    ) return;
			if ( mNSingular2L.find(NodesId[i]) != mNSingular2L.end() ) return;
        }
		
		double maxSingValue = 3.00;
		
		for(i=0; i<4; i++)
        {
            if ( mNSingular3L.find(NodesId[i]) != mNSingular3L.end() )
            {	
				if (mNSingular3L[NodesId[i]] > maxSingValue) 
				{
				    maxSingValue = mNSingular3L[NodesId[i]];
				}
            }
        }

        if ( maxSingValue > 3.00 )
        {    
            for (i=0; i<4; i++) 
            {
                if ( mNSingular3L.find(NodesId[i]) == mNSingular3L.end() ) 
				{
				    mNSingular4L[NodesId[i]] = maxSingValue; 
				}
            }
        } 	
	}

	//***************************************
    //* - Fill mNSingular5L map
    //***************************************
	void Modeler::FindSing5L(int* NodesId)
	{
		int i;
		
		for(i=0; i<4; i++)
        {
            if ( mSingular.find(NodesId[i])    != mSingular.end()    ) return;
			if ( mNSingular2L.find(NodesId[i]) != mNSingular2L.end() ) return;
			if ( mNSingular3L.find(NodesId[i]) != mNSingular3L.end() ) return;
        }
		
		double maxSingValue = 4.00;
		
		for(i=0; i<4; i++)
        {
            if ( mNSingular4L.find(NodesId[i]) != mNSingular4L.end() )
            {	
				if (mNSingular4L[NodesId[i]] > maxSingValue) 
				{
				    maxSingValue = mNSingular4L[NodesId[i]];
				}
            }
        }

        if ( maxSingValue > 4.00 )
        {    
            for (i=0; i<4; i++) 
            {
                if ( mNSingular4L.find(NodesId[i]) == mNSingular4L.end() ) 
				{
				    mNSingular5L[NodesId[i]] = maxSingValue; 
				}
            }
        } 	
	}

	//***************************************
    //* - Fill mNSingular6L map
    //***************************************
	void Modeler::FindSing6L(int* NodesId)
	{
		int i;
		
		for(i=0; i<4; i++)
        {
            if ( mSingular.find(NodesId[i])    != mSingular.end()    ) return;
			if ( mNSingular2L.find(NodesId[i]) != mNSingular2L.end() ) return;
			if ( mNSingular3L.find(NodesId[i]) != mNSingular3L.end() ) return;
			if ( mNSingular4L.find(NodesId[i]) != mNSingular4L.end() ) return;
        }
		
		double maxSingValue = 5.00;
		
		for(i=0; i<4; i++)
        {
            if ( mNSingular5L.find(NodesId[i]) != mNSingular5L.end() )
            {	
				if (mNSingular5L[NodesId[i]] > maxSingValue) 
				{
				    maxSingValue = mNSingular5L[NodesId[i]];
				}
            }
        }

        if ( maxSingValue > 5.00 )
        {    
            for (i=0; i<4; i++) 
            {
                if ( mNSingular5L.find(NodesId[i]) == mNSingular5L.end() ) 
				{
				    mNSingular6L[NodesId[i]] = maxSingValue; 
				}
            }
        } 	
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - High Order Elements
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//************************************************************************************************* 
    //* - Add high order elements to the HONodesId vector in volume elements
    //*************************************************************************************************
	void Modeler::PushHONodesOnVolume(int* NodesId, std::vector<int>& HONodesId)
	{
		HONodesId.push_back(NodesId[0]);
        HONodesId.push_back(NodesId[1]);
        HONodesId.push_back(NodesId[2]);
		HONodesId.push_back(NodesId[3]);

		if ( mElementOrder == 1 ) return;

        HONodesOnEdge(NodesId[0], NodesId[1], HONodesId);
        HONodesOnEdge(NodesId[1], NodesId[2], HONodesId);
        HONodesOnEdge(NodesId[2], NodesId[0], HONodesId);
		HONodesOnEdge(NodesId[0], NodesId[3], HONodesId);
        HONodesOnEdge(NodesId[1], NodesId[3], HONodesId);
        HONodesOnEdge(NodesId[2], NodesId[3], HONodesId);

		if ( mElementOrder == 2 ) return;

		HONodesOnFace(NodesId[0], NodesId[1], NodesId[2], HONodesId);
		HONodesOnFace(NodesId[0], NodesId[1], NodesId[3], HONodesId);
		HONodesOnFace(NodesId[1], NodesId[2], NodesId[3], HONodesId);
		HONodesOnFace(NodesId[0], NodesId[2], NodesId[3], HONodesId);

		if ( mElementOrder == 3 ) return;
			
		HONodesOnVolum(NodesId, HONodesId);


		/* 4th order /////////////////////////////////////////////////////////////////////////////////////////
	    std::cout<<"Points:"<<std::endl;
	    std::cout<<HONodesId[0]<<" , "<<HONodesId[1]<<" , "<<HONodesId[2]<<" , "<<HONodesId[3]<<std::endl;

	    std::cout<<"Edges:"<<std::endl;
	    std::cout<<"[0 1] : "<<HONodesId[4] <<" , "<<HONodesId[5] <<" , "<<HONodesId[6] <<std::endl;
	    std::cout<<"[1 2] : "<<HONodesId[7] <<" , "<<HONodesId[8] <<" , "<<HONodesId[9] <<std::endl;
	    std::cout<<"[2 0] : "<<HONodesId[10]<<" , "<<HONodesId[11]<<" , "<<HONodesId[12]<<std::endl;
	    std::cout<<"[0 3] : "<<HONodesId[13]<<" , "<<HONodesId[14]<<" , "<<HONodesId[15]<<std::endl;
	    std::cout<<"[1 3] : "<<HONodesId[16]<<" , "<<HONodesId[17]<<" , "<<HONodesId[18]<<std::endl;
	    std::cout<<"[2 3] : "<<HONodesId[19]<<" , "<<HONodesId[20]<<" , "<<HONodesId[21]<<std::endl;

	    std::cout<<"Faces:"<<std::endl;
	    std::cout<<"[0 1 2] : "<<HONodesId[22]<<" , "<<HONodesId[23]<<" , "<<HONodesId[24]<<std::endl;
	    std::cout<<"[0 1 3] : "<<HONodesId[25]<<" , "<<HONodesId[26]<<" , "<<HONodesId[27]<<std::endl;
	    std::cout<<"[1 2 3] : "<<HONodesId[28]<<" , "<<HONodesId[29]<<" , "<<HONodesId[30]<<std::endl;
	    std::cout<<"[0 2 3] : "<<HONodesId[31]<<" , "<<HONodesId[32]<<" , "<<HONodesId[33]<<std::endl;

	    std::cout<<"Volume:"<<std::endl;
	    std::cout<<"[0 1 2 3] : "<<HONodesId[34]<<std::endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////*/

		/* 3th order /////////////////////////////////////////////////////////////////////////////////////////
	    std::cout<<"Points:"<<std::endl;
	    std::cout<<HONodesId[0]<<" , "<<HONodesId[1]<<" , "<<HONodesId[2]<<" , "<<HONodesId[3]<<std::endl;

	    std::cout<<"Edges:"<<std::endl;
	    std::cout<<"[0 1] : "<<HONodesId[4] <<" , "<<HONodesId[5] <<std::endl;
	    std::cout<<"[1 2] : "<<HONodesId[6] <<" , "<<HONodesId[7] <<std::endl;
	    std::cout<<"[2 0] : "<<HONodesId[8] <<" , "<<HONodesId[9] <<std::endl;
	    std::cout<<"[0 3] : "<<HONodesId[10]<<" , "<<HONodesId[11]<<std::endl;
	    std::cout<<"[1 3] : "<<HONodesId[12]<<" , "<<HONodesId[13]<<std::endl;
	    std::cout<<"[2 3] : "<<HONodesId[14]<<" , "<<HONodesId[15]<<std::endl;

	    std::cout<<"Faces:"<<std::endl;
	    std::cout<<"[0 1 2] : "<<HONodesId[16]<<std::endl;
	    std::cout<<"[0 1 3] : "<<HONodesId[17]<<std::endl;
	    std::cout<<"[1 2 3] : "<<HONodesId[18]<<std::endl;
	    std::cout<<"[0 2 3] : "<<HONodesId[19]<<std::endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////*/

		/* 2nd order ////////////////////////////////////////////////////////////////////////////////////////
	    std::cout<<"Points:"<<std::endl;
	    std::cout<<HONodesId[0]<<" , "<<HONodesId[1]<<" , "<<HONodesId[2]<<" , "<<HONodesId[3]<<std::endl;

	    std::cout<<"Edges:"<<std::endl;
	    std::cout<<"[0 1] : "<<HONodesId[4]<<std::endl;
	    std::cout<<"[1 2] : "<<HONodesId[5]<<std::endl;
	    std::cout<<"[2 0] : "<<HONodesId[6]<<std::endl;
	    std::cout<<"[0 3] : "<<HONodesId[7]<<std::endl;
	    std::cout<<"[1 3] : "<<HONodesId[8]<<std::endl;
	    std::cout<<"[2 3] : "<<HONodesId[9]<<std::endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}

	//************************************************************************************************* 
    //* - Add high order elements to the HONodesId vector in surface elements
    //*************************************************************************************************
	void Modeler::PushHONodesOnSurface(int* NodesId, std::vector<int>& HONodesId)
	{
		HONodesId.push_back(NodesId[0]);
        HONodesId.push_back(NodesId[1]);
        HONodesId.push_back(NodesId[2]);

		if ( mElementOrder == 1 ) return;
  
        HONodesOnEdge(NodesId[0], NodesId[1], HONodesId);
        HONodesOnEdge(NodesId[1], NodesId[2], HONodesId);
        HONodesOnEdge(NodesId[2], NodesId[0], HONodesId);

        if ( mElementOrder == 2 ) return;
			
		HONodesOnFace(NodesId[0], NodesId[1], NodesId[2], HONodesId);
	}

	//************************************************************************************************* 
    //* - Add high order elements to the HONodesId vector in surface elements
    //*************************************************************************************************
	void Modeler::PushHONodesOnSurface(std::vector<int>& NodesId, std::vector<int>& HONodesId)
	{
		HONodesId.push_back(NodesId[0]);
        HONodesId.push_back(NodesId[1]);
        HONodesId.push_back(NodesId[2]);

		if ( mElementOrder == 1 ) return;
  
        HONodesOnEdge(NodesId[0], NodesId[1], HONodesId);
        HONodesOnEdge(NodesId[1], NodesId[2], HONodesId);
        HONodesOnEdge(NodesId[2], NodesId[0], HONodesId);

        if ( mElementOrder == 2 ) return;
			
		HONodesOnFace(NodesId[0], NodesId[1], NodesId[2], HONodesId);
	}

	//************************************************************************************************* 
    //* - Looking for high order nodes on edge
    //************************************************************************************************* 
    void Modeler::HONodesOnEdge(int IdNode1, int IdNode2, std::vector<int>& HONodesId)
    {
         bool inverseOrder;
         unsigned int minId, maxId;
         
         if (IdNode1 < IdNode2) 
         {
            inverseOrder = false;
            minId = IdNode1;
            maxId = IdNode2; 
         }
         else
         {
            inverseOrder = true;
            minId = IdNode2;
            maxId = IdNode1; 
         }

         std::vector<std::vector<unsigned int> > SameKeyEdges(mEdgeHONodes[minId*maxId]);
         std::vector<std::vector<unsigned int> >::iterator it_edges;

         for(it_edges = SameKeyEdges.begin(); it_edges != SameKeyEdges.end(); ++it_edges)
         {
             std::vector<unsigned int> edge(*it_edges);
     
             if ( (edge[0] == minId) && (edge[1] == maxId) )
             {
                  if (!inverseOrder)
                  {
                      std::vector<unsigned int>::iterator it_hon = edge.begin()+2;

                      while( it_hon != edge.end() )
                      {
                           HONodesId.push_back(*it_hon);
                           it_hon++;
                      }
                  }
                  else
                  {
                      std::vector<unsigned int>::reverse_iterator rit_hon = edge.rbegin();
                      
                      while( rit_hon != edge.rend()-2 )
                      {  
                          HONodesId.push_back(*rit_hon);
                          rit_hon++;
                      }
                  }
             }
         }
    }

	//************************************************************************************************* 
    //* - Looking for high order nodes on face
    //************************************************************************************************* 
	void Modeler::HONodesOnFace (int IdNode1, int IdNode2, int IdNode3, std::vector<int>& HONodesId)
	{
         unsigned int minId, medId, maxId;
		 int orderingCase;
		
		 if (IdNode1 < IdNode2) 
		 {
		      if (IdNode2 < IdNode3)
			  { 
				   //[n1,n2,n3]
			       orderingCase = 1;
				   minId = IdNode1;
				   medId = IdNode2; 
				   maxId = IdNode3; 
			  }
			  else
			  {
			       if (IdNode1 < IdNode3)
				   {
					    //[n1,n3,n2]
			            orderingCase = 2;
					    minId = IdNode1;
					    medId = IdNode3; 
					    maxId = IdNode2; 
				   }
				   else
				   {
					    //[n3,n1,n2]
			            orderingCase = 3;
					    minId = IdNode3;
				        medId = IdNode1; 
				        maxId = IdNode2; 
				   }
			  }
		 }
		 else
		 {
			  if (IdNode1 < IdNode3)
			  {
				   //[n2,n1,n3]
			       orderingCase = 4;
				   minId = IdNode2;
				   medId = IdNode1; 
				   maxId = IdNode3; 
			  }
			  else
			  {
				   if (IdNode2 < IdNode3)
				   {
					    //[n2,n3,n1]
			            orderingCase = 5;
					    minId = IdNode2;
					    medId = IdNode3; 
					    maxId = IdNode1; 
				   }
				   else
				   {
					    //[n3,n2,n1]
			            orderingCase = 6;
					    minId = IdNode3;
					    medId = IdNode2; 
					    maxId = IdNode1; 
				   }
			  }
		 }

		 std::vector<std::vector<unsigned int> > SameKeyFaces(mFaceHONodes[minId*medId*maxId]);
         std::vector<std::vector<unsigned int> >::iterator it_faces;

		 for(it_faces = SameKeyFaces.begin(); it_faces != SameKeyFaces.end(); ++it_faces)
         {
             std::vector<unsigned int> face(*it_faces);
     
             if ( (face[0]==minId) && (face[1]==medId) && (face[2]==maxId) )
             {
				 if      ( mElementOrder == 3 ) HONodesId.push_back(face[3]);
				 else if ( mElementOrder == 4 ) ReorderHONodesOnFace_4th(orderingCase, face, HONodesId);
             }
         }
	}

	//************************************************************************************************* 
    //* - Reordering nodes for a giving surface order
    //************************************************************************************************* 
	void Modeler::ReorderHONodesOnFace_4th(int orderingCase, std::vector<unsigned int>& face, std::vector<int>& HONodesId)
	{
		unsigned int minf = face[3];
		unsigned int medf = face[4];
		unsigned int maxf = face[5];

		 //[n1,n2,n3] = [minId, medId, maxId] 
		 if ( orderingCase == 1 )
		 {
			 HONodesId.push_back(minf);
			 HONodesId.push_back(medf);
			 HONodesId.push_back(maxf);
		 }
		 //[n1,n3,n2] = [minId, maxId, medId] 
		 else if ( orderingCase == 2 )
		 {
			 HONodesId.push_back(minf);
			 HONodesId.push_back(maxf);
			 HONodesId.push_back(medf);
		 }
		 //[n3,n1,n2] = [medId, maxId, minId] 
		 else if ( orderingCase == 3 )
		 {
			 HONodesId.push_back(medf);
			 HONodesId.push_back(maxf);
			 HONodesId.push_back(minf);
		 }
		 //[n2,n1,n3] = [medId, minId, maxId] 
		 else if ( orderingCase == 4 )
		 {
			 HONodesId.push_back(medf);
			 HONodesId.push_back(minf);
			 HONodesId.push_back(maxf);
		 }
		 //[n2,n3,n1] = [maxId, minId, medId] 
		 else if ( orderingCase == 5 )
		 {
			 HONodesId.push_back(maxf);
			 HONodesId.push_back(minf);
			 HONodesId.push_back(medf);
		 }
		 //[n3,n2,n1] = [maxId, medId, minId] 
		 else if ( orderingCase == 6 )
		 {
			 HONodesId.push_back(maxf);
			 HONodesId.push_back(medf);
			 HONodesId.push_back(minf);
		 }
	}

	//************************************************************************************************* 
    //* - Looking for high order nodes on volume
    //************************************************************************************************* 
	void Modeler::HONodesOnVolum(int* NodesId, std::vector<int>& HONodesId)
	{
	     unsigned int key = NodesId[0]*NodesId[1]*NodesId[2]*NodesId[3];

         std::vector<std::vector<unsigned int> > SameKeyElement(mVolmHONodes[key]);
         std::vector<std::vector<unsigned int> >::iterator it_elements;

         for(it_elements = SameKeyElement.begin(); it_elements != SameKeyElement.end(); ++it_elements)
         {
             std::vector<unsigned int> nInElement(*it_elements);
     
             if ( (nInElement[0]==NodesId[0]) && (nInElement[1]==NodesId[1]) && (nInElement[2]==NodesId[2]) && (nInElement[3]==NodesId[3]) )
             {
				 if ( mElementOrder == 4 ) HONodesId.push_back(nInElement[4]);
             }
         }
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - High order nodes creation
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//************************************************************************************************* 
    //* - Create high order nodes
    //*************************************************************************************************
	void Modeler::CreateHONodes(int* NodesId)
	{
		if ( mElementOrder == 1 ) return; 
		
		CreateHONodes_OnEdges(NodesId);

		if ( mElementOrder == 2 ) return; 

		CreateHONodes_OnFaces(NodesId);
		    
		if ( mElementOrder == 3 ) return;
			
		CreateHONodes_OnVolum(NodesId);
	}

    //************************************************************************************************* 
    //* - Create high order nodes if quadratic meshing (only for 2nd order elements)
    //*************************************************************************************************
	void Modeler::GenerateVEQElement(std::vector<int>& QHONodesId)
	{
		int i, j;
		
		unsigned int IdNode1,IdNode2;
        unsigned int key;

		unsigned int middleNode[6] = {QHONodesId[4],QHONodesId[6],QHONodesId[7],
		                              QHONodesId[5],QHONodesId[8],QHONodesId[9]};

		// If PropertyId == 0 throw error, material not assigned
		if (QHONodesId[10] == 0) { throw 0; }

		// Save material properties for contact elements
		for ( i=0; i<4; i++ ) mTempNodeProperties[QHONodesId[i]] = QHONodesId[10];

        int midNodeIndex = 0;

		for ( i=0; i<3; i++ )
		{
			for ( j=i+1; j<4; j++ )
			{
				if (QHONodesId[i] < QHONodesId[j]) 
				{
					IdNode1 = QHONodesId[i];
					IdNode2 = QHONodesId[j]; 
				}
				else
				{
					IdNode1 = QHONodesId[j];
					IdNode2 = QHONodesId[i]; 
				}

				key = IdNode1*IdNode2;

				if ( mEdgeHONodes.find(key) != mEdgeHONodes.end() )
				{
					std::vector<std::vector<unsigned int> > ExistantNodesOnEdges(mEdgeHONodes[key]);
					std::vector<std::vector<unsigned int> >::iterator it_ee;

					int repe = 0;

					for( it_ee = ExistantNodesOnEdges.begin(); it_ee != ExistantNodesOnEdges.end(); ++it_ee )
					{
						std::vector<unsigned int> NodesOnEdge(*it_ee);
                
						if ( (IdNode1 == NodesOnEdge[0]) && (IdNode2 == NodesOnEdge[1]) ) repe++;
					}

					if ( repe == 0 )
					{
						std::vector<unsigned int> edgeNodes(3);

						edgeNodes[0] = IdNode1;
						edgeNodes[1] = IdNode2;
						edgeNodes[2] = middleNode[midNodeIndex];

						mEdgeHONodes[IdNode1*IdNode2].push_back(edgeNodes);  
					}
				}
				else
				{
					std::vector<unsigned int> edgeNodes(3);

					edgeNodes[0] = IdNode1;
					edgeNodes[1] = IdNode2;
					edgeNodes[2] = middleNode[midNodeIndex];

					mEdgeHONodes[IdNode1*IdNode2].push_back(edgeNodes);  
				}
				
				midNodeIndex++;
			}
		}
	}

	//************************************************************************************************* 
    //* - Create High order nodes on edges
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnEdges(int* NodesId)
	{
		unsigned int IdNode1,IdNode2;
        unsigned int key;
		
		for (int i=0; i<3; i++)
		{
			for (int j=i+1; j<4; j++)
			{
				if (NodesId[i] < NodesId[j]) 
				{
					IdNode1 = NodesId[i];
					IdNode2 = NodesId[j]; 
				}
				else
				{
					IdNode1 = NodesId[j];
					IdNode2 = NodesId[i]; 
				}

				key = IdNode1*IdNode2;

				if ( mEdgeHONodes.find(key) != mEdgeHONodes.end() )
				{
					std::vector<std::vector<unsigned int> > ExistantNodesOnEdges(mEdgeHONodes[key]);
					std::vector<std::vector<unsigned int> >::iterator it_ee;

					int repe = 0;

					for( it_ee = ExistantNodesOnEdges.begin(); it_ee != ExistantNodesOnEdges.end(); ++it_ee )
					{
						std::vector<unsigned int> NodesOnEdge(*it_ee);
                
						if ( (IdNode1 == NodesOnEdge[0]) && (IdNode2 == NodesOnEdge[1]) ) repe++;
					}

					if ( repe == 0 )
					{
						if      ( mElementOrder == 2 ) CreateHONodes_OnEdge_2nd(IdNode1, IdNode2);
						else if ( mElementOrder == 3 ) CreateHONodes_OnEdge_3th(IdNode1, IdNode2);
						else if ( mElementOrder == 4 ) CreateHONodes_OnEdge_4th(IdNode1, IdNode2);
					}
				}
				else
				{
					if      ( mElementOrder == 2 ) CreateHONodes_OnEdge_2nd(IdNode1, IdNode2);
					else if ( mElementOrder == 3 ) CreateHONodes_OnEdge_3th(IdNode1, IdNode2);
					else if ( mElementOrder == 4 ) CreateHONodes_OnEdge_4th(IdNode1, IdNode2);
				}	
			}
		}
	}

	//************************************************************************************************* 
    //* - Create High order nodes on faces
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnFaces(int* NodesId)
	{
		unsigned int IdNode1,IdNode2,IdNode3;
        unsigned int key;
		
		for (int i=0; i<2; i++)
		{
			for (int j=i+1; j<3; j++)
			{
				for (int m=j+1; m<4; m++)
				{
					if (NodesId[i] < NodesId[j]) 
					{
						if (NodesId[j] < NodesId[m])
						{
							IdNode1 = NodesId[i];
						    IdNode2 = NodesId[j]; 
							IdNode3 = NodesId[m]; 
						}
						else
						{
							if (NodesId[i] < NodesId[m])
							{
								IdNode1 = NodesId[i];
								IdNode2 = NodesId[m]; 
								IdNode3 = NodesId[j]; 
							}
							else
							{
								IdNode1 = NodesId[m];
								IdNode2 = NodesId[i]; 
								IdNode3 = NodesId[j]; 
							}
						}
					}
					else
					{
						if (NodesId[i] < NodesId[m])
						{
							IdNode1 = NodesId[j];
						    IdNode2 = NodesId[i]; 
							IdNode3 = NodesId[m]; 
						}
						else
						{
							if (NodesId[j] < NodesId[m])
							{
								IdNode1 = NodesId[j];
								IdNode2 = NodesId[m]; 
								IdNode3 = NodesId[i]; 
							}
							else
							{
								IdNode1 = NodesId[m];
								IdNode2 = NodesId[j]; 
								IdNode3 = NodesId[i]; 
							}
						}
					}

					key = IdNode1*IdNode2*IdNode3;

					if ( mFaceHONodes.find(key) != mFaceHONodes.end() )
					{
						std::vector<std::vector<unsigned int> > ExistantNodesOnFaces(mFaceHONodes[key]);
						std::vector<std::vector<unsigned int> >::iterator it_ef;

						int repe = 0;

						for( it_ef = ExistantNodesOnFaces.begin(); it_ef != ExistantNodesOnFaces.end(); ++it_ef )
						{
							std::vector<unsigned int> NodesOnFace(*it_ef);
                
							if ((IdNode1==NodesOnFace[0]) && (IdNode2==NodesOnFace[1]) && (IdNode3==NodesOnFace[2])) repe++;
						}

						if ( repe == 0 )
						{
							if      ( mElementOrder == 3 ) CreateHONodes_OnFace_3th(IdNode1, IdNode2, IdNode3);
							else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
						}
					}
					else
					{
						if      ( mElementOrder == 3 ) CreateHONodes_OnFace_3th(IdNode1, IdNode2, IdNode3);
						else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
					}
				}
			}
		}	
	}

	//************************************************************************************************* 
    //* - Create High order nodes on volums
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnVolum(int* NodesId)
	{
		if ( mElementOrder == 4 ) CreateHONodes_OnVolum_4th(NodesId);	
	}

	//************************************************************************************************* 
    //* - Creating new higher order nodes on edge (2nd order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_2nd(int IdNode1, int IdNode2)
    {
		std::vector<unsigned int> edgeNodes(3);
        
        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back(edgeNodes);  

        Node::Pointer pNode1 = mpModel->GetNode(IdNode1-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNode2-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 0.5;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(edgeNodes[2], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);
	}

	//************************************************************************************************* 
    //* - Creating new higher order nodes on edge (3th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_3th(int IdNode1, int IdNode2)
    {
		std::vector<unsigned int> edgeNodes(4);
        
        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;
        edgeNodes[3] = edgeNodes[2]+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back(edgeNodes);  

        Node::Pointer pNode1 = mpModel->GetNode(IdNode1-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNode2-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 1.00/3.00;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(edgeNodes[2], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);

		//Creating new high order node on edge number 2
        alpha = 2.00/3.00;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2(new Node(edgeNodes[3], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node2);
        new_ho_node2->AddDofs(mDofSet);
	}

	//************************************************************************************************* 
    //* - Creating new higher order nodes on face (3th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnFace_3th(int IdNode1, int IdNode2, int IdNode3)
    {
		std::vector<unsigned int> faceNodes(4);
        
        faceNodes[0] = IdNode1;
        faceNodes[1] = IdNode2;
        faceNodes[2] = IdNode3;

        faceNodes[3] = (mpModel->GetNodesArray()).size()+1;

        mFaceHONodes[IdNode1*IdNode2*IdNode3].push_back(faceNodes);
		
		Node::Pointer pNode1 = mpModel->GetNode(IdNode1-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNode2-1);
        Node::Pointer pNode3 = mpModel->GetNode(IdNode3-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();

		double hoX, hoY, hoZ;
        double alpha, beta;
		
        //Creating new high order node on face number 1
        alpha = 1.00/3.00;
        beta  = 1.00/3.00;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(faceNodes[3], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);   
	}

	//************************************************************************************************* 
    //* - Creating new higher order nodes on edge (4th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_4th(int IdNode1, int IdNode2)
    {
        std::vector<unsigned int> edgeNodes(5);
        
        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;
        edgeNodes[3] = edgeNodes[2]+1;
        edgeNodes[4] = edgeNodes[3]+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back(edgeNodes);  

        Node::Pointer pNode1 = mpModel->GetNode(IdNode1-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNode2-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 0.25;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(edgeNodes[2], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);

		//Creating new high order node on edge number 2
        alpha = 0.50;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2(new Node(edgeNodes[3], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node2);
        new_ho_node2->AddDofs(mDofSet);

		//Creating new high order node on edge number 3
        alpha = 0.75;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node3(new Node(edgeNodes[4], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node3);
        new_ho_node3->AddDofs(mDofSet);    
    }

	//************************************************************************************************* 
    //* - Creating new higher order nodes on face (4th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnFace_4th(int IdNode1, int IdNode2, int IdNode3)
    {
        std::vector<unsigned int> faceNodes(6);
        
        faceNodes[0] = IdNode1;
        faceNodes[1] = IdNode2;
        faceNodes[2] = IdNode3;

        faceNodes[3] = (mpModel->GetNodesArray()).size()+1;
        faceNodes[4] = faceNodes[3]+1;
        faceNodes[5] = faceNodes[4]+1;

        mFaceHONodes[IdNode1*IdNode2*IdNode3].push_back(faceNodes);
		
		Node::Pointer pNode1 = mpModel->GetNode(IdNode1-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNode2-1);
        Node::Pointer pNode3 = mpModel->GetNode(IdNode3-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();

		double hoX, hoY, hoZ;
        double alpha, beta;
		
        //Creating new high order node on face number 1
        alpha = 0.25;
        beta  = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(faceNodes[3], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);   

		//Creating new high order node on face number 2
        alpha = 0.50;
        beta  = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node2(new Node(faceNodes[4], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node2);
        new_ho_node2->AddDofs(mDofSet);  
		
		//Creating new high order node on face number 3
        alpha = 0.25;
        beta  = 0.50;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node3(new Node(faceNodes[5], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node3);
        new_ho_node3->AddDofs(mDofSet);   
    }

	//************************************************************************************************* 
    //* - Creating new higher order nodes on volume (4th order)
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnVolum_4th(int* IdNodes)
    {
		std::vector<unsigned int> volumNodes(5);
        
        volumNodes[0] = IdNodes[0];
        volumNodes[1] = IdNodes[1];
        volumNodes[2] = IdNodes[2];
		volumNodes[3] = IdNodes[3];
        
        volumNodes[4] = (mpModel->GetNodesArray()).size()+1;
        
        mVolmHONodes[IdNodes[0]*IdNodes[1]*IdNodes[2]*IdNodes[3]].push_back(volumNodes);
		
		Node::Pointer pNode1 = mpModel->GetNode(IdNodes[0]-1);
        Node::Pointer pNode2 = mpModel->GetNode(IdNodes[1]-1);
        Node::Pointer pNode3 = mpModel->GetNode(IdNodes[2]-1);
		Node::Pointer pNode4 = mpModel->GetNode(IdNodes[3]-1);

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();
		double X4 = pNode4->X(), Y4 = pNode4->Y(), Z4 = pNode4->Z();

		double hoX, hoY, hoZ;
        double alpha, beta, gamma;
		
        //Creating new high order node on volume number 1
        alpha = 0.25;
        beta  = 0.25;
		gamma = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + gamma * ( X4 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + gamma * ( Y4 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + gamma * ( Z4 - Z1 ) + Z1;

        Node::Pointer new_ho_node1(new Node(volumNodes[4], hoX, hoY, hoZ));
        mpModel->AddNode(new_ho_node1);
        new_ho_node1->AddDofs(mDofSet);   
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - Building
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//************************************************************************************************* 
    //* - Assembling high order volume elements
    //*************************************************************************************************
    void Modeler::GenerateVolumeElement(int* NodesId, unsigned int PropertiesId)
    { 
		// If PropertyId == 0 throw error, material not assigned
		if (PropertiesId == 0) { throw 0; }
		
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));
    
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new VolumeElement_2ndQ(nodes,properties));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new VolumeElement_1st(nodes,properties));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new VolumeElement_2nd(nodes,properties));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new VolumeElement_3th(nodes,properties));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new VolumeElement_4th(nodes,properties));
		}

		bool IsUgLayer = false;

		for(int i=0; i<4; i++)
		{
			if ( ( mSingular.find   (NodesId[i]) != mSingular.end   () ) || ( mNSingular2L.find(NodesId[i]) != mNSingular2L.end() ) ||
				 ( mNSingular3L.find(NodesId[i]) != mNSingular3L.end() ) || ( mNSingular4L.find(NodesId[i]) != mNSingular4L.end() ) ||
				 ( mNSingular5L.find(NodesId[i]) != mNSingular5L.end() ) || ( mNSingular6L.find(NodesId[i]) != mNSingular6L.end() )  )
			{
				IsUgLayer = true;
			}
		}

		//pElement->SetPeso(0.00);

		if ( IsUgLayer == true ) 
		{
			//CalculateDistToSing( NodesId, nodes );
			pElement->SetPeso(0.00);	
		}
		else
		{
			pElement->SetPeso(1.00);	
		}

        Vector<int> EleIdVector; 
        Matrix<std::complex<double> > EleStiffMatrix;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);

		//Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
	    PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);

		//Trans(T)*A*T for contact discontinuities
		ContactHO_TKT(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for PEC, PMC conditions 
        HO_TKT(HONodesId, EleStiffMatrix);

		//Trans(T)*A*T for axisymmetric problems
		Axisym_TKT(HONodesId, EleStiffMatrix);

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for (int i_local = 0; i_local < eleSize; i_local++)
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize ) 
			{
				for (int j_local = 0; j_local < eleSize; j_local++)
				{
					int j_global = EleIdVector[j_local];

					if ( ( i_global <= j_global ) && (j_global < mSystemSize) ) 
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( cKij != czero ) A_matrix(i_global,j_global) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//************************************************************************************************* 
    //* - Calculates the maximum and the minimum distance to the singularity
    //*************************************************************************************************
	void Modeler::CalculateDistToSing( int* NodesId, std::vector<Node::Pointer> nodes )
	{
		double X1 = nodes[0]->X(), Y1 = nodes[0]->Y(), Z1 = nodes[0]->Z();
		double X2 = nodes[1]->X(), Y2 = nodes[1]->Y(), Z2 = nodes[1]->Z();
		double X3 = nodes[2]->X(), Y3 = nodes[2]->Y(), Z3 = nodes[2]->Z();
		double X4 = nodes[3]->X(), Y4 = nodes[3]->Y(), Z4 = nodes[3]->Z();

		double distance; 

		//distance [0 1]
		distance = sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2) );

		if ( mMinDistToSing  == 0.00 )
		{ 
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[0];
			mMinDistNode2  = NodesId[1];
		}

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[0];
			mMaxDistNode2  = NodesId[1];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[0];
			mMinDistNode2  = NodesId[1];
		}

		//distance [0 2]
		distance = sqrt( (X1-X3)*(X1-X3) + (Y1-Y3)*(Y1-Y3) + (Z1-Z3)*(Z1-Z3) );

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[0];
			mMaxDistNode2  = NodesId[2];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[0];
			mMinDistNode2  = NodesId[2];
		}

		//distance [0 3]
		distance = sqrt( (X1-X4)*(X1-X4) + (Y1-Y4)*(Y1-Y4) + (Z1-Z4)*(Z1-Z4) );

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[0];
			mMaxDistNode2  = NodesId[3];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[0];
			mMinDistNode2  = NodesId[3];
		}

		//distance [1 2]
		distance = sqrt( (X2-X3)*(X2-X3) + (Y2-Y3)*(Y2-Y3) + (Z2-Z3)*(Z2-Z3) );

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[1];
			mMaxDistNode2  = NodesId[2];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[1];
			mMinDistNode2  = NodesId[2];
		}

		//distance [1 3]
		distance = sqrt( (X2-X4)*(X2-X4) + (Y2-Y4)*(Y2-Y4) + (Z2-Z4)*(Z2-Z4) );

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[1];
			mMaxDistNode2  = NodesId[3];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[1];
			mMinDistNode2  = NodesId[3];
		}

		//distance [2 3]
		distance = sqrt( (X3-X4)*(X3-X4) + (Y3-Y4)*(Y3-Y4) + (Z3-Z4)*(Z3-Z4) );

		if ( distance > mMaxDistToSing )
		{
			mMaxDistToSing = distance;
			mMaxDistNode1  = NodesId[2];
			mMaxDistNode2  = NodesId[3];
		}
		else if ( distance < mMinDistToSing )
		{
			mMinDistToSing = distance;
			mMinDistNode1  = NodesId[2];
			mMinDistNode2  = NodesId[3];
		}	 
	}

	//************************************************************************************************* 
    //* - Assembling high order source elements
    //*************************************************************************************************
    void Modeler::GenerateSourceElement(int* NodesId, unsigned int PropertiesId)
    { 
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) nodes.push_back(mpModel->GetNode((*hoit)-1));

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new JSource_2ndQ(nodes,properties));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new JSource_1st(nodes,properties));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new JSource_2nd(nodes,properties));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new JSource_3th(nodes,properties));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new JSource_4th(nodes,properties));
		}
        
        Vector<int> EleIdVector; 
        Vector<std::complex<double> > EleResVector;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetResidualVector(EleResVector);

		//Inv(T)*b for PBC condition
		PBC_FrontBack_TR(HONodesId, EleIdVector, EleResVector);
	    PBC_RightLeft_TR(HONodesId, EleIdVector, EleResVector);
		PBC_FrontBack_TR(HONodesId, EleIdVector, EleResVector);

		//Inv(T)*b in dielectric contact surfaces
		ContactHO_TR(HONodesId, EleIdVector, EleResVector);

        //Inv(T)*b
        HO_TR(HONodesId,EleResVector);

		//Inv(T)*b for axisymmetric problems
		Axisym_TR(HONodesId,EleResVector);

        //Assembling to global source vector
        int eleSize = EleIdVector.size();

		for (int i_local=0; i_local < eleSize; i_local++)
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize ) 
			{	
				b_vector[i_global] += EleResVector[i_local];
			}
		}

        EleIdVector.FreeData(); 
        EleResVector.FreeData(); 
    }

	//************************************************************************************************* 
    //* - Assembling far field elements
    //*************************************************************************************************
    void Modeler::GenerateFarFieldElement(int* NodesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));
        
        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new FarField_2ndQ(nodes));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new FarField_1st(nodes));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new FarField_2nd(nodes));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new FarField_3th(nodes));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new FarField_4th(nodes));
		}

		pElement->SetFrequency(mProblemFrequency);

		Vector<int> EleIdVector; 
        Matrix<std::complex<double> > EleStiffMatrix;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);
		
        //Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
	    PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		
		//Trans(T)*A*T
        HO_TKT(HONodesId, EleStiffMatrix);

		//Trans(T)*A*T for axisymmetric problems
		Axisym_TKT(HONodesId,EleStiffMatrix);

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for (int i_local = 0; i_local < eleSize; i_local++)
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize ) 
			{
				for (int j_local = 0; j_local < eleSize; j_local++)
				{
					int j_global = EleIdVector[j_local];

					if ( ( i_global <= j_global ) && (j_global < mSystemSize) ) 
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( cKij != czero ) A_matrix(i_global,j_global) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
	}

	//************************************************************************************************* 
    //* - Assembling generic Robin elements
    //*************************************************************************************************
    void Modeler::GenerateGenericRobinElement(int* NodesId, unsigned int PropertiesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);
        
        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new GenericRobin_2ndQ(nodes,properties));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new GenericRobin_1st(nodes,properties));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new GenericRobin_2nd(nodes,properties));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new GenericRobin_3th(nodes,properties));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new GenericRobin_4th(nodes,properties));
		}

		Vector<int> EleIdVector; 
        Matrix<std::complex<double> > EleStiffMatrix;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);

		//Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
	    PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		
        //Trans(T)*A*T
        HO_TKT(HONodesId, EleStiffMatrix);

		//Trans(T)*A*T for axisymmetric problems
		Axisym_TKT(HONodesId,EleStiffMatrix);

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for (int i_local = 0; i_local < eleSize; i_local++)
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize ) 
			{
				for (int j_local = 0; j_local < eleSize; j_local++)
				{
					int j_global = EleIdVector[j_local];

					if ( ( i_global <= j_global ) && (j_global < mSystemSize) ) 
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( cKij != czero ) A_matrix(i_global,j_global) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
	}

	//************************************************************************************************* 
    //* - Assembling boundary elements
    //*************************************************************************************************
    void Modeler::GenerateRWPortTE10Element(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if      ( mElementOrder == 1 ) pElement = Element::Pointer(new RWPortTE10_1st(nodes,properties));
        else if ( mElementOrder == 2 ) pElement = Element::Pointer(new RWPortTE10_2nd(nodes,properties));
        else if ( mElementOrder == 3 ) pElement = Element::Pointer(new RWPortTE10_3th(nodes,properties));
		else if ( mElementOrder == 4 ) pElement = Element::Pointer(new RWPortTE10_4th(nodes,properties));

		Vector<int> EleIdVector; 
		Vector<std::complex<double> > EleResVector;
        Matrix<std::complex<double> > EleStiffMatrix;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);
        pElement->GetResidualVector(EleResVector);
   
        //Trans(T)*A*T
        HO_TKT(HONodesId,EleStiffMatrix);

		//Inv(T)*b
        HO_TR(HONodesId,EleResVector);

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for (int i_Klocal = 0; i_Klocal < eleSize; i_Klocal++)
		{
			int i_Kglobal = EleIdVector[i_Klocal];

			if ( i_Kglobal < mSystemSize ) 
			{
				for (int j_Klocal = 0; j_Klocal < eleSize; j_Klocal++)
				{
					int j_Kglobal = EleIdVector[j_Klocal];

					if ( ( i_Kglobal <= j_Kglobal ) && (j_Kglobal < mSystemSize) ) 
                    {
                        cKij = EleStiffMatrix[i_Klocal][j_Klocal];
                        if ( cKij != czero ) A_matrix(i_Kglobal,j_Kglobal) += cKij;
                    }
				}
			}
        }

        //Assembling to global source vector
		for (int i_blocal=0; i_blocal < eleSize; i_blocal++)
		{
			int i_bglobal = EleIdVector[i_blocal];

			if ( i_bglobal < mSystemSize ) 
			{	
				b_vector[i_bglobal] += EleResVector[i_blocal];
			}
		}

        EleIdVector.FreeData(); 
        EleResVector.FreeData(); 
        EleStiffMatrix.FreeData();
    }

	//************************************************************************************************* 
    //* - Assembling boundary elements
    //*************************************************************************************************
    void Modeler::GenerateCoaxPortTEMElement(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new CoaxialPortTEM_2ndQ(nodes,properties));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new CoaxialPortTEM_1st(nodes,properties));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new CoaxialPortTEM_2nd(nodes,properties));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new CoaxialPortTEM_3th(nodes,properties));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new CoaxialPortTEM_4th(nodes,properties));
		}

		Vector<int> EleIdVector; 
		Vector<std::complex<double> > EleResVector;
        Matrix<std::complex<double> > EleStiffMatrix;
        
        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);
        pElement->GetResidualVector(EleResVector);
   
        //Trans(T)*A*T
        HO_TKT(HONodesId,EleStiffMatrix);

		//Inv(T)*b
        HO_TR(HONodesId,EleResVector);

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for (int i_Klocal = 0; i_Klocal < eleSize; i_Klocal++)
		{
			int i_Kglobal = EleIdVector[i_Klocal];

			if ( i_Kglobal < mSystemSize ) 
			{
				for (int j_Klocal = 0; j_Klocal < eleSize; j_Klocal++)
				{
					int j_Kglobal = EleIdVector[j_Klocal];

					if ( ( i_Kglobal <= j_Kglobal ) && (j_Kglobal < mSystemSize) ) 
                    {
                        cKij = EleStiffMatrix[i_Klocal][j_Klocal];
                        if ( cKij != czero ) A_matrix(i_Kglobal,j_Kglobal) += cKij;
                    }
				}
			}
        }

        //Assembling to global source vector
		for (int i_blocal=0; i_blocal < eleSize; i_blocal++)
		{
			int i_bglobal = EleIdVector[i_blocal];

			if ( i_bglobal < mSystemSize ) 
			{	
				b_vector[i_bglobal] += EleResVector[i_blocal];
			}
		}

        EleIdVector.FreeData(); 
        EleResVector.FreeData(); 
        EleStiffMatrix.FreeData();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - Matrix rotations
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//************************************************************************************************* 
    //* - Inv(T)*StiffnesMatrix*T, high order 
    //*************************************************************************************************
    void Modeler::HO_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix)
    {
		if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;
		
		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find(*hoit) != mNormals.end() )
            {
                //////////////  Local coordinate system  //////////////////
                std::vector<double> n(mNormals[*hoit]);
                std::vector<double> t(3), b(3);
                TangencialCoordinates(n,t,b);

                //////////////  Inv(T)*StiffnesMatrix*T  //////////////////
                int i;

                int step    = HONodesId.size();
				int step2   = 2*step;
                int dofSize = 3*step;

                std::complex<double> dof_x, dof_y, dof_z;

                //StiffnesMatrix*T
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[i][local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[i][local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

                //Inv(T)*StiffnesMatrix
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

   					StiffMatrix[local_i        ][i] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[local_i + step ][i] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[local_i + step2][i] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

            }//end if( mNormals.find(*hoit) != mNormals.end() )
			
            local_i++; 

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 
    }

	//************************************************************************************************* 
    //* - Inv(T)*StiffnesMatrix*T, for axisymmetric problems
    //*************************************************************************************************
    void Modeler::Axisym_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix)
    {
        if ( !mAxisymmetric ) return;
		
		int local_i = 0;

		double distToAxis;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
			if ( ( mNormals.find(*hoit) == mNormals.end() ) && ( DistanceToAxis(*hoit) > 1e-12 ) )
            {
				//////////////  Local coordinate system  //////////////////
				std::vector<double> n(3), t(3), b(3);

				AxisymmetricCoordinates(*hoit,n,t,b);

                //////////////  Inv(T)*StiffnesMatrix*T  //////////////////
                int i;

                int step    = HONodesId.size();
				int step2   = 2*step;
                int dofSize = 3*step;

                std::complex<double> dof_x, dof_y, dof_z;

                //StiffnesMatrix*T
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[i][local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[i][local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

                //Inv(T)*StiffnesMatrix
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

   					StiffMatrix[local_i        ][i] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[local_i + step ][i] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[local_i + step2][i] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

            }//end if( mNormals.find(*hoit) != mNormals.end() )
			
            local_i++; 

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 
    }

	//************************************************************************************************* 
    //* - Inv(T)*b
    //*************************************************************************************************
    void Modeler::HO_TR(std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector)
    {
        if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;
		
		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find(*hoit) != mNormals.end() )
            {
                //////////////  Local coordinate system  //////////////////
                std::vector<double> n(mNormals[*hoit]);
                std::vector<double> t(3), b(3);
                TangencialCoordinates(n,t,b);

                ///////////////////  Inv(T)*b  ////////////////////
                int step  = HONodesId.size();
				int step2 = 2*step;
              
                std::complex<double> dof_x, dof_y, dof_z;
               
				dof_x = ResVector[local_i        ];
                dof_y = ResVector[local_i + step ];
                dof_z = ResVector[local_i + step2];

                ResVector[local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                ResVector[local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                ResVector[local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
            }

            local_i++; 
        }   
    }

	//************************************************************************************************* 
    //* - Inv(T)*b, for axisymmetric problems.
    //*************************************************************************************************
    void Modeler::Axisym_TR(std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector)
    {
        if ( !mAxisymmetric ) return;
		
		int local_i = 0;

        double distToAxis;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
			if ( ( mNormals.find(*hoit) == mNormals.end() ) && ( DistanceToAxis(*hoit) > 1e-12 ) )
            {
				//////////////  Local coordinate system  //////////////////
				std::vector<double> n(3), t(3), b(3);

				AxisymmetricCoordinates(*hoit,n,t,b);

                ///////////////////  Inv(T)*b  ////////////////////
                int step  = HONodesId.size();
				int step2 = 2*step;
              
                std::complex<double> dof_x, dof_y, dof_z;
               
				dof_x = ResVector[local_i        ];
                dof_y = ResVector[local_i + step ];
                dof_z = ResVector[local_i + step2];

                ResVector[local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                ResVector[local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                ResVector[local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
            }

            local_i++; 
        }   
    }

	//************************************************************************************************* 
    //* -  Local system of coordinates in PEC surface for n 
    //*************************************************************************************************
    void Modeler::TangencialCoordinates(std::vector<double>& n, 
                                        std::vector<double>& t, 
                                        std::vector<double>& b)
    {
         double tnorm = sqrt(n[0]*n[0]+n[1]*n[1]);

         if (tnorm > 0.00)
         {
             t[0] = n[1]/tnorm; 
             t[1] =-n[0]/tnorm; 
             t[2] = 0.00;
         }
         else
         {
             t[0]  = n[2]; 
             t[1]  = 0.00; 
             t[2]  = 0.00;
         }

         b[0] = n[1]*t[2] - n[2]*t[1];
         b[1] = n[2]*t[0] - n[0]*t[2];
         b[2] = n[0]*t[1] - n[1]*t[0];
    }

	//************************************************************************************************* 
    //* - Distance to axisymmetric axis.
    //*************************************************************************************************
	double Modeler::DistanceToAxis(unsigned int NodeId)
	{
		Node::Pointer pNode = mpModel->GetNode(NodeId-1);

		double X_coord = pNode->X();
		double Z_coord = pNode->Z();

		double dist = sqrt( X_coord*X_coord + Z_coord*Z_coord );

		return dist;
	}

	//************************************************************************************************* 
    //* -  Axisymmetric unitary vectors
    //*************************************************************************************************
    void Modeler::AxisymmetricCoordinates(unsigned int NodeId,
		                                  std::vector<double>& n, 
                                          std::vector<double>& t, 
                                          std::vector<double>& b)
    {   
		Node::Pointer pNode = mpModel->GetNode(NodeId-1);

		t[0] = pNode->X();
		t[1] = 0.00;
		t[2] = pNode->Z();

		double tnorm = sqrt( t[0]*t[0] + t[2]*t[2] );

		if (tnorm > 0.00)
		{
			t[0] = t[0]/tnorm; 
			t[2] = t[2]/tnorm; 
		}

		n[0] =-t[2]; 
		n[1] = 0.00; 
		n[2] = t[0];

		b[0] = 0.00;
		b[1] = 1.00;
		b[2] = 0.00;
    }

    //************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC nodes
    //*************************************************************************************************
    void Modeler::PBC_FrontBack_TKT(std::vector<int>& HONodesId, 
		                            Vector<int>& IdVector,
							        Matrix<std::complex<double> >& StiffMatrix)
   {   
	    if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 ) return;
	    
		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );
		
		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Front.find(*ohoit) != mPBC_NodeNodePairs_Front.end() )
            {
			    PBC_FB_NodeToNode_cID( *ohoit, local_i, HONodesId, IdVector ); 
            }
			else if ( mPBC_NodeEdgePairs_Front.find(*ohoit) != mPBC_NodeEdgePairs_Front.end() )
			{
			    PBC_FB_NodeToEdge_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 			
			}
			else if ( mPBC_NodeElementPairs_Front.find(*ohoit) != mPBC_NodeElementPairs_Front.end() )
			{
			    PBC_FB_NodeToElement_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 					
			}
            local_i++; 
        }
    }

	//************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC node to node
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToNode_cID(int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector)
	{
		int backNode = mPBC_NodeNodePairs_Front[NodeId];
		int step     = IdVector.size()/3;
	
		Node::Pointer pbNode = mpModel->GetNode( backNode - 1 );

        HONodesId[local_i         ] = pbNode->Id();

		IdVector [local_i         ] = pbNode->pDofcEx()->EquationId();
		IdVector [local_i +   step] = pbNode->pDofcEy()->EquationId();
		IdVector [local_i + 2*step] = pbNode->pDofcEz()->EquationId();        
	}

	//************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC node to edge
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToEdge_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<unsigned int> bNodes( mPBC_NodeEdgePairs_Front[NodeId] );	

		std::vector<Node::Pointer> pbNodes; 
		
		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		double L2 = mPBC_NodeEdgeCoord_Front[NodeId];
        double L1 = 1.00 - L2;

		std::vector<double> N(3);

		N[0] = L1 * (2.00*L1 - 1.00);
		N[1] = L2 * (2.00*L2 - 1.00);
		N[2] = 4.00 * L1 * L2;

		PBC_StiffMatrixReDef(local_i, oldSize, newSize, N, StiffMatrix);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Redefinition of IdVector for PBC boundary conditions
    //*************************************************************************************************
	void Modeler::PBC_IdVectorReDef(int local_i, std::vector<Node::Pointer>& pNodes, std::vector<int>& HONodesId, Vector<int>& IdVector)
	{
		std::vector<int> newHONodesId( HONodesId.size() + ( ( pNodes.size() - 1 )     ) );

		Vector<int>      newIdVector ( IdVector.size()  + ( ( pNodes.size() - 1 ) * 3 ) );

		int step     = IdVector.size()/3;
		int step2    = 2*step;

		int newStep  = newIdVector.size()/3;
		int newStep2 = 2*newStep;

		int i, j;

		for ( i = 0; i < local_i; i++ ) 
		{ 
            newHONodesId[i           ] = HONodesId[i        ];

			newIdVector [i           ] = IdVector [i        ]; 
			newIdVector [i + newStep ] = IdVector [i + step ]; 
			newIdVector [i + newStep2] = IdVector [i + step2]; 
		}
		for ( i = local_i, j = 0; j < pNodes.size(); i++, j++ ) 
		{ 
			newHONodesId[i           ] = pNodes[j]->Id();
			
			newIdVector [i           ] = pNodes[j]->pDofcEx()->EquationId(); 
			newIdVector [i + newStep ] = pNodes[j]->pDofcEy()->EquationId(); 
			newIdVector [i + newStep2] = pNodes[j]->pDofcEz()->EquationId();
		}
		for ( i = local_i + pNodes.size(), j = local_i + 1; j < step; i++, j++ ) 
		{ 
			newHONodesId[i           ] = HONodesId[j        ];
			
			newIdVector [i           ] = IdVector [j        ]; 
			newIdVector [i + newStep ] = IdVector [j + step ]; 
			newIdVector [i + newStep2] = IdVector [j + step2]; 
		}

		HONodesId.resize( newHONodesId.size() ); 
		
		for ( i = 0; i < newHONodesId.size(); i++ ) { HONodesId[i] = newHONodesId[i]; }
		
		IdVector.resize( newIdVector.size() ); 
		
		for ( i = 0; i < newIdVector.size() ; i++ ) { IdVector [i] = newIdVector [i]; }
	}

	//************************************************************************************************* 
    //* - Inv(T)*StiffnesMatrix*T, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<unsigned int> bNodes( mPBC_NodeElementPairs_Front[NodeId] );	

		std::vector<Node::Pointer> pbNodes; 
		
		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		std::vector<double> natCoord(mPBC_NodeElementCoord_Front[NodeId]);
		
		double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;
		
		std::vector<double> N(6);

        N[0] = L1 * (2.00*L1 - 1.00);
        N[1] = L2 * (2.00*L2 - 1.00);
        N[2] = L3 * (2.00*L3 - 1.00);
		N[3] = 4.00 * L1 * L2;
		N[4] = 4.00 * L2 * L3;
        N[5] = 4.00 * L3 * L1;

		PBC_StiffMatrixReDef(local_i, oldSize, newSize, N, StiffMatrix);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Redefinition of the stiffness matrix for PBC boundary conditions
    //*************************************************************************************************
	void Modeler::PBC_StiffMatrixReDef(int local_i, int oldSize, int newSize,
		                               std::vector<double>& N, 
									   Matrix<std::complex<double> >& StiffMatrix)
	{
	    int inew, jnew;
		int iold, jold;
		int ib  , ibf ;

		int oldBlockSize = oldSize/3;
		int newBlockSize = newSize/3;

		Matrix<std::complex<double> > newStiffMatrix( newSize, newSize );

		for ( iold = 0, inew = 0; iold < oldSize; iold++ )
		{
			if ( ( iold != local_i ) && ( iold != local_i + oldBlockSize ) && ( iold != local_i + 2*oldBlockSize ) )
			{
				for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
				{
					if ( jold != local_i )
					{
						newStiffMatrix[inew][jnew                 ] = StiffMatrix[iold][jold                 ]; 
						newStiffMatrix[inew][jnew +   newBlockSize] = StiffMatrix[iold][jold +   oldBlockSize]; 
						newStiffMatrix[inew][jnew + 2*newBlockSize] = StiffMatrix[iold][jold + 2*oldBlockSize]; 
						
						jnew++;
						continue;
					}
						
					for ( ib = 0; ib < N.size(); ib++ ) 
					{
						newStiffMatrix[inew][jnew +                  ib] = N[ib] * StiffMatrix[iold][jold                 ];
						newStiffMatrix[inew][jnew +   newBlockSize + ib] = N[ib] * StiffMatrix[iold][jold +   oldBlockSize];
						newStiffMatrix[inew][jnew + 2*newBlockSize + ib] = N[ib] * StiffMatrix[iold][jold + 2*oldBlockSize]; 
					}

					jnew = jnew + N.size();
				}

				inew++;
				continue;
			}
			
			for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
			{
				if ( jold != local_i )
				{
					for ( ibf = 0; ibf < N.size(); ibf++ ) 
					{
						newStiffMatrix[inew + ibf][jnew                 ] = N[ibf] * StiffMatrix[iold][jold                 ]; 
						newStiffMatrix[inew + ibf][jnew +   newBlockSize] = N[ibf] * StiffMatrix[iold][jold +   oldBlockSize];
						newStiffMatrix[inew + ibf][jnew + 2*newBlockSize] = N[ibf] * StiffMatrix[iold][jold + 2*oldBlockSize]; 		
					}

					jnew++;
					continue;
				}

				for ( ib = 0; ib < N.size(); ib++ ) 
				{
					for ( ibf = 0; ibf < N.size(); ibf++ ) 
					{
						newStiffMatrix[inew + ibf][jnew +                  ib] = N[ibf] * N[ib] * StiffMatrix[iold][jold                 ];
						newStiffMatrix[inew + ibf][jnew +   newBlockSize + ib] = N[ibf] * N[ib] * StiffMatrix[iold][jold +   oldBlockSize];
						newStiffMatrix[inew + ibf][jnew + 2*newBlockSize + ib] = N[ibf] * N[ib] * StiffMatrix[iold][jold + 2*oldBlockSize]; 
					}
				}

				jnew = jnew + N.size();
			}

			inew = inew + N.size();
		}

		StiffMatrix.Resize( newSize, newSize );
		
		for (int i = 0; i < newSize; i++ )
		{
			for (int j = 0; j < newSize; j++ )
			{
		    	StiffMatrix[i][j] = newStiffMatrix[i][j];
			}	
		}
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC nodes
    //*************************************************************************************************
	void Modeler::PBC_FrontBack_TR(std::vector<int>& HONodesId, 
		                           Vector<int>& IdVector, 
							       Vector<std::complex<double> >& ResVector)
	{
        if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 ) return;
		
		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );
		
		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Front.find(*ohoit) != mPBC_NodeNodePairs_Front.end() )
            {
			    PBC_FB_NodeToNode_cID( *ohoit, local_i, HONodesId, IdVector ); 
            }
			else if ( mPBC_NodeEdgePairs_Front.find(*ohoit) != mPBC_NodeEdgePairs_Front.end() )
			{
			    PBC_FB_NodeToEdge_TR( *ohoit, local_i, HONodesId, IdVector, ResVector ); 			
			}
			else if ( mPBC_NodeElementPairs_Front.find(*ohoit) != mPBC_NodeElementPairs_Front.end() )
			{
			    PBC_FB_NodeToElement_TR( *ohoit, local_i, HONodesId, IdVector, ResVector ); 					
			}
            local_i++; 
        }
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to edge
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToEdge_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<unsigned int> bNodes( mPBC_NodeEdgePairs_Front[NodeId] );	

		std::vector<Node::Pointer> pbNodes; 
		
		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		double L2 = mPBC_NodeEdgeCoord_Front[NodeId];
        double L1 = 1.00 - L2;

		std::vector<double> N(3);

		N[0] = L1 * (2.00*L1 - 1.00);
		N[1] = L2 * (2.00*L2 - 1.00);
		N[2] = 4.00 * L1 * L2;

		PBC_ResVectorReDef(local_i, oldSize, newSize, N, ResVector);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<unsigned int> bNodes( mPBC_NodeElementPairs_Front[NodeId] );	

		std::vector<Node::Pointer> pbNodes; 
		
		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		std::vector<double> natCoord(mPBC_NodeElementCoord_Front[NodeId]);
		
		double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;
		
		std::vector<double> N(6);

        N[0] = L1 * (2.00*L1 - 1.00);
        N[1] = L2 * (2.00*L2 - 1.00);
        N[2] = L3 * (2.00*L3 - 1.00);
		N[3] = 4.00 * L1 * L2;
		N[4] = 4.00 * L2 * L3;
        N[5] = 4.00 * L3 * L1;

		PBC_ResVectorReDef(local_i, oldSize, newSize, N, ResVector);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Redefinition of the residual vector for PBC boundary conditions
    //*************************************************************************************************
	void Modeler::PBC_ResVectorReDef(int local_i, int oldSize, int newSize, std::vector<double>& N, Vector<std::complex<double> >& ResVector)
	{
		int inew, iold, ib;

		int oldBlockSize = oldSize/3;
		int newBlockSize = newSize/3;

		Vector<std::complex<double> > newResVector( newSize );

		for ( iold = 0, inew = 0; iold < oldBlockSize; iold++ )
		{
			if ( iold != local_i )
			{
				newResVector[inew                 ] = ResVector[iold                 ]; 
			    newResVector[inew +   newBlockSize] = ResVector[iold +   oldBlockSize]; 
			    newResVector[inew + 2*newBlockSize] = ResVector[iold + 2*oldBlockSize]; 
				
				inew++;
				continue;
			}

			for ( ib = 0; ib < N.size(); ib++ ) 
			{
				newResVector[inew + ib                 ] = N[ib] * ResVector[iold                 ]; 
				newResVector[inew + ib +   newBlockSize] = N[ib] * ResVector[iold +   oldBlockSize];
				newResVector[inew + ib + 2*newBlockSize] = N[ib] * ResVector[iold + 2*oldBlockSize]; 		
			}

			inew = inew + N.size();
		}

		ResVector.resize( newSize ); 
		
		for (int i = 0; i < newSize; i++ )
		{
			ResVector[i] = newResVector[i];	
		}
	}

	//************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC nodes
    //*************************************************************************************************
    void Modeler::PBC_RightLeft_TKT(std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
   {
	    if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 ) return;
	    
		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );
		
		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if      ( mPBC_NodeNodePairs_Right.find(*ohoit) != mPBC_NodeNodePairs_Right.end() )
            {
				PBC_RL_NodeToNode_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 
            }
			else if ( mPBC_NodeEdgePairs_Right.find(*ohoit) != mPBC_NodeEdgePairs_Right.end() )
			{
				PBC_RL_NodeToEdge_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 
			}
			else if ( mPBC_NodeElementPairs_Right.find(*ohoit) != mPBC_NodeElementPairs_Right.end() )
			{
				PBC_RL_NodeToElement_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 
			}
            local_i++; 
        }
    }

	//************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC node to node
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToNode_TKT(int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		int lNode = mPBC_NodeNodePairs_Right[NodeId];
		int size  = IdVector.size();
		int step  = size/3;
	
		Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

		HONodesId[local_i         ] = plNode->Id();

		IdVector [local_i         ] = plNode->pDofcEx()->EquationId();
		IdVector [local_i +   step] = plNode->pDofcEy()->EquationId();
		IdVector [local_i + 2*step] = plNode->pDofcEz()->EquationId();   

		std::vector<double> N(1); 
		
		N[0] = 1.00;
		
		PBC_StiffMatrixReDef_Rot(local_i, size, size, N, NodeId, StiffMatrix);
	}

    //************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC node to edge
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToEdge_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<unsigned int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );	

		std::vector<Node::Pointer> plNodes; 
		
		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		double L2 = mPBC_NodeEdgeCoord_Right[NodeId];
        double L1 = 1.00 - L2;

		std::vector<double> N(3);

		N[0] = L1 * (2.00*L1 - 1.00);
		N[1] = L2 * (2.00*L2 - 1.00);
		N[2] = 4.00 * L1 * L2;

		PBC_StiffMatrixReDef_Rot(local_i, oldSize, newSize, N, NodeId, StiffMatrix);

		local_i = local_i + ( lNodes.size() - 1 );
	}

    //************************************************************************************************* 
    //* - Trans(T)*StiffnesMatrix*T, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<unsigned int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );	

		std::vector<Node::Pointer> plNodes; 
		
		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		std::vector<double> natCoord(mPBC_NodeElementCoord_Right[NodeId]);
		
		double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;
		
		std::vector<double> N(6);

        N[0] = L1 * (2.00*L1 - 1.00);
        N[1] = L2 * (2.00*L2 - 1.00);
        N[2] = L3 * (2.00*L3 - 1.00);
		N[3] = 4.00 * L1 * L2;
		N[4] = 4.00 * L2 * L3;
        N[5] = 4.00 * L3 * L1;

		PBC_StiffMatrixReDef_Rot(local_i, oldSize, newSize, N, NodeId, StiffMatrix);

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Redefinition of the stiffness matrix for PBC boundary conditions
    //*************************************************************************************************
	void Modeler::PBC_StiffMatrixReDef_Rot(int local_i, int oldSize, int newSize,
		                                   std::vector<double>& N, int NodeId,
									       Matrix<std::complex<double> >& StiffMatrix)
	{      
		double cos_A = 1.0;
		double sin_A = 0.0;

		if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

		int inew, jnew;
		int iold, jold;
		int ib  , ibf ;

		int oldBlockSize = oldSize / 3;
		int newBlockSize = newSize / 3;

		Matrix<std::complex<double> > newStiffMatrix( newSize, newSize );

		for ( iold = 0, inew = 0; iold < oldSize; iold++ )
		{
			if ( iold == local_i )
			{
			    for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
				{
					if ( jold != local_i )
					{   
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew                 ] = N[ibf] * ( cos_A * StiffMatrix[iold][jold                 ] - sin_A * StiffMatrix[iold + oldBlockSize][jold                 ] );
							newStiffMatrix[inew + ibf][jnew +   newBlockSize] = N[ibf] * ( cos_A * StiffMatrix[iold][jold +   oldBlockSize] - sin_A * StiffMatrix[iold + oldBlockSize][jold +   oldBlockSize] );
							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize] = N[ibf] * ( cos_A * StiffMatrix[iold][jold + 2*oldBlockSize] - sin_A * StiffMatrix[iold + oldBlockSize][jold + 2*oldBlockSize] );
						}

						jnew++;
						continue;
					}

					for ( ib = 0; ib < N.size(); ib++ ) 
					{
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew +                  ib] = N[ibf] * N[ib] * cos_A * ( cos_A * StiffMatrix[iold               ][jold] - sin_A * StiffMatrix[iold               ][jold + oldBlockSize] ) -
								                                                     N[ibf] * N[ib] * sin_A * ( cos_A * StiffMatrix[iold + oldBlockSize][jold] - sin_A * StiffMatrix[iold + oldBlockSize][jold + oldBlockSize] );

							newStiffMatrix[inew + ibf][jnew +   newBlockSize + ib] = N[ibf] * N[ib] * cos_A * ( sin_A * StiffMatrix[iold               ][jold] + cos_A * StiffMatrix[iold               ][jold + oldBlockSize] ) -
								                                                     N[ibf] * N[ib] * sin_A * ( sin_A * StiffMatrix[iold + oldBlockSize][jold] + cos_A * StiffMatrix[iold + oldBlockSize][jold + oldBlockSize] );

							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize + ib] = N[ibf] * N[ib] * ( cos_A * StiffMatrix[iold][jold + 2*oldBlockSize] - sin_A * StiffMatrix[iold + oldBlockSize][jold + 2*oldBlockSize] ); 
						}
					}

					jnew = jnew + N.size();
				}

				inew = inew + N.size();
			}
			else if ( iold == (local_i + oldBlockSize) )
			{
				for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
				{
					if ( jold != local_i )
					{   
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew                 ] = N[ibf] * ( sin_A * StiffMatrix[iold - oldBlockSize][jold                 ] + cos_A * StiffMatrix[iold][jold                 ] );
							newStiffMatrix[inew + ibf][jnew +   newBlockSize] = N[ibf] * ( sin_A * StiffMatrix[iold - oldBlockSize][jold +   oldBlockSize] + cos_A * StiffMatrix[iold][jold +   oldBlockSize] );
							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize] = N[ibf] * ( sin_A * StiffMatrix[iold - oldBlockSize][jold + 2*oldBlockSize] + cos_A * StiffMatrix[iold][jold + 2*oldBlockSize] );
						}

						jnew++;
						continue;
					}

					for ( ib = 0; ib < N.size(); ib++ ) 
					{
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew +                  ib] = N[ibf] * N[ib] * sin_A * ( cos_A * StiffMatrix[iold - oldBlockSize][jold] - sin_A * StiffMatrix[iold - oldBlockSize][jold + oldBlockSize] ) +
								                                                     N[ibf] * N[ib] * cos_A * ( cos_A * StiffMatrix[iold               ][jold] - sin_A * StiffMatrix[iold               ][jold + oldBlockSize] );

							newStiffMatrix[inew + ibf][jnew +   newBlockSize + ib] = N[ibf] * N[ib] * sin_A * ( sin_A * StiffMatrix[iold - oldBlockSize][jold] + cos_A * StiffMatrix[iold - oldBlockSize][jold + oldBlockSize] ) +
								                                                     N[ibf] * N[ib] * cos_A * ( sin_A * StiffMatrix[iold               ][jold] + cos_A * StiffMatrix[iold               ][jold + oldBlockSize] );

							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize + ib] = N[ibf] * N[ib] * ( sin_A * StiffMatrix[iold - oldBlockSize][jold + 2*oldBlockSize] + cos_A * StiffMatrix[iold][jold + 2*oldBlockSize] ); 
						}
					}

					jnew = jnew + N.size();
				}

				inew = inew + N.size();
			}
			else if ( iold == (local_i + 2*oldBlockSize) )
			{
				for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
				{
					if ( jold != local_i )
					{   
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew                 ] = N[ibf] * StiffMatrix[iold][jold                 ];
							newStiffMatrix[inew + ibf][jnew +   newBlockSize] = N[ibf] * StiffMatrix[iold][jold +   oldBlockSize];
							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize] = N[ibf] * StiffMatrix[iold][jold + 2*oldBlockSize]; 		
						}

						jnew++;
						continue;
					}

					for ( ib = 0; ib < N.size(); ib++ ) 
					{
						for ( ibf = 0; ibf < N.size(); ibf++ ) 
						{
							newStiffMatrix[inew + ibf][jnew +                  ib] = N[ibf] * N[ib] * ( cos_A * StiffMatrix[iold][jold] - sin_A * StiffMatrix[iold][jold + oldBlockSize] );
							newStiffMatrix[inew + ibf][jnew +   newBlockSize + ib] = N[ibf] * N[ib] * ( sin_A * StiffMatrix[iold][jold] + cos_A * StiffMatrix[iold][jold + oldBlockSize] );
							newStiffMatrix[inew + ibf][jnew + 2*newBlockSize + ib] = N[ibf] * N[ib] * StiffMatrix[iold][jold + 2*oldBlockSize]; 
						}
					}

					jnew = jnew + N.size();
				}

				inew = inew + N.size();
			}
			else
			{
				for ( jold = 0, jnew = 0; jold < oldBlockSize; jold++ )
				{   
					if ( jold != local_i )
					{
						newStiffMatrix[inew][jnew                 ] = StiffMatrix[iold][jold                 ]; 
						newStiffMatrix[inew][jnew +   newBlockSize] = StiffMatrix[iold][jold +   oldBlockSize]; 
						newStiffMatrix[inew][jnew + 2*newBlockSize] = StiffMatrix[iold][jold + 2*oldBlockSize]; 
						
						jnew++;
						continue;
					}
						
					for ( ib = 0; ib < N.size(); ib++ ) 
					{
						newStiffMatrix[inew][jnew +                  ib] = N[ib] * ( cos_A * StiffMatrix[iold][jold] - sin_A * StiffMatrix[iold][jold + oldBlockSize] );
						newStiffMatrix[inew][jnew +   newBlockSize + ib] = N[ib] * ( sin_A * StiffMatrix[iold][jold] + cos_A * StiffMatrix[iold][jold + oldBlockSize] );
						newStiffMatrix[inew][jnew + 2*newBlockSize + ib] = N[ib] * StiffMatrix[iold][jold + 2*oldBlockSize]; 
					}

					jnew = jnew + N.size();
				}

				inew++;
			}
		}

		StiffMatrix.Resize( newSize, newSize );
		
		for (int i = 0; i < newSize; i++ )
		{
			for (int j = 0; j < newSize; j++ )
			{
		    	StiffMatrix[i][j] = newStiffMatrix[i][j];
			}	
		}
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC nodes
    //*************************************************************************************************
	void Modeler::PBC_RightLeft_TR(std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
	    if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 ) return;
	    
		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );
		
		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if      ( mPBC_NodeNodePairs_Right.find(*ohoit) != mPBC_NodeNodePairs_Right.end() )
            {
				PBC_RL_NodeToNode_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );
            }
			else if ( mPBC_NodeEdgePairs_Right.find(*ohoit) != mPBC_NodeEdgePairs_Right.end() )
			{
				PBC_RL_NodeToEdge_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );			
			}
			else if ( mPBC_NodeElementPairs_Right.find(*ohoit) != mPBC_NodeElementPairs_Right.end() )
			{
				PBC_RL_NodeToElement_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );			
			}
            local_i++; 
        }
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to node
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToNode_TR(int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		int lNode = mPBC_NodeNodePairs_Right[NodeId];
		int size  = IdVector.size();
		int step  = size/3;
	
		Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

		HONodesId[local_i         ] = plNode->Id();

		IdVector [local_i         ] = plNode->pDofcEx()->EquationId();
		IdVector [local_i +   step] = plNode->pDofcEy()->EquationId();
		IdVector [local_i + 2*step] = plNode->pDofcEz()->EquationId();   

		std::vector<double> N(1); 
		
		N[0] = 1.00;
		
		PBC_ResVectorReDef_Rot(local_i, size, size, N, NodeId, ResVector);
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to edge
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToEdge_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<unsigned int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );	

		std::vector<Node::Pointer> plNodes; 
		
		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		double L2 = mPBC_NodeEdgeCoord_Right[NodeId];
        double L1 = 1.00 - L2;

		std::vector<double> N(3);

		N[0] = L1 * (2.00*L1 - 1.00);
		N[1] = L2 * (2.00*L2 - 1.00);
		N[2] = 4.00 * L1 * L2;

		PBC_ResVectorReDef_Rot(local_i, oldSize, newSize, N, NodeId, ResVector);

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<unsigned int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );	

		std::vector<Node::Pointer> plNodes; 
		
		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

		std::vector<double> natCoord(mPBC_NodeElementCoord_Right[NodeId]);
		
		double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;
		
		std::vector<double> N(6);

        N[0] = L1 * (2.00*L1 - 1.00);
        N[1] = L2 * (2.00*L2 - 1.00);
        N[2] = L3 * (2.00*L3 - 1.00);
		N[3] = 4.00 * L1 * L2;
		N[4] = 4.00 * L2 * L3;
        N[5] = 4.00 * L3 * L1;

		PBC_ResVectorReDef_Rot(local_i, oldSize, newSize, N, NodeId, ResVector);

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//************************************************************************************************* 
    //* - Trans(T)*b, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_ResVectorReDef_Rot(int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Vector<std::complex<double> >& ResVector)
	{
		double cos_A = 1.0;
		double sin_A = 0.0;

		if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );
	
		int inew, iold, ib;

		int oldBlockSize = oldSize / 3;
		int newBlockSize = newSize / 3;

		Vector<std::complex<double> > newResVector( newSize );

		for ( iold = 0, inew = 0; iold < oldSize; iold++ )
		{
			if ( iold == local_i )
			{
				for ( ib = 0; ib < N.size(); ib++ ) 
				{
					newResVector[inew + ib] = N[ib] * ( cos_A * ResVector[iold] - sin_A * ResVector[iold + oldBlockSize] ); 
				}
				inew = inew + N.size();
			}
			else if ( iold == (local_i + oldBlockSize) )
			{
				for ( ib = 0; ib < N.size(); ib++ ) 
				{
					newResVector[inew + ib] = N[ib] * ( sin_A * ResVector[iold - oldBlockSize] + cos_A * ResVector[iold] ); 
				}
				inew = inew + N.size();
			}
			else if ( iold == (local_i + 2*oldBlockSize) )
			{
				for ( ib = 0; ib < N.size(); ib++ ) 
				{
					newResVector[inew + ib] = N[ib] * ResVector[iold]; 
				}
				inew = inew + N.size();
			}
			else
			{
				newResVector[inew] = ResVector[iold]; 
				inew++;
			}
		}

		ResVector.resize( newSize ); 
		
		for (int i = 0; i < newSize; i++ )
		{
			ResVector[i] = newResVector[i];	
		}
	}

    //************************************************************************************************* 
    //* - Inv(T)*StiffnesMatrix*T, contact nodes
    //*************************************************************************************************
    void Modeler::ContactHO_TKT(std::vector<int>& HONodesId,
		                        Vector<int>& IdVector,
								Matrix<std::complex<double> >& StiffMatrix)
    {
		if ( mContactPairs.size() == 0 ) return;
		
		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mContactPairs.find(*hoit) != mContactPairs.end() )
            {
                int i;

                int step    = HONodesId.size();
				int step2   = 2*step;
                int dofSize = 3*step;

				ContactPairData& pContactPair = mContactPairs[*hoit];

				Node::Pointer pNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

				HONodesId[local_i        ] = pNode->Id();

				IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
				IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
				IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();
				
				///////////////  Material properties matrix  //////////////////
				Properties::Pointer Properties;

				double eo = 8.8541878176e-12;
        
				double freq = mProblemFrequency;

                double sigma, eps_real, eps_imag;

				//Material R
				Properties = mpModel->GetProperties(pContactPair.matR);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

				//Material L
				Properties = mpModel->GetProperties(pContactPair.matL);

				sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
                eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

				//Properties matrix
				std::complex<double> ContactMatrix[3][3];

                std::complex<double> rmat = (cEpsL/cEpsR) - 1.00;

				Vector<double> n = pContactPair.ctcNormal;

				ContactMatrix[0][0] = (n[0]*n[0]) * rmat + 1.00; ContactMatrix[0][1] = (n[0]*n[1]) * rmat       ; ContactMatrix[0][2] = (n[0]*n[2]) * rmat;
				ContactMatrix[1][0] = (n[1]*n[0]) * rmat       ; ContactMatrix[1][1] = (n[1]*n[1]) * rmat + 1.00; ContactMatrix[1][2] = (n[1]*n[2]) * rmat;
				ContactMatrix[2][0] = (n[2]*n[0]) * rmat       ; ContactMatrix[2][1] = (n[2]*n[1]) * rmat       ; ContactMatrix[2][2] = (n[2]*n[2]) * rmat + 1.00;

                ////////////////  Inv(T)*StiffnesMatrix*T  //////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                //StiffnesMatrix*T
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[i][local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[i][local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

                //Inv(T)*StiffnesMatrix
                for (i=0; i<dofSize; i++)
                { 
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

   					StiffMatrix[local_i        ][i] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[local_i + step ][i] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[local_i + step2][i] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

            }//end if( mNormals.find(*hoit) != mNormals.end() )

            local_i++; 
        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 	
    }

	//************************************************************************************************* 
    //* - Inv(T)*b, contact nodes
    //*************************************************************************************************
    void Modeler::ContactHO_TR(std::vector<int>& HONodesId,
		                       Vector<int>& IdVector, 
							   Vector<std::complex<double> >& ResVector)
    {
        if ( mContactPairs.size() == 0 ) return;
		
		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mContactPairs.find(*hoit) != mContactPairs.end() )
            {
                int step    = HONodesId.size();
				int step2   = 2*step;

				ContactPairData& pContactPair = mContactPairs[*hoit];

				Node::Pointer pNode = mpModel->GetNode( (pContactPair.ctcNode) - 1 );

				HONodesId[local_i        ] = pNode->Id();

				IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
				IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
				IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();
				
				///////////////  Material properties matrix  //////////////////
				Properties::Pointer Properties;

				double eo = 8.8541878176e-12;

				double freq = mProblemFrequency;

                double sigma, eps_real, eps_imag;

				//Material R
				Properties = mpModel->GetProperties(pContactPair.matR);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

				//Material L
				Properties = mpModel->GetProperties(pContactPair.matL);

				sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
                eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

				//Properties matrix
				std::complex<double> ContactMatrix[3][3];

                std::complex<double> rmat = (cEpsL/cEpsR) - 1.00;

				Vector<double> n = pContactPair.ctcNormal;

				ContactMatrix[0][0] = (n[0]*n[0]) * rmat + 1.00; ContactMatrix[0][1] = (n[0]*n[1]) * rmat       ; ContactMatrix[0][2] = (n[0]*n[2]) * rmat;
				ContactMatrix[1][0] = (n[1]*n[0]) * rmat       ; ContactMatrix[1][1] = (n[1]*n[1]) * rmat + 1.00; ContactMatrix[1][2] = (n[1]*n[2]) * rmat;
				ContactMatrix[2][0] = (n[2]*n[0]) * rmat       ; ContactMatrix[2][1] = (n[2]*n[1]) * rmat       ; ContactMatrix[2][2] = (n[2]*n[2]) * rmat + 1.00;

                //////////////  Inv(T)*b  //////////////////
                std::complex<double> dof_x, dof_y, dof_z;
               
				dof_x = ResVector[local_i        ];
                dof_y = ResVector[local_i + step ];
                dof_z = ResVector[local_i + step2];

                ResVector[local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                ResVector[local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                ResVector[local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
            }

            local_i++; 
        }  
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - PEC, PMC, PBC and projection elements
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns each PBC node to its corresponding std::map
    ////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::GeneratePBCElement(int* NodesId, unsigned int PropertiesId)
    {
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);
        
		int FaceTypeId = (*properties)(COMPLEX_IBC)[0];
        int MaterialId = NodesId[3];

		mPBCTolerance = (*properties)(COMPLEX_IBC)[1];

		std::vector<int> HONodesId;
		PushHONodesOnSurface(NodesId, HONodesId);

        //It gives an error if we assign directly HONodesId to std::map with push_back (¿?)
		int numberOfNodes = HONodesId.size();
		std::vector<unsigned int> PBCNodes(numberOfNodes);
		for (int id = 0; id < numberOfNodes; id++) PBCNodes[id] = HONodesId[id];
		
        //Front face
		if      (FaceTypeId == 11)
		{
		    mPBCFrontElements[MaterialId].push_back(PBCNodes);
		}
		//Back face
		else if (FaceTypeId == 12)
		{
			mPBCBackElements[MaterialId].push_back(PBCNodes);
		}
		//Right face 
		else if (FaceTypeId == 21)
		{
		    mPBCRightElements[MaterialId].push_back(PBCNodes);
		}
		//Left face 
		else if (FaceTypeId == 22)
		{
			mPBCLeftElements[MaterialId].push_back(PBCNodes);
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns the nodes in a periodic surface to the elements in the other
	//   corresponding periodic surface.
	// - Only implemented for 2nd order elements.
	// - Geometry must be placed carefully: 
	//   1) Front-Back surfaces in the XY-plane.
	//   2) Left surface in the 00-YZ-Plane. 
	//	 3) For cylindrical symmetry the central axis must be placed along the Z axis.
    ///////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::SetPBC()
	{
		if ( mElementOrder != 2 )
		{
			//std::cout << std::endl; 
			//std::cout << "WARNING: PBC condition is only implemented for 2nd order elements." << std::endl; 
			//std::cout << "!!!!!!!!!!!!!!!!!!!!!! Change element order !!!!!!!!!!!!!!!!!!!!!!" << std::endl;
			//std::cout << std::endl;
			mPBCFrontElements.clear(); mPBCFrontElements.swap(std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCFrontElements));
			mPBCBackElements.clear (); mPBCBackElements.swap (std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCBackElements) );
			mPBCRightElements.clear(); mPBCRightElements.swap(std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCRightElements));
			mPBCLeftElements.clear (); mPBCLeftElements.swap (std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCLeftElements) );
			return;
		}
		
		if ( mPBCFrontElements.size() > 0 ) Set_FrontBack_PBC();

        if ( mPBCRightElements.size() > 0 ) Set_RightLeft_PBC();

        ////Debugging /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//std::map<unsigned int, unsigned int>::iterator it_nn;
		//for(it_nn = mPBC_NodeNodePairs_Front.begin(); it_nn != mPBC_NodeNodePairs_Front.end(); ++it_nn) 
		//{
		//	std::cout << "Frontal: " << it_nn->first << " , " << it_nn->second << std::endl;
		//}
        //      
		//std::map<unsigned int, std::vector<unsigned int> >::iterator it_neg;
		//std::map<unsigned int, double>::iterator it_edco = mPBC_NodeEdgeCoord_Front.begin();
		//for(it_neg = mPBC_NodeEdgePairs_Front.begin(); it_neg != mPBC_NodeEdgePairs_Front.end(); ++it_neg, ++it_edco)
		//{	
		//	std::cout << it_neg->first << " , "; 
		//	std::cout << (it_neg->second)[0] << " , " << (it_neg->second)[1] << " , " << (it_neg->second )[2] << "   ";
		//	std::cout << "FNatX: " << std::setprecision(12) << (it_edco->second) << "   " << std::endl;
		//}
		//	
		//std::map<unsigned int, std::vector<unsigned int> >::iterator it_nel;
		//std::map<unsigned int, std::vector<double> >::iterator it_elco =  mPBC_NodeElementCoord_Front.begin();
		//for(it_nel = mPBC_NodeElementPairs_Front.begin(); it_nel != mPBC_NodeElementPairs_Front.end(); ++it_nel, ++it_elco) 
		//{			
		//	std::vector<Node::Pointer> pNodes;
		//	std::vector<unsigned int> NodesIds(it_nel->second);
		//	std::vector<unsigned int>::iterator it_nId;
		//	for(it_nId = NodesIds.begin(); it_nId != NodesIds.end(); ++it_nId) pNodes.push_back( mpModel->GetNode( *it_nId - 1 ) );
		//	double Xi, Yi; nCoordToXY( (it_elco->second)[0], (it_elco->second)[1], Xi, Yi, pNodes ); 
		//	
		//	std::cout << it_nel->first << " , "; 
		//	std::cout << (it_nel->second)[0] << " , " << (it_nel->second)[1] << " , " << (it_nel->second)[2] << " , ";
		//	std::cout << (it_nel->second)[3] << " , " << (it_nel->second)[4] << " , " << (it_nel->second)[5] << "   ";
		//	std::cout << "FNatXY: " << std::setprecision(12) << (it_elco->second)[0] << " , " << (it_elco->second)[1];	
		//	std::cout << "   FXY: " << std::setprecision(12) << Xi << " , " << Yi << std::endl;	
		//}
		//
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		//std::map<unsigned int, unsigned int>::iterator it_nnr;
		//for(it_nnr = mPBC_NodeNodePairs_Right.begin(); it_nnr != mPBC_NodeNodePairs_Right.end(); ++it_nnr) 
		//{
		//    std::cout << it_nnr->first << " , " << it_nnr->second << "   ";
		//	std::cout << "YZrot: " << std::setprecision(12) << PBC_RotateY(mpModel->GetNode((it_nnr->first)-1)) << " , " 
		//		                                            << mpModel->GetNode((it_nnr->first)-1)->Z()         << std::endl;
		//}
        //       
		//std::map<unsigned int, std::vector<unsigned int> >::iterator it_negr;
		//std::map<unsigned int, double>::iterator it_edcor = mPBC_NodeEdgeCoord_Right.begin();
		//for(it_negr = mPBC_NodeEdgePairs_Right.begin(); it_negr != mPBC_NodeEdgePairs_Right.end(); ++it_negr, ++it_edcor)
		//{			
		//	std::cout << it_negr->first << " , "; 
		//	std::cout << (it_negr->second)[0] << " , " << (it_negr->second)[1] << " , " << (it_negr->second )[2] << "   ";
		//	std::cout << "RLNatX: " << std::setprecision(12) << (it_edcor->second) << "   " ;
		//	std::cout << "YZrot: "  << std::setprecision(12) << PBC_RotateY(mpModel->GetNode((it_negr->first)-1)) << " , " 
		//		                                             << mpModel->GetNode((it_negr->first)-1)->Z()         << std::endl;
		//}
		//	
		//std::map<unsigned int, std::vector<unsigned int> >::iterator it_nelr;
		//std::map<unsigned int, std::vector<double> >::iterator it_elcor =  mPBC_NodeElementCoord_Right.begin();
		//for(it_nelr = mPBC_NodeElementPairs_Right.begin(); it_nelr != mPBC_NodeElementPairs_Right.end(); ++it_nelr, ++it_elcor) 
		//{
		//	std::vector<Node::Pointer> pNodes;
		//	std::vector<unsigned int> NodesIds(it_nelr->second);
		//	std::vector<unsigned int>::iterator it_nId;
		//	for(it_nId = NodesIds.begin(); it_nId != NodesIds.end(); ++it_nId) pNodes.push_back( mpModel->GetNode( *it_nId - 1 ) );
		//	double Zi, Yi; nCoordToZY( (it_elcor->second)[0], (it_elcor->second)[1], Zi, Yi, pNodes ); 
		//	
		//	std::cout << it_nelr->first << " , "; 
		//	std::cout << (it_nelr->second)[0] << " , " << (it_nelr->second)[1] << " , " << (it_nelr->second)[2] << " , ";
		//	std::cout << (it_nelr->second)[3] << " , " << (it_nelr->second)[4] << " , " << (it_nelr->second)[5] << "   ";
		//	std::cout << "RLNatXY: " << std::setprecision(12) << (it_elcor->second)[0] << " , " << (it_elcor->second)[1] << "   ";
		//	std::cout << "YZele: " << std::setprecision(12) << Yi << " , " << Zi << "   "; 
		//	std::cout << "YZrot: " << PBC_RotateY( mpModel->GetNode((it_nelr->first)-1) ) << " , " << ( mpModel->GetNode((it_nelr->first)-1) )->Z() << std::endl;
		//}
		//
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns the nodes in the front periodic surface to the elements in the back periodic surface.
	// - Front-back surfaces must be placed in the XY-plane.
	// - Only for 2nd order elements.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Set_FrontBack_PBC()
    {
		//Loop through frontal elements for each material
		std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it_Front;
        for(it_Front = mPBCFrontElements.begin(); it_Front != mPBCFrontElements.end(); ++it_Front) 
        {
			int pbcMaterial = it_Front->first;

			std::vector<std::vector<unsigned int> > elemSetFront(it_Front->second);
			
			//Loop through frontal elements with pbcMaterial
			std::vector<std::vector<unsigned int> >::iterator it_nodesF;
			for(it_nodesF = elemSetFront.begin(); it_nodesF != elemSetFront.end(); ++it_nodesF) 
            {
				std::vector<unsigned int> eleFront(*it_nodesF);
				
				//Loop through nodes in one frontal element
				std::vector<unsigned int>::iterator it_nF;
				for(it_nF = eleFront.begin(); it_nF != eleFront.end(); ++it_nF) 
				{
					if (mPBC_NodeNodePairs_Front.find   (*it_nF) != mPBC_NodeNodePairs_Front.end   ()) continue;
					if (mPBC_NodeEdgePairs_Front.find   (*it_nF) != mPBC_NodeEdgePairs_Front.end   ()) continue;
					if (mPBC_NodeElementPairs_Front.find(*it_nF) != mPBC_NodeElementPairs_Front.end()) continue;

					//Front node to be compared
					Node::Pointer pfNode = mpModel->GetNode( (*it_nF) - 1 );

					//Frontal node FixDof
					(*mpModel)(cEx, *pfNode) = czero; pfNode->pDofcEx()->FixDof();
					(*mpModel)(cEy, *pfNode) = czero; pfNode->pDofcEy()->FixDof();
					(*mpModel)(cEz, *pfNode) = czero; pfNode->pDofcEz()->FixDof();

					//Back elements with pbcMaterial
					std::vector<std::vector<unsigned int> > elemSetBack( mPBCBackElements[pbcMaterial] );

					//Nearest node to Xo,Yo
					unsigned int centralNodeId = PBC_FB_SearchNearestNode( pfNode, elemSetBack );

					//Check if the nearest node is the frontal node Xo,Yo
					if ( FrontNodeIsBackNode( pfNode, centralNodeId ) == true ) continue;

					//Search elements around nearest node 
					std::vector<std::vector<unsigned int> > elemSetNearFNode;
                    
					FillElemSetNearNode( centralNodeId, elemSetBack, elemSetNearFNode );

					//Find nearest point inside each element of elemSet
					std::vector<std::vector<double> > minNatCoordSet;
					std::vector<unsigned int>         minEleNodesIds;
					std::vector<double>               minNatCoord;

					PBC_FB_NearestPointInElementSet( pfNode, elemSetNearFNode, minNatCoordSet, minEleNodesIds, minNatCoord );
					
					while ( 1 )
					{
						if ( NearPointIsInExtEdge( centralNodeId, minEleNodesIds, minNatCoord ) == true )
						{
							FillElemSetNearNode( centralNodeId, elemSetBack, elemSetNearFNode );
							PBC_FB_NearestPointInElementSet( pfNode, elemSetNearFNode, minNatCoordSet, minEleNodesIds, minNatCoord );
						}
						else
						{
							break;
						}
					};

					double prec = 0.01;

					while ( prec >= mPBCTolerance )
					{
						if ( ( abs( minNatCoord[0] ) < mPBCTolerance ) || 
							 ( abs( minNatCoord[1] ) < mPBCTolerance ) ||
					         ( abs( minNatCoord[0] + minNatCoord[1] - 1.0 ) < mPBCTolerance ) )
						{
							PBC_FB_SearchMinPointInEdges( pfNode, prec, elemSetNearFNode, minNatCoordSet, minEleNodesIds, minNatCoord );
						}
						else
						{
							PBC_FB_SearchMinPointInElementToTol( pfNode, prec, minEleNodesIds, minNatCoord );
							break;
						}

						prec = prec * 0.1;
					}

					PBC_FB_SetPointInElement( pfNode, minEleNodesIds, minNatCoord );

				}//End of loop through nodes in one frontal element

			}//End of loop through frontal elements with pbcMaterial

		}//End of loop through frontal elements 

	    mPBCFrontElements.clear(); mPBCFrontElements.swap( std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCFrontElements) );
		mPBCBackElements.clear (); mPBCBackElements.swap ( std::map<unsigned int, std::vector<std::vector<unsigned int> > >(mPBCBackElements ) );
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Check if the near point is in an exterior edge.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool Modeler::NearPointIsInExtEdge( unsigned int& centralNodeId, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord )
	{
		if ( ( abs( minNatCoord[0] ) < mPBCTolerance ) && 
			 ( minEleNodesIds[0] != centralNodeId    ) && 
			 ( minEleNodesIds[2] != centralNodeId    ) && 
			 ( minEleNodesIds[5] != centralNodeId    )   )
		{
		    centralNodeId = minEleNodesIds[5];
			return true;
		}
		else if ( ( abs( minNatCoord[1] ) < mPBCTolerance ) && 
			      ( minEleNodesIds[0] != centralNodeId    ) && 
			      ( minEleNodesIds[1] != centralNodeId    ) && 
			      ( minEleNodesIds[3] != centralNodeId    )   )
		{
		    centralNodeId = minEleNodesIds[3];
			return true;
		}
		else if ( ( abs( minNatCoord[0] + minNatCoord[1] - 1.00 ) < mPBCTolerance ) && 
			      ( minEleNodesIds[1] != centralNodeId                            ) && 
			      ( minEleNodesIds[2] != centralNodeId                            ) && 
			      ( minEleNodesIds[4] != centralNodeId                            )   )
		{
			centralNodeId = minEleNodesIds[4];
			return true;
		}
		else 
		{
		    return false;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Set the point Xo,Yo in the element
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_FB_SetPointInElement( Node::Pointer pfNode, std::vector<unsigned int>& bNodesIds, std::vector<double>& natCoord)
	{
		if      ( abs( natCoord[0] ) < mPBCTolerance )  
		{ 
			if      ( abs( natCoord[1]       ) < mPBCTolerance ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[0];
			else if ( abs( natCoord[1] - 0.5 ) < mPBCTolerance ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[5];
			else if ( abs( natCoord[1] - 1.0 ) < mPBCTolerance ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[2];
			else                                                 PBC_FB_SetPointInEdge( pfNode, bNodesIds, 0.0, natCoord[1] );	
		}
		else if ( abs( natCoord[1] ) < mPBCTolerance )  
		{ 
			if      ( abs( natCoord[0] - 0.5 ) < mPBCTolerance ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[3];
			else if ( abs( natCoord[0] - 1.0 ) < mPBCTolerance ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[1];
			else                                                 PBC_FB_SetPointInEdge( pfNode, bNodesIds, natCoord[0], 0.0 );
		}
		else if ( abs( natCoord[0] + natCoord[1] - 1.0 ) < mPBCTolerance )  
		{ 
			if ( ( abs( natCoord[0] - 0.5 ) < mPBCTolerance ) && 
				 ( abs( natCoord[1] - 0.5 ) < mPBCTolerance )    ) mPBC_NodeNodePairs_Front[pfNode->Id()] = bNodesIds[4];
			else                                                   PBC_FB_SetPointInEdge( pfNode, bNodesIds, natCoord[0], natCoord[1] );
		}
		else
		{
			mPBC_NodeElementPairs_Front[pfNode->Id()] = bNodesIds;
			mPBC_NodeElementCoord_Front[pfNode->Id()] = natCoord;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Set the point of natural coordinates nX,nY in the edge of the element 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_FB_SetPointInEdge( Node::Pointer pfNode, std::vector<unsigned int>& bNodesIds, double nx, double ny)
	{
        int    n1,n2,n3;
		double edgeCoord; 

		if      ( nx == 0.0 ) { n1 = 2; n2 = 0; n3 = 5; edgeCoord = 1 - ny; }
		else if ( ny == 0.0 ) { n1 = 0; n2 = 1; n3 = 3; edgeCoord = nx    ; }
		else                  { n1 = 1; n2 = 2; n3 = 4; edgeCoord = 1 - nx; }

		std::vector<unsigned int> edgeNodesIds(3);
		edgeNodesIds[0] = bNodesIds[n1];
		edgeNodesIds[1] = bNodesIds[n2];
		edgeNodesIds[2] = bNodesIds[n3];

		mPBC_NodeEdgePairs_Front[pfNode->Id()] = edgeNodesIds;
        mPBC_NodeEdgeCoord_Front[pfNode->Id()] = edgeCoord;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Transforms natural coordinates nx, ny of element pNodes in cartesian coordinates X,Y
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::nCoordToXY (double nx, double ny, double& X, double& Y, std::vector<Node::Pointer>& pNodes)
	{
		double L2 = nx; 
		double L3 = ny; 
		double L1 = 1.00 - L2 - L3;

		double N0 = L1 * (2.00*L1 - 1.00);
        double N1 = L2 * (2.00*L2 - 1.00);
        double N2 = L3 * (2.00*L3 - 1.00);
		double N3 = 4.00 * L1 * L2;
		double N4 = 4.00 * L2 * L3;
		double N5 = 4.00 * L3 * L1;

		X = N0 * pNodes[0]->X() + N1 * pNodes[1]->X() + N2 * pNodes[2]->X() +
			N3 * pNodes[3]->X() + N4 * pNodes[4]->X() + N5 * pNodes[5]->X() ;

		Y = N0 * pNodes[0]->Y() + N1 * pNodes[1]->Y() + N2 * pNodes[2]->Y() +
			N3 * pNodes[3]->Y() + N4 * pNodes[4]->Y() + N5 * pNodes[5]->Y() ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Xo,Yo in each element of elemSet
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_FB_NearestPointInElementSet( Node::Pointer pfNode, std::vector<std::vector<unsigned int> >& elemSet, 
												   std::vector<std::vector<double> >& minNatCoordSet,
												   std::vector<unsigned int>&         minEleNodesIds,
					                               std::vector<double>&               minNatCoord )
	{
		minNatCoordSet.clear(); minNatCoordSet.swap( std::vector<std::vector<double> > ( minNatCoordSet ) );
		minEleNodesIds.clear(); minEleNodesIds.swap( std::vector<unsigned int>         ( minEleNodesIds ) );
		minNatCoord.clear()   ; minNatCoord.swap   ( std::vector<double>               ( minNatCoord    ) );

		std::vector<unsigned int> iniNodesIds( *elemSet.begin() );

		for( int i = 0; i < iniNodesIds.size(); i++ ) 
		{ 
			minEleNodesIds.push_back( iniNodesIds[i] );
		}

		minNatCoord.push_back( 0.0 ); 
		minNatCoord.push_back( 0.0 );

		Node::Pointer iniNode = mpModel->GetNode( iniNodesIds[0] - 1 );

        double distToEle;
		double minDistToEle = abs( pfNode->X() - iniNode->X() ) + abs( pfNode->Y() - iniNode->Y() );

		double nx_ini = 0.0; double nx_end = 1.0; 
		double ny_ini = 0.0; double ny_end = 1.0; 
		double step   = 0.1;

		std::vector<std::vector<unsigned int> >::iterator it_eSet;

		for(it_eSet = elemSet.begin(); it_eSet != elemSet.end(); ++it_eSet) 
		{
			std::vector<unsigned int> eNodesIds(*it_eSet); 
			std::vector<double>       mdnCoord (2, 0.00);

			distToEle = PBC_FB_MinDistPointInEle( pfNode, eNodesIds, nx_ini, nx_end, ny_ini, ny_end, step, mdnCoord );

			minNatCoordSet.push_back( mdnCoord );

			if ( distToEle <= minDistToEle )
			{
				minDistToEle = distToEle;				
				
				minNatCoord[0] = mdnCoord[0]; 
				minNatCoord[1] = mdnCoord[1];
				
				for( int i = 0; i < eNodesIds.size(); i++ ) minEleNodesIds[i] = eNodesIds[i];
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Xo,Yo in each element of elemSet when the point is in an edge
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_FB_SearchMinPointInEdges( Node::Pointer pfNode, double prec, 
		                                        std::vector<std::vector<unsigned int> >& elemSet, 
												std::vector<std::vector<double> >& minNatCoordSet,
												std::vector<unsigned int>&         minEleNodesIds,
					                            std::vector<double>&               minNatCoord )
	{
		Node::Pointer iniNode = mpModel->GetNode( minEleNodesIds[0] - 1 );
		
		double distToEdge;
		double minDistToEdge = abs( pfNode->X() - iniNode->X() ) + abs( pfNode->Y() - iniNode->Y() );

        double incInt = prec * 10.0;

		std::vector<std::vector<unsigned int> >::iterator it_eSet;

		std::vector<std::vector<double> >::iterator it_nC = minNatCoordSet.begin();

		for( it_eSet = elemSet.begin(); it_eSet != elemSet.end(); ++it_eSet, ++it_nC ) 
		{
			std::vector<unsigned int> eNodesIds( *it_eSet );
			std::vector<double>       minCoord ( *it_nC   );

            double nx_ini = minCoord[0] - incInt; 
			double nx_end = minCoord[0] + incInt; 

			double ny_ini = minCoord[1] - incInt; 
			double ny_end = minCoord[1] + incInt; 

			distToEdge = PBC_FB_MinDistPointInEle( pfNode, eNodesIds, nx_ini, nx_end, ny_ini, ny_end, prec, minCoord );

			(*it_nC)[0] = minCoord[0]; 
			(*it_nC)[1] = minCoord[1];

			if ( distToEdge <= minDistToEdge )
			{
				minDistToEdge = distToEdge;				
				
				minNatCoord[0] = minCoord[0]; 
				minNatCoord[1] = minCoord[1];
				
				for( int i = 0; i < eNodesIds.size(); i++ ) minEleNodesIds[i] = eNodesIds[i];
			}
		}	
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Xo,Yo in an element up to tolerance
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_FB_SearchMinPointInElementToTol( Node::Pointer pfNode, double prec, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord )
	{
		do 
		{
			double incInt = prec * 10.0;

			double nx_ini = minNatCoord[0] - incInt; 
			double nx_end = minNatCoord[0] + incInt; 

			double ny_ini = minNatCoord[1] - incInt; 
			double ny_end = minNatCoord[1] + incInt; 
				
			PBC_FB_MinDistPointInEle( pfNode, minEleNodesIds, nx_ini, nx_end, ny_ini, ny_end, prec, minNatCoord );

			prec = prec * 0.1;	

		}while ( prec >= mPBCTolerance );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Xo,Yo in an element
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Modeler::PBC_FB_MinDistPointInEle( Node::Pointer pfNode, std::vector<unsigned int>& NodesIds, 
		                                      double nx_ini, double nx_end, 
											  double ny_ini, double ny_end, 
											  double step, 
											  std::vector<double>& natCoord )
	{
		std::vector<Node::Pointer> peNodes;

		std::vector<unsigned int>::iterator it_nId;

		for(it_nId = NodesIds.begin(); it_nId != NodesIds.end(); ++it_nId) 
		{
			peNodes.push_back( mpModel->GetNode( *it_nId - 1 ) );
		}

		double Xo = pfNode->X();
		double Yo = pfNode->Y();
	
		double nx, ny;
		double Xi, Yi;
		
		double dist;
		double min_dist = abs( Xo - peNodes[0]->X() ) + abs ( Yo - peNodes[0]->Y() );

		for ( nx = nx_ini; nx <= nx_end; nx = nx + step ) 
		{ 
			for ( ny = ny_ini; ny <= ny_end; ny = ny + step ) 
			{ 
				if ( ( ( nx + ny ) > 1.0 ) || ( nx < 0.0 ) || ( ny < 0.0 ) ) continue;

				nCoordToXY(nx, ny, Xi, Yi, peNodes); 
				
				dist = abs( Xo - Xi ) + abs ( Yo - Yi );

				if ( dist <= min_dist )
				{
					min_dist    = dist;
					natCoord[0] = nx;
					natCoord[1] = ny;
				}
			}
		}

		return min_dist;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search the nearest back node to Xo,Yo
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int Modeler::PBC_FB_SearchNearestNode(Node::Pointer pfNode, std::vector<std::vector<unsigned int> >& elemSetBack)
	{
		double Xo = pfNode->X();
		double Yo = pfNode->Y();
		
		Node::Pointer iniNode = mpModel->GetNode( *(*elemSetBack.begin()).begin() - 1 );

		unsigned int  nearestNode = iniNode->Id();
		double        minDistance = abs( Xo - iniNode->X() ) + abs( Yo - iniNode->Y() );

		std::vector<std::vector<unsigned int> >::iterator it_elemB;

		for(it_elemB = elemSetBack.begin(); it_elemB != elemSetBack.end(); ++it_elemB) 
        {
            std::vector<unsigned int> elemNodesB(*it_elemB);
            
			std::vector<unsigned int>::iterator it_pbn;
			for(it_pbn = elemNodesB.begin(); it_pbn != elemNodesB.end(); ++it_pbn) 
			{
				Node::Pointer pbNode = mpModel->GetNode( (*it_pbn) - 1 );

				double eDist = abs( Xo - pbNode->X() ) + abs( Yo - pbNode->Y() );

				if ( eDist <= minDistance ) 
				{ 
					nearestNode = *it_pbn; 
					minDistance = eDist;
				} 
			}
		}

		return nearestNode;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Check if the nearest node is the frontal node Xo,Yo
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool Modeler::FrontNodeIsBackNode(Node::Pointer pfNode, unsigned int backNodeId)
	{
		Node::Pointer pbNode = mpModel->GetNode( backNodeId - 1 );

		double Xi = pbNode->X();
		double Yi = pbNode->Y();

		double Xo = pfNode->X();
		double Yo = pfNode->Y();

		bool FNodeIsBNode = ( abs(Xo - Xi) <= ( abs(Xo) * mPBCTolerance ) ) && 
			                ( abs(Yo - Yi) <= ( abs(Yo) * mPBCTolerance ) );

		if ( ( FNodeIsBNode == true ) || ( mPBCTolerance > 0.1 ) )
		{
			mPBC_NodeNodePairs_Front[pfNode->Id()] = backNodeId;
			return true;
		}
		else
		{
			return false;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Fill a subset of elements from another set of elements which have the node NodeId
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::FillElemSetNearNode( unsigned int NodeId, 
		                               std::vector<std::vector<unsigned int> >& elemSet, 
									   std::vector<std::vector<unsigned int> >& elemSetNearNode )
	{
		elemSetNearNode.clear(); elemSetNearNode.swap( std::vector<std::vector<unsigned int> > (elemSetNearNode) );

		std::vector<std::vector<unsigned int> >::iterator it_elemSet;

		for(it_elemSet = elemSet.begin(); it_elemSet != elemSet.end(); ++it_elemSet) 
        {
            std::vector<unsigned int> elemNodesIds(*it_elemSet);
            
			std::vector<unsigned int>::iterator it_eleNodId;

			for(it_eleNodId = elemNodesIds.begin(); it_eleNodId != elemNodesIds.end(); ++it_eleNodId) 
			{
				if ( NodeId == (*it_eleNodId) ) 
				{ 
					elemSetNearNode.push_back( elemNodesIds ); 
					break; 
				}
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns the nodes in the right periodic surface to the elements in the left periodic surface.
	// - Only implemented for 2nd order elements.
	// - Left surface must be placed in the 00-YZ-Plane. 
	// - For cylindrical symmetry the central axis must be placed along the Z axis.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::Set_RightLeft_PBC()
	{
		//Loop through right elements for each material
		std::map<unsigned int, std::vector<std::vector<unsigned int> > >::iterator it_Right;
        for(it_Right = mPBCRightElements.begin(); it_Right != mPBCRightElements.end(); ++it_Right) 
        {
			int pbcMaterial = it_Right->first;

			std::vector<std::vector<unsigned int> > elemSetRight(it_Right->second);
			
			//Loop through right elements with pbcMaterial
			std::vector<std::vector<unsigned int> >::iterator it_nodesR;
			for(it_nodesR = elemSetRight.begin(); it_nodesR != elemSetRight.end(); ++it_nodesR) 
            {
				std::vector<unsigned int> eleRight(*it_nodesR);

				double tXn0   = ( mpModel->GetNode( eleRight[0]-1 ) )->X();
				double tXn1   = ( mpModel->GetNode( eleRight[1]-1 ) )->X();
				double tXn2   = ( mpModel->GetNode( eleRight[2]-1 ) )->X();

				double gidTol = 1e-14;

				mIsRightPBCTilted = ( abs(tXn0-tXn1) > gidTol ) || ( abs(tXn1-tXn2) > gidTol ) || ( abs(tXn2-tXn0) > gidTol );
								
				//Loop through nodes in one right element
                std::vector<unsigned int>::iterator it_nR;
				for(it_nR = eleRight.begin(); it_nR != eleRight.end(); ++it_nR) 
				{
					if ( mPBC_NodeNodePairs_Right.find   (*it_nR) != mPBC_NodeNodePairs_Right.end   () ) continue;
					if ( mPBC_NodeEdgePairs_Right.find   (*it_nR) != mPBC_NodeEdgePairs_Right.end   () ) continue;
					if ( mPBC_NodeElementPairs_Right.find(*it_nR) != mPBC_NodeElementPairs_Right.end() ) continue;

					//Right node to be compared
					Node::Pointer prNode = mpModel->GetNode( (*it_nR) - 1 );

					//If mIsRightPBCTilted == true and node is in the Z-axis then Ex = Ey = 0 and continue
					if ( ( mIsRightPBCTilted == true ) && ( abs(prNode->X()) < gidTol ) && ( abs(prNode->Y()) < gidTol ) ) 
				    {   
						if (mNormals.find(*it_nR) == mNormals.end())
						{   
						    (*mpModel)(cEx, *prNode) = czero; prNode->pDofcEx()->FixDof();
						    (*mpModel)(cEy, *prNode) = czero; prNode->pDofcEy()->FixDof();
						}
						continue;
					}
	
					//Right node FixDof
					(*mpModel)(cEx, *prNode) = czero; prNode->pDofcEx()->FixDof();
					(*mpModel)(cEy, *prNode) = czero; prNode->pDofcEy()->FixDof();
					(*mpModel)(cEz, *prNode) = czero; prNode->pDofcEz()->FixDof();

					//Back elements with pbcMaterial
					std::vector<std::vector<unsigned int> > elemSetLeft( mPBCLeftElements[pbcMaterial] );

                    //Nearest node to Zo,Yo
					unsigned int centralNodeId = PBC_RL_SearchNearestNode( prNode, elemSetLeft );

					//Check if the nearest node is the right node Zo,Yo
					if ( RightNodeIsLeftNode( prNode, centralNodeId ) == true ) continue;
                    
					//Search elements around nearest node 
					std::vector<std::vector<unsigned int> > elemSetNearRNode;

					FillElemSetNearNode( centralNodeId, elemSetLeft, elemSetNearRNode );

					//Find nearest point inside each element of elemSet
					std::vector<std::vector<double> > minNatCoordSet;
					std::vector<unsigned int>         minEleNodesIds;
					std::vector<double>               minNatCoord;
					
					PBC_RL_NearestPointInElementSet( prNode, elemSetNearRNode, minNatCoordSet, minEleNodesIds, minNatCoord );
			
					while ( 1 )
					{
						if ( NearPointIsInExtEdge( centralNodeId, minEleNodesIds, minNatCoord ) == true )
						{
							FillElemSetNearNode( centralNodeId, elemSetLeft, elemSetNearRNode );
							PBC_RL_NearestPointInElementSet( prNode, elemSetNearRNode, minNatCoordSet, minEleNodesIds, minNatCoord );
						}
						else
						{
							break;
						}
					};

					double prec = 0.01;

					while ( prec >= mPBCTolerance )
					{
						if ( ( abs( minNatCoord[0] ) < mPBCTolerance ) || 
							 ( abs( minNatCoord[1] ) < mPBCTolerance ) ||
					         ( abs( minNatCoord[0] + minNatCoord[1] - 1.0 ) < mPBCTolerance ) )
						{
							PBC_RL_SearchMinPointInEdges( prNode, prec, elemSetNearRNode, minNatCoordSet, minEleNodesIds, minNatCoord );
						}
						else
						{
							PBC_RL_SearchMinPointInElementToTol( prNode, prec, minEleNodesIds, minNatCoord );
							break;
						}

						prec = prec * 0.1;
					}

					PBC_RL_SetPointInElement( prNode, minEleNodesIds, minNatCoord );
                    
				}//End of loop through nodes in one frontal element

			}//End of loop through frontal elements with pbcMaterial

		}//End of loop through frontal elements 
		
		mPBCRightElements.clear(); mPBCRightElements.swap( std::map<unsigned int, std::vector<std::vector<unsigned int> > >( mPBCRightElements ) );
		mPBCLeftElements.clear (); mPBCLeftElements.swap ( std::map<unsigned int, std::vector<std::vector<unsigned int> > >( mPBCLeftElements  ) );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search the nearest back node to Zo,Yo
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	unsigned int Modeler::PBC_RL_SearchNearestNode(Node::Pointer prNode, std::vector<std::vector<unsigned int> >& elemSetLeft)
	{
	    double Zo = prNode->Z();
		double Yo = PBC_RotateY(prNode);
		
		Node::Pointer iniNode = mpModel->GetNode( *(*elemSetLeft.begin()).begin() - 1 );

		unsigned int  nearestNode = iniNode->Id();
		double        minDistance = abs( Zo - iniNode->Z() ) + abs( Yo - iniNode->Y() );

		std::vector<std::vector<unsigned int> >::iterator it_elemL;

		for(it_elemL = elemSetLeft.begin(); it_elemL != elemSetLeft.end(); ++it_elemL) 
        {
            std::vector<unsigned int> elemNodesL(*it_elemL);
            
			std::vector<unsigned int>::iterator it_pln;

			for(it_pln = elemNodesL.begin(); it_pln != elemNodesL.end(); ++it_pln) 
			{
				Node::Pointer plNode = mpModel->GetNode( (*it_pln) - 1 );

				double eDist = abs( Zo - plNode->Z() ) + abs( Yo - plNode->Y() );

				if ( eDist <= minDistance ) 
				{ 
					nearestNode = *it_pln; 
					minDistance = eDist;
				} 
			}
		}

		return nearestNode;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Rotate right node around Z-axis from right plane to left plane (00-YZ-plane)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Modeler::PBC_RotateY(Node::Pointer prNode)
	{
		if ( mIsRightPBCTilted == false ) return ( prNode->Y() );

		double X = prNode->X();
		double Y = prNode->Y();

		double modXY2 = X*X + Y*Y; 

		return ( modXY2 / sqrt(modXY2) );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Check if the nearest node is the right node Zo,Yo
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool Modeler::RightNodeIsLeftNode(Node::Pointer prNode, unsigned int leftNodeId)
	{
		Node::Pointer plNode = mpModel->GetNode( leftNodeId - 1 );

		double Zi = plNode->Z();
		double Yi = plNode->Y();

	    double Zo = prNode->Z();
		double Yo = PBC_RotateY(prNode);

		bool RNodeIsLNode = ( abs(Zo - Zi) <= ( abs(Zo) * mPBCTolerance ) ) && 
			                ( abs(Yo - Yi) <= ( abs(Yo) * mPBCTolerance ) );

		if ( ( RNodeIsLNode == true ) || ( mPBCTolerance > 0.1 ) )
		{
			mPBC_NodeNodePairs_Right[prNode->Id()] = leftNodeId;
			return true;
		}
		else
		{
			return false;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Zo,Yo in each element of elemSet
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_RL_NearestPointInElementSet( Node::Pointer prNode, std::vector<std::vector<unsigned int> >& elemSet, 
												   std::vector<std::vector<double> >& minNatCoordSet,
												   std::vector<unsigned int>&         minEleNodesIds,
					                               std::vector<double>&               minNatCoord )
	{
		minNatCoordSet.clear(); minNatCoordSet.swap( std::vector<std::vector<double> > ( minNatCoordSet ) );
		minEleNodesIds.clear(); minEleNodesIds.swap( std::vector<unsigned int>         ( minEleNodesIds ) );
		minNatCoord.clear()   ; minNatCoord.swap   ( std::vector<double>               ( minNatCoord    ) );

		std::vector<unsigned int> iniNodesIds( *elemSet.begin() );

		for( int i = 0; i < iniNodesIds.size(); i++ ) 
		{ 
			minEleNodesIds.push_back( iniNodesIds[i] );
		}

		minNatCoord.push_back( 0.0 ); 
		minNatCoord.push_back( 0.0 );

		Node::Pointer iniNode = mpModel->GetNode( iniNodesIds[0] - 1 );

		double Zo = prNode->Z();
		double Yo = PBC_RotateY(prNode);

		double minDistToEle = abs( Zo - iniNode->Z() ) + abs( Yo - iniNode->Y() );
		double distToEle;

		double nx_ini = 0.0; double nx_end = 1.0; 
		double ny_ini = 0.0; double ny_end = 1.0; 
		double step   = 0.1;

		std::vector<std::vector<unsigned int> >::iterator it_eSet;

		for(it_eSet = elemSet.begin(); it_eSet != elemSet.end(); ++it_eSet) 
		{
			std::vector<unsigned int> eNodesIds(*it_eSet); 
			std::vector<double>       mdnCoord (2, 0.00);

			distToEle = PBC_RL_MinDistPointInEle( prNode, eNodesIds, nx_ini, nx_end, ny_ini, ny_end, step, mdnCoord );

			minNatCoordSet.push_back( mdnCoord );

			if ( distToEle <= minDistToEle )
			{
				minDistToEle = distToEle;				
				
				minNatCoord[0] = mdnCoord[0]; 
				minNatCoord[1] = mdnCoord[1];
				
				for( int i = 0; i < eNodesIds.size(); i++ ) minEleNodesIds[i] = eNodesIds[i];
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Zo,Yo in an element
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Modeler::PBC_RL_MinDistPointInEle( Node::Pointer prNode, std::vector<unsigned int>& NodesIds, 
		                                      double nx_ini, double nx_end, 
											  double ny_ini, double ny_end, 
											  double step, 
											  std::vector<double>& natCoord )
	{
		std::vector<Node::Pointer> peNodes;

		std::vector<unsigned int>::iterator it_nId;

		for(it_nId = NodesIds.begin(); it_nId != NodesIds.end(); ++it_nId) 
		{
			peNodes.push_back( mpModel->GetNode( *it_nId - 1 ) );
		}

		double Zo = prNode->Z();
		double Yo = PBC_RotateY(prNode);
	
		double nx, ny;
		double Zi, Yi;
		
		double dist;
		double min_dist = abs( Zo - peNodes[0]->Z() ) + abs ( Yo - peNodes[0]->Y() );

		for ( nx = nx_ini; nx <= nx_end; nx = nx + step ) 
		{ 
			for ( ny = ny_ini; ny <= ny_end; ny = ny + step ) 
			{ 
				if ( ( ( nx + ny ) > 1.0 ) || ( nx < 0.0 ) || ( ny < 0.0 ) ) continue;

				nCoordToZY(nx, ny, Zi, Yi, peNodes); 
				
				dist = abs( Zo - Zi ) + abs ( Yo - Yi );

				if ( dist <= min_dist )
				{
					min_dist    = dist;
					natCoord[0] = nx;
					natCoord[1] = ny;
				}
			}
		}

		return min_dist;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Transforms natural coordinates nx, ny of element pNodes in cartesian coordinates Z,Y
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::nCoordToZY (double nx, double ny, double& Z, double& Y, std::vector<Node::Pointer>& pNodes)
	{
		double L2 = nx; 
		double L3 = ny; 
		double L1 = 1.00 - L2 - L3;

		double N0 = L1 * (2.00*L1 - 1.00);
        double N1 = L2 * (2.00*L2 - 1.00);
        double N2 = L3 * (2.00*L3 - 1.00);
		double N3 = 4.00 * L1 * L2;
		double N4 = 4.00 * L2 * L3;
		double N5 = 4.00 * L3 * L1;

		Z = N0 * pNodes[0]->Z() + N1 * pNodes[1]->Z() + N2 * pNodes[2]->Z() +
			N3 * pNodes[3]->Z() + N4 * pNodes[4]->Z() + N5 * pNodes[5]->Z() ;

		Y = N0 * pNodes[0]->Y() + N1 * pNodes[1]->Y() + N2 * pNodes[2]->Y() +
			N3 * pNodes[3]->Y() + N4 * pNodes[4]->Y() + N5 * pNodes[5]->Y() ;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Zo,Yo in each element of elemSet when the point is in an edge
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_RL_SearchMinPointInEdges( Node::Pointer prNode, double prec, 
		                                        std::vector<std::vector<unsigned int> >& elemSet, 
												std::vector<std::vector<double> >& minNatCoordSet,
												std::vector<unsigned int>&         minEleNodesIds,
					                            std::vector<double>&               minNatCoord )
	{
		Node::Pointer iniNode = mpModel->GetNode( minEleNodesIds[0] - 1 );

		double Zo = prNode->Z();
		double Yo = PBC_RotateY(prNode);
		
		double minDistToEdge = abs( Zo - iniNode->X() ) + abs( Yo - iniNode->Y() );
        double distToEdge;

        double incInt = prec * 10.0;

		std::vector<std::vector<unsigned int> >::iterator it_eSet;

		std::vector<std::vector<double> >::iterator it_nC = minNatCoordSet.begin();

		for( it_eSet = elemSet.begin(); it_eSet != elemSet.end(); ++it_eSet, ++it_nC ) 
		{
			std::vector<unsigned int> eNodesIds( *it_eSet );
			std::vector<double>       minCoord ( *it_nC   );

            double nx_ini = minCoord[0] - incInt; 
			double nx_end = minCoord[0] + incInt; 

			double ny_ini = minCoord[1] - incInt; 
			double ny_end = minCoord[1] + incInt; 

			distToEdge = PBC_RL_MinDistPointInEle( prNode, eNodesIds, nx_ini, nx_end, ny_ini, ny_end, prec, minCoord );

			(*it_nC)[0] = minCoord[0]; 
			(*it_nC)[1] = minCoord[1];

			if ( distToEdge <= minDistToEdge )
			{
				minDistToEdge = distToEdge;				
				
				minNatCoord[0] = minCoord[0]; 
				minNatCoord[1] = minCoord[1];
				
				for( int i = 0; i < eNodesIds.size(); i++ ) minEleNodesIds[i] = eNodesIds[i];
			}
		}	
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Search for the nearest point to Zo,Yo in an element up to tolerance
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_RL_SearchMinPointInElementToTol( Node::Pointer prNode, double prec, std::vector<unsigned int>& minEleNodesIds, std::vector<double>& minNatCoord )
	{
		do 
		{
			double incInt = prec * 10.0;

			double nx_ini = minNatCoord[0] - incInt; 
			double nx_end = minNatCoord[0] + incInt; 

			double ny_ini = minNatCoord[1] - incInt; 
			double ny_end = minNatCoord[1] + incInt; 
				
			PBC_RL_MinDistPointInEle( prNode, minEleNodesIds, nx_ini, nx_end, ny_ini, ny_end, prec, minNatCoord );

			prec = prec * 0.1;	

		}while ( prec >= mPBCTolerance );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Set the point Zo,Yo in the element
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_RL_SetPointInElement( Node::Pointer prNode, std::vector<unsigned int>& lNodesIds, std::vector<double>& natCoord)
	{
		if      ( abs( natCoord[0] ) < mPBCTolerance )  
		{ 
			if      ( abs( natCoord[1]       ) < mPBCTolerance ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[0];
			else if ( abs( natCoord[1] - 0.5 ) < mPBCTolerance ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[5];
			else if ( abs( natCoord[1] - 1.0 ) < mPBCTolerance ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[2];
			else                                                 PBC_RL_SetPointInEdge( prNode, lNodesIds, 0.0, natCoord[1] );	
		}
		else if ( abs( natCoord[1] ) < mPBCTolerance )  
		{ 
			if      ( abs( natCoord[0] - 0.5 ) < mPBCTolerance ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[3];
			else if ( abs( natCoord[0] - 1.0 ) < mPBCTolerance ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[1];
			else                                                 PBC_RL_SetPointInEdge( prNode, lNodesIds, natCoord[0], 0.0 );
		}
		else if ( abs( natCoord[0] + natCoord[1] - 1.0 ) < mPBCTolerance )  
		{ 
			if ( ( abs( natCoord[0] - 0.5 ) < mPBCTolerance ) && 
				 ( abs( natCoord[1] - 0.5 ) < mPBCTolerance )    ) mPBC_NodeNodePairs_Right[prNode->Id()] = lNodesIds[4];
			else                                                   PBC_RL_SetPointInEdge( prNode, lNodesIds, natCoord[0], natCoord[1] );
		}
		else
		{
			mPBC_NodeElementPairs_Right[prNode->Id()] = lNodesIds;
			mPBC_NodeElementCoord_Right[prNode->Id()] = natCoord;
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Set the point of natural coordinates nX,nY in the edge of the element 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::PBC_RL_SetPointInEdge( Node::Pointer prNode, std::vector<unsigned int>& lNodesIds, double nx, double ny)
	{
        int    n1,n2,n3;
		double edgeCoord; 

		if      ( nx == 0.0 ) { n1 = 2; n2 = 0; n3 = 5; edgeCoord = 1 - ny; }
		else if ( ny == 0.0 ) { n1 = 0; n2 = 1; n3 = 3; edgeCoord = nx    ; }
		else                  { n1 = 1; n2 = 2; n3 = 4; edgeCoord = 1 - nx; }

		std::vector<unsigned int> edgeNodesIds(3);
		edgeNodesIds[0] = lNodesIds[n1];
		edgeNodesIds[1] = lNodesIds[n2];
		edgeNodesIds[2] = lNodesIds[n3];

		mPBC_NodeEdgePairs_Right[prNode->Id()] = edgeNodesIds;
        mPBC_NodeEdgeCoord_Right[prNode->Id()] = edgeCoord;
	}

	//***********************************************************************************
    //* - Calculates the volume integral of the fields.
    //***********************************************************************************
    void Modeler::GenerateVolIntElement(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId, HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double volumeID = cIBC[0];

		int gp;

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order6(cX, cY, cZ, W);

		Matrix<double> N;

		if      ( mElementOrder == 1 ) Lagrange3D_Ni_1st(N, cX, cY, cZ);
        else if ( mElementOrder == 2 ) Lagrange3D_Ni_2nd(N, cX, cY, cZ);
        else if ( mElementOrder == 3 ) Lagrange3D_Ni_3th(N, cX, cY, cZ);
		else if ( mElementOrder == 4 ) Lagrange3D_Ni_4th(N, cX, cY, cZ);

		double elementVolume = 0.00;
		
		if ( mQuadraticGeometry )
		{
			std::vector<double> detJ;
			Calculate_detJ(detJ, cX, cY, cZ, nodes);

			for ( gp=0; gp<nGaussPoints; gp++ ) W[gp]         *= detJ[gp]; 
			for ( gp=0; gp<nGaussPoints; gp++ ) elementVolume +=    W[gp]; 
		}
		else
		{
			elementVolume = VolumeOfTetrahedra(nodes);
			double jacob  = 6.00*elementVolume;

			for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob; 
		}

		std::vector<std::complex<double> > volumeIntegralE(5);

		volumeIntegralE[0] = std::complex<double>(0.00,0.00);
		volumeIntegralE[1] = std::complex<double>(0.00,0.00);
		volumeIntegralE[2] = std::complex<double>(0.00,0.00);
		volumeIntegralE[3] = std::complex<double>(0.00,0.00);
		volumeIntegralE[4] = std::complex<double>(0.00,0.00);

		std::vector<std::complex<double> > volumeIntegralH(5);

		volumeIntegralH[0] = std::complex<double>(0.00,0.00);
		volumeIntegralH[1] = std::complex<double>(0.00,0.00);
		volumeIntegralH[2] = std::complex<double>(0.00,0.00);
		volumeIntegralH[3] = std::complex<double>(0.00,0.00);
		volumeIntegralH[4] = std::complex<double>(0.00,0.00);

		std::vector<std::complex<double> > volumeIntegralJ(6);

		volumeIntegralJ[0] = std::complex<double>(0.00,0.00);
		volumeIntegralJ[1] = std::complex<double>(0.00,0.00);
		volumeIntegralJ[2] = std::complex<double>(0.00,0.00);
		volumeIntegralJ[3] = std::complex<double>(0.00,0.00);
		volumeIntegralJ[4] = std::complex<double>(0.00,0.00);
		volumeIntegralJ[5] = std::complex<double>(0.00,0.00);

        std::vector<Node::Pointer>::iterator node_it;

		int in = 0;

        for (node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
		    double vIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ ) vIntNi += W[gp]*N[in][gp]; 
			
			// Electric field
			std::complex<double> cExi = (*properties)(cEx, **node_it);
			std::complex<double> cEyi = (*properties)(cEy, **node_it);
			std::complex<double> cEzi = (*properties)(cEz, **node_it);

			double modE2 = std::real(cExi)*std::real(cExi) + std::imag(cExi)*std::imag(cExi) +
                           std::real(cEyi)*std::real(cEyi) + std::imag(cEyi)*std::imag(cEyi) +
                           std::real(cEzi)*std::real(cEzi) + std::imag(cEzi)*std::imag(cEzi) ;

			double modE  = sqrt(modE2);

			volumeIntegralE[0] += cExi*vIntNi;
			volumeIntegralE[1] += cEyi*vIntNi;
			volumeIntegralE[2] += cEzi*vIntNi;

            volumeIntegralE[3] += modE2 * vIntNi;
			volumeIntegralE[4] += modE  * vIntNi;

            // Magnetic field
			std::vector<double> H_real = (*mpModel)(REAL_H, **node_it); 
            std::vector<double> H_imag = (*mpModel)(IMAG_H, **node_it);

			std::complex<double> cHxi (H_real[0], H_imag[0]);
			std::complex<double> cHyi (H_real[1], H_imag[1]);
			std::complex<double> cHzi (H_real[2], H_imag[2]);

			double modH2 = H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                           H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                           H_real[2]*H_real[2] + H_imag[2]*H_imag[2] ;

			double modH  = sqrt(modH2);

			volumeIntegralH[0] += cHxi*vIntNi;
			volumeIntegralH[1] += cHyi*vIntNi;
			volumeIntegralH[2] += cHzi*vIntNi;

            volumeIntegralH[3] += modH2 * vIntNi;
			volumeIntegralH[4] += modH  * vIntNi;

			// Current density
			std::vector<double> J_real = (*mpModel)(REAL_J, **node_it); 
            std::vector<double> J_imag = (*mpModel)(IMAG_J, **node_it);

			std::complex<double> cJxi (J_real[0], J_imag[0]);
			std::complex<double> cJyi (J_real[1], J_imag[1]);
			std::complex<double> cJzi (J_real[2], J_imag[2]);

			double modJ2 = J_real[0]*J_real[0] + J_imag[0]*J_imag[0] +
                           J_real[1]*J_real[1] + J_imag[1]*J_imag[1] +
                           J_real[2]*J_real[2] + J_imag[2]*J_imag[2] ;

            double modJ  = sqrt(modJ2);

            std::complex<double> conjJpE = std::conj(cJxi)*cExi + 
				                           std::conj(cJyi)*cEyi + 
										   std::conj(cJzi)*cEzi;

			volumeIntegralJ[0] += cJxi*vIntNi;
			volumeIntegralJ[1] += cJyi*vIntNi;
			volumeIntegralJ[2] += cJzi*vIntNi;

            volumeIntegralJ[3] += modJ2 * vIntNi;
			volumeIntegralJ[4] += modJ  * vIntNi;

			volumeIntegralJ[5] += conjJpE * vIntNi;

			in++;
        }

		if ( mVolumeIntegralE.find(volumeID) != mVolumeIntegralE.end() )
		{
			mVolumeValue    [volumeID] += elementVolume;

			mVolumeIntegralE[volumeID] += volumeIntegralE;
			mVolumeIntegralH[volumeID] += volumeIntegralH;
			mVolumeIntegralJ[volumeID] += volumeIntegralJ;
		}
        else
		{
			mVolumeValue    [volumeID] = elementVolume;

			mVolumeIntegralE[volumeID] = volumeIntegralE;
			mVolumeIntegralH[volumeID] = volumeIntegralH;
			mVolumeIntegralJ[volumeID] = volumeIntegralJ;
		}
	}

   /***********************************************************************************
    * - Calculation of dN/da, dN/db, dN/du on cX,cY,cZ points
    ***********************************************************************************/
    void Modeler::NaturalDerivatives3D_2ndQ(Matrix<double>& dNda,
		                                    Matrix<double>& dNdb,
											Matrix<double>& dNdu,
											std::vector<double>& cX, 
											std::vector<double>& cY, 
											std::vector<double>& cZ)
    {
		int p;
		int numNodes  = 10;
		int numPoints = cX.size();

		dNda.Resize(numNodes, numPoints);
		dNdb.Resize(numNodes, numPoints);
		dNdu.Resize(numNodes, numPoints);

		double L1,L2,L3,L4;

		//dNi/da ,dNi/db, dNi/du
        for ( p=0; p<numPoints; p++ ) 
        {
            L2 = cX[p];
            L3 = cY[p];
			L4 = cZ[p];
            L1 = 1.00 - L2 - L3 - L4;

			dNda[0][p] =  1.00 - 4.00 * L1;
			dNdb[0][p] =  1.00 - 4.00 * L1;
			dNdu[0][p] =  1.00 - 4.00 * L1;

			dNda[1][p] =  4.00 * L2 - 1.00;
			dNdb[1][p] =  0.00;
			dNdu[1][p] =  0.00;

			dNda[2][p] =  0.00;
			dNdb[2][p] =  4.00 * L3 - 1.00;
			dNdu[2][p] =  0.00;

			dNda[3][p] =  0.00;
			dNdb[3][p] =  0.00;
			dNdu[3][p] =  4.00 * L4 - 1.00;

			dNda[4][p] = -4.00 * L2 + 4.00 * L1;
			dNdb[4][p] = -4.00 * L2;
			dNdu[4][p] = -4.00 * L2;

			dNda[5][p] =  4.00 * L3;
			dNdb[5][p] =  4.00 * L2;
			dNdu[5][p] =  0.00;

			dNda[6][p] = -4.00 * L3 ;
			dNdb[6][p] = -4.00 * L3 + 4.00 * L1;
			dNdu[6][p] = -4.00 * L3 ;

			dNda[7][p] = -4.00 * L4;
			dNdb[7][p] = -4.00 * L4;
			dNdu[7][p] = -4.00 * L4 + 4.00 * L1;

			dNda[8][p] =  4.00 * L4;
			dNdb[8][p] =  0.00;
			dNdu[8][p] =  4.00 * L2;

			dNda[9][p] =  0.00;
			dNdb[9][p] =  4.00 * L4;
			dNdu[9][p] =  4.00 * L3;
		}
	}

   /***********************************************************************************
    * - Calculation of the Jacobian determinant
    ***********************************************************************************/
    void Modeler::Calculate_detJ(std::vector<double>& detJ,
		                         std::vector<double>& cX, 
								 std::vector<double>& cY, 
								 std::vector<double>& cZ,
								 std::vector<Node::Pointer>& cnodes)
    {
		int i, p;
		int numNodes  = 10;
		int numPoints = cX.size();

		detJ.resize(numPoints,0.00);

		Matrix<double> dNda, dNdb, dNdu;

		NaturalDerivatives3D_2ndQ(dNda, dNdb, dNdu, cX, cY, cZ);

		for ( p=0; p<numPoints; p++ ) 
        {
			double dXda = 0.00;	double dXdb = 0.00; double dXdu = 0.00;
			double dYda = 0.00;	double dYdb = 0.00; double dYdu = 0.00;
			double dZda = 0.00;	double dZdb = 0.00; double dZdu = 0.00;
			
			for ( i=0; i<numNodes; i++ ) 
			{
				dXda += dNda[i][p] * cnodes[i]->X();  dXdb += dNdb[i][p] * cnodes[i]->X();  dXdu += dNdu[i][p] * cnodes[i]->X();
				dYda += dNda[i][p] * cnodes[i]->Y();  dYdb += dNdb[i][p] * cnodes[i]->Y();  dYdu += dNdu[i][p] * cnodes[i]->Y();
				dZda += dNda[i][p] * cnodes[i]->Z();  dZdb += dNdb[i][p] * cnodes[i]->Z();  dZdu += dNdu[i][p] * cnodes[i]->Z();
			}
	
			detJ[p] = dXda*dYdb*dZdu + dXdb*dZda*dYdu + dXdu*dYda*dZdb - dXda*dZdb*dYdu - dXdb*dYda*dZdu - dXdu*dZda*dYdb;
		}
	}

    //***************************************************************************************
    //* - Calculate volume of a linear tetrahedra
    //***************************************************************************************
    double Modeler::VolumeOfTetrahedra(std::vector<Node::Pointer>& cnodes)
    {
        double X1 = cnodes[0]->X(), Y1 = cnodes[0]->Y(), Z1 = cnodes[0]->Z();
		double X2 = cnodes[1]->X(), Y2 = cnodes[1]->Y(), Z2 = cnodes[1]->Z();
		double X3 = cnodes[2]->X(), Y3 = cnodes[2]->Y(), Z3 = cnodes[2]->Z();
		double X4 = cnodes[3]->X(), Y4 = cnodes[3]->Y(), Z4 = cnodes[3]->Z();

		double det;

        det = + X2*Y3*Z4 + X4*Y2*Z3 + X3*Y4*Z2 - X4*Y3*Z2 - X2*Y4*Z3 - X3*Y2*Z4 
              - X1*Y3*Z4 - X4*Y1*Z3 - X3*Y4*Z1 + X4*Y3*Z1 + X1*Y4*Z3 + X3*Y1*Z4 
              + X1*Y2*Z4 + X4*Y1*Z2 + X2*Y4*Z1 - X4*Y2*Z1 - X1*Y4*Z2 - X2*Y1*Z4 
              - X1*Y2*Z3 - X3*Y1*Z2 - X2*Y3*Z1 + X3*Y2*Z1 + X1*Y3*Z2 + X2*Y1*Z3;
         
        return fabs(det/6.00);
    }

	//***************************************************************************************
    //* - Calculate surface of a linear tetrahedra
    //***************************************************************************************
    double Modeler::AreaOfTriangle(std::vector<Node::Pointer>& cnodes)
    {
		double n[3],v1[3],v2[3];

        v2[0] = cnodes[2]->X() - cnodes[0]->X();
        v2[1] = cnodes[2]->Y() - cnodes[0]->Y();
        v2[2] = cnodes[2]->Z() - cnodes[0]->Z();

        v1[0] = cnodes[1]->X() - cnodes[0]->X();
        v1[1] = cnodes[1]->Y() - cnodes[0]->Y();
        v1[2] = cnodes[1]->Z() - cnodes[0]->Z();

        //area = 0.5*(v2 x v1)
        n[0] = v2[1]*v1[2] - v2[2]*v1[1];
        n[1] = v2[2]*v1[0] - v2[0]*v1[2];
        n[2] = v2[0]*v1[1] - v2[1]*v1[0];

		return ( 0.5*sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] ) );
    }

   /***********************************************************************************
    * - Calculation of the Jacobian determinant
    ***********************************************************************************/
    void Modeler::Calculate_detJ(std::vector<double>& detJ,
		                         std::vector<double>& cX, 
							     std::vector<double>& cY,
								 std::vector<Node::Pointer>& cnodes)
    {
		int i, p;
		int numNodes  = 6;
		int numPoints = cX.size();

		std::vector<double> dNda(numNodes);
		std::vector<double> dNdb(numNodes);

		detJ.resize(numPoints);

        for (p=0; p<numPoints; p++) 
        {
			// dNi/da ,dNi/db
			dNda[0] =   4.00*cX[p] + 4.00*cY[p] - 3.00;
			dNdb[0] =   4.00*cX[p] + 4.00*cY[p] - 3.00;

			dNda[1] =   4.00*cX[p] - 1.00;
			dNdb[1] =   0.00;

			dNda[2] =   0.00;
			dNdb[2] =   4.00*cY[p] - 1.00;

			dNda[3] = - 8.00*cX[p] - 4.00*cY[p] + 4.00;
			dNdb[3] = - 4.00*cX[p];

			dNda[4] =   4.00*cY[p];
			dNdb[4] =   4.00*cX[p];

			dNda[5] = - 4.00*cY[p];
			dNdb[5] = - 8.00*cY[p] - 4.00*cX[p] + 4.00;

			double dXda = 0.00;	double dXdb = 0.00;
			double dYda = 0.00;	double dYdb = 0.00;
			double dZda = 0.00;	double dZdb = 0.00;

			for (i=0; i<numNodes; i++) 
			{
				dXda += dNda[i] * cnodes[i]->X();  dXdb += dNdb[i] * cnodes[i]->X();
				dYda += dNda[i] * cnodes[i]->Y();  dYdb += dNdb[i] * cnodes[i]->Y();         
				dZda += dNda[i] * cnodes[i]->Z();  dZdb += dNdb[i] * cnodes[i]->Z();
			}

			Vector<double> normal(3);

			normal[0] = dZda*dYdb - dYda*dZdb;
			normal[1] = dXda*dZdb - dZda*dXdb;
			normal[2] = dYda*dXdb - dXda*dYdb;

			detJ[p] = sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00));
		}
	}

	//***********************************************************************************
    //* - Calculates the surface integral of the fields.
    //***********************************************************************************
    void Modeler::GenerateSurfIntElement(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId, HONodesId);

		std::vector<int>::iterator hoit;
        
        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double surfaceID = cIBC[0];

		int gp;

		std::vector<double> cX; 
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order8(cX, cY, W);

		Matrix<double> N;

		if      ( mElementOrder == 1 ) Lagrange2D_Ni_1st(N, cX, cY);
        else if ( mElementOrder == 2 ) Lagrange2D_Ni_2nd(N, cX, cY);
        else if ( mElementOrder == 3 ) Lagrange2D_Ni_3th(N, cX, cY);
		else if ( mElementOrder == 4 ) Lagrange2D_Ni_4th(N, cX, cY);

		double elementArea = 0.00;
		
		if ( mQuadraticGeometry )
		{
			std::vector<double> detJ;
			Calculate_detJ(detJ, cX, cY, nodes);

			for ( gp=0; gp<nGaussPoints; gp++ ) W[gp]       *= detJ[gp]; 
			for ( gp=0; gp<nGaussPoints; gp++ ) elementArea += W[gp]; 
		}
		else
		{
			elementArea  = AreaOfTriangle(nodes);
			double jacob = 2.00*elementArea;

			for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob; 
		}

		std::vector<std::complex<double> > surfaceIntegralE(5);

		surfaceIntegralE[0] = std::complex<double>(0.00,0.00);
		surfaceIntegralE[1] = std::complex<double>(0.00,0.00);
		surfaceIntegralE[2] = std::complex<double>(0.00,0.00);
		surfaceIntegralE[3] = std::complex<double>(0.00,0.00);
		surfaceIntegralE[4] = std::complex<double>(0.00,0.00);

		std::vector<std::complex<double> > surfaceIntegralH(5);

		surfaceIntegralH[0] = std::complex<double>(0.00,0.00);
		surfaceIntegralH[1] = std::complex<double>(0.00,0.00);
		surfaceIntegralH[2] = std::complex<double>(0.00,0.00);
		surfaceIntegralH[3] = std::complex<double>(0.00,0.00);
		surfaceIntegralH[4] = std::complex<double>(0.00,0.00);

		std::vector<std::complex<double> > surfaceIntegralJ(5);

		surfaceIntegralJ[0] = std::complex<double>(0.00,0.00);
		surfaceIntegralJ[1] = std::complex<double>(0.00,0.00);
		surfaceIntegralJ[2] = std::complex<double>(0.00,0.00);
		surfaceIntegralJ[3] = std::complex<double>(0.00,0.00);
		surfaceIntegralJ[4] = std::complex<double>(0.00,0.00);

        std::vector<Node::Pointer>::iterator node_it;

		int in = 0;

        for ( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
		    double sIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ ) sIntNi += W[gp]*N[in][gp]; 
			
			// Electric field
			std::complex<double> cExi = (*properties)(cEx, **node_it);
			std::complex<double> cEyi = (*properties)(cEy, **node_it);
			std::complex<double> cEzi = (*properties)(cEz, **node_it);

			surfaceIntegralE[0] += cExi*sIntNi;
			surfaceIntegralE[1] += cEyi*sIntNi;
			surfaceIntegralE[2] += cEzi*sIntNi;

			double modE2 = std::real(cExi)*std::real(cExi) + std::imag(cExi)*std::imag(cExi) +
                           std::real(cEyi)*std::real(cEyi) + std::imag(cEyi)*std::imag(cEyi) +
                           std::real(cEzi)*std::real(cEzi) + std::imag(cEzi)*std::imag(cEzi) ;

			double modE  = sqrt(modE2);

			surfaceIntegralE[3] += modE2 * sIntNi;
			surfaceIntegralE[4] += modE  * sIntNi;

            // Magnetic field
			std::vector<double> H_real = (*mpModel)(REAL_H, **node_it); 
            std::vector<double> H_imag = (*mpModel)(IMAG_H, **node_it);

			std::complex<double> cHxi (H_real[0], H_imag[0]);
			std::complex<double> cHyi (H_real[1], H_imag[1]);
			std::complex<double> cHzi (H_real[2], H_imag[2]);

			surfaceIntegralH[0] += cHxi*sIntNi;
			surfaceIntegralH[1] += cHyi*sIntNi;
			surfaceIntegralH[2] += cHzi*sIntNi;

			double modH2 = H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                           H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                           H_real[2]*H_real[2] + H_imag[2]*H_imag[2] ;
			
			double modH  = sqrt(modH2);

			surfaceIntegralH[3] += modH2 * sIntNi;
			surfaceIntegralH[4] += modH  * sIntNi;

			// Current density
			std::vector<double> J_real = (*mpModel)(REAL_J, **node_it); 
            std::vector<double> J_imag = (*mpModel)(IMAG_J, **node_it);

			std::complex<double> cJxi (J_real[0], J_imag[0]);
			std::complex<double> cJyi (J_real[1], J_imag[1]);
			std::complex<double> cJzi (J_real[2], J_imag[2]);

			surfaceIntegralJ[0] += cJxi*sIntNi;
			surfaceIntegralJ[1] += cJyi*sIntNi;
			surfaceIntegralJ[2] += cJzi*sIntNi;

			double modJ2 = J_real[0]*J_real[0] + J_imag[0]*J_imag[0] +
                           J_real[1]*J_real[1] + J_imag[1]*J_imag[1] +
                           J_real[2]*J_real[2] + J_imag[2]*J_imag[2] ;

			double modJ  = sqrt(modJ2);

			surfaceIntegralJ[3] += modJ2 * sIntNi;
            surfaceIntegralJ[4] += modJ  * sIntNi;

			in++;
        }

		if ( mSurfaceIntegralE.find(surfaceID) != mSurfaceIntegralE.end() )
		{
			mSurfaceValue    [surfaceID] += elementArea;
			mSurfaceIntegralE[surfaceID] += surfaceIntegralE;
			mSurfaceIntegralH[surfaceID] += surfaceIntegralH;
			mSurfaceIntegralJ[surfaceID] += surfaceIntegralJ;
		}
        else
		{
			mSurfaceValue    [surfaceID] = elementArea;
			mSurfaceIntegralE[surfaceID] = surfaceIntegralE;
			mSurfaceIntegralH[surfaceID] = surfaceIntegralH;
			mSurfaceIntegralJ[surfaceID] = surfaceIntegralJ;
		}
	}
    
	//***********************************************************************************
    //* - Projecting resulting field on plane
    //***********************************************************************************
    void Modeler::GenerateProjectionRWPortTE10(int*  NodesId, unsigned int PropertiesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if      ( mElementOrder == 1 ) pElement = Element::Pointer(new RWPortTE10_1st(nodes,properties));
        else if ( mElementOrder == 2 ) pElement = Element::Pointer(new RWPortTE10_2nd(nodes,properties));
        else if ( mElementOrder == 3 ) pElement = Element::Pointer(new RWPortTE10_3th(nodes,properties));
		else if ( mElementOrder == 4 ) pElement = Element::Pointer(new RWPortTE10_4th(nodes,properties));

        double portNumber = pElement->getPortNumber();
		
		std::complex<double> TE10Proyec    = pElement->getProjection();
		std::complex<double> Normalization = pElement->getNormalization();
		
		//Elemental projection
		if (mProjection.find(portNumber) != mProjection.end()) mProjection[portNumber] += TE10Proyec;
        else                                                   mProjection[portNumber]  = TE10Proyec;

		//Elemental normalization
		if (mNormalization.find(portNumber) != mNormalization.end()) mNormalization[portNumber] += Normalization;
        else                                                         mNormalization[portNumber]  = Normalization;
    }

	//***********************************************************************************
    //* - Projecting resulting field on plane
    //***********************************************************************************
    void Modeler::GenerateProjectionCoaxPortTEM(int*  NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;
        
        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if ( mQuadraticGeometry )
		{
			pElement = Element::Pointer(new CoaxialPortTEM_2ndQ(nodes,properties));
		}
		else
		{
			if      ( mElementOrder == 1 ) pElement = Element::Pointer(new CoaxialPortTEM_1st(nodes,properties));
			else if ( mElementOrder == 2 ) pElement = Element::Pointer(new CoaxialPortTEM_2nd(nodes,properties));
			else if ( mElementOrder == 3 ) pElement = Element::Pointer(new CoaxialPortTEM_3th(nodes,properties));
			else if ( mElementOrder == 4 ) pElement = Element::Pointer(new CoaxialPortTEM_4th(nodes,properties));
		}

        double portNumber = pElement->getPortNumber();
		
		std::complex<double> TEMProyec     = pElement->getProjection();
		std::complex<double> Normalization = pElement->getNormalization();
		
		//Elemental projection
		if (mProjection.find(portNumber) != mProjection.end()) mProjection[portNumber] += TEMProyec;
        else                                                   mProjection[portNumber]  = TEMProyec;

		//Elemental normalization
		if (mNormalization.find(portNumber) != mNormalization.end()) mNormalization[portNumber] += Normalization;
        else                                                         mNormalization[portNumber]  = Normalization; 
    }

    //************************************************************************************************* 
    //* - Calculates average normal in PEC nodes
    //*************************************************************************************************
    void Modeler::GeneratePECElement(std::vector<int>& LONodesId)
    {
         if ( mQuadraticGeometry )
		 {  
			 GenerateQPECElement(LONodesId);
			 return;
		 }

		 Node::Pointer pNode;
         std::vector<Node::Pointer> nodes;
         std::vector<int>::const_iterator id_it;

         std::vector<int> NodesId;

		 PushHONodesOnSurface(LONodesId, NodesId);

         for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
         {
             pNode = mpModel->GetNode((*id_it)-1);
             nodes.push_back(pNode);
         }

         // Exterior normal of the element  
         Vector<double> normal(3);

         // Normal = (v2 x v1)/mod_normal
         normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
                     (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
         normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
                     (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
         normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
                     (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

		 // Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				 mSetPECnormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				mSetPECnormals[*id_it].push_back(normal);
			}
		}
	}

	//************************************************************************************************* 
    //* - Calculates average normal in PEC nodes for quadratic meshing
    //*************************************************************************************************
    void Modeler::GenerateQPECElement(std::vector<int>& LONodesId)
	{
		std::vector<int> NodesId;
		std::vector<int>::const_iterator id_it;

		PushHONodesOnSurface(LONodesId, NodesId);

		int numNodes = 6;
		int nN       = 0;
		int pt       = 0;

		std::vector<double> Xi(numNodes), Yi(numNodes), Zi(numNodes);

		for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
		{
			Node::Pointer pNode = mpModel->GetNode((*id_it)-1);

			Xi[nN] = pNode->X();
			Yi[nN] = pNode->Y();
			Zi[nN] = pNode->Z();

			nN++;
		}

		Matrix<double> dNda, dNdb;

		SurfDerivativesInNodes_2ndQ(dNda, dNdb);

		for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
		{
			// Exterior normal of the element  
		    Vector<double> normal(3);

			double dXda = 0.00;	double dXdb = 0.00;
			double dYda = 0.00;	double dYdb = 0.00;
			double dZda = 0.00;	double dZdb = 0.00;

			//dNi/da ,dNi/db
			for (nN=0; nN<numNodes; nN++) 
			{
				dXda += dNda[nN][pt]*Xi[nN];  dXdb += dNdb[nN][pt]*Xi[nN];
				dYda += dNda[nN][pt]*Yi[nN];  dYdb += dNdb[nN][pt]*Yi[nN];         
				dZda += dNda[nN][pt]*Zi[nN];  dZdb += dNdb[nN][pt]*Zi[nN];
			}

			normal[0] = dZda*dYdb - dYda*dZdb;
			normal[1] = dXda*dZdb - dZda*dXdb;
			normal[2] = dYda*dXdb - dXda*dYdb;

			if ( mNormalsAreaWtd )
			{
				mSetPECnormals[*id_it].push_back(normal);
			}
			else
			{
				normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));
				mSetPECnormals[*id_it].push_back(normal);
			}

			pt++;

		}//end for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
	}

   /*******************************************************************
    * - Calculate derivatives in nodes for normals calculations
    ********************************************************************/
    void Modeler::SurfDerivativesInNodes_2ndQ(Matrix<double>& dNda, 
		                                      Matrix<double>& dNdb)
    {
		const int numNodes = 6;

		std::vector<double> alfa(numNodes); 
		std::vector<double>	beta(numNodes);

		//Points 0,1,2
		alfa[0] = 0.00;  beta[0] = 0.00;
		alfa[1] = 1.00;  beta[1] = 0.00;
		alfa[2] = 0.00;  beta[2] = 1.00;
		//Edge [0 1]
		alfa[3] = 0.50;  beta[3] = 0.00;
		//Edge [1 2]
		alfa[4] = 0.50;  beta[4] = 0.50;
		//Edge [2 0]
		alfa[5] = 0.00;  beta[5] = 0.50;
		
		dNda.Resize(numNodes, numNodes);
		dNdb.Resize(numNodes, numNodes);

		//dNi/da ,dNi/db
        for (int p=0; p<numNodes; p++) 
        {
			dNda[0][p] =   4.00*alfa[p] + 4.00*beta[p] - 3.00;
			dNdb[0][p] =   4.00*alfa[p] + 4.00*beta[p] - 3.00;

			dNda[1][p] =   4.00*alfa[p] - 1.00;
			dNdb[1][p] =   0.00;

			dNda[2][p] =   0.00;
			dNdb[2][p] =   4.00*beta[p] - 1.00;

			dNda[3][p] = - 8.00*alfa[p] - 4.00*beta[p] + 4.00;
			dNdb[3][p] = - 4.00*alfa[p];

			dNda[4][p] =   4.00*beta[p];
			dNdb[4][p] =   4.00*alfa[p];

			dNda[5][p] = - 4.00*beta[p];
			dNdb[5][p] = - 8.00*beta[p] - 4.00*alfa[p] + 4.00;
		}
	}

   //************************************************************************************************* 
   //* - Calculates average normal in PMC nodes
   //*************************************************************************************************
   void Modeler::GeneratePMCElement(std::vector<int>& LONodesId)
   {
	     Node::Pointer pNode;
         std::vector<Node::Pointer> nodes;
         std::vector<int>::const_iterator id_it;

		 std::vector<int> NodesId;

		 PushHONodesOnSurface(LONodesId, NodesId);

         for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
         {
             pNode = mpModel->GetNode((*id_it)-1);
             nodes.push_back(pNode);
         }

         // Exterior normal of the element  
         Vector<double> normal(3,0.00);

         // Normal = (v2 x v1)/mod_normal
         normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
                     (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
         normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
                     (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
         normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
                     (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

		 // Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				 mSetPMCnormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				mSetPMCnormals[*id_it].push_back(normal);
			}
		}
    }

    //************************************************************************************************* 
    //* - Calculates average normal in symmetric PEC nodes
    //*************************************************************************************************
    void Modeler::GenerateSPECElement(std::vector<int>& LONodesId)
    {
		 Node::Pointer pNode;
         std::vector<Node::Pointer> nodes;
         std::vector<int>::const_iterator id_it;

         std::vector<int> NodesId;

		 PushHONodesOnSurface(LONodesId, NodesId);

         for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
         {
             pNode = mpModel->GetNode((*id_it)-1);
             nodes.push_back(pNode);
         }

         // Exterior normal of the element  
         Vector<double> normal(3);

         // Normal = (v2 x v1)/mod_normal
         normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
                     (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
         normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
                     (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
         normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
                     (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

		// Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				 mSetSPECnormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				mSetSPECnormals[*id_it].push_back(normal);
			}
		}
	}

	//************************************************************************************************* 
    //* - Calculates average normal in TE PMC nodes
    //*************************************************************************************************
    void Modeler::GenerateTEPMCElement(std::vector<int>& LONodesId)
    {
		 Node::Pointer pNode;
         std::vector<Node::Pointer> nodes;
         std::vector<int>::const_iterator id_it;

		 std::vector<int> NodesId;

		 PushHONodesOnSurface(LONodesId, NodesId);

         for(id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it) 
         {
             pNode = mpModel->GetNode((*id_it)-1);
             nodes.push_back(pNode);
         }

         // Exterior normal of the element  
         Vector<double> normal(3,0.00);

         // Normal = (v2 x v1)/mod_normal
         normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) - 
                     (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());
         normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) - 
                     (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());
         normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) - 
                     (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

		 // Area weighted normal
		if ( mNormalsAreaWtd )
		{
			normal *= (1.00/6.00); 

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				 mSetTEPMCnormals[*id_it].push_back(normal);
			}
		}
		// Geometric average normal
		else 
		{
			normal *= (1.00/sqrt(std::inner_product(normal.begin(),normal.end(),normal.begin(),0.00)));

			for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it) 
			{
				mSetTEPMCnormals[*id_it].push_back(normal);
			}
		}	
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	// - Setting, fixing, solving, updating, steps.
    //
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //************************************************************************************************* 
    //* - Finishing the first step
    //*************************************************************************************************
	void Modeler::FinishFirstStep()
	{	
		if ( mProjection.size() > 0 )
        {
            CalculateSijParameters();

	        PrintSijParameters (mProjectionNormalized, mProblemFrequency);
            WriteSijHeaderFiles(mProjectionNormalized);
            WriteSijInFiles    (mProjectionNormalized, mProblemFrequency);

		    mProjectionNormalized.clear();
            mProjectionNormalized.swap(std::map<double, std::complex<double> >(mProjectionNormalized));
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
	        PrintViParameters (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ, 
				               mVolumeValue    , mProblemFrequency);

            WriteViHeaderFiles(mVolumeIntegralE);

            WriteViInFiles    (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ, 
						       mProblemFrequency);

		    mVolumeIntegralE.clear();
            mVolumeIntegralE.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralE));

			mVolumeIntegralH.clear();
            mVolumeIntegralH.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralH));

			mVolumeIntegralJ.clear();
            mVolumeIntegralJ.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralJ));

			mVolumeValue.clear();
            mVolumeValue.swap(std::map<double, double>(mVolumeValue));
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters (mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ, 
				                  mSurfaceValue    , mProblemFrequency);

            WriteSurfiHeaderFiles(mSurfaceIntegralE);

            WriteSurfiInFiles    (mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ, 
				                  mProblemFrequency);

		    mSurfaceIntegralE.clear();
            mSurfaceIntegralE.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralE));

			mSurfaceIntegralH.clear();
            mSurfaceIntegralH.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralH));

			mSurfaceIntegralJ.clear();
            mSurfaceIntegralJ.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralJ));

			mSurfaceValue.clear();
            mSurfaceValue.swap(std::map<double, double>(mSurfaceValue));
        }

		//x_vector.resize(mSystemSize,0.00); 
        b_vector.resize(mSystemSize,0.00); 
        A_matrix.Resize(mSystemSize); 
	}

    //************************************************************************************************* 
    //* - Preparing the next step in frequency sweep mode
    //*************************************************************************************************
	void Modeler::FinishStep()
	{
		if ( mProjection.size() > 0 )
        {
            CalculateSijParameters();

	        PrintSijParameters(mProjectionNormalized, mProblemFrequency);
            WriteSijInFiles   (mProjectionNormalized, mProblemFrequency);

		    mProjectionNormalized.clear();
            mProjectionNormalized.swap(std::map<double, std::complex<double> >(mProjectionNormalized));
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
            PrintViParameters(mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ, 
				              mVolumeValue    , mProblemFrequency);

            WriteViInFiles   (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ, 
							  mProblemFrequency);

		    mVolumeIntegralE.clear();
            mVolumeIntegralE.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralE));

			mVolumeIntegralH.clear();
            mVolumeIntegralH.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralH));

			mVolumeIntegralJ.clear();
            mVolumeIntegralJ.swap(std::map<double, Vector<std::complex<double> > >(mVolumeIntegralJ));

			mVolumeValue.clear();
            mVolumeValue.swap(std::map<double, double>(mVolumeValue));
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters(mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ, 
				                 mSurfaceValue    , mProblemFrequency);

            WriteSurfiInFiles(mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ, 
				              mProblemFrequency);

		    mSurfaceIntegralE.clear();
            mSurfaceIntegralE.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralE));

			mSurfaceIntegralH.clear();
            mSurfaceIntegralH.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralH));

			mSurfaceIntegralJ.clear();
            mSurfaceIntegralJ.swap(std::map<double, Vector<std::complex<double> > >(mSurfaceIntegralJ));

			mSurfaceValue.clear();
            mSurfaceValue.swap(std::map<double, double>(mSurfaceValue));
        }

		//x_vector.resize(mSystemSize,0.00); 
        b_vector.resize(mSystemSize,0.00); 
        A_matrix.Resize(mSystemSize); 
	}

    //*******************************************************************************
    //* - Solving linear problem
    //*******************************************************************************
    void Modeler::SolveProblem()
    {
        SolveLinearSystem();

		//Setting to cero all the variables
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nit;
		for (nit = nodes.begin(); nit != nodes.end(); ++nit)
		{
			(*mpModel)(cEx, **nit) = czero;
			(*mpModel)(cEy, **nit) = czero;
			(*mpModel)(cEz, **nit) = czero;
		}

        SaveSystemSolution(cEx);
        SaveSystemSolution(cEy);
        SaveSystemSolution(cEz);

		UpdateAxisToCartesian();
		UpdateCoordEPEC();
		UpdateContactCoord();
		UpdateCoordEPBC();
		
		//x_vector.FreeData();

        //Keep mDofSet for the next iteration 
        //mDofSet.clear();
        //mDofSet.swap(std::set<Dof::Pointer,ComparePDof>(mDofSet));  
    }

	//*****************************************************************************
    //* - Updates coordinate system in axisymmetric problems
    //***************************************************************************** 
    void Modeler::UpdateAxisToCartesian()
    {
        if ( !mAxisymmetric ) return;
		
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;
		unsigned int nodeId;

        for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
        {
			nodeId = (*nodes_it)->Id();

            if ( ( mNormals.find(nodeId) == mNormals.end() )  && ( DistanceToAxis(nodeId) > 1e-12 ) )
            {
                std::vector<double> n(3), t(3), b(3);

				AxisymmetricCoordinates(nodeId,n,t,b);

                std::complex<double> En = (*mpModel)(cEx, **nodes_it),
					                 Et = (*mpModel)(cEy, **nodes_it),
                                     Eb = (*mpModel)(cEz, **nodes_it);

                (*mpModel)(cEx, **nodes_it) = n[0]*En + t[0]*Et + b[0]*Eb;
                (*mpModel)(cEy, **nodes_it) = n[1]*En + t[1]*Et + b[1]*Eb;
                (*mpModel)(cEz, **nodes_it) = n[2]*En + t[2]*Et + b[2]*Eb;
            }
        }
    }

	//*****************************************************************************
    //* - Change local coordinates in PEC for cartesian coordinates E3D
    //***************************************************************************** 
    void Modeler::UpdateCoordEPEC()
    {
        if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;
		
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;

        for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
        {
            if ( mNormals.find((*nodes_it)->Id()) != mNormals.end() )
            {
                std::vector<double> n(mNormals[(*nodes_it)->Id()]);
				std::vector<double> t(3), b(3);
                TangencialCoordinates(n,t,b);

                std::complex<double> En = (*mpModel)(cEx, **nodes_it),
					                 Et = (*mpModel)(cEy, **nodes_it),
                                     Eb = (*mpModel)(cEz, **nodes_it);

                (*mpModel)(cEx, **nodes_it) = n[0]*En + t[0]*Et + b[0]*Eb;
                (*mpModel)(cEy, **nodes_it) = n[1]*En + t[1]*Et + b[1]*Eb;
                (*mpModel)(cEz, **nodes_it) = n[2]*En + t[2]*Et + b[2]*Eb;
            }
        }

	    //Keep normal vectors for the next iteration
        //mNormals.clear();
        //mNormals.swap(std::map<unsigned int,  Vector<double> >(mNormals));  
    }

	//*****************************************************************************
    //* - Update electric field values in PBC nodes 
    //***************************************************************************** 
    void Modeler::UpdateCoordEPBC()
	{
		std::vector<unsigned int> FrontLeftNodesId;

		UpdateCoord_RL_PBC(FrontLeftNodesId);
		UpdateCoord_FB_PBC();
		UpdateCoord_LF_PBC(FrontLeftNodesId);
	}

	//*****************************************************************************
    //* - Check if the left node is also a front node 
    //***************************************************************************** 
	bool Modeler::PBC_LeftNodeIsFrontNode(unsigned int rightNodeId, std::vector<unsigned int>& leftNodesId, std::vector<unsigned int>& FrontLeftNodesId)
	{
		bool LeftNodeIsFrontNode = false;

		for ( int i = 0; i < leftNodesId.size(); i++ )
		{
			if ( ( mPBC_NodeNodePairs_Front.find   (leftNodesId[i]) != mPBC_NodeNodePairs_Front.end   () ) ||
				 ( mPBC_NodeEdgePairs_Front.find   (leftNodesId[i]) != mPBC_NodeEdgePairs_Front.end   () ) ||
				 ( mPBC_NodeElementPairs_Front.find(leftNodesId[i]) != mPBC_NodeElementPairs_Front.end() )  )
			{
				LeftNodeIsFrontNode = true;
				FrontLeftNodesId.push_back(rightNodeId);
				break;
			}
		}

		return LeftNodeIsFrontNode;
	}

	//*****************************************************************************
    //* - Update electric field values in PBC nodes (Right-Left faces)
    //***************************************************************************** 
    void Modeler::UpdateCoord_RL_PBC(std::vector<unsigned int>& FrontLeftNodesId)
	{
		if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 ) { return; }
		
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nit;

		for ( nit = nodes.begin(); nit!= nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();
			
			if ( mPBC_NodeNodePairs_Right.find(NodeId) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNode = mPBC_NodeNodePairs_Right[NodeId];

				std::vector<unsigned int> leftNodesId(1); leftNodesId[0] = lNode;

				if ( PBC_LeftNodeIsFrontNode( NodeId, leftNodesId, FrontLeftNodesId ) == true ) continue;
		
		        Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

				std::complex<double> lnEx = (*mpModel)(cEx, *plNode); 
				std::complex<double> lnEy = (*mpModel)(cEy, *plNode);
				std::complex<double> lnEz = (*mpModel)(cEz, *plNode);

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

				(*mpModel)(cEx, **nit) =  cos_A * lnEx + sin_A * lnEy;
                (*mpModel)(cEy, **nit) = -sin_A * lnEx + cos_A * lnEy;
                (*mpModel)(cEz, **nit) = lnEz;	  

				continue;
            }

			if ( mPBC_NodeEdgePairs_Right.find(NodeId) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<unsigned int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );	

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) continue;

				std::vector<std::complex<double> > lnEx(lNodes.size()); 
				std::vector<std::complex<double> > lnEy(lNodes.size()); 
				std::vector<std::complex<double> > lnEz(lNodes.size());

				for ( int i = 0; i < lNodes.size(); i++ )
				{ 
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );
					lnEx[i] = (*mpModel)(cEx, *plNode); 
					lnEy[i] = (*mpModel)(cEy, *plNode); 
					lnEz[i] = (*mpModel)(cEz, *plNode); 
				}

                std::vector<double> N(3);

                double L2 = mPBC_NodeEdgeCoord_Right[NodeId];
				double L1 = 1.00 - L2;

				N[0] = L1 * (2.00*L1 - 1.00); 
				N[1] = L2 * (2.00*L2 - 1.00); 
				N[2] = 4.00 * L1 * L2;

				std::complex<double> clEx = N[0] * lnEx[0] + N[1] * lnEx[1] + N[2] * lnEx[2];
				std::complex<double> clEy = N[0] * lnEy[0] + N[1] * lnEy[1] + N[2] * lnEy[2];
				std::complex<double> clEz = N[0] * lnEz[0] + N[1] * lnEz[1] + N[2] * lnEz[2];

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

				(*mpModel)(cEx, **nit) =  cos_A * clEx + sin_A * clEy;
                (*mpModel)(cEy, **nit) = -sin_A * clEx + cos_A * clEy;
                (*mpModel)(cEz, **nit) =  clEz;	  

				continue;
			}

			if ( mPBC_NodeElementPairs_Right.find(NodeId) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<unsigned int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );	

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) continue;

				std::vector<std::complex<double> > lnEx(lNodes.size()); 
				std::vector<std::complex<double> > lnEy(lNodes.size()); 
				std::vector<std::complex<double> > lnEz(lNodes.size());

				for ( int i = 0; i < lNodes.size(); i++ )
				{ 
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );
					lnEx[i] = (*mpModel)(cEx, *plNode); 
					lnEy[i] = (*mpModel)(cEy, *plNode); 
					lnEz[i] = (*mpModel)(cEz, *plNode); 
				}

 				std::vector<double> N(6);

				std::vector<double> natCoord(mPBC_NodeElementCoord_Right[NodeId]);
				
				double L2 = natCoord[0];
				double L3 = natCoord[1];
				double L1 = 1.00 - L2 - L3;
				
				N[0] = L1 * (2.00*L1 - 1.00);
				N[1] = L2 * (2.00*L2 - 1.00);
				N[2] = L3 * (2.00*L3 - 1.00);
				N[3] = 4.00 * L1 * L2;
				N[4] = 4.00 * L2 * L3;
				N[5] = 4.00 * L3 * L1;

				std::complex<double> clEx = N[0] * lnEx[0] + N[1] * lnEx[1] + N[2] * lnEx[2] + N[3] * lnEx[3] + N[4] * lnEx[4] + N[5] * lnEx[5];
				std::complex<double> clEy = N[0] * lnEy[0] + N[1] * lnEy[1] + N[2] * lnEy[2] + N[3] * lnEy[3] + N[4] * lnEy[4] + N[5] * lnEy[5];
				std::complex<double> clEz = N[0] * lnEz[0] + N[1] * lnEz[1] + N[2] * lnEz[2] + N[3] * lnEz[3] + N[4] * lnEz[4] + N[5] * lnEz[5];

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );
				
				(*mpModel)(cEx, **nit) =  cos_A * clEx + sin_A * clEy;
                (*mpModel)(cEy, **nit) = -sin_A * clEx + cos_A * clEy;
                (*mpModel)(cEz, **nit) =  clEz;	  

				continue;
			}
		}
	}

	//*****************************************************************************
    //* - Computes inverse rotation of a point
    //***************************************************************************** 
    void Modeler::PBC_InvRot( unsigned int NodeId, double& cos_A, double& sin_A )
    {
        double X = ( mpModel->GetNode( NodeId - 1 ) )->X();
		double Y = ( mpModel->GetNode( NodeId - 1 ) )->Y();

		double mod_r = sqrt(X*X + Y*Y); 
					
		cos_A = Y / mod_r;
		sin_A = X / mod_r;
	}

	//*****************************************************************************
    //* - Update electric field values in PBC nodes (Front-Back surfaces)
    //***************************************************************************** 
    void Modeler::UpdateCoord_FB_PBC()
	{
		if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 ) { return; }
		
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        
		std::vector<Node::Pointer>::iterator nit;

		for ( nit = nodes.begin(); nit!= nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();

			if ( mPBC_NodeNodePairs_Front.find(NodeId) != mPBC_NodeNodePairs_Front.end() )
            {
				int bNode = mPBC_NodeNodePairs_Front[NodeId];
		
		        Node::Pointer pbNode = mpModel->GetNode( bNode - 1 );

				std::complex<double> bnEx = (*mpModel)(cEx, *pbNode); 
				std::complex<double> bnEy = (*mpModel)(cEy, *pbNode);
				std::complex<double> bnEz = (*mpModel)(cEz, *pbNode);

				(*mpModel)(cEx, **nit) = bnEx;
                (*mpModel)(cEy, **nit) = bnEy;
                (*mpModel)(cEz, **nit) = bnEz;	

				continue;
            }

			if ( mPBC_NodeEdgePairs_Front.find(NodeId) != mPBC_NodeEdgePairs_Front.end() )
			{
				std::vector<unsigned int> bNodes( mPBC_NodeEdgePairs_Front[NodeId] );	

				std::vector<std::complex<double> > bnEx(bNodes.size()); 
				std::vector<std::complex<double> > bnEy(bNodes.size()); 
				std::vector<std::complex<double> > bnEz(bNodes.size());

				for ( int i = 0; i < bNodes.size(); i++ )
				{ 
					Node::Pointer pbNode = mpModel->GetNode( bNodes[i] - 1 );
					bnEx[i] = (*mpModel)(cEx, *pbNode); 
					bnEy[i] = (*mpModel)(cEy, *pbNode); 
					bnEz[i] = (*mpModel)(cEz, *pbNode); 
				}

                std::vector<double> N(3);

                double L2 = mPBC_NodeEdgeCoord_Front[NodeId];
				double L1 = 1.00 - L2;

				N[0] = L1 * (2.00*L1 - 1.00); 
				N[1] = L2 * (2.00*L2 - 1.00); 
				N[2] = 4.00 * L1 * L2;

				(*mpModel)(cEx, **nit) = N[0] * bnEx[0] + N[1] * bnEx[1] + N[2] * bnEx[2];
				(*mpModel)(cEy, **nit) = N[0] * bnEy[0] + N[1] * bnEy[1] + N[2] * bnEy[2];
				(*mpModel)(cEz, **nit) = N[0] * bnEz[0] + N[1] * bnEz[1] + N[2] * bnEz[2];

				continue;
			}

			if ( mPBC_NodeElementPairs_Front.find(NodeId) != mPBC_NodeElementPairs_Front.end() )
			{
				std::vector<unsigned int> bNodes( mPBC_NodeElementPairs_Front[NodeId] );	

				std::vector<std::complex<double> > bnEx(bNodes.size()); 
				std::vector<std::complex<double> > bnEy(bNodes.size()); 
				std::vector<std::complex<double> > bnEz(bNodes.size());

				for ( int i = 0; i < bNodes.size(); i++ )
				{ 
					Node::Pointer pbNode = mpModel->GetNode( bNodes[i] - 1 );		
					bnEx[i] = (*mpModel)(cEx, *pbNode); 
					bnEy[i] = (*mpModel)(cEy, *pbNode); 
					bnEz[i] = (*mpModel)(cEz, *pbNode); 
				}
				
				std::vector<double> N(6);

				std::vector<double> natCoord(mPBC_NodeElementCoord_Front[NodeId]);
				
				double L2 = natCoord[0];
				double L3 = natCoord[1];
				double L1 = 1.00 - L2 - L3;
				
				N[0] = L1 * (2.00*L1 - 1.00);
				N[1] = L2 * (2.00*L2 - 1.00);
				N[2] = L3 * (2.00*L3 - 1.00);
				N[3] = 4.00 * L1 * L2;
				N[4] = 4.00 * L2 * L3;
				N[5] = 4.00 * L3 * L1;

				(*mpModel)(cEx, **nit) = N[0] * bnEx[0] + N[1] * bnEx[1] + N[2] * bnEx[2] + N[3] * bnEx[3] + N[4] * bnEx[4] + N[5] * bnEx[5];
				(*mpModel)(cEy, **nit) = N[0] * bnEy[0] + N[1] * bnEy[1] + N[2] * bnEy[2] + N[3] * bnEy[3] + N[4] * bnEy[4] + N[5] * bnEy[5];
				(*mpModel)(cEz, **nit) = N[0] * bnEz[0] + N[1] * bnEz[1] + N[2] * bnEz[2] + N[3] * bnEz[3] + N[4] * bnEz[4] + N[5] * bnEz[5];

				continue;
			}
        }
	}

	//*****************************************************************************
    //* - Update electric field values in PBC nodes (Left-Front nodes)
    //***************************************************************************** 
    void Modeler::UpdateCoord_LF_PBC(std::vector<unsigned int>& FrontLeftNodes)
	{
		if ( FrontLeftNodes.size() == 0 ) { return; }
		
		std::vector<unsigned int>::iterator it_fln;

		for ( it_fln = FrontLeftNodes.begin(); it_fln!= FrontLeftNodes.end(); ++it_fln )
        {
			unsigned int NodeId = *it_fln;
			
			if ( mPBC_NodeNodePairs_Right.find(NodeId) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNode = mPBC_NodeNodePairs_Right[NodeId];

				Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

				std::complex<double> lnEx = (*mpModel)(cEx, *plNode); 
				std::complex<double> lnEy = (*mpModel)(cEy, *plNode);
				std::complex<double> lnEz = (*mpModel)(cEz, *plNode);

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

				Node::Pointer prNode = mpModel->GetNode( NodeId - 1 );

				(*mpModel)(cEx, *prNode) =  cos_A * lnEx + sin_A * lnEy;
                (*mpModel)(cEy, *prNode) = -sin_A * lnEx + cos_A * lnEy;
                (*mpModel)(cEz, *prNode) = lnEz;	  

				continue;
            }

			if ( mPBC_NodeEdgePairs_Right.find(NodeId) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<unsigned int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );	

				std::vector<std::complex<double> > lnEx(lNodes.size()); 
				std::vector<std::complex<double> > lnEy(lNodes.size()); 
				std::vector<std::complex<double> > lnEz(lNodes.size());

				for ( int i = 0; i < lNodes.size(); i++ )
				{ 
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );
					lnEx[i] = (*mpModel)(cEx, *plNode); 
					lnEy[i] = (*mpModel)(cEy, *plNode); 
					lnEz[i] = (*mpModel)(cEz, *plNode); 
				}

                std::vector<double> N(3);

                double L2 = mPBC_NodeEdgeCoord_Right[NodeId];
				double L1 = 1.00 - L2;

				N[0] = L1 * (2.00*L1 - 1.00); 
				N[1] = L2 * (2.00*L2 - 1.00); 
				N[2] = 4.00 * L1 * L2;

				std::complex<double> clEx = N[0] * lnEx[0] + N[1] * lnEx[1] + N[2] * lnEx[2];
				std::complex<double> clEy = N[0] * lnEy[0] + N[1] * lnEy[1] + N[2] * lnEy[2];
				std::complex<double> clEz = N[0] * lnEz[0] + N[1] * lnEz[1] + N[2] * lnEz[2];

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

				Node::Pointer prNode = mpModel->GetNode( NodeId - 1 );

				(*mpModel)(cEx, *prNode) =  cos_A * clEx + sin_A * clEy;
                (*mpModel)(cEy, *prNode) = -sin_A * clEx + cos_A * clEy;
                (*mpModel)(cEz, *prNode) =  clEz;	  

				continue;
			}

			if ( mPBC_NodeElementPairs_Right.find(NodeId) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<unsigned int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );	

				std::vector<std::complex<double> > lnEx(lNodes.size()); 
				std::vector<std::complex<double> > lnEy(lNodes.size()); 
				std::vector<std::complex<double> > lnEz(lNodes.size());

				for ( int i = 0; i < lNodes.size(); i++ )
				{ 
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );
					lnEx[i] = (*mpModel)(cEx, *plNode); 
					lnEy[i] = (*mpModel)(cEy, *plNode); 
					lnEz[i] = (*mpModel)(cEz, *plNode); 
				}

 				std::vector<double> N(6);

				std::vector<double> natCoord(mPBC_NodeElementCoord_Right[NodeId]);
				
				double L2 = natCoord[0];
				double L3 = natCoord[1];
				double L1 = 1.00 - L2 - L3;
				
				N[0] = L1 * (2.00*L1 - 1.00);
				N[1] = L2 * (2.00*L2 - 1.00);
				N[2] = L3 * (2.00*L3 - 1.00);
				N[3] = 4.00 * L1 * L2;
				N[4] = 4.00 * L2 * L3;
				N[5] = 4.00 * L3 * L1;

				std::complex<double> clEx = N[0] * lnEx[0] + N[1] * lnEx[1] + N[2] * lnEx[2] + N[3] * lnEx[3] + N[4] * lnEx[4] + N[5] * lnEx[5];
				std::complex<double> clEy = N[0] * lnEy[0] + N[1] * lnEy[1] + N[2] * lnEy[2] + N[3] * lnEy[3] + N[4] * lnEy[4] + N[5] * lnEy[5];
				std::complex<double> clEz = N[0] * lnEz[0] + N[1] * lnEz[1] + N[2] * lnEz[2] + N[3] * lnEz[3] + N[4] * lnEz[4] + N[5] * lnEz[5];

				double cos_A = 1.0; 
				double sin_A = 0.0; 
				
				if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

				Node::Pointer prNode = mpModel->GetNode( NodeId - 1 );
				
				(*mpModel)(cEx, *prNode) =  cos_A * clEx + sin_A * clEy;
                (*mpModel)(cEy, *prNode) = -sin_A * clEx + cos_A * clEy;
                (*mpModel)(cEz, *prNode) =  clEz;	  

				continue;
			}
		}
	}

	//*****************************************************************************
    //* - Updating contact coordinates
    //***************************************************************************** 
    void Modeler::UpdateContactCoord()
    {
        if ( mContactPairs.size() == 0 ) return;
		
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;

        for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
        {
			if ( mContactPairs.find((*nodes_it)->Id()) != mContactPairs.end() )
            {
				ContactPairData& pContactPair = mContactPairs[(*nodes_it)->Id()];

				Node::Pointer pLNode = mpModel->GetNode((pContactPair.ctcNode) - 1);

                std::complex<double> ExL = (*mpModel)(cEx, *pLNode),
					                 EyL = (*mpModel)(cEy, *pLNode),
				                     EzL = (*mpModel)(cEz, *pLNode);
                //Properties
				Properties::Pointer Properties;

				double eo = 8.8541878176e-12;

				double freq = mProblemFrequency;

                double sigma, eps_real, eps_imag;

				//Material R
				Properties = mpModel->GetProperties(pContactPair.matR);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

				//Material L
				Properties = mpModel->GetProperties(pContactPair.matL);

				sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
                eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
               
		        std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

				//Properties matrix
				std::complex<double> ContactMatrix[3][3];

                std::complex<double> rmat = (cEpsL/cEpsR) - 1.00;

				Vector<double> n = pContactPair.ctcNormal;

				ContactMatrix[0][0] = (n[0]*n[0]) * rmat + 1.00; ContactMatrix[0][1] = (n[0]*n[1]) * rmat       ; ContactMatrix[0][2] = (n[0]*n[2]) * rmat;
				ContactMatrix[1][0] = (n[1]*n[0]) * rmat       ; ContactMatrix[1][1] = (n[1]*n[1]) * rmat + 1.00; ContactMatrix[1][2] = (n[1]*n[2]) * rmat;
				ContactMatrix[2][0] = (n[2]*n[0]) * rmat       ; ContactMatrix[2][1] = (n[2]*n[1]) * rmat       ; ContactMatrix[2][2] = (n[2]*n[2]) * rmat + 1.00;

				//Update
                (*mpModel)(cEx, **nodes_it) = ContactMatrix[0][0]*ExL + ContactMatrix[0][1]*EyL + ContactMatrix[0][2]*EzL;
                (*mpModel)(cEy, **nodes_it) = ContactMatrix[1][0]*ExL + ContactMatrix[1][1]*EyL + ContactMatrix[1][2]*EzL;
				(*mpModel)(cEz, **nodes_it) = ContactMatrix[2][0]*ExL + ContactMatrix[2][1]*EyL + ContactMatrix[2][2]*EzL;
            }
        }
    }

	//*******************************************************************************
    //* - Setting normals in a PEC 
    //*******************************************************************************
	void Modeler::SettingPECNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetPECnormals.begin(); it_Set != mSetPECnormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                Vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

			//Setting normal and deegrees of freedom
            mNormals[it_Set->first] = ac_normal;
			
			Node::Pointer pNode = mpModel->GetNode((it_Set->first)-1);
                
			(*mpModel)(cEy, *pNode) = czero;
            pNode->pDofcEy()->FixDof();

			(*mpModel)(cEz, *pNode) = czero;
            pNode->pDofcEz()->FixDof();
        } 

		mSetPECnormals.clear();
        mSetPECnormals.swap(std::map<unsigned int, std::vector<std::vector<double> > >(mSetPECnormals));
	}

	//*******************************************************************************
    //* - Setting normals in a PMC 
    //*******************************************************************************
	void Modeler::SettingSPECNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetSPECnormals.begin(); it_Set != mSetSPECnormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                Vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vToBeOrto(mNormals[it_Set->first]);
				Vector<double> vOrto(3,0.00);

				Ortogonalization(vToBeOrto, ac_normal, vOrto);

				mNormals[it_Set->first] = vOrto;
			}
			else
			{
				mNormals[it_Set->first] = ac_normal;
			
				Node::Pointer pNode = mpModel->GetNode((it_Set->first)-1);
                
				(*mpModel)(cEy, *pNode) = czero;
				pNode->pDofcEy()->FixDof();

				(*mpModel)(cEz, *pNode) = czero;
				pNode->pDofcEz()->FixDof();
			}
        } 

		mSetSPECnormals.clear();
        mSetSPECnormals.swap(std::map<unsigned int, std::vector<std::vector<double> > >(mSetSPECnormals));
	}

	//*******************************************************************************
    //* - Setting normals in a PMC 
    //*******************************************************************************
	void Modeler::SettingPMCNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

		for(it_Set = mSetPMCnormals.begin(); it_Set != mSetPMCnormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                std::vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vToBeOrto(mNormals[it_Set->first]);
				Vector<double> vOrto(3,0.00);

				Ortogonalization(vToBeOrto, ac_normal, vOrto);

				mNormals[it_Set->first] = vOrto;
			}
			else
			{
				mNormals[it_Set->first] = ac_normal;
			
				Node::Pointer pNode = mpModel->GetNode((it_Set->first)-1);
                
				(*mpModel)(cEx, *pNode) = czero;
				pNode->pDofcEx()->FixDof();
			}
        } 

        mSetPMCnormals.clear();
        mSetPMCnormals.swap(std::map<unsigned int, std::vector<std::vector<double> > >(mSetPMCnormals));
	}

	//*******************************************************************************
    //* - Setting normals in a TE port 
    //*******************************************************************************
	void Modeler::SettingTEPMCNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

		for(it_Set = mSetTEPMCnormals.begin(); it_Set != mSetTEPMCnormals.end(); ++it_Set) 
        {
            std::vector<std::vector<double> > vectorSet(it_Set->second);
            std::vector<std::vector<double> >::iterator it_vs;
            Vector<double> ac_normal(3,0.00);

            for(it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs)
            {
                std::vector<double> normal_in_set(*it_vs);
                ac_normal += normal_in_set;
            }

            ac_normal *= (1.00/sqrt(std::inner_product(ac_normal.begin(),ac_normal.end(),ac_normal.begin(),0.00)));

			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vToBeOrto(mNormals[it_Set->first]);
				Vector<double> vOrto(3,0.00);

				Ortogonalization(vToBeOrto, ac_normal, vOrto);

				mNormals[it_Set->first] = vOrto;
			}
			else
			{
				mNormals[it_Set->first] = ac_normal;
			
				Node::Pointer pNode = mpModel->GetNode((it_Set->first)-1);
                
				(*mpModel)(cEx, *pNode) = czero;
				pNode->pDofcEx()->FixDof();
			}
        } 

        mSetTEPMCnormals.clear();
        mSetTEPMCnormals.swap(std::map<unsigned int, std::vector<std::vector<double> > >(mSetTEPMCnormals));
	}

	//*******************************************************************************
    //* - Making the vector vToBeOrto ortogonal to the vector vAxis 
    //*******************************************************************************
	void Modeler::Ortogonalization(Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto)
    {
		//vOrto == (vAxis X vToBeOrto) X vAxis, 
		//vOrto is ortogonal to vAxis and pointing in the same direction as vToBeOrto

		Vector<double> nDummy(3,0.00);

		//nDummy  = (vAxis X vToBeOrto)
		nDummy[0] = vAxis[1]*vToBeOrto[2] - vAxis[2]*vToBeOrto[1];
		nDummy[1] = vAxis[2]*vToBeOrto[0] - vAxis[0]*vToBeOrto[2];
		nDummy[2] = vAxis[0]*vToBeOrto[1] - vAxis[1]*vToBeOrto[0];

		//vOrto  = (nDummy X vAxis) 
		vOrto[0] = nDummy[1]*vAxis[2] - nDummy[2]*vAxis[1];
		vOrto[1] = nDummy[2]*vAxis[0] - nDummy[0]*vAxis[2];
		vOrto[2] = nDummy[0]*vAxis[1] - nDummy[1]*vAxis[0];

		vOrto *= (1.00/sqrt(std::inner_product(vOrto.begin(), vOrto.end(), vOrto.begin(),0.00)));
	}

	//*******************************************************************************
    //* - Ortogonalization of the contact normals respect to PEC and PMC normals
    //*******************************************************************************
	void Modeler::OrtogonalizeContactNormals()
    {
		std::map<unsigned int, ContactPairData>::iterator itCtc;

        for(itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc) 
        {
			unsigned int fCtcNode  = itCtc->first;
			
			if ( (mNormals.find(fCtcNode) != mNormals.end()) )
            {
				Vector<double> sNormal(mNormals[fCtcNode]);
				Vector<double> vToBeOrto(mContactPairs[fCtcNode].ctcNormal);
				Vector<double> vOrto(3,0.00);

				Ortogonalization(vToBeOrto, sNormal, vOrto);

				mContactPairs[fCtcNode].ctcNormal = vOrto;
			}
		}
	}

    //*******************************************************************************
    //* - Setting normals in PEC and PMC surfaces
    //*******************************************************************************
    void Modeler::SetNormals()
    {
        //We add to mSetPECnormals the contacts nodes in edges and PEC corners 
		AddContactPECNormals();

		//Setting PEC normals
		SettingPECNormals();
		SettingSPECNormals();

		//PEC normals are forced in dielectric-PEC intersections.
		ForcePECNormalsInContact();

		//Setting PMC normals
		SettingPMCNormals();
		SettingTEPMCNormals();

		//Contact normals must be ortogonal to any normal
		OrtogonalizeContactNormals();
    }

    //*******************************************************************************
    //* - Initialize Building process
    //*******************************************************************************
    void Modeler::IniBuilding()
    {
        std::cout << "Nodes read." << std::endl << std::endl;
        
		std::cout << "Setting normals and DOFs..." << std::endl;
		
		SetContactDiscontinuity();
		SetNormals();
		SetPBC    ();

		if (mSetAllEzToZero ) SetAllcEzToCero ();
		if (mSetAllExyToZero) SetAllcExyToCero();
		if (mAxisymmetric)    SetAllcEyzToCero();

        SetDofSet();

		std::cout << "Done." << std::endl << std::endl;

		if (mDebugSolutionMode) return; 
    }

	//*******************************************************************************
    //* - Set all cEz to cero (3D-Exy problems)
    //*******************************************************************************
    void Modeler::SetAllcEzToCero()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator it;
		std::complex<double> czero(0.00,0.00);

		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			Node::Pointer pNode = *it;

			(*mpModel)(cEz, *pNode) = czero;
			pNode->pDofcEz()->FixDof();
		}
	}

	//*******************************************************************************
    //* - Set all cEx cEy to cero (3D-Ez problems)
    //*******************************************************************************
    void Modeler::SetAllcExyToCero()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator it;
		std::complex<double> czero(0.00,0.00);

		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			Node::Pointer pNode = *it;

			(*mpModel)(cEx, *pNode) = czero;
			pNode->pDofcEx()->FixDof();

			(*mpModel)(cEy, *pNode) = czero;
			pNode->pDofcEy()->FixDof();
		}
	}

	//*******************************************************************************
    //* - Set all cEy cEz to cero (Axisimmetric problems, actually Er and Ephs)
    //*******************************************************************************
    void Modeler::SetAllcEyzToCero()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator it;
		std::complex<double> czero(0.00,0.00);

		Node::Pointer pNode;
		unsigned int nodeId;
		
		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			pNode  = *it;
			nodeId = pNode->Id();

			(*mpModel)(cEy, *pNode) = czero;
			pNode->pDofcEy()->FixDof();

			(*mpModel)(cEz, *pNode) = czero;
			pNode->pDofcEz()->FixDof();

			if ( ( mNormals.find(nodeId)  != mNormals.end() ) || 
				 ( DistanceToAxis(nodeId) <= 1e-12 ) )
            {	
				(*mpModel)(cEx, *pNode) = czero;
				pNode->pDofcEx()->FixDof();
			}
		}
	}

    //*******************************************************************************
    //* - Solve linear system
    //*******************************************************************************
    void Modeler::SolveLinearSystem()
    {    
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::cout << "Max. distance to singularity: " << mMaxDistToSing << "  in nodes: "<< mMaxDistNode1 << " , " << mMaxDistNode2 << std::endl;
		//std::cout << "Min. distance to singularity: " << mMinDistToSing << "  in nodes: "<< mMinDistNode1 << " , " << mMinDistNode2 << std::endl;
		//std::cout << std::endl;
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		if (mReadSolutionMode)
		{
			ReadSolutionVector(x_vector);
			A_matrix.FreeData();
			b_vector.FreeData();
			return;
		}

		if (mSolveWithExternal)
		{
			WriteMatrixInFile(A_matrix);
			WriteVectorInFile(b_vector);
			WriteInSolvedFile(0);

			A_matrix.FreeData();
			b_vector.FreeData();

			if (mExternalSolverPath == "")
			{
				std::cout << "ERMES matrix generation finished - solving with external solver." << std::endl;
				std::cout << "Read solution when external solver finishes.                    " << std::endl << std::endl;
				std::cout << "ERMES finished.                                                 " << std::endl << std::endl;
				exit(EXIT_SUCCESS);
		    }

			// This fragment calls the external solver an waits for its solution
			std::cout << "ERMES matrix generation finished - solving with external solver." << std::endl << std::endl;
			std::cout << mExternalSolverPath                                                << std::endl << std::endl;
			std::cout << "Waiting until external solver finishes...                       " << std::endl << std::endl;
            
			const char* MSDOSOrder = mExternalSolverPath;

			system(MSDOSOrder);

			std::cout << std::endl << "External solver has finished." << std::endl << std::endl;

			ReadSolutionVector(x_vector);
		}
		else
		{
			std::cout << "Solving linear system..." << std::endl << std::endl;

			if (mReadFileInitialGuess) ReadInitialGuess(x_vector);

			mpLinearSolver->SetWriteStepResult(mWriteResultEveryStep);

            mpLinearSolver->Solve(A_matrix, b_vector, x_vector);

            if (mWriteResultFinalStep) WriteResultsInFile(x_vector);

            std::cout << "Solver finished." << std::endl << std::endl;
		}
		
        A_matrix.FreeData();
        b_vector.FreeData(); 
    }

    //*************************************************
    //* - Set singular nodes
    //*************************************************
    void Modeler::FixDof(unsigned int NodeId, 
		                 const Variable<double>& rVariable, 
		                 const double& Value)
    {
		int varKey = rVariable.getKey();

		if(varKey == Sg.getKey()) mSingular[NodeId] = Value;
    }

    //*************************************************
    //* - Set complex degrees of freedom
    //*************************************************
    void Modeler::FixCDof(unsigned int NodeId, 
		                  const Variable<double>& rVariable, 
                          const std::vector<double>& vValue)
    {
    }

	//*************************************************
    //* - Free complex degrees of freedom
    //*************************************************
    void Modeler::FreeCDof(unsigned int NodeId, 
		                   const Variable<double>& rVariable)
    {
          Node::Pointer pnode = mpModel->GetNode(NodeId-1);

          int varKey = rVariable.getKey();
          
          if     (varKey == Ex.getKey()) pnode->pDofcEx()->FreeDof();
		  else if(varKey == Ey.getKey()) pnode->pDofcEy()->FreeDof();
		  else if(varKey == Ez.getKey()) pnode->pDofcEz()->FreeDof();
    }

    //*******************************************************************************
    //* - Set Linear Solver
    //*******************************************************************************
    void Modeler::GenerateLinearSolver(String Name, 
				                       String SolverType, 
				                       int MaxNumberOfIterations,
									   int StepIterations,
				                       double Tolerance)
    {
        mSolveWithBiCG      = false;
		mSolveWithQMR       = false;
		mSolveWithCG        = false;
		mSolveWithExternal  = false;
		mExternalSolverPath = "";

		int NumTreads = mItSolverNumThreads;

		if (SolverType == "Bi_Conjugate_Gradient")
		{
			mSolveWithBiCG = true;
			mpLinearSolver = ComplexSolver::Pointer(new BiCGComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}
		else if (SolverType == "Quasi_Minimal_Residual")
		{
			mSolveWithQMR  = true;
			mpLinearSolver = ComplexSolver::Pointer(new QMRComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}	
		else if (SolverType == "Conjugate_Gradient")
		{
			mSolveWithCG   = true;
			mpLinearSolver = ComplexSolver::Pointer(new CGComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}	
		else if (SolverType == "External_solver")
		{
			mSolveWithExternal = true;
			ReadExternalPathFromFile(mBaseFileName, mExternalSolverPath);
		}
		else 
		{
			mSolveWithQMR  = true;
			mpLinearSolver = ComplexSolver::Pointer(new QMRComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}
    }

    //*****************************************************************************************
    //* - Set Dofs
    //*****************************************************************************************
    void Modeler::SetDofSet()
    {
        std::set<Dof::Pointer,ComparePDof>::iterator it3;

        int i=0;
        int j=mDofSet.size();

        for (it3 = mDofSet.begin(); it3 != mDofSet.end(); ++it3)
        {
            if ( (*it3)->IsFixed() )
            {
	            (*it3)->SetEquationId(--j);
            }
            else
            {
	            (*it3)->SetEquationId(i++);
            }
        }

        mSystemSize = j;

        x_vector.resize(mSystemSize,0.00); 
        b_vector.resize(mSystemSize,0.00); 
        A_matrix.Resize(mSystemSize);  
    }
	
	//*******************************************************
    //* - Set projection
    //*******************************************************
	void Modeler::PrintProjection()
    {
		if ( mProjection.size() > 0 ) 
		{
		     CalculateSijParameters();
			 PrintSijParameters(mProjectionNormalized, mProblemFrequency);
		}

		if ( mVolumeIntegralE.size() > 0 )
        {
            PrintViParameters(mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ, 
				              mVolumeValue    , mProblemFrequency);
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters(mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ, 
				                 mSurfaceValue    , mProblemFrequency);
		}
    }

    //*******************************************************
    //* - Setting problem type
    //*******************************************************
    void Modeler::SetProblemType(const String& ProblemType)
    {
		if (ProblemType == "E3D") 
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = false;
			mAxisymmetric    = false;
			mPlasmaMode      = false;
		}
		else if (ProblemType == "Ez3D") 
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = true;
			mAxisymmetric    = false;
			mPlasmaMode      = false;
		}
		else if (ProblemType == "Exy3D") 
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = true;
			mSetAllExyToZero = false;
			mAxisymmetric    = false;
			mPlasmaMode      = false;
		}
		else if (ProblemType == "Ea3D") 
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = false;
			mAxisymmetric    = true;
			mPlasmaMode      = false;
		}
		else if (ProblemType == "Plasma")
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = false;
			mAxisymmetric    = false;
			mPlasmaMode      = true;
		}
		else if (ProblemType == "RELSSOL")
		{
			mReleaseSolutionMode = true;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = false;		
		}
		else if (ProblemType == "DEBGSOL")
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = true;
			mReadSolutionMode    = false;
		}
		else if (ProblemType == "READSOL")
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = true;
		}
	    else if (ProblemType == "Project") 
		{
			PrintProjection();
		}
		else if (ProblemType == "GAv") 
		{
			mNormalsGeomAvg = true;
			mNormalsAreaWtd = false;
		}
		else if (ProblemType == "AWg")
		{
			mNormalsGeomAvg = false;
			mNormalsAreaWtd = true;
		}
		else if (ProblemType == "RFIG")
		{
			mReadFileInitialGuess = true;
		}
		else if (ProblemType == "NRFIG")
		{
			mReadFileInitialGuess = false;
		}
		else if (ProblemType == "NWRIF")
		{
			mWriteResultEveryStep = false;
			mWriteResultFinalStep = false;
		}
		else if (ProblemType == "ESWRIF")
		{
			mWriteResultEveryStep = true;
			mWriteResultFinalStep = true;
		}
		else if (ProblemType == "FSWRIF")
		{
			mWriteResultEveryStep = false;
			mWriteResultFinalStep = true;
		}
		else if (ProblemType == "OFFBIN")
		{
			mOutputFileFormatIsBIN   = true;
			mOutputFileFormatIsASCII = false;
		}
		else if (ProblemType == "OFFASCII")
		{
			mOutputFileFormatIsBIN   = false;
			mOutputFileFormatIsASCII = true;
		}
		else if (ProblemType == "QE")
		{
			mQuadraticGeometry = true;
			mElementOrder      = 2;
		}
		else if (ProblemType == "LE")
		{
			mQuadraticGeometry = false;
		}
	    else if (ProblemType == "1st") mElementOrder = 1;
        else if (ProblemType == "2nd") mElementOrder = 2;
        else if (ProblemType == "3rd") mElementOrder = 3;
        else if (ProblemType == "4th") mElementOrder = 4;
		else if (ProblemType == "1pr" ) mItSolverNumThreads = 1 ;
        else if (ProblemType == "2pr" ) mItSolverNumThreads = 2 ;
        else if (ProblemType == "3pr" ) mItSolverNumThreads = 3 ;
        else if (ProblemType == "4pr" ) mItSolverNumThreads = 4 ;
		else if (ProblemType == "5pr" ) mItSolverNumThreads = 5 ;
        else if (ProblemType == "6pr" ) mItSolverNumThreads = 6 ;
        else if (ProblemType == "7pr" ) mItSolverNumThreads = 7 ;
        else if (ProblemType == "8pr" ) mItSolverNumThreads = 8 ;
		else if (ProblemType == "9pr" ) mItSolverNumThreads = 9 ;
        else if (ProblemType == "10pr") mItSolverNumThreads = 10;
		else if (ProblemType == "11pr") mItSolverNumThreads = 11;
		else if (ProblemType == "12pr") mItSolverNumThreads = 12;
		else if (ProblemType == "13pr") mItSolverNumThreads = 13;
		else if (ProblemType == "14pr") mItSolverNumThreads = 14;
		else if (ProblemType == "15pr") mItSolverNumThreads = 15;
		else if (ProblemType == "16pr") mItSolverNumThreads = 16;
    }

	//*******************************************************
    //* - Giving element order
    //*******************************************************
    int Modeler::ElementOrder()
    {
        return mElementOrder;
    }

	//*******************************************************
	//* - Returns true if the output file format is binary.
	//*******************************************************
	bool Modeler::Is_OutputFileFormat_BIN()
	{
		return mOutputFileFormatIsBIN;
	}

	//*******************************************************
	//* - Returns true if the output file format is ASCII.
	//*******************************************************
	bool Modeler::Is_OutputFileFormat_ASCII()
	{
		return mOutputFileFormatIsASCII;
	}

    //*******************************************************
    //* - Setting frequency
    //*******************************************************
    void Modeler::SetFrequency(double ProblemFrequency)
    {
        mProblemFrequency = ProblemFrequency;
    }

    //*******************************************************
    //* - Calculates of E 3D
    //*******************************************************
    void Modeler::CalculateNodalE3D(int variableKey)
    {
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        //E(t)
        if (variableKey == E.getKey())
        {	
             std::vector<double> E_total(3),
                                 E_real, E_imag;

             double time = mpModel->CurrentTime();

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 E_real = (*mpModel)(REAL_E, **it);
                 E_imag = (*mpModel)(IMAG_E, **it);

                 E_total[0] = E_real[0]*cos(mProblemFrequency*time) + 
                              E_imag[0]*sin(mProblemFrequency*time);

                 E_total[1] = E_real[1]*cos(mProblemFrequency*time) + 
                              E_imag[1]*sin(mProblemFrequency*time); 

                 E_total[2] = E_real[2]*cos(mProblemFrequency*time) + 
                              E_imag[2]*sin(mProblemFrequency*time); 

                 (*mpModel)(E, **it) = E_total;
             }    
        } 
		//H(t)
        if (variableKey == H.getKey())
        {	
             std::vector<double> H_total(3),
                                 H_real, H_imag;

             double time = mpModel->CurrentTime();

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 H_real = (*mpModel)(REAL_H, **it);
                 H_imag = (*mpModel)(IMAG_H, **it);

                 H_total[0] = H_real[0]*cos(mProblemFrequency*time) + 
                              H_imag[0]*sin(mProblemFrequency*time);

                 H_total[1] = H_real[1]*cos(mProblemFrequency*time) + 
                              H_imag[1]*sin(mProblemFrequency*time); 

                 H_total[2] = H_real[2]*cos(mProblemFrequency*time) + 
                              H_imag[2]*sin(mProblemFrequency*time); 

                 (*mpModel)(H, **it) = H_total;
             }    
        } 
		//J(t)
        if (variableKey == J.getKey())
        {	
             std::vector<double> J_total(3),
                                 J_real, J_imag;

             double time = mpModel->CurrentTime();

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 J_real = (*mpModel)(REAL_J, **it);
                 J_imag = (*mpModel)(IMAG_J, **it);

                 J_total[0] = J_real[0]*cos(mProblemFrequency*time) + 
                              J_imag[0]*sin(mProblemFrequency*time);

                 J_total[1] = J_real[1]*cos(mProblemFrequency*time) + 
                              J_imag[1]*sin(mProblemFrequency*time); 

                 J_total[2] = J_real[2]*cos(mProblemFrequency*time) + 
                              J_imag[2]*sin(mProblemFrequency*time); 

                 (*mpModel)(J, **it) = J_total;
             }    
        } 
        //REAL_E
        else if (variableKey == REAL_E.getKey())
        {
			 std::vector<double> E_real(3);

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
				 E_real[0] = std::real((*mpModel)(cEx, **it));
                 E_real[1] = std::real((*mpModel)(cEy, **it));
                 E_real[2] = std::real((*mpModel)(cEz, **it));

                 (*mpModel)(REAL_E, **it) = E_real;
             }
        } 
        //IMAG_E 
        else if (variableKey == IMAG_E.getKey())
        {   
             std::vector<double> E_imag(3);

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
				 E_imag[0] = std::imag((*mpModel)(cEx, **it));
                 E_imag[1] = std::imag((*mpModel)(cEy, **it));
                 E_imag[2] = std::imag((*mpModel)(cEz, **it));

                 (*mpModel)(IMAG_E, **it) = E_imag;
             } 
        }
        //MOD_E
        else if (variableKey == MOD_E.getKey())
        {    
             std::vector<double> E_real, E_imag;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 E_real = (*mpModel)(REAL_E, **it); 
                 E_imag = (*mpModel)(IMAG_E, **it);

	             (*mpModel)(MOD_E, **it) = sqrt( E_real[0]*E_real[0] + E_imag[0]*E_imag[0] +
                                                 E_real[1]*E_real[1] + E_imag[1]*E_imag[1] +
                                                 E_real[2]*E_real[2] + E_imag[2]*E_imag[2] );
             }
        }
		//LOG_MOD_E
		else if (variableKey == LOG_MOD_E.getKey())
        {    
             std::vector<double> E_real, E_imag;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 E_real = (*mpModel)(REAL_E, **it); 
                 E_imag = (*mpModel)(IMAG_E, **it);

				 (*mpModel)(LOG_MOD_E, **it) = log10( sqrt( E_real[0]*E_real[0] + E_imag[0]*E_imag[0] +
                                                            E_real[1]*E_real[1] + E_imag[1]*E_imag[1] +
                                                            E_real[2]*E_real[2] + E_imag[2]*E_imag[2] ) );
             }
        }
		//REAL_H 
		else if (variableKey == REAL_H.getKey())
        {
			 Vector<double> H_real;

             int NumNeighbours;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 NumNeighbours = mNeigbourElements[(*it)->Id()-1];
             
				 H_real  = (*mpModel)(REAL_H, **it);

				 if (NumNeighbours > 0)
				 { 
				     H_real /= NumNeighbours;
				 }

                 (*mpModel)(REAL_H, **it) = H_real;
             } 
        }
		//IMAG_H 
        else if (variableKey == IMAG_H.getKey())
        {   
             Vector<double> H_imag;

             int NumNeighbours;
                
             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 NumNeighbours = mNeigbourElements[(*it)->Id()-1];

                 H_imag  = (*mpModel)(IMAG_H, **it);

				 if (NumNeighbours > 0)
				 {
					 H_imag /= NumNeighbours;
				 }

                 (*mpModel)(IMAG_H, **it) = H_imag;
             } 
        }
		//MOD_H
        else if (variableKey == MOD_H.getKey())
        {    
             std::vector<double> H_real, H_imag;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 H_real = (*mpModel)(REAL_H, **it); 
                 H_imag = (*mpModel)(IMAG_H, **it);

	             (*mpModel)(MOD_H, **it) = sqrt( H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                                                 H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                                                 H_real[2]*H_real[2] + H_imag[2]*H_imag[2] );
             }
        }
		//LOG_MOD_H
		else if (variableKey == LOG_MOD_H.getKey())
        {    
             std::vector<double> H_real, H_imag;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 H_real = (*mpModel)(REAL_H, **it); 
                 H_imag = (*mpModel)(IMAG_H, **it);

				 (*mpModel)(LOG_MOD_H, **it) = log10( sqrt( H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                                                            H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                                                            H_real[2]*H_real[2] + H_imag[2]*H_imag[2] ) );
             }
        }
		//BOUNDARY_NORMALS
        else if (variableKey == BOUNDARY_NORMALS.getKey())
        {
             for (it=nodes.begin(); it!= nodes.end(); ++it)
             {
				 int nodoid=(*it)->Id();

				 std::vector<double> ncero(3,0.00);
				 std::vector<double> n(3,0.00);

				 if ( mNormals.find(nodoid) != mNormals.end() )
				 {   
				     n[0] = mNormals[nodoid][0];
					 n[1] = mNormals[nodoid][1];
					 n[2] = mNormals[nodoid][2];

					 (*mpModel)(BOUNDARY_NORMALS, **it) = n;
				 }
				 else
				 {
					 (*mpModel)(BOUNDARY_NORMALS, **it) = ncero;
				 }
             }
        } 
        //CONTACT_NORMALS
        else if (variableKey == CONTACT_NORMALS.getKey())
        {   
             for (it=nodes.begin(); it!= nodes.end(); ++it)
             {
				  int nodoid=(*it)->Id();

				  std::vector<double> ncero(3,0.00);
				  std::vector<double> n(3,0.00);

				  if ( mContactPairs.find(nodoid) != mContactPairs.end() )
                  {   
                     n[0] = mContactPairs[nodoid].ctcNormal[0];
					 n[1] = mContactPairs[nodoid].ctcNormal[1];
					 n[2] = mContactPairs[nodoid].ctcNormal[2];

                     (*mpModel)(CONTACT_NORMALS, **it) = n;
                  }
                  else
                  {
                     (*mpModel)(CONTACT_NORMALS, **it) = ncero;
                  }
             } 
        }
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Singular elements surface - Debugger mode
        else if (variableKey == MOD_H.getKey())
        {    
             std::vector<double> E_real, E_imag;

             for (it=nodes.begin(); it!= nodes.end(); ++it)
             {
				 int nodoid=(*it)->Id();
				 
				 if      (    mSingular.find(nodoid) != mSingular.end()    ) (*mpModel)(MOD_H, **it) = 1.00;
				 else if ( mNSingular2L.find(nodoid) != mNSingular2L.end() ) (*mpModel)(MOD_H, **it) = 2.00;
				 else if ( mNSingular3L.find(nodoid) != mNSingular3L.end() ) (*mpModel)(MOD_H, **it) = 3.00;
				 else if ( mNSingular4L.find(nodoid) != mNSingular4L.end() ) (*mpModel)(MOD_H, **it) = 4.00;
				 else if ( mNSingular5L.find(nodoid) != mNSingular5L.end() ) (*mpModel)(MOD_H, **it) = 5.00;
				 else if ( mNSingular6L.find(nodoid) != mNSingular6L.end() ) (*mpModel)(MOD_H, **it) = 6.00;
				 else                                                        (*mpModel)(MOD_H, **it) = 0.00;
             }
        }
		//////////////////////////////////////////////////////////////////////////////////////////////////////////*/
    }

	//*******************************************************
    //* - Initial set up to calculate H field.
    //*******************************************************
	void Modeler::IniHDerivation()
	{
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		int i;
        
        int size = nodes.size();

        mNeigbourElements.resize(size);

		for (i=0; i<size; i++) mNeigbourElements[i] = 0;
		
		std::vector<Node::Pointer>::iterator it; 
        
        Vector<double> cero(3);

		for (i=0; i<3; i++) cero[i] = 0.00;

        for (it=nodes.begin(); it!=nodes.end(); ++it)
        {
			 (*mpModel)(REAL_H, **it) = cero; 
             (*mpModel)(IMAG_H, **it) = cero; 
        }  
	}

	//*******************************************************
    //* - Initial set up to calculate J.
    //*******************************************************
	void Modeler::IniJCalculation()
	{
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
       
        Vector<double> cero(3);

		for (int i=0; i<3; i++) cero[i] = 0.00;

	    std::vector<Node::Pointer>::iterator it; 

        for (it=nodes.begin(); it!=nodes.end(); ++it)
        {
			 (*mpModel)(REAL_J, **it) = cero; 
             (*mpModel)(IMAG_J, **it) = cero;
        }  
	}

	//*******************************************************
    //* - Calculate H field in elements.
    //*******************************************************
	void Modeler::CalculateHElement(int* NodesId, unsigned int PropertiesId)
	{
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
		{
			mNeigbourElements[(*hoit)-1] += 1;
			nodes.push_back(mpModel->GetNode((*hoit)-1));
		}
    
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);
        
        Element::Pointer pElement;

		if      ( mElementOrder == 1 ) pElement = Element::Pointer(new VolumeElement_1st(nodes,properties));
		else if ( mElementOrder == 2 ) pElement = Element::Pointer(new VolumeElement_2nd(nodes,properties));
		else if ( mElementOrder == 3 ) pElement = Element::Pointer(new VolumeElement_3th(nodes,properties));
		else if ( mElementOrder == 4 ) pElement = Element::Pointer(new VolumeElement_4th(nodes,properties));

        Vector<std::complex<double> > ncHx, ncHy, ncHz; 
    
        pElement->CalculateDerivatives(ncHx, ncHy, ncHz);

        Vector<Node::Pointer>::iterator node_it; 

		Vector<double> rH_vector(3);
		Vector<double> iH_vector(3);
		Vector<double> v_temp(3); 

		int n = 0;
   
        for(node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
			rH_vector[0] = std::real(ncHx[n]);
			rH_vector[1] = std::real(ncHy[n]);
			rH_vector[2] = std::real(ncHz[n]);

			iH_vector[0] = std::imag(ncHx[n]);
			iH_vector[1] = std::imag(ncHy[n]);
			iH_vector[2] = std::imag(ncHz[n]);

			v_temp  = (*mpModel)(REAL_H, **node_it);
            v_temp += rH_vector;
            (*mpModel)(REAL_H, **node_it) = v_temp;

            v_temp  = (*mpModel)(IMAG_H, **node_it); 
            v_temp += iH_vector; 
            (*mpModel)(IMAG_H, **node_it) = v_temp; 

			n++;
        }
	}

	//*******************************************************
    //* - Calculate Joule heating in elements.
    //*******************************************************
	void Modeler::CalculateJouleHElement(int* NodesId, unsigned int PropertiesId)
	{
		//Detecting if the element is in a medium with losses.
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);
		
		double eo   = 8.8541878176e-12;
		double freq = mProblemFrequency;

		double sigma    = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);
		double eps_imag = (*properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		double eff_cond = freq*eps_imag + sigma;

		bool HasLosses; 

		if (eff_cond > 0) HasLosses = true;
		else              HasLosses = false;

		//Only medium with losses will have Joule heating
		if (!HasLosses) return;
		
		//Recovering all the nodes of the element
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));
    
		std::vector<Node::Pointer>::iterator it;

		std::vector<double> E_real, E_imag;

		for (it=nodes.begin(); it!= nodes.end(); ++it)
		{
			E_real = (*mpModel)(REAL_E, **it); 
			E_imag = (*mpModel)(IMAG_E, **it);

			double modE2 = E_real[0]*E_real[0] + E_imag[0]*E_imag[0] +
						   E_real[1]*E_real[1] + E_imag[1]*E_imag[1] +
						   E_real[2]*E_real[2] + E_imag[2]*E_imag[2] ;

			(*mpModel)(JOULE_HEATING, **it) = 0.5*eff_cond*modE2;
		}	
	}

    //****************************************************************************
    //* - Calculate J with VE elements.
    //****************************************************************************
	void Modeler::CalculateJVEElement(int* NodesId, unsigned int PropertiesId)
	{
		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

		double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

		if ( sigma == 0.00 ) return;
					
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));
    
		Vector<Node::Pointer>::iterator node_it; 

		Vector<double> J_vector(3);

		for(node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
		{
			J_vector[0] = sigma * std::real((*mpModel)(cEx, **node_it));
			J_vector[1] = sigma * std::real((*mpModel)(cEy, **node_it));
			J_vector[2] = sigma * std::real((*mpModel)(cEz, **node_it));

			(*mpModel)(REAL_J, **node_it) = J_vector;
			
			J_vector[0] = sigma * std::imag((*mpModel)(cEx, **node_it));
			J_vector[1] = sigma * std::imag((*mpModel)(cEy, **node_it));
			J_vector[2] = sigma * std::imag((*mpModel)(cEz, **node_it));

			(*mpModel)(IMAG_J, **node_it) = J_vector; 
		}
	}

	//*****************************************************************************
    //* - Calculate J with JE elements.
    //*****************************************************************************
	void Modeler::CalculateJJEElement(int* NodesId, unsigned int PropertiesId)
	{
		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        std::vector<double> vJ = (*properties)(SINUSOIDAL_SURFACE_CURRENT);

		std::complex<double> cJx( vJ[0]*cos(vJ[1]), vJ[0]*sin(vJ[1]) );
		std::complex<double> cJy( vJ[2]*cos(vJ[3]), vJ[2]*sin(vJ[3]) );
		std::complex<double> cJz( vJ[4]*cos(vJ[5]), vJ[4]*sin(vJ[5]) );

		std::vector<double> vJa = (*properties)(COMPLEX_IBC);

		std::complex<double> cJa( vJa[0]*cos(vJa[1]), vJa[0]*sin(vJa[1]) );

		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));
    
		Vector<Node::Pointer>::iterator node_it; 

		Vector<double> J_vector(3);

		double rx, rz, norm;

		int nodeId;

		for(node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
			rx   = -((*node_it)->Z());
			rz   =  ((*node_it)->X());

		    norm = sqrt( rx*rx + rz*rz );

			if ( norm > 0.00 )
			{
				rx = rx/norm;
				rz = rz/norm;
			}
			else
			{
				rx = 0.00;
				rz = 0.00;
			}

			nodeId = (*node_it)->Id();

			J_vector[0] = std::real(cJx) + std::real(cJa)*rx;
            J_vector[1] = std::real(cJy);
            J_vector[2] = std::real(cJz) + std::real(cJa)*rz;

			tmp_rJ[nodeId] = J_vector;

			J_vector[0] = std::imag(cJx) + std::imag(cJa)*rx;
            J_vector[1] = std::imag(cJy);
            J_vector[2] = std::imag(cJz) + std::imag(cJa)*rz;

			tmp_iJ[nodeId] = J_vector;
        }
	}

	//*******************************************************
    //* - End calculate J.
    //*******************************************************
	void Modeler::EndJCalculation()
	{
        int i;

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        
		std::vector<Node::Pointer>::iterator node_it; 

		Vector<double> Eddy_vector(3, 0.00);

		Vector<double> rJtot(3, 0.00);
		Vector<double> iJtot(3, 0.00);

        for (node_it=nodes.begin(); node_it!=nodes.end(); ++node_it)
        {
			int nodeId = (*node_it)->Id();
			
			// REAL_J
			Eddy_vector = (*mpModel)(REAL_J, **node_it);

			if ( tmp_rJ.find(nodeId) != tmp_rJ.end() )
			{
				for (i=0; i<3; i++) rJtot[i] = Eddy_vector[i] + tmp_rJ[nodeId][i];
			}
			else
			{
				for (i=0; i<3; i++) rJtot[i] = Eddy_vector[i];
			}

			(*mpModel)(REAL_J, **node_it) = rJtot;
			
			// IMAG_J
			Eddy_vector  = (*mpModel)(IMAG_J, **node_it);

			if ( tmp_iJ.find(nodeId) != tmp_iJ.end() ) 
			{
				for (i=0; i<3; i++) iJtot[i] = Eddy_vector[i] + tmp_iJ[nodeId][i];
			}
			else 
			{
				for (i=0; i<3; i++) iJtot[i] = Eddy_vector[i];
			}

			(*mpModel)(IMAG_J, **node_it) = iJtot;
			
			// MOD_J
	        (*mpModel)(MOD_J , **node_it) = sqrt( rJtot[0]*rJtot[0] + iJtot[0]*iJtot[0] +
                                                  rJtot[1]*rJtot[1] + iJtot[1]*iJtot[1] +
                                                  rJtot[2]*rJtot[2] + iJtot[2]*iJtot[2] );
        }  

		tmp_rJ.clear();
        tmp_rJ.swap(std::map<unsigned int, std::vector<double> >(tmp_rJ));

		tmp_iJ.clear();
        tmp_iJ.swap(std::map<unsigned int, std::vector<double> >(tmp_iJ));
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	// - Projections
    //
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*******************************************************
    //* - Calculate Sij parameters
    //*******************************************************
    void Modeler::CalculateSijParameters()
	{
		std::map<double,  std::complex<double> >::iterator proyec_it;

		for(proyec_it = mProjection.begin(); proyec_it != mProjection.end(); proyec_it++)
		{
		    std::complex<double> p = proyec_it->second/(mNormalization[proyec_it->first]);
           
			if ( proyec_it->first != 0.00 ) 
            {
				mProjectionNormalized[proyec_it->first] = p;
            }
			else
            {
				mProjectionNormalized[proyec_it->first] = p - 1.00;   
            }
		}
        
        mNormalization.clear();
        mNormalization.swap(std::map<double, std::complex<double> >(mNormalization));
	    mProjection.clear();
        mProjection.swap(std::map<double, std::complex<double> >(mProjection));
	}
 
} // Namespace Kratos
