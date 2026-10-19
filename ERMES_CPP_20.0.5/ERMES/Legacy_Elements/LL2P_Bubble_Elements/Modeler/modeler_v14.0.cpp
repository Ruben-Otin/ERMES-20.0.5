
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

#include "../ERMES/VolumeElement_1st_Plasma.h"
#include "../ERMES/FarField_1st_Plasma.h"
#include "../ERMES/GenericRobin_1st_Plasma.h"

#include "../ERMES/VolumeElement_2ndQ_DivD_Plasma.h"
#include "../ERMES/FarField_2ndQ_DivD_Plasma.h"
#include "../ERMES/GenericRobin_2ndQ_Plasma.h"

#include "../ERMES/VolumeElement_1st.h"
#include "../ERMES/JSource_1st.h"
#include "../ERMES/FarField_1st.h"
#include "../ERMES/RWPortTE10_1st.h"
#include "../ERMES/CoaxialPortTEM_1st.h"
#include "../ERMES/GenericRobin_1st.h"

#include "../ERMES/VolumeElement_1bb.h"
#include "../ERMES/JSource_1bb.h"
#include "../ERMES/FarField_1bb.h"

#include "../ERMES/VolumeElement_1st_Static.h"
#include "../ERMES/GenericRobin_1st_Static.h"

#include "../ERMES/Utils/MatrixFileWritting.h"
#include "../ERMES/Utils/LagrangeBases.h"
#include "../ERMES/Utils/GaussIntegrationTables.h"

#include "../includes/modeler.h"
#include "../includes/model.h"
#include "../linear_solvers/kratos_linear_solvers.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#ifdef _WIN64
    #include <direct.h>
    #include "dirent/include/dirent.h"
#else
    #include <dirent.h>
    #include <sys/stat.h>
#endif


namespace Kratos
{
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Print error messages and exit ERMES
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Send_Error_Msg( String ErrorMsg, String CheckMsg, int VerticalSpaces )
    {
         for ( int i=0; i<VerticalSpaces; i++ )
         {
             std::cout << std::endl;
         }
         std::cout << "!!!!!!!!!! " << ErrorMsg << " !!!!!!!!!!" << std::endl;
         std::cout << std::endl;
         std::cout << CheckMsg << std::endl;
         std::cout << std::endl;
         std::cout << "ERMES analysis finished." << std::endl;
         exit(EXIT_SUCCESS);
    }

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
		 for ( int i=0; i<4; i++ ) mTempNodeProperties[NodesId[i]] = PropertiesId;
	}

	//*************************************************************************************************
    //* - Generate pairs of nodes in a contact surface.
    //*************************************************************************************************
	void Modeler::GenerateContactPairs(Vector<int>& CtCNodesId)
	{
        if (mStaticMode == true) return;
 
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

		double sigma_Last3    = (*Last3_prop)(IHL_ELECTRIC_CONDUCTIVITY);
		double real_eps_Last3 = (*Last3_prop)(REAL_ELECTRIC_PERMITTIVITY) * eo;
		double imag_eps_Last3 = (*Last3_prop)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

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
    //* - Calculates geometric averaged normal in contact surfaces
    //*************************************************************************************************
	void Modeler::CalculateContactNormals(Vector<int>& NodesId)
	{
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

		// Normal = V2 x V1 = 2 * element_area * n
		normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
				    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

		normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
				    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

		normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
				    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        //if ( mElementOrder == 0 )
        //{
      //      std::vector<double> cX; 
		    //std::vector<double>	cY;
		    //std::vector<double>  W;
      //      
		    //int nGaussPoints = GaussPoints2D_Order8(cX, cY, W);
      //      
		    //Matrix<double> N;
      //      
		    //Lagrange2D_Ni_2bb(N, cX, cY);

      //      double intS_Ni_X = 0;
      //      double intS_Ni_Y = 0;
      //      double intS_Ni_Z = 0;

      //      Vector<double> Ns(normal);

      //      Ns /= sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
      //  
      //      for (int i=0; i<3; i++)
      //      {
      //          for (int gp=0; gp<nGaussPoints; gp++) 
      //          {
      //              intS_Ni_X += W[gp] * N[i][gp] * Ns[0];
      //              intS_Ni_Y += W[gp] * N[i][gp] * Ns[1];
      //              intS_Ni_Z += W[gp] * N[i][gp] * Ns[2];
      //          }
      //      }

      //      normal[0] = intS_Ni_X;
      //      normal[1] = intS_Ni_Y;
      //      normal[2] = intS_Ni_Z;
        //}

        // Pushing back the normal vector of the element to the contact normal set at each node
        for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
		{
			mSetContactNormals[*id_it].push_back(normal);
		}
	}

	//*************************************************************************************************
    //* - Set average normals in contact surfaces
    //*************************************************************************************************
	void Modeler::SetAvgNormalsInContact()
	{
		// Setting average normals in surface dilelectric nodes
		std::map<unsigned int, std::vector<std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetContactNormals.begin(); it_Set != mSetContactNormals.end(); ++it_Set)
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet(it_Set->second);

            // Averaged normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd == true ) normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            else                           normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );

            // Including the calculated normal in the mNormal vector
            mContactNormals[it_Set->first] = normal_atNode;
        }

		mSetContactNormals.clear();
        std::map<unsigned int,std::vector<std::vector<double> > >().swap(mSetContactNormals);
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

            std::complex<double> cZero(0.00,0.00);

			(*mpModel)(cEx, *pRNode) = cZero;
			(*mpModel)(cEy, *pRNode) = cZero;
			(*mpModel)(cEz, *pRNode) = cZero;

			pRNode->pDofcEx()->FixDof();
			pRNode->pDofcEy()->FixDof();
			pRNode->pDofcEz()->FixDof();

            //pRNode->pDofcEx()->FreeDof();
			//pRNode->pDofcEy()->FreeDof();
			//pRNode->pDofcEz()->FreeDof();

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

            mType_Of_BC_Normal[ pRNode->Id() ] = 'R';
            mType_Of_BC_Normal[ pLNode->Id() ] = 'L';

            mLNormal[ pLNode->Id() ] = pContactPair.ctcNormal;
            mRNormal[ pRNode->Id() ] = pContactPair.ctcNormal;
		}
        
		DelTempContactVectors();

		////////////////////////////////////////////////////////////////////////////////////////
		//std::map<unsigned int, ContactPairData>::iterator it;
		//for(it = mContactPairs.begin(); it != mContactPairs.end(); ++it)
		//{
		//	std::cout<<"Node1: "<<it->first            <<"  mat1: "<<(it->second).matR<<"   "
		//			 <<"Node2: "<<(it->second).ctcNode <<"  mat2: "<<(it->second).matL<<"   "
		//			 <<"n :"<<(it->second).ctcNormal
		//			 <<std::endl;
		//}
		///////////////////////////////////////////////////////////////////////////////////////
	}

	//*************************************************************************************************
    //* - Delete temporal vectors generated along contact discontinuity making.
    //*************************************************************************************************
	void Modeler::DelTempContactVectors()
	{
		mTempNodeProperties.clear();
		std::map<unsigned int, unsigned short int>().swap(mTempNodeProperties);

		mContactNormals.clear();
		std::map<unsigned int, std::vector<double> >().swap(mContactNormals);

		mTempContactPairs.clear();
		std::map<unsigned int, unsigned int>().swap(mTempContactPairs);
	}

	//*************************************************************************************************
    //* - Joining normals to second node from the first node
    //*************************************************************************************************
	void Modeler::JoinContactNormalsInDirichletS()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		// Unify normal in Dirichlet surface with second node in dielectric pair
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

            if ( mSetPSBCnormals.find(FirstNode) != mSetPSBCnormals.end() )
			{
				unsigned int SecondNode = (itCtc->second).ctcNode;

				std::vector<std::vector<double> > FirstVectorSet (mSetPSBCnormals[FirstNode]);
				std::vector<std::vector<double> >::iterator it_Fvs;

				for(it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs)
				{
					Vector<double> vectorF(*it_Fvs);
					mSetPSBCnormals[SecondNode].push_back(vectorF);
				}
			}
		}
	}

	//*************************************************************************************************
    //* - Force Dirichlet normals in all contact nodes that are in a Dirichlet surface
    //*************************************************************************************************
	void Modeler::ForceDirichletNormalsInContact()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		// New dielectric normal. All nodes in a Dirichlet surface have the same normal.
		// PEC/PSBC/PMC/TEPMC dominates normal direction.
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

                std::complex<double> cZero(0.00,0.00);

				(*mpModel)(cEx, *pNode) = cZero;
				(*mpModel)(cEy, *pNode) = cZero;
				(*mpModel)(cEz, *pNode) = cZero;

				pNode->pDofcEx()->FixDof();
				pNode->pDofcEy()->FixDof();
				pNode->pDofcEz()->FixDof();
			}
		}

		mUsefulNodes.clear(); std::vector<bool>().swap(mUsefulNodes);
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
            if ( mSingular.find   (NodesId[i]) != mSingular.end   () ) return;
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
            if ( mSingular.find   (NodesId[i]) != mSingular.end   () ) return;
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
            if ( mSingular.find   (NodesId[i]) != mSingular.end   () ) return;
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
	void Modeler::PushHONodesOnVolume( int* NodesId, std::vector<int>& HONodesId )
	{
		HONodesId.push_back( NodesId[0] );
        HONodesId.push_back( NodesId[1] );
        HONodesId.push_back( NodesId[2] );
		HONodesId.push_back( NodesId[3] );

        if ( mElementOrder == 0 ) 
        {
            HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[2], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[1], NodesId[2], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[2], NodesId[3], HONodesId );
            HONodesOnVolume( NodesId   ,                         HONodesId );
        }            
             
		if ( mElementOrder <= 1 ) return;

        HONodesOnEdge( NodesId[0], NodesId[1], HONodesId );
        HONodesOnEdge( NodesId[1], NodesId[2], HONodesId );
        HONodesOnEdge( NodesId[2], NodesId[0], HONodesId );
		HONodesOnEdge( NodesId[0], NodesId[3], HONodesId );
        HONodesOnEdge( NodesId[1], NodesId[3], HONodesId );
        HONodesOnEdge( NodesId[2], NodesId[3], HONodesId );

		if ( mElementOrder == 2 ) return;

		HONodesOnFace( NodesId[0], NodesId[1], NodesId[2], HONodesId );
		HONodesOnFace( NodesId[0], NodesId[1], NodesId[3], HONodesId );
		HONodesOnFace( NodesId[1], NodesId[2], NodesId[3], HONodesId );
		HONodesOnFace( NodesId[0], NodesId[2], NodesId[3], HONodesId );

		if ( mElementOrder == 3 ) return;

		HONodesOnVolume( NodesId, HONodesId );
	}

	//*************************************************************************************************
    //* - Add high order elements to the HONodesId vector in surface elements
    //*************************************************************************************************
	void Modeler::PushHONodesOnSurface( int* NodesId, std::vector<int>& HONodesId )
	{
		HONodesId.push_back( NodesId[0] );
        HONodesId.push_back( NodesId[1] );
        HONodesId.push_back( NodesId[2] );

        if ( mElementOrder == 0 ) 
        {
            HONodesOnFace( NodesId[0], NodesId[1], NodesId[2], HONodesId );
        }

		if ( mElementOrder <= 1 ) return;

        HONodesOnEdge( NodesId[0], NodesId[1], HONodesId );
        HONodesOnEdge( NodesId[1], NodesId[2], HONodesId );
        HONodesOnEdge( NodesId[2], NodesId[0], HONodesId );

        if ( mElementOrder == 2 ) return;

		HONodesOnFace( NodesId[0], NodesId[1], NodesId[2], HONodesId );
	}

	//*************************************************************************************************
    //* - Add high order elements to the HONodesId vector in surface elements
    //*************************************************************************************************
	void Modeler::PushHONodesOnSurface( std::vector<int>& NodesId, std::vector<int>& HONodesId )
	{
		HONodesId.push_back( NodesId[0] );
        HONodesId.push_back( NodesId[1] );
        HONodesId.push_back( NodesId[2] );

        if ( mElementOrder == 0 ) 
        {
            HONodesOnFace( NodesId[0], NodesId[1], NodesId[2], HONodesId );
        }

		if ( mElementOrder <= 1 ) return;

        HONodesOnEdge( NodesId[0], NodesId[1], HONodesId );
        HONodesOnEdge( NodesId[1], NodesId[2], HONodesId );
        HONodesOnEdge( NodesId[2], NodesId[0], HONodesId );

        if ( mElementOrder == 2 ) return;

		HONodesOnFace( NodesId[0], NodesId[1], NodesId[2], HONodesId );
	}

	//*************************************************************************************************
    //* - Looking for high order nodes on edge
    //*************************************************************************************************
    void Modeler::HONodesOnEdge( int IdNode1, int IdNode2, std::vector<int>& HONodesId )
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
	void Modeler::HONodesOnFace( int IdNode1, int IdNode2, int IdNode3, std::vector<int>& HONodesId )
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
				 if      ( mElementOrder == 0 ) ReorderHONodesOnFace_4th(orderingCase, face, HONodesId);
                 else if ( mElementOrder == 3 ) HONodesId.push_back(face[3]);
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
	void Modeler::HONodesOnVolume(int* NodesId, std::vector<int>& HONodesId)
	{
	    unsigned int key = NodesId[0]*NodesId[1]*NodesId[2]*NodesId[3];
        
        std::vector<std::vector<unsigned int> > SameKeyElement(mVolmHONodes[key]);
        std::vector<std::vector<unsigned int> >::iterator it_elements;
        
        for(it_elements = SameKeyElement.begin(); it_elements != SameKeyElement.end(); ++it_elements)
        {
            std::vector<unsigned int> nInElement(*it_elements);
        
            if ( (nInElement[0] == NodesId[0]) && (nInElement[1] == NodesId[1]) && 
                 (nInElement[2] == NodesId[2]) && (nInElement[3] == NodesId[3]) )
            {
                HONodesId.push_back(nInElement[4]);
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
        if ( mElementOrder == 0 ) 
        {
            CreateHONodes_OnFaces ( NodesId );
            CreateHONodes_OnVolume( NodesId ); 
        }

        if ( mElementOrder <= 1 ) return; 

		CreateHONodes_OnEdges( NodesId );

		if ( mElementOrder == 2 ) return; 
       
		CreateHONodes_OnFaces( NodesId );

		if ( mElementOrder == 3 ) return; 
      
		CreateHONodes_OnVolume( NodesId );
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
    //* - Creates high order nodes on edges
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
    //* - Creates high order nodes on faces
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
							if      ( mElementOrder == 0 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
                            else if ( mElementOrder == 3 ) CreateHONodes_OnFace    (IdNode1, IdNode2, IdNode3);
							else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
						}
					}
					else
					{
					    if      ( mElementOrder == 0 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
                        else if ( mElementOrder == 3 ) CreateHONodes_OnFace    (IdNode1, IdNode2, IdNode3);
					    else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th(IdNode1, IdNode2, IdNode3);
					}
				}
			}
		}
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
    //* - Creating a new higher order node on face 
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnFace(int IdNode1, int IdNode2, int IdNode3)
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
    //* - Creating a new higher order node inside volume 
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnVolume(int* IdNodes)
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
    //* - Checks if the problem set-up is consistent
    //*************************************************************************************************
    void Modeler::CheckConsistency()
    {
        // Checking mode consistency
        if    ( (mProblemFrequency == 0.0) && (mStaticMode == false) )
        {
            Send_Error_Msg ( "Problem frequency equal to 0",
                             "Please, check problem frequency.", 1 );
        }
        else if ( (mFrequencySweep == true) && (mStaticMode == true) )
        {
            Send_Error_Msg ( "Sweep frequency mode not supported in static analysis",
                             "Please, change frequency mode or problem mode.", 1 );
        }
        else if ( (mPlasmaMode == true) && (mElementOrder > 2) )
        {
            Send_Error_Msg( "Plasma mode not implemented for element orders higher than 2",
                            "Please, change element order or problem mode.", 1 );
        }
        else if ( (mPlasmaMode == true) && (mElementOrder == 0) )
        {
            Send_Error_Msg( "Plasma mode not implemented for bubble elements",
                            "Please, change element type or problem mode.", 1 );
        }
        else if ( (mStaticMode == true) && (mElementOrder != 1) )
        {
            Send_Error_Msg( "Static mode just implemented for 1st order elements",
                            "Please, change element order or problem mode.", 1 );
        }
        else if ( (mStaticMode == true) && (mDebugSolutionMode == true) )
        {
            Send_Error_Msg( "Debug solution mode not supported in static analysis",
                            "Please, change solution mode or problem mode.", 1 );
        }
        else if ( (mQuadraticGeometry == true) && (mElementOrder != 2) )
        {
            Send_Error_Msg( "Quadratic meshing just supported for 2nd order elements",
                            "Please, change element order or meshing mode.", 1 );
        }

        // Static mode do not import currents
        if ( mStaticMode == true )  mImportCurrents = false;
    }

    //*************************************************************************************************
    //* - Assembling volume elements in electro-static mode
    //*************************************************************************************************
    void Modeler::GenerateVolumeElement_Static( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertyId == 0 throw error, material not assigned
        if ( PropertiesId == 0 ) { throw 0; }

        // Getting local nodes
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

        // Getting material properties
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        // Creating element
        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Static (nodes, properties) );

        // Getting DOFs Id vector
        Vector<int> EleIdVector;

        pElement->GetEquationIdVector( EleIdVector );

        // Getting elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffMatrix( EleStiffMatrix );

        // Including Dirichlet BC in the residual vector
        Vector<double> EleResVector;

        pElement->GetResidualVector_Dirichlet( mFixVoltage, EleStiffMatrix, EleResVector );

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

        // Assembling matrix
        for ( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdVector[i_local];

            if ( i_global < mSystemSize )
            {
                for ( int j_local = 0; j_local < eleSize; j_local++ )
                {
                    int j_global = EleIdVector[j_local];

                    if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( abs( cKij ) > 0.0 )  A_matrix( i_global, j_global ) += cKij;
                    }
                }
            }
        }

        //Assembling global source vector
        for ( int i_blocal=0; i_blocal < eleSize; i_blocal++ )
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
    //* - Assembling high order volume elements in plasma mode
    //*************************************************************************************************
    void Modeler::GenerateVolumeElement_Plasma( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertyId == 0 throw error, material not assigned
        if ( PropertiesId == 0 ) { throw 0; }

        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Plasma element definition
        Element::Pointer pElement;

        if   ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_Plasma      ( nodes, properties ) );
		else                        pElement = Element::Pointer( new VolumeElement_2ndQ_DivD_Plasma( nodes, properties ) );

        pElement->Set_Plasma_Model( mpColdPlasma );

        bool Is_IHL = pElement->Set_Plasma_or_IHL();

        // Checking if the element is inside a sigularity layer
        bool IsUgLayer = false;

        for ( int i=0; i<4; i++ )
        {
            if ( ( mSingular.find   ( NodesId[i] ) != mSingular.end   () ) || ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) ||
                 ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() ) || ( mNSingular4L.find( NodesId[i] ) != mNSingular4L.end() ) ||
                 ( mNSingular5L.find( NodesId[i] ) != mNSingular5L.end() ) || ( mNSingular6L.find( NodesId[i] ) != mNSingular6L.end() )  )
            {
                IsUgLayer = true;
            }
        }

        if ( IsUgLayer == true ) pElement->SetPeso( 0.00 );
        else                     pElement->SetPeso( 1.00 );

        // Getting equation id vector
        Vector<int> EleIdVector;       

        pElement->GetEquationIdVector( EleIdVector );

        // Getting element stiffness matrix
        Matrix<std::complex<double> > EleStiffMatrix; 

        pElement->GetStiffMatrix( EleStiffMatrix );

        // Trans(T)*A*T for PBC conditions
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Trans(T)*A*T for contact discontinuities
        ContactHO_TKT_Plasma( HONodesId, EleIdVector, EleStiffMatrix );

        // Trans(T)*A*T for PEC, PMC conditions
        HO_TKT( HONodesId, EleStiffMatrix );

        // Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        Plasma_TKT( HONodesId, EleStiffMatrix );

        // Assembling local matrix to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

        double Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        if ( mpColdPlasma->Is_Full_Matrix() == true )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];

                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];

                        if ( ( j_global < mSystemSize ) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];
                            if ( abs( cKij ) > Kij_Tol ) A_matrix( i_global, j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else if ( Is_IHL == false )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];

                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];

                        if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];
                            if ( abs( cKij ) > Kij_Tol ) A_matrix( i_global, j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else if ( mpColdPlasma->Is_HermSymm_Matrix() == true )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];

                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];

                        if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];
                            if ( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global, j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else if ( mpColdPlasma->Is_HermFull_Matrix() == true )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];

                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];

                        if ( ( j_global < mSystemSize ) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];
                            if( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global, j_global ) += cKij;
                        }
                    }
                }
            }
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************
    //* - Assembling high order volume elements in standard mode
    //*************************************************************************************************
    void Modeler::GenerateVolumeElement( int* NodesId, unsigned int PropertiesId )
    {
		// If PropertyId == 0 throw error, material not assigned
		if ( PropertiesId == 0 ) { throw 0; }

		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		//PushHONodesOnVolume( NodesId, HONodesId );

        if ( mElementOrder == 4 ) 
        {
		    HONodesId.push_back( NodesId[0] );
            HONodesId.push_back( NodesId[1] );
            HONodesId.push_back( NodesId[2] );
		    HONodesId.push_back( NodesId[3] );

            HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[2], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[1], NodesId[2], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[2], NodesId[3], HONodesId );
            HONodesOnVolume( NodesId   ,                         HONodesId );               
        }
        else
        {
            PushHONodesOnVolume( NodesId, HONodesId );
        }

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;
        
		if      ( mQuadraticGeometry ) pElement = Element::Pointer( new VolumeElement_2ndQ( nodes, properties ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_1bb ( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st ( nodes, properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd ( nodes, properties ) );
		else if ( mElementOrder == 3 ) pElement = Element::Pointer( new VolumeElement_3th ( nodes, properties ) );
		else if ( mElementOrder == 4 ) pElement = Element::Pointer( new VolumeElement_1bb ( nodes, properties ) );//pElement = Element::Pointer( new VolumeElement_4th ( nodes, properties ) );

		// Checking if the element is inside a sigularity layer
		bool IsUgLayer = false;

		for ( int i=0; i<4; i++ )
		{
			if ( ( mSingular.find   ( NodesId[i] ) != mSingular.end   () ) || ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) ||
				 ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() ) || ( mNSingular4L.find( NodesId[i] ) != mNSingular4L.end() ) ||
				 ( mNSingular5L.find( NodesId[i] ) != mNSingular5L.end() ) || ( mNSingular6L.find( NodesId[i] ) != mNSingular6L.end() )  )
			{
				IsUgLayer = true;
			}
		}
        
        // Set weighted regularization in elements near a singularity
		if ( IsUgLayer == true ) pElement->SetPeso( 0.00 );
		else                     pElement->SetPeso( 1.00 );
        
        // Setting connectivities for bubble elements
        if ( (mElementOrder == 0) || (mElementOrder == 4) ) 
        {    
            /*
            Vector< Vector< Node::Pointer > > CtC_pNodes    ( 4 );
            Vector<  std::complex<double>   > CtC_Properties( 4,  std::complex<double>( 0.00, 0.00 ) );

            Vector< Vector<int> > Faces; Faces.resize( 4 );

		    Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
            Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		    Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		    Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;
            
            for (int iFace = 0; iFace<4 ; iFace++ )
            {
                int ctcNodeId = HONodesId[ Faces[iFace][3] ];
                
                if ( mContactPairs.find( ctcNodeId ) != mContactPairs.end() ) 
                {
                    Properties::Pointer Properties = mpModel->GetProperties( mContactPairs[ ctcNodeId ].matL );

                    double eo   = 8.8541878176e-12;
				    double freq = mProblemFrequency;
                    
				    double sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
                    double eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			        double eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
                    
		            std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

                    CtC_Properties[ iFace ] = cEpsL;

                    CtC_pNodes[ iFace ].resize( 6 );

                    for ( int pnf = 0; pnf<6; pnf++ ) CtC_pNodes[iFace][pnf] = mpModel->GetNode( mContactPairs[ HONodesId[ Faces[iFace][pnf] ] ].ctcNode - 1 );
                }
            }
            */
            pElement->Set_Mesh_Connectivities( mPropertiesMap, mNodesConnectivity, mElementsVector, mType_Of_BC_Normal );
            pElement->Set_Global_Derivatives ( mbVolumes, mb_dNk_dX, mb_dNk_dY, mb_dNk_dZ );
            //pElement->Set_Contact_Info       ( CtC_pNodes, CtC_Properties );
        }
       
        // Getting equation id vector
        Vector<int> EleIdVector; 
        
        pElement->GetEquationIdVector( EleIdVector );

        // Getting element stiffness matrix
        Matrix<std::complex<double> > EleStiffMatrix; 
        
        pElement->GetStiffMatrix( EleStiffMatrix );

        //if ( mElementOrder == 0 ) 
        //{ 
        //    HONodesId.resize(4); 
        //    for(int i=0; i<4; i++ ) HONodesId[i] = NodesId[i]; 
        //}
        
		//Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
	    PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

		//Trans(T)*A*T for contact discontinuities
		ContactHO_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        //Trans(T)*A*T for PEC, PMC conditions
        HO_TKT( HONodesId, EleStiffMatrix );

		//Trans(T)*A*T for axisymmetric problems
		Axisym_TKT( HONodesId, EleStiffMatrix );

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

        // Assembling symmetric matrix
		for ( int i_local = 0; i_local < eleSize; i_local++ )
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize )
			{
				for ( int j_local = 0; j_local < eleSize; j_local++ )
				{
					int j_global = EleIdVector[j_local];

					if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( abs( cKij ) > 0.0 ) A_matrix( i_global, j_global ) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

    //*************************************************************************************************
    //* - 
    //*************************************************************************************************
    void Modeler::Apply_Jump_Integral( int* NodesId, unsigned int PropertiesId )
    {
        return;
        
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );        

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Vector< std::vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = HONodesId[0]; Faces[0][1] = HONodesId[1]; Faces[0][2] = HONodesId[2]; Faces[0][3] = HONodesId[ 4]; Faces[0][4] = HONodesId[ 5]; Faces[0][5] = HONodesId[ 6];
        Faces[1].resize( 6 ); Faces[1][0] = HONodesId[3]; Faces[1][1] = HONodesId[1]; Faces[1][2] = HONodesId[0]; Faces[1][3] = HONodesId[ 7]; Faces[1][4] = HONodesId[ 8]; Faces[1][5] = HONodesId[ 9];
		Faces[2].resize( 6 ); Faces[2][0] = HONodesId[1]; Faces[2][1] = HONodesId[3]; Faces[2][2] = HONodesId[2]; Faces[2][3] = HONodesId[10]; Faces[2][4] = HONodesId[11]; Faces[2][5] = HONodesId[12];
		Faces[3].resize( 6 ); Faces[3][0] = HONodesId[2]; Faces[3][1] = HONodesId[3]; Faces[3][2] = HONodesId[0]; Faces[3][3] = HONodesId[13]; Faces[3][4] = HONodesId[14]; Faces[3][5] = HONodesId[15];
         
        for (int iface=0; iface<4 ; iface++ )
        {
            if ( mContactPairs.find( Faces[ iface ][ 3 ] ) != mContactPairs.end() ) 
            {
                int ctcNumNodes = 12;
                int ctcNumDofs  = ctcNumNodes*3; 

                std::complex<double> cZero( 0.00, 0.00 );
                
                std::vector<int> FHONodesId( ctcNumNodes );

                for ( int pnf=0; pnf<6; pnf++ ) 
                {
                    FHONodesId[ pnf     ] =                Faces[ iface ][ pnf ]           ;// R nodes
                    FHONodesId[ pnf + 6 ] = mContactPairs[ Faces[ iface ][ pnf ] ].ctcNode ;// L nodes
                }

                Vector< Node::Pointer > pctcNodes;

                std::vector<int>::iterator fhoit;

                for ( fhoit = FHONodesId.begin(); fhoit != FHONodesId.end(); ++fhoit ) pctcNodes.push_back( mpModel->GetNode( (*fhoit)-1 ) );      

                // Material properties matrix
				Properties::Pointer Properties;

				double eo = 8.8541878176e-12;
                double pi = 3.1415926535897932384626433832795;
                double mo = (4.00e-7) * pi;   

				double freq = mProblemFrequency;

                double sigma, eps_real, eps_imag;

                ContactPairData& pContactPair = mContactPairs[ Faces[ iface ][ 3 ] ];

				// Material R
				Properties = mpModel->GetProperties(pContactPair.matR);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

				// Material L
				Properties = mpModel->GetProperties(pContactPair.matL);

				sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
                eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		        std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

                double mu_real  = (*Properties)(REAL_MAGNETIC_PERMEABILITY) * mo;
                double mu_imag  = (*Properties)(IMAG_MAGNETIC_PERMEABILITY) * mo;

                std::complex<double> cMu( mu_real , mu_imag );

                // Getting equation id vector
                Vector<int> EleIdVector; EleIdVector.resize( ctcNumDofs );

                for( int i = 0; i < ctcNumNodes; ++i ) 
                {  
                    EleIdVector[ i                 ] = pctcNodes[ i ]->pDofcEx()->EquationId();
                    EleIdVector[ i + ctcNumNodes   ] = pctcNodes[ i ]->pDofcEy()->EquationId();
                	EleIdVector[ i + ctcNumNodes*2 ] = pctcNodes[ i ]->pDofcEz()->EquationId();
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////
 
                // Loop indexes
		        int i, j, gp;

                // Gauss points and weights
		        Vector<double> cX; 
		        Vector<double> cY;
		        Vector<double>  W;
		        
		        int nGaussPoints = GaussPoints2D_Order10( cX, cY, W );

                double narea[3],v1[3],v2[3];

                v2[0] = pctcNodes[ 2 ]->X() - pctcNodes[ 0 ]->X();
                v2[1] = pctcNodes[ 2 ]->Y() - pctcNodes[ 0 ]->Y();
                v2[2] = pctcNodes[ 2 ]->Z() - pctcNodes[ 0 ]->Z();
                
                v1[1] = pctcNodes[ 1 ]->Y() - pctcNodes[ 0 ]->Y();
                v1[2] = pctcNodes[ 1 ]->Z() - pctcNodes[ 0 ]->Z();
                v1[0] = pctcNodes[ 1 ]->X() - pctcNodes[ 0 ]->X();
                
                //area = 0.5*(v2 x v1)
                narea[0] = v2[1]*v1[2] - v2[2]*v1[1];
                narea[1] = v2[2]*v1[0] - v2[0]*v1[2];
                narea[2] = v2[0]*v1[1] - v2[1]*v1[0];
                
                double Area = 0.5*sqrt( narea[0]*narea[0] + narea[1]*narea[1] + narea[2]*narea[2] );

                ///////////////////////////////////////////////////////////////

                //Vector<double> ctcN = pContactPair.ctcNormal;

                //Vector< Vector<double> > v_ctcN; v_ctcN.resize(6); for ( j= 0; j<6; j++ ) v_ctcN[j] = mContactPairs[ Faces[ iface ][ j ] ].ctcNormal;
                
                Vector<double> ctcN; ctcN.resize( 3 );

                ctcN[0] = v2[1]*v1[2] - v2[2]*v1[1];
                ctcN[1] = v2[2]*v1[0] - v2[0]*v1[2];
                ctcN[2] = v2[0]*v1[1] - v2[1]*v1[0];

                double inv_extN_norm = 1.00 / sqrt( ctcN[0]*ctcN[0] + ctcN[1]*ctcN[1] + ctcN[2]*ctcN[2] );

		        ctcN *= inv_extN_norm;

                //////////////////////////////////////////////////////

                // Jacobian
		        double detJ = 2.00 * Area;
		        
		        for( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= detJ; 
		        
                // Lagrangian base on Gauss points
		        Matrix<double> N;
		        
		        Lagrange2D_Ni_2bb( N, cX, cY );

                ///////////////////////////////////////////////////////////////////

                std::vector<Node::Pointer> cpNodes; 
                
                for(int i=0; i<4; i++) cpNodes.push_back( mpModel->GetNode( HONodesId[i]-1 ) );

                double X1 = cpNodes[0]->X(); double X2 = cpNodes[1]->X(); double X3 = cpNodes[2]->X(); double X4 = cpNodes[3]->X(); 
                double Y1 = cpNodes[0]->Y(); double Y2 = cpNodes[1]->Y(); double Y3 = cpNodes[2]->Y(); double Y4 = cpNodes[3]->Y(); 
                double Z1 = cpNodes[0]->Z(); double Z2 = cpNodes[1]->Z(); double Z3 = cpNodes[2]->Z(); double Z4 = cpNodes[3]->Z(); 

                double cVolume = + X2*Y3*Z4 + X4*Y2*Z3 + X3*Y4*Z2 - X4*Y3*Z2 - X2*Y4*Z3 - X3*Y2*Z4 
                                 - X1*Y3*Z4 - X4*Y1*Z3 - X3*Y4*Z1 + X4*Y3*Z1 + X1*Y4*Z3 + X3*Y1*Z4 
                                 + X1*Y2*Z4 + X4*Y1*Z2 + X2*Y4*Z1 - X4*Y2*Z1 - X1*Y4*Z2 - X2*Y1*Z4 
                                 - X1*Y2*Z3 - X3*Y1*Z2 - X2*Y3*Z1 + X3*Y2*Z1 + X1*Y3*Z2 + X2*Y1*Z3;
                 
                cVolume = fabs( cVolume / 6.00 );

                Matrix<double> DN; 
                
                // [ Coord: X,Y,Z ][ Local_node: 0,1,2,3 ]
                DN.Resize( 3 , 4 );
                
                double cte = 1.00 / ( 6.00 * cVolume ); 
                
                // dN/dx
                DN[0][0] = cte * ( Y3*Z2 + Y4*Z3 + Y2*Z4 - Y3*Z4 - Y2*Z3 - Y4*Z2 );
                DN[0][1] = cte * ( Y3*Z4 + Y4*Z1 + Y1*Z3 - Y3*Z1 - Y1*Z4 - Y4*Z3 );
                DN[0][2] = cte * ( Y2*Z1 + Y1*Z4 + Y4*Z2 - Y2*Z4 - Y4*Z1 - Y1*Z2 );
                DN[0][3] = cte * ( Y2*Z3 + Y3*Z1 + Y1*Z2 - Y2*Z1 - Y3*Z2 - Y1*Z3 );
                                 
                // dN/dy         
                DN[1][0] = cte * ( X3*Z4 + X4*Z2 + X2*Z3 - X3*Z2 - X2*Z4 - X4*Z3 ); 
                DN[1][1] = cte * ( X3*Z1 + X1*Z4 + X4*Z3 - X1*Z3 - X3*Z4 - X4*Z1 );
                DN[1][2] = cte * ( X2*Z4 + X4*Z1 + X1*Z2 - X2*Z1 - X4*Z2 - X1*Z4 );
                DN[1][3] = cte * ( X1*Z3 + X3*Z2 + X2*Z1 - X3*Z1 - X1*Z2 - X2*Z3 );
                                 
                // dN/dz         
                DN[2][0] = cte * ( X4*Y3 + X2*Y4 + X3*Y2 - X2*Y3 - X4*Y2 - X3*Y4 );
                DN[2][1] = cte * ( X3*Y4 + X4*Y1 + X1*Y3 - X3*Y1 - X4*Y3 - X1*Y4 ); 
                DN[2][2] = cte * ( X2*Y1 + X4*Y2 + X1*Y4 - X4*Y1 - X2*Y4 - X1*Y2 ); 
                DN[2][3] = cte * ( X2*Y3 + X3*Y1 + X1*Y2 - X2*Y1 - X3*Y2 - X1*Y3 );  


                Vector<double> dNjdx( ctcNumNodes, 0.00 ); 
                Vector<double> dNjdy( ctcNumNodes, 0.00 ); 
                Vector<double> dNjdz( ctcNumNodes, 0.00 );

                Vector< Vector<int> > cFaces; cFaces.resize( 4 );

                cFaces[0].resize( 6 ); cFaces[0][0] = 0; cFaces[0][1] = 1; cFaces[0][2] = 2; cFaces[0][3] =  4; cFaces[0][4] =  5; cFaces[0][5] =  6;
                cFaces[1].resize( 6 ); cFaces[1][0] = 3; cFaces[1][1] = 1; cFaces[1][2] = 0; cFaces[1][3] =  7; cFaces[1][4] =  8; cFaces[1][5] =  9;
		        cFaces[2].resize( 6 ); cFaces[2][0] = 1; cFaces[2][1] = 3; cFaces[2][2] = 2; cFaces[2][3] = 10; cFaces[2][4] = 11; cFaces[2][5] = 12;
		        cFaces[3].resize( 6 ); cFaces[3][0] = 2; cFaces[3][1] = 3; cFaces[3][2] = 0; cFaces[3][3] = 13; cFaces[3][4] = 14; cFaces[3][5] = 15;

                for( int dj=0; dj<4; dj++ ) 
                {
                    dNjdx[dj    ] = DN[0][cFaces[iface][dj]];
                    dNjdy[dj    ] = DN[1][cFaces[iface][dj]];
                    dNjdz[dj    ] = DN[2][cFaces[iface][dj]];

                    dNjdx[dj + 6] = DN[0][cFaces[iface][dj]];
                    dNjdy[dj + 6] = DN[1][cFaces[iface][dj]];
                    dNjdz[dj + 6] = DN[2][cFaces[iface][dj]]; 
                }

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////

                Matrix< std::complex<double> > C_div; C_div.Resize( 3, ctcNumDofs, cZero );
                Matrix< std::complex<double> > C_cur; C_cur.Resize( 9, ctcNumDofs, cZero );
                Matrix< std::complex<double> > D_cur; D_cur.Resize( 9, ctcNumDofs, cZero );

                for( i=0; i<3; i++ )
		        {
		            for( j=0; j<6; j++ )
		            {
		            	std::complex<double> Sintg_NiNj( 0.00, 0.00 );
		                
                        for( gp=0; gp<nGaussPoints; gp++ ) 
			    		{
			    			Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
			    		}
		            		
                        C_div[ i     ][ j                     ] += ( cEpsR * ctcN[0] * Sintg_NiNj ); 
                        C_div[ i     ][ j + ctcNumNodes       ] += ( cEpsR * ctcN[1] * Sintg_NiNj ); 
                        C_div[ i     ][ j + ctcNumNodes*2     ] += ( cEpsR * ctcN[2] * Sintg_NiNj );
                                     
                        C_div[ i     ][ j +                 6 ] -= ( cEpsL * ctcN[0] * Sintg_NiNj ); 	  
                        C_div[ i     ][ j + ctcNumNodes   + 6 ] -= ( cEpsL * ctcN[1] * Sintg_NiNj ); 	  
                        C_div[ i     ][ j + ctcNumNodes*2 + 6 ] -= ( cEpsL * ctcN[2] * Sintg_NiNj ); 	  

                        C_cur[ i     ][ j + ctcNumNodes       ] += (  ctcN[2] * Sintg_NiNj ); 
                        C_cur[ i     ][ j + ctcNumNodes*2     ] -= (  ctcN[1] * Sintg_NiNj );  	 					 
		        	    C_cur[ i + 3 ][ j                     ] -= (  ctcN[2] * Sintg_NiNj ); 
                        C_cur[ i + 3 ][ j + ctcNumNodes*2     ] += (  ctcN[0] * Sintg_NiNj );		 	 						
		        	    C_cur[ i + 6 ][ j                     ] += (  ctcN[1] * Sintg_NiNj ); 
                        C_cur[ i + 6 ][ j + ctcNumNodes       ] -= (  ctcN[0] * Sintg_NiNj ); 

                        C_cur[ i     ][ j + ctcNumNodes   + 6 ] += ( -ctcN[2] * Sintg_NiNj ); 
                        C_cur[ i     ][ j + ctcNumNodes*2 + 6 ] -= ( -ctcN[1] * Sintg_NiNj );  	 					 
		        	    C_cur[ i + 3 ][ j                 + 6 ] -= ( -ctcN[2] * Sintg_NiNj ); 
                        C_cur[ i + 3 ][ j + ctcNumNodes*2 + 6 ] += ( -ctcN[0] * Sintg_NiNj );		 	 						
		        	    C_cur[ i + 6 ][ j                 + 6 ] += ( -ctcN[1] * Sintg_NiNj ); 
                        C_cur[ i + 6 ][ j + ctcNumNodes   + 6 ] -= ( -ctcN[0] * Sintg_NiNj ); 

                        //////////////////////////////////////////////////////////////////////////////////////

 		    	        double sintg_Ni_dNjdX = 0.00;
		                double sintg_Ni_dNjdY = 0.00;
		                double sintg_Ni_dNjdZ = 0.00;
		                
                        for( gp=0; gp<nGaussPoints; gp++ )
                        {
                        	sintg_Ni_dNjdX += W[gp] * N[i][gp] * dNjdx[j]; 
		    	        	sintg_Ni_dNjdY += W[gp] * N[i][gp] * dNjdy[j];
		    	        	sintg_Ni_dNjdZ += W[gp] * N[i][gp] * dNjdz[j];
                        }      

                        D_cur[ i     ][ j + ctcNumNodes       ] += sintg_Ni_dNjdZ; 
                        D_cur[ i     ][ j + ctcNumNodes*2     ] -= sintg_Ni_dNjdY; 
		        	    D_cur[ i + 3 ][ j                     ] -= sintg_Ni_dNjdZ; 
                        D_cur[ i + 3 ][ j + ctcNumNodes*2     ] += sintg_Ni_dNjdX;
		        	    D_cur[ i + 6 ][ j                     ] += sintg_Ni_dNjdY; 
                        D_cur[ i + 6 ][ j + ctcNumNodes       ] -= sintg_Ni_dNjdX; 
                       
                        D_cur[ i     ][ j + ctcNumNodes   + 6 ] += sintg_Ni_dNjdZ;  
                        D_cur[ i     ][ j + ctcNumNodes*2 + 6 ] -= sintg_Ni_dNjdY;  
		        	    D_cur[ i + 3 ][ j                 + 6 ] -= sintg_Ni_dNjdZ;  
                        D_cur[ i + 3 ][ j + ctcNumNodes*2 + 6 ] += sintg_Ni_dNjdX;
		        	    D_cur[ i + 6 ][ j                 + 6 ] += sintg_Ni_dNjdY;  
                        D_cur[ i + 6 ][ j + ctcNumNodes   + 6 ] -= sintg_Ni_dNjdX;  
		            }
		        }

                Matrix< std::complex<double> > iA_div; iA_div.Resize( 3, 3, cZero );
                Matrix< std::complex<double> > iA_cur; iA_cur.Resize( 9, 9, cZero );

                for( i=0; i<3; i++ )
                {
                    for( j=0; j<3; j++ )
                    {
                		if ( i != j ) iA_div[ i ][ j ] = -3.00 / Area;
                		else          iA_div[ i ][ j ] = +9.00 / Area;                          
                    
                        iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                        iA_cur[ i + 3 ][ j + 3 ] = iA_div[ i ][ j ];
                        iA_cur[ i + 6 ][ j + 6 ] = iA_div[ i ][ j ];
                    }
                }

                iA_div /= ( cMu*cEpsR*cEpsL  );
                iA_cur /= ( cMu   );

                // Getting element stiffness matrix
                Matrix< std::complex<double> > EleStiffMatrix; 
                EleStiffMatrix.Resize( ctcNumDofs, ctcNumDofs, cZero );

                Matrix< std::complex<double> > iAd_x_Ddiv = iA_div * C_div;
		        Matrix< std::complex<double> > Cdiv_tran  = C_div.Transpose();
                EleStiffMatrix += ( Cdiv_tran * iAd_x_Ddiv );

                Matrix< std::complex<double> > iAc_x_Ccur = iA_cur * C_cur;
		        Matrix< std::complex<double> > Ccur_tran  = C_cur.Transpose();
		        EleStiffMatrix += ( Ccur_tran * iAc_x_Ccur );

                Matrix< std::complex<double> > iAc_x_Dcur = iA_cur * D_cur;
                //EleStiffMatrix -= ( Ccur_tran * iAc_x_Dcur );

                Matrix< std::complex<double> > Dcur_tran  = D_cur.Transpose();
                //EleStiffMatrix -= ( Dcur_tran * iAc_x_Ccur );


                /////////////////////////////////////////////////////////////////////////////////////////////////

                double hX1 = pctcNodes[0]->X(), hY1 = pctcNodes[0]->Y(), hZ1 = pctcNodes[0]->Z();
		        double hX2 = pctcNodes[1]->X(), hY2 = pctcNodes[1]->Y(), hZ2 = pctcNodes[1]->Z();
		        double hX3 = pctcNodes[2]->X(), hY3 = pctcNodes[2]->Y(), hZ3 = pctcNodes[2]->Z();
		        
		        // Element edge lenghts 
		        double Lenght_12 = sqrt( (hX1-hX2)*(hX1-hX2) + (hY1-hY2)*(hY1-hY2) + (hZ1-hZ2)*(hZ1-hZ2) );
		        double Lenght_13 = sqrt( (hX1-hX3)*(hX1-hX3) + (hY1-hY3)*(hY1-hY3) + (hZ1-hZ3)*(hZ1-hZ3) );
		        double Lenght_23 = sqrt( (hX2-hX3)*(hX2-hX3) + (hY2-hY3)*(hY2-hY3) + (hZ2-hZ3)*(hZ2-hZ3) );
		        
                // Maximum lenght 
		        double hF = Lenght_12;
		        
		        if ( Lenght_13 > hF ) hF = Lenght_13;
		        if ( Lenght_23 > hF ) hF = Lenght_23;
                
			    std::complex<double> cteNxN; 

                ////////////////////////////////////////////////////////////////////////////////////////////////

                Matrix< std::complex<double> > D_dif; D_dif.Resize( ctcNumDofs, ctcNumDofs, cZero );

                cteNxN = hF / ( cMu*cEpsR*cEpsL  );

                for( i=0; i<6; i++ )
		        {
		            for( j=0; j<6; j++ )
		            {
		            	std::complex<double> Sintg_NiNj( 0.00, 0.00 );
		                
                        for( gp=0; gp<nGaussPoints; gp++ ) 
			    		{
			    			Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
			    		}

                        Sintg_NiNj *= cteNxN;
		            	
                        D_dif[ i      ][ j     ] = +ctcN[0] * ctcN[0] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j     ] = -ctcN[0] * ctcN[0] * cEpsR * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j     ] = +ctcN[1] * ctcN[0] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j     ] = -ctcN[1] * ctcN[0] * cEpsR * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j     ] = +ctcN[2] * ctcN[0] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j     ] = -ctcN[2] * ctcN[0] * cEpsR * cEpsL * Sintg_NiNj; 


                        D_dif[ i      ][ j + 6 ] = -ctcN[0] * ctcN[0] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j + 6 ] = +ctcN[0] * ctcN[0] * cEpsL * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j + 6 ] = -ctcN[1] * ctcN[0] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j + 6 ] = +ctcN[1] * ctcN[0] * cEpsL * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j + 6 ] = -ctcN[2] * ctcN[0] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j + 6 ] = +ctcN[2] * ctcN[0] * cEpsL * cEpsL * Sintg_NiNj; 


                        D_dif[ i      ][ j + 12 ] = +ctcN[0] * ctcN[1] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j + 12 ] = -ctcN[0] * ctcN[1] * cEpsR * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j + 12 ] = +ctcN[1] * ctcN[1] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j + 12 ] = -ctcN[1] * ctcN[1] * cEpsR * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j + 12 ] = +ctcN[2] * ctcN[1] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j + 12 ] = -ctcN[2] * ctcN[1] * cEpsR * cEpsL * Sintg_NiNj; 


                        D_dif[ i      ][ j + 18 ] = -ctcN[0] * ctcN[1] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j + 18 ] = +ctcN[0] * ctcN[1] * cEpsL * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j + 18 ] = -ctcN[1] * ctcN[1] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j + 18 ] = +ctcN[1] * ctcN[1] * cEpsL * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j + 18 ] = -ctcN[2] * ctcN[1] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j + 18 ] = +ctcN[2] * ctcN[1] * cEpsL * cEpsL * Sintg_NiNj; 


                        D_dif[ i      ][ j + 24 ] = +ctcN[0] * ctcN[2] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j + 24 ] = -ctcN[0] * ctcN[2] * cEpsR * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j + 24 ] = +ctcN[1] * ctcN[2] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j + 24 ] = -ctcN[1] * ctcN[2] * cEpsR * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j + 24 ] = +ctcN[2] * ctcN[2] * cEpsR * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j + 24 ] = -ctcN[2] * ctcN[2] * cEpsR * cEpsL * Sintg_NiNj; 


                        D_dif[ i      ][ j + 30 ] = -ctcN[0] * ctcN[2] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i +  6 ][ j + 30 ] = +ctcN[0] * ctcN[2] * cEpsL * cEpsL * Sintg_NiNj;  
                        D_dif[ i + 12 ][ j + 30 ] = -ctcN[1] * ctcN[2] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 18 ][ j + 30 ] = +ctcN[1] * ctcN[2] * cEpsL * cEpsL * Sintg_NiNj; 
                        D_dif[ i + 24 ][ j + 30 ] = -ctcN[2] * ctcN[2] * cEpsL * cEpsR * Sintg_NiNj; 
                        D_dif[ i + 30 ][ j + 30 ] = +ctcN[2] * ctcN[2] * cEpsL * cEpsL * Sintg_NiNj; 
		            }
		        }

                EleStiffMatrix += D_dif;

                /////////////////////////////////////////////////////////////////////////////////////////

                Matrix< std::complex<double> > C_dif; C_dif.Resize( ctcNumDofs, ctcNumDofs, cZero );

                cteNxN = hF / cMu;

                for( i=0; i<6; i++ )
		        {
		            for( j=0; j<6; j++ )
		            {
		            	std::complex<double> Sintg_NiNj( 0.00, 0.00 );
		                
                        for( gp=0; gp<nGaussPoints; gp++ ) 
			    		{
			    			Sintg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
			    		}

                        Sintg_NiNj *= cteNxN;
		            	
                        C_dif[ i      ][ j     ] = +( ctcN[2] * ctcN[2] + ctcN[1] * ctcN[1] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j     ] = -( ctcN[2] * ctcN[2] + ctcN[1] * ctcN[1] ) * Sintg_NiNj; 
                        C_dif[ i + 12 ][ j     ] = -(                     ctcN[1] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j     ] = +(                     ctcN[1] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j     ] = -( ctcN[2] * ctcN[0]                     ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j     ] = +( ctcN[2] * ctcN[0]                     ) * Sintg_NiNj; 
                        
                        C_dif[ i      ][ j + 6 ] = -( ctcN[2] * ctcN[2] + ctcN[1] * ctcN[1] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j + 6 ] = +( ctcN[2] * ctcN[2] + ctcN[1] * ctcN[1] ) * Sintg_NiNj;  
                        C_dif[ i + 12 ][ j + 6 ] = +(                     ctcN[1] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j + 6 ] = -(                     ctcN[1] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j + 6 ] = +( ctcN[2] * ctcN[0]                     ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j + 6 ] = -( ctcN[2] * ctcN[0]                     ) * Sintg_NiNj; 
                        
                        C_dif[ i      ][ j + 12 ] = -(                     ctcN[0] * ctcN[1] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j + 12 ] = +(                     ctcN[0] * ctcN[1] ) * Sintg_NiNj;  
                        C_dif[ i + 12 ][ j + 12 ] = +( ctcN[2] * ctcN[2] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j + 12 ] = -( ctcN[2] * ctcN[2] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j + 12 ] = -( ctcN[2] * ctcN[1]                     ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j + 12 ] = +( ctcN[2] * ctcN[1]                     ) * Sintg_NiNj; 
                        
                        C_dif[ i      ][ j + 18 ] = +(                     ctcN[0] * ctcN[1] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j + 18 ] = -(                     ctcN[0] * ctcN[1] ) * Sintg_NiNj;  
                        C_dif[ i + 12 ][ j + 18 ] = -( ctcN[2] * ctcN[2] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j + 18 ] = +( ctcN[2] * ctcN[2] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j + 18 ] = +( ctcN[2] * ctcN[1]                     ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j + 18 ] = -( ctcN[2] * ctcN[1]                     ) * Sintg_NiNj; 
                       
                        C_dif[ i      ][ j + 24 ] = -(                     ctcN[0] * ctcN[2] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j + 24 ] = +(                     ctcN[0] * ctcN[2] ) * Sintg_NiNj; 
                        C_dif[ i + 12 ][ j + 24 ] = -( ctcN[1] * ctcN[2]                     ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j + 24 ] = +( ctcN[1] * ctcN[2]                     ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j + 24 ] = +( ctcN[1] * ctcN[1] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j + 24 ] = -( ctcN[1] * ctcN[1] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 

                        C_dif[ i      ][ j + 30 ] = +(                     ctcN[0] * ctcN[2] ) * Sintg_NiNj; 
                        C_dif[ i +  6 ][ j + 30 ] = -(                     ctcN[0] * ctcN[2] ) * Sintg_NiNj;  
                        C_dif[ i + 12 ][ j + 30 ] = +( ctcN[1] * ctcN[2]                     ) * Sintg_NiNj; 
                        C_dif[ i + 18 ][ j + 30 ] = -( ctcN[1] * ctcN[2]                     ) * Sintg_NiNj; 
                        C_dif[ i + 24 ][ j + 30 ] = -( ctcN[1] * ctcN[1] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
                        C_dif[ i + 30 ][ j + 30 ] = +( ctcN[1] * ctcN[1] + ctcN[0] * ctcN[0] ) * Sintg_NiNj; 
		            }
		        }

                EleStiffMatrix += C_dif;

                /////////////////////////////////////////////////////////////////////////////////////////

                //Trans(T)*A*T for PBC conditions
		        PBC_FrontBack_TKT( FHONodesId, EleIdVector, EleStiffMatrix );
	            PBC_RightLeft_TKT( FHONodesId, EleIdVector, EleStiffMatrix );
		        PBC_FrontBack_TKT( FHONodesId, EleIdVector, EleStiffMatrix );
                
		        //Trans(T)*A*T for contact discontinuities
		        ContactHO_TKT( FHONodesId, EleIdVector, EleStiffMatrix );
                
                //Trans(T)*A*T for PEC, PMC conditions
                HO_TKT( FHONodesId, EleStiffMatrix );
              
		        //Trans(T)*A*T for axisymmetric problems
		        Axisym_TKT( FHONodesId, EleStiffMatrix );

                //Assembling to global matrix
                int eleSize = EleIdVector.size();
                
                std::complex<double> cKij;
                
                // Assembling symmetric matrix
		        for ( int i_local = 0; i_local < eleSize; i_local++ )
		        {
		        	int i_global = EleIdVector[i_local];
                
		        	if ( i_global < mSystemSize )
		        	{
		        		for ( int j_local = 0; j_local < eleSize; j_local++ )
		        		{
		        			int j_global = EleIdVector[j_local];
                
		        			if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                            {
                                cKij = EleStiffMatrix[i_local][j_local];
                                if ( abs( cKij ) > 0.0 ) A_matrix( i_global, j_global ) += cKij;
                            }
		        		}
		        	}
                }

                EleIdVector.FreeData();
                EleStiffMatrix.FreeData();
            }
        }
    }

    //*************************************************************************************************
    //* - Setting element vector
    //*************************************************************************************************
    void Modeler::Set_Elements_Vector( int* NodesId, unsigned int PropertiesId )
    {
        return;

  //      if ( mElementOrder != 0 ) return;
  //      
  //      std::vector<int> HONodesId;

  //      PushHONodesOnVolume( NodesId, HONodesId );

  //      Node::Pointer pNode;

  //      std::complex<double> cZero(0.00,0.00);

  //      pNode = mpModel->GetNode( HONodesId[16]-1 );
  //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();


  //      pNode = mpModel->GetNode( HONodesId[4]-1 );
  //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //      pNode = mpModel->GetNode( HONodesId[7]-1 );
  //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //      pNode = mpModel->GetNode( HONodesId[10]-1 );
  //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //      pNode = mpModel->GetNode( HONodesId[13]-1 );
  //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

 /*       for (int i=4;i<22;i++)
        {
            pNode = mpModel->GetNode( HONodesId[i]-1 );
            (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	        (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
        } */

        //return;       

        /*
        std::vector<int> HONodesId; 
        
        PushHONodesOnVolume( NodesId, HONodesId );

        Vector< Vector<int> > Faces; Faces.resize( 4 );

		// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
        Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

        for (int iFace = 0; iFace<4 ; iFace++ )
        {
            int ctcNodeId = HONodesId[ Faces[iFace][3] ];
            
            if ( mContactPairs.find( ctcNodeId ) != mContactPairs.end() ) 
            {
                int ctcNodeId1 = HONodesId[ Faces[iFace][3] ];

                mType_Of_BC_Normal[ ctcNodeId1 ] = 'C';
                mType_Of_BC_Normal[ mContactPairs[ ctcNodeId1 ].ctcNode ] = 'R';

                int ctcNodeId2 = HONodesId[ Faces[iFace][4] ];

                mType_Of_BC_Normal[ ctcNodeId2 ] = 'C';
                mType_Of_BC_Normal[ mContactPairs[ ctcNodeId2 ].ctcNode ] = 'R';

                int ctcNodeId3 = HONodesId[ Faces[iFace][5] ];

                mType_Of_BC_Normal[ ctcNodeId3 ] = 'C';
                mType_Of_BC_Normal[ mContactPairs[ ctcNodeId3 ].ctcNode ] = 'R';
            }
        }
        */

        //return;
        
        //std::cout<<"Eliminating bubbles..."<<std::endl;

        ///////////////////////////////////////////////////////////
        //// Element info vector
        //Vector<unsigned int> ElementInfo( 5 );
        //
        //// Adding node Ids and material Id to the ElementInfo vector 
        //ElementInfo[0] = NodesId[0]-1;
        //ElementInfo[1] = NodesId[1]-1;
        //ElementInfo[2] = NodesId[2]-1;
        //ElementInfo[3] = NodesId[3]-1;
        //ElementInfo[4] = PropertiesId;
        //
        //// Adding element info to elements vector 
        //mElementsVector.push_back( ElementInfo );
        ////////////////////////////////////////////////////////

  //      std::vector<int> HONodesId; 
  //      
  //      PushHONodesOnVolume( NodesId, HONodesId );

  //      Vector< Vector<int> > Faces; Faces.resize( 4 );

		//// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		//Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
  //      Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		//Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		//Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

  //      Node::Pointer pNode;

  //      for (int iFace = 0; iFace<4 ; iFace++ )
  //      {
  //          int ctcNodeId = HONodesId[ Faces[iFace][3] ];
  //          
  //          if ( (mContactPairs.find( ctcNodeId ) != mContactPairs.end()) || (mType_Of_BC_Normal[ mContactPairs[ ctcNodeId ].ctcNode ] == 'L') ) 
  //          {
  //              Node::Pointer pNode;
  //              std::complex<double> cZero(0.00,0.00);

  //              for (int i=4;i<17;i++)
  //              {
  //                  pNode = mpModel->GetNode(HONodesId[i]-1);
  //                  (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //               (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//            (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
  //              } 
  //          }
  //      }

  //      return;

  //      std::vector<int> HONodesId; 
  //      
  //      PushHONodesOnVolume( NodesId, HONodesId );

  //      std::complex<double> cZero(0.00,0.00);

  //      Vector< Vector<int> > Faces; Faces.resize( 4 );

		//// Tetrahedra: [1-2-3] [2-4-3] [3-4-1] [4-2-1]
		//Faces[0].resize( 6 ); Faces[0][0] = 0; Faces[0][1] = 1; Faces[0][2] = 2; Faces[0][3] =  4; Faces[0][4] =  5; Faces[0][5] =  6;
  //      Faces[1].resize( 6 ); Faces[1][0] = 3; Faces[1][1] = 1; Faces[1][2] = 0; Faces[1][3] =  7; Faces[1][4] =  8; Faces[1][5] =  9;
		//Faces[2].resize( 6 ); Faces[2][0] = 1; Faces[2][1] = 3; Faces[2][2] = 2; Faces[2][3] = 10; Faces[2][4] = 11; Faces[2][5] = 12;
		//Faces[3].resize( 6 ); Faces[3][0] = 2; Faces[3][1] = 3; Faces[3][2] = 0; Faces[3][3] = 13; Faces[3][4] = 14; Faces[3][5] = 15;

  //      Node::Pointer pNode;

  //      for (int iFace = 0; iFace<4 ; iFace++ )
  //      {
  //          int ctcNodeId = HONodesId[ Faces[iFace][3] ];
  //          
  //          if ( mContactPairs.find( ctcNodeId ) != mContactPairs.end() ) 
  //          {
  //              pNode = mpModel->GetNode(HONodesId[ Faces[iFace][3] ]-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode(HONodesId[ Faces[iFace][4] ]-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode(HONodesId[ Faces[iFace][5] ]-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode( mContactPairs[ HONodesId[ Faces[iFace][3] ] ].ctcNode-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode( mContactPairs[ HONodesId[ Faces[iFace][4] ] ].ctcNode-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode( mContactPairs[ HONodesId[ Faces[iFace][5] ] ].ctcNode-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();

  //              pNode = mpModel->GetNode( HONodesId[ 16 ]-1 );
  //              (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	 //           (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		//        (*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
  //          }
  //      }

    }

    //*************************************************************************************************
    //* - Setting nodal connectivity vector
    //*************************************************************************************************
    void Modeler::Set_Nodal_Connectivities()
    {
        mNodesConnectivity.resize( mpModel->GetNodesArray().size() );        
               
        for( int ElementId = 0; ElementId < mElementsVector.size(); ++ElementId ) 
		{
            Vector<unsigned int> ElementInfo = mElementsVector[ ElementId ];
            
            mNodesConnectivity[ ElementInfo[0] ].push_back( ElementId );
            mNodesConnectivity[ ElementInfo[1] ].push_back( ElementId );
            mNodesConnectivity[ ElementInfo[2] ].push_back( ElementId );
            mNodesConnectivity[ ElementInfo[3] ].push_back( ElementId );

            mPropertiesMap[ ElementInfo[4] ] = mpModel->GetProperties( ElementInfo[4] );
        }

        int CornerNodesCount = 0;

        for( int NodeId = 0; NodeId < mNodesConnectivity.size(); ++NodeId ) 
        {
            if ( mNodesConnectivity[ NodeId ].size() > 0 ) CornerNodesCount++;
        }

        mNodesConnectivity.resize( CornerNodesCount );
    }

    //*************************************************************************************************
    //* - Calculates 1st order derivatives and volumes for all elements in bubble mode
    //*************************************************************************************************
    void Modeler::Set_Bubble_Derivatives()
    {
        int numElements = mElementsVector.size();

        mbVolumes.resize( numElements );
        mb_dNk_dX.resize( numElements );
        mb_dNk_dY.resize( numElements );
        mb_dNk_dZ.resize( numElements ); 

        // Loop over all the elements in the mesh
		for ( int ElementId = 0; ElementId < numElements; ++ElementId )
		{
            Vector<unsigned int> ElementInfo = mElementsVector[ ElementId ];

            // Element nodes 
            Node::Pointer pNode_1 = mpModel->GetNode( ElementInfo[0] );
            Node::Pointer pNode_2 = mpModel->GetNode( ElementInfo[1] );
            Node::Pointer pNode_3 = mpModel->GetNode( ElementInfo[2] );
            Node::Pointer pNode_4 = mpModel->GetNode( ElementInfo[3] );

			// Nodes coordinates
		    double X1 = pNode_1->X(), Y1 = pNode_1->Y(), Z1 = pNode_1->Z();
		    double X2 = pNode_2->X(), Y2 = pNode_2->Y(), Z2 = pNode_2->Z();
		    double X3 = pNode_3->X(), Y3 = pNode_3->Y(), Z3 = pNode_3->Z();
		    double X4 = pNode_4->X(), Y4 = pNode_4->Y(), Z4 = pNode_4->Z();

			// Computing element volume
		    double det = + X2*Y3*Z4 + X4*Y2*Z3 + X3*Y4*Z2 - X4*Y3*Z2 - X2*Y4*Z3 - X3*Y2*Z4
                         - X1*Y3*Z4 - X4*Y1*Z3 - X3*Y4*Z1 + X4*Y3*Z1 + X1*Y4*Z3 + X3*Y1*Z4
                         + X1*Y2*Z4 + X4*Y1*Z2 + X2*Y4*Z1 - X4*Y2*Z1 - X1*Y4*Z2 - X2*Y1*Z4
                         - X1*Y2*Z3 - X3*Y1*Z2 - X2*Y3*Z1 + X3*Y2*Z1 + X1*Y3*Z2 + X2*Y1*Z3;

            double eVol = fabs( det / 6.00 );

			mbVolumes[ ElementId ] = eVol;

			// Computing derivatives
            mb_dNk_dX[ ElementId ].resize( 4 );
            mb_dNk_dY[ ElementId ].resize( 4 );
            mb_dNk_dZ[ ElementId ].resize( 4 );

            double cte = 1.00 / ( 6.00 * eVol ); 
		    
            // dNk / dx
            mb_dNk_dX[ ElementId ][ 0 ] = cte * ( Y3*Z2 + Y4*Z3 + Y2*Z4 - Y3*Z4 - Y2*Z3 - Y4*Z2 );
            mb_dNk_dX[ ElementId ][ 1 ] = cte * ( Y3*Z4 + Y4*Z1 + Y1*Z3 - Y3*Z1 - Y1*Z4 - Y4*Z3 );
            mb_dNk_dX[ ElementId ][ 2 ] = cte * ( Y2*Z1 + Y1*Z4 + Y4*Z2 - Y2*Z4 - Y4*Z1 - Y1*Z2 );
            mb_dNk_dX[ ElementId ][ 3 ] = cte * ( Y2*Z3 + Y3*Z1 + Y1*Z2 - Y2*Z1 - Y3*Z2 - Y1*Z3 );
                           		 
            // dNk / dy          
            mb_dNk_dY[ ElementId ][ 0 ] = cte * ( X3*Z4 + X4*Z2 + X2*Z3 - X3*Z2 - X2*Z4 - X4*Z3 ); 
            mb_dNk_dY[ ElementId ][ 1 ] = cte * ( X3*Z1 + X1*Z4 + X4*Z3 - X1*Z3 - X3*Z4 - X4*Z1 );
            mb_dNk_dY[ ElementId ][ 2 ] = cte * ( X2*Z4 + X4*Z1 + X1*Z2 - X2*Z1 - X4*Z2 - X1*Z4 );
            mb_dNk_dY[ ElementId ][ 3 ] = cte * ( X1*Z3 + X3*Z2 + X2*Z1 - X3*Z1 - X1*Z2 - X2*Z3 );
                              	 
            // dNk / dz          
            mb_dNk_dZ[ ElementId ][ 0 ] = cte * ( X4*Y3 + X2*Y4 + X3*Y2 - X2*Y3 - X4*Y2 - X3*Y4 );
            mb_dNk_dZ[ ElementId ][ 1 ] = cte * ( X3*Y4 + X4*Y1 + X1*Y3 - X3*Y1 - X4*Y3 - X1*Y4 ); 
            mb_dNk_dZ[ ElementId ][ 2 ] = cte * ( X2*Y1 + X4*Y2 + X1*Y4 - X4*Y1 - X2*Y4 - X1*Y2 ); 
            mb_dNk_dZ[ ElementId ][ 3 ] = cte * ( X2*Y3 + X3*Y1 + X1*Y2 - X2*Y1 - X3*Y2 - X1*Y3 );    
        }
    }

    //*************************************************************************************************
    //* - Cleaning derivatives and volumes vectors in bubble mode
    //*************************************************************************************************
    void Modeler::Clear_Bubble_Derivatives()
    {
        mbVolumes.clear(); Vector<        double  >().swap( mbVolumes );
        mb_dNk_dX.clear(); Vector< Vector<double> >().swap( mb_dNk_dX );
        mb_dNk_dY.clear(); Vector< Vector<double> >().swap( mb_dNk_dY ); 
        mb_dNk_dZ.clear(); Vector< Vector<double> >().swap( mb_dNk_dZ );
    }

    //*************************************************************************************************
    //* - Cleaning nodal connectivity objects
    //*************************************************************************************************
    void Modeler::Clear_Nodal_Connectivities()
    {
        // Cleaning mNodesConnectivity     
        mNodesConnectivity.clear(); 
        Vector< Vector<unsigned int> >().swap( mNodesConnectivity );    

        // Cleaning mElementsVector     
        mElementsVector.clear();
        Vector< Vector<unsigned int> >().swap( mElementsVector );    

        // Cleaning mPropertiesMap   
        mPropertiesMap.clear();
        std::map< unsigned int, Properties::Pointer >().swap( mPropertiesMap );    
    }

	//*************************************************************************************************
    //* - Calculates the maximum and the minimum distance to the singularity
    //*************************************************************************************************
	void Modeler::CalculateDistToSing( int* NodesId, std::vector<Node::Pointer>& nodes )
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
    void Modeler::GenerateSourceElement( int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        if ( mElementOrder == 4 ) 
        {
		    HONodesId.push_back( NodesId[0] );
            HONodesId.push_back( NodesId[1] );
            HONodesId.push_back( NodesId[2] );
		    HONodesId.push_back( NodesId[3] );

            HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[2], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[1], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[1], NodesId[2], NodesId[3], HONodesId );
		    HONodesOnFace  ( NodesId[0], NodesId[2], NodesId[3], HONodesId );
            HONodesOnVolume( NodesId   ,                         HONodesId );               
        }
        else
        {
            PushHONodesOnVolume( NodesId, HONodesId );
        }

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mQuadraticGeometry ) pElement = Element::Pointer( new JSource_2ndQ( nodes, properties ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new JSource_1bb ( nodes, properties ) );
	    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new JSource_1st ( nodes, properties ) );
	    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new JSource_2nd ( nodes, properties ) );
	    else if ( mElementOrder == 3 ) pElement = Element::Pointer( new JSource_3th ( nodes, properties ) );
	    else if ( mElementOrder == 4 ) pElement = Element::Pointer( new JSource_1bb ( nodes, properties ) );//pElement = Element::Pointer( new JSource_4th ( nodes, properties ) );
	    
        // Get DOF Ids vector
        Vector<int> EleIdVector;
        
        pElement->GetEquationIdVector( EleIdVector );

        // Get residual vector
        Vector<std::complex<double> > EleResVector;

        pElement->GetResidualVector( EleResVector );

		// Inv(T)*b for PBC condition
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );
	    PBC_RightLeft_TR( HONodesId, EleIdVector, EleResVector );
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );

		// Inv(T)*b for contact nodes
        if ( mPlasmaMode ) ContactHO_TR_Plasma( HONodesId, EleIdVector, EleResVector );
        else               ContactHO_TR       ( HONodesId, EleIdVector, EleResVector );

        // Inv(T)*b for PEC, PMC conditions
        HO_TR( HONodesId, EleResVector );

		// Inv(T)*b for axisymmetric mode
        if ( mAxisymmetric ) Axisym_TR( HONodesId, EleResVector );

        // Inv(T)*b for plasma mode
        if ( mPlasmaMode   ) Plasma_TR( HONodesId, EleResVector );

        //Assembling to global source vector
        int eleSize = EleIdVector.size();

		for( int i_local=0; i_local < eleSize; i_local++ )
		{
			int i_global = EleIdVector[i_local];

			if( i_global < mSystemSize )
			{
				b_vector[i_global] += EleResVector[i_local];
			}
		}

        EleIdVector.FreeData ();
        EleResVector.FreeData();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assembling boundary plasma element matrix on auxiliar global matrix  
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix )
    {
         int eleSize = EleIdVector.size();
         
         std::complex<double> cKij;
         
         double Kij_Tol = mpColdPlasma->Get_Kij_Tol();
         
         if ( mpColdPlasma->Is_HermSymm_Matrix() )
         {
             for (int i_local = 0; i_local < eleSize; i_local++)
             {
                 int i_global = EleIdVector[i_local];
         
                 if ( i_global < mSystemSize )
                 {
                     for (int j_local = 0; j_local < eleSize; j_local++)
                     {
                         int j_global = EleIdVector[j_local];
         
                         if ( (i_global <= j_global) && (j_global < mSystemSize) )
                         {
                             cKij = EleStiffMatrix[i_local][j_local];
                             if ( abs(cKij) > Kij_Tol ) A_matrix_aux(i_global,j_global) += cKij;
                         }
                     }
                 }
             }
         }
         else if ( mpColdPlasma->Is_HermFull_Matrix() )
         {
             for (int i_local = 0; i_local < eleSize; i_local++)
             {
                 int i_global = EleIdVector[i_local];
         
                 if ( i_global < mSystemSize )
                 {
                     for (int j_local = 0; j_local < eleSize; j_local++)
                     {
                         int j_global = EleIdVector[j_local];
         
                         if ( (j_global < mSystemSize) )
                         {
                             cKij = EleStiffMatrix[i_local][j_local];
                             if ( abs(cKij) > Kij_Tol ) A_matrix_aux(i_global,j_global) += cKij;
                         }
                     }
                 }
             }
         }
         else
         {
             for (int i_local = 0; i_local < eleSize; i_local++)
             {
                 int i_global = EleIdVector[i_local];
         
                 if ( i_global < mSystemSize )
                 {
                     for (int j_local = 0; j_local < eleSize; j_local++)
                     {
                         int j_global = EleIdVector[j_local];
         
                         if ( (j_global < mSystemSize) )
                         {
                             cKij = EleStiffMatrix[i_local][j_local];
                             if ( abs(cKij) > Kij_Tol ) A_matrix(i_global,j_global) += cKij;
                         }
                     }
                 }
             }
         }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assembling generic Robin elements in plasma mode
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::GenerateGenericRobinElement_Plasma( int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnSurface(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) nodes.push_back(mpModel->GetNode((*hoit)-1));

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;

        if ( mElementOrder == 1 ) pElement = Element::Pointer( new GenericRobin_1st_Plasma (nodes,properties) );
        else                      pElement = Element::Pointer( new GenericRobin_2ndQ_Plasma(nodes,properties) );

        pElement->Set_Plasma_Model( mpColdPlasma      );
        pElement->SetFrequency    ( mProblemFrequency );

        Vector<int> EleIdVector;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);

        //Trans(T)*A*T for PBC conditions
        PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
        PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
        PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for discontinuities
        ContactHO_TKT_Plasma(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for PEC and PMC
        HO_TKT(HONodesId, EleStiffMatrix);

        //Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        Plasma_TKT(HONodesId, EleStiffMatrix);

        //Trans(T)*A*T for axisymmetric problems
        Axisym_TKT(HONodesId, EleStiffMatrix);

        //Assembling to global matrix
        Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );

        //Cleaning elemental matrix
        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    // - Assembling far field elements in plasma mode
    ////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::GenerateFarFieldElement_Plasma( int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnSurface(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) nodes.push_back(mpModel->GetNode((*hoit)-1));

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;

        if ( mElementOrder == 1 ) pElement = Element::Pointer( new FarField_1st_Plasma      (nodes,properties) );
        else                      pElement = Element::Pointer( new FarField_2ndQ_DivD_Plasma(nodes,properties) );
        
        pElement->Set_Plasma_Model( mpColdPlasma      );
        pElement->SetFrequency    ( mProblemFrequency );
 
        Vector<int> EleIdVector;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);

        //Trans(T)*A*T for PBC conditions
        PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
        PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
        PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for discontinuities
        ContactHO_TKT_Plasma(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for PEC and PMC
        HO_TKT(HONodesId, EleStiffMatrix);

        //Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        Plasma_TKT(HONodesId, EleStiffMatrix);

        //Trans(T)*A*T for axisymmetric problems
        Axisym_TKT(HONodesId, EleStiffMatrix);

        //Assembling to global matrix
        Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );

        // Cleaning elemental matrix
        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//***************************************************************************************************************
    //* - Assembling far field elements
    //***************************************************************************************************************
    void Modeler::GenerateFarFieldElement( int* NodesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Element::Pointer pElement;
 
        if      ( mQuadraticGeometry ) pElement = Element::Pointer( new FarField_2ndQ( nodes ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new FarField_1bb ( nodes ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new FarField_1st ( nodes ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new FarField_2nd ( nodes ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new FarField_3th ( nodes ) );
        else if ( mElementOrder == 4 ) pElement = Element::Pointer( new FarField_4th ( nodes ) );

		pElement->SetFrequency(mProblemFrequency);

		Vector<int> EleIdVector;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector(EleIdVector);
        pElement->GetStiffMatrix(EleStiffMatrix);

        //Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);
	    PBC_RightLeft_TKT(HONodesId, EleIdVector, EleStiffMatrix);
		PBC_FrontBack_TKT(HONodesId, EleIdVector, EleStiffMatrix);

        //Trans(T)*A*T for contact discontinuities
        if ( mPlasmaMode == true ) ContactHO_TKT_Plasma(HONodesId, EleIdVector, EleStiffMatrix);
        else                       ContactHO_TKT       (HONodesId, EleIdVector, EleStiffMatrix);

		//Trans(T)*A*T
        HO_TKT(HONodesId, EleStiffMatrix);

        //Trans(T)*A*T for axisymmetric problems
        if ( mAxisymmetric == true ) Axisym_TKT(HONodesId, EleStiffMatrix);

		//Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        if ( mPlasmaMode   == true ) Plasma_TKT(HONodesId, EleStiffMatrix);

        //Assembling to global matrix 
        if ( mPlasmaMode == true )
        {
             Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
             int eleSize = EleIdVector.size();
             
             std::complex<double> cKij;
             
             //Assembling to symmetric matrix
             for (int i_local = 0; i_local < eleSize; i_local++)
             {
                  int i_global = EleIdVector[i_local];
                  
                  if ( i_global < mSystemSize )
                  {
                       for (int j_local = 0; j_local < eleSize; j_local++)
                       {
                       	    int j_global = EleIdVector[j_local];
                            
                       	    if ( (i_global <= j_global) && (j_global < mSystemSize) )
                            {
                                 cKij = EleStiffMatrix[i_local][j_local];
                                 if ( abs(cKij) > 0.0 ) A_matrix(i_global,j_global) += cKij;
                            }
                       }
                  }
             }
        }

        //Cleaning elemental matrix
        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
	}

	//***********************************************************************************************************
    //* - Assembling generic Robin elements
    //***********************************************************************************************************
    void Modeler::GenerateGenericRobinElement( int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;
		
        if      ( mQuadraticGeometry ) pElement = Element::Pointer( new GenericRobin_2ndQ( nodes, properties ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new GenericRobin_1st ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new GenericRobin_1st ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new GenericRobin_2nd ( nodes, properties ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new GenericRobin_3th ( nodes, properties ) );
        else if ( mElementOrder == 4 ) pElement = Element::Pointer( new GenericRobin_4th ( nodes, properties ) );

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
                        if ( abs(cKij) > 0.0 ) A_matrix(i_global,j_global) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
	}

    //*************************************************************************************************
    //* - Assembling generic Robin elements (Static mode)
    //*************************************************************************************************
    void Modeler::GenerateGenericRobinElement_Static( int* NodesId, unsigned int PropertiesId )
    {
        // Getting local nodes
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnSurface(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

        // Getting material properties
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        // Creating element
        Element::Pointer pElement;

        pElement = Element::Pointer( new GenericRobin_1st_Static(nodes, properties) );

        // Getting DOFs Id vector
        Vector<int> EleIdVector;

        pElement->GetEquationIdVector( EleIdVector );

        // Getting elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffMatrix( EleStiffMatrix );

        // Including Dirichlet BC in the residual vector
        Vector<double> EleResVector;

        pElement->GetResidualVector( EleResVector );

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

        // Assembling matrix
        for (int i_local = 0; i_local < eleSize; i_local++)
        {
            int i_global = EleIdVector[i_local];

            if ( i_global < mSystemSize )
            {
                for (int j_local = 0; j_local < eleSize; j_local++)
                {
                    int j_global = EleIdVector[j_local];

                    if ( (i_global <= j_global) && (j_global < mSystemSize) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
                        if ( abs(cKij) > 0.0 )  A_matrix(i_global,j_global) += cKij;
                    }
                }
            }
        }

        //Assembling global source vector
        for ( int i_blocal=0; i_blocal < eleSize; i_blocal++ )
        {
            int i_bglobal = EleIdVector[i_blocal];

            if ( i_bglobal < mSystemSize )
            {
                b_vector[i_bglobal] += EleResVector[i_blocal];
            }
        }

        EleIdVector.FreeData   ();
        EleResVector.FreeData  ();
        EleStiffMatrix.FreeData();
    }

	//********************************************************************************************************
    //* - Assembling boundary elements
    //********************************************************************************************************
    void Modeler::GenerateRWPortTE10Element(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new RWPortTE10_1st( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new RWPortTE10_1st( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new RWPortTE10_2nd( nodes, properties ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new RWPortTE10_3th( nodes, properties ) );
		else if ( mElementOrder == 4 ) pElement = Element::Pointer( new RWPortTE10_4th( nodes, properties ) );

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

        if ( mPlasmaMode == true )
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
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
                            if ( abs(cKij) > 0.0 ) A_matrix(i_Kglobal,j_Kglobal) += cKij;
                        }
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

        //Cleaning element matrix and vectors
        EleIdVector.FreeData   ();
        EleResVector.FreeData  ();
        EleStiffMatrix.FreeData();
    }

	//***********************************************************************************************************
    //* - Assembling boundary elements
    //***********************************************************************************************************
    void Modeler::GenerateCoaxPortTEMElement(int* NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

        if      ( mQuadraticGeometry ) pElement = Element::Pointer( new CoaxialPortTEM_2ndQ( nodes, properties ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new CoaxialPortTEM_1st ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new CoaxialPortTEM_2nd ( nodes, properties ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new CoaxialPortTEM_3th ( nodes, properties ) );
        else if ( mElementOrder == 4 ) pElement = Element::Pointer( new CoaxialPortTEM_4th ( nodes, properties ) );

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

        std::cout<<"before assembly"<<std::endl;

        // Assembling to global matrix 
        int eleSize = EleIdVector.size();

        if ( mPlasmaMode == true )
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
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
                            if ( abs(cKij) > 0.0 ) A_matrix(i_Kglobal,j_Kglobal) += cKij;
                        }
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

        //Cleaning element matrix and vectors
        EleIdVector.FreeData   ();
        EleResVector.FreeData  ();
        EleStiffMatrix.FreeData();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - Matrix rotations
	//
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////
    // - Trans(T) * StiffMatrix * T to rotate matrix to the n,t,b coordinate system
    /////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::HO_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix)
    {
		if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;

		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find(*hoit) != mNormals.end() )
            {
                ////////////////  Local coordinate system  /////////////////////
                std::vector<double> n(mNormals[*hoit]);
                std::vector<double> t(3), b(3);
                TangencialCoordinates(n,t,b);

                //////////////  Inv(T) * StiffnessMatrix * T  //////////////////
                int step    = HONodesId.size();
				int step2   = 2*step;
                int dofSize = 3*step;

                std::complex<double> dof_x, dof_y, dof_z;

                // StiffnessMatrix * T
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[i][local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[i][local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

                // Inv(T) * StiffnessMatrix
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

   					StiffMatrix[local_i        ][i] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                    StiffMatrix[local_i + step ][i] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                    StiffMatrix[local_i + step2][i] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
                }

                // Applying {En} = [ R * invEp * Trans(R) ] * {Dn} on plasma sheath BC
                if ( ( mPlasmaMode == true ) && ( mType_Of_BC_Normal[*hoit] == 'D' ) )
                {
                    // {En} = [ R * invEp * Trans(R) ] * {Dn}
                    Matrix< std::complex<double> > RiEpRt;

                    Calculate_En_R_invEp_Rt_Dn( RiEpRt, *hoit );

                    // StiffnessMatrix * [R * invEp * Trans(R)]
                    for (int i=0; i<dofSize; i++)
                    {
                        dof_x = StiffMatrix[i][local_i        ];
                        dof_y = StiffMatrix[i][local_i + step ];
                        dof_z = StiffMatrix[i][local_i + step2];

                        StiffMatrix[i][local_i        ] = RiEpRt[0][0]*dof_x + RiEpRt[0][1]*dof_y + RiEpRt[0][2]*dof_z;
                        StiffMatrix[i][local_i + step ] = RiEpRt[1][0]*dof_x + RiEpRt[1][1]*dof_y + RiEpRt[1][2]*dof_z;
                        StiffMatrix[i][local_i + step2] = RiEpRt[2][0]*dof_x + RiEpRt[2][1]*dof_y + RiEpRt[2][2]*dof_z;
                    }

                    // Trans( [R * invEp * Trans(R)] ) * StiffnessMatrix
                    for (int i=0; i<dofSize; i++)
                    {
                        dof_x = StiffMatrix[local_i        ][i];
                        dof_y = StiffMatrix[local_i + step ][i];
                        dof_z = StiffMatrix[local_i + step2][i];

                    	StiffMatrix[local_i        ][i] = RiEpRt[0][0]*dof_x + RiEpRt[0][1]*dof_y + RiEpRt[0][2]*dof_z;
                        StiffMatrix[local_i + step ][i] = RiEpRt[1][0]*dof_x + RiEpRt[1][1]*dof_y + RiEpRt[1][2]*dof_z;
                        StiffMatrix[local_i + step2][i] = RiEpRt[2][0]*dof_x + RiEpRt[2][1]*dof_y + RiEpRt[2][2]*dof_z;
                    }
                }

            }//end if( mNormals.find(*hoit) != mNormals.end() )

            local_i++;

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    }

	//////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Calculates the matrix that relates En with Dn -> {En} = [ R * invEp * Trans(R) ] * {Dn}
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Calculate_En_R_invEp_Rt_Dn( Matrix< std::complex<double> >& R_invEp_Rt, int NodeId )
    {
        // Permittivity tensor on node ( {Dxyz} = [TEp] * {Exyz} )
        Matrix< std::complex<double> > TEp;

        Node::Pointer pNode = mpModel->GetNode( NodeId - 1 );

        mpColdPlasma->Get_PermittivityTensor_InNode( TEp, pNode, mProblemFrequency );

        // Inverse of TEp
        Matrix< std::complex<double> > invTEp(3,3);

        std::complex<double> detTEp = Determinant( TEp, 3 );

        invTEp[0][0] = ( TEp[1][1] * TEp[2][2] - TEp[1][2] * TEp[2][1] ) / detTEp;
        invTEp[0][1] = ( TEp[0][2] * TEp[2][1] - TEp[0][1] * TEp[2][2] ) / detTEp;
        invTEp[0][2] = ( TEp[0][1] * TEp[1][2] - TEp[0][2] * TEp[1][1] ) / detTEp;

        invTEp[1][0] = ( TEp[1][2] * TEp[2][0] - TEp[1][0] * TEp[2][2] ) / detTEp;
        invTEp[1][1] = ( TEp[0][0] * TEp[2][2] - TEp[0][2] * TEp[2][0] ) / detTEp;
        invTEp[1][2] = ( TEp[0][2] * TEp[1][0] - TEp[0][0] * TEp[1][2] ) / detTEp;

        invTEp[2][0] = ( TEp[1][0] * TEp[2][1] - TEp[1][1] * TEp[2][0] ) / detTEp;
        invTEp[2][1] = ( TEp[0][1] * TEp[2][0] - TEp[0][0] * TEp[2][1] ) / detTEp;
        invTEp[2][2] = ( TEp[0][0] * TEp[1][1] - TEp[0][1] * TEp[1][0] ) / detTEp;

        // Normalizing permittivity tensor with eo 
        invTEp *= 8.8541878176e-12;

        // Surface normal on the node and associated coordinate system
        std::vector<double> n( mNormals[NodeId] );

        std::vector<double> t(3), b(3);

        TangencialCoordinates(n,t,b);

        // Rotation matrix {Entb} = [R] * {Exyz}
        Matrix< std::complex<double> > R(3,3);

        R[0][0] = n[0];  R[0][1] = n[1];  R[0][2] = n[2];
        R[1][0] = t[0];  R[1][1] = t[1];  R[1][2] = t[2];
        R[2][0] = b[0];  R[2][1] = b[1];  R[2][2] = b[2];

        // En = [ R * invEp * Trans(R) ] * Dn
        // Carefull!!, T.Transpose() changes T permanently
        R_invEp_Rt.Resize(3,3);

		R_invEp_Rt = R * invTEp;

		R_invEp_Rt = R_invEp_Rt * R.Transpose();
    }

	//*************************************************************************************************
    //* - Inv(T)*StiffMatrix*T to rotate stiffness matrix to the axisymmetric coordinate system
    //*************************************************************************************************
    void Modeler::Axisym_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix)
    {
        if ( !mAxisymmetric ) return;

		int local_i = 0;

		double distToAxis;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
			if ( ( mNormals.find(*hoit) == mNormals.end() ) && ( DistanceToAxis(*hoit) > mGiDTolerance ) )
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
    //* - Inv(T)*StiffMatrix*T to rotate matrix to the diagonal permittivity tensor base (RLP)
    //* - Being:
    //* - K_RPL = [U][T] * K_CC * [Trans(T)][Herm(U)]
    //* - E_SDP = [T] * E_CC  //  E_RLP = [U] * E_SDP
    //*************************************************************************************************
    void Modeler::Plasma_TKT(std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix)
    {
        // Checking if E parallel tolerance is activated
        if ( mpColdPlasma->Get_Eparallel_Tolerance() == "Off" ) return;

        int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find(*hoit) == mNormals.end() )
            {
                // Complex constants
                std::complex<double> cZero    ( 0.0        , 0.0         );
                std::complex<double> cUnit    ( 0.0        , 1.0         );
                std::complex<double> rUnit    ( 1.0        , 0.0         );
                std::complex<double> rInvSqrt2( 1.0/sqrt(2), 0.0         );
                std::complex<double> cInvSqrt2( 0.0        , 1.0/sqrt(2) );

                // E_RLP = [U] * E_SDP
                Matrix<std::complex<double> > U(3,3);

                U[0][0] = rInvSqrt2; U[0][1] =  cInvSqrt2; U[0][2] = cZero;
                U[1][0] = rInvSqrt2; U[1][1] = -cInvSqrt2; U[1][2] = cZero;
                U[2][0] = cZero    ; U[2][1] =  cZero    ; U[2][2] = rUnit;

                // E_SDP = [T] * E_CC
                Matrix<std::complex<double> > T(3,3);

                Rotation_Matrix_SDP_T_CC( T, *hoit );

                // inv[R] = [U][T]
                Matrix<std::complex<double> > invR(3,3);

                invR = U * T;

                // [R] = [Trans(T)][Herm(U)]
                Matrix<std::complex<double> > R(3,3);

                for (int i=0; i<3; i++)
                {
                    for (int j=0; j<3; j++)
                    {
                        U[i][j] = std::conj(U[i][j]);
                    }
                }

                R = T.Transpose() * U.Transpose();

                //////////////  Inv(R) * StiffMatrix * R  //////////////////
                int step    = HONodesId.size();
                int step2   = 2*step;
                int dofSize = 3*step;

                std::complex<double> dof_x, dof_y, dof_z;

                // StiffMatrix * R
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = R[0][0]*dof_x + R[1][0]*dof_y + R[2][0]*dof_z;
                    StiffMatrix[i][local_i + step ] = R[0][1]*dof_x + R[1][1]*dof_y + R[2][1]*dof_z;
                    StiffMatrix[i][local_i + step2] = R[0][2]*dof_x + R[1][2]*dof_y + R[2][2]*dof_z;
                }

                // Inv(R) * StiffMatrix
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

                    StiffMatrix[local_i        ][i] = invR[0][0]*dof_x + invR[0][1]*dof_y + invR[0][2]*dof_z;
                    StiffMatrix[local_i + step ][i] = invR[1][0]*dof_x + invR[1][1]*dof_y + invR[1][2]*dof_z;
                    StiffMatrix[local_i + step2][i] = invR[2][0]*dof_x + invR[2][1]*dof_y + invR[2][2]*dof_z;
                }

            }//end if( mNormals.find(*hoit) == mNormals.end() && ... )

            local_i++;

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////
    // - Inv(T)*b for PEC and PMC boundary conditions
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::HO_TR(std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector)
    {
        if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;

		int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find( *hoit ) != mNormals.end() )
            {
                //////////////  Local coordinate system  /////////////////
                std::vector<double> n(mNormals[*hoit]);
                std::vector<double> t(3), b(3);
                TangencialCoordinates(n,t,b);

                ///////////////////  Inv(T)*b  ///////////////////////////
                int step  = HONodesId.size();
				int step2 = 2*step;

                std::complex<double> dof_x, dof_y, dof_z;

				dof_x = ResVector[local_i        ];
                dof_y = ResVector[local_i + step ];
                dof_z = ResVector[local_i + step2];

                ResVector[local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
                ResVector[local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
                ResVector[local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;

                // Applying {En} = [ R * invEp * Trans(R) ] * {Dn} on plasma sheath BC
                if ( ( mPlasmaMode == true ) && ( mType_Of_BC_Normal[*hoit] == 'D' ) )
                {
                    Matrix< std::complex<double> > RiEpRt;

                    Calculate_En_R_invEp_Rt_Dn( RiEpRt, *hoit );

                    dof_x = ResVector[local_i        ];
                    dof_y = ResVector[local_i + step ];
                    dof_z = ResVector[local_i + step2];

                    ResVector[local_i        ] = RiEpRt[0][0]*dof_x + RiEpRt[1][0]*dof_y + RiEpRt[2][0]*dof_z;
                    ResVector[local_i + step ] = RiEpRt[0][1]*dof_x + RiEpRt[1][1]*dof_y + RiEpRt[2][1]*dof_z;
                    ResVector[local_i + step2] = RiEpRt[0][2]*dof_x + RiEpRt[1][2]*dof_y + RiEpRt[2][2]*dof_z;
                }

            }// End if( mNormals.find(*hoit) != mNormals.end() )

            local_i++;

        }// End for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    }

	//*************************************************************************************************
    //* - Inv(T)*b for axisymmetric problems.
    //*************************************************************************************************
    void Modeler::Axisym_TR(std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector)
    {
        int local_i = 0;

        double distToAxis;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
			if ( ( mNormals.find(*hoit) == mNormals.end() ) && ( DistanceToAxis(*hoit) > mGiDTolerance ) )
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
    //* - Inv(T)*b for plasma mode
    //*************************************************************************************************
    void Modeler::Plasma_TR(std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector)
    {
        // Checking if E parallel tolerance is activated
        if ( mpColdPlasma->Get_Eparallel_Tolerance() == "Off" ) return;

        int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mNormals.find(*hoit) == mNormals.end() )
            {
                // Complex constants
                std::complex<double> cZero    ( 0.0        , 0.0         );
                std::complex<double> cUnit    ( 0.0        , 1.0         );
                std::complex<double> rUnit    ( 1.0        , 0.0         );
                std::complex<double> rInvSqrt2( 1.0/sqrt(2), 0.0         );
                std::complex<double> cInvSqrt2( 0.0        , 1.0/sqrt(2) );

                // E_RLP = [U] * E_SDP
                Matrix<std::complex<double> > U(3,3);

                U[0][0] = rInvSqrt2; U[0][1] =  cInvSqrt2; U[0][2] = cZero;
                U[1][0] = rInvSqrt2; U[1][1] = -cInvSqrt2; U[1][2] = cZero;
                U[2][0] = cZero    ; U[2][1] =  cZero    ; U[2][2] = rUnit;

                // E_SDP = [T] * E_CC
                Matrix<std::complex<double> > T(3,3);

                Rotation_Matrix_SDP_T_CC( T, *hoit );

                // inv[R] = [U][T]
                Matrix<std::complex<double> > invR(3,3);

                invR = U * T;

                ///////////////////  Inv(T)*b  ////////////////////
                int step  = HONodesId.size();
                int step2 = 2*step;

                std::complex<double> dof_x, dof_y, dof_z;

                dof_x = ResVector[local_i        ];
                dof_y = ResVector[local_i + step ];
                dof_z = ResVector[local_i + step2];

                ResVector[local_i        ] = invR[0][0]*dof_x + invR[0][1]*dof_y + invR[0][2]*dof_z;
                ResVector[local_i + step ] = invR[1][0]*dof_x + invR[1][1]*dof_y + invR[1][2]*dof_z;
                ResVector[local_i + step2] = invR[2][0]*dof_x + invR[2][1]*dof_y + invR[2][2]*dof_z;
            }

            local_i++;
        }
    }

	//*************************************************************************************************
    //* - Local system of coordinates in PEC surface for n
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
    //* -  Calculates the rotation matrix from CC to SDP (E_SDP = [T] * E_CC)
    //*************************************************************************************************
    void Modeler::Rotation_Matrix_SDP_T_CC( Matrix<std::complex<double> >& T, unsigned int NodeId )
    {
        // Unit vectors
        Vector<double> n(3), t(3), b(3);

        // Get node
        Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

        // Plasma parameters on node
        Vector<double> Bext = mpColdPlasma->Get_Bexternal_InNode( pNode );

        double Bx = Bext[0];
        double By = Bext[1];
        double Bz = Bext[2];

        // b = B / |B|
        double Bnorm = std::sqrt( Bx*Bx + By*By + Bz*Bz );

        b[0] = Bx / Bnorm;
        b[1] = By / Bnorm;
        b[2] = Bz / Bnorm;

        // n _|_ b
        double nnorm = std::sqrt( b[0]*b[0] + b[1]*b[1] );

        if (nnorm > 0.00)
        {
            n[0] = b[1] / nnorm;
            n[1] =-b[0] / nnorm;
            n[2] = 0.00;
        }
        else
        {
            n[0] = b[2];
            n[1] = 0.00;
            n[2] = 0.00;
        }

        // t = b x n
        t[0] = b[1]*n[2] - b[2]*n[1];
        t[1] = b[2]*n[0] - b[0]*n[2];
        t[2] = b[0]*n[1] - b[1]*n[0];

        // Rotation matrix from CC to SDP (E_SDP = [T] * E_CC)
        T[0][0] = n[0]; T[0][1] = n[1]; T[0][2] = n[2];
        T[1][0] = t[0]; T[1][1] = t[1]; T[1][2] = t[2];
        T[2][0] = b[0]; T[2][1] = b[1]; T[2][2] = b[2];
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
    //* - Axisymmetric unitary vectors
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
    //* - Trans(T)*StiffMatrix*T, PBC nodes
    //*************************************************************************************************
    void Modeler::PBC_FrontBack_TKT(std::vector<int>& HONodesId,
		                            Vector<int>& IdVector,
							        Matrix<std::complex<double> >& StiffMatrix)
   {
	    if ( ( mPBC_NodeNodePairs_Front.size   () +
               mPBC_NodeEdgePairs_Front.size   () +
               mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

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
    //* - Trans(T)*StiffMatrix*T, PBC node to node
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
    //* - Trans(T)*StiffMatrix*T, PBC node to edge
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToEdge_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<int> bNodes( mPBC_NodeEdgePairs_Front[NodeId] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element edge
		double natCoord = mPBC_NodeEdgeCoord_Front[NodeId];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

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
    //* - Inv(T)*StiffMatrix*T, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TKT(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix)
	{
		std::vector<int> bNodes( mPBC_NodeElementPairs_Front[NodeId] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Front[NodeId] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_StiffMatrixReDef(local_i, oldSize, newSize, N, StiffMatrix);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//*************************************************************************************************
    //* - Values of the Back element base N for the proyected location of the Front node
    //*************************************************************************************************
    void Modeler::PBC_NaturalBase_N_Surface( std::vector<double>& N, std::vector<double>& natCoord )
    {
        double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;

        if ( mElementOrder == 0 )
        {
		    N.resize(3);

            N[0] = L1;
            N[1] = L2;
            N[2] = L3;
			N[3] = 27.0 * L1 * L2 * L3;
        }
        else if ( mElementOrder == 1 )
        {
		    N.resize(3);

            N[0] = L1;
            N[1] = L2;
            N[2] = L3;
        }
        else if ( mElementOrder == 2 )
        {
		    N.resize(6);

            N[0] = L1 * ( 2.00*L1 - 1.00 );
            N[1] = L2 * ( 2.00*L2 - 1.00 );
            N[2] = L3 * ( 2.00*L3 - 1.00 );
		    N[3] = 4.00 * L1 * L2;
		    N[4] = 4.00 * L2 * L3;
            N[5] = 4.00 * L3 * L1;
        }
        else if ( mElementOrder == 3 )
        {
            N.resize(10);

            N[0] = 0.50 * L1 * ( 3.00*L1 - 1.00 ) * ( 3.00*L1 - 2.00 );
            N[1] = 0.50 * L2 * ( 3.00*L2 - 1.00 ) * ( 3.00*L2 - 2.00 );
            N[2] = 0.50 * L3 * ( 3.00*L3 - 1.00 ) * ( 3.00*L3 - 2.00 );
            N[3] = 4.50 * L1 * L2 * ( 3.00*L1 - 1.00 );
            N[4] = 4.50 * L1 * L2 * ( 3.00*L2 - 1.00 );
            N[5] = 4.50 * L2 * L3 * ( 3.00*L2 - 1.00 );
            N[6] = 4.50 * L2 * L3 * ( 3.00*L3 - 1.00 );
            N[7] = 4.50 * L3 * L1 * ( 3.00*L3 - 1.00 );
            N[8] = 4.50 * L3 * L1 * ( 3.00*L1 - 1.00 );
            N[9] = 27.0 * L1 * L2 * L3;
        }
        else if ( mElementOrder == 4 )
        {
            N.resize(15);

            N[0 ] = (1.00/6.00)  * L1 * ( 4.00*L1 - 1.00 ) * ( 4.00*L1 - 2.00 ) * ( 4.00*L1 - 3.00 );
            N[1 ] = (1.00/6.00)  * L2 * ( 4.00*L2 - 1.00 ) * ( 4.00*L2 - 2.00 ) * ( 4.00*L2 - 3.00 );
            N[2 ] = (1.00/6.00)  * L3 * ( 4.00*L3 - 1.00 ) * ( 4.00*L3 - 2.00 ) * ( 4.00*L3 - 3.00 );
            N[3 ] = (8.00/3.00)  * L1 * L2 * ( 4.00*L1 - 1.00 ) * ( 4.00*L1 - 2.00 );
            N[4 ] = (4.00)       * L1 * L2 * ( 4.00*L1 - 1.00 ) * ( 4.00*L2 - 1.00 );
            N[5 ] = (8.00/3.00)  * L1 * L2 * ( 4.00*L2 - 1.00 ) * ( 4.00*L2 - 2.00 );
            N[6 ] = (8.00/3.00)  * L2 * L3 * ( 4.00*L2 - 1.00 ) * ( 4.00*L2 - 2.00 );
            N[7 ] = (4.00)       * L2 * L3 * ( 4.00*L2 - 1.00 ) * ( 4.00*L3 - 1.00 );
            N[8 ] = (8.00/3.00)  * L2 * L3 * ( 4.00*L3 - 1.00 ) * ( 4.00*L3 - 2.00 );
            N[9 ] = (8.00/3.00)  * L1 * L3 * ( 4.00*L3 - 1.00 ) * ( 4.00*L3 - 2.00 );
            N[10] = (4.00)       * L1 * L3 * ( 4.00*L1 - 1.00 ) * ( 4.00*L3 - 1.00 );
            N[11] = (8.00/3.00)  * L1 * L3 * ( 4.00*L1 - 1.00 ) * ( 4.00*L1 - 2.00 );
            N[12] = (32.00) * L1 * L2 * L3 * ( 4.00*L1 - 1.00 );
            N[13] = (32.00) * L1 * L2 * L3 * ( 4.00*L2 - 1.00 );
            N[14] = (32.00) * L1 * L2 * L3 * ( 4.00*L3 - 1.00 );
        }
    }

	//*************************************************************************************************
    //* - Values of the Back element edge base N for the proyected location of the Front node
    //*************************************************************************************************
    void Modeler::PBC_NaturalBase_N_Line( std::vector<double>& N, double natCoord )
    {
        double L2 = natCoord;
        double L1 = 1.00 - L2;

        if      ( mElementOrder <= 1 )
        {
		    N.resize(2);

            N[0] = L1;
            N[1] = L2;
        }
        else if ( mElementOrder == 2 )
        {
		    N.resize(3);

		    N[0] = L1 * ( 2.00*L1 - 1.00 );
		    N[1] = L2 * ( 2.00*L2 - 1.00 );
		    N[2] = 4.00 * L1 * L2;
        }
        else if ( mElementOrder == 3 )
        {
            N.resize(4);

            N[0] = 0.50 * L1 * ( 3.00*L1 - 1.00 ) * ( 3.00*L1 - 2.00 );
            N[1] = 0.50 * L2 * ( 3.00*L2 - 1.00 ) * ( 3.00*L2 - 2.00 );
            N[2] = 4.50 * L1 * L2 * ( 3.00*L1 - 1.00 );
            N[3] = 4.50 * L1 * L2 * ( 3.00*L2 - 1.00 );
        }
        else if ( mElementOrder == 4 )
        {
            N.resize(5);

            N[0] = (1.00/6.00) * L1 * ( 4.00*L1 - 1.00 ) * ( 4.00*L1 - 2.00 ) * ( 4.00*L1 - 3.00 );
            N[1] = (1.00/6.00) * L2 * ( 4.00*L2 - 1.00 ) * ( 4.00*L2 - 2.00 ) * ( 4.00*L2 - 3.00 );
            N[2] = (8.00/3.00) * L1 * L2 * ( 4.00*L1 - 1.00 ) * ( 4.00*L1 - 2.00 );
            N[3] = (4.00)      * L1 * L2 * ( 4.00*L1 - 1.00 ) * ( 4.00*L2 - 1.00 );
            N[4] = (8.00/3.00) * L1 * L2 * ( 4.00*L2 - 1.00 ) * ( 4.00*L2 - 2.00 );
        }
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
        if ( ( mPBC_NodeNodePairs_Front.size   () +
               mPBC_NodeEdgePairs_Front.size   () +
               mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

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
		std::vector<int> bNodes( mPBC_NodeEdgePairs_Front[NodeId] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element edge
		double natCoord = mPBC_NodeEdgeCoord_Front[NodeId];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_ResVectorReDef(local_i, oldSize, newSize, N, ResVector);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//*************************************************************************************************
    //* - Trans(T)*b, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<int> bNodes( mPBC_NodeElementPairs_Front[NodeId] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) pbNodes.push_back( mpModel->GetNode( bNodes[in] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, pbNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Front[NodeId] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

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
		std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Right natural base values of the proyected Left node on the Right element edge
		double natCoord = mPBC_NodeEdgeCoord_Right[NodeId];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_StiffMatrixReDef_Rot(local_i, oldSize, newSize, N, NodeId, StiffMatrix);

		local_i = local_i + ( lNodes.size() - 1 );
	}

    //*************************************************************************************************
    //* - Trans(T)*StiffnesMatrix*T, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TKT( int NodeId, int& local_i   ,
                                            std::vector<int>& HONodesId,
                                            Vector<int>& IdVector      ,
                                            Matrix<std::complex<double> >& StiffMatrix )
	{
		std::vector<int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Right natural base values of the proyected Left node on the Right element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_StiffMatrixReDef_Rot(local_i, oldSize, newSize, N, NodeId, StiffMatrix);

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//*************************************************************************************************
    //* - Redefinition of the stiffness matrix for PBC boundary conditions
    //*************************************************************************************************
	void Modeler::PBC_StiffMatrixReDef_Rot(int local_i, int oldSize, int newSize,
		                                   std::vector<double>& N, int NodeId   ,
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
		std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Left natural base values of the proyected Right node on the Left element edge
		double natCoord = mPBC_NodeEdgeCoord_Right[NodeId];

        std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_ResVectorReDef_Rot(local_i, oldSize, newSize, N, NodeId, ResVector);

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//*************************************************************************************************
    //* - Trans(T)*b, PBC node to element
    //*************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TR(int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector)
	{
		std::vector<int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) plNodes.push_back( mpModel->GetNode( lNodes[i] - 1 ) );

		int oldSize = IdVector.size();

		PBC_IdVectorReDef(local_i, plNodes, HONodesId, IdVector);

		int newSize = IdVector.size();

        // Left natural base values of the proyected Right node on the Left element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

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

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Trans(T)*StiffMatrix*T to eliminate one of the nodes in a contact pair
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::ContactHO_TKT( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
    {
  //      if ( mContactPairs.size() == 0 ) return;

		//int local_i = 0;

  //      std::vector<int>::iterator hoit;

  //      for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
  //      {
  //          if ( mType_Of_BC_Normal[*hoit] == 'R' )
  //          {
  //              int i;

  //              int step    = HONodesId.size();
		//		int step2   = 2*step;
  //              int dofSize = 3*step;

		//		ContactPairData& pContactPair = mContactPairs[*hoit];

		//		Node::Pointer pNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

		//		HONodesId[local_i        ] = pNode->Id();

		//		IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
		//		IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
		//		IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

  //              ////////////////  Local coordinate system  /////////////////////
  //              std::vector<double> n = pContactPair.ctcNormal; //std::vector<double> n (mRNormal[*hoit]);
  //              std::vector<double> t(3), b(3);
  //              TangencialCoordinates(n,t,b);

  //              /////////////////////////////////////////////////////////////////////////////

		//		// Material properties matrix
		//		Properties::Pointer Properties;

		//		double eo = 8.8541878176e-12;

		//		double freq = mProblemFrequency;

  //              double sigma, eps_real, eps_imag;

		//		// Material R
		//		Properties = mpModel->GetProperties(pContactPair.matR);

  //              sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
		//		eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
		//	    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		//        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

  //              //cEpsR /= eo;

  //              std::complex<double> invEp = (1.00 / cEpsR);

  //              //////////////////////////////////////////////////////////////////////////////

  //              std::complex<double> dof_x, dof_y, dof_z;

  //              // StiffnessMatrix * T
  //              for (int i=0; i<dofSize; i++)
  //              {
  //                  dof_x = StiffMatrix[i][local_i        ];
  //                  dof_y = StiffMatrix[i][local_i + step ];
  //                  dof_z = StiffMatrix[i][local_i + step2];

  //                  StiffMatrix[i][local_i        ] = invEp*n[0]*dof_x + invEp*n[1]*dof_y + invEp*n[2]*dof_z;
  //                  StiffMatrix[i][local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
  //                  StiffMatrix[i][local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
  //              }

  //              // Inv(T) * StiffnessMatrix
  //              for (int i=0; i<dofSize; i++)
  //              {
  //                  dof_x = StiffMatrix[local_i        ][i];
  //                  dof_y = StiffMatrix[local_i + step ][i];
  //                  dof_z = StiffMatrix[local_i + step2][i];

  // 					StiffMatrix[local_i        ][i] = invEp*n[0]*dof_x + invEp*n[1]*dof_y + invEp*n[2]*dof_z;
  //                  StiffMatrix[local_i + step ][i] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
  //                  StiffMatrix[local_i + step2][i] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
  //              }

  //          }//end if( mContactPairs.find(*hoit) != mContactPairs.end() )

  //          local_i++;

  //      }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)

  //      ///////////////////////////////////////////////////////////////////////////////////////////////////////////

		//local_i = 0;

  //      for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
  //      {
  //          if ( mType_Of_BC_Normal[*hoit] == 'L' )
  //          {
  //              int i;

  //              int step    = HONodesId.size();
		//		int step2   = 2*step;
  //              int dofSize = 3*step;

  //              ////////////////  Local coordinate system  /////////////////////
  //              std::vector<double> n(mLNormal[*hoit]);
  //              std::vector<double> t(3), b(3);
  //              TangencialCoordinates(n,t,b);

  //              ///////////////////////////////////////////////////////////////

  //              double eo = 8.8541878176e-12;

  //              std::complex<double> cEpsL( eo, 0.0 );

  //              //cEpsL /= eo;

  //              std::complex<double> invEp = (1.00 / cEpsL);

  //              ///////////////////////////////////////////////////////////////


  //              std::complex<double> dof_x, dof_y, dof_z;

  //              // StiffnessMatrix * T
  //              for (int i=0; i<dofSize; i++)
  //              {
  //                  dof_x = StiffMatrix[i][local_i        ];
  //                  dof_y = StiffMatrix[i][local_i + step ];
  //                  dof_z = StiffMatrix[i][local_i + step2];

  //                  StiffMatrix[i][local_i        ] = invEp*n[0]*dof_x + invEp*n[1]*dof_y + invEp*n[2]*dof_z;
  //                  StiffMatrix[i][local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
  //                  StiffMatrix[i][local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
  //              }

  //              // Inv(T) * StiffnessMatrix
  //              for (int i=0; i<dofSize; i++)
  //              {
  //                  dof_x = StiffMatrix[local_i        ][i];
  //                  dof_y = StiffMatrix[local_i + step ][i];
  //                  dof_z = StiffMatrix[local_i + step2][i];

  // 					StiffMatrix[local_i        ][i] = invEp*n[0]*dof_x + invEp*n[1]*dof_y + invEp*n[2]*dof_z;
  //                  StiffMatrix[local_i + step ][i] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
  //                  StiffMatrix[local_i + step2][i] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
  //              }

  //          }

  //          local_i++;

  //      }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
  //      
  //      return;

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

				// Material properties matrix
				Properties::Pointer Properties;

				double eo = 8.8541878176e-12;

				double freq = mProblemFrequency;

                double sigma, eps_real, eps_imag;

				// Material R
				Properties = mpModel->GetProperties(pContactPair.matR);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );

				// Material L
				Properties = mpModel->GetProperties(pContactPair.matL);

                sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		        std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );
            
				///////////////////////////////////////////////////    Material properties matrix    ////////////////////////////////////////////////////////////
				std::complex<double> ContactMatrix[3][3];

                std::complex<double> rmat = (cEpsL/cEpsR) - 1.00;

				Vector<double> n = pContactPair.ctcNormal;

				ContactMatrix[0][0] = (n[0]*n[0]) * rmat + 1.00; ContactMatrix[0][1] = (n[0]*n[1]) * rmat       ; ContactMatrix[0][2] = (n[0]*n[2]) * rmat;
				ContactMatrix[1][0] = (n[1]*n[0]) * rmat       ; ContactMatrix[1][1] = (n[1]*n[1]) * rmat + 1.00; ContactMatrix[1][2] = (n[1]*n[2]) * rmat;
				ContactMatrix[2][0] = (n[2]*n[0]) * rmat       ; ContactMatrix[2][1] = (n[2]*n[1]) * rmat       ; ContactMatrix[2][2] = (n[2]*n[2]) * rmat + 1.00;

                //////////////////////////////////////////////////     Inv(T)*StiffnesMatrix*T     //////////////////////////////////////////////////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                // StiffnesMatrix*T
                for (i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[i][local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[i][local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

                // Inv(T)*StiffnesMatrix
                for (i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

   					StiffMatrix[local_i        ][i] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[local_i + step ][i] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[local_i + step2][i] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

            }//end if( mContactPairs.find(*hoit) != mContactPairs.end() )

            local_i++;

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        
    }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Inv(T)*b in contact nodes
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::ContactHO_TR( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
    {
  //      if ( mContactPairs.size() == 0 ) return;

		//int local_i = 0;

  //      std::vector<int>::iterator hoit;

  //      for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
  //      {
  //          if ( mContactPairs.find(*hoit) != mContactPairs.end() )
  //          {
  //              int step    = HONodesId.size();
		//		int step2   = 2*step;

		//		ContactPairData& pContactPair = mContactPairs[*hoit];

		//		Node::Pointer pNode = mpModel->GetNode( (pContactPair.ctcNode) - 1 );

		//		HONodesId[local_i        ] = pNode->Id();

		//		//IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
		//		IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
		//		IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

  //              //////////////  Local coordinate system  /////////////////
  //              std::vector<double> n(pContactPair.ctcNormal);
  //              std::vector<double> t(3), b(3);
  //              TangencialCoordinates(n,t,b);

  //              ///////////////////  Inv(T)*b  ///////////////////////////
  //              std::complex<double> dof_x, dof_y, dof_z;

		//		dof_x = ResVector[local_i        ];
  //              dof_y = ResVector[local_i + step ];
  //              dof_z = ResVector[local_i + step2];

  //              ResVector[local_i        ] = n[0]*dof_x + n[1]*dof_y + n[2]*dof_z;
  //              ResVector[local_i + step ] = t[0]*dof_x + t[1]*dof_y + t[2]*dof_z;
  //              ResVector[local_i + step2] = b[0]*dof_x + b[1]*dof_y + b[2]*dof_z;
  //          }

  //          local_i++;
  //      }

  //      return;
 //////////////////////////////////////////////////////////////////////////////////////////////       
        
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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Trans(T)*StiffMatrix*T to eliminate one of the nodes in a contact pair for Dielectric-Plasma interfaces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::ContactHO_TKT_Plasma( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
    {
    	if ( mContactPairs.size() == 0 ) return;

    	int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mContactPairs.find(*hoit) != mContactPairs.end() )
            {
                int RNodeId = (*hoit);
                int LNodeId = mContactPairs[RNodeId].ctcNode;

                ///////////////////////////////////////////////////   Modify Id vector    ////////////////////////////////////////////////
                int step    = HONodesId.size();
    			int step2   = 2*step;
                int dofSize = 3*step;

                Node::Pointer pNode = mpModel->GetNode( LNodeId - 1 );

                // Carefull!!!, this changes the content of (*hoit)
                HONodesId[local_i        ] = pNode->Id();

    			IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
    			IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
    			IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

                //////////////////////////////////////////////    Inv(T)*StiffnesMatrix*T    /////////////////////////////////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                Matrix< std::complex<double> > ContactMatrix(3,3);

                Get_Contact_Matrix_Plasma( ContactMatrix, RNodeId );

                // StiffnesMatrix*T
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[i][local_i + step ] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[i][local_i + step2] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

                // Inv(T)*StiffnesMatrix
                for (int i=0; i<dofSize; i++)
                {
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

    				StiffMatrix[local_i        ][i] = ContactMatrix[0][0]*dof_x + ContactMatrix[1][0]*dof_y + ContactMatrix[2][0]*dof_z;
                    StiffMatrix[local_i + step ][i] = ContactMatrix[0][1]*dof_x + ContactMatrix[1][1]*dof_y + ContactMatrix[2][1]*dof_z;
                    StiffMatrix[local_i + step2][i] = ContactMatrix[0][2]*dof_x + ContactMatrix[1][2]*dof_y + ContactMatrix[2][2]*dof_z;
                }

            }//end if( mContactPairs.find(*hoit) != mContactPairs.end() )

            local_i++;

        }//end for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Calculates the contact matrix for Dielectric-Plasma interfaces
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Get_Contact_Matrix_Plasma( Matrix< std::complex<double> >& ContactMatrix, int RNodesId )
    {
        ContactPairData& pContactPair = mContactPairs[RNodesId];

        Node::Pointer prNode = mpModel->GetNode( RNodesId             - 1 );
        Node::Pointer plNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

        std::vector<double> n = pContactPair.ctcNormal;
        std::vector<double> t(3), b(3);

        TangencialCoordinates(n,t,b);

        // Material properties constants
        Properties::Pointer Properties;

        double eo   = 8.8541878176e-12;
        double freq = mProblemFrequency;

        std::complex<double> cZero( 0.0, 0.0 );
        std::complex<double> cUnit( 0.0, 1.0 );

        double sigma;
        double epr_real, epr_imag;
        double mur_real, mur_imag;

        Matrix< std::complex<double> > TEn(3,3);

        // Material R
        Properties = mpModel->GetProperties( pContactPair.matR );

        sigma      = (*Properties)( IHL_ELECTRIC_CONDUCTIVITY  );
        epr_real   = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
        epr_imag   = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );
        mur_real   = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
        mur_imag   = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );

        // Plasma permittivity tensor
        if ( (sigma == 0.0) && (epr_real == 1.0) && (epr_imag == 0.0) && (mur_real == 1.0) && (mur_imag == 0.0) )
        {
            mpColdPlasma->Get_PermittivityTensor_InNode( TEn, prNode, freq );
        }
        // IHL permittivity tensor
        else
        {
            std::complex<double> cEps( eo*epr_real , eo*epr_imag + ( sigma / freq ) );

            TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero;
            TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero;
            TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
        }
        // [TR] * {Er} = [TL] * {El}
        Matrix< std::complex<double> > TR(3,3);

        TR[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0];
        TR[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1];
        TR[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2];

        TR[1][0] = t[0]; TR[1][1] = t[1]; TR[1][2] = t[2];
        TR[2][0] = b[0]; TR[2][1] = b[1]; TR[2][2] = b[2];

        // Material L
        Properties = mpModel->GetProperties( pContactPair.matL );

        sigma      = (*Properties)( IHL_ELECTRIC_CONDUCTIVITY  );
        epr_real   = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
        epr_imag   = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );
        mur_real   = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
        mur_imag   = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );

        // Plasma permittivity tensor
        if ( (sigma == 0.0) && (epr_real == 1.0) && (epr_imag == 0.0) && (mur_real == 1.0) && (mur_imag == 0.0) )
        {
            mpColdPlasma->Get_PermittivityTensor_InNode( TEn, plNode, freq );
        }
        // IHL permittivity tensor
        else
        {
            std::complex<double> cEps( eo*epr_real , eo*epr_imag + ( sigma / freq ) );

            TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero;
            TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero;
            TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
        }
        // [TR] * {Er} = [TL] * {El}
        Matrix< std::complex<double> > TL(3,3);

        TL[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0];
        TL[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1];
        TL[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2];

        TL[1][0] = t[0]; TL[1][1] = t[1]; TL[1][2] = t[2];
        TL[2][0] = b[0]; TL[2][1] = b[1]; TL[2][2] = b[2];

        // {Er} = ( [invTR] * [TL] ) * {El}
        Matrix< std::complex<double> > invTR ( 3, 3 );

        std::complex<double> detTR = Determinant( TR, 3 );

        invTR[0][0] = ( TR[1][1] * TR[2][2] - TR[1][2] * TR[2][1] ) / detTR;
        invTR[0][1] = ( TR[0][2] * TR[2][1] - TR[0][1] * TR[2][2] ) / detTR;
        invTR[0][2] = ( TR[0][1] * TR[1][2] - TR[0][2] * TR[1][1] ) / detTR;

        invTR[1][0] = ( TR[1][2] * TR[2][0] - TR[1][0] * TR[2][2] ) / detTR;
        invTR[1][1] = ( TR[0][0] * TR[2][2] - TR[0][2] * TR[2][0] ) / detTR;
        invTR[1][2] = ( TR[0][2] * TR[1][0] - TR[0][0] * TR[1][2] ) / detTR;

        invTR[2][0] = ( TR[1][0] * TR[2][1] - TR[1][1] * TR[2][0] ) / detTR;
        invTR[2][1] = ( TR[0][1] * TR[2][0] - TR[0][0] * TR[2][1] ) / detTR;
        invTR[2][2] = ( TR[0][0] * TR[1][1] - TR[0][1] * TR[1][0] ) / detTR;

        // Contact matrix = [invTR] * [TL]
        ContactMatrix.Resize(3,3);

        ContactMatrix = invTR * TL;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Inv(T)*b in contact nodes for Dielectric-Plasma interfaces
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::ContactHO_TR_Plasma( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
    {
        if ( mContactPairs.size() == 0 ) return;

    	int local_i = 0;

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            if ( mContactPairs.find(*hoit) != mContactPairs.end() )
            {
                int RNodeId = (*hoit);
                int LNodeId = mContactPairs[RNodeId].ctcNode;

                ///////////////////////////////////////////////   Modify Id vector    ///////////////////////////////////////////
                int step  = HONodesId.size();
    			int step2 = 2*step;

                Node::Pointer pNode = mpModel->GetNode( LNodeId - 1 );

                // Carefull!!!, this changes the content of (*hoit)
    			HONodesId[local_i        ] = pNode->Id();

    			IdVector [local_i        ] = pNode->pDofcEx()->EquationId();
    			IdVector [local_i + step ] = pNode->pDofcEy()->EquationId();
    			IdVector [local_i + step2] = pNode->pDofcEz()->EquationId();

                ///////////////////////////////////////////////////  Inv(T)*b  //////////////////////////////////////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                Matrix< std::complex<double> > ContactMatrix;

                Get_Contact_Matrix_Plasma( ContactMatrix, RNodeId );

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

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Returns the determinant of a square matrix M of size n (for complex matrices)
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::complex<double> Modeler::Determinant( Matrix< std::complex<double> >& M, int n )
    {
        std::complex<double> det( 0.0, 0.0 );

        Matrix< std::complex<double> > subM(n,n);

        if (n == 1)
        {
            return M[0][0];
        }

        if (n == 2)
        {
            return ( (M[0][0] * M[1][1]) - (M[1][0] * M[0][1]) );
        }

        for ( int x = 0; x < n; x++ )
        {
            int subi = 0;
            for ( int i = 1; i < n; i++ )
            {
                int subj = 0;
                for ( int j = 0; j < n; j++ )
                {
                    if (j == x) continue;

                    subM[subi][subj] = M[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + ( pow( -1, x ) * M[0][x] * Determinant( subM, n-1 ) );
        }

        return det;
    }

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// - PEC, PMC, PBC and projection elements
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns each PBC node to its corresponding std::map
    ////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::GeneratePBCElement(int* NodesId, unsigned int PropertiesId)
    {
        // Each material is assigned to a different group to speed the PBC set-up
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

		int FaceTypeId = (*properties)(COMPLEX_IBC)[0];
        int MaterialId = NodesId[3];

		mPBCTolerance = (*properties)(COMPLEX_IBC)[1];

		std::vector<int> HONodesId;
		PushHONodesOnSurface(NodesId, HONodesId);

        //Front face
		if      (FaceTypeId == 11)
		{
		    mPBCFrontElements[MaterialId].push_back(HONodesId);
		}
		//Back face
		else if (FaceTypeId == 12)
		{
			mPBCBackElements[MaterialId].push_back(HONodesId);
		}
		//Right face
		else if (FaceTypeId == 21)
		{
		    mPBCRightElements[MaterialId].push_back(HONodesId);
		}
		//Left face
		else if (FaceTypeId == 22)
		{
			mPBCLeftElements[MaterialId].push_back(HONodesId);
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
		if ( mPBCFrontElements.size() > 0 ) Set_FrontBack_PBC();
        if ( mPBCRightElements.size() > 0 ) Set_RightLeft_PBC();
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns the nodes in the front periodic surface to the elements in the back periodic surface.
	// - Front-back surfaces must be placed in the XY-plane.
	// - Only for 2nd order elements.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::Set_FrontBack_PBC()
    {
		// Loop through materials in Back PBC
        std::map<unsigned int, std::vector<std::vector<int> > >::iterator it_Back;
        for( it_Back = mPBCBackElements.begin(); it_Back != mPBCBackElements.end(); ++it_Back )
        {
            // Material ID
			int pbcMaterial = it_Back->first;

            // Group of Back elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetBack( it_Back->second );

            // Group of Frontal elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetFront( mPBCFrontElements[pbcMaterial] );

            // List of nodes (sorted and unique) in the Front PBC
            std::list<int> FrontNodesList;

             // Filling the list with the nodes on each Front element
            std::vector<std::vector<int> >::iterator it_eF;
            for( it_eF = elemSetFront.begin(); it_eF != elemSetFront.end(); ++it_eF )
            {
                std::vector<int>::iterator it_nF;
                for( it_nF = (*it_eF).begin(); it_nF != (*it_eF).end(); ++it_nF )
				{
                    FrontNodesList.push_back( *it_nF );
                }
            }

            // Sorting list
            FrontNodesList.sort();

            // Collapsing equal nodes Id
            FrontNodesList.unique();

            // Loop through the Back elements to calculate the Frontal nodes projections
			std::vector<std::vector<int> >::iterator it_eB;
			for( it_eB = elemSetBack.begin(); it_eB != elemSetBack.end() && FrontNodesList.size() > 0 ; ++it_eB )
            {
                // Recovering Back element nodes coordinates on the XY plane
                std::vector<double> eX, eY;

                for ( int i=0; i<(*it_eB).size(); i++ )
                {
                    eX.push_back( mpModel->GetNode( (*it_eB)[i]-1 )->X() );
                    eY.push_back( mpModel->GetNode( (*it_eB)[i]-1 )->Y() );
                }

                // Calculating cross section (only nodes inside this area will be consider)
                double eCross_Section;
                double edge01_Lenght = sqrt( pow(eX[0]-eX[1], 2) + pow(eY[0]-eY[1], 2) );
                double edge02_Lenght = sqrt( pow(eX[0]-eX[2], 2) + pow(eY[0]-eY[2], 2) );

                if ( edge02_Lenght > edge01_Lenght ) eCross_Section = 1.1 * edge02_Lenght;
                else                                 eCross_Section = 1.1 * edge01_Lenght;

                // Parameters of the geometrical mapping interpolation
                std::vector<double> a, b;

                // Computing the parameters a,b of the geometrical mapping interpolation
                PBC_Mapping_Parameters( a, b, eX, eY );

                // Front nodes projected on this Back element
                std::vector<int> AssignedFrontNodes;

                // Loop over the list of nodes in the Front PBC
                std::list<int>::iterator it_nFl;
                for( it_nFl = FrontNodesList.begin(); it_nFl != FrontNodesList.end(); ++it_nFl )
                {
                    // Recovering the node pointer from its Id
                    Node::Pointer pfNode = mpModel->GetNode( (*it_nFl) - 1 );

                    // Cartesina coordinates of the node
                    double carX = pfNode->X();
                    double carY = pfNode->Y();

                    // Checking if the Front node projection is close enough to the Back element
                    if ( sqrt( pow(eX[0]-carX, 2) + pow(eY[0]-carY, 2) ) > eCross_Section ) continue;

                    // Natural coordinates of the Frontal node projection on the Back element
                    double natX;
                    double natY;

                    if ( a.size() == 6 )
                    {
                        natX = a[0] + ( a[1] * carX ) + ( a[2] * carX * carX ) + ( a[3] * carX * carY ) + ( a[4] * carY ) + ( a[5] * carY * carY );
                        natY = b[0] + ( b[1] * carX ) + ( b[2] * carX * carX ) + ( b[3] * carX * carY ) + ( b[4] * carY ) + ( b[5] * carY * carY );
                    }
                    else
                    {
                        natX = a[0] + ( a[1] * carX ) + ( a[2] * carY );
                        natY = b[0] + ( b[1] * carX ) + ( b[2] * carY );
                    }

                    // If the Front node is outside the Back element -> continue
                    if ( ( natX < (-mPBCTolerance) ) || ( natY < (-mPBCTolerance) ) || ( (natX + natY) > (1.00 + sqrt(2)*mPBCTolerance) ) ) continue;

                    // This node is going to be assigned
                    AssignedFrontNodes.push_back( *it_nFl );

                    // Setting DOFs for this node
                    std::complex<double> cZero(0.00,0.00);

                    (*mpModel)(cEx, *pfNode) = cZero; pfNode->pDofcEx()->FixDof();
                    (*mpModel)(cEy, *pfNode) = cZero; pfNode->pDofcEy()->FixDof();
                    (*mpModel)(cEz, *pfNode) = cZero; pfNode->pDofcEz()->FixDof();

                    // Cheking if Frontal node is projected onto a Back node
                    int BackNodeId = PBC_Node_To_Node( natX, natY, *it_eB );

                    if ( BackNodeId > 0 )
                    {
                        mPBC_NodeNodePairs_Front[pfNode->Id()] = BackNodeId;
                        continue;
                    }

                    // Cheking if Frontal node is projected onto a Back element edge
                    std::vector<int> BackEdgeNodeIds;

                    double NatCoord_InEdge = PBC_Node_To_Edge( natX, natY, *it_eB, BackEdgeNodeIds );

                    if ( BackEdgeNodeIds.size() > 0 )
                    {
                        mPBC_NodeEdgePairs_Front[pfNode->Id()] = BackEdgeNodeIds;
                        mPBC_NodeEdgeCoord_Front[pfNode->Id()] = NatCoord_InEdge;
                        continue;
                    }

                    // Frontal node is projected onto a Back element
                    std::vector<double> NatCoord_InElement(2);
                    NatCoord_InElement[0] = natX;
                    NatCoord_InElement[1] = natY;

                    mPBC_NodeElementPairs_Front[pfNode->Id()] = (*it_eB);
                    mPBC_NodeElementCoord_Front[pfNode->Id()] = NatCoord_InElement;

                }// End loop over the list of nodes in the Front PBC

                // Removing from the FrontNodesList the assigned Front nodes
                std::vector<int>::iterator ita;
                for(ita = AssignedFrontNodes.begin(); ita != AssignedFrontNodes.end(); ++ita) FrontNodesList.remove(*ita);

            }// End loop through Back elements

            // Checking if all the Front nodes in a pcbMaterial have been assigned
            if ( FrontNodesList.size() > 0 )
            {
                std::cout << FrontNodesList.size() << " Front nodes are not assigned in PBC material " << pbcMaterial << std::endl;
            }

        }// End loop through materials in Back PBC

	    mPBCFrontElements.clear();
	    std::map<unsigned int, std::vector<std::vector<int> > >().swap(mPBCFrontElements);

		mPBCBackElements.clear();
		std::map<unsigned int, std::vector<std::vector<int> > >().swap(mPBCBackElements);

        ///////////////////////////////////////////////    Debugging     ///////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::map<unsigned int, int>::iterator it_nn;
		//for(it_nn = mPBC_NodeNodePairs_Front.begin(); it_nn != mPBC_NodeNodePairs_Front.end(); ++it_nn)
		//{
		//    std::cout << "Frontal: " << it_nn->first << " -> ";
        //    std::cout << it_nn->second << std::endl;
		//}
        //
		//std::map<unsigned int, std::vector<int> >::iterator it_neg;
		//std::map<unsigned int, double>::iterator it_edco = mPBC_NodeEdgeCoord_Front.begin();
		//for(it_neg = mPBC_NodeEdgePairs_Front.begin(); it_neg != mPBC_NodeEdgePairs_Front.end(); ++it_neg, ++it_edco)
		//{
		//    std::cout << "Frontal: " << it_neg->first;
		//
		//    std::cout << " -> " << (it_neg->second)[0] << " , " << (it_neg->second)[1];
        //    if      (mElementOrder == 2)
        //    std::cout << " , "  << (it_neg->second)[2];
        //    else if (mElementOrder == 3)
		//    std::cout << " , "  << (it_neg->second)[2]
        //              << " , "  << (it_neg->second)[3];
        //    else if (mElementOrder == 4)
		//    std::cout << " , "  << (it_neg->second)[2]
        //              << " , "  << (it_neg->second)[3]
        //              << " , "  << (it_neg->second)[4];
        //
		//    std::cout << " - FNatX: " << std::setprecision(12) << (it_edco->second) << "   " << std::endl;
		//}
        //
		//std::map<unsigned int, std::vector< int  > >::iterator it_nel;
		//std::map<unsigned int, std::vector<double> >::iterator it_elco =  mPBC_NodeElementCoord_Front.begin();
		//for(it_nel = mPBC_NodeElementPairs_Front.begin(); it_nel != mPBC_NodeElementPairs_Front.end(); ++it_nel, ++it_elco)
		//{
		//    std::cout << "Frontal: " << it_nel->first;
        //
		//    std::cout << " -> " << (it_nel->second)[0] << " , " << (it_nel->second)[1] << " , " << (it_nel->second)[2];
        //    if      (mElementOrder == 2)
        //    std::cout << " , "  << (it_nel->second)[3] << " , " << (it_nel->second)[4] << " , " << (it_nel->second)[5];
        //    else if (mElementOrder == 3)
		//    std::cout << " , " << (it_nel->second)[3] << " , " << (it_nel->second)[4]
        //              << " , " << (it_nel->second)[5] << " , " << (it_nel->second)[6]
        //              << " , " << (it_nel->second)[7] << " , " << (it_nel->second)[8];
        //    else if (mElementOrder == 4)
		//    std::cout << " , " << (it_nel->second)[3] << " , " << (it_nel->second)[ 4] << " , " << (it_nel->second)[ 5]
        //              << " , " << (it_nel->second)[6] << " , " << (it_nel->second)[ 7] << " , " << (it_nel->second)[ 8]
        //              << " , " << (it_nel->second)[9] << " , " << (it_nel->second)[10] << " , " << (it_nel->second)[11];
        //
		//    std::cout << " - FNatXY: " << std::setprecision(12) << (it_elco->second)[0] << " , " << (it_elco->second)[1]
        //              << std::endl;
		//}
        //Send_Error_Msg ( "Front-Back PBC debugging finished", "Please, check PBC pairs.", 1 );
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Calculates the parameters for the geometrical mapping
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::PBC_Mapping_Parameters( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY )
    {
        if  ( mQuadraticGeometry == true )
        {
            PBC_Mapping_Parameters_Quadrt( a, b, eX, eY );
        }
        else
        {
            PBC_Mapping_Parameters_Linear( a, b, eX, eY );
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Calculates the parameters for the geometrical mapping of a linear interpolation
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::PBC_Mapping_Parameters_Linear( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY )
    {
        // Mapping parameters (linear mapping)
        a.resize(3);
        b.resize(3);

        double cte = 1.0 / ( (eY[1]-eY[2]) * (eX[0]-eX[2]) - (eY[2]-eY[0]) * (eX[2]-eX[1]) );

        a[0] = cte * ( eX[2] * eY[0] - eY[2] * eX[0] );
        a[1] = cte * ( eY[2] - eY[0] );
        a[2] = cte * ( eX[0] - eX[2] );

        b[0] = cte * ( eX[0] * eY[1] - eY[0] * eX[1] );
        b[1] = cte * ( eY[0] - eY[1] );
        b[2] = cte * ( eX[1] - eX[0] );
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Calculates the parameters for the geometrical mapping of a quadratic interpolation
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void Modeler::PBC_Mapping_Parameters_Quadrt( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY )
    {
        // Coeficient matrix for a second order mapping (see Jin pag. 157)
        Matrix<double> Mpp(6,6);

        for ( int i=0; i<6 ; i++ )
        {
            Mpp(i,0) = 1.00       ;
            Mpp(i,1) = eX[i]      ;
            Mpp(i,2) = eX[i]*eX[i];
            Mpp(i,3) = eX[i]*eY[i];
            Mpp(i,4) = eY[i]      ;
            Mpp(i,5) = eY[i]*eY[i];
        }

        double detMpp = Determinant( Mpp, 6 );

        // Residual vectors
        double fa[6];

        fa[0] = 0.0;
        fa[1] = 1.0;
        fa[2] = 0.0;
        fa[3] = 0.5;
        fa[4] = 0.5;
        fa[5] = 0.0;

        double fb[6];

        fb[0] = 0.0;
        fb[1] = 0.0;
        fb[2] = 1.0;
        fb[3] = 0.0;
        fb[4] = 0.5;
        fb[5] = 0.5;

        // Mapping parameters (second order mapping)
        a.resize(6);
        b.resize(6);

        // Adjoint matrix
        Matrix<double> Adj(6,6);

        // Solving the linear system for the mapping parameters
        for ( int j=0; j<6 ; j++ )
        {
            // Adjoint matrix for a[j]
            Adj = Mpp; for (int i=0; i<6 ; i++) Adj(i,j) = fa[i];

            a[j] = Determinant( Adj, 6 ) / detMpp;

            // Adjoint matrix for b[j]
            Adj = Mpp; for (int i=0; i<6 ; i++) Adj(i,j) = fb[i];

            b[j] = Determinant( Adj, 6 ) / detMpp;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Returns the determinant of a square matrix M of size n
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double Modeler::Determinant( Matrix<double>& M, int n )
    {
        double det = 0.0;

        Matrix<double> subM(n,n);

        if (n == 1)
        {
            return M[0][0];
        }

        if (n == 2)
        {
            return ( (M[0][0] * M[1][1]) - (M[1][0] * M[0][1]) );
        }

        for ( int x = 0; x < n; x++ )
        {
            int subi = 0;
            for ( int i = 1; i < n; i++ )
            {
                int subj = 0;
                for ( int j = 0; j < n; j++ )
                {
                    if (j == x) continue;

                    subM[subi][subj] = M[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + ( pow( -1, x ) * M[0][x] * Determinant( subM, n-1 ) );
        }

        return det;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Checks if the projected node with natural coordinates (natX, natY) is on a node of the element ElementNodeIds.
    // - If the projected node projects on a node it returns its Node_Id > 0.
    // - If the projected node is not any of the element nodes it returns 0.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int Modeler::PBC_Node_To_Node( double natX, double natY, std::vector<int>& ElementNodeIds )
    {
        // Nodes on the vertices (valid for all element order)
        if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) )     { return ElementNodeIds[0]; }
        if ( ( abs(natX - 1.00) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) )     { return ElementNodeIds[1]; }
        if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY - 1.00) < mPBCTolerance ) )     { return ElementNodeIds[2]; }

        // Bubble node in face
        if      ( mElementOrder == 0 )
        {
            if ( ( abs(natX - 1.00/3.00) < mPBCTolerance ) && ( abs(natY - 1.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[3]; }
        }
        // 2nd order elements middle nodes
        else if ( mElementOrder == 2 )
        {
            if ( ( abs(natX - 0.50) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) ) { return ElementNodeIds[3]; }
            if ( ( abs(natX - 0.50) < mPBCTolerance ) && ( abs(natY - 0.50) < mPBCTolerance ) ) { return ElementNodeIds[4]; }
            if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY - 0.50) < mPBCTolerance ) ) { return ElementNodeIds[5]; }
        }
        // 3rd order elements middle nodes
        else if ( mElementOrder == 3 )
        {
            if ( ( abs(natX - 1.00/3.00) < mPBCTolerance ) && ( abs(natY            ) < mPBCTolerance ) ) { return ElementNodeIds[3]; }
            if ( ( abs(natX - 2.00/3.00) < mPBCTolerance ) && ( abs(natY            ) < mPBCTolerance ) ) { return ElementNodeIds[4]; }

            if ( ( abs(natX - 2.00/3.00) < mPBCTolerance ) && ( abs(natY - 1.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[5]; }
            if ( ( abs(natX - 1.00/3.00) < mPBCTolerance ) && ( abs(natY - 2.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[6]; }

            if ( ( abs(natX            ) < mPBCTolerance ) && ( abs(natY - 2.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[7]; }
            if ( ( abs(natX            ) < mPBCTolerance ) && ( abs(natY - 1.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[8]; }

            if ( ( abs(natX - 1.00/3.00) < mPBCTolerance ) && ( abs(natY - 1.00/3.00) < mPBCTolerance ) ) { return ElementNodeIds[9]; }
        }
        // 4th order elements middle nodes
        else if ( mElementOrder == 4 )
        {
            if ( ( abs(natX - 0.25) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) ) { return ElementNodeIds[3 ]; }
            if ( ( abs(natX - 0.50) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) ) { return ElementNodeIds[4 ]; }
            if ( ( abs(natX - 0.75) < mPBCTolerance ) && ( abs(natY       ) < mPBCTolerance ) ) { return ElementNodeIds[5 ]; }

            if ( ( abs(natX - 0.75) < mPBCTolerance ) && ( abs(natY - 0.25) < mPBCTolerance ) ) { return ElementNodeIds[6 ]; }
            if ( ( abs(natX - 0.50) < mPBCTolerance ) && ( abs(natY - 0.50) < mPBCTolerance ) ) { return ElementNodeIds[7 ]; }
            if ( ( abs(natX - 0.25) < mPBCTolerance ) && ( abs(natY - 0.75) < mPBCTolerance ) ) { return ElementNodeIds[8 ]; }

            if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY - 0.75) < mPBCTolerance ) ) { return ElementNodeIds[9 ]; }
            if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY - 0.50) < mPBCTolerance ) ) { return ElementNodeIds[10]; }
            if ( ( abs(natX       ) < mPBCTolerance ) && ( abs(natY - 0.25) < mPBCTolerance ) ) { return ElementNodeIds[11]; }

            if ( ( abs(natX - 0.25) < mPBCTolerance ) && ( abs(natY - 0.25) < mPBCTolerance ) ) { return ElementNodeIds[12]; }
            if ( ( abs(natX - 0.50) < mPBCTolerance ) && ( abs(natY - 0.25) < mPBCTolerance ) ) { return ElementNodeIds[13]; }
            if ( ( abs(natX - 0.25) < mPBCTolerance ) && ( abs(natY - 0.50) < mPBCTolerance ) ) { return ElementNodeIds[14]; }
        }

        return 0;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Checks if the projected node with natural coordinates (natX, natY) is on an edge of the element ElementNodeIds.
    // - If the projected node is on an edge then it returns its natural coordinate to respect that edge.
    // - If the projected node is not on an edge then it returns 0.0.
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double Modeler::PBC_Node_To_Edge( double natX, double natY, std::vector<int>& ElementNodeIds, std::vector<int>& EdgeNodeIds )
    {
        // Projected node is on the edge [0 1]
        if ( abs( natY ) < mPBCTolerance )
        {
            EdgeNodeIds.push_back(ElementNodeIds[0]);
            EdgeNodeIds.push_back(ElementNodeIds[1]);

            // 2nd order element middle nodes
            if      ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[3]);
            }
            // 3rd order element middle nodes
            else if ( mElementOrder == 3 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[3]);
                EdgeNodeIds.push_back(ElementNodeIds[4]);
            }
            // 4th order element middle nodes
            else if ( mElementOrder == 4 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[3]);
                EdgeNodeIds.push_back(ElementNodeIds[4]);
                EdgeNodeIds.push_back(ElementNodeIds[5]);
            }

            return natX;
        }
        // Projected node is on the edge [1 2]
        else if ( abs( natX + natY - 1.00 ) < ( sqrt(2) * mPBCTolerance ) )
        {
            EdgeNodeIds.push_back(ElementNodeIds[1]);
            EdgeNodeIds.push_back(ElementNodeIds[2]);

            // 2nd order element middle nodes
            if      ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[4]);
            }
            // 3rd order element middle nodes
            else if ( mElementOrder == 3 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[5]);
                EdgeNodeIds.push_back(ElementNodeIds[6]);
            }
            // 4th order element middle nodes
            else if ( mElementOrder == 4 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[6]);
                EdgeNodeIds.push_back(ElementNodeIds[7]);
                EdgeNodeIds.push_back(ElementNodeIds[8]);
            }

            return ( 1.00 - natX );
        }
        // Projected node is on the edge [2 0]
        else if ( abs( natX ) < mPBCTolerance )
        {
            EdgeNodeIds.push_back(ElementNodeIds[2]);
            EdgeNodeIds.push_back(ElementNodeIds[0]);

            // 2nd order element middle nodes
            if      ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[5]);
            }
            // 3rd order element middle nodes
            else if ( mElementOrder == 3 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[7]);
                EdgeNodeIds.push_back(ElementNodeIds[8]);
            }
            // 4th order element middle nodes
            else if ( mElementOrder == 4 )
            {
                EdgeNodeIds.push_back(ElementNodeIds[ 9]);
                EdgeNodeIds.push_back(ElementNodeIds[10]);
                EdgeNodeIds.push_back(ElementNodeIds[11]);
            }

            return ( 1.00 - natY );
        }
        else
        {
            return 0.0;
        }
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
    // - Assigns the nodes on the right periodic surface to the elements on the left periodic surface.
	// - Left surface must be placed in the 00-YZ-Plane.
	// - For cylindrical symmetry (cyclic periodic) the central axis must be placed along the Z axis.
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Modeler::Set_RightLeft_PBC()
	{
         // Loop through materials in Left PBC
		std::map<unsigned int, std::vector<std::vector<int> > >::iterator it_Left;
        for(it_Left = mPBCLeftElements.begin(); it_Left != mPBCLeftElements.end(); ++it_Left)
        {
            // Material ID
			int pbcMaterial = it_Left->first;

            // Group of Left elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetLeft( it_Left->second );

            // Group of Right elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetRight( mPBCRightElements[pbcMaterial] );

            // List of nodes (sorted and unique) in the Right PBC
            std::list<int> RightNodesList;

            // Filling the list with the Right elements nodes and checking if the Right surface is tilted
            std::vector<std::vector<int> >::iterator it_eR;
            for( it_eR = elemSetRight.begin(); it_eR != elemSetRight.end(); ++it_eR )
            {
                // Filling the list with the Right elements nodes
                std::vector<int>::iterator it_nR;
                for( it_nR = (*it_eR).begin(); it_nR != (*it_eR).end(); ++it_nR )
				{
                    RightNodesList.push_back( *it_nR );
                }

                // Checks if the Right surface is tilted or not (Periodic/Cyclic boundary conditions)
                double x[3]; for ( int i=0; i<3; i++ )
                {
                    x[i] = ( mpModel->GetNode( (*it_eR)[i]-1 ) )->X();
                }
                mIsRightPBCTilted = ( abs(x[0]-x[1]) > mPBCTolerance ) || ( abs(x[1]-x[2]) > mPBCTolerance ) || ( abs(x[2]-x[0]) > mPBCTolerance );
            }

            // Sorting list
            RightNodesList.sort();

            // Collapsing equal nodes Id
            RightNodesList.unique();

            // Loop through the Left elements to calculate the Right nodes projections
			std::vector<std::vector<int> >::iterator it_eL;
			for( it_eL = elemSetLeft.begin(); it_eL != elemSetLeft.end() && RightNodesList.size() > 0 ; ++it_eL )
            {
                // Recovering Left element nodes coordinates on the ZY plane
                std::vector<double> eX, eY;

                for ( int i=0; i<(*it_eL).size(); i++ )
                {
                    eX.push_back( mpModel->GetNode( (*it_eL)[i]-1 )->Z() );
                    eY.push_back( mpModel->GetNode( (*it_eL)[i]-1 )->Y() );
                }

                // Calculating cross section (only nodes inside this area will be consider)
                double eCross_Section;
                double edge01_Lenght = sqrt( pow(eX[0]-eX[1], 2) + pow(eY[0]-eY[1], 2) );
                double edge02_Lenght = sqrt( pow(eX[0]-eX[2], 2) + pow(eY[0]-eY[2], 2) );

                if ( edge02_Lenght > edge01_Lenght ) eCross_Section = 1.1 * edge02_Lenght;
                else                                 eCross_Section = 1.1 * edge01_Lenght;

                // Parameters of the geometrical mapping interpolation
                std::vector<double> a, b;

                // Computing the parameters a,b of the geometrical mapping interpolation
                PBC_Mapping_Parameters( a, b, eX, eY );

                // Right nodes projected on this Left element
                std::vector<int> AssignedRightNodes;

                // Loop over the list of nodes in the Right PBC
                std::list<int>::iterator it_nRl;
                for( it_nRl = RightNodesList.begin(); it_nRl != RightNodesList.end(); ++it_nRl )
                {
                    // Recovering the node pointer from its Id
                    Node::Pointer prNode = mpModel->GetNode( (*it_nRl) - 1 );

                    // Cartesina coordinates of the node (2D projection)
                    double carX = prNode->Z();
                    double carY = PBC_RotateY(prNode);

                    // Checking if the Right node projection is close enough to the Left element
                    if ( sqrt( pow(eX[0]-carX, 2) + pow(eY[0]-carY, 2) ) > eCross_Section ) continue;

                    // Natural coordinates of the Frontal node projection on the Back element
                    double natX;
                    double natY;

                    if ( a.size() == 6 )
                    {
                        natX = a[0] + ( a[1] * carX ) + ( a[2] * carX * carX ) + ( a[3] * carX * carY ) + ( a[4] * carY ) + ( a[5] * carY * carY );
                        natY = b[0] + ( b[1] * carX ) + ( b[2] * carX * carX ) + ( b[3] * carX * carY ) + ( b[4] * carY ) + ( b[5] * carY * carY );
                    }
                    else
                    {
                        natX = a[0] + ( a[1] * carX ) + ( a[2] * carY );
                        natY = b[0] + ( b[1] * carX ) + ( b[2] * carY );
                    }

                    // If the Right node is outside the Left element -> continue
                    if ( ( natX < (-mPBCTolerance) ) || ( natY < (-mPBCTolerance) ) || ( (natX + natY) > (1.00 + sqrt(2)*mPBCTolerance) ) ) continue;

                    // This node is going to be assigned
                    AssignedRightNodes.push_back( *it_nRl );

                    // Setting DOFs for this node
                    std::complex<double> cZero(0.00,0.00);

                    // If the Right surface is tilted and the node is in the Z-axis then make Ex=Ey=0 and continue
                    if ( ( mIsRightPBCTilted == true ) && ( abs(prNode->X()) < mPBCTolerance ) && ( abs(prNode->Y()) < mPBCTolerance ) )
                    {
                        if (mNormals.find(*it_nRl) == mNormals.end())
                        {
                            (*mpModel)(cEx, *prNode) = cZero; prNode->pDofcEx()->FixDof();
                            (*mpModel)(cEy, *prNode) = cZero; prNode->pDofcEy()->FixDof();
                        }
                        continue;
                    }

                    // Right node FixDof
                    (*mpModel)(cEx, *prNode) = cZero; prNode->pDofcEx()->FixDof();
                    (*mpModel)(cEy, *prNode) = cZero; prNode->pDofcEy()->FixDof();
                    (*mpModel)(cEz, *prNode) = cZero; prNode->pDofcEz()->FixDof();

                    // Cheking if the Right node is projected onto a Left node
                    int LeftNodeId = PBC_Node_To_Node( natX, natY, *it_eL );

                    if ( LeftNodeId > 0 )
                    {
                        mPBC_NodeNodePairs_Right[prNode->Id()] = LeftNodeId;
                        continue;
                    }

                    // Cheking if Right node is projected onto a Left element edge
                    std::vector<int> LeftEdgeNodeIds;

                    double NatCoord_InEdge = PBC_Node_To_Edge( natX, natY, *it_eL, LeftEdgeNodeIds );

                    if ( LeftEdgeNodeIds.size() > 0 )
                    {
                        mPBC_NodeEdgePairs_Right[prNode->Id()] = LeftEdgeNodeIds;
                        mPBC_NodeEdgeCoord_Right[prNode->Id()] = NatCoord_InEdge;
                        continue;
                    }

                    // Right node is projected onto a Left element
                    std::vector<double> NatCoord_InElement(2);
                    NatCoord_InElement[0] = natX;
                    NatCoord_InElement[1] = natY;

                    mPBC_NodeElementPairs_Right[prNode->Id()] = (*it_eL);
                    mPBC_NodeElementCoord_Right[prNode->Id()] = NatCoord_InElement;

                }// End loop over the list of nodes in the Right PBC

                // Removing from the RightNodesList the assigned Right nodes
                std::vector<int>::iterator ita;
                for(ita = AssignedRightNodes.begin(); ita != AssignedRightNodes.end(); ++ita) RightNodesList.remove(*ita);

            }// End loop through Left elements

            // Checking if all the Right nodes in a pcbMaterial have been assigned
            if ( RightNodesList.size() > 0 )
            {
                std::cout << RightNodesList.size() << " Right nodes are not assigned in PBC material " << pbcMaterial << std::endl;
            }

        }// End loop through materials in Left PBC

		mPBCRightElements.clear();
		std::map<unsigned int, std::vector<std::vector<int> > >().swap(mPBCRightElements);

		mPBCLeftElements.clear();
		std::map<unsigned int, std::vector<std::vector<int> > >().swap(mPBCLeftElements);

        ///////////////////////////////////////////////    Debugging     ///////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::map<unsigned int, int>::iterator it_nn;
		//for(it_nn = mPBC_NodeNodePairs_Right.begin(); it_nn != mPBC_NodeNodePairs_Right.end(); ++it_nn)
		//{
		//    std::cout << "Right: " << it_nn->first << " -> ";
        //    std::cout << it_nn->second << std::endl;
		//}
        //
		//std::map<unsigned int, std::vector<int> >::iterator it_neg;
		//std::map<unsigned int, double>::iterator it_edco = mPBC_NodeEdgeCoord_Right.begin();
		//for(it_neg = mPBC_NodeEdgePairs_Right.begin(); it_neg != mPBC_NodeEdgePairs_Right.end(); ++it_neg, ++it_edco)
		//{
		//    std::cout << "Right: " << it_neg->first;
		//
		//    std::cout << " -> "<< (it_neg->second)[0] << " , " << (it_neg->second)[1];
        //    if      (mElementOrder == 2)
        //    std::cout << " , " << (it_neg->second)[2];
        //    else if (mElementOrder == 3)
		//    std::cout << " , " << (it_neg->second)[2]
        //              << " , " << (it_neg->second)[3];
        //    else if (mElementOrder == 4)
		//    std::cout << " , " << (it_neg->second)[2]
        //              << " , " << (it_neg->second)[3]
        //              << " , " << (it_neg->second)[4];
        //
		//    std::cout << " - FNatX: " << std::setprecision(12) << (it_edco->second) << "   " << std::endl;
		//}
        //
		//std::map<unsigned int, std::vector< int  > >::iterator it_nel;
		//std::map<unsigned int, std::vector<double> >::iterator it_elco =  mPBC_NodeElementCoord_Right.begin();
		//for(it_nel = mPBC_NodeElementPairs_Right.begin(); it_nel != mPBC_NodeElementPairs_Right.end(); ++it_nel, ++it_elco)
		//{
		//    std::cout << "Right: " << it_nel->first;
        //
		//    std::cout << " -> "<< (it_nel->second)[0] << " , " << (it_nel->second)[1] << " , " << (it_nel->second)[2];
        //    if      (mElementOrder == 2)
        //    std::cout << " , " << (it_nel->second)[3] << " , " << (it_nel->second)[4] << " , " << (it_nel->second)[5];
        //    else if (mElementOrder == 3)
		//    std::cout << " , " << (it_nel->second)[3] << " , " << (it_nel->second)[4]
        //              << " , " << (it_nel->second)[5] << " , " << (it_nel->second)[6]
        //              << " , " << (it_nel->second)[7] << " , " << (it_nel->second)[8];
        //    else if (mElementOrder == 4)
		//    std::cout << " , " << (it_nel->second)[3] << " , " << (it_nel->second)[ 4] << " , " << (it_nel->second)[ 5]
        //              << " , " << (it_nel->second)[6] << " , " << (it_nel->second)[ 7] << " , " << (it_nel->second)[ 8]
        //              << " , " << (it_nel->second)[9] << " , " << (it_nel->second)[10] << " , " << (it_nel->second)[11];
        //
		//    std::cout << " - FNatXY: " << std::setprecision(12) << (it_elco->second)[0] << " , " << (it_elco->second)[1]
        //              << std::endl;
		//}
        //Send_Error_Msg ( "Left-Right PBC debugging finished", "Please, check PBC pairs.", 1 );
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Rotates Right node around Z-axis from Right plane to Left plane (00-YZ-plane)
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	double Modeler::PBC_RotateY(Node::Pointer prNode)
	{
		if ( mIsRightPBCTilted == false ) return ( prNode->Y() );

		double X = prNode->X();
		double Y = prNode->Y();

		double modXY2 = X*X + Y*Y;

		return ( modXY2 / sqrt(modXY2) );
	}

	//****************************************************************************************************************
    //* - Calculates the volume integral of the fields
    //****************************************************************************************************************
    void Modeler::GenerateVolIntElement( int* NodesId, unsigned int PropertiesId )
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume( NodesId, HONodesId );

		std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double volumeID = cIBC[0];

		int gp;

		std::vector<double> cX;
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order6( cX, cY, cZ, W );

		Matrix<double> N;

        if      ( mElementOrder == 0 ) Lagrange3D_Ni_2bb( N, cX, cY, cZ );
		else if ( mElementOrder == 1 ) Lagrange3D_Ni_1st( N, cX, cY, cZ );
        else if ( mElementOrder == 2 ) Lagrange3D_Ni_2nd( N, cX, cY, cZ );
        else if ( mElementOrder == 3 ) Lagrange3D_Ni_3th( N, cX, cY, cZ );
		else if ( mElementOrder == 4 ) Lagrange3D_Ni_4th( N, cX, cY, cZ );

		double elementVolume = 0.00;

		if ( mQuadraticGeometry )
		{
			std::vector<double> detJ;

			Calculate_detJ( detJ, cX, cY, cZ, nodes );

			for ( gp=0; gp<nGaussPoints; ++gp ) W[gp]         *= detJ[gp];
			for ( gp=0; gp<nGaussPoints; ++gp ) elementVolume +=    W[gp];
		}
		else
		{
			elementVolume = VolumeOfTetrahedra(nodes);

			double jacob  = 6.00 * elementVolume;

			for ( gp=0; gp<nGaussPoints; ++gp ) W[gp] *= jacob;
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

        for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
		    double vIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; ++gp ) vIntNi += W[gp] * N[in][gp];

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

    //***********************************************************************************
    //* - Calculation of dN/da, dN/db, dN/du on cX,cY,cZ points
    //***********************************************************************************
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

    //**********************************************************************************
    //* - Calculation of the Jacobian determinant
    //**********************************************************************************
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

		PushHONodesOnSurface( NodesId, HONodesId );

		std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double surfaceID = cIBC[0];

		int gp;

		std::vector<double> cX;
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order8( cX, cY, W );

		Matrix<double> N;

        if      ( mElementOrder == 0 ) Lagrange2D_Ni_2bb( N, cX, cY );
		else if ( mElementOrder == 1 ) Lagrange2D_Ni_1st( N, cX, cY );
        else if ( mElementOrder == 2 ) Lagrange2D_Ni_2nd( N, cX, cY );
        else if ( mElementOrder == 3 ) Lagrange2D_Ni_3th( N, cX, cY );
		else if ( mElementOrder == 4 ) Lagrange2D_Ni_4th( N, cX, cY );

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

    //*************************************************************************************
    //* - Calculates volumetric integrals of the fields in static mode.
    //*************************************************************************************
    void Modeler::GenerateVolIntElement_Static( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId = 0 then the volume element is used for getting field gradients.
        // These gradients are used in the surface integrals of non-smoothed fields.
        if ( PropertiesId == 0 )
        {
            GetGrad_SurfInt_NonSmooth_Static( NodesId );
            return;
        }

        // Getting all element nodes
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume(NodesId, HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        // Getting volume Id
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Vector<double> cIBC = (*properties)(COMPLEX_IBC);

        double volumeID = cIBC[0];

        // Integration Gauss points
        int gp; std::vector<double> cX, cY, cZ, W;

        int nGaussPoints = GaussPoints3D_Order6( cX, cY, cZ, W );

        Matrix<double> N;

        Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // Volume of the element
        double elementVolume = VolumeOfTetrahedra(nodes);
        double jacob         = 6.00 * elementVolume;

        for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob;

        // Volumetric integrals vectors
        Vector<double> volIntegralE  (5, 0.00);
        Vector<double> volIntegralJ  (5, 0.00);
        Vector<double> volIntegralEgp(5, 0.00);

        // Getting gradients to calculate electric field non-smoothed
        Vector< Vector<double> > gradOnNodes;

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        pElement->CalculateGradientOnNodes( gradOnNodes );

        // Volumetric nodal integration of the fields
        std::vector<Node::Pointer>::iterator node_it;

        // Null vector
        Vector<double> vCero(3, 0.00);

        int in = 0;

        for (node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
            double vIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ ) vIntNi += W[gp] * N[in][gp];

            // Electric field
            Vector<double> nodalE(3, 0.00);

            // Catch exception when ELECTRIC_FIELD_SMOOTHED is not calculated
            try       { nodalE = (*properties)(ELECTRIC_FIELD_SMOOTHED, **node_it); }
            catch(...){ nodalE = vCero; }

            double modE2 = nodalE[0]*nodalE[0] + nodalE[1]*nodalE[1] + nodalE[2]*nodalE[2];
            double modE  = sqrt(modE2);

            // Volume integral of Ex, Ey, Ez
            volIntegralE[0] += nodalE[0] * vIntNi;
            volIntegralE[1] += nodalE[1] * vIntNi;
            volIntegralE[2] += nodalE[2] * vIntNi;

            // Volume integral of |E|^2, |E|
            volIntegralE[3] += modE2 * vIntNi;
            volIntegralE[4] += modE  * vIntNi;

            // Current density
            Vector<double> nodalJ(3, 0.00);

            // Catch exception when CURRENT_DENSITY_SMOOTHED is not calculated
            try       { nodalJ = (*properties)(CURRENT_DENSITY_SMOOTHED, **node_it); }
            catch(...){ nodalJ = vCero; }

            double modJ2 = nodalJ[0]*nodalJ[0] + nodalJ[1]*nodalJ[1] + nodalJ[2]*nodalJ[2];
            double modJ  = sqrt(modJ2);

            // Volume integral of Jx, Jy, Jz
            volIntegralJ[0] += nodalJ[0] * vIntNi;
            volIntegralJ[1] += nodalJ[1] * vIntNi;
            volIntegralJ[2] += nodalJ[2] * vIntNi;

            // Volume integral of |J|^2, |J|
            volIntegralJ[3] += modJ2 * vIntNi;
            volIntegralJ[4] += modJ  * vIntNi;

            // Electric field non-smoothed
            Vector<double> nodalEgp = gradOnNodes[in];

            // E = -grad(V)
            nodalEgp *= -1.00;

            double modEgp2 = nodalEgp[0]*nodalEgp[0] + nodalEgp[1]*nodalEgp[1] + nodalEgp[2]*nodalEgp[2];
            double modEgp  = sqrt(modEgp2);

            // Volume integral of Ex, Ey, Ez non-smoothed
            volIntegralEgp[0] += nodalEgp[0] * vIntNi;
            volIntegralEgp[1] += nodalEgp[1] * vIntNi;
            volIntegralEgp[2] += nodalEgp[2] * vIntNi;

            // Volume integral of |E|^2, |E| non-smoothed
            volIntegralEgp[3] += modEgp2 * vIntNi;
            volIntegralEgp[4] += modEgp  * vIntNi;

            in++;
        }

        // Adding volume integrals of this element to the total
        if ( mVolumeValue.find(volumeID) != mVolumeValue.end() )
        {
            mVolumeValue           [volumeID] += elementVolume;
            mVolIntegralE_Static   [volumeID] += volIntegralE;
            mVolIntegralJ_Static   [volumeID] += volIntegralJ;
            mVolIntegralE_Static_GP[volumeID] += volIntegralEgp;
        }
        else
        {
            mVolumeValue           [volumeID] = elementVolume;
            mVolIntegralE_Static   [volumeID] = volIntegralE;
            mVolIntegralJ_Static   [volumeID] = volIntegralJ;
            mVolIntegralE_Static_GP[volumeID] = volIntegralEgp;
        }
    }

    //***********************************************************************************************
    //* - Calculates gradients for the non-smoothed fields surface integrals.
    //*   If PropertiesId = 0 then the PVIE element is used for getting field gradients.
    //*   These gradients are used in the surface integrals of non-smoothed fields.
    //***********************************************************************************************
    void Modeler::GetGrad_SurfInt_NonSmooth_Static( int* NodesId )
    {
        // Getting all element nodes
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume(NodesId, HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        // Dummy property, no materials are used for the gradients
        Properties::Pointer properties = mpModel->GetProperties(1);

        // Setting element
        Element::Pointer pElement = Element::Pointer(new VolumeElement_1st_Static(nodes,properties));

        // Calculate gradients on nodes (non-smoothed gradients are constant in 1st order elements)
        Vector< Vector<double> > gradOnNodes;

        pElement->CalculateGradientOnNodes( gradOnNodes );

        // Add gradient of first node to the mSurfInt_NonSmooth_E_Static vector
        mSurfInt_NonSmooth_E_Static.push_back( gradOnNodes[0] );
    }

    //***********************************************************************************
    //* - Calculates surface integral of the fields for the static case.
    //***********************************************************************************
    void Modeler::GenerateSurfIntElement_Static( int* NodesId, unsigned int PropertiesId )
    {
        // Getting all element nodes
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnSurface(NodesId, HONodesId);

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        // Getting surface Id
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Vector<double> cIBC = (*properties)(COMPLEX_IBC);

        double surfaceID = cIBC[0];

        // Integration Gauss points
        int gp; std::vector<double> cX, cY, W;

        int nGaussPoints = GaussPoints2D_Order8(cX, cY, W);

        Matrix<double> N;

        Lagrange2D_Ni_1st(N, cX, cY);

        // Surface of the element
        double elementArea = AreaOfTriangle(nodes);
        double jacob       = 2.00 * elementArea;

        for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob;

        // Surface integrals vectors
        Vector<double> surfIntegralE  (5, 0.00);
        Vector<double> surfIntegralJ  (5, 0.00);
        Vector<double> surfIntegralEgp(5, 0.00);

        // Getting gradients to calculate electric field non-smoothed
        Vector<double> gradOnElement = mSurfInt_NonSmooth_E_Static[mNonSmooth_Surf_Id];

        // Surface nodal integration of the fields
        std::vector<Node::Pointer>::iterator node_it;

        // Null vector
        Vector<double> vCero(3, 0.00);

        int in = 0;

        for ( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            double sIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ ) sIntNi += W[gp]*N[in][gp];

            // Electric field
            Vector<double> nodalE(3, 0.00);

            // Catch exception when ELECTRIC_FIELD_SMOOTHED is not calculated
            try       { nodalE = (*properties)(ELECTRIC_FIELD_SMOOTHED, **node_it); }
            catch(...){ nodalE = vCero; }

            double modE2 = nodalE[0]*nodalE[0] + nodalE[1]*nodalE[1] + nodalE[2]*nodalE[2];
            double modE  = sqrt(modE2);

            // Surface integral of Ex, Ey, Ez
            surfIntegralE[0] += nodalE[0] * sIntNi;
            surfIntegralE[1] += nodalE[1] * sIntNi;
            surfIntegralE[2] += nodalE[2] * sIntNi;

            // Surface integral of |E|^2, |E|
            surfIntegralE[3] += modE2 * sIntNi;
            surfIntegralE[4] += modE  * sIntNi;

            // Current density
            Vector<double> nodalJ(3, 0.00);

            // Catch exception when CURRENT_DENSITY_SMOOTHED is not calculated
            try       { nodalJ = (*properties)(CURRENT_DENSITY_SMOOTHED, **node_it); }
            catch(...){ nodalJ = vCero; }

            double modJ2 = nodalJ[0]*nodalJ[0] + nodalJ[1]*nodalJ[1] + nodalJ[2]*nodalJ[2];
            double modJ  = sqrt(modJ2);

            // Volume integral of Jx, Jy, Jz
            surfIntegralJ[0] += nodalJ[0] * sIntNi;
            surfIntegralJ[1] += nodalJ[1] * sIntNi;
            surfIntegralJ[2] += nodalJ[2] * sIntNi;

            // Volume integral of |J|^2, |J|
            surfIntegralJ[3] += modJ2 * sIntNi;
            surfIntegralJ[4] += modJ  * sIntNi;

            // Electric field non-smoothed
            Vector<double> nodalEgp = gradOnElement;

            // E = -grad(V)
            nodalEgp *= -1.00;

            double modEgp2 = nodalEgp[0]*nodalEgp[0] + nodalEgp[1]*nodalEgp[1] + nodalEgp[2]*nodalEgp[2];
            double modEgp  = sqrt(modEgp2);

            // Surface integral of Ex, Ey, Ez non-smoothed
            surfIntegralEgp[0] += nodalEgp[0] * sIntNi;
            surfIntegralEgp[1] += nodalEgp[1] * sIntNi;
            surfIntegralEgp[2] += nodalEgp[2] * sIntNi;

            // Surface integral of |E|^2, |E| non-smoothed
            surfIntegralEgp[3] += modEgp2 * sIntNi;
            surfIntegralEgp[4] += modEgp  * sIntNi;

            in++;
        }

        // Adding surface integrals of this element to the total
        if ( mSurfaceValue.find(surfaceID) != mSurfaceValue.end() )
        {
            mSurfaceValue           [surfaceID] += elementArea;
            mSurfIntegralE_Static   [surfaceID] += surfIntegralE;
            mSurfIntegralJ_Static   [surfaceID] += surfIntegralJ;
            mSurfIntegralE_Static_GP[surfaceID] += surfIntegralEgp;
        }
        else
        {
            mSurfaceValue           [surfaceID] = elementArea;
            mSurfIntegralE_Static   [surfaceID] = surfIntegralE;
            mSurfIntegralJ_Static   [surfaceID] = surfIntegralJ;
            mSurfIntegralE_Static_GP[surfaceID] = surfIntegralEgp;
        }

        // Increase surface_id counter
        mNonSmooth_Surf_Id++;
    }

	//*******************************************************************************************************
    //* - Projecting resulting field on plane
    //*******************************************************************************************************
    void Modeler::GenerateProjectionRWPortTE10(int*  NodesId, unsigned int PropertiesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new RWPortTE10_1st( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new RWPortTE10_1st( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new RWPortTE10_2nd( nodes, properties ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new RWPortTE10_3th( nodes, properties ) );
		else if ( mElementOrder == 4 ) pElement = Element::Pointer( new RWPortTE10_4th( nodes, properties ) );

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

	//***********************************************************************************************************
    //* - Projecting resulting field on plane
    //***********************************************************************************************************
    void Modeler::GenerateProjectionCoaxPortTEM(int*  NodesId, unsigned int PropertiesId)
    {
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

		Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

        if      ( mQuadraticGeometry ) pElement = Element::Pointer( new CoaxialPortTEM_2ndQ( nodes, properties ) );
        else if ( mElementOrder == 0 ) pElement = Element::Pointer( new CoaxialPortTEM_1st ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new CoaxialPortTEM_2nd ( nodes, properties ) );
        else if ( mElementOrder == 3 ) pElement = Element::Pointer( new CoaxialPortTEM_3th ( nodes, properties ) );
        else if ( mElementOrder == 4 ) pElement = Element::Pointer( new CoaxialPortTEM_4th ( nodes, properties ) );
		
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
        if ( mStaticMode == true ) return;
 
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

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PEC normal set at each node
        for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
		{
			mSetPECnormals[*id_it].push_back(normal);
		}

      //// ////////////////////////////////////////////////////////////////////
      //  if (mElementOrder == 0)
      //  {
      //      std::cout<<"FB - PEC -> zero..."<<std::endl;
      //      
      //      std::complex<double> cZero(0.00,0.00);
      //      
      //      pNode = mpModel->GetNode(NodesId[3]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //      
      //      pNode = mpModel->GetNode(NodesId[4]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //      
      //      pNode = mpModel->GetNode(NodesId[5]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //  }
      ////////////////////////////////////////////////////////////
	}

    //*************************************************************************************************
    //* - Calculates average normal in Plasma Sheath Boundary Condition (PSBC) nodes
    //*************************************************************************************************
    void Modeler::GeneratePSBCElement( std::vector<int>& LONodesId )
    {
        if ( mPlasmaMode == false ) return;

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

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PSBC normal set at each node
        for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
		{
			mSetPSBCnormals[*id_it].push_back(normal);
		}
    }

    //*************************************************************************************************
    //* - Calculates normal in PMC nodes (symmetry planes)
    //*************************************************************************************************
    void Modeler::GeneratePMCElement(std::vector<int>& LONodesId)
    {
        if ( mStaticMode == true ) return;

        if ( mPlasmaMode == true )
        {
            String sEparTol = mpColdPlasma->Get_Eparallel_Tolerance();

            if (sEparTol != "Off")
            {
                Send_Error_Msg( "PMC boundary condition incompatible with active E par tolerance",
                                "Please, eliminate PMC or change E par tolerance to Off.", 1 );
            }
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

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PMC normal set at each node
        for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
		{
			mSetPMCnormals[*id_it].push_back(normal);
		}

      //// ////////////////////////////////////////////////////////////////////
      //  if (mElementOrder == 0)
      //  {
      //      std::cout<<"FB - PMC -> zero..."<<std::endl;
      //      
      //      std::complex<double> cZero(0.00,0.00);
      //      
      //      pNode = mpModel->GetNode(NodesId[3]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //      
      //      pNode = mpModel->GetNode(NodesId[4]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //      
      //      pNode = mpModel->GetNode(NodesId[5]-1);
      //      (*mpModel)(cEx, *pNode) = cZero; pNode->pDofcEx()->FixDof();
	     //   (*mpModel)(cEy, *pNode) = cZero; pNode->pDofcEy()->FixDof();
		    //(*mpModel)(cEz, *pNode) = cZero; pNode->pDofcEz()->FixDof();
      //  }
      ////////////////////////////////////////////////////////////
    }

	//*************************************************************************************************
    //* - Calculates average normal in TE PMC nodes (symmetry plane ortogonal to a PMC plane)
    //*************************************************************************************************
    void Modeler::GenerateTEPMCElement(std::vector<int>& LONodesId)
    {
        if ( mStaticMode == true ) return;

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

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the TEPMC normal set at each node
        for(id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it)
		{
			mSetTEPMCnormals[*id_it].push_back(normal);
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
		    std::map<double, std::complex<double> >().swap(mProjectionNormalized);
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
	        PrintViParameters (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
				               mVolumeValue    , mProblemFrequency);

            WriteViHeaderFiles(mVolumeIntegralE);

            WriteViInFiles    (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
						       mProblemFrequency);

		    mVolumeIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralE);

			mVolumeIntegralH.clear();
            std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralH);

			mVolumeIntegralJ.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralJ);

			mVolumeValue.clear();
			std::map<double, double>().swap(mVolumeValue);
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters (mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                  mSurfaceValue    , mProblemFrequency);

            WriteSurfiHeaderFiles(mSurfaceIntegralE);

            WriteSurfiInFiles    (mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                  mProblemFrequency);

		    mSurfaceIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralE);

			mSurfaceIntegralH.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralH);

			mSurfaceIntegralJ.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralJ);

			mSurfaceValue.clear();
			std::map<double, double>().swap(mSurfaceValue);
        }

		//x_vector.resize(mSystemSize,0.00);
        b_vector.resize(mSystemSize,0.00);
        A_matrix.Resize(mSystemSize);

        if ( mA_matrix_aux_Required == true )
        {
            A_matrix_aux.Resize(mSystemSize);
        }
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
		    std::map<double, std::complex<double> >().swap(mProjectionNormalized);
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
            PrintViParameters(mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
				              mVolumeValue    , mProblemFrequency);

            WriteViInFiles   (mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
							  mProblemFrequency);

		    mVolumeIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralE);

			mVolumeIntegralH.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralH);

			mVolumeIntegralJ.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mVolumeIntegralJ);

			mVolumeValue.clear();
			std::map<double, double>().swap(mVolumeValue);
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters(mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                 mSurfaceValue    , mProblemFrequency);

            WriteSurfiInFiles(mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				              mProblemFrequency);

		    mSurfaceIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralE);

			mSurfaceIntegralH.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralH);

			mSurfaceIntegralJ.clear();
			std::map<double, Vector<std::complex<double> > >().swap(mSurfaceIntegralJ);

			mSurfaceValue.clear();
			std::map<double, double>().swap(mSurfaceValue);
        }

		//x_vector.resize(mSystemSize,0.00);
        b_vector.resize(mSystemSize,0.00);
        A_matrix.Resize(mSystemSize);

        if ( mA_matrix_aux_Required == true )
        {
            A_matrix_aux.Resize(mSystemSize);
        }
	}

    //*******************************************************************************
    //* - Solving linear problem
    //*******************************************************************************
    void Modeler::SolveProblem()
    {
        // Solving linear system Ax=b
        SolveLinearSystem();

		// Setting to cero the fundamental variables
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nit;
        std::complex<double> cZero( 0.00 , 0.00 );

        for (nit = nodes.begin(); nit != nodes.end(); ++nit)
		{
			(*mpModel)(cEx, **nit) = cZero;
			(*mpModel)(cEy, **nit) = cZero;
			(*mpModel)(cEz, **nit) = cZero;
		}

        // Storing x_vector solution on cExyz variables
        SaveSystemSolution(cEx);
        SaveSystemSolution(cEy);
        SaveSystemSolution(cEz);

        // Updating fixed DOFs
        if ( mPlasmaMode   ) UpdateCoordPlasmaRLP ();
        if ( mAxisymmetric ) UpdateAxisToCartesian();

		UpdateCoordEPEC   ();
		UpdateContactCoord();
		UpdateCoordEPBC   ();
    }

    //*******************************************************************************
    //* - Solving static linear problem
    //*******************************************************************************
    void Modeler::SolveStaticProblem()
    {
        // Solving linear system Ax=b
        SolveLinearSystem();

        // Saving solution vector on VOLTAGE variable
        SaveSystemSolution_Static(VOLTAGE);

        // Updating fixed DOFs
        UpdateFixedVOLTAGE();
    }

    //*******************************************************************************
    //* - Updating fixed VOLTAGE DOFs
    //*******************************************************************************
    void Modeler::UpdateFixedVOLTAGE()
    {
        // Updating fixed DOFs
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            if ( (*it)->pDofV()->IsFixed() )
            {
                (*mpModel)(VOLTAGE, **it) = mFixVoltage[(*it)->Id()];
            }
        }
    }

    //*****************************************************************************
    //* - Updates coordinate system in plasma mode
    //*****************************************************************************
    void Modeler::UpdateCoordPlasmaRLP()
    {
        // Checking if E parallel tolerance is activated
        if ( mpColdPlasma->Get_Eparallel_Tolerance() == "Off" ) return;

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;

        for ( nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it )
        {
            unsigned int nodeId = (*nodes_it)->Id();

            if ( mNormals.find(nodeId) == mNormals.end() )
            {
                // Complex constants
                std::complex<double> cZero    ( 0.0        , 0.0         );
                std::complex<double> cUnit    ( 0.0        , 1.0         );
                std::complex<double> rUnit    ( 1.0        , 0.0         );
                std::complex<double> rInvSqrt2( 1.0/sqrt(2), 0.0         );
                std::complex<double> cInvSqrt2( 0.0        , 1.0/sqrt(2) );

                // E_RLP = [U] * E_SDP
                Matrix<std::complex<double> > U_Herm(3,3);

                U_Herm[0][0] =  rInvSqrt2; U_Herm[0][1] = rInvSqrt2; U_Herm[0][2] = cZero;
                U_Herm[1][0] = -cInvSqrt2; U_Herm[1][1] = cInvSqrt2; U_Herm[1][2] = cZero;
                U_Herm[2][0] =  cZero    ; U_Herm[2][1] = cZero    ; U_Herm[2][2] = rUnit;

                // E_SDP = [T] * E_CC
                Matrix<std::complex<double> > T(3,3);

                Rotation_Matrix_SDP_T_CC( T, nodeId );

                // [R] = [Trans(T)][Herm(U)]
                Matrix<std::complex<double> > R(3,3);

                R = T.Transpose() * U_Herm;

                std::complex<double> Er = (*mpModel)(cEx, **nodes_it);
                std::complex<double> El = (*mpModel)(cEy, **nodes_it);
                std::complex<double> Ep = (*mpModel)(cEz, **nodes_it);

                (*mpModel)(cEx, **nodes_it) = R[0][0]*Er + R[0][1]*El + R[0][2]*Ep;
                (*mpModel)(cEy, **nodes_it) = R[1][0]*Er + R[1][1]*El + R[1][2]*Ep;
                (*mpModel)(cEz, **nodes_it) = R[2][0]*Er + R[2][1]*El + R[2][2]*Ep;
            }
        }
    }

	//*****************************************************************************
    //* - Updates coordinate system in axisymmetric problems
    //*****************************************************************************
    void Modeler::UpdateAxisToCartesian()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;
		unsigned int nodeId;

        for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
        {
			nodeId = (*nodes_it)->Id();

            if ( ( mNormals.find(nodeId) == mNormals.end() )  && ( DistanceToAxis(nodeId) > mGiDTolerance ) )
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
            int NodeId = (*nodes_it)->Id();

            if ( mNormals.find( NodeId ) != mNormals.end() )
            {
                // Transform from D_ntb to E_ntb in plasma sheath BC: {E_ntb} = [R * invEp * Rt] * {D_ntb}
                if ( ( mPlasmaMode == true ) && ( mType_Of_BC_Normal[NodeId] == 'D' ) )
                {
                    Matrix< std::complex<double> > RiEpRt;

                    Calculate_En_R_invEp_Rt_Dn( RiEpRt, NodeId );

                    std::complex<double> Dn = (*mpModel)(cEx, **nodes_it),
                                         Dt = (*mpModel)(cEy, **nodes_it),
                                         Db = (*mpModel)(cEz, **nodes_it);

                    (*mpModel)(cEx, **nodes_it) = RiEpRt[0][0]*Dn + RiEpRt[0][1]*Dt + RiEpRt[0][2]*Db;
                    (*mpModel)(cEy, **nodes_it) = RiEpRt[1][0]*Dn + RiEpRt[1][1]*Dt + RiEpRt[1][2]*Db;
                    (*mpModel)(cEz, **nodes_it) = RiEpRt[2][0]*Dn + RiEpRt[2][1]*Dt + RiEpRt[2][2]*Db;
                }

                // Transform from E_ntb to E_xyz in PEC/PMC BC: {E_xyz} = [Rt] * {E_ntb}
                std::vector<double> n(mNormals[NodeId]);
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
        //std::map<unsigned int, Vector<double> >().swap(mNormals);
    }

	//*****************************************************************************
    //* - Update electric field values in PBC nodes
    //*****************************************************************************
    void Modeler::UpdateCoordEPBC()
	{
		std::vector<int> FrontLeftNodesId;

		UpdateCoord_RL_PBC(FrontLeftNodesId);
		UpdateCoord_FB_PBC();
		UpdateCoord_LF_PBC(FrontLeftNodesId);
	}

	//*****************************************************************************
    //* - Check if the left node is also a front node
    //*****************************************************************************
	bool Modeler::PBC_LeftNodeIsFrontNode( int rightNodeId,
                                           std::vector<int>& leftNodesId,
                                           std::vector<int>& FrontLeftNodesId )
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
    void Modeler::UpdateCoord_RL_PBC(std::vector<int>& FrontLeftNodesId)
	{
		if ( ( mPBC_NodeNodePairs_Right.size   () +
               mPBC_NodeEdgePairs_Right.size   () +
               mPBC_NodeElementPairs_Right.size() ) == 0 )
        {
            return;
        }

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nit;

		for ( nit = nodes.begin(); nit!= nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();

            std::complex<double> lcEx(0.00, 0.00);
			std::complex<double> lcEy(0.00, 0.00);
			std::complex<double> lcEz(0.00, 0.00);

			if      ( mPBC_NodeNodePairs_Right.find(NodeId) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNode = mPBC_NodeNodePairs_Right[NodeId];

				std::vector<int> leftNodesId(1); leftNodesId[0] = lNode;

				if ( PBC_LeftNodeIsFrontNode( NodeId, leftNodesId, FrontLeftNodesId ) == true ) continue;

		        Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

				lcEx = (*mpModel)(cEx, *plNode);
				lcEy = (*mpModel)(cEy, *plNode);
				lcEz = (*mpModel)(cEz, *plNode);
            }
			else if ( mPBC_NodeEdgePairs_Right.find(NodeId) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) continue;

                double natCoord = mPBC_NodeEdgeCoord_Right[NodeId];

		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i = 0; i < lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    lcEx += N[i] * (*mpModel)(cEx, *plNode);
				    lcEy += N[i] * (*mpModel)(cEy, *plNode);
				    lcEz += N[i] * (*mpModel)(cEz, *plNode);
                }
			}
			else if ( mPBC_NodeElementPairs_Right.find(NodeId) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<int> lNodes( mPBC_NodeElementPairs_Right[NodeId] );

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) continue;

		        std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );

                std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i = 0; i < lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    lcEx += N[i] * (*mpModel)(cEx, *plNode);
				    lcEy += N[i] * (*mpModel)(cEy, *plNode);
				    lcEz += N[i] * (*mpModel)(cEz, *plNode);
                }
			}
            else
            {
                continue;
            }

            double cos_A = 1.0;
			double sin_A = 0.0;

			if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

			(*mpModel)(cEx, **nit) =  cos_A * lcEx + sin_A * lcEy;
            (*mpModel)(cEy, **nit) = -sin_A * lcEx + cos_A * lcEy;
            (*mpModel)(cEz, **nit) =  lcEz;
		}
	}

	//*****************************************************************************
    //* - Computes inverse rotation of a point
    //*****************************************************************************
    void Modeler::PBC_InvRot( int NodeId, double& cos_A, double& sin_A )
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
		if ( ( mPBC_NodeNodePairs_Front.size   () +
               mPBC_NodeEdgePairs_Front.size   () +
               mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator nit;

		for ( nit = nodes.begin(); nit!= nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();

            std::complex<double> bcEx(0.00, 0.00);
			std::complex<double> bcEy(0.00, 0.00);
			std::complex<double> bcEz(0.00, 0.00);

			if      ( mPBC_NodeNodePairs_Front.find(NodeId) != mPBC_NodeNodePairs_Front.end() )
            {
				int bNode = mPBC_NodeNodePairs_Front[NodeId];

		        Node::Pointer pbNode = mpModel->GetNode( bNode - 1 );

				bcEx = (*mpModel)(cEx, *pbNode);
				bcEy = (*mpModel)(cEy, *pbNode);
				bcEz = (*mpModel)(cEz, *pbNode);
            }
			else if ( mPBC_NodeEdgePairs_Front.find(NodeId) != mPBC_NodeEdgePairs_Front.end() )
			{
				std::vector<int>  bNodes( mPBC_NodeEdgePairs_Front[NodeId] );

                double natCoord = mPBC_NodeEdgeCoord_Front[NodeId];

		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i = 0; i < bNodes.size(); i++ )
                {
					Node::Pointer pbNode = mpModel->GetNode( bNodes[i] - 1 );

                    bcEx += N[i] * (*mpModel)(cEx, *pbNode);
				    bcEy += N[i] * (*mpModel)(cEy, *pbNode);
				    bcEz += N[i] * (*mpModel)(cEz, *pbNode);
                }
			}
			else if ( mPBC_NodeElementPairs_Front.find(NodeId) != mPBC_NodeElementPairs_Front.end() )
			{
				std::vector<int>    bNodes  ( mPBC_NodeElementPairs_Front[NodeId] );

                std::vector<double> natCoord( mPBC_NodeElementCoord_Front[NodeId] );

		        std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i = 0; i < bNodes.size(); i++ )
                {
					Node::Pointer pbNode = mpModel->GetNode( bNodes[i] - 1 );

                    bcEx += N[i] * (*mpModel)(cEx, *pbNode);
				    bcEy += N[i] * (*mpModel)(cEy, *pbNode);
				    bcEz += N[i] * (*mpModel)(cEz, *pbNode);
                }
			}
            else
            {
                 continue;
            }

            (*mpModel)(cEx, **nit) = bcEx;
            (*mpModel)(cEy, **nit) = bcEy;
            (*mpModel)(cEz, **nit) = bcEz;
        }
	}

	//*****************************************************************************
    //* - Update electric field values in PBC nodes (Left-Front nodes)
    //*****************************************************************************
    void Modeler::UpdateCoord_LF_PBC(std::vector<int>& FrontLeftNodes)
	{
		if ( FrontLeftNodes.size() == 0 ) { return; }

		std::vector<int>::iterator it_fln;

		for ( it_fln = FrontLeftNodes.begin(); it_fln != FrontLeftNodes.end(); ++it_fln )
        {
			int NodeId = *it_fln;

            std::complex<double> lcEx(0.00, 0.00);
			std::complex<double> lcEy(0.00, 0.00);
			std::complex<double> lcEz(0.00, 0.00);

			if      ( mPBC_NodeNodePairs_Right.find(NodeId) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNode = mPBC_NodeNodePairs_Right[NodeId];

				Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

				lcEx = (*mpModel)(cEx, *plNode);
				lcEy = (*mpModel)(cEy, *plNode);
				lcEz = (*mpModel)(cEz, *plNode);
            }
			else if ( mPBC_NodeEdgePairs_Right.find(NodeId) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<int>  lNodes( mPBC_NodeEdgePairs_Right[NodeId] );

                double natCoord = mPBC_NodeEdgeCoord_Right[NodeId];

		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i = 0; i < lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    lcEx += N[i] * (*mpModel)(cEx, *plNode);
				    lcEy += N[i] * (*mpModel)(cEy, *plNode);
				    lcEz += N[i] * (*mpModel)(cEz, *plNode);
                }
			}
			else if ( mPBC_NodeElementPairs_Right.find(NodeId) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<int>    lNodes  ( mPBC_NodeElementPairs_Right[NodeId] );

		        std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );

                std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i = 0; i < lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    lcEx += N[i] * (*mpModel)(cEx, *plNode);
				    lcEy += N[i] * (*mpModel)(cEy, *plNode);
				    lcEz += N[i] * (*mpModel)(cEz, *plNode);
                }
     		}
            else
            {
                continue;
            }

			double cos_A = 1.0;
			double sin_A = 0.0;

			if ( mIsRightPBCTilted == true ) PBC_InvRot( NodeId, cos_A, sin_A );

			Node::Pointer prNode = mpModel->GetNode( NodeId - 1 );

			(*mpModel)(cEx, *prNode) =  cos_A * lcEx + sin_A * lcEy;
            (*mpModel)(cEy, *prNode) = -sin_A * lcEx + cos_A * lcEy;
            (*mpModel)(cEz, *prNode) =  lcEz;
		}
	}

	//*****************************************************************************
    //* - Updating contact coordinates
    //*****************************************************************************
    void Modeler::UpdateContactCoord()
    {
 /*
  //      if ( mContactPairs.size() == 0 ) return;

		//std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

  //      std::vector<Node::Pointer>::iterator nodes_it;

  //      for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
  //      {
		//	if ( mContactPairs.find( (*nodes_it)->Id() ) != mContactPairs.end() )
  //          {
		//		ContactPairData& pContactPair = mContactPairs[ (*nodes_it)->Id() ];

  //              Node::Pointer pRNode = mpModel->GetNode( (*nodes_it)->Id()    - 1 );
		//		Node::Pointer pLNode = mpModel->GetNode( pContactPair.ctcNode - 1 );
  //              
  //              // Transform from E_ntb to E_xyz in PEC/PMC BC: {E_xyz} = [Rt] * {E_ntb}
  //              std::vector<double> n( pContactPair.ctcNormal );
		//		std::vector<double> t(3), b(3);
  //              TangencialCoordinates(n,t,b);

  //              // Material properties matrix
		//		Properties::Pointer Properties;
		//		double eo = 8.8541878176e-12;
		//		double freq = mProblemFrequency;
  //              double sigma, eps_real, eps_imag;

		//		// Material R
		//		Properties = mpModel->GetProperties(pContactPair.matR);

  //              sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
		//		eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
		//	    eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;

		//        std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );
  //              std::complex<double> cEpsL(       eo,                     0.0 );

  //              std::complex<double> invEpR = ( 1.00 / cEpsR );
  //              std::complex<double> invEpL = ( 1.00 / cEpsL );
  //             
  //              std::complex<double> EnL = invEpL * (*mpModel)(cEx, *pLNode);
  //              std::complex<double> EtL = (*mpModel)(cEy, *pLNode);
  //              std::complex<double> EbL = (*mpModel)(cEz, *pLNode);

  //              std::complex<double> EnR = invEpR * (*mpModel)(cEx, *pLNode);
  //              std::complex<double> EtR = EtL;
  //              std::complex<double> EbR = EbL;

  //              (*mpModel)(cEx, *pLNode) = n[0]*EnL + t[0]*EtL + b[0]*EbL;
  //              (*mpModel)(cEy, *pLNode) = n[1]*EnL + t[1]*EtL + b[1]*EbL;
  //              (*mpModel)(cEz, *pLNode) = n[2]*EnL + t[2]*EtL + b[2]*EbL;

  //              (*mpModel)(cEx, *pRNode) = n[0]*EnR + t[0]*EtR + b[0]*EbR;
  //              (*mpModel)(cEy, *pRNode) = n[1]*EnR + t[1]*EtR + b[1]*EbR;
  //              (*mpModel)(cEz, *pRNode) = n[2]*EnR + t[2]*EtR + b[2]*EbR;
  //          }
  //      }
*/
//        return;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
        if ( mContactPairs.size() == 0 ) return;

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator nodes_it;

        for (nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it)
        {
			if ( mContactPairs.find( (*nodes_it)->Id() ) != mContactPairs.end() )
            {
				ContactPairData& pContactPair = mContactPairs[ (*nodes_it)->Id() ];

				Node::Pointer pLNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

                // Contact matrix
                Matrix< std::complex<double> > ContactMatrix( 3, 3 );

                if ( mPlasmaMode )
                {
                    Get_Contact_Matrix_Plasma( ContactMatrix, (*nodes_it)->Id() );
                }
                else
                {
                    // Properties
                    Properties::Pointer Properties;

                    double eo   = 8.8541878176e-12;
                    double freq = mProblemFrequency;

                    double sigma;
                    double eps_real, eps_imag;

                    // Material R
                    Properties = mpModel->GetProperties( pContactPair.matR );
                    sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				    eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			        eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
                    
                    std::complex<double> cEpsR( eps_real, eps_imag + (sigma/freq) );
                    
                    // Material L
                    Properties = mpModel->GetProperties( pContactPair.matL );
                    sigma    = (*Properties)(IHL_ELECTRIC_CONDUCTIVITY);
				    eps_real = (*Properties)(REAL_ELECTRIC_PERMITTIVITY) * eo;
			        eps_imag = (*Properties)(IMAG_ELECTRIC_PERMITTIVITY) * eo;
                    
                    std::complex<double> cEpsL( eps_real, eps_imag + (sigma/freq) );

                    // Contact matrix materials properties
                    std::complex<double> rmat = (cEpsL/cEpsR) - 1.00;

                    Vector<double> n = pContactPair.ctcNormal;

                    ContactMatrix[0][0] = (n[0]*n[0]) * rmat + 1.00; ContactMatrix[0][1] = (n[0]*n[1]) * rmat       ; ContactMatrix[0][2] = (n[0]*n[2]) * rmat;
                    ContactMatrix[1][0] = (n[1]*n[0]) * rmat       ; ContactMatrix[1][1] = (n[1]*n[1]) * rmat + 1.00; ContactMatrix[1][2] = (n[1]*n[2]) * rmat;
                    ContactMatrix[2][0] = (n[2]*n[0]) * rmat       ; ContactMatrix[2][1] = (n[2]*n[1]) * rmat       ; ContactMatrix[2][2] = (n[2]*n[2]) * rmat + 1.00;
                }

                // Updating Right nodes
                std::complex<double> ExL = (*mpModel)(cEx, *pLNode),
                                     EyL = (*mpModel)(cEy, *pLNode),
                                     EzL = (*mpModel)(cEz, *pLNode);

                (*mpModel)(cEx, **nodes_it) = ContactMatrix[0][0]*ExL + ContactMatrix[0][1]*EyL + ContactMatrix[0][2]*EzL;
                (*mpModel)(cEy, **nodes_it) = ContactMatrix[1][0]*ExL + ContactMatrix[1][1]*EyL + ContactMatrix[1][2]*EzL;
				(*mpModel)(cEz, **nodes_it) = ContactMatrix[2][0]*ExL + ContactMatrix[2][1]*EyL + ContactMatrix[2][2]*EzL;
            }
        } 
    }

	//******************************************************************************************************
    //* - Returns the area weighted normal of a set of normals.
    //* - The module of the normals must be proportional to the element area.
    //******************************************************************************************************
    Vector<double> Modeler::Calculate_Area_Weighted_Normal( std::vector< std::vector<double> >& vectorSet )
    {
        // Area weighted normal
        Vector<double> aw_normal(3, 0.00);

        // Iterator over the set of normals
        std::vector< std::vector<double> >::iterator it_vs;

        // Adding all the normals in a node
        for( it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs )
        {
            Vector<double> nInSet(*it_vs);
            aw_normal += nInSet;
        }

        // Normalizing the total normal
        aw_normal *= (1.00/sqrt(std::inner_product(aw_normal.begin(),aw_normal.end(),aw_normal.begin(),0.00)));

        return aw_normal;
    }

    //******************************************************************************************************
    //* - Returns the geometric averaged normal of a set of normals.
    //******************************************************************************************************
    Vector<double> Modeler::Calculate_Geom_Averaged_Normal( std::vector< std::vector<double> >& vectorSet )
    {
        // Geometric average of a set of vectors
        Vector<double> ga_normal(3, 0.00);

        // Number of vectors in the set
        int setSize = vectorSet.size();

        // Group of vectors with a separation angle between them lower than 45 degrees
        Vector< Vector<double> > vector_group;

        // True is the vector has been assigned to a group
        Vector<bool> used_vector(setSize, 0);

        // Vectors to compare
        Vector<double> v1, v2;

        // Dot product of v1 and v2
        double v1_dot_v2;

        // Iteration over vector set
        for ( int ivs = 0; ivs<setSize; ivs++ )
        {
            // If the vector is assigned -> next vector
            if ( used_vector[ivs] == true ) continue;

            // Take vector from set and normalize
            v1 = vectorSet[ivs];
            v1 *= ( 1.00 / sqrt( std::inner_product(v1.begin(),v1.end(),v1.begin(),0.00) ) );

            // Initiate group of vectors
            Vector<double> v_group(v1);

            // Look for the vectors with a separation angle with v1 lower than 45 degrees
            for ( int jvs = ivs+1; jvs<setSize; jvs++ )
            {
                // If the vector has been assigned -> next vector
                if (used_vector[jvs] == true) continue;

                // Take vector from set and normalize
                v2 = vectorSet[jvs];
                v2 *= ( 1.00 / sqrt( std::inner_product(v2.begin(),v2.end(),v2.begin(),0.00) ) );

                // Dot product v1*v2
                v1_dot_v2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

                // If v1 and v2 angle separation is lower than 45 -> add to n1 group
                if (v1_dot_v2 > 0.5)
                {
                    v_group += v2;
                    used_vector[jvs] = true;
                }
            }

            // Normalize the group resultant
            v_group *= ( 1.00 / sqrt( std::inner_product(v_group.begin(), v_group.end(), v_group.begin(), 0.00) ) );

            // Add group resultant to the geometric average normal
            ga_normal += v_group;
        }

        // Normalize geometric average normal
        ga_normal *= (1.00/sqrt(std::inner_product(ga_normal.begin(),ga_normal.end(),ga_normal.begin(),0.00)));

        return ga_normal;
    }

	//*******************************************************************************
    //* - Setting normals in a PEC node ( nxE = 0 )
    //*******************************************************************************
	void Modeler::SettingPECNormals()
    {
		std::map< unsigned int, std::vector< std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetPECnormals.begin(); it_Set != mSetPECnormals.end(); ++it_Set)
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet(it_Set->second);

            // Resultant normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd == true ) normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            else                           normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );

			// Including the calculated normal in the mNormal vector
            mNormals[it_Set->first] = normal_atNode;

            // Setting degrees of freedom for a PEC node ( nxE = 0 )
			Node::Pointer pNode = mpModel->GetNode((it_Set->first)-1);

            // Setting DOFs in a PEC node
            std::complex<double> cZero(0.00,0.00);

			(*mpModel)(cEy, *pNode) = cZero;
            pNode->pDofcEy()->FixDof();

			(*mpModel)(cEz, *pNode) = cZero;
            pNode->pDofcEz()->FixDof();

            mType_Of_BC_Normal[pNode->Id()] = 'E';
        }

		mSetPECnormals.clear();
		std::map<unsigned int, std::vector<std::vector<double> > >().swap(mSetPECnormals);
	}

    //*******************************************************************************
    //* - Setting normals in a PSBC node ( n*D = 0 )
    //*******************************************************************************
	void Modeler::SettingPSBCNormals()
    {
		std::map< unsigned int, std::vector< std::vector<double> > >::iterator it_Set;

        for(it_Set = mSetPSBCnormals.begin(); it_Set != mSetPSBCnormals.end(); ++it_Set)
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet(it_Set->second);

            // Resultant normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd == true ) normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            else                           normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );

            // If the node is a PEC -> calculate the average normal
			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vPEC ( mNormals[it_Set->first] );
				Vector<double> vPSBC( normal_atNode           );

				vPSBC += vPEC;
                vPSBC *= (1.00/sqrt(std::inner_product(vPSBC.begin(),vPSBC.end(),vPSBC.begin(),0.00)));

				mNormals[it_Set->first] = vPSBC;
			}
            // if the node is not PEC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[it_Set->first] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first)-1 );

                std::complex<double> cZero(0.00,0.00);

				(*mpModel)(cEx, *pNode) = cZero;
				pNode->pDofcEx()->FixDof();

                mType_Of_BC_Normal[pNode->Id()] = 'D';
			}
        }

		mSetPSBCnormals.clear();
		std::map<unsigned int, std::vector<std::vector<double> > >().swap(mSetPSBCnormals);
	}

	//*******************************************************************************
    //* - Setting normals in a PMC node ( n*E = 0 )
    //*******************************************************************************
	void Modeler::SettingPMCNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

		for(it_Set = mSetPMCnormals.begin(); it_Set != mSetPMCnormals.end(); ++it_Set)
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC -> make the normals ortogonal
			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[it_Set->first] );
				Vector<double> vOrto(3, 0.00);

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[it_Set->first] = vOrto;
			}
            // if the node is not PEC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[it_Set->first] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first)-1 );

                std::complex<double> cZero(0.00,0.00);

				(*mpModel)(cEx, *pNode) = cZero;
				pNode->pDofcEx()->FixDof();

                mType_Of_BC_Normal[pNode->Id()] = 'H';
			}
        }

        mSetPMCnormals.clear();
        std::map<unsigned int, std::vector<std::vector<double> > >().swap(mSetPMCnormals);
	}

	//*******************************************************************************
    //* - Setting normals in a TE port ( ortogonal to a PMC surface )
    //*******************************************************************************
	void Modeler::SettingTEPMCNormals()
    {
		std::map<unsigned int,  std::vector<std::vector<double> > >::iterator it_Set;

		for(it_Set = mSetTEPMCnormals.begin(); it_Set != mSetTEPMCnormals.end(); ++it_Set)
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // TE_PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC or PMC -> make the normals ortogonal
			if ( mNormals.find(it_Set->first) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[it_Set->first] );
				Vector<double> vOrto(3, 0.00);

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[it_Set->first] = vOrto;
			}
            // If the node is not PEC/PMC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[it_Set->first] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first)-1 );

                std::complex<double> cZero(0.00,0.00);

				(*mpModel)(cEx, *pNode) = cZero;
				pNode->pDofcEx()->FixDof();

                mType_Of_BC_Normal[pNode->Id()] = 'T';
			}
        }

        mSetTEPMCnormals.clear();
        std::map<unsigned int, std::vector<std::vector<double> > >().swap(mSetTEPMCnormals);
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
        //We join to the Dirichlet sets the contacts nodes normals
		JoinContactNormalsInDirichletS();

		//Setting Dirichlet normals
		SettingPECNormals ();
        SettingPSBCNormals();

		//Dirichlet normals are forced in dielectric-Dirichlet intersections.
		ForceDirichletNormalsInContact();

		//Setting simmetry normals
		SettingPMCNormals  ();
		SettingTEPMCNormals();

		//Contact normals must be ortogonal to simmetry normals
		OrtogonalizeContactNormals();
    }

    //*******************************************************************************
    //* - Initialize Building process
    //*******************************************************************************
    void Modeler::IniBuilding()
    {
        std::cout << "Nodes read."    << std::endl << std::endl;

        std::cout << "Setting normals and DOFs..." << std::endl;

		SetContactDiscontinuity();
		SetNormals();
		SetPBC    ();

		if (mSetAllEzToZero ) SetAllcEzToCero ();
		if (mSetAllExyToZero) SetAllcExyToCero();
		if (mAxisymmetric   ) SetAllcEyzToCero();

        SetDofSet();

		std::cout << "Done." << std::endl << std::endl;
    }

	//*******************************************************************************
    //* - Set all cEz to cero (3D-Exy problems)
    //*******************************************************************************
    void Modeler::SetAllcEzToCero()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator it;
		std::complex<double> cZero(0.00,0.00);

		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			Node::Pointer pNode = *it;

			(*mpModel)(cEz, *pNode) = cZero;
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
		std::complex<double> cZero(0.00,0.00);

		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			Node::Pointer pNode = *it;

			(*mpModel)(cEx, *pNode) = cZero;
			pNode->pDofcEx()->FixDof();

			(*mpModel)(cEy, *pNode) = cZero;
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
		std::complex<double> cZero(0.00,0.00);

		Node::Pointer pNode;
		unsigned int nodeId;

		for (it = nodes.begin(); it != nodes.end(); ++it)
		{
			pNode  = *it;
			nodeId = pNode->Id();

			(*mpModel)(cEy, *pNode) = cZero;
			pNode->pDofcEy()->FixDof();

			(*mpModel)(cEz, *pNode) = cZero;
			pNode->pDofcEz()->FixDof();

			if ( ( mNormals.find(nodeId)  != mNormals.end() ) ||
				 ( DistanceToAxis(nodeId) <= mGiDTolerance  ) )
            {
				(*mpModel)(cEx, *pNode) = cZero;
				pNode->pDofcEx()->FixDof();
			}
		}
	}

    //*******************************************************************************
    //* - Show global matrix
    //*******************************************************************************
    void Modeler::Show_Global_Matrix()
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        for (int i = 0; i < A_matrix.RowsNumber(); i++)
        {
            data_iterator     = A_matrix[i].begin();
            end_data_iterator = A_matrix[i].end();

            while(data_iterator != end_data_iterator)
            {
                int         i_col = data_iterator->first;
                ComplexType aij   = data_iterator->second;

                std::cout << i << " , " << i_col << " : " << aij << std::endl;

                data_iterator++;
            }
        }

        if ( mA_matrix_aux_Required == false ) return;

        for (int i = 0; i < A_matrix_aux.RowsNumber(); i++)
        {
            data_iterator     = A_matrix_aux[i].begin();
            end_data_iterator = A_matrix_aux[i].end();

            while(data_iterator != end_data_iterator)
            {
                int         i_col = data_iterator->first;
                ComplexType aij   = data_iterator->second;

                std::cout << i << " , " << i_col << " (aux): " << aij << std::endl;

                data_iterator++;
            }
        }
    }

    //*******************************************************************************
    //* - Solve linear system
    //*******************************************************************************
    void Modeler::SolveLinearSystem()
    {
        if (mReadSolutionMode)
		{
			ReadSolutionVector(x_vector);
			A_matrix.FreeData();
			b_vector.FreeData();
            A_matrix_aux.FreeData();
            return;
		}

		if (mSolveWithExternal)
		{
			WriteMatrixInFile(A_matrix);
			WriteVectorInFile(b_vector);

			A_matrix.FreeData();
			b_vector.FreeData();

            if ( mA_matrix_aux_Required == true )
            {
                WriteAuxMatrixInFile(A_matrix_aux);
            }

            A_matrix_aux.FreeData();

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

            if ( mA_matrix_aux_Required == false )
            {
                mpLinearSolver->Solve( A_matrix, b_vector, x_vector );
            }
            else if ( mpColdPlasma->Is_HermSymm_Matrix() )
            {
                mpLinearSolver->Solve_HermSymm( A_matrix, A_matrix_aux, b_vector, x_vector );
            }
            else
            {
                mpLinearSolver->Solve_HermFull( A_matrix, A_matrix_aux, b_vector, x_vector );
            }

            if (mWriteResultFinalStep) WriteResultsInFile(x_vector);

            std::cout << "Solver finished." << std::endl << std::endl;
		}

        A_matrix.FreeData();
        b_vector.FreeData();
        A_matrix_aux.FreeData();
    }

    //*************************************************
    //* - Set singular nodes
    //*************************************************
    void Modeler::FixDof(unsigned int NodeId,
		                 const Variable<double>& rVariable,
		                 const double& Value)
    {
		int varKey = rVariable.getKey();

		if ( varKey == Sg.getKey() )
        {
            mSingular[NodeId] = Value;
        }
        else if ( varKey == VOLTAGE.getKey() )
        {
            Node::Pointer pnode = mpModel->GetNode(NodeId-1);

            pnode->pDofV()->FixDof();

            mFixVoltage[NodeId] = Value;
        }
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
    //* - Set Linear solver in plasma mode
    //*******************************************************************************
    void Modeler::GenerateLinearSolver_Plasma( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
    {
        mSolveWithExternal  = false;
        mExternalSolverPath = ""   ;

        int NumTreads = mItSolverNumThreads;

        if      ( SolverType == "Bi_Conjugate_Gradient"  )
        {
            mpLinearSolver = ComplexSolver::Pointer(new BiCG_CP_ComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
        }
        else if ( SolverType == "Quasi_Minimal_Residual" )
        {
            mpLinearSolver = ComplexSolver::Pointer(new QMR_CP_ComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
        }
        else if ( SolverType == "External_solver"        )
        {
            mSolveWithExternal = true;
            ReadExternalPathFromFile(mBaseFileName, mExternalSolverPath);
        }
        else
        {
            mpLinearSolver = ComplexSolver::Pointer(new BiCG_CP_ComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
        }
    }

    //*******************************************************************************
    //* - Set linear solver
    //*******************************************************************************
    void Modeler::GenerateLinearSolver( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
    {
        mSolveWithExternal  = false;
		mExternalSolverPath = ""   ;

		int NumTreads = mItSolverNumThreads;

		if      ( SolverType == "Bi_Conjugate_Gradient" )
		{
			mpLinearSolver = ComplexSolver::Pointer(new BiCGComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}
		else if ( SolverType == "Quasi_Minimal_Residual" )
		{
			mpLinearSolver = ComplexSolver::Pointer(new QMRComplexSolver(NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name));
		}
		else if ( SolverType == "External_solver"        )
		{
			mSolveWithExternal = true;
			ReadExternalPathFromFile(mBaseFileName, mExternalSolverPath);
		}
		else
		{
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

        x_vector.resize(mSystemSize, 0.00);
        b_vector.resize(mSystemSize, 0.00);
        A_matrix.Resize(mSystemSize);

        mA_matrix_aux_Required = false;

        if ( mPlasmaMode == true )
        {
            if ( !mpColdPlasma->Is_Full_Matrix() )
            {
                mA_matrix_aux_Required = true;
                A_matrix_aux.Resize(mSystemSize);
            }
        }
    }

	//*******************************************************
    //* - Print projections
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

    //*************************************************************************
    //* - Print the results of the field integrals for the static mode.
    //*************************************************************************
    void Modeler::PrintProjection_Static()
    {
        if ( mVolumeValue.size() > 0 )
        {
            PrintViParameters_Static( mVolIntegralE_Static   ,
                                      mVolIntegralJ_Static   ,
                                      mVolIntegralE_Static_GP,
                                      mVolumeValue           );
        }

        if ( mSurfaceValue.size() > 0 )
        {
            PrintSurfiParameters_Static( mSurfIntegralE_Static   ,
                                         mSurfIntegralJ_Static   ,
                                         mSurfIntegralE_Static_GP,
                                         mSurfaceValue           );
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
            mStaticMode      = false;
		}
		else if (ProblemType == "Ez3D")
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = true;
			mAxisymmetric    = false;
			mPlasmaMode      = false;
            mStaticMode      = false;
		}
		else if (ProblemType == "Exy3D")
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = true;
			mSetAllExyToZero = false;
			mAxisymmetric    = false;
			mPlasmaMode      = false;
            mStaticMode      = false;
		}
		else if (ProblemType == "Ea3D")
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = false;
			mAxisymmetric    = true;
			mPlasmaMode      = false;
            mStaticMode      = false;
		}
		else if (ProblemType == "Plasma")
		{
			mProblemType     = E3D;
			mSetAllEzToZero  = false;
			mSetAllExyToZero = false;
			mAxisymmetric    = false;
			mPlasmaMode      = true;
            mStaticMode      = false;
		}
        else if (ProblemType == "Static")
        {
            mProblemType     = E3D;
            mSetAllEzToZero  = false;
            mSetAllExyToZero = false;
            mAxisymmetric    = false;
            mPlasmaMode      = false;
            mStaticMode      = true;
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
        else if (ProblemType == "Project_Static_Fields")
        {
            PrintProjection_Static();
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
        else if (ProblemType == "IMPJON")
		{
			mImportCurrents = true;
		}
        else if (ProblemType == "IMPJOFF")
		{
			mImportCurrents = false;
		}
		else if (ProblemType == "QE")
		{
			mQuadraticGeometry = true;
		}
		else if (ProblemType == "LE")
		{
			mQuadraticGeometry = false;
		}
        else if (ProblemType == "CheckConsistency")
        {
            CheckConsistency();
        }

        else if (ProblemType == "Bubble") { mElementOrder = 0; }
        else if (ProblemType == "1st"   ) { mElementOrder = 1; }
        else if (ProblemType == "2nd"   ) { mElementOrder = 2; }
        else if (ProblemType == "3rd"   ) { mElementOrder = 3; }
        else if (ProblemType == "4th"   ) { mElementOrder = 4; }
        
        else if (ProblemType == "GiDTol3" ) mGiDTolerance = 1e-3 ;
        else if (ProblemType == "GiDTol6" ) mGiDTolerance = 1e-6 ;
        else if (ProblemType == "GiDTol9" ) mGiDTolerance = 1e-9 ;
        else if (ProblemType == "GiDTol12") mGiDTolerance = 1e-12;

		else if (ProblemType ==   "1pr") mItSolverNumThreads =   1;
        else if (ProblemType ==   "2pr") mItSolverNumThreads =   2;
        else if (ProblemType ==   "4pr") mItSolverNumThreads =   4;
        else if (ProblemType ==   "8pr") mItSolverNumThreads =   8;
		else if (ProblemType ==  "16pr") mItSolverNumThreads =  16;
        else if (ProblemType ==  "32pr") mItSolverNumThreads =  32;
        else if (ProblemType ==  "64pr") mItSolverNumThreads =  64;
        else if (ProblemType == "128pr") mItSolverNumThreads = 128;
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
    //* - Returns true if the problem mode is Static.
    //*******************************************************
    bool Modeler::Is_StaticMode()
    {
        return mStaticMode;
    }

    //*******************************************************
    //* - Returns true if the problem mode is Plasma.
    //*******************************************************
    bool Modeler::Is_PlasmaMode()
    {
        return mPlasmaMode;
    }

    //*******************************************************
    //* - Setting frequency mode
    //*******************************************************
    void Modeler::SetFrequencyMode(bool SweepingFreq)
    {
        mFrequencySweep = SweepingFreq;
    }

    //*******************************************************
    //* - Setting frequency
    //*******************************************************
    void Modeler::SetFrequency(double ProblemFrequency)
    {
        mProblemFrequency = ProblemFrequency;
    }

	//************************************************************
	//* - Creates plasma object and reads plasma data from files
	//************************************************************
	void Modeler::LoadPlasmaParameters(String plasma_file_name)
	{
		mpColdPlasma = ColdPlasma::Pointer(new ColdPlasma());
		mpColdPlasma->Load_Data(plasma_file_name);
	}

	//************************************************************
	//* - Obtains electron density at every node
	//************************************************************
	void Modeler::Extract_ColdPlasma_ElectronDensity()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			(*mpModel)(ELECTRON_DENSITY, **it) = mpColdPlasma->Get_ElectronDensity_InNode( *it );
		}
	}

	//************************************************************
	//* - Obtains applied external B field at every node
	//************************************************************
	void Modeler::Extract_ColdPlasma_Bexternal()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			(*mpModel)(B_EXT, **it) = mpColdPlasma->Get_Bexternal_InNode( *it );
		}
	}

	//*******************************************************************************************
	//* - Extracts S,D,P,R,L components of the plasma permittivity tensor at every node
	//*******************************************************************************************
	void Modeler::Extract_ColdPlasma_PermittivityTensor()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			Vector<std::complex<double> >  SDPRL = mpColdPlasma->Get_SDPRL_InNode( *it, mProblemFrequency );

			(*mpModel)(PLASMA_PERMITTIVITY_S, **it) = std::real( SDPRL[0] );
			(*mpModel)(PLASMA_PERMITTIVITY_D, **it) = std::real( SDPRL[1] );
			(*mpModel)(PLASMA_PERMITTIVITY_P, **it) = std::real( SDPRL[2] );
			(*mpModel)(PLASMA_PERMITTIVITY_R, **it) = std::real( SDPRL[3] );
			(*mpModel)(PLASMA_PERMITTIVITY_L, **it) = std::real( SDPRL[4] );
		}
	}

    //********************************************************************************
    //* - Obtains E field component parallel to B_ext
    //********************************************************************************
    void Modeler::Extract_ColdPlasma_Eparallel()
    {
        Model::NodesArrayType::iterator it;

        for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
        {
            // External magnetic flux density
            Vector<double> Bext = mpColdPlasma->Get_Bexternal_InNode(*it);

            // b = B / |B|
            double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

            Vector<double> b(3);
            b[0] = Bext[0] / Bnorm;
            b[1] = Bext[1] / Bnorm;
            b[2] = Bext[2] / Bnorm;

            // E field
            Vector<double> E_real = (*mpModel)(REAL_E, **it);
            Vector<double> E_imag = (*mpModel)(IMAG_E, **it);

            // E field proyection onto B_ext
            double E_par_real = E_real[0]*b[0] + E_real[1]*b[1] + E_real[2]*b[2];
            double E_par_imag = E_imag[0]*b[0] + E_imag[1]*b[1] + E_imag[2]*b[2];

            double mod_E_par  = std::sqrt( E_par_real*E_par_real + E_par_imag*E_par_imag );

            (*mpModel)(MOD_E_PARALLEL, **it) = mod_E_par;
        }
    }

    //********************************************************************************
    //* - Obtains E field component perpendicular to B_ext
    //********************************************************************************
    void Modeler::Extract_ColdPlasma_Eperpendicular()
    {
        Model::NodesArrayType::iterator it;

        for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
        {
            // External magnetic flux density
            Vector<double> Bext = mpColdPlasma->Get_Bexternal_InNode(*it);

            // b = B / |B|
            double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

            Vector<double> b(3);
            b[0] = Bext[0] / Bnorm;
            b[1] = Bext[1] / Bnorm;
            b[2] = Bext[2] / Bnorm;

            // E field
            Vector<double> E_real = (*mpModel)(REAL_E, **it);
            Vector<double> E_imag = (*mpModel)(IMAG_E, **it);

            // E field proyection onto B_ext
            double E_par_real = E_real[0]*b[0] + E_real[1]*b[1] + E_real[2]*b[2];
            double E_par_imag = E_imag[0]*b[0] + E_imag[1]*b[1] + E_imag[2]*b[2];

            // Real part of the perpendicular component (Eper = E - Epar)
            Vector<double> E_per_real(3);

            E_per_real[0] = E_real[0] - E_par_real*b[0];
            E_per_real[1] = E_real[1] - E_par_real*b[1];
            E_per_real[2] = E_real[2] - E_par_real*b[2];

            // Imaginary part of the perpendicular component (Eper = E - Epar)
            Vector<double> E_per_imag(3);

            E_per_imag[0] = E_imag[0] - E_par_imag*b[0];
            E_per_imag[1] = E_imag[1] - E_par_imag*b[1];
            E_per_imag[2] = E_imag[2] - E_par_imag*b[2];

            (*mpModel)(MOD_E_PERPENDICULAR, **it) = std::sqrt( E_per_real[0]*E_per_real[0] + E_per_imag[0]*E_per_imag[0] +
                                                               E_per_real[1]*E_per_real[1] + E_per_imag[1]*E_per_imag[1] +
                                                               E_per_real[2]*E_per_real[2] + E_per_imag[2]*E_per_imag[2] );
        }
    }

	//*********************************************************************
	//* - Extracts all the plasma parameters at every node
	//*********************************************************************
	void Modeler::Extract_ColdPlasma_AllParameters()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			double eDensity;
            Vector<double> Bext;
            Vector<std::complex<double> > SDPRL;

            mpColdPlasma->Get_AllPlasmaParameters_InNode( eDensity, Bext, SDPRL, *it, mProblemFrequency );

			(*mpModel)(ELECTRON_DENSITY,      **it) = eDensity;
			(*mpModel)(B_EXT,                 **it) = Bext;
            (*mpModel)(PLASMA_PERMITTIVITY_S, **it) = std::real( SDPRL[0] );
            (*mpModel)(PLASMA_PERMITTIVITY_D, **it) = std::real( SDPRL[1] );
            (*mpModel)(PLASMA_PERMITTIVITY_P, **it) = std::real( SDPRL[2] );
            (*mpModel)(PLASMA_PERMITTIVITY_R, **it) = std::real( SDPRL[3] );
            (*mpModel)(PLASMA_PERMITTIVITY_L, **it) = std::real( SDPRL[4] );
		}
	}

    //**************************************************************************************
    //* - Set to zero de parallel component of E when tolerance criteria is met
    //**************************************************************************************
    void Modeler::Set_E_parallel_zeros()
    {
        String sEparTol = mpColdPlasma->Get_Eparallel_Tolerance();

        // Option deactivated
        if (sEparTol == "Off")
        {
            return;
        }

        Model::NodesArrayType::iterator it;

        // E parallel equal to zero when P < 0 and |P| > Tol
        double Tol = atof(sEparTol);

        for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
        {
            Vector<std::complex<double> > SDPRL = mpColdPlasma->Get_SDPRL_InNode( *it, mProblemFrequency );

            double P = std::real( SDPRL[2] );

            if ( ( P < 0 ) && ( abs(P) > Tol ) )
            {
                std::complex<double> cZero(0.00,0.00);

                (*mpModel)(cEz, **it) = cZero;
                (*it)->pDofcEz()->FixDof();
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // - Makes the average of H in contact nodes
    ////////////////////////////////////////////////////////////////////////////
    void Modeler::Average_H_In_Contact_Nodes()
    {
        if ( mContactPairs.size() == 0 ) return;

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        for ( it = nodes.begin(); it != nodes.end(); ++it )
        {
            Node::Pointer pRNode = (*it);

            int RNodeId = pRNode->Id();

            if ( mContactPairs.find( RNodeId ) != mContactPairs.end() )
            {
                ContactPairData& pContactPair = mContactPairs[ RNodeId ];

                Node::Pointer pLNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

                int LNodeId = pLNode->Id();

                int Total_Neighbours = mNeigbourElements[ RNodeId-1 ] +
                                       mNeigbourElements[ LNodeId-1 ] ;

                mNeigbourElements[ RNodeId-1 ] = Total_Neighbours;
                mNeigbourElements[ LNodeId-1 ] = Total_Neighbours;

                Vector<double> RH_real = (*mpModel)(REAL_H, *pRNode);
                Vector<double> RH_imag = (*mpModel)(IMAG_H, *pRNode);

                Vector<double> LH_real = (*mpModel)(REAL_H, *pLNode);
                Vector<double> LH_imag = (*mpModel)(IMAG_H, *pLNode);

                Vector<double> Total_H_real(3,0.00);
                Total_H_real += RH_real;
                Total_H_real += LH_real;

                Vector<double> Total_H_imag(3,0.00);
                Total_H_imag += RH_imag;
                Total_H_imag += LH_imag;

                (*mpModel)(REAL_H, *pRNode) = Total_H_real;
                (*mpModel)(IMAG_H, *pRNode) = Total_H_imag;

                (*mpModel)(REAL_H, *pLNode) = Total_H_real;
                (*mpModel)(IMAG_H, *pLNode) = Total_H_imag;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // - Applies PEC and PMC conditions to H
    ////////////////////////////////////////////////////////////////////////////
    void Modeler::Apply_PEC_and_PMC_NBC_To_H()
    {
        if ( mNormals.size() == 0 ) return;

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        for ( it = nodes.begin(); it != nodes.end(); ++it )
        {
            int NodeId = (*it)->Id();

            if ( mNormals.find( NodeId ) != mNormals.end() )
            {
                std::vector<double> n( mNormals[NodeId] );
                std::vector<double> t( 3 ), b( 3 );

                TangencialCoordinates(n, t, b);

                Vector<double> free_H_real = (*mpModel)(REAL_H, **it);
                Vector<double> free_H_imag = (*mpModel)(IMAG_H, **it);

                Vector<double> fixd_H_real(3, 0.00);
                Vector<double> fixd_H_imag(3, 0.00);

                // If PEC make n*H = 0
                if ( mType_Of_BC_Normal[NodeId] == 'H' )//'E' )
                {
                    double Ht_real = t[0] * free_H_real[0] + t[1] * free_H_real[1] + t[2] * free_H_real[2];
                    double Hb_real = b[0] * free_H_real[0] + b[1] * free_H_real[1] + b[2] * free_H_real[2];

                    double Ht_imag = t[0] * free_H_imag[0] + t[1] * free_H_imag[1] + t[2] * free_H_imag[2];
                    double Hb_imag = b[0] * free_H_imag[0] + b[1] * free_H_imag[1] + b[2] * free_H_imag[2];

                    fixd_H_real[0] = t[0] * Ht_real + b[0] * Hb_real;
                    fixd_H_real[1] = t[1] * Ht_real + b[1] * Hb_real;
                    fixd_H_real[2] = t[2] * Ht_real + b[2] * Hb_real;

                    fixd_H_imag[0] = t[0] * Ht_imag + b[0] * Hb_imag;
                    fixd_H_imag[1] = t[1] * Ht_imag + b[1] * Hb_imag;
                    fixd_H_imag[2] = t[2] * Ht_imag + b[2] * Hb_imag;
                }
                // If PMC make nxH = 0
                else if ( mType_Of_BC_Normal[NodeId] == 'E' )//'H' )
                {
                    double Hn_real = n[0] * free_H_real[0] + n[1] * free_H_real[1] + n[2] * free_H_real[2];
                    double Hn_imag = n[0] * free_H_imag[0] + n[1] * free_H_imag[1] + n[2] * free_H_imag[2];

                    fixd_H_real[0] = n[0] * Hn_real;
                    fixd_H_real[1] = n[1] * Hn_real;
                    fixd_H_real[2] = n[2] * Hn_real;

                    fixd_H_imag[0] = n[0] * Hn_imag;
                    fixd_H_imag[1] = n[1] * Hn_imag;
                    fixd_H_imag[2] = n[2] * Hn_imag;
                }
                // If not PEC or PMC let H as calculated
                else
                {
                    fixd_H_real[0] = free_H_real[0];
                    fixd_H_real[1] = free_H_real[1];
                    fixd_H_real[2] = free_H_real[2];

                    fixd_H_imag[0] = free_H_imag[0];
                    fixd_H_imag[1] = free_H_imag[1];
                    fixd_H_imag[2] = free_H_imag[2];                
                }

                (*mpModel)(REAL_H, **it) = fixd_H_real;
                (*mpModel)(IMAG_H, **it) = fixd_H_imag;
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // - Calculates of nodal magnitudes
    ////////////////////////////////////////////////////////////////////////////
    void Modeler::CalculateNodalE3D(int variableKey)
    {
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        // E(t)
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
		// H(t)
        else if (variableKey == H.getKey())
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
		// J(t)
        else if (variableKey == J.getKey())
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
        // REAL_E
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
        // IMAG_E
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
        // MOD_E
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
		// LOG_MOD_E
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
		// REAL_H
		else if (variableKey == REAL_H.getKey())
        {
             int NumNeighbours;

             Vector<double> H_real;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 NumNeighbours = mNeigbourElements[(*it)->Id()-1];

                 H_real = (*mpModel)(REAL_H, **it);

                 if (NumNeighbours > 0)
				 {
				     H_real /= NumNeighbours;
				 }

                 (*mpModel)(REAL_H, **it) = H_real;
             }
        }
		// IMAG_H
        else if (variableKey == IMAG_H.getKey())
        {
             int NumNeighbours;

             Vector<double> H_imag;

             for (it=nodes.begin(); it!=nodes.end(); ++it)
             {
                 NumNeighbours = mNeigbourElements[(*it)->Id()-1];

                 H_imag = (*mpModel)(IMAG_H, **it);

				 if (NumNeighbours > 0)
				 {
					 H_imag /= NumNeighbours;
				 }

                 (*mpModel)(IMAG_H, **it) = H_imag;
             }
        }
		// MOD_H
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
		// LOG_MOD_H
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
		// BOUNDARY_NORMALS
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
        // CONTACT_NORMALS
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
        // E field parallel to B_ext
        else if (variableKey == E_PARALLEL_T.getKey())
        {
            double time = mpModel->CurrentTime();

            for (it=nodes.begin(); it!=nodes.end(); ++it)
            {
                // External magnetic flux density
                Vector<double> Bext = mpColdPlasma->Get_Bexternal_InNode(*it);

                // b = B / |B|
                double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

                Vector<double> b(3);
                b[0] = Bext[0] / Bnorm;
                b[1] = Bext[1] / Bnorm;
                b[2] = Bext[2] / Bnorm;

                // E field
                Vector<double> E_real = (*mpModel)(REAL_E, **it);
                Vector<double> E_imag = (*mpModel)(IMAG_E, **it);

                std::vector<double> E_total(3);

                E_total[0] = E_real[0]*cos(mProblemFrequency*time) + E_imag[0]*sin(mProblemFrequency*time);
                E_total[1] = E_real[1]*cos(mProblemFrequency*time) + E_imag[1]*sin(mProblemFrequency*time);
                E_total[2] = E_real[2]*cos(mProblemFrequency*time) + E_imag[2]*sin(mProblemFrequency*time);

                // E field proyection on B_ext
                double EparDotb = E_total[0]*b[0] + E_total[1]*b[1] + E_total[2]*b[2];

                std::vector<double> E_parallel_t(3);

                E_parallel_t[0] = EparDotb * b[0];
                E_parallel_t[1] = EparDotb * b[1];
                E_parallel_t[2] = EparDotb * b[2];

                (*mpModel)(E_PARALLEL_T, **it) = E_parallel_t;
            }
        }
        // E field perpendicular to B_ext
        else if (variableKey == E_PERPENDICULAR_T.getKey())
        {
            double time = mpModel->CurrentTime();

            for (it=nodes.begin(); it!=nodes.end(); ++it)
            {
                // External magnetic flux density
                Vector<double> Bext = mpColdPlasma->Get_Bexternal_InNode(*it);

                // b = B / |B|
                double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

                Vector<double> b(3);
                b[0] = Bext[0] / Bnorm;
                b[1] = Bext[1] / Bnorm;
                b[2] = Bext[2] / Bnorm;

                // E field
                Vector<double> E_real = (*mpModel)(REAL_E, **it);
                Vector<double> E_imag = (*mpModel)(IMAG_E, **it);

                std::vector<double> E_total(3);

                E_total[0] = E_real[0]*cos(mProblemFrequency*time) + E_imag[0]*sin(mProblemFrequency*time);
                E_total[1] = E_real[1]*cos(mProblemFrequency*time) + E_imag[1]*sin(mProblemFrequency*time);
                E_total[2] = E_real[2]*cos(mProblemFrequency*time) + E_imag[2]*sin(mProblemFrequency*time);

                // E field proyection on the perpendicular direction to B_ext (Eper = E - Epar)
                double EparDotb = E_total[0]*b[0] + E_total[1]*b[1] + E_total[2]*b[2];

                std::vector<double> E_perpendicular_t(3);

                E_perpendicular_t[0] = E_total[0] - EparDotb * b[0];
                E_perpendicular_t[1] = E_total[1] - EparDotb * b[1];
                E_perpendicular_t[2] = E_total[2] - EparDotb * b[2];

                (*mpModel)(E_PERPENDICULAR_T, **it) = E_perpendicular_t;
            }
        }
        // Poynting vector
        else if (variableKey == POYNTING_VECTOR.getKey())
        {
            // Imaginary unit
            std::complex<double> cUnit(0.00, 1.00);

            for (it=nodes.begin(); it!=nodes.end(); ++it)
            {
                // E field
                Vector<double> E_real = (*mpModel)(REAL_E, **it);
                Vector<double> E_imag = (*mpModel)(IMAG_E, **it);

                Vector<std::complex<double> > cEv(3);
                cEv[0] = E_real[0] + cUnit * E_imag[0];
                cEv[1] = E_real[1] + cUnit * E_imag[1];
                cEv[2] = E_real[2] + cUnit * E_imag[2];

                // H field
                Vector<double> H_real = (*mpModel)(REAL_H, **it);
                Vector<double> H_imag = (*mpModel)(IMAG_H, **it);

                Vector<std::complex<double> > cHv(3);
                cHv[0] = H_real[0] + cUnit * H_imag[0];
                cHv[1] = H_real[1] + cUnit * H_imag[1];
                cHv[2] = H_real[2] + cUnit * H_imag[2];

                //  E X conj(H)
                Vector<std::complex<double> > E_X_conj_H(3);
                E_X_conj_H[0] = cEv[1] * std::conj(cHv[2]) - cEv[2] * std::conj(cHv[1]);
                E_X_conj_H[1] = cEv[2] * std::conj(cHv[0]) - cEv[0] * std::conj(cHv[2]);
                E_X_conj_H[2] = cEv[0] * std::conj(cHv[1]) - cEv[1] * std::conj(cHv[0]);

                // < S > = 0.5 * real( E X conj(H) )
                Vector<double> avgS(3);
                avgS[0] = 0.5 * std::real( E_X_conj_H[0] );
                avgS[1] = 0.5 * std::real( E_X_conj_H[1] );
                avgS[2] = 0.5 * std::real( E_X_conj_H[2] );

                (*mpModel)(POYNTING_VECTOR, **it) = avgS;
            }
        }
        // Lorentz force
        else if (variableKey == LORENTZ_FORCE.getKey())
        {
            // Imaginary unit
            std::complex<double> cUnit(0.00, 1.00);

            for (it=nodes.begin(); it!=nodes.end(); ++it)
            {
                // J total
                Vector<double> J_real = (*mpModel)(REAL_J, **it);
                Vector<double> J_imag = (*mpModel)(IMAG_J, **it);

                Vector< std::complex<double> > cJv(3);
                cJv[0] = J_real[0] + cUnit * J_imag[0];
                cJv[1] = J_real[1] + cUnit * J_imag[1];
                cJv[2] = J_real[2] + cUnit * J_imag[2];

                // H field
                Vector<double> H_real = (*mpModel)(REAL_H, **it);
                Vector<double> H_imag = (*mpModel)(IMAG_H, **it);

		        double mu0 = 4.00e-7 * 3.1415926535897932384626433832795;

                Vector< std::complex<double> > cBv(3);
                cBv[0] = mu0 * ( H_real[0] + cUnit * H_imag[0] );
                cBv[1] = mu0 * ( H_real[1] + cUnit * H_imag[1] );
                cBv[2] = mu0 * ( H_real[2] + cUnit * H_imag[2] );

                //  J X conj(B)
                Vector< std::complex<double> > J_X_conj_B(3);
                J_X_conj_B[0] = cJv[1] * std::conj(cBv[2]) - cJv[2] * std::conj(cBv[1]);
                J_X_conj_B[1] = cJv[2] * std::conj(cBv[0]) - cJv[0] * std::conj(cBv[2]);
                J_X_conj_B[2] = cJv[0] * std::conj(cBv[1]) - cJv[1] * std::conj(cBv[0]);

                // < F > = 0.5 * real( J X conj(B) )
                Vector<double> avgF(3);
                avgF[0] = 0.5 * std::real( J_X_conj_B[0] );
                avgF[1] = 0.5 * std::real( J_X_conj_B[1] );
                avgF[2] = 0.5 * std::real( J_X_conj_B[2] );

                (*mpModel)(LORENTZ_FORCE, **it) = avgF;
            }
        }
    }

    //*****************************************************************
    //* - Initial set up to calculate statict E field.
    //*****************************************************************
    void Modeler::IniStatic_E_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        int size = nodes.size();

        mNeigbourElements.resize(size);

        for (int i=0; i<size; i++) mNeigbourElements[i] = 0;

        std::vector<Node::Pointer>::iterator it;

        Vector<double> vCero(3, 0.00);

        for (it=nodes.begin(); it!=nodes.end(); ++it)
        {
            (*mpModel)(ELECTRIC_FIELD_SMOOTHED, **it) = vCero;
        }
    }

    //*****************************************************************
    //* - Ending calculation of statict E field.
    //*****************************************************************
    void Modeler::EndStatic_E_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        Vector<double> static_E;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeigbourElements[(*it)->Id()-1];

            static_E = (*mpModel)(ELECTRIC_FIELD_SMOOTHED, **it);

            if (NumNeighbours > 0)
            {
                static_E /= NumNeighbours;
            }

            (*mpModel)(ELECTRIC_FIELD_SMOOTHED, **it) = static_E;
        }
    }

    //*****************************************************************
    //* - Initial set up to calculate static J current density.
    //*****************************************************************
    void Modeler::IniStatic_J_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        int size = nodes.size();

        mNeigbourElements.resize(size);

        for (int i=0; i<size; i++) mNeigbourElements[i] = 0;

        std::vector<Node::Pointer>::iterator it;

        Vector<double> vCero(3, 0.00);

        for (it=nodes.begin(); it!=nodes.end(); ++it)
        {
            (*mpModel)(CURRENT_DENSITY_SMOOTHED, **it) = vCero;
        }
    }

    //*****************************************************************
    //* - Ending calculation of static J current density.
    //*****************************************************************
    void Modeler::EndStatic_J_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        Vector<double> static_J;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeigbourElements[(*it)->Id()-1];

            static_J = (*mpModel)(CURRENT_DENSITY_SMOOTHED, **it);

            if (NumNeighbours > 0)
            {
                static_J /= NumNeighbours;
            }

            (*mpModel)(CURRENT_DENSITY_SMOOTHED, **it) = static_J;
        }
    }

    //*****************************************************************
    //* - Initial set up to calculate static Joule heating.
    //*****************************************************************
    void Modeler::IniStatic_JouleH_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        int size = nodes.size();

        mNeigbourElements.resize(size);

        for (int i=0; i<size; i++) mNeigbourElements[i] = 0;

        std::vector<Node::Pointer>::iterator it;

        for (it=nodes.begin(); it!=nodes.end(); ++it)
        {
            (*mpModel)(JOULE_HEATING_SMOOTHED, **it) = 0.0;
        }
    }

    //*****************************************************************
    //* - Ending calculation of static Joule heating.
    //*****************************************************************
    void Modeler::EndStatic_JouleH_Derivation()
    {
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        double static_JouleH;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeigbourElements[(*it)->Id()-1];

            static_JouleH = (*mpModel)(JOULE_HEATING_SMOOTHED, **it);

            if (NumNeighbours > 0)
            {
                static_JouleH /= NumNeighbours;
            }

            (*mpModel)(JOULE_HEATING_SMOOTHED, **it) = static_JouleH;
        }
    }

	//*****************************************************************
    //* - Initial set up to calculate H field.
    //*****************************************************************
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

    //******************************************************************************
    //* - Initial set up to calculate J.
    //******************************************************************************
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

	//******************************************************************************
    //* - Calculate H field in elements.
    //******************************************************************************
	void Modeler::CalculateHElement(int* NodesId, unsigned int PropertiesId)
	{
		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
		{
			mNeigbourElements[(*hoit)-1] += 1;
            //mNeigbourElements[(*hoit)-1] = 0;
			nodes.push_back(mpModel->GetNode((*hoit)-1));
		}

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        this->SetProperties(PropertiesId, FREQUENCY, mProblemFrequency);

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_1bb( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st( nodes, properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd( nodes, properties ) );
		else if ( mElementOrder == 3 ) pElement = Element::Pointer( new VolumeElement_3th( nodes, properties ) );
		else if ( mElementOrder == 4 ) pElement = Element::Pointer( new VolumeElement_1bb( nodes, properties ) );//pElement = Element::Pointer( new VolumeElement_4th( nodes, properties ) );

        if ( (mElementOrder == 0) || (mElementOrder == 4) ) 
        {
            pElement->Set_Mesh_Connectivities( mPropertiesMap, mNodesConnectivity, mElementsVector, mType_Of_BC_Normal);
            pElement->Set_Global_Derivatives ( mbVolumes, mb_dNk_dX, mb_dNk_dY, mb_dNk_dZ );
        }
                                                                                                 
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

			v_temp = (*mpModel)(REAL_H, **node_it);

            //double normrH = std::sqrt( rH_vector[0]*rH_vector[0] + rH_vector[1]*rH_vector[1] + rH_vector[2]*rH_vector[2] );
            //double normrV = std::sqrt( v_temp[0]*v_temp[0] + v_temp[1]*v_temp[1] + v_temp[2]*v_temp[2] );
            //if ( normrH > normrV ) (*mpModel)(REAL_H, **node_it) = rH_vector;
            //if ( (normrV == 0.00)  || ( normrH < normrV ) ) (*mpModel)(REAL_H, **node_it) = rH_vector;
           
            v_temp += rH_vector;
            (*mpModel)(REAL_H, **node_it) = v_temp;

            v_temp  = (*mpModel)(IMAG_H, **node_it);

            //double normiH = std::sqrt( iH_vector[0]*iH_vector[0] + iH_vector[1]*iH_vector[1] + iH_vector[2]*iH_vector[2] );
            //double normiV = std::sqrt( v_temp[0]*v_temp[0] + v_temp[1]*v_temp[1] + v_temp[2]*v_temp[2] );
            //if ( normiH > normiV ) (*mpModel)(IMAG_H, **node_it) = iH_vector;
            //if ( (normiV == 0.00)  || ( normiH < normiV ) ) (*mpModel)(IMAG_H, **node_it) = iH_vector;

            v_temp += iH_vector;
            (*mpModel)(IMAG_H, **node_it) = v_temp;

			n++;
        }
	}

    //*****************************************************************************************
    //* - Calculate static E field E=-grad(V) in the Gauss points of the element.
    //*****************************************************************************************
    void Modeler::Calculate_Static_E_Element_GP( Vector<Vector<double> >& Static_E_OnElementGPs,
                                                 int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector<Vector<double> > gradOnGaussPs;

        pElement->CalculateGradientOnGaussPs( gradOnGaussPs );

        // E = -gradV in the central Gauss point
        Vector<double> Static_E_OnCentralGP = gradOnGaussPs[0];
        Static_E_OnCentralGP *= -1.00;

        // Adding gradient to the element list
        Static_E_OnElementGPs.push_back( Static_E_OnCentralGP );
    }

    //*****************************************************************************************
    //* - Calculates static J field J=-sigma*grad(V) in the Gauss points of the element.
    //*****************************************************************************************
    void Modeler::Calculate_Static_J_Element_GP( Vector<Vector<double> >& Static_J_OnElementGPs,
                                                 int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector< Vector<double> > gradOnGaussPs;

        pElement->CalculateGradientOnGaussPs( gradOnGaussPs );

        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        // J = -sigma*gradV in the central Gauss point
        Vector<double> Static_J_OnCentralGP = gradOnGaussPs[0];

        Static_J_OnCentralGP *= -sigma;

        // Adding gradient to the element list
        Static_J_OnElementGPs.push_back( Static_J_OnCentralGP );
    }

    //******************************************************************************************************************
    //* - It writes element nodes Ids and the calculated static current on the export currents file.
    //******************************************************************************************************************
    void Modeler::Export_Static_J_Element_GP( int* NodesId, unsigned int PropertiesId )
    {
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        // Checking if the element has an electrical conductivity > 0.0
        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);
        
        if ( sigma == 0.0 ) return;
        
        // Continue the computation of J only in elements of conductive materials
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back( mpModel->GetNode((*hoit)-1) );
        }

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector< Vector<double> > gradOnGaussPs;

        pElement->CalculateGradientOnGaussPs( gradOnGaussPs );

        // J = -sigma*gradV in the central Gauss point
        Vector<double> Static_J_OnCentralGP = gradOnGaussPs[0];

        Static_J_OnCentralGP *= -sigma;

        double modJ2 = Static_J_OnCentralGP[0]*Static_J_OnCentralGP[0] + 
                       Static_J_OnCentralGP[1]*Static_J_OnCentralGP[1] +
                       Static_J_OnCentralGP[2]*Static_J_OnCentralGP[2] ;

        // If |J| == 0.0 do not write the result in the export file
        if ( modJ2 <= 0.0 ) return;

        // Adding element nodes to the export file
        mExportCurrentFile << NodesId[0] << "  " << NodesId[1] << "  " 
                           << NodesId[2] << "  " << NodesId[3] << "  "; 
          
        // Adding current source components to the export file
        mExportCurrentFile << Static_J_OnCentralGP[0] << "  " 
                           << Static_J_OnCentralGP[1] << "  " 
                           << Static_J_OnCentralGP[2] ;

        // End of line
        mExportCurrentFile << std::endl; 
    }

    //*********************************************************************************************
    //* - It opens the export current file and writes the phase of the current.
    //*********************************************************************************************
    void Modeler::Ini_Exporting_Static_Current()
    {
        String ExpCurrent_FileName; 
        double ExpCurrent_Phase;

        // Reading info from *-9.dat to export currents
        ReadExpCurrentInfoFromFile( mBaseFileName, ExpCurrent_FileName, ExpCurrent_Phase );

        std::cout << "Writing current density in file " << ExpCurrent_FileName << "..." << std::endl;

        // Creating directory to store exported currents
        String ExpJs_Dir_Name = "Exp_J_Sources/";

        #ifdef _WIN64
            mkdir( ExpJs_Dir_Name );
        #else
            mkdir( ExpJs_Dir_Name, 0755 );
        #endif

        // Open file to write exported currents
        mExportCurrentFile.open( ExpJs_Dir_Name + ExpCurrent_FileName, std::fstream::out );

        // Writing phase
        mExportCurrentFile << ExpCurrent_Phase << std::endl;        
    }

    //*********************************************************************************************
    //* - It closes the export current file.
    //*********************************************************************************************
    void Modeler::End_Exporting_Static_Current()
    {
        mExportCurrentFile.close();
    }

    //*********************************************************************************************
    //* - Add to the b_vector the J currents from the folder "Exp_J_Sources/"
    //*********************************************************************************************
    void Modeler::Build_Imported_Static_Currents()
    {        
        // Check if import currents mode is activated
        if ( mImportCurrents == false ) return; 

        // Objectes from "dirent.h" 
        DIR           *ExpJs_Dir;
        struct dirent *Dir_Ent  ;

        // Name of the folder where the exported J currents are stored
        String ExpJs_Dir_Name = "Exp_J_Sources/";

        // If the directory "Exp_J_Sources/" does not exist then exit function
        if ( ( ExpJs_Dir = opendir( ExpJs_Dir_Name ) ) == NULL ) return;
       
        // Read all the files inside the directoy "Exp_J_Sources/"
        while ( ( Dir_Ent = readdir( ExpJs_Dir ) ) != NULL ) 
        {
            // If the entity is not a file go to the next entity
            if ( Dir_Ent->d_type != DT_REG ) continue; 
                      
            // Saving file name as a string
            String File_Name( Dir_Ent->d_name );

            std::cout << "Importing J sources from " << File_Name << "..." << std::endl;

            // Open J source file
            std::fstream ImpCurrentFile( ExpJs_Dir_Name + File_Name, std::fstream::in );  

            int    NodesId [4];
            double J_vector[3]; 
            double J_phase;

            // Reading phase from file
            ImpCurrentFile >> J_phase;

            while( !ImpCurrentFile.eof() )
            {
                // Reading Node Ids from file
                ImpCurrentFile >> NodesId[0]; 
                ImpCurrentFile >> NodesId[1];
                ImpCurrentFile >> NodesId[2];
                ImpCurrentFile >> NodesId[3];
                
                // Reading J source from file
                ImpCurrentFile >> J_vector[0]; 
                ImpCurrentFile >> J_vector[1];
                ImpCurrentFile >> J_vector[2];

                // Assembling imported source element to b_vector
                Assemble_Imported_Current_Element( NodesId, J_vector, J_phase );
            }

            // Closing J source file
            ImpCurrentFile.close(); 
        }

        // Close directory Exp_J_Sources/
        closedir( ExpJs_Dir );
    }

	//*************************************************************************************************
    //* - Assembling imported current element
    //*************************************************************************************************
    void Modeler::Assemble_Imported_Current_Element( int* NodesId, double* Jvec, double Jpha )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume( NodesId, HONodesId );

        for( std::vector<int>::iterator hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) 
        {
            nodes.push_back( mpModel->GetNode( (*hoit) - 1 ) );
        }

        int numNodes;

        if      ( mElementOrder == 0 ) numNodes = 17;
        else if ( mElementOrder == 1 ) numNodes =  4;
		else if ( mElementOrder == 2 ) numNodes = 10;
		else if ( mElementOrder == 3 ) numNodes = 20;
		else if ( mElementOrder == 4 ) numNodes = 35;

        int numDofs = numNodes * 3;

        // DOFs Id vector
        Vector<int> EleIdVector( numDofs );
       
        for ( int i=0; i<numNodes; i++ ) 
        { 
            EleIdVector[i              ] = nodes[i]->pDofcEx()->EquationId();
            EleIdVector[i + numNodes   ] = nodes[i]->pDofcEy()->EquationId();
			EleIdVector[i +(numNodes*2)] = nodes[i]->pDofcEz()->EquationId();
        }

        // Residual vector
        Vector< std::complex<double> > EleResVector( numDofs );
        
        // Initializing residual vector
        for ( int i = 0; i<numDofs; i++ ) EleResVector = std::complex<double>( 0.00, 0.00 );

        // j*w*J
        double freq = mProblemFrequency;

        std::complex<double> cJx( -freq * Jvec[0] * sin( Jpha ), freq * Jvec[0] * cos( Jpha ) );
        std::complex<double> cJy( -freq * Jvec[1] * sin( Jpha ), freq * Jvec[1] * cos( Jpha ) );
        std::complex<double> cJz( -freq * Jvec[2] * sin( Jpha ), freq * Jvec[2] * cos( Jpha ) );
 
        // Gauss points and FEM base
        int nGaussPoints;

        std::vector<double> cX, cY, cZ, W;

        Matrix<double> N;

        if ( mElementOrder == 0 ) 
        { 
            nGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W ); 
            Lagrange3D_Ni_2bb( N, cX, cY, cZ );
        } 
        else if ( mElementOrder == 1 ) 
        { 
            nGaussPoints = GaussPoints3D_Order3( cX, cY, cZ, W ); 
            Lagrange3D_Ni_1st( N, cX, cY, cZ );
        } 
		else if ( mElementOrder == 2 ) 
        { 
            nGaussPoints = GaussPoints3D_Order4( cX, cY, cZ, W ); 
            Lagrange3D_Ni_2nd( N, cX, cY, cZ );
        } 
		else if ( mElementOrder == 3 )
        { 
            nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W ); 
            Lagrange3D_Ni_3th( N, cX, cY, cZ );
        } 
		else if ( mElementOrder == 4 )
        { 
            nGaussPoints = GaussPoints3D_Order6( cX, cY, cZ, W ); 
            Lagrange3D_Ni_4th( N, cX, cY, cZ );
        } 

        // Jacobian
        double X[4], Y[4], Z[4];
        double vol, jacob;
        
        // Define X,Y,Z vector from 1 (instead of 0) to keep the same vol formula as in the elements 
        X[1] = nodes[0]->X(); X[2] = nodes[1]->X(); X[3] = nodes[2]->X(); X[4] = nodes[3]->X();
        Y[1] = nodes[0]->Y(); Y[2] = nodes[1]->Y(); Y[3] = nodes[2]->Y(); Y[4] = nodes[3]->Y();
        Z[1] = nodes[0]->Z(); Z[2] = nodes[1]->Z(); Z[3] = nodes[2]->Z(); Z[4] = nodes[3]->Z();

        // Volume formula (copy&paste from the source elements)
        vol = + X[2]*Y[3]*Z[4] + X[4]*Y[2]*Z[3] + X[3]*Y[4]*Z[2] - X[4]*Y[3]*Z[2] - X[2]*Y[4]*Z[3] - X[3]*Y[2]*Z[4] 
              - X[1]*Y[3]*Z[4] - X[4]*Y[1]*Z[3] - X[3]*Y[4]*Z[1] + X[4]*Y[3]*Z[1] + X[1]*Y[4]*Z[3] + X[3]*Y[1]*Z[4] 
              + X[1]*Y[2]*Z[4] + X[4]*Y[1]*Z[2] + X[2]*Y[4]*Z[1] - X[4]*Y[2]*Z[1] - X[1]*Y[4]*Z[2] - X[2]*Y[1]*Z[4] 
              - X[1]*Y[2]*Z[3] - X[3]*Y[1]*Z[2] - X[2]*Y[3]*Z[1] + X[3]*Y[2]*Z[1] + X[1]*Y[3]*Z[2] + X[2]*Y[1]*Z[3];
         
        jacob = fabs( vol );

        for ( int i = 0; i < numNodes; i++ )
        {
        	double intNi = 0.00;
        
        	for ( int gp = 0; gp < nGaussPoints; gp++ ) intNi += W[gp] * N[i][gp];
        
        	EleResVector[i              ] += cJx * ( jacob * intNi );
        	EleResVector[i + numNodes   ] += cJy * ( jacob * intNi );
        	EleResVector[i +(numNodes*2)] += cJz * ( jacob * intNi );
        }

		// Inv(T)*b for PBC condition
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );
	    PBC_RightLeft_TR( HONodesId, EleIdVector, EleResVector );
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );

		// Inv(T)*b for contact nodes
        if ( mPlasmaMode ) ContactHO_TR_Plasma( HONodesId, EleIdVector, EleResVector );
        else               ContactHO_TR       ( HONodesId, EleIdVector, EleResVector );

        // Inv(T)*b for PEC, PMC conditions
        HO_TR( HONodesId, EleResVector );

		// Inv(T)*b for axisymmetric mode
        if ( mAxisymmetric ) Axisym_TR( HONodesId, EleResVector );

        // Inv(T)*b for plasma mode
        if ( mPlasmaMode   ) Plasma_TR( HONodesId, EleResVector );

        //Assembling to global source vector
        int eleSize = EleIdVector.size();

		for ( int i_local=0; i_local < eleSize; i_local++ )
		{
			int i_global = EleIdVector[i_local];

			if ( i_global < mSystemSize )
			{
				b_vector[i_global] += EleResVector[i_local];
			}
		}

        EleIdVector.FreeData ();
        EleResVector.FreeData();
    }

    //*********************************************************************************************
    //* - Calculates Joule heating Q = sigma * ||grad(V)||^2 in the Gauss points of the element.
    //*********************************************************************************************
    void Modeler::Calculate_Static_JouleH_Element_GP( Vector<double>& Static_JouleH_OnElementGPs,
                                                      int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        // Grad(V) on Gauss points
        Vector< Vector<double> > gradOnGaussPs;

        pElement->CalculateGradientOnGaussPs( gradOnGaussPs );

        // Square of the module of grad(V) on the central Gauss point [0]
        double modGradV2 = gradOnGaussPs[0][0]*gradOnGaussPs[0][0] +
                           gradOnGaussPs[0][1]*gradOnGaussPs[0][1] +
                           gradOnGaussPs[0][2]*gradOnGaussPs[0][2] ;

        // Electrical conductivity of the element
        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        // Joule heating
        double JouleHeating = sigma * modGradV2;

        // Adding gradient to the element list
        Static_JouleH_OnElementGPs.push_back( JouleHeating );
    }

    //*****************************************************************************************
    //* - Calculate static E field E = -grad(V) in elements.
    //*****************************************************************************************
    void Modeler::Calculate_Static_E_Element( int* NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            mNeigbourElements[(*hoit)-1] += 1;
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector<Vector<double> > gradOnNodes;

        pElement->CalculateGradientOnNodes( gradOnNodes );

        Vector<Node::Pointer>::iterator node_it;

        int n = 0;

        for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            (*mpModel)(ELECTRIC_FIELD_SMOOTHED, **node_it) -= gradOnNodes[n];

            n++;
        }
    }

    //*****************************************************************************************
    //* - Calculate static J current density J = sigma * (-grad(V)).
    //*****************************************************************************************
    void Modeler::Calculate_Static_J_Element(int* NodesId, unsigned int PropertiesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            mNeigbourElements[(*hoit)-1] += 1;
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector<Vector<double> > gradOnNodes;

        pElement->CalculateGradientOnNodes( gradOnNodes );

        Vector<Node::Pointer>::iterator node_it;

        int n = 0;

        for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            Vector<double> temp_J = gradOnNodes[n];

            temp_J *= sigma;

            (*mpModel)(CURRENT_DENSITY_SMOOTHED, **node_it) -= temp_J;

            n++;
        }
    }

    //*****************************************************************************************
    //* - Calculate static Joule heating Q = sigma * ||grad(V)||^2.
    //*****************************************************************************************
    void Modeler::Calculate_Static_JouleH_Element(int* NodesId, unsigned int PropertiesId)
    {
        std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            mNeigbourElements[(*hoit)-1] += 1;
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Static(nodes,properties) );

        Vector< Vector<double> > gradOnNodes;

        pElement->CalculateGradientOnNodes( gradOnNodes );

        Vector<Node::Pointer>::iterator node_it;

        double temp_JouleH;

        int n = 0;

        for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            temp_JouleH = sigma * ( gradOnNodes[n][0] * gradOnNodes[n][0] +
                                    gradOnNodes[n][1] * gradOnNodes[n][1] +
                                    gradOnNodes[n][2] * gradOnNodes[n][2]  );

            (*mpModel)(JOULE_HEATING_SMOOTHED, **node_it) += temp_JouleH;

            n++;
        }
    }

	//***************************************************************************
    //* - Calculate Joule heating in elements.
    //***************************************************************************
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

			(*mpModel)(JOULE_HEATING, **it) = 0.5 * eff_cond * modE2;
		}
	}

    //***********************************************************************************
    //* - Calculate Joule heating in plasma elements.
    //***********************************************************************************
	void Modeler::CalculateJouleHElement_Plasma(int* NodesId, unsigned int PropertiesId)
	{
		// Material properties
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        double sigma    = (*properties)( IHL_ELECTRIC_CONDUCTIVITY  );

        double epr_real = (*properties)( REAL_ELECTRIC_PERMITTIVITY );
        double epr_imag = (*properties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mur_real = (*properties)( REAL_MAGNETIC_PERMEABILITY );
        double mur_imag = (*properties)( IMAG_MAGNETIC_PERMEABILITY );

        // In Plasma_Mode if the material is not vacuum then it is considered as an IHL material.
        if ( (sigma    != 0.0) ||
             (epr_real != 1.0) || (epr_imag != 0.0) ||
             (mur_real != 1.0) || (mur_imag != 0.0)  )
        {
            CalculateJouleHElement(NodesId, PropertiesId);
            return;
        }

        // In Plasma_Mode if the material is vacuum then it is considered as a plasma.
        std::vector<Node::Pointer> nodes;

		std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId, HONodesId);

		std::vector<int>::iterator hoit;

		for (hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
		{
			nodes.push_back(mpModel->GetNode((*hoit) - 1));
		}

		Vector<Node::Pointer>::iterator node_it;

        Matrix<std::complex<double> > TEn(3,3);

		Vector<std::complex<double> > cDv(3);
		Vector<std::complex<double> > cEv(3);

		std::complex<double> EconjD;

		for (node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
		{
			cEv[0] = (*mpModel)(cEx, **node_it);
			cEv[1] = (*mpModel)(cEy, **node_it);
			cEv[2] = (*mpModel)(cEz, **node_it);

            mpColdPlasma->Get_PermittivityTensor_InNode(TEn, *node_it, mProblemFrequency);

			cDv[0] = TEn[0][0]*cEv[0] + TEn[0][1]*cEv[1] + TEn[0][2]*cEv[2];
			cDv[1] = TEn[1][0]*cEv[0] + TEn[1][1]*cEv[1] + TEn[1][2]*cEv[2];
			cDv[2] = TEn[2][0]*cEv[0] + TEn[2][1]*cEv[1] + TEn[2][2]*cEv[2];

            EconjD = cEv[0] * std::conj( cDv[0] ) +
                     cEv[1] * std::conj( cDv[1] ) +
                     cEv[2] * std::conj( cDv[2] ) ;

            (*mpModel)(JOULE_HEATING, **node_it) = - 0.5 * mProblemFrequency * std::imag(EconjD);
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

	//****************************************************************************
    //* - Calculate J with VE elements in plasma mode
    //****************************************************************************
    void Modeler::CalculateJVEElement_Plasma(int* NodesId, unsigned int PropertiesId)
    {
        // Material properties
        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        double sigma    = (*properties)( IHL_ELECTRIC_CONDUCTIVITY  );

        double epr_real = (*properties)( REAL_ELECTRIC_PERMITTIVITY );
        double epr_imag = (*properties)( IMAG_ELECTRIC_PERMITTIVITY );

        double mur_real = (*properties)( REAL_MAGNETIC_PERMEABILITY );
        double mur_imag = (*properties)( IMAG_MAGNETIC_PERMEABILITY );

        // In Plasma_Mode if the material is not vacuum then it is considered as an IHL material.
        if ( (sigma    != 0.0) ||
             (epr_real != 1.0) || (epr_imag != 0.0) ||
             (mur_real != 1.0) || (mur_imag != 0.0)  )
        {
            CalculateJVEElement(NodesId, PropertiesId);
            return;
        }

        // In Plasma_Mode if the material is vacuum then it is considered as a plasma.
        std::vector<Node::Pointer> nodes;

		std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId, HONodesId);

		std::vector<int>::iterator hoit;

		for (hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)
		{
			nodes.push_back(mpModel->GetNode((*hoit) - 1));
		}

		Vector<Node::Pointer>::iterator node_it;

        Matrix<std::complex<double> > TCn(3,3);

		Vector<std::complex<double> > cJv(3);
		Vector<std::complex<double> > cEv(3);

		Vector<double> rJv(3);
		Vector<double> iJv(3);

		for (node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
		{
			cEv[0] = (*mpModel)(cEx, **node_it);
			cEv[1] = (*mpModel)(cEy, **node_it);
			cEv[2] = (*mpModel)(cEz, **node_it);

			mpColdPlasma->Get_ConductivityTensor_InNode(TCn, *node_it, mProblemFrequency);

			cJv[0] = TCn[0][0]*cEv[0] + TCn[0][1]*cEv[1] + TCn[0][2]*cEv[2];
			cJv[1] = TCn[1][0]*cEv[0] + TCn[1][1]*cEv[1] + TCn[1][2]*cEv[2];
			cJv[2] = TCn[2][0]*cEv[0] + TCn[2][1]*cEv[1] + TCn[2][2]*cEv[2];

			rJv[0] = std::real( cJv[0] );
			rJv[1] = std::real( cJv[1] );
			rJv[2] = std::real( cJv[2] );

			(*mpModel)(REAL_J, **node_it) = rJv;

			iJv[0] = std::imag( cJv[0] );
			iJv[1] = std::imag( cJv[1] );
			iJv[2] = std::imag( cJv[2] );

			(*mpModel)(IMAG_J, **node_it) = iJv;
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

		std::vector<double> vJa   = (*properties)(COMPLEX_IBC);
		std::vector<double> Jaxis = (*properties)(COMPLEX_NEUMANN_FLOW);

		std::complex<double> cJa( vJa[0]*cos(vJa[1]), vJa[0]*sin(vJa[1]) );

		std::vector<Node::Pointer> nodes;

        std::vector<int> HONodesId;

		PushHONodesOnVolume(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  nodes.push_back(mpModel->GetNode((*hoit)-1));

		Vector<Node::Pointer>::iterator node_it;

		Vector<double> J_vector(3);

		int nodeId;

		for(node_it = nodes.begin(); node_it != nodes.end(); ++node_it)
        {
			// Position vector
			double px = (*node_it)->X() - Jaxis[0];
			double py = (*node_it)->Y() - Jaxis[1];
			double pz = (*node_it)->Z() - Jaxis[2];

			// Rotation axis
			double ax = Jaxis[3];
			double ay = Jaxis[4];
			double az = Jaxis[5];

			// r = a x p
			double rx = ay*pz - az*py;
			double ry = az*px - ax*pz;
			double rz = ax*py - ay*px;

			double norm = sqrt( rx*rx + ry*ry + rz*rz );

			if (norm > 0.00)
			{
				rx = rx / norm;
				ry = ry / norm;
				rz = rz / norm;
			}
			else
			{
				rx = 0.00;
				ry = 0.00;
				rz = 0.00;
			}

			nodeId = (*node_it)->Id();

			J_vector[0] = std::real(cJx) + std::real(cJa)*rx;
            J_vector[1] = std::real(cJy) + std::real(cJa)*ry;
            J_vector[2] = std::real(cJz) + std::real(cJa)*rz;

			tmp_rJ[nodeId] = J_vector;

			J_vector[0] = std::imag(cJx) + std::imag(cJa)*rx;
            J_vector[1] = std::imag(cJy) + std::imag(cJa)*ry;
            J_vector[2] = std::imag(cJz) + std::imag(cJa)*rz;

			tmp_iJ[nodeId] = J_vector;
        }
	}

	//*********************************************************************************************************
    //* - Add imported currents to element nodes
    //*********************************************************************************************************
	void Modeler::Calculate_Imported_Static_Currents()
	{
        // Check if import currents mode is activated
        if ( mImportCurrents == false ) return;

        // Objectes from "dirent.h" 
        DIR           *ExpJs_Dir;
        struct dirent *Dir_Ent  ;

        // Name of the folder where the exported J currents are stored
        String ExpJs_Dir_Name = "Exp_J_Sources/";

        // If the directory "Exp_J_Sources/" does not exist then exit function
        if ( ( ExpJs_Dir = opendir( ExpJs_Dir_Name ) ) == NULL ) return;
       
        // Read all the files inside the directoy "Exp_J_Sources/"
        while ( ( Dir_Ent = readdir( ExpJs_Dir ) ) != NULL ) 
        {
            // If the entity is not a file go to the next entity
            if ( Dir_Ent->d_type != DT_REG ) continue; 
                      
            // Saving file name as a string
            String File_Name( Dir_Ent->d_name );

            // Open J source file
            std::fstream ImpCurrentFile( ExpJs_Dir_Name + File_Name, std::fstream::in );  

            // Info to read from file
            double         J_phase    ;
            Vector<int>    NodesId (4);
            Vector<double> J_vector(3);

            // Dummy vectors
            Vector<double> rtmp(3);
            Vector<double> itmp(3);

            // Counter of neighbords elements for smoothing
            std::map<unsigned int, int> numNeighbors;

            // Counters 
            int i, j;
           
            // First line of the file (phase)
            ImpCurrentFile >> J_phase;

            while( !ImpCurrentFile.eof() )
            {
                // Reading Node Ids from file
                ImpCurrentFile >> NodesId[0]; 
                ImpCurrentFile >> NodesId[1];
                ImpCurrentFile >> NodesId[2];
                ImpCurrentFile >> NodesId[3];
                
                // Reading J source from fila
                ImpCurrentFile >> J_vector[0]; 
                ImpCurrentFile >> J_vector[1];
                ImpCurrentFile >> J_vector[2];

                for ( i=0; i<4; i++ )
                {
                    int nodeId = NodesId[i];

                    // Check for the first occurrence of the node
                    if ( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() )
                    {
                        for ( j=0; j<3; j++ )
                        {
                            tmp_ist_rJ[nodeId][j] += J_vector[j] * cos( J_phase ); 
                            tmp_ist_iJ[nodeId][j] += J_vector[j] * sin( J_phase ); 
                        }

                        numNeighbors[nodeId]++;
                    }
                    else
                    {
                        for ( j=0; j<3; j++ )
                        {
                            rtmp[j] = J_vector[j] * cos( J_phase ); 
                            itmp[j] = J_vector[j] * sin( J_phase ); 
                        }

                        tmp_ist_rJ[nodeId] = rtmp;
                        tmp_ist_iJ[nodeId] = itmp;

                        numNeighbors[nodeId] = 1;
                    }
                }
            }

            // Smoothing
            std::map<unsigned int, int>::iterator nit;
            for ( nit = numNeighbors.begin(); nit != numNeighbors.end(); nit++ )
            {
                rtmp = tmp_ist_rJ[nit->first]; 
                itmp = tmp_ist_iJ[nit->first];

                for ( j=0; j<3; j++ )
                {
                    rtmp[j] /= nit->second; 
                    itmp[j] /= nit->second; 
                }

                tmp_ist_rJ[nit->first] = rtmp; 
                tmp_ist_iJ[nit->first] = itmp; 
            }

            // Close file
            ImpCurrentFile.close();

            // Cleaning neighbords counting map
            numNeighbors.clear(); std::map<unsigned int, int>().swap( numNeighbors );
        }

        // Close directory Exp_J_Sources/
        closedir( ExpJs_Dir );
	}

	//***********************************************************************************
    //* - Adding all the current contributions ( eddy, imposed and imported )
    //***********************************************************************************
	void Modeler::EndJCalculation()
	{
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		std::vector<Node::Pointer>::iterator node_it;

		Vector<double> rEddy(3, 0.00); Vector<double> iEddy(3, 0.00);
		Vector<double> rJtot(3, 0.00); Vector<double> iJtot(3, 0.00);

        for ( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            // Loop index
			int i;

            // Nodes Id
            int nodeId = (*node_it)->Id();

            // Eddy currents contribution
			rEddy = (*mpModel)(REAL_J, **node_it);
            iEddy = (*mpModel)(IMAG_J, **node_it);

            for ( i=0; i<3; i++ )
            {
                rJtot[i] = rEddy[i];
                iJtot[i] = iEddy[i]; 
            }

            // Adding imposed currents contribution
			if ( tmp_rJ.find( nodeId ) != tmp_rJ.end() )
			{
                for ( i=0; i<3; i++ )
                {
                    rJtot[i] += tmp_rJ[nodeId][i];
                    iJtot[i] += tmp_iJ[nodeId][i];
                }
			}

            // Adding imported static currents contribution
            if ( ( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() ) && 
                 ( tmp_rJ.find    ( nodeId ) == tmp_rJ.end    () )  )
			{
                for ( i=0; i<3; i++ )
                {
                    rJtot[i] += tmp_ist_rJ[nodeId][i];
                    iJtot[i] += tmp_ist_iJ[nodeId][i];
                }
			}

            // Total current distribution
            (*mpModel)(REAL_J, **node_it) = rJtot;
			(*mpModel)(IMAG_J, **node_it) = iJtot;

			// Module of the total current
	        (*mpModel)(MOD_J , **node_it) = sqrt( rJtot[0]*rJtot[0] + rJtot[1]*rJtot[1] + rJtot[2]*rJtot[2] +
                                                  iJtot[0]*iJtot[0] + iJtot[1]*iJtot[1] + iJtot[2]*iJtot[2] );
        }

		tmp_rJ.clear()    ; std::map<unsigned int, std::vector<double> >().swap(tmp_rJ    );
		tmp_iJ.clear()    ; std::map<unsigned int, std::vector<double> >().swap(tmp_iJ    );

        tmp_ist_rJ.clear(); std::map<unsigned int, std::vector<double> >().swap(tmp_ist_rJ);
		tmp_ist_iJ.clear(); std::map<unsigned int, std::vector<double> >().swap(tmp_ist_iJ);
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
        std::map<double, std::complex<double> >().swap(mNormalization);

	    mProjection.clear();
        std::map<double, std::complex<double> >().swap(mProjection);
	}

} // Namespace Kratos
