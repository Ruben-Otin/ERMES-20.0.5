
#include "../ERMES/VolumeElement_1st_FullWave.h"
#include "../ERMES/JSource_1st_FullWave.h"
#include "../ERMES/GenericRobin_1st_FullWave.h"
#include "../ERMES/FarField_1st_FullWave.h"
#include "../ERMES/CoaxialPortTEM_1st_FullWave.h"
#include "../ERMES/RWPortTE10_1st_FullWave.h"

#include "../ERMES/VolumeElement_2nd_FullWave.h"
#include "../ERMES/JSource_2nd_FullWave.h"
#include "../ERMES/GenericRobin_2ndQ_FullWave.h"
#include "../ERMES/FarField_2ndQ_FullWave.h"
#include "../ERMES/CoaxialPortTEM_2ndQ_FullWave.h"
#include "../ERMES/RWPortTE10_2ndQ_FullWave.h"

#include "../ERMES/VolumeElement_3sb_FullWave.h"
#include "../ERMES/JSource_3sb_FullWave.h"
#include "../ERMES/GenericRobin_3sb_FullWave.h"
#include "../ERMES/FarField_3sb_FullWave.h"
#include "../ERMES/CoaxialPortTEM_3sb_FullWave.h"
#include "../ERMES/RWPortTE10_3sb_FullWave.h"

#include "../ERMES/VolumeElement_1st_ColdPlasma.h"
#include "../ERMES/FarField_1st_ColdPlasma.h"
#include "../ERMES/GenericRobin_1st_ColdPlasma.h"

#include "../ERMES/VolumeElement_2ndQ_ColdPlasma.h"
#include "../ERMES/FarField_2ndQ_ColdPlasma.h"
#include "../ERMES/GenericRobin_2ndQ_ColdPlasma.h"

#include "../ERMES/VolumeElement_3sb_ColdPlasma.h"
#include "../ERMES/FarField_3sb_ColdPlasma.h"
#include "../ERMES/GenericRobin_3sb_ColdPlasma.h"

#include "../ERMES/VolumeElement_1st_Electrostatic.h"
#include "../ERMES/GenericRobin_1st_Electrostatic.h"

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
    //*****************************************************************************
    //* - Generate nodes
    //*****************************************************************************
    void Modeler::GenerateNode( IndexType NodeId, double X, double Y, double Z )
    {
        Node::Pointer new_node( new Node( NodeId, X, Y, Z ) );
    
        mpModel->AddNode( new_node );
    
        if ( mElectrostaticMode )
        {
            new_node->Add_V_Dof( mDofSet );
        }
        else if ( mPotentials_On ) 
        {
            new_node->Add_cAV_Dofs( mDofSet );
        }
        else   
        {
            new_node->Add_cE_Dofs( mDofSet );
        }
    }
    
    //*************************************************************************************************
    //* - Print error messages and exit ERMES
    //*************************************************************************************************
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
         exit( EXIT_SUCCESS );
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// - Contact discontinuity in dielectrics surface
	///////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************
    //* - Fill nodes properties vector map
    //*************************************************************************************************
	void Modeler::ContactProperties( int* NodesId, unsigned int PropertiesId )
	{
		 for ( int i=0; i<4; i++ ) 
         {
             mTempNodeProperties[ NodesId[i] ] = PropertiesId;
         }
	}

	//*************************************************************************************************
    //* - Generate pairs of nodes in a contact surface.
    //*************************************************************************************************
	void Modeler::GenerateContactPairs( Vector<int>& CtCNodesId )
	{
        if ( mElectrostaticMode == true ) 
        {
            return;
        }
 
        Vector<int> sR_NodesId;
		Vector<int> sL_NodesId;

		unsigned short int R_PropertyId;
		unsigned short int L_PropertyId;

        //sR_NodesId will have the nodes with higher mod(eps)
		//sL_NodesId will have the nodes with lower  mod(eps)
		Properties::Pointer First3_prop = mpModel->GetProperties( mTempNodeProperties[ CtCNodesId[0] ] );
		Properties::Pointer Last3_prop  = mpModel->GetProperties( mTempNodeProperties[ CtCNodesId[3] ] );

		double freq = mProblemFrequency;
		double eo   = 8.8541878176e-12;

		double sigma_First3    = (*First3_prop)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
		double real_eps_First3 = (*First3_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double imag_eps_First3 = (*First3_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> cEps_First3( real_eps_First3, imag_eps_First3 + ( sigma_First3 / freq ) );

		double mod_epc_First3 = std::abs(cEps_First3);

		double sigma_Last3    = (*Last3_prop)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
		double real_eps_Last3 = (*Last3_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double imag_eps_Last3 = (*Last3_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> cEps_Last3( real_eps_Last3, imag_eps_Last3 + ( sigma_Last3 / freq ) );

		double mod_epc_Last3 = std::abs( cEps_Last3 );

		double mod_epc_R;
        double mod_epc_L;

		if ( mod_epc_First3 > mod_epc_Last3 )
		{
			Vector<int> OneSideNodesId( 3 );

			for ( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i ];
            }

            PushHONodesOnSurface( OneSideNodesId, sR_NodesId );

			for ( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i + 3 ];
            }

            PushHONodesOnSurface( OneSideNodesId, sL_NodesId );

			CalculateContactNormals( sR_NodesId );
			CalculateContactNormals( sL_NodesId );

			R_PropertyId = mTempNodeProperties[ CtCNodesId[0] ];
			L_PropertyId = mTempNodeProperties[ CtCNodesId[3] ];

			mod_epc_R = mod_epc_First3;
			mod_epc_L = mod_epc_Last3 ;
		}
		else
		{
			Vector<int> OneSideNodesId( 3 );

			for ( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i ];
            }
			
            PushHONodesOnSurface( OneSideNodesId, sL_NodesId );

			for ( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i + 3 ];
            }
			
            PushHONodesOnSurface( OneSideNodesId, sR_NodesId );

			CalculateContactNormals( sL_NodesId );
			CalculateContactNormals( sR_NodesId );

			L_PropertyId = mTempNodeProperties[ CtCNodesId[0] ];
			R_PropertyId = mTempNodeProperties[ CtCNodesId[3] ];

			mod_epc_L = mod_epc_First3;
			mod_epc_R = mod_epc_Last3;
		}

        // Setting temporal contact pairs
        if ( mPotentials_On )
        {
            Set_TempContactPairs_A( sR_NodesId, sL_NodesId, mod_epc_L );
        }
        else
        {
            Set_TempContactPairs_E( sR_NodesId, sL_NodesId, mod_epc_L );
        }

        // Setting temporal nodal properties
        Vector<int>::const_iterator R_Idit = sR_NodesId.begin();
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for ( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			mTempNodeProperties[ (*R_Idit) ] = R_PropertyId;
		    mTempNodeProperties[ (*L_Idit) ] = L_PropertyId;
        }
	}

	//**********************************************************************************************************************
    //* - Set temporal contact pairs in E field formulation
    //**********************************************************************************************************************
	void Modeler::Set_TempContactPairs_E( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L )
    {
		Vector<int>::const_iterator R_Idit = sR_NodesId.begin();
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for ( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			int RNodeId = (*R_Idit);
			int LNodeId = (*L_Idit);

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

			if ( pRNode->pDofcEx()->IsFixed() )
			{
				if ( mTempContactPairs[ RNodeId ] != LNodeId )
				{
					int currLNodeId = mTempContactPairs[ RNodeId ];

					Properties::Pointer currL_prop = mpModel->GetProperties( mTempNodeProperties[ currLNodeId ] );

                    double freq = mProblemFrequency;
		            double eo   = 8.8541878176e-12 ;

					double sigma_currL    = (*currL_prop)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
					double real_eps_currL = (*currL_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
					double imag_eps_currL = (*currL_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

					std::complex<double> cEps_currL( real_eps_currL, imag_eps_currL + ( sigma_currL / freq ) );

					double mod_epc_currL = std::abs( cEps_currL );

					if ( mod_epc_currL > mod_epc_L )
					{
						Node::Pointer pcurrLNode = mpModel->GetNode( currLNodeId-1 );

						pcurrLNode->pDofcEx()->FixDof();

						pcurrLNode->pDofcEy()->FreeDof();

						pLNode->pDofcEy()->FixDof();

						mTempContactPairs[ RNodeId     ] = LNodeId;
						mTempContactPairs[ currLNodeId ] = LNodeId;

						std::map<unsigned int, unsigned int>::iterator cPair_it;

						for ( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
						{
							if ( cPair_it->second == currLNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
						}
					}
					else
					{
						pLNode->pDofcEx()->FixDof();

						mTempContactPairs[ LNodeId ] = currLNodeId;
					}
				}
			}
			else if ( pLNode->pDofcEx()->IsFixed() )
			{
				if ( !pRNode->pDofcEy()->IsFixed() )
				{
					pRNode->pDofcEx()->FixDof();

					mTempContactPairs[ RNodeId ] = mTempContactPairs[ LNodeId ];
				}
			}
			else if ( pRNode->pDofcEy()->IsFixed() )
			{
				pRNode->pDofcEx()->FixDof();

			    pRNode->pDofcEy()->FreeDof();

				pLNode->pDofcEy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;

				std::map<unsigned int, unsigned int>::iterator cPair_it;

		        for ( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
				{
					if ( cPair_it->second == RNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
				}
			}
			else
			{
				pRNode->pDofcEx()->FixDof();

				pLNode->pDofcEy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;
			}
		}    
    }

	//**********************************************************************************************************************
    //* - Set temporal contact pairs in potential A formulation
    //**********************************************************************************************************************
	void Modeler::Set_TempContactPairs_A( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L )
    {
		Vector<int>::const_iterator R_Idit = sR_NodesId.begin();
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for ( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			int RNodeId = (*R_Idit);
			int LNodeId = (*L_Idit);

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

			if ( pRNode->pDofcAx()->IsFixed() )
			{
				if ( mTempContactPairs[ RNodeId ] != LNodeId )
				{
					int currLNodeId = mTempContactPairs[ RNodeId ];

					Properties::Pointer currL_prop = mpModel->GetProperties( mTempNodeProperties[ currLNodeId ] );

                    double freq = mProblemFrequency;
		            double eo   = 8.8541878176e-12 ;

					double sigma_currL    = (*currL_prop)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
					double real_eps_currL = (*currL_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
					double imag_eps_currL = (*currL_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

					std::complex<double> cEps_currL( real_eps_currL, imag_eps_currL + ( sigma_currL / freq ) );

					double mod_epc_currL = std::abs( cEps_currL );

					if ( mod_epc_currL > mod_epc_L )
					{
						Node::Pointer pcurrLNode = mpModel->GetNode( currLNodeId-1 );

						pcurrLNode->pDofcAx()->FixDof();

						pcurrLNode->pDofcAy()->FreeDof();

						pLNode->pDofcAy()->FixDof();

						mTempContactPairs[ RNodeId     ] = LNodeId;
						mTempContactPairs[ currLNodeId ] = LNodeId;

						std::map<unsigned int, unsigned int>::iterator cPair_it;

						for ( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
						{
							if ( cPair_it->second == currLNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
						}
					}
					else
					{
						pLNode->pDofcAx()->FixDof();

						mTempContactPairs[ LNodeId ] = currLNodeId;
					}
				}
			}
			else if ( pLNode->pDofcAx()->IsFixed() )
			{
				if ( !pRNode->pDofcAy()->IsFixed() )
				{
					pRNode->pDofcAx()->FixDof();

					mTempContactPairs[ RNodeId ] = mTempContactPairs[ LNodeId ];
				}
			}
			else if ( pRNode->pDofcAy()->IsFixed() )
			{
				pRNode->pDofcAx()->FixDof();

			    pRNode->pDofcAy()->FreeDof();

				pLNode->pDofcAy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;

				std::map<unsigned int, unsigned int>::iterator cPair_it;

		        for ( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
				{
					if ( cPair_it->second == RNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
				}
			}
			else
			{
				pRNode->pDofcAx()->FixDof();

				pLNode->pDofcAy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;
			}
		}    
    }

	//**********************************************************************************************************************
    //* - Calculates geometric averaged normal in contact surfaces
    //**********************************************************************************************************************
	void Modeler::CalculateContactNormals( Vector<int>& NodesId )
	{
		Node::Pointer pNode;

		std::vector<Node::Pointer> nodes;

		std::vector<int>::const_iterator id_it;

		for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			pNode = mpModel->GetNode( (*id_it)-1 );

			nodes.push_back( pNode );
		}

		// Exterior normal of the element
		Vector<double> normal( 3 );

		// Normal = V2 x V1 = 2 * element_area * n
		normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
				    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

		normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
				    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

		normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
				    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the contact normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetContactNormals[*id_it].push_back( normal );
		}
	}

	//**********************************************************************************************************************
    //* - Set average normals in contact surfaces
    //**********************************************************************************************************************
	void Modeler::SetAvgNormalsInContact()
	{
		// Setting average normals in surface dilelectric nodes
		std::map< unsigned int, std::vector<std::vector<double> > >::iterator it_Set;

        for ( it_Set = mSetContactNormals.begin(); it_Set != mSetContactNormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Averaged normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd == true ) 
            {
                normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            }
            else       
            {
                normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );
            }

            // Including the calculated normal in the mNormal vector
            mContactNormals[ it_Set->first ] = normal_atNode;
        }

		mSetContactNormals.clear();
        std::map< unsigned int,std::vector<std::vector<double> > >().swap( mSetContactNormals );
	}

	//**********************************************************************************************************************
    //* - Setting discontinuity on dielectric surface.
    //**********************************************************************************************************************
	void Modeler::SetContacts()
	{
        SetAvgNormalsInContact();

		std::map<unsigned int, unsigned int>::iterator cPair_it;

		for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
        {
			unsigned int RNodeId = cPair_it->first;
			unsigned int LNodeId = mTempContactPairs[ RNodeId ];

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

            std::complex<double> cZero( 0.00, 0.00 );

            if ( mPotentials_On )
            {
                (*mpModel)( cAx, *pRNode) = cZero;
			    (*mpModel)( cAy, *pRNode) = cZero;
			    (*mpModel)( cAz, *pRNode) = cZero;
                (*mpModel)( cVs, *pRNode) = cZero;

                pRNode->pDofcAx()->FixDof();
                pRNode->pDofcAy()->FixDof();
			    pRNode->pDofcAz()->FixDof();

                if ( pRNode->pDof( cVs ) != NULL ) 
                {
                    pRNode->pDofcVs()->FixDof(); 
                }

			    pLNode->pDofcAx()->FreeDof();
			    pLNode->pDofcAy()->FreeDof();
			    pLNode->pDofcAz()->FreeDof();

                if ( pLNode->pDof( cVs ) != NULL ) 
                {
                    pLNode->pDofcVs()->FreeDof(); 
                }

			    if ( !mAVContinuity_On )
                {
                    //pRNode->pDofcAx()->FreeDof();
                }
            }
            else
            {
			    (*mpModel)( cEx, *pRNode ) = cZero;
			    (*mpModel)( cEy, *pRNode ) = cZero;
			    (*mpModel)( cEz, *pRNode ) = cZero;
                
			    pRNode->pDofcEx()->FixDof();
			    pRNode->pDofcEy()->FixDof();
			    pRNode->pDofcEz()->FixDof();
                
			    pLNode->pDofcEx()->FreeDof();
			    pLNode->pDofcEy()->FreeDof();
			    pLNode->pDofcEz()->FreeDof();            
            }

			Vector<double> AvgNormal( 3, 0.00 );

			Vector<double> RNormal( mContactNormals[ RNodeId ] );
			Vector<double> LNormal( mContactNormals[ LNodeId ] );

            // n or -n produces the same transformation matrix.
			AvgNormal[0] = -RNormal[0];
			AvgNormal[1] = -RNormal[1];
			AvgNormal[2] = -RNormal[2];

			// Making the contact pair
			ContactPairData& pContactPair = mContactPairs[ RNodeId ];

			pContactPair.ctcNormal = AvgNormal;
			pContactPair.ctcNode   = LNodeId;
			pContactPair.matR      = mTempNodeProperties[ RNodeId ];
			pContactPair.matL      = mTempNodeProperties[ LNodeId ];
		}

		DelTempContactVectors();
	}

	//**********************************************************************************************************************
    //* - Delete temporal vectors generated along contact discontinuity making.
    //**********************************************************************************************************************
	void Modeler::DelTempContactVectors()
	{
		mTempNodeProperties.clear();
		std::map<unsigned int, unsigned short int>().swap( mTempNodeProperties );

		mContactNormals.clear();
		std::map<unsigned int, std::vector<double> >().swap( mContactNormals );

		mTempContactPairs.clear();
		std::map<unsigned int, unsigned int>().swap( mTempContactPairs );
	}

	//**********************************************************************************************************************
    //* - Calculates "surface smoothing" on discontinuity surfaces (LL2P formulation)
    //**********************************************************************************************************************
    void Modeler::Calculate_Contact_Smoothing_Srf( Vector<int>& CtCNodesId )
    {
        // Reading left-right nodes on discontinuity surface
        Vector<int> ctcNodeIds_1( 3 );
        Vector<int> ctcNodeIds_2( 3 );

        for ( int i=0; i<3; i++ )
        {
            ctcNodeIds_1[ i ] = CtCNodesId[ i     ];
            ctcNodeIds_2[ i ] = CtCNodesId[ i + 3 ];
        }

        // Bubble elements on contact surface
        Vector<int> HONodesId_1; PushHONodesOnSurface( ctcNodeIds_1, HONodesId_1 ); 
        Vector<int> HONodesId_2; PushHONodesOnSurface( ctcNodeIds_2, HONodesId_2 );

        // Pointers to contact nodes 
        Vector<Node::Pointer> ctcNodes_1;
        Vector<Node::Pointer> ctcNodes_2;

        Vector<int>::iterator hoit;

        for ( hoit = HONodesId_1.begin(); hoit != HONodesId_1.end(); hoit++ ) 
        {
            ctcNodes_1.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

        for ( hoit = HONodesId_2.begin(); hoit != HONodesId_2.end(); hoit++ ) 
        {
            ctcNodes_2.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

        // Number of nodes and dofs on contact surface
        int ctcNumNodes = 6;
        int ctcNumDofs  = ctcNumNodes * 3;

		// Electric field on contact nodes
        Matrix< std::complex<double> > cE_1( ctcNumDofs, 1 );
        Matrix< std::complex<double> > cE_2( ctcNumDofs, 1 );
		     
		for( int n=0; n<ctcNumNodes; n++ )
		{
			cE_1[ n                 ][0] = (*mpModel)( cEx, *ctcNodes_1[n] );
			cE_1[ n + ctcNumNodes   ][0] = (*mpModel)( cEy, *ctcNodes_1[n] );
			cE_1[ n + ctcNumNodes*2 ][0] = (*mpModel)( cEz, *ctcNodes_1[n] );
                                                                           
            cE_2[ n                 ][0] = (*mpModel)( cEx, *ctcNodes_2[n] );
			cE_2[ n + ctcNumNodes   ][0] = (*mpModel)( cEy, *ctcNodes_2[n] );
			cE_2[ n + ctcNumNodes*2 ][0] = (*mpModel)( cEz, *ctcNodes_2[n] );
		}

        // Computing surface area
		double na[3], v1[3], v2[3];

        v2[0] = ctcNodes_1[2]->X() - ctcNodes_1[0]->X();
        v2[1] = ctcNodes_1[2]->Y() - ctcNodes_1[0]->Y();
        v2[2] = ctcNodes_1[2]->Z() - ctcNodes_1[0]->Z();

        v1[0] = ctcNodes_1[1]->X() - ctcNodes_1[0]->X();
        v1[1] = ctcNodes_1[1]->Y() - ctcNodes_1[0]->Y();
        v1[2] = ctcNodes_1[1]->Z() - ctcNodes_1[0]->Z();

        // Area = 0.5 * ( v2 x v1 )
        na[0] = v2[1]*v1[2] - v2[2]*v1[1];
        na[1] = v2[2]*v1[0] - v2[0]*v1[2];
        na[2] = v2[0]*v1[1] - v2[1]*v1[0];

		double Area = 0.5 * sqrt( na[0]*na[0] + na[1]*na[1] + na[2]*na[2] );

        // Gauss points and weights
		std::vector<double> cX, cY, W;

		int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		double detJ = 2.00 * Area;

		for( int gp=0; gp<numGaussPoints; gp++ ) W[gp] *= detJ; 

        // Basis functions on Gauss points
		Matrix<double> N;

		Lagrange2D_Ni_3sb( N, cX, cY );

		// Building projection matrices 
		Matrix< std::complex<double> > matrix_NiNj( 9, ctcNumDofs );

		for( int i=0; i<3; i++ )
		{
		    for( int j=0; j<ctcNumNodes; j++ )
		    {
		    	double intg_NiNj = 0.00;
		    
                for( int gp=0; gp<numGaussPoints; gp++ )
				{
					intg_NiNj += W[gp] * N[i][gp] * N[j][gp]; 
                }      
		    	
				matrix_NiNj[ i     ][ j                 ] = intg_NiNj; 
                matrix_NiNj[ i + 3 ][ j + ctcNumNodes   ] = intg_NiNj;  						 
                matrix_NiNj[ i + 6 ][ j + ctcNumNodes*2 ] = intg_NiNj; 
		    }
		}

		// P1 projection matrices
		Matrix< std::complex<double> > iA_div( 3, 3 );
		Matrix< std::complex<double> > iA_cur( 9, 9 );
		
		for( int i=0; i<3; i++ )
		{
            for( int j=0; j<3; j++ )
            {
        		if ( i != j ) iA_div[ i ][ j ] = -3.00 / Area;
        		else          iA_div[ i ][ j ] = +9.00 / Area;                          
				
				iA_cur[ i     ][ j     ] = iA_div[ i ][ j ]; 
                iA_cur[ i + 3 ][ j + 3 ] = iA_div[ i ][ j ];
                iA_cur[ i + 6 ][ j + 6 ] = iA_div[ i ][ j ];
            }									   
		}

		// LL2 projection on nodes 
		Matrix< std::complex<double> > EProj_1( 9, 1 ); 
        Matrix< std::complex<double> > EProj_2( 9, 1 ); 
		
		EProj_1 = ( iA_cur * matrix_NiNj ) * cE_1;
        EProj_2 = ( iA_cur * matrix_NiNj ) * cE_2;

        // Average E field on nodes
		Vector< Vector< std::complex<double> > > avgE_1( 3 );
        Vector< Vector< std::complex<double> > > avgE_2( 3 );

		for( int i=0; i<3; i++ )
		{
            avgE_1[ i ].resize( 3 );
            avgE_1[ i ][ 0 ] = ( 1.0 / 3.0 ) * ( EProj_1[ 0 ][ 0 ] + EProj_1[ 1 ][ 0 ] + EProj_1[ 2 ][ 0 ] ); 
            avgE_1[ i ][ 1 ] = ( 1.0 / 3.0 ) * ( EProj_1[ 3 ][ 0 ] + EProj_1[ 4 ][ 0 ] + EProj_1[ 5 ][ 0 ] ); 
            avgE_1[ i ][ 2 ] = ( 1.0 / 3.0 ) * ( EProj_1[ 6 ][ 0 ] + EProj_1[ 7 ][ 0 ] + EProj_1[ 8 ][ 0 ] ); 

            avgE_2[ i ].resize( 3 );
            avgE_2[ i ][ 0 ] = ( 1.0 / 3.0 ) * ( EProj_2[ 0 ][ 0 ] + EProj_2[ 1 ][ 0 ] + EProj_2[ 2 ][ 0 ] ); 
            avgE_2[ i ][ 1 ] = ( 1.0 / 3.0 ) * ( EProj_2[ 3 ][ 0 ] + EProj_2[ 4 ][ 0 ] + EProj_2[ 5 ][ 0 ] ); 
            avgE_2[ i ][ 2 ] = ( 1.0 / 3.0 ) * ( EProj_2[ 6 ][ 0 ] + EProj_2[ 7 ][ 0 ] + EProj_2[ 8 ][ 0 ] ); 
		}

        // Adding E field values to contact surfaces
        for ( int i=0; i<3; i++ )
        {
            // Surface 1
            if ( mContactEfield.find( HONodesId_1[i] ) != mContactEfield.end() )
            {
                mContactEfield    [ HONodesId_1[i] ] += avgE_1[ i ];
                mContactNeighbours[ HONodesId_1[i] ] += 1;
            }
            else
            {
                mContactEfield    [ HONodesId_1[i] ] = avgE_1[ i ];
                mContactNeighbours[ HONodesId_1[i] ] = 1;
            }

            // Surface 2
            if ( mContactEfield.find( HONodesId_2[i] ) != mContactEfield.end() )
            {
                mContactEfield    [ HONodesId_2[i] ] += avgE_2[ i ];
                mContactNeighbours[ HONodesId_2[i] ] += 1;            
            }
            else
            {
                mContactEfield    [ HONodesId_2[i] ] = avgE_2[ i ];
                mContactNeighbours[ HONodesId_2[i] ] = 1;            
            }
        }
    }

	//**********************************************************************************************************************
    //* - Initiates "volume smoothing" (LL2P formulation)
    //**********************************************************************************************************************
    void Modeler::Ini_Contact_Smoothing_Vol( Vector<int>& CtCNodesId )
    {
        for ( int i=0; i<6; i++ )
        {
            mContactNode[ CtCNodesId[ i ] ] = true;
        }
    }

	//**********************************************************************************************************************
    //* - Calculates "volume smoothing" on discontinuity surfaces (LL2P formulation)
    //**********************************************************************************************************************
    void Modeler::Calculate_Contact_Smoothing_Vol( int* NodesId, unsigned int PropertiesId )
    {
        // Check if element has a node in a contact surface
		bool ThisElementIsInContact = false;

        for ( int i=0; i<4; i++ )
		{
            if ( mContactNode.find( NodesId[i] ) !=  mContactNode.end() )
			{
				ThisElementIsInContact = true;
                break;
			}
		}

        if ( ThisElementIsInContact == false ) return;

        // Getting bubble nodes
        std::vector<int> HONodesId;

		PushHONodesOnVolume( NodesId, HONodesId );

        // Node pointers vector
        std::vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        // Building 3sb element for computing P1 projection on Gauss points
        Element::Pointer pElement;
        
        pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, mpModel->GetProperties(PropertiesId) ) );

        // Electric field P1 projection on central Gauss point
        cVector2Type P1P_Efield;

        pElement->Calculate_E_field_OnGaussPoints( P1P_Efield, 1 );
        
        for ( int i=0; i<4; i++ )
        {
            if ( mContactEfield.find( HONodesId[i] ) != mContactEfield.end() )
            {
                mContactEfield    [ HONodesId[i] ] += P1P_Efield[ 0 ];
                mContactNeighbours[ HONodesId[i] ] += 1;
            }
            else
            {
                mContactEfield    [ HONodesId[i] ] = P1P_Efield[ 0 ];
                mContactNeighbours[ HONodesId[i] ] = 1;
            }
        } 
    }

	//**********************************************************************************************************************
    //* - Updating electric field values on discontinuity surfaces ( LL2P formulation )
    //**********************************************************************************************************************    
    void Modeler::End_Contact_Smooting()
    {
        std::map<unsigned int, Vector<std::complex<double> > >::iterator cef_it;

        for ( cef_it = mContactEfield.begin(); cef_it != mContactEfield.end(); cef_it++ )
        {
            Vector< std::complex<double> > cEfield = cef_it->second;

            int numNeighbours = mContactNeighbours[ cef_it->first ];
        
            cEfield /= numNeighbours;

            Node::Pointer pNode = mpModel->GetNode( ( cef_it->first ) - 1 );

            (*mpModel)( cEx, *pNode ) = cEfield[0];
            (*mpModel)( cEy, *pNode ) = cEfield[1];
            (*mpModel)( cEz, *pNode ) = cEfield[2];

            Vector<double> E_real( 3 );
            Vector<double> E_imag( 3 );

            E_real[0] = std::real( (*mpModel)( cEx, *pNode ) );
            E_real[1] = std::real( (*mpModel)( cEy, *pNode ) );
            E_real[2] = std::real( (*mpModel)( cEz, *pNode ) );
                                                           
            E_imag[0] = std::imag( (*mpModel)( cEx, *pNode ) );
            E_imag[1] = std::imag( (*mpModel)( cEy, *pNode ) );
            E_imag[2] = std::imag( (*mpModel)( cEz, *pNode ) );
                                                           
            (*mpModel)( REAL_E, *pNode ) = E_real;
            (*mpModel)( IMAG_E, *pNode ) = E_imag;
        }

        mContactNode.clear();
        std::map<unsigned int, bool>().swap( mContactNode );

        mContactNeighbours.clear();
        std::map<unsigned int, int>().swap( mContactNeighbours );
	  
        mContactEfield.clear();
        std::map<unsigned int, Vector<std::complex<double> > >().swap( mContactEfield );
    }

	//**********************************************************************************************************************
    //* - Joining normals to second node from the first node
    //**********************************************************************************************************************
	void Modeler::JoinContactNormalsInDirichletS()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		// Unify normal in Dirichlet surface with second node in dielectric pair
		for ( itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc )
        {
			unsigned int FirstNode  = itCtc->first;

			if ( mSetPECnormals.find( FirstNode ) != mSetPECnormals.end() )
			{
				unsigned int SecondNode = ( itCtc->second ).ctcNode;

				std::vector<std::vector<double> > FirstVectorSet ( mSetPECnormals[ FirstNode ] );

				std::vector<std::vector<double> >::iterator it_Fvs;

				for ( it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs )
				{
					Vector<double> vectorF( *it_Fvs );

					mSetPECnormals[ SecondNode ].push_back( vectorF );
				}
			}

            if ( mSetPSBCnormals.find( FirstNode ) != mSetPSBCnormals.end() )
			{
				unsigned int SecondNode = ( itCtc->second ).ctcNode;

				std::vector<std::vector<double> > FirstVectorSet ( mSetPSBCnormals[ FirstNode ] );

				std::vector<std::vector<double> >::iterator it_Fvs;

				for ( it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs )
				{
					Vector<double> vectorF( *it_Fvs );

					mSetPSBCnormals[ SecondNode ].push_back( vectorF );
				}
			}
		}
	}

	//**********************************************************************************************************************
    //* - Force Dirichlet normals in all contact nodes that are in a Dirichlet surface
    //**********************************************************************************************************************
	void Modeler::ForceDirichletNormalsInContact()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		// New dielectric normal. All nodes in a Dirichlet surface have the same normal.
		// PEC/PSBC/PMC/TEPMC dominates normal direction.
        for ( itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc )
        {
			unsigned int FirstNode  =   itCtc->first           ;
			unsigned int SecondNode = ( itCtc->second ).ctcNode;

			if ( ( mNormals.find( SecondNode ) != mNormals.end() ) )
            {
				Vector<double> PECNormal   ( mNormals[ SecondNode ] );
				Vector<double> NEWctcNormal( 3, 0.00 );
				Vector<double> Dummy       ( 3, 0.00 );

				TangencialCoordinates( PECNormal, NEWctcNormal, Dummy );

				mNormals     [ FirstNode ]           = PECNormal   ;
				mContactPairs[ FirstNode ].ctcNormal = NEWctcNormal;
			}
		}
	}

    ////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Fill singularity maps
	////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //**************************************************************************************************
    //* - Fill mNSingular2L map
    //**************************************************************************************************
	void Modeler::FindSing2L( int* NodesId )
	{
		double maxSingValue = 1.00;

		for ( int i=0; i<4; i++ )
        {
            if ( mSingular.find( NodesId[i] ) != mSingular.end() )
            {
				if ( mSingular[ NodesId[i] ] > maxSingValue )
				{
				    maxSingValue = mSingular[ NodesId[i] ];
				}
            }
        }

        if ( maxSingValue > 1.00 )
        {
            for ( int i=0; i<4; i++ )
            {
                if ( mSingular.find( NodesId[i] ) == mSingular.end() )
				{
				    mNSingular2L[ NodesId[i] ] = maxSingValue;
				}
            }
        }
	}

	//**************************************************************************************************
    //* - Fill mNSingular3L map
    //**************************************************************************************************
	void Modeler::FindSing3L( int* NodesId )
	{
		for ( int i=0; i<4; i++ )
        {
            if ( mSingular.find( NodesId[i] ) != mSingular.end() ) return;
        }

		double maxSingValue = 2.00;

		for ( int i=0; i<4; i++ )
        {
            if ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() )
            {
				if ( mNSingular2L[ NodesId[i] ] > maxSingValue )
				{
				    maxSingValue = mNSingular2L[ NodesId[i] ];
				}
            }
        }

        if ( maxSingValue > 2.00 )
        {
            for ( int i=0; i<4; i++ )
            {
                if ( mNSingular2L.find( NodesId[i] ) == mNSingular2L.end() )
				{
				    mNSingular3L[ NodesId[i] ] = maxSingValue;
				}
            }
        }
	}

	//**************************************************************************************************
    //* - Fill mNSingular4L map
    //**************************************************************************************************
	void Modeler::FindSing4L( int* NodesId )
	{
		for( int i=0; i<4; i++ )
        {
            if ( mSingular.find   ( NodesId[i] ) != mSingular.end   () ) return;
			if ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) return;
        }

		double maxSingValue = 3.00;

		for ( int i=0; i<4; i++ )
        {
            if ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() )
            {
				if ( mNSingular3L[ NodesId[i] ] > maxSingValue )
				{
				    maxSingValue = mNSingular3L[ NodesId[i] ];
				}
            }
        }

        if ( maxSingValue > 3.00 )
        {
            for ( int i=0; i<4; i++ )
            {
                if ( mNSingular3L.find( NodesId[i] ) == mNSingular3L.end() )
				{
				    mNSingular4L[ NodesId[i] ] = maxSingValue;
				}
            }
        }
	}

	//**************************************************************************************************
    //* - Fill mNSingular5L map
    //**************************************************************************************************
	void Modeler::FindSing5L( int* NodesId )
	{
		for ( int i=0; i<4; i++ )
        {
            if ( mSingular.find   ( NodesId[i] ) != mSingular.end   () ) return;
			if ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) return;
			if ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() ) return;
        }

		double maxSingValue = 4.00;

		for ( int i=0; i<4; i++ )
        {
            if ( mNSingular4L.find( NodesId[i] ) != mNSingular4L.end() )
            {
				if ( mNSingular4L[ NodesId[i] ] > maxSingValue )
				{
				    maxSingValue = mNSingular4L[ NodesId[i] ];
				}
            }
        }

        if ( maxSingValue > 4.00 )
        {
            for ( int i=0; i<4; i++ )
            {
                if ( mNSingular4L.find( NodesId[i] ) == mNSingular4L.end() )
				{
				    mNSingular5L[ NodesId[i] ] = maxSingValue;
				}
            }
        }
	}

	//**************************************************************************************************
    //* - Fill mNSingular6L map
    //**************************************************************************************************
	void Modeler::FindSing6L( int* NodesId )
	{
		for ( int i=0; i<4; i++ )
        {
            if ( mSingular.find   ( NodesId[i] ) != mSingular.end   () ) return;
			if ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) return;
			if ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() ) return;
			if ( mNSingular4L.find( NodesId[i] ) != mNSingular4L.end() ) return;
        }

		double maxSingValue = 5.00;

		for ( int i=0; i<4; i++ )
        {
            if ( mNSingular5L.find( NodesId[i] ) != mNSingular5L.end() )
            {
				if (mNSingular5L[ NodesId[i] ] > maxSingValue)
				{
				    maxSingValue = mNSingular5L[ NodesId[i] ];
				}
            }
        }

        if ( maxSingValue > 5.00 )
        {
            for ( int i=0; i<4; i++ )
            {
                if ( mNSingular5L.find( NodesId[i] ) == mNSingular5L.end() )
				{
				    mNSingular6L[ NodesId[i] ] = maxSingValue;
				}
            }
        }
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////
	// - High Order Elements
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//**************************************************************************************************
    //* - Add high order elements to the HONodesId vector in volume elements
    //**************************************************************************************************
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

		if( IdNode1 < IdNode2 )
		{
		    if( IdNode2 < IdNode3 )
		    {
			    //[n1,n2,n3]
		        orderingCase = 1;
			    minId = IdNode1;
			    medId = IdNode2;
			    maxId = IdNode3;
		    }
		    else
		    {
		        if( IdNode1 < IdNode3 )
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
		    if( IdNode1 < IdNode3 )
		    {
			    //[n2,n1,n3]
		        orderingCase = 4;
			    minId = IdNode2;
			    medId = IdNode1;
			    maxId = IdNode3;
		    }
		    else
		    {
			    if( IdNode2 < IdNode3 )
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

		std::vector<std::vector<unsigned int> > SameKeyFaces( mFaceHONodes[minId*medId*maxId] );
        std::vector<std::vector<unsigned int> >::iterator it_faces;

		for( it_faces = SameKeyFaces.begin(); it_faces != SameKeyFaces.end(); ++it_faces )
        {
            std::vector<unsigned int> face(*it_faces);

            if ( ( face[0]==minId ) && ( face[1]==medId ) && ( face[2]==maxId ) )
            {
			    if      ( mElementOrder == 0 ) ReorderHONodesOnFace_4th( orderingCase, face, HONodesId );
                else if ( mElementOrder == 3 ) HONodesId.push_back     ( face[3]                       );
			    else if ( mElementOrder == 4 ) ReorderHONodesOnFace_4th( orderingCase, face, HONodesId );
            }
        }
	}

	//*************************************************************************************************
    //* - Reordering nodes for a giving surface order
    //*************************************************************************************************
	void Modeler::ReorderHONodesOnFace_4th( int orderingCase, std::vector<unsigned int>& face, std::vector<int>& HONodesId )
	{
		unsigned int minf = face[3];
		unsigned int medf = face[4];
		unsigned int maxf = face[5];
        
		//[n1,n2,n3] = [minId, medId, maxId]
		if( orderingCase == 1 )
		{
		    HONodesId.push_back(minf);
		    HONodesId.push_back(medf);
		    HONodesId.push_back(maxf);
		}
		//[n1,n3,n2] = [minId, maxId, medId]
		else if( orderingCase == 2 )
		{
		    HONodesId.push_back(minf);
		    HONodesId.push_back(maxf);
		    HONodesId.push_back(medf);
		}
		//[n3,n1,n2] = [medId, maxId, minId]
		else if( orderingCase == 3 )
		{
		    HONodesId.push_back(medf);
		    HONodesId.push_back(maxf);
		    HONodesId.push_back(minf);
		}
		//[n2,n1,n3] = [medId, minId, maxId]
		else if( orderingCase == 4 )
		{
		    HONodesId.push_back(medf);
		    HONodesId.push_back(minf);
		    HONodesId.push_back(maxf);
		}
		//[n2,n3,n1] = [maxId, minId, medId]
		else if( orderingCase == 5 )
		{
		    HONodesId.push_back(maxf);
		    HONodesId.push_back(minf);
		    HONodesId.push_back(medf);
		}
		//[n3,n2,n1] = [maxId, medId, minId]
		else if( orderingCase == 6 )
		{
		    HONodesId.push_back(maxf);
		    HONodesId.push_back(medf);
		    HONodesId.push_back(minf);
	    }
	}

	//*************************************************************************************************
    //* - Looking for high order nodes on volume
    //*************************************************************************************************
	void Modeler::HONodesOnVolume( int* NodesId, std::vector<int>& HONodesId )
	{
	    unsigned int key = NodesId[0]*NodesId[1]*NodesId[2]*NodesId[3];
        
        std::vector<std::vector<unsigned int> > SameKeyElement( mVolmHONodes[key] );

        std::vector<std::vector<unsigned int> >::iterator it_elements;
        
        for( it_elements = SameKeyElement.begin(); it_elements != SameKeyElement.end(); ++it_elements )
        {
            std::vector<unsigned int> nInElement( *it_elements );
        
            if ( ( nInElement[0] == NodesId[0] ) && ( nInElement[1] == NodesId[1] ) && 
                 ( nInElement[2] == NodesId[2] ) && ( nInElement[3] == NodesId[3] ) )
            {
                HONodesId.push_back( nInElement[4] );
            }
        }
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	// - High order nodes creation
	///////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************
    //* - Create high order nodes
    //*************************************************************************************************
	void Modeler::CreateHONodes( int* NodesId )
	{
        if( mElementOrder == 0 ) 
        {
            CreateHONodes_OnFaces ( NodesId );
            CreateHONodes_OnVolume( NodesId ); 
        }

        if( mElementOrder <= 1 ) return; 

		CreateHONodes_OnEdges( NodesId );

		if( mElementOrder == 2 ) return; 
       
		CreateHONodes_OnFaces( NodesId );

		if( mElementOrder == 3 ) return; 
      
		CreateHONodes_OnVolume( NodesId );
	}

    //*************************************************************************************************
    //* - Create high order nodes if quadratic meshing (only for 2nd order elements)
    //*************************************************************************************************
	void Modeler::GenerateVEQElement( std::vector<int>& QHONodesId )
	{
		int i, j;

		unsigned int IdNode1,IdNode2;
        unsigned int key;

		unsigned int middleNode[6] = { QHONodesId[4], QHONodesId[6], QHONodesId[7], QHONodesId[5], QHONodesId[8], QHONodesId[9] };

		// If PropertyId == 0 throw error, material not assigned
		if( QHONodesId[10] == 0 ) { throw 0; }

		// Save material properties for contact elements
		for( i=0; i<4; i++ ) mTempNodeProperties[QHONodesId[i]] = QHONodesId[10];

        int midNodeIndex = 0;

		for( i=0; i<3; i++ )
		{
			for( j=i+1; j<4; j++ )
			{
				if( QHONodesId[i] < QHONodesId[j] )
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

				if( mEdgeHONodes.find(key) != mEdgeHONodes.end() )
				{
					std::vector<std::vector<unsigned int> > ExistantNodesOnEdges(mEdgeHONodes[key]);
					std::vector<std::vector<unsigned int> >::iterator it_ee;

					int repe = 0;

					for( it_ee = ExistantNodesOnEdges.begin(); it_ee != ExistantNodesOnEdges.end(); ++it_ee )
					{
						std::vector<unsigned int> NodesOnEdge(*it_ee);

						if ( (IdNode1 == NodesOnEdge[0]) && (IdNode2 == NodesOnEdge[1]) ) repe++;
					}

					if( repe == 0 )
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
	void Modeler::CreateHONodes_OnEdges( int* NodesId )
	{
		unsigned int IdNode1,IdNode2;
        unsigned int key;

		for( int i=0; i<3; i++ )
		{
			for( int j=i+1; j<4; j++ )
			{
				if( NodesId[i] < NodesId[j] )
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

				if( mEdgeHONodes.find(key) != mEdgeHONodes.end() )
				{
					std::vector<std::vector<unsigned int> > ExistantNodesOnEdges(mEdgeHONodes[key]);
					std::vector<std::vector<unsigned int> >::iterator it_ee;

					int repe = 0;

					for( it_ee = ExistantNodesOnEdges.begin(); it_ee != ExistantNodesOnEdges.end(); ++it_ee )
					{
						std::vector<unsigned int> NodesOnEdge(*it_ee);

						if ( (IdNode1 == NodesOnEdge[0]) && (IdNode2 == NodesOnEdge[1]) ) repe++;
					}

					if( repe == 0 )
					{
						if      ( mElementOrder == 2 ) CreateHONodes_OnEdge_2nd( IdNode1, IdNode2 );
						else if ( mElementOrder == 3 ) CreateHONodes_OnEdge_3th( IdNode1, IdNode2 );
						else if ( mElementOrder == 4 ) CreateHONodes_OnEdge_4th( IdNode1, IdNode2 );
					}
				}
				else
				{
					if      ( mElementOrder == 2 ) CreateHONodes_OnEdge_2nd( IdNode1, IdNode2 );
					else if ( mElementOrder == 3 ) CreateHONodes_OnEdge_3th( IdNode1, IdNode2 );
					else if ( mElementOrder == 4 ) CreateHONodes_OnEdge_4th( IdNode1, IdNode2 );
				}
			}
		}
	}

	//*************************************************************************************************
    //* - Creates high order nodes on faces
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnFaces( int* NodesId )
	{
		unsigned int IdNode1,IdNode2,IdNode3;
        unsigned int key;

		for( int i=0; i<2; i++ )
		{
			for( int j=i+1; j<3; j++ )
			{
				for( int m=j+1; m<4; m++ )
				{
					if( NodesId[i] < NodesId[j] )
					{
						if( NodesId[j] < NodesId[m] )
						{
							IdNode1 = NodesId[i];
						    IdNode2 = NodesId[j];
							IdNode3 = NodesId[m];
						}
						else
						{
							if( NodesId[i] < NodesId[m] )
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
						if( NodesId[i] < NodesId[m] )
						{
							IdNode1 = NodesId[j];
						    IdNode2 = NodesId[i];
							IdNode3 = NodesId[m];
						}
						else
						{
							if( NodesId[j] < NodesId[m] )
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

					if( mFaceHONodes.find(key) != mFaceHONodes.end() )
					{
						std::vector<std::vector<unsigned int> > ExistantNodesOnFaces( mFaceHONodes[key] );
						std::vector<std::vector<unsigned int> >::iterator it_ef;

						int repe = 0;

						for( it_ef = ExistantNodesOnFaces.begin(); it_ef != ExistantNodesOnFaces.end(); ++it_ef )
						{
							std::vector<unsigned int> NodesOnFace(*it_ef);

							if( ( IdNode1==NodesOnFace[0] ) && ( IdNode2==NodesOnFace[1] ) && ( IdNode3==NodesOnFace[2] ) ) repe++;
						}

						if( repe == 0 )
						{
							if      ( mElementOrder == 0 ) CreateHONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
                            else if ( mElementOrder == 3 ) CreateHONodes_OnFace    ( IdNode1, IdNode2, IdNode3 );
							else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
						}
					}
					else
					{
					    if      ( mElementOrder == 0 ) CreateHONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
                        else if ( mElementOrder == 3 ) CreateHONodes_OnFace    ( IdNode1, IdNode2, IdNode3 );
					    else if ( mElementOrder == 4 ) CreateHONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
					}
				}
			}
		}
	}

	//*************************************************************************************************
    //* - Creating new higher order nodes on edge (2nd order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_2nd( int IdNode1, int IdNode2 )
    {
		std::vector<unsigned int> edgeNodes(3);

        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2-1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 0.5;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[2], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node1 );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );
	}

	//*************************************************************************************************
    //* - Creating new higher order nodes on edge (3th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_3th( int IdNode1, int IdNode2 )
    {
		std::vector<unsigned int> edgeNodes(4);

        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;
        edgeNodes[3] = edgeNodes[2]+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2-1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 1.00/3.00;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[2], hoX, hoY, hoZ ) );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );

		//Creating new high order node on edge number 2
        alpha = 2.00/3.00;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( edgeNodes[3], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node2 );

        if ( mPotentials_On ) new_ho_node2->Add_cA_Dofs( mDofSet );
        else                  new_ho_node2->Add_cE_Dofs( mDofSet );
	}

	//*************************************************************************************************
    //* - Creating a new higher order node on face 
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnFace( int IdNode1, int IdNode2, int IdNode3 )
    {
		std::vector<unsigned int> faceNodes(4);

        faceNodes[0] = IdNode1;
        faceNodes[1] = IdNode2;
        faceNodes[2] = IdNode3;

        faceNodes[3] = (mpModel->GetNodesArray()).size()+1;

        mFaceHONodes[IdNode1*IdNode2*IdNode3].push_back( faceNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNode1-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2-1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNode3-1 );

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

        Node::Pointer new_ho_node1( new Node( faceNodes[3], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node1 );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );
	}

	//*************************************************************************************************
    //* - Creating new higher order nodes on edge (4th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnEdge_4th( int IdNode1, int IdNode2 )
    {
        std::vector<unsigned int> edgeNodes(5);

        edgeNodes[0] = IdNode1;
        edgeNodes[1] = IdNode2;

        edgeNodes[2] = (mpModel->GetNodesArray()).size()+1;
        edgeNodes[3] = edgeNodes[2]+1;
        edgeNodes[4] = edgeNodes[3]+1;

        mEdgeHONodes[IdNode1*IdNode2].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2-1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        //Creating new high order node on edge number 1
        alpha = 0.25;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[2], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node1 );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );

		//Creating new high order node on edge number 2
        alpha = 0.50;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( edgeNodes[3], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node2 );

        if ( mPotentials_On ) new_ho_node2->Add_cA_Dofs( mDofSet );
        else                  new_ho_node2->Add_cE_Dofs( mDofSet );

		//Creating new high order node on edge number 3
        alpha = 0.75;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node3( new Node( edgeNodes[4], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node3 );

        if ( mPotentials_On ) new_ho_node3->Add_cA_Dofs( mDofSet );
        else                  new_ho_node3->Add_cE_Dofs( mDofSet );
    }

	//*************************************************************************************************
    //* - Creating new higher order nodes on face (4th order)
    //*************************************************************************************************
    void Modeler::CreateHONodes_OnFace_4th( int IdNode1, int IdNode2, int IdNode3 )
    {
        std::vector<unsigned int> faceNodes(6);

        faceNodes[0] = IdNode1;
        faceNodes[1] = IdNode2;
        faceNodes[2] = IdNode3;

        faceNodes[3] = (mpModel->GetNodesArray()).size()+1;
        faceNodes[4] = faceNodes[3]+1;
        faceNodes[5] = faceNodes[4]+1;

        mFaceHONodes[IdNode1*IdNode2*IdNode3].push_back( faceNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNode1-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2-1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNode3-1 );

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

        Node::Pointer new_ho_node1( new Node( faceNodes[3], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node1 );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );

		//Creating new high order node on face number 2
        alpha = 0.50;
        beta  = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( faceNodes[4], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node2 );

        if ( mPotentials_On ) new_ho_node2->Add_cA_Dofs( mDofSet );
        else                  new_ho_node2->Add_cE_Dofs( mDofSet );

		//Creating new high order node on face number 3
        alpha = 0.25;
        beta  = 0.50;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node3( new Node( faceNodes[5], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node3 );

        if ( mPotentials_On ) new_ho_node3->Add_cA_Dofs( mDofSet );
        else                  new_ho_node3->Add_cE_Dofs( mDofSet );
    }

	//*************************************************************************************************
    //* - Creating a new higher order node inside volume 
    //*************************************************************************************************
	void Modeler::CreateHONodes_OnVolume( int* IdNodes )
    {
		std::vector<unsigned int> volumNodes(5);

        volumNodes[0] = IdNodes[0];
        volumNodes[1] = IdNodes[1];
        volumNodes[2] = IdNodes[2];
		volumNodes[3] = IdNodes[3];

        volumNodes[4] = (mpModel->GetNodesArray()).size()+1;

        mVolmHONodes[IdNodes[0]*IdNodes[1]*IdNodes[2]*IdNodes[3]].push_back( volumNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNodes[0]-1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNodes[1]-1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNodes[2]-1 );
		Node::Pointer pNode4 = mpModel->GetNode( IdNodes[3]-1 );

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

        Node::Pointer new_ho_node1( new Node( volumNodes[4], hoX, hoY, hoZ ) );

        mpModel->AddNode( new_ho_node1 );

        if ( mPotentials_On ) new_ho_node1->Add_cA_Dofs( mDofSet );
        else                  new_ho_node1->Add_cE_Dofs( mDofSet );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Building
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*****************************************************************************************************************
    //* - Checks if the problem set-up is consistent
    //*****************************************************************************************************************
    void Modeler::CheckConsistency()
    {
        mQuadraticGeometry = false;
        
        // Full wave mode
        if ( mFullWaveMode == true )
        {
            if ( mProblemFrequency == 0.0 ) 
            {
                Send_Error_Msg( "Problem frequency equal to 0.0", "Please, change problem frequency.", 1 );
            }
        }

        // Cold plasma mode
        if ( mColdPlasmaMode == true )
        {
            mPotentials_On = false;

            if ( mElementOrder == 0 ) 
            {   
                Send_Error_Msg( "LL2P-3sb elements not implemented for cold plasma mode", "Please, change element type or problem mode.", 1 );
            }

            if ( mProblemFrequency == 0.0 ) 
            {
                Send_Error_Msg( "Problem frequency equal to 0.0", "Please, change problem frequency.", 1 );
            }
        }

        // Electrostatic mode 
        if ( mElectrostaticMode == true )  
        {
            mProblemFrequency  = 0.0  ;
            mPotentials_On     = true ;
            mFrequencySweep    = false;
            
            mImportCurrents    = false;

            mQuadraticGeometry = false;
            mElementOrder      = 1    ;
        }

        // Quadratic geometry only applies to 2nd order elements 
        if ( ( mQuadraticGeometry == true ) && ( mElementOrder != 2 ) ) 
        {
            mQuadraticGeometry = false;
        }
    }

    //*****************************************************************************************************************
    //* - Assembling volume elements in electro-static mode
    //*****************************************************************************************************************
    void Modeler::GenerateVolumeElement_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertyId == 0 throw error, material not assigned
        if ( PropertiesId == 0 ) 
        { 
            throw 0; 
        }

        // Get local nodes
        std::vector<Node::Pointer> pNodes   ;
        std::vector<int>           HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get material properties
        Properties::Pointer Properties = mpModel->GetProperties( PropertiesId );

        // Create element
        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, Properties ) );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector
        Vector<double> EleResVector;

        pElement->GetResidualVectorDirichlet( mFixVoltage, EleStiffMatrix, EleResVector );

        // Assemble elemental stiffness matrix to symmetric global matrix
        Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );

        // Assemble elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Clean objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

    //*****************************************************************************************************************
    //* - Sets element parameters ( singularity weight, normal type, LL2P factor )
    //*****************************************************************************************************************
    void Modeler::Set_Volume_Element_Parameters( int* NodesId, Element::Pointer pElement )
    {
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

        // Setting regularization weight 
		if ( IsUgLayer ) 
        {
            pElement->SetPeso( 0.00 );
        }
		else       
        {
            pElement->SetPeso( 1.00 );
        }

        // Setting parameters for LL2P bubble elements
        if ( mElementOrder == 0 ) 
        {
            pElement->SetNormalType( mType_Of_BC_Normal );
            pElement->SetFactor    ( mLL2P_hk_factor    );
        }

        // Setting plasma model
        if ( mColdPlasmaMode )        
        {
            pElement->SetPlasmaModel( mpColdPlasma );
        }

        // Activates potentials on element
        pElement->SetPotentials( mPotentials_On );
    }

    //**********************************************************************************************************************************************************
    //* - Applies boundary conditions to element stiffness matrix
    //**********************************************************************************************************************************************************
    void Modeler::Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        // PBC conditions ( Trans(T)*A*T )
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Contact discontinuities ( Trans(T)*A*T )
        if ( mColdPlasmaMode == true ) 
        {
            Apply_Contact_Matrix_TKT_ColdPlasma( HONodesId, EleIdVector, EleStiffMatrix );
        }
        else  
        {
            Apply_Contact_Matrix_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        }
       
        // PEC, PMC conditions ( Trans(T)*A*T )
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

        // Axisymmetric mode ( Trans(T)*A*T ) 
		Axisym_TKT( HONodesId, EleStiffMatrix );

        // Cold palsma RLP base ( Inv(T)*A*T )
        if ( mColdPlasmaMode == true ) 
        {
            ColdPlasma_TKT( HONodesId, EleStiffMatrix );
        }
    }

    //**********************************************************************************************************************************************************
    //* - Applies boundary conditions to element force vector
    //**********************************************************************************************************************************************************
    void Modeler::Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Vector< std::complex<double> >& EleResVector )
    {
	    // PBC conditions ( Inv(T)*b ) 
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );
	    PBC_RightLeft_TR( HONodesId, EleIdVector, EleResVector );
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );

		// Contact discontinuities ( Inv(T)*b )
        if ( mColdPlasmaMode == true ) 
        {
            Apply_Contact_Vector_TR_ColdPlasma( HONodesId, EleIdVector, EleResVector );
        }
        else   
        {
            Apply_Contact_Vector_TR( HONodesId, EleIdVector, EleResVector );
        }

        // PEC, PMC conditions ( Inv(T)*B ) 
        Apply_Rotation_To_Vector( HONodesId, EleResVector );

		// Axisymmetric mode ( Inv(T)*B )
        if ( mAxisymmetric == true ) 
        {
            Axisym_TR( HONodesId, EleResVector );
        }

        // Cold plasma RLP base ( Inv(T)*B )
        if ( mColdPlasmaMode == true ) 
        {
            ColdPlasma_TR( HONodesId, EleResVector );
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembling elemental stiffness matrix (real) to global symmetric matrix
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Symmetric( Vector<int>& EleIdVector, Matrix<double>& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

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
    }

    //*****************************************************************************************************************************************
    //* - Assembling elemental stiffness matrix (complex) to global symmetric matrix
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Symmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if ( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;

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

    //*****************************************************************************************************************************************
    //* - Assembling elemental stiffness matrix (complex) to global non-symmetric matrix
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_NonSymmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if ( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;

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

    //*****************************************************************************************************************************************
    //* - Assembling elemental stiffness matrix (complex) to auxiliar symmetric global matrix
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Aux_Symmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if ( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;        
        
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

    //*****************************************************************************************************************************************
    //* - Assembling elemental stiffness matrix (complex) to auxiliar non-symmetric global matrix
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if ( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;        

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

    //*****************************************************************************************************************************************
    //* - Assembling elemental residual vector (complex) to global residual vector
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ResidualVector_to_Global( Vector<int>& EleIdsVector, Vector< std::complex<double> >& EleResVector )
    {
        int eleSize = EleIdsVector.size();
        
        for ( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdsVector[ i_local ];

            if ( i_global < mSystemSize )
            {
                b_vector[ i_global ] += EleResVector[ i_local ];
            }
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembling elemental residual vector (real) to global residual vector
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ResidualVector_to_Global( Vector<int>& EleIdsVector, Vector<double>& EleResVector )
    {
        int eleSize = EleIdsVector.size();
        
        for ( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdsVector[ i_local ];

            if ( i_global < mSystemSize )
            {
                b_vector[ i_global ] += EleResVector[ i_local ];
            }
        }
    }

    //*************************************************************************************************************************
    //* - Assembling volume elements in cold plasma mode 
    //*************************************************************************************************************************
    void Modeler::GenerateVolumeElement_ColdPlasma( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertyId == 0 throw error, material not assigned
        if ( PropertiesId == 0 ) 
        { 
            throw 0; 
        }

        // Get volume element nodes
		Vector<Node::Pointer> pNodes;

        Vector<int> HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer Properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Plasma element definition
        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_ColdPlasma ( pNodes, Properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_ColdPlasma ( pNodes, Properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2ndQ_ColdPlasma( pNodes, Properties ) );

        // Set element parameters
        Set_Volume_Element_Parameters( NodesId, pElement );

        // Get DOFs Id vector
        Vector<int> EleIdsVector; 

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix; 

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Apply boundary conditions
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Assemble elemental stiffness matrix to global matrix 
        if      ( mpColdPlasma->Is_Full_Matrix    () ) Assemble_ElementMatrix_to_Global_NonSymmetric    ( EleIdsVector, EleStiffMatrix );
        else if (    !pElement->GetIsIHL          () ) Assemble_ElementMatrix_to_Global_Symmetric       ( EleIdsVector, EleStiffMatrix );
        else if ( mpColdPlasma->Is_HermSymm_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_Symmetric   ( EleIdsVector, EleStiffMatrix );           
        else if ( mpColdPlasma->Is_HermFull_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( EleIdsVector, EleStiffMatrix );
     
        // Clean objects
        EleIdsVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*****************************************************************************************************************
    //* - Assembling volume elements in Full-Wave mode 
    //*****************************************************************************************************************
    void Modeler::GenerateVolumeElement_FullWave( int* NodesId, unsigned int PropertiesId )
    {
		// If PropertyId == 0 then throw error ( material not assigned )
		if ( PropertiesId == 0 ) 
        { 
            throw 0; 
        }

        // Get volume element nodes
		Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, nodes );

        // Get element properties
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId ); 
        
        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Volume element definition
        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

        // Set element parameters
        Set_Volume_Element_Parameters( NodesId, pElement );

        // Get DOFs Id vector
        Vector<int> EleIdsVector; 

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix; 

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Apply boundary conditions
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Assemble elemental stiffness matrix to symmetric global matrix
        Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );

        // Clean objects
        EleIdsVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembling high order source elements
    //*************************************************************************************************************************************
    void Modeler::GenerateSourceElement( int* NodesId, unsigned int PropertiesId )
    {
        // Get volume element nodes
		Vector<Node::Pointer> pNodes;

        Vector<int> HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer Properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Volumetric source element definition
        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, Properties ) );
	    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, Properties ) );
	    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, Properties ) );
	    
        // Activates potentials on element
        pElement->SetPotentials( mPotentials_On );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental residual vector
        Vector<std::complex<double> > EleResVector; 

        pElement->GetResidualVector( EleResVector );

        // Apply boundary conditions to force vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Cleaning vectors 
        EleIdsVector.FreeData();
        EleResVector.FreeData();
    }

    //*************************************************************************************************************************************
    //* - Assembling boundary plasma element matrix on auxiliar global matrix  
    //*************************************************************************************************************************************
    void Modeler::Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();
         
        std::complex<double> cKij;
        
        double Kij_Tol = mpColdPlasma->Get_Kij_Tol();
        
        if ( mpColdPlasma->Is_HermSymm_Matrix() )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];
        
                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];
        
                        if ( (i_global <= j_global) && (j_global < mSystemSize) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];

                            if ( abs(cKij) > Kij_Tol ) A_matrix_aux( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else if ( mpColdPlasma->Is_HermFull_Matrix() )
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];
        
                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];
        
                        if ( (j_global < mSystemSize) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];

                            if ( abs(cKij) > Kij_Tol ) A_matrix_aux( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else
        {
            for ( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[i_local];
        
                if ( i_global < mSystemSize )
                {
                    for ( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[j_local];
        
                        if ( (j_global < mSystemSize) )
                        {
                            cKij = EleStiffMatrix[i_local][j_local];

                            if ( abs(cKij) > Kij_Tol ) A_matrix( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
    }

    //*************************************************************************************************************************************
    //* - Assembling generic Robin elements in plasma mode
    //*************************************************************************************************************************************
    void Modeler::GenerateGenericRobinElement_ColdPlasma( int* NodesId, unsigned int PropertiesId )
    {
        Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

        PushHONodesOnSurface(NodesId,HONodesId);

        std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) nodes.push_back(mpModel->GetNode((*hoit)-1));

        Properties::Pointer properties = mpModel->GetProperties(PropertiesId);

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new GenericRobin_3sb_ColdPlasma ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new GenericRobin_1st_ColdPlasma ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new GenericRobin_2ndQ_ColdPlasma( nodes, properties ) );

        pElement->SetPlasmaModel( mpColdPlasma      );
        pElement->SetFrequency  ( mProblemFrequency );

        Vector<int>                   EleIdVector   ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );

        //Trans(T)*A*T for PBC conditions
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        //Trans(T)*A*T for discontinuities
        Apply_Contact_Matrix_TKT_ColdPlasma( HONodesId, EleIdVector, EleStiffMatrix );

        //Trans(T)*A*T for PEC and PMC
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

        //Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        ColdPlasma_TKT( HONodesId, EleStiffMatrix );

        //Trans(T)*A*T for axisymmetric problems
        Axisym_TKT( HONodesId, EleStiffMatrix );

        //Assembling to global matrix
        Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );

        //Cleaning elemental matrix
        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

    //*************************************************************************************************************************************
    //* - Assembling far field elements in plasma mode
    //*************************************************************************************************************************************
    void Modeler::GenerateFarFieldElement_ColdPlasma( int* NodesId, unsigned int PropertiesId )
    {
        Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

        PushHONodesOnSurface( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit ) nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new FarField_3sb_ColdPlasma ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new FarField_1st_ColdPlasma ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new FarField_2ndQ_ColdPlasma( nodes, properties ) );

        pElement->SetPlasmaModel( mpColdPlasma      );
        pElement->SetFrequency  ( mProblemFrequency );
 
        Vector<int>                   EleIdVector   ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );

        //Trans(T)*A*T for PBC conditions
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        //Trans(T)*A*T for discontinuities
        Apply_Contact_Matrix_TKT_ColdPlasma( HONodesId, EleIdVector, EleStiffMatrix );

        //Trans(T)*A*T for PEC and PMC
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

        //Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        ColdPlasma_TKT( HONodesId, EleStiffMatrix );

        //Trans(T)*A*T for axisymmetric problems
        Axisym_TKT( HONodesId, EleStiffMatrix );

        //Assembling to global matrix
        Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );

        // Cleaning elemental matrix
        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembling far field elements
    //*************************************************************************************************************************************
    void Modeler::GenerateFarFieldElement( int* NodesId )
    {
        Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

		PushHONodesOnSurface(NodesId,HONodesId);

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit) 
        {
            nodes.push_back(mpModel->GetNode((*hoit)-1));
        }

        Element::Pointer pElement;
 
        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new FarField_3sb_FullWave ( nodes ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new FarField_1st_FullWave ( nodes ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new FarField_2ndQ_FullWave( nodes ) );

		pElement->SetFrequency( mProblemFrequency );

		Vector<int>                   EleIdVector   ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );

        // Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
	    PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Trans(T)*A*T for contact discontinuities
        if ( mColdPlasmaMode == true ) 
        {
            Apply_Contact_Matrix_TKT_ColdPlasma( HONodesId, EleIdVector, EleStiffMatrix );
        }
        else   
        {
            Apply_Contact_Matrix_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        }
        
		// Trans(T)*A*T
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

        //Trans(T)*A*T for axisymmetric problems
        if ( mAxisymmetric   == true ) Axisym_TKT( HONodesId, EleStiffMatrix );

		//Inv(T)*A*T to rotate matrix to the diagonal permittivity tensor base (RLP)
        if ( mColdPlasmaMode == true ) ColdPlasma_TKT( HONodesId, EleStiffMatrix );

        //Assembling to global matrix 
        if ( mColdPlasmaMode == true )
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
             int eleSize = EleIdVector.size();
             
             std::complex<double> cKij;
             
             //Assembling to symmetric matrix
             for ( int i_local = 0; i_local < eleSize; i_local++ )
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
    void Modeler::GenerateGenericRobinElement_FullWave( int* NodesId, unsigned int PropertiesId )
    {
        Vector<Node::Pointer> nodes;
        
        Vector<int> HONodesId;

		PushHONodesOnSurface( NodesId, HONodesId );

		std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )  
        {
            nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement;
		
        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new GenericRobin_3sb_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new GenericRobin_1st_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new GenericRobin_2ndQ_FullWave( nodes, properties ) );

		Vector<int>                   EleIdVector   ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );

		// Trans(T)*A*T for PBC conditions
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
	    PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
		PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Trans(T)*A*T for contact discontinuities
        Apply_Contact_Matrix_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Trans(T)*A*T
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

		// Trans(T)*A*T for axisymmetric problems
		Axisym_TKT( HONodesId,EleStiffMatrix );

        //Assembling to global matrix
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

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

                        if ( abs(cKij) > 0.0 ) A_matrix( i_global, j_global ) += cKij;
                    }
				}
			}
        }

        EleIdVector.FreeData();
        EleStiffMatrix.FreeData();
	}

    //*************************************************************************************************
    //* - Assembling generic Robin elements (Electrostatic mode)
    //*************************************************************************************************
    void Modeler::GenerateGenericRobinElement_Electrostatic( int* NodesId, unsigned int PropertiesId )
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

        pElement = Element::Pointer( new GenericRobin_1st_Electrostatic( nodes, properties ) );

        // Getting DOFs Id vector
        Vector<int> EleIdVector;

        pElement->GetEquationIdVector( EleIdVector );

        // Getting elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

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

        EleIdVector.FreeData();
        EleResVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//********************************************************************************************************
    //* - Assembling boundary elements
    //********************************************************************************************************
    void Modeler::GenerateRWPortTE10Element( int* NodesId, unsigned int PropertiesId )
    {
		Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

		PushHONodesOnSurface( NodesId, HONodesId );

		std::vector<int>::iterator hoit;

        for(hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit)  
        {
            nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new RWPortTE10_3sb_FullWave ( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new RWPortTE10_1st_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new RWPortTE10_2ndQ_FullWave( nodes, properties ) );

		Vector<int>                   EleIdVector   ;
		Vector<std::complex<double> > EleResVector  ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );
        pElement->GetResidualVector  ( EleResVector   );

        // Trans(T)*A*T
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

		// Inv(T)*b
        Apply_Rotation_To_Vector( HONodesId, EleResVector );

        // Assembling to global matrix
        int eleSize = EleIdVector.size();

        if ( mColdPlasmaMode == true )
        {
             Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
             std::complex<double> cKij;
             
             for ( int i_Klocal = 0; i_Klocal < eleSize; i_Klocal++ )
             {
                  int i_Kglobal = EleIdVector[i_Klocal];
                  
                  if ( i_Kglobal < mSystemSize )
                  {
                       for ( int j_Klocal = 0; j_Klocal < eleSize; j_Klocal++ )
                       {
                       	    int j_Kglobal = EleIdVector[j_Klocal];
                            
                       	    if ( ( i_Kglobal <= j_Kglobal ) && ( j_Kglobal < mSystemSize ) )
                            {
                                 cKij = EleStiffMatrix[i_Klocal][j_Klocal];

                                 if ( abs(cKij) > 0.0 ) A_matrix( i_Kglobal, j_Kglobal ) += cKij;
                            }
                       }    
                  }
             }
        }

        //Assembling to global source vector
		for ( int i_blocal=0; i_blocal < eleSize; i_blocal++ )
		{
		     int i_bglobal = EleIdVector[i_blocal];
             
		     if ( i_bglobal < mSystemSize )
		     {
		          b_vector[i_bglobal] += EleResVector[i_blocal];
		     }
		}

        //Cleaning element matrix and vectors
        EleIdVector.FreeData();
        EleResVector.FreeData();
        EleStiffMatrix.FreeData();
    }

	//***********************************************************************************************************
    //* - Assembling boundary elements
    //***********************************************************************************************************
    void Modeler::GenerateCoaxPortTEMElement( int* NodesId, unsigned int PropertiesId )
    {
		Vector<Node::Pointer> nodes;

        Vector<int> HONodesId;

		PushHONodesOnSurface( NodesId, HONodesId );

		std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )  
        {
            nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new CoaxialPortTEM_3sb_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new CoaxialPortTEM_2ndQ_FullWave( nodes, properties ) );

		Vector<int>                   EleIdVector   ;
		Vector<std::complex<double> > EleResVector  ;
        Matrix<std::complex<double> > EleStiffMatrix;

        pElement->GetEquationIdVector( EleIdVector    );
        pElement->GetStiffnessMatrix ( EleStiffMatrix );
        pElement->GetResidualVector  ( EleResVector   );

        // Trans(T)*A*T
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

		// Inv(T)*b
        Apply_Rotation_To_Vector( HONodesId, EleResVector );

        // Assembling to global matrix 
        int eleSize = EleIdVector.size();

        if ( mColdPlasmaMode == true )
        {
             Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdVector, EleStiffMatrix );    
        }
        else
        {
             std::complex<double> cKij;
             
             for ( int i_Klocal = 0; i_Klocal < eleSize; i_Klocal++ )
             {
                  int i_Kglobal = EleIdVector[i_Klocal];
                  
                  if ( i_Kglobal < mSystemSize )
                  {
                       for ( int j_Klocal = 0; j_Klocal < eleSize; j_Klocal++ )
                       {
                            int j_Kglobal = EleIdVector[j_Klocal];
                            
                            if ( ( i_Kglobal <= j_Kglobal ) && (j_Kglobal < mSystemSize) )
                            {
                                 cKij = EleStiffMatrix[i_Klocal][j_Klocal];

                                 if ( abs(cKij) > 0.0 ) A_matrix( i_Kglobal, j_Kglobal ) += cKij;
                            }
                       }
                  }
             }
        }

        //Assembling to global source vector
		for ( int i_blocal=0; i_blocal < eleSize; i_blocal++ )
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

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - Matrix rotations
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*****************************************************************************************************************
    // - Trans(T) * K * T 
    //*****************************************************************************************************************
    void Modeler::TransT_K_T( Matrix< std::complex<double> >& K, 
                              Matrix< std::complex<double> >& T, 
                              Vector<int>& NodesIds, int NodeIndex )
    {
        int NumDofs = 3 * NodesIds.size();

        if ( mPotentials_On ) 
        {
            for ( int i=0; i<NodesIds.size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( NodesIds[i]-1 );
                
                if ( pNode->pDof( cVs ) != NULL ) 
                {
                    NumDofs++;
                }
            }
        }
       
        int Step = NodesIds.size();

        std::complex<double> Dof_X;
        std::complex<double> Dof_Y;
        std::complex<double> Dof_Z;
        
        // K * T
        for ( int i=0; i<NumDofs; i++ )
        {
            Dof_X = K[ i ][ NodeIndex          ];
            Dof_Y = K[ i ][ NodeIndex +   Step ];
            Dof_Z = K[ i ][ NodeIndex + 2*Step ];
        
            K[ i ][ NodeIndex          ] = T[0][0]*Dof_X + T[0][1]*Dof_Y + T[0][2]*Dof_Z;
            K[ i ][ NodeIndex +   Step ] = T[1][0]*Dof_X + T[1][1]*Dof_Y + T[1][2]*Dof_Z;
            K[ i ][ NodeIndex + 2*Step ] = T[2][0]*Dof_X + T[2][1]*Dof_Y + T[2][2]*Dof_Z;
        }
        
        // Transpose(T) * K
        for ( int j=0; j<NumDofs; j++ )
        {
            Dof_X = K[ NodeIndex          ][ j ];
            Dof_Y = K[ NodeIndex +   Step ][ j ];
            Dof_Z = K[ NodeIndex + 2*Step ][ j ];
        
        	K[ NodeIndex          ][ j ] = T[0][0]*Dof_X + T[0][1]*Dof_Y + T[0][2]*Dof_Z;
            K[ NodeIndex +   Step ][ j ] = T[1][0]*Dof_X + T[1][1]*Dof_Y + T[1][2]*Dof_Z;
            K[ NodeIndex + 2*Step ][ j ] = T[2][0]*Dof_X + T[2][1]*Dof_Y + T[2][2]*Dof_Z;
        }
    }

    //*****************************************************************************************************************
    // - Trans(T) * B 
    //*****************************************************************************************************************
    void Modeler::TransT_B( Vector< std::complex<double> >& B, 
                            Matrix< std::complex<double> >& T, 
                            Vector<int>& NodesIds, int NodeIndex )
    {
        int Step = NodesIds.size();
   
        std::complex<double> Dof_X;
        std::complex<double> Dof_Y;
        std::complex<double> Dof_Z;
        
	    Dof_X = B[ NodeIndex          ];
        Dof_Y = B[ NodeIndex +   Step ];
        Dof_Z = B[ NodeIndex + 2*Step ];
        
        B[ NodeIndex          ] = T[0][0]*Dof_X + T[0][1]*Dof_Y + T[0][2]*Dof_Z;
        B[ NodeIndex +   Step ] = T[1][0]*Dof_X + T[1][1]*Dof_Y + T[1][2]*Dof_Z;
        B[ NodeIndex + 2*Step ] = T[2][0]*Dof_X + T[2][1]*Dof_Y + T[2][2]*Dof_Z;
    }

	//*****************************************************************************************************************
    // - Rotates the stiffness matrix to the n, t, b coordinate system ( for PEC/PMC boundary conditions )
    //*****************************************************************************************************************
    void Modeler::Apply_Rotation_To_Matrix( Vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix )
    {
		if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;

		int local_i = 0;

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit, ++local_i )
        {
            if ( mNormals.find( *hoit ) != mNormals.end() )
            {
                // Rotation matrix to the n, t, b coordinate system
                Matrix< std::complex<double> > RotMatrix;

                Vector<double> n( mNormals[ *hoit ] );
        	    
                Get_Rotation_Matrix( RotMatrix, n );
                
                TransT_K_T( StiffMatrix, RotMatrix, HONodesId, local_i );

                // Applying {En} = [ R * invEp * Trans(R) ] * {Dn} on plasma sheath BC
                if ( ( mColdPlasmaMode == true ) && ( mType_Of_BC_Normal[*hoit] == 'D' ) )
                {
                    Matrix< std::complex<double> > RiEpRt;

                    Get_En_R_invEp_Rt_Dn( RiEpRt, *hoit );

                    TransT_K_T( StiffMatrix, RiEpRt, HONodesId, local_i );
                }
            }
        }
    }

	//*****************************************************************************************************************
    // - Calculates rotation matrix to n, t, b coordinate system
    //*****************************************************************************************************************
    void Modeler::Get_Rotation_Matrix( Matrix< std::complex<double> >& RotMatrix, Vector<double>& n )
    {
        RotMatrix.Resize( 3, 3 );
        
        Vector<double> t( 3 );
        Vector<double> b( 3 );

        TangencialCoordinates( n, t, b );

        RotMatrix[0][0] = n[0]; RotMatrix[0][1] = n[1]; RotMatrix[0][2] = n[2];
        RotMatrix[1][0] = t[0]; RotMatrix[1][1] = t[1]; RotMatrix[1][2] = t[2];
        RotMatrix[2][0] = b[0]; RotMatrix[2][1] = b[1]; RotMatrix[2][2] = b[2];
    }

	//*****************************************************************************************************************
    // - Calculates the matrix that relates En with Dn -> {En} = [ R * invEp * Trans(R) ] * {Dn}
    //*****************************************************************************************************************
    void Modeler::Get_En_R_invEp_Rt_Dn( Matrix< std::complex<double> >& R_invEp_Rt, int NodeId )
    {
        // Permittivity tensor on node ( {Dxyz} = [TEp] * {Exyz} )
        Matrix< std::complex<double> > TEp;

        Node::Pointer pNode = mpModel->GetNode( NodeId - 1 );

        mpColdPlasma->Get_PermittivityTensor_InNode( TEp, pNode, mProblemFrequency );

        // Inverse of TEp
        Matrix< std::complex<double> > invTEp( 3, 3 );

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
        Vector<double> n( mNormals[NodeId] );
        Vector<double> t( 3 );
        Vector<double> b( 3 );

        TangencialCoordinates( n, t, b );

        // Rotation matrix {Entb} = [R] * {Exyz}
        Matrix< std::complex<double> > R( 3, 3 );

        R[0][0] = n[0];  R[0][1] = n[1];  R[0][2] = n[2];
        R[1][0] = t[0];  R[1][1] = t[1];  R[1][2] = t[2];
        R[2][0] = b[0];  R[2][1] = b[1];  R[2][2] = b[2];

        // En = [ R * invEp * Trans(R) ] * Dn
        // Carefull!!, T.Transpose() changes T permanently
        R_invEp_Rt.Resize( 3, 3 );

		R_invEp_Rt = R * invTEp;

		R_invEp_Rt = R_invEp_Rt * R.Transpose();
    }

	//*************************************************************************************************
    //* - Inv(T)*StiffMatrix*T to rotate stiffness matrix to the axisymmetric coordinate system
    //*************************************************************************************************
    void Modeler::Axisym_TKT( std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix )
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
				Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( *hoit, n, t, b );

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

    //******************************************************************************************************
    //* - Inv(T)*StiffMatrix*T to rotate matrix to the diagonal permittivity tensor base (RLP)
    //* - Being:
    //* - K_RPL = [U][T] * K_CC * [Trans(T)][Herm(U)]
    //* - E_SDP = [T] * E_CC  //  E_RLP = [U] * E_SDP
    //******************************************************************************************************
    void Modeler::ColdPlasma_TKT( std::vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix )
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

	//*************************************************************************************************************
    // - Rotates the force vector to the n, t, b coordinate system ( for PEC/PMC boundary conditions )
    //*************************************************************************************************************
    void Modeler::Apply_Rotation_To_Vector( Vector<int>& HONodesId, Vector<std::complex<double> >& ResVector )
    {
        if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) return;

		int local_i = 0;

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit, ++local_i )
        {
            if ( mNormals.find( *hoit ) != mNormals.end() )
            {
                // Rotate vector to n, t, b coordinate system
                Matrix< std::complex<double> > RotMatrix;

                Vector<double> n( mNormals[ *hoit ] );
        	    
                Get_Rotation_Matrix( RotMatrix, n );

                TransT_B( ResVector, RotMatrix, HONodesId, local_i );

                // Applying {En} = [ R * invEp * Trans(R) ] * {Dn} on plasma sheath BC
                if ( ( mColdPlasmaMode == true ) && ( mType_Of_BC_Normal[*hoit] == 'D' ) )
                {
                    Matrix< std::complex<double> > RiEpRt;

                    Get_En_R_invEp_Rt_Dn( RiEpRt, *hoit );

                    TransT_B( ResVector, RiEpRt.Transpose(), HONodesId, local_i );
                }
            }
        }
    }

	//*************************************************************************************************
    //* - Inv(T)*b for axisymmetric problems.
    //*************************************************************************************************
    void Modeler::Axisym_TR( std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector )
    {
        int local_i = 0;

        double distToAxis;

        std::vector<int>::iterator hoit;

        for ( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
			if ( ( mNormals.find(*hoit) == mNormals.end() ) && ( DistanceToAxis(*hoit) > mGiDTolerance ) )
            {
				//////////////  Local coordinate system  //////////////////
				Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( *hoit, n, t, b );

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
    void Modeler::ColdPlasma_TR( std::vector<int>& HONodesId, Vector<std::complex<double> >& ResVector )
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

	//***********************************************************************************************************
    //* - Local system of coordinates in PEC surface for n
    //***********************************************************************************************************
    void Modeler::TangencialCoordinates( Vector<double>& n, Vector<double>& t, Vector<double>& b )
    {
         double tnorm = sqrt( n[0]*n[0] + n[1]*n[1] );

         if ( tnorm > 0.00 )
         {
             t[0] = n[1] / tnorm;
             t[1] =-n[0] / tnorm;
             t[2] = 0.00;
         }
         else
         {
             t[0] = n[2];
             t[1] = 0.00;
             t[2] = 0.00;
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
	double Modeler::DistanceToAxis( unsigned int NodeId )
	{
		Node::Pointer pNode = mpModel->GetNode(NodeId-1);

		double X_coord = pNode->X();
		double Z_coord = pNode->Z();

		double dist = sqrt( X_coord*X_coord + Z_coord*Z_coord );

		return dist;
	}

	//**********************************************************************************************************************************
    //* - Axisymmetric unitary vectors
    //**********************************************************************************************************************************
    void Modeler::AxisymmetricCoordinates( unsigned int NodeId, Vector<double>& n, Vector<double>& t, Vector<double>& b )
    {
		Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

		t[0] = pNode->X();
		t[1] = 0.00;
		t[2] = pNode->Z();

		double tnorm = sqrt( t[0]*t[0] + t[2]*t[2] );

		if ( tnorm > 0.00 )
		{
			t[0] = t[0] / tnorm;
			t[2] = t[2] / tnorm;
		}

		n[0] =-t[2];
		n[1] = 0.00;
		n[2] = t[0];

		b[0] = 0.00;
		b[1] = 1.00;
		b[2] = 0.00;
    }

    //**************************************************************************************************************************************************************
    //* - Trans(T)*StiffMatrix*T, PBC nodes
    //**************************************************************************************************************************************************************
    void Modeler::PBC_FrontBack_TKT( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
    {
	    // Checking if the PBC Front-Back boundary condition is in use
        if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );

		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Front.find( *ohoit ) != mPBC_NodeNodePairs_Front.end() )
            {
			    PBC_FB_NodeToNode_cID( *ohoit, local_i, HONodesId, IdVector );            
            }
			else if ( mPBC_NodeEdgePairs_Front.find( *ohoit ) != mPBC_NodeEdgePairs_Front.end() )
			{
			    PBC_FB_NodeToEdge_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix );      
			}
			else if ( mPBC_NodeElementPairs_Front.find( *ohoit ) != mPBC_NodeElementPairs_Front.end() )
			{
			    PBC_FB_NodeToElement_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix ); 
			}

            local_i++;
        }
    }

	//**************************************************************************************************************************************************************
    //* - Trans(T)*StiffMatrix*T, PBC node to node
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FB_NodeToNode_cID( int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector )
	{
		// Back node assigned to NodeId Front node
        int bNodeId = mPBC_NodeNodePairs_Front[ NodeId ];
		
		Node::Pointer pbNode = mpModel->GetNode( bNodeId-1 );

        // Updating DOFs Id Vector
        if ( mPotentials_On )
        {
            int numVdofs = 0;

            for ( int i=0; i<HONodesId.size(); i++ )
            {
                Node::Pointer pfNode = mpModel->GetNode( HONodesId[ i ]-1 );
                
                if ( pfNode->pDof( cVs ) != NULL ) numVdofs++;
            }
            
            int Step = ( IdVector.size() - numVdofs ) / 3;
            
            IdVector[ local_i          ] = pbNode->pDofcAx()->EquationId();
		    IdVector[ local_i +   Step ] = pbNode->pDofcAy()->EquationId();
		    IdVector[ local_i + 2*Step ] = pbNode->pDofcAz()->EquationId();   

            if ( ( pbNode->pDof( cVs ) != NULL ) && ( local_i < numVdofs ) )
            {
                IdVector[ local_i + 3*Step ] = pbNode->pDofcVs()->EquationId();  
            }
        }
        else
        {
            int Step = IdVector.size() / 3;
            
            IdVector[ local_i          ] = pbNode->pDofcEx()->EquationId();
		    IdVector[ local_i +   Step ] = pbNode->pDofcEy()->EquationId();
		    IdVector[ local_i + 2*Step ] = pbNode->pDofcEz()->EquationId();        
        }

        // Updating HO nodes Id vector
        HONodesId[ local_i ] = pbNode->Id();  
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*StiffMatrix*T, PBC node to edge
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FB_NodeToEdge_TKT( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
	{
		std::vector<int> bNodeIds( mPBC_NodeEdgePairs_Front[ NodeId ] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodeIds.size(); in++ ) 
        {
            pbNodes.push_back( mpModel->GetNode( bNodeIds[ in ]-1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, pbNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element edge
		double natCoord = mPBC_NodeEdgeCoord_Front[ NodeId ];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_StiffMatrixReDef( local_i, oldSize, newSize, N, StiffMatrix );

		local_i = local_i + ( bNodeIds.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Redefinition of IdVector for PBC boundary conditions
    //**************************************************************************************************************************************************************
	void Modeler::PBC_IdVectorReDef( int local_i, std::vector<Node::Pointer>& pNodes, std::vector<int>& HONodesId, Vector<int>& IdVector )
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

	//**************************************************************************************************************************************************************
    //* - Inv(T)*StiffMatrix*T, PBC node to element
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TKT( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
	{
		std::vector<int> bNodes( mPBC_NodeElementPairs_Front[ NodeId ] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) 
        {
            pbNodes.push_back( mpModel->GetNode( bNodes[ in ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, pbNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Front[ NodeId ] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_StiffMatrixReDef( local_i, oldSize, newSize, N, StiffMatrix );

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Values of the Back element base N for the proyected location of the Front node
    //**************************************************************************************************************************************************************
    void Modeler::PBC_NaturalBase_N_Surface( std::vector<double>& N, std::vector<double>& natCoord )
    {
        double L2 = natCoord[0];
        double L3 = natCoord[1];
        double L1 = 1.00 - L2 - L3;

        if ( mElementOrder == 0 )
        {
		    N.resize(6);

            N[0] = L1;
            N[1] = L2;
            N[2] = L3;
            N[3] = 32.00 * L1 * L2 * L3 * ( 4.00*L1 - 1.00 );
            N[4] = 32.00 * L1 * L2 * L3 * ( 4.00*L2 - 1.00 );
            N[5] = 32.00 * L1 * L2 * L3 * ( 4.00*L3 - 1.00 ); 
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
    }

	//**************************************************************************************************************************************************************
    //* - Values of the Back element edge base N for the proyected location of the Front node
    //**************************************************************************************************************************************************************
    void Modeler::PBC_NaturalBase_N_Line( std::vector<double>& N, double natCoord )
    {
        double L2 = natCoord;
        double L1 = 1.00 - L2;

        if ( mElementOrder == 0 )
        {
		    N.resize(2);

            N[0] = L1;
            N[1] = L2;
        }
        else if ( mElementOrder == 1 )
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
    }

	//**************************************************************************************************************************************************************
    //* - Redefinition of the stiffness matrix for PBC boundary conditions
    //**************************************************************************************************************************************************************
	void Modeler::PBC_StiffMatrixReDef( int local_i, int oldSize, int newSize, std::vector<double>& N, Matrix<std::complex<double> >& StiffMatrix )
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

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC nodes
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FrontBack_TR( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
        // Checking if the PBC Front-Back boundary condition is in use
        if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );

		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Front.find( *ohoit ) != mPBC_NodeNodePairs_Front.end() )
            {
			    PBC_FB_NodeToNode_cID( *ohoit, local_i, HONodesId, IdVector );                
            }
			else if ( mPBC_NodeEdgePairs_Front.find( *ohoit ) != mPBC_NodeEdgePairs_Front.end() )
			{
			    PBC_FB_NodeToEdge_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );       
			}
			else if ( mPBC_NodeElementPairs_Front.find( *ohoit ) != mPBC_NodeElementPairs_Front.end() )
			{
			    PBC_FB_NodeToElement_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );  
			}
            local_i++;
        }
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to edge
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FB_NodeToEdge_TR( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
		std::vector<int> bNodes( mPBC_NodeEdgePairs_Front[ NodeId ] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) 
        {
            pbNodes.push_back( mpModel->GetNode( bNodes[ in ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, pbNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element edge
		double natCoord = mPBC_NodeEdgeCoord_Front[ NodeId ];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_ResVectorReDef( local_i, oldSize, newSize, N, ResVector );

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to element
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FB_NodeToElement_TR( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
		std::vector<int> bNodes( mPBC_NodeElementPairs_Front[ NodeId ] );

		std::vector<Node::Pointer> pbNodes;

		for ( int in = 0; in < bNodes.size(); in++ ) 
        {
            pbNodes.push_back( mpModel->GetNode( bNodes[ in ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, pbNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Back natural base values of the proyected Front node on the Back element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Front[ NodeId ] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_ResVectorReDef(local_i, oldSize, newSize, N, ResVector);

		local_i = local_i + ( bNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Redefinition of the residual vector for PBC boundary conditions
    //**************************************************************************************************************************************************************
	void Modeler::PBC_ResVectorReDef( int local_i, int oldSize, int newSize, std::vector<double>& N, Vector<std::complex<double> >& ResVector )
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

	//**************************************************************************************************************************************************************
    //* - Trans(T)*StiffnesMatrix*T, PBC nodes
    //**************************************************************************************************************************************************************
    void Modeler::PBC_RightLeft_TKT( std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
   {
	    // Checking if the PBC Right-Left boundary condition is in use
        if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 ) 
        {
            return;
        }

		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );

		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Right.find( *ohoit ) != mPBC_NodeNodePairs_Right.end() )
            {
				PBC_RL_NodeToNode_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix );
            }
			else if ( mPBC_NodeEdgePairs_Right.find( *ohoit ) != mPBC_NodeEdgePairs_Right.end() )
			{
				PBC_RL_NodeToEdge_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix );
			}
			else if ( mPBC_NodeElementPairs_Right.find( *ohoit ) != mPBC_NodeElementPairs_Right.end() )
			{
				PBC_RL_NodeToElement_TKT( *ohoit, local_i, HONodesId, IdVector, StiffMatrix );
			}
            local_i++;
        }
    }

	//**************************************************************************************************************************************************************
    //* - Trans(T)*StiffnesMatrix*T, PBC node to node
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToNode_TKT( int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
	{
		int lNode = mPBC_NodeNodePairs_Right[ NodeId ];
		int size  = IdVector.size();
		int step  = size / 3;

		Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

		HONodesId[local_i         ] = plNode->Id();

		IdVector [local_i         ] = plNode->pDofcEx()->EquationId();
		IdVector [local_i +   step] = plNode->pDofcEy()->EquationId();
		IdVector [local_i + 2*step] = plNode->pDofcEz()->EquationId();

		std::vector<double> N(1);

		N[0] = 1.00;

		PBC_StiffMatrixReDef_Rot( local_i, size, size, N, NodeId, StiffMatrix );
	}

    //**************************************************************************************************************************************************************
    //* - Trans(T)*StiffnesMatrix*T, PBC node to edge
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToEdge_TKT( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
	{
		std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[NodeId] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) 
        {
            plNodes.push_back( mpModel->GetNode( lNodes[ i ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, plNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Right natural base values of the proyected Left node on the Right element edge
		double natCoord = mPBC_NodeEdgeCoord_Right[NodeId];

		std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_StiffMatrixReDef_Rot( local_i, oldSize, newSize, N, NodeId, StiffMatrix );

		local_i = local_i + ( lNodes.size() - 1 );
	}

    //**************************************************************************************************************************************************************
    //* - Trans(T)*StiffnesMatrix*T, PBC node to element
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TKT( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Matrix<std::complex<double> >& StiffMatrix )
	{
		std::vector<int> lNodes( mPBC_NodeElementPairs_Right[ NodeId ] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) 
        {
            plNodes.push_back( mpModel->GetNode( lNodes[ i ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, plNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Right natural base values of the proyected Left node on the Right element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Right[ NodeId ] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_StiffMatrixReDef_Rot( local_i, oldSize, newSize, N, NodeId, StiffMatrix );

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Redefinition of the stiffness matrix for PBC boundary conditions
    //**************************************************************************************************************************************************************
	void Modeler::PBC_StiffMatrixReDef_Rot( int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Matrix<std::complex<double> >& StiffMatrix )
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

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC nodes
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RightLeft_TR( std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
	    if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 ) 
        {
            return;
        }

		int local_i = 0;

		std::vector<int> OldHONodesId( HONodesId );

		std::vector<int>::iterator ohoit;

        for( ohoit = OldHONodesId.begin(); ohoit != OldHONodesId.end(); ++ohoit )
        {
            if ( mPBC_NodeNodePairs_Right.find( *ohoit ) != mPBC_NodeNodePairs_Right.end() )
            {
				PBC_RL_NodeToNode_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );
            }
			else if ( mPBC_NodeEdgePairs_Right.find( *ohoit ) != mPBC_NodeEdgePairs_Right.end() )
			{
				PBC_RL_NodeToEdge_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );
			}
			else if ( mPBC_NodeElementPairs_Right.find( *ohoit ) != mPBC_NodeElementPairs_Right.end() )
			{
				PBC_RL_NodeToElement_TR( *ohoit, local_i, HONodesId, IdVector, ResVector );
			}
            local_i++;
        }
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to node
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToNode_TR( int NodeId, int local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
		int lNode = mPBC_NodeNodePairs_Right[NodeId];
		int size  = IdVector.size();
		int step  = size / 3;

		Node::Pointer plNode = mpModel->GetNode( lNode - 1 );

		HONodesId[local_i         ] = plNode->Id();

		IdVector [local_i         ] = plNode->pDofcEx()->EquationId();
		IdVector [local_i +   step] = plNode->pDofcEy()->EquationId();
		IdVector [local_i + 2*step] = plNode->pDofcEz()->EquationId();

		std::vector<double> N(1);

		N[0] = 1.00;

		PBC_ResVectorReDef_Rot( local_i, size, size, N, NodeId, ResVector );
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to edge
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToEdge_TR( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
		std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[ NodeId ] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) 
        {
            plNodes.push_back( mpModel->GetNode( lNodes[ i ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, plNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Left natural base values of the proyected Right node on the Left element edge
		double natCoord = mPBC_NodeEdgeCoord_Right[ NodeId ];

        std::vector<double> N;

        PBC_NaturalBase_N_Line( N, natCoord );

		PBC_ResVectorReDef_Rot( local_i, oldSize, newSize, N, NodeId, ResVector );

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to element
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RL_NodeToElement_TR( int NodeId, int& local_i, std::vector<int>& HONodesId, Vector<int>& IdVector, Vector<std::complex<double> >& ResVector )
	{
		std::vector<int> lNodes( mPBC_NodeElementPairs_Right[ NodeId ] );

		std::vector<Node::Pointer> plNodes;

		for ( int i = 0; i < lNodes.size(); i++ ) 
        {
            plNodes.push_back( mpModel->GetNode( lNodes[ i ] - 1 ) );
        }

		int oldSize = IdVector.size();

		PBC_IdVectorReDef( local_i, plNodes, HONodesId, IdVector );

		int newSize = IdVector.size();

        // Left natural base values of the proyected Right node on the Left element
		std::vector<double> natCoord( mPBC_NodeElementCoord_Right[ NodeId ] );

        std::vector<double> N;

        PBC_NaturalBase_N_Surface( N, natCoord );

		PBC_ResVectorReDef_Rot( local_i, oldSize, newSize, N, NodeId, ResVector );

		local_i = local_i + ( lNodes.size() - 1 );
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T)*b, PBC node to element
    //**************************************************************************************************************************************************************
	void Modeler::PBC_ResVectorReDef_Rot( int local_i, int oldSize, int newSize, std::vector<double>& N, int NodeId, Vector<std::complex<double> >& ResVector )
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

    //****************************************************************************************************************************************************
    // - Apply Trans(T) * StiffMatrix * T to eliminate one of the nodes in a contact pair
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Matrix_TKT( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix )
    {
		if ( mContactPairs.size() == 0 ) 
        {
            return;
        }

        for ( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if ( mContactPairs.find( HONodesId[ local_i ] ) != mContactPairs.end() )
            {
                if ( mPotentials_On )
                {
                    Apply_Contact_TKT_AV( local_i, HONodesId, DIdVector, StiffMatrix );
                }
                else
                {
                    Apply_Contact_TKT_Ef( local_i, HONodesId, DIdVector, StiffMatrix );
                }
            }
        }
    }

    //****************************************************************************************************************************************************
    // - Apply Trans(T) * StiffMatrix * T to eliminate one of the nodes in a contact pair ( E field fromulation )
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_TKT_Ef( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix )
    {   
        // Get contact matrix
        Matrix< std::complex<double> > CtCMatrix( 3, 3 );

        Get_Contact_Matrix_FullWave( CtCMatrix, HONodesId[ local_i ] );

        // Step size and DOF size
        int Step    =   HONodesId.size();
        int DofSize = 3*HONodesId.size();

        //  Inv(T) * StiffnesMatrix * T     
        std::complex<double> dof_x;
        std::complex<double> dof_y;
        std::complex<double> dof_z;
        
        // StiffnesMatrix * T
        for ( int i=0; i<DofSize; i++ )
        {
            dof_x = StiffMatrix[ i ][ local_i          ];
            dof_y = StiffMatrix[ i ][ local_i +   Step ];
            dof_z = StiffMatrix[ i ][ local_i + 2*Step ];
        
            StiffMatrix[ i ][ local_i          ] = CtCMatrix[0][0]*dof_x + CtCMatrix[1][0]*dof_y + CtCMatrix[2][0]*dof_z;
            StiffMatrix[ i ][ local_i +   Step ] = CtCMatrix[0][1]*dof_x + CtCMatrix[1][1]*dof_y + CtCMatrix[2][1]*dof_z;
            StiffMatrix[ i ][ local_i + 2*Step ] = CtCMatrix[0][2]*dof_x + CtCMatrix[1][2]*dof_y + CtCMatrix[2][2]*dof_z;
        }
        
        // Inv(T) * StiffnesMatrix
        for ( int i=0; i<DofSize; i++ )
        {
            dof_x = StiffMatrix[ local_i          ][ i ];
            dof_y = StiffMatrix[ local_i +   Step ][ i ];
            dof_z = StiffMatrix[ local_i + 2*Step ][ i ];
        
        	StiffMatrix[ local_i          ][ i ] = CtCMatrix[0][0]*dof_x + CtCMatrix[1][0]*dof_y + CtCMatrix[2][0]*dof_z;
            StiffMatrix[ local_i +   Step ][ i ] = CtCMatrix[0][1]*dof_x + CtCMatrix[1][1]*dof_y + CtCMatrix[2][1]*dof_z;
            StiffMatrix[ local_i + 2*Step ][ i ] = CtCMatrix[0][2]*dof_x + CtCMatrix[1][2]*dof_y + CtCMatrix[2][2]*dof_z;
        }   

        // Get Left node
        Node::Pointer pNode = mpModel->GetNode( mContactPairs[ HONodesId[ local_i ] ].ctcNode - 1 );

        // Change node DOFs Ids
        DIdVector[ local_i          ] = pNode->pDofcEx()->EquationId();
        DIdVector[ local_i +   Step ] = pNode->pDofcEy()->EquationId();
        DIdVector[ local_i + 2*Step ] = pNode->pDofcEz()->EquationId();

        // Change node Id
        HONodesId[ local_i ] = pNode->Id();
    }

    //****************************************************************************************************************************************************
    // - Calculates contact matrix for dielectric interfaces
    //****************************************************************************************************************************************************
    void Modeler::Get_Contact_Matrix_FullWave( Matrix< std::complex<double> >& CtCMatrix, int RNodesId )
    {
        ContactPairData& pContactPair = mContactPairs[ RNodesId ];
 
        Properties::Pointer Properties;
        
        double eo   = 8.8541878176e-12 ;
        double freq = mProblemFrequency;
        
        double sigma   ;
        double eps_real;
        double eps_imag;
        
        // Material R
        Properties = mpModel->GetProperties( pContactPair.matR );
        
        sigma    = (*Properties)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
        eps_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
        eps_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
        
        std::complex<double> cEpsR( eps_real, eps_imag + ( sigma / freq ) );
        
        // Material L
        Properties = mpModel->GetProperties( pContactPair.matL );
        
        sigma    = (*Properties)(  IHL_ELECTRIC_CONDUCTIVITY )     ;
        eps_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
        eps_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
        
        std::complex<double> cEpsL( eps_real, eps_imag + ( sigma / freq ) );
                
        std::complex<double> rmat = ( cEpsL / cEpsR ) - 1.00;
        
        Vector<double> n = pContactPair.ctcNormal;
                                                                                                                                    
        CtCMatrix[0][0] = ( n[0]*n[0] ) * rmat + 1.00; CtCMatrix[0][1] = ( n[0]*n[1] ) * rmat       ; CtCMatrix[0][2] = ( n[0]*n[2] ) * rmat       ;
        CtCMatrix[1][0] = ( n[1]*n[0] ) * rmat       ; CtCMatrix[1][1] = ( n[1]*n[1] ) * rmat + 1.00; CtCMatrix[1][2] = ( n[1]*n[2] ) * rmat       ;
        CtCMatrix[2][0] = ( n[2]*n[0] ) * rmat       ; CtCMatrix[2][1] = ( n[2]*n[1] ) * rmat       ; CtCMatrix[2][2] = ( n[2]*n[2] ) * rmat + 1.00;
    }

    //****************************************************************************************************************************************************
    // - Eliminate one of the nodes in a contact pair ( AV potentials fromulation )
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_TKT_AV( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix )
    {
        // Get Left node
        Node::Pointer pNode = mpModel->GetNode( mContactPairs[ HONodesId[ local_i ] ].ctcNode - 1 );

        // Step size and DOF size
        int Step    = HONodesId.size();
        int DofSize = DIdVector.size();

        // Change node DOFs Ids
        DIdVector[ local_i          ] = pNode->pDofcAx()->EquationId();
        DIdVector[ local_i +   Step ] = pNode->pDofcAy()->EquationId();
        DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();

        if ( ( local_i + 3*Step ) < DofSize ) 
        {
            DIdVector[ local_i + 3*Step ] = pNode->pDofcVs()->EquationId();
        }   

        // Change node Id
        HONodesId[ local_i ] = pNode->Id();
    }

	//****************************************************************************************************************************************************
    // - Inv(T) * b in contact nodes
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Vector_TR( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector<std::complex<double> >& ResVector )
    {
        if ( mContactPairs.size() == 0 ) 
        {
            return;
        }

		for ( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if ( mContactPairs.find( HONodesId[ local_i ] ) != mContactPairs.end() )
            {
                // Get Left node
				Node::Pointer pNode = mpModel->GetNode( mContactPairs[ HONodesId[ local_i ] ].ctcNode - 1 );
                    
                // Step size 
                int Step = HONodesId.size();
 
                if ( mPotentials_On )
                {
                    // Change node DOFs Ids
				    DIdVector[ local_i          ] = pNode->pDofcAx()->EquationId();
				    DIdVector[ local_i +   Step ] = pNode->pDofcAy()->EquationId();
				    DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();                
                }
                else
                {
                    // Get contact matrix
                    Matrix< std::complex<double> > CtCMatrix( 3, 3 );
                    
                    Get_Contact_Matrix_FullWave( CtCMatrix, HONodesId[ local_i ] );     
                    
                    // Inv(T) * b  
                    std::complex<double> dof_x;
                    std::complex<double> dof_y;
                    std::complex<double> dof_z;
                    
				    dof_x = ResVector[ local_i          ];
                    dof_y = ResVector[ local_i +   Step ];
                    dof_z = ResVector[ local_i + 2*Step ];
                    
                    ResVector[ local_i          ] = CtCMatrix[0][0]*dof_x + CtCMatrix[1][0]*dof_y + CtCMatrix[2][0]*dof_z;
                    ResVector[ local_i +   Step ] = CtCMatrix[0][1]*dof_x + CtCMatrix[1][1]*dof_y + CtCMatrix[2][1]*dof_z;
                    ResVector[ local_i + 2*Step ] = CtCMatrix[0][2]*dof_x + CtCMatrix[1][2]*dof_y + CtCMatrix[2][2]*dof_z;
                    
                    // Change node DOFs Ids
				    DIdVector[ local_i          ] = pNode->pDofcEx()->EquationId();
				    DIdVector[ local_i +   Step ] = pNode->pDofcEy()->EquationId();
				    DIdVector[ local_i + 2*Step ] = pNode->pDofcEz()->EquationId();
                }

                // Change node Id
				HONodesId[ local_i ] = pNode->Id();
            }
        }
    }

    //****************************************************************************************************************************************************
    // - Trans(T)*StiffMatrix*T to eliminate one of the nodes in a contact pair for Dielectric-Plasma interfaces
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Matrix_TKT_ColdPlasma( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix<std::complex<double> >& StiffMatrix )
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

    			DIdVector[local_i        ] = pNode->pDofcEx()->EquationId();
    			DIdVector[local_i + step ] = pNode->pDofcEy()->EquationId();
    			DIdVector[local_i + step2] = pNode->pDofcEz()->EquationId();

                //////////////////////////////////////////////    Inv(T)*StiffnesMatrix*T    /////////////////////////////////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                Matrix< std::complex<double> > ContactMatrix(3,3);

                Get_Contact_Matrix_ColdPlasma( ContactMatrix, RNodeId );

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
    void Modeler::Get_Contact_Matrix_ColdPlasma( Matrix< std::complex<double> >& CtCMatrix, int RNodesId )
    {
        ContactPairData& pContactPair = mContactPairs[RNodesId];

        Node::Pointer prNode = mpModel->GetNode( RNodesId             - 1 );
        Node::Pointer plNode = mpModel->GetNode( pContactPair.ctcNode - 1 );

        Vector<double> n = pContactPair.ctcNormal;
        Vector<double> t( 3 );
        Vector<double> b( 3 );

        TangencialCoordinates( n, t, b );

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
        CtCMatrix.Resize(3,3);

        CtCMatrix = invTR * TL;
    }

    //****************************************************************************************************************************************************
    // - Inv(T)*b in contact nodes for Dielectric-Plasma interfaces
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Vector_TR_ColdPlasma( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector<std::complex<double> >& ResVector )
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

    			DIdVector[local_i        ] = pNode->pDofcEx()->EquationId();
    			DIdVector[local_i + step ] = pNode->pDofcEy()->EquationId();
    			DIdVector[local_i + step2] = pNode->pDofcEz()->EquationId();

                ///////////////////////////////////////////////////  Inv(T)*b  //////////////////////////////////////////////////
                std::complex<double> dof_x, dof_y, dof_z;

                Matrix< std::complex<double> > ContactMatrix;

                Get_Contact_Matrix_ColdPlasma( ContactMatrix, RNodeId );

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

    //*************************************************************************************************************************************
    // - Returns the determinant of a square matrix M of size n (for complex matrices)
    //*************************************************************************************************************************************
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

    //**********************************************************************************************************************
    // - Fix all DOFs in a node to zero  
    //**********************************************************************************************************************
    void Modeler::Fix_All_DOFs_In_Node( Node::Pointer pNode )
    {
        std::complex<double> cZero( 0.00, 0.00 );
        
        if ( mPotentials_On )
        {
            (*mpModel)( cAx, *pNode ) = cZero; 
            (*mpModel)( cAy, *pNode ) = cZero; 
            (*mpModel)( cAz, *pNode ) = cZero;  

            pNode->pDofcAx()->FixDof();
            pNode->pDofcAy()->FixDof();
            pNode->pDofcAz()->FixDof();
        
            if ( pNode->pDof( cVs ) != NULL ) 
            {
                (*mpModel)( cVs, *pNode ) = cZero; 
                
                pNode->pDofcVs()->FixDof(); 
            }
        }                            
        else                         
        {                            
            (*mpModel)( cEx, *pNode ) = cZero; 
            (*mpModel)( cEy, *pNode ) = cZero; 
            (*mpModel)( cEz, *pNode ) = cZero; 

            pNode->pDofcEx()->FixDof();
            pNode->pDofcEy()->FixDof();
            pNode->pDofcEz()->FixDof();
        }                                
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - PEC, PMC, PBC and projection elements
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//**********************************************************************************************************************
    // - Assigns each PBC node to its corresponding std::map
    //**********************************************************************************************************************
    void Modeler::GeneratePBCElement( int* NodesId, unsigned int PropertiesId )
    {
        // Each material is assigned to a different group to speed the PBC set-up
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        int MaterialId = NodesId[3];

		int FaceTypeId = (*properties)( COMPLEX_IBC )[0];

		mPBCTolerance  = (*properties)( COMPLEX_IBC )[1];

		std::vector<int> HONodesId;

        PushHONodesOnSurface( NodesId, HONodesId );

		if      ( FaceTypeId == 11 ) mPBCFrontElements[ MaterialId ].push_back( HONodesId );
		else if ( FaceTypeId == 12 ) mPBCBackElements [ MaterialId ].push_back( HONodesId );
		else if ( FaceTypeId == 21 ) mPBCRightElements[ MaterialId ].push_back( HONodesId );
		else if ( FaceTypeId == 22 ) mPBCLeftElements [ MaterialId ].push_back( HONodesId );
	}

	//**********************************************************************************************************************
    // - Assigns the nodes in a periodic surface to the elements in the other
	//   corresponding periodic surface.
	// - Geometry must be placed carefully:
	//   1) Front-Back surfaces in the XY-plane.
	//   2) Left surface in the 00-YZ-Plane.
	//	 3) For cylindrical symmetry the central axis must be placed along the Z axis.
    //**********************************************************************************************************************
	void Modeler::SetPBC()
	{
        // Setting PBC for Front-Back surfaces
		if ( mPBCFrontElements.size() > 0 ) Set_FrontBack_PBC();

        //// Check if Front-Back PBC is compatible with LL2P 3sb elements
        //if ( ( mElementOrder == 0 ) && ( ( mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) > 0 ) )
        //{
        //    Send_Error_Msg( "PBC for LL2P 3sb elements is only valid when Front-Back surface meshes are equal", 
        //                    "Please, re-mesh or change element type.", 1 );            
        //}
	
        // Setting PBC for Right-Left surfaces
        if ( mPBCRightElements.size() > 0 ) Set_RightLeft_PBC();

        //// Check if Right-Left PBC is compatible with LL2P 3sb elements
        //if ( ( mElementOrder == 0 ) && ( ( mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) > 0 ) )
        //{
        //    Send_Error_Msg( "PBC for LL2P 3sb elements is only valid when Right-Left surface meshes are equal", 
        //                    "Please, re-mesh or change element type.", 1 );                         
        //}
	}

	//**********************************************************************************************************************
    // - Assigns the nodes in the front periodic surface to the elements in the back periodic surface.
	// - Front-back surfaces must be placed in the XY-plane.
    //**********************************************************************************************************************
    void Modeler::Set_FrontBack_PBC()
    {
	    // Loop through materials in Back PBC
        std::map<unsigned int, std::vector<std::vector<int> > >::iterator it_Back;

        for( it_Back = mPBCBackElements.begin(); it_Back != mPBCBackElements.end(); ++it_Back )
        {
            // Material ID
			int pbcMaterial = it_Back->first;

            // Group of Back elements made of the same material (pbcMaterial)
			std::vector< std::vector<int> > elemSetBack( it_Back->second );

            // Group of Frontal elements made of the same material (pbcMaterial)
			std::vector< std::vector<int> > elemSetFront( mPBCFrontElements[ pbcMaterial ] );

            // List of nodes (sorted and unique) in the Front PBC
            std::list<int> FrontNodesList;

            // Boolean map for the nodes in element vertices 
            std::map< unsigned int, bool > fNodeInVertice;

             // Filling the list with the nodes on each Front element
            std::vector< std::vector<int> >::iterator it_eF;

            for( it_eF = elemSetFront.begin(); it_eF != elemSetFront.end(); ++it_eF )
            {
                std::vector<int>::iterator it_nF;

                int iLocal = 0;

                for( it_nF = (*it_eF).begin(); it_nF != (*it_eF).end(); ++it_nF )
				{
                    FrontNodesList.push_back( *it_nF );

                    if ( iLocal < 3 ) fNodeInVertice[ *it_nF ] = true ;
                    else              fNodeInVertice[ *it_nF ] = false;

                    iLocal++;
                }
            }

            // Sorting list
            FrontNodesList.sort();

            // Collapsing equal nodes Id
            FrontNodesList.unique();

            // Loop through the Back elements to calculate the Frontal nodes projections
			std::vector< std::vector<int> >::iterator it_eB;

            // Loop through Back elements for one material
			for( it_eB = elemSetBack.begin(); it_eB != elemSetBack.end() && FrontNodesList.size() > 0 ; ++it_eB )
            {
                // Recovering Back element nodes coordinates on the XY plane
                std::vector<double> eX, eY;

                for( int i=0; i<(*it_eB).size(); i++ )
                {
                    eX.push_back( mpModel->GetNode( (*it_eB)[i]-1 )->X() );
                    eY.push_back( mpModel->GetNode( (*it_eB)[i]-1 )->Y() );
                }

                // Calculating cross section (only nodes inside this area will be consider)
                double eCross_Section;

                double edge01_Lenght = sqrt( pow( eX[0]-eX[1], 2 ) + pow( eY[0]-eY[1], 2 ) ); 
                double edge02_Lenght = sqrt( pow( eX[0]-eX[2], 2 ) + pow( eY[0]-eY[2], 2 ) );

                if ( edge02_Lenght > edge01_Lenght ) 
                {
                    eCross_Section = ( 1.0 + mPBCTolerance ) * edge02_Lenght;
                }
                else   
                {
                    eCross_Section = ( 1.0 + mPBCTolerance ) * edge01_Lenght;
                }

                // Parameters of the geometrical mapping interpolation
                std::vector<double> a, b;

                // Computing the parameters a, b of the geometrical mapping interpolation
                PBC_Mapping_Parameters( a, b, eX, eY );

                // Front nodes projected on this Back element
                std::vector<int> AssignedFrontNodes;

                // Loop over the list of nodes in the Front PBC
                std::list<int>::iterator it_nFl;

                // Loop over the list of nodes in the Front PBC
                for( it_nFl = FrontNodesList.begin(); it_nFl != FrontNodesList.end(); ++it_nFl )
                {
                    // Recovering front node pointer 
                    Node::Pointer pfNode = mpModel->GetNode( (*it_nFl)-1 );

                    // Cartesina coordinates of the node
                    double carX = pfNode->X();
                    double carY = pfNode->Y();

                    // Checking if the Front node projection is close enough to the Back element
                    if ( sqrt( pow( eX[0]-carX, 2 ) + pow( eY[0]-carY, 2 ) ) > eCross_Section ) continue;

                    // Natural coordinates of the Frontal node projection on the Back element
                    double natX = a[0] + ( a[1] * carX ) + ( a[2] * carY );
                    double natY = b[0] + ( b[1] * carX ) + ( b[2] * carY );

                    // If the Front node is outside the Back element -> continue
                    if ( ( natX < -mPBCTolerance ) || ( natX > ( 1.0 + mPBCTolerance ) ) ) continue; 
                    if ( ( natY < -mPBCTolerance ) || ( natY > ( 1.0 + mPBCTolerance ) ) ) continue; 
                    if ( ( natX + natY )           >           ( 1.0 + mPBCTolerance )   ) continue;

                    // Push assigned node in AssignedFrontNodes list 
                    AssignedFrontNodes.push_back( *it_nFl );

                    // Fix all DOFs of the assigned Front node
                    Fix_All_DOFs_In_Node( pfNode );

                    // Cheking if Frontal node is projected onto a Back node 
                    int  bLocalId = PBC_Node_To_Node( natX, natY );

                    // Only nodes with the same DOFs are assigned
                    if ( bLocalId > -1 )
                    {
                        if ( ( ( fNodeInVertice[ pfNode->Id() ] == true  ) && ( bLocalId < 3 ) ) ||
                             ( ( fNodeInVertice[ pfNode->Id() ] == false ) && ( bLocalId > 2 ) )  )
                        {
                            mPBC_NodeNodePairs_Front[ pfNode->Id() ] = (*it_eB)[ bLocalId ]; 

                            continue;                         
                        }
                    }

                    // Cheking if Frontal node is projected onto a Back element edge
                    std::vector<int> BackEdgeNodeIds;

                    double NatCoord_InEdge = PBC_Node_To_Edge( natX, natY, *it_eB, BackEdgeNodeIds );

                    if ( BackEdgeNodeIds.size() > 0 )
                    {
                        mPBC_NodeEdgePairs_Front[ pfNode->Id() ] = BackEdgeNodeIds;
                        mPBC_NodeEdgeCoord_Front[ pfNode->Id() ] = NatCoord_InEdge;

                        continue;
                    }

                    // If Frontal node is not projected onto a node/edge then it is projected onto a Back element
                    std::vector<double> NatCoord_InElement( 2 );

                    NatCoord_InElement[0] = natX;
                    NatCoord_InElement[1] = natY;

                    mPBC_NodeElementPairs_Front[ pfNode->Id() ] = (*it_eB);
                    mPBC_NodeElementCoord_Front[ pfNode->Id() ] = NatCoord_InElement;

                }// End loop over the list of nodes in the Front PBC

                // Removing assigned Front nodes from the FrontNodesList 
                std::vector<int>::iterator ita;

                for( ita = AssignedFrontNodes.begin(); ita != AssignedFrontNodes.end(); ++ita ) 
                {
                    FrontNodesList.remove( *ita );
                }

            }// End loop through Back elements

            // Assign the remaining nodes to the closest Back element edge 
            std::list<int>::iterator it_nFl;

            for( it_nFl = FrontNodesList.begin(); it_nFl != FrontNodesList.end(); ++it_nFl )
            {
                Assign_Front_Node_To_Closest_Edge( *it_nFl, elemSetBack );
            }

        }// End loop through materials in Back PBC

        // Clearing PBC element maps
	    mPBCFrontElements.clear(); 
        std::map< unsigned int, std::vector< std::vector<int> > >().swap( mPBCFrontElements );

		mPBCBackElements.clear(); 
        std::map< unsigned int, std::vector< std::vector<int> > >().swap( mPBCBackElements );
    }

    //***********************************************************************************************************************************************
	// - Assign a Front node to the closest element edge
	//***********************************************************************************************************************************************
    void Modeler::Assign_Front_Node_To_Closest_Edge( int NodeId, std::vector< std::vector<int> >& ElementSet )
    {
        Node::Pointer pNode = mpModel->GetNode( NodeId-1 );
       
        double MinDistance = 1e12;

        std::vector< std::vector<int> >::iterator it_ES;
       
        for( it_ES = ElementSet.begin(); it_ES != ElementSet.end(); ++it_ES )
        {
            // Loop over element edges 
            for ( int edg=0; edg<3; edg++ )
            {
                int iniEdg = edg    ;
                int endEdg = edg + 1;

                if ( endEdg > 2 ) endEdg = 0;

                double iniX = mpModel->GetNode( (*it_ES)[ iniEdg ]-1 )->X();
                double iniY = mpModel->GetNode( (*it_ES)[ iniEdg ]-1 )->Y();
                                                          
                double endX = mpModel->GetNode( (*it_ES)[ endEdg ]-1 )->X();
                double endY = mpModel->GetNode( (*it_ES)[ endEdg ]-1 )->Y();
                
                double edgX = endX - iniX;
                double edgY = endY - iniY;

                double nodX = pNode->X() - iniX;
                double nodY = pNode->Y() - iniY;
                
                double ModEdge = sqrt( edgX*edgX + edgY*edgY );
       
                double Project = ( nodX*edgX + nodY*edgY ) / ModEdge;
       
                if ( ( Project < 0.0 ) || ( Project > ModEdge ) ) 
                {
                    continue;
                }
       
                double Distance = abs( ( nodY*edgX - nodX*edgY  ) / ModEdge ); 
       
                if ( Distance > MinDistance ) 
                {
                    continue;
                }
       
                MinDistance = Distance;
       
                std::vector<int> EdgeNodeIds;
       
                EdgeNodeIds.push_back( (*it_ES)[ iniEdg ] );
                EdgeNodeIds.push_back( (*it_ES)[ endEdg ] );

                if ( mElementOrder == 2 ) 
                {
                    EdgeNodeIds.push_back( (*it_ES)[ iniEdg + 3 ] );
                }
              
                mPBC_NodeEdgePairs_Front[ NodeId ] = EdgeNodeIds;
                mPBC_NodeEdgeCoord_Front[ NodeId ] = Project / ModEdge;
            }
        }
    }

    //***********************************************************************************************************************************************
	// - Assign a Right node to the closest element edge
	//***********************************************************************************************************************************************
    void Modeler::Assign_Right_Node_To_Closest_Edge( int NodeId, std::vector< std::vector<int> >& ElementSet )
    {
        Node::Pointer pNode = mpModel->GetNode( NodeId-1 );
       
        double MinDistance = 1e12;

        std::vector< std::vector<int> >::iterator it_ES;
       
        for( it_ES = ElementSet.begin(); it_ES != ElementSet.end(); ++it_ES )
        {
            // Loop over element edges 
            for ( int edg=0; edg<3; edg++ )
            {
                int iniEdg = edg    ;
                int endEdg = edg + 1;

                if ( endEdg > 2 ) endEdg = 0;

                double iniX = mpModel->GetNode( (*it_ES)[ iniEdg ]-1 )->Z();
                double iniY = mpModel->GetNode( (*it_ES)[ iniEdg ]-1 )->Y();
                                                          
                double endX = mpModel->GetNode( (*it_ES)[ endEdg ]-1 )->Z();
                double endY = mpModel->GetNode( (*it_ES)[ endEdg ]-1 )->Y();
                
                double edgX = endX - iniX;
                double edgY = endY - iniY;

                double nodX = pNode->Z()           - iniX;
                double nodY = PBC_RotateY( pNode ) - iniY;

                double ModEdge = sqrt( edgX*edgX + edgY*edgY );
       
                double Project = ( nodX*edgX + nodY*edgY ) / ModEdge;
       
                if ( ( Project < 0.0 ) || ( Project > ModEdge ) ) 
                {
                    continue;
                }
       
                double Distance = abs( ( nodY*edgX - nodX*edgY  ) / ModEdge ); 
       
                if ( Distance > MinDistance ) 
                {
                    continue;
                }
       
                MinDistance = Distance;
       
                std::vector<int> EdgeNodeIds;
       
                EdgeNodeIds.push_back( (*it_ES)[ iniEdg ] );
                EdgeNodeIds.push_back( (*it_ES)[ endEdg ] );

                if ( mElementOrder == 2 ) 
                {
                    EdgeNodeIds.push_back( (*it_ES)[ iniEdg + 3 ] );
                }
              
                mPBC_NodeEdgePairs_Right[ NodeId ] = EdgeNodeIds;
                mPBC_NodeEdgeCoord_Right[ NodeId ] = Project / ModEdge;
            }
        }
    }

    //***********************************************************************************************************************************************
	// - Calculates mapping parameters (linear interpolation)
	//***********************************************************************************************************************************************
    void Modeler::PBC_Mapping_Parameters( std::vector<double>& a, std::vector<double>& b, std::vector<double>& eX, std::vector<double>& eY )
    {
        a.resize( 3 );
        b.resize( 3 );

        double cte = 1.0 / ( (eY[1]-eY[2]) * (eX[0]-eX[2]) - (eY[2]-eY[0]) * (eX[2]-eX[1]) );

        a[0] = cte * ( eX[2] * eY[0] - eY[2] * eX[0] );
        a[1] = cte * ( eY[2] - eY[0] );
        a[2] = cte * ( eX[0] - eX[2] );

        b[0] = cte * ( eX[0] * eY[1] - eY[0] * eX[1] );
        b[1] = cte * ( eY[0] - eY[1] );
        b[2] = cte * ( eX[1] - eX[0] );
    }

    //***********************************************************************************************************************************************
	// - Checks if the projected node with natural coordinates ( natX, natY ) is a node of the natural element.
    // - If the projected node projects onto a node it returns the local node Id.
    // - If the projected node is not any of the element nodes then it returns -1.
	//***********************************************************************************************************************************************
    int Modeler::PBC_Node_To_Node( double natX, double natY )
    {
        // Nodes on the vertices (valid for any element order)
        if ( ( abs( natX        ) < mPBCTolerance ) && ( abs( natY        ) < mPBCTolerance ) ) { return 0; }
        if ( ( abs( natX - 1.00 ) < mPBCTolerance ) && ( abs( natY        ) < mPBCTolerance ) ) { return 1; }
        if ( ( abs( natX        ) < mPBCTolerance ) && ( abs( natY - 1.00 ) < mPBCTolerance ) ) { return 2; }

        // 0th (3sb) order elements middle nodes
        if ( mElementOrder == 0 )
        {
            if ( ( abs( natX - 0.25 ) < mPBCTolerance ) && ( abs( natY - 0.25 ) < mPBCTolerance ) ) { return 3; }
            if ( ( abs( natX - 0.50 ) < mPBCTolerance ) && ( abs( natY - 0.25 ) < mPBCTolerance ) ) { return 4; }
            if ( ( abs( natX - 0.25 ) < mPBCTolerance ) && ( abs( natY - 0.50 ) < mPBCTolerance ) ) { return 5; }
        }
        // 2nd order elements middle nodes
        else if ( mElementOrder == 2 )
        {
            if ( ( abs( natX - 0.50 ) < mPBCTolerance ) && ( abs( natY        ) < mPBCTolerance ) ) { return 3; }
            if ( ( abs( natX - 0.50 ) < mPBCTolerance ) && ( abs( natY - 0.50 ) < mPBCTolerance ) ) { return 4; }
            if ( ( abs( natX        ) < mPBCTolerance ) && ( abs( natY - 0.50 ) < mPBCTolerance ) ) { return 5; }
        }

        return -1;
    }

    //***********************************************************************************************************************************************
	// - Checks if the projected node with natural coordinates ( natX, natY ) is on an edge of the element ElementNodeIds.
    // - If the projected node is on an edge then it returns its natural coordinate to respect that edge.
    // - If the projected node is not on an edge then it returns -1.0.
	//***********************************************************************************************************************************************
    double Modeler::PBC_Node_To_Edge( double natX, double natY, std::vector<int>& ElementNodeIds, std::vector<int>& EdgeNodeIds )
    {
        // Projected node is on the edge [0 1]
        if ( abs( natY ) < mPBCTolerance )
        {
            EdgeNodeIds.push_back( ElementNodeIds[0] );
            EdgeNodeIds.push_back( ElementNodeIds[1] );

            // 2nd order element middle nodes
            if ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back( ElementNodeIds[3] );
            }

            return natX;
        }
        // Projected node is on the edge [1 2]
        else if ( abs( natX + natY - 1.00 ) < mPBCTolerance )
        {
            EdgeNodeIds.push_back( ElementNodeIds[1] );
            EdgeNodeIds.push_back( ElementNodeIds[2] );

            // 2nd order element middle nodes
            if ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back( ElementNodeIds[4] );
            }

            return ( 1.00 - natX );
        }
        // Projected node is on the edge [2 0]
        else if ( abs( natX ) < mPBCTolerance )
        {
            EdgeNodeIds.push_back( ElementNodeIds[2] );
            EdgeNodeIds.push_back( ElementNodeIds[0] );

            // 2nd order element middle nodes
            if ( mElementOrder == 2 )
            {
                EdgeNodeIds.push_back( ElementNodeIds[5] );
            }

            return ( 1.00 - natY );
        }
        else
        {
            return -1.0;
        }
    }

	//***********************************************************************************************************************************************
    // - Assigns the nodes on the right periodic surface to the elements on the left periodic surface.
	// - Left surface must be placed in the 00-YZ-Plane.
	// - For cylindrical symmetry (cyclic periodic) the central axis must be placed along the Z axis.
    //***********************************************************************************************************************************************
	void Modeler::Set_RightLeft_PBC()
	{
         // Loop through materials in Left PBC
		std::map<unsigned int, std::vector<std::vector<int> > >::iterator it_Left;

        for( it_Left = mPBCLeftElements.begin(); it_Left != mPBCLeftElements.end(); ++it_Left )
        {
            // Material ID
			int pbcMaterial = it_Left->first;

            // Group of Left elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetLeft( it_Left->second );

            // Group of Right elements made of the same material (pbcMaterial)
			std::vector<std::vector<int> > elemSetRight( mPBCRightElements[ pbcMaterial ] );

            // List of nodes (sorted and unique) in the Right PBC
            std::list<int> RightNodesList;

            // Boolean map for the nodes in element vertices 
            std::map< unsigned int, bool > rNodeInVertice;

            // Filling the list with the Right elements nodes and checking if the Right surface is tilted
            std::vector<std::vector<int> >::iterator it_eR;

            for( it_eR = elemSetRight.begin(); it_eR != elemSetRight.end(); ++it_eR )
            {
                std::vector<int>::iterator it_nR;

                int iLocal = 0;

                for( it_nR = (*it_eR).begin(); it_nR != (*it_eR).end(); ++it_nR )
				{
                    RightNodesList.push_back( *it_nR );

                    if ( iLocal < 3 ) rNodeInVertice[ *it_nR ] = true ;
                    else              rNodeInVertice[ *it_nR ] = false;

                    iLocal++;
                }

                // Checks if the Right surface is tilted or not (Periodic/Cyclic boundary conditions)
                double x[3]; 
                
                for( int i=0; i<3; i++ )
                {
                    x[i] = ( mpModel->GetNode( (*it_eR)[i]-1 ) )->X();
                }

                mIsRightPBCTilted = ( abs( x[0] - x[1] ) > mPBCTolerance ) || 
                                    ( abs( x[1] - x[2] ) > mPBCTolerance ) || 
                                    ( abs( x[2] - x[0] ) > mPBCTolerance )  ;
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

                for( int i=0; i<(*it_eL).size(); i++ )
                {
                    eX.push_back( mpModel->GetNode( (*it_eL)[i]-1 )->Z() );
                    eY.push_back( mpModel->GetNode( (*it_eL)[i]-1 )->Y() );
                }

                // Calculating cross section (only nodes inside this area will be consider)
                double eCross_Section;

                double edge01_Lenght = sqrt( pow( eX[0]-eX[1], 2 ) + pow( eY[0]-eY[1], 2 ) );
                double edge02_Lenght = sqrt( pow( eX[0]-eX[2], 2 ) + pow( eY[0]-eY[2], 2 ) );

                if ( edge02_Lenght > edge01_Lenght ) 
                {
                    eCross_Section = ( 1.0 + mPBCTolerance ) * edge02_Lenght;
                }
                else
                {
                    eCross_Section = ( 1.0 + mPBCTolerance ) * edge01_Lenght;
                }

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
                    // Recovering right node pointer
                    Node::Pointer prNode = mpModel->GetNode( (*it_nRl)-1 );

                    // Cartesina coordinates of the node (2D projection)
                    double carX = prNode->Z();
                    double carY = PBC_RotateY( prNode );

                    // Checking if the Right node projection is close enough to the Left element
                    if ( sqrt( pow( eX[0]-carX, 2 ) + pow( eY[0]-carY, 2 ) ) > eCross_Section ) continue;

                    // Natural coordinates of the Right node projection on the Left element
                    double natX = a[0] + ( a[1] * carX ) + ( a[2] * carY );
                    double natY = b[0] + ( b[1] * carX ) + ( b[2] * carY );

                    // If the Right node is outside the Left element -> continue
                    if ( ( natX < -mPBCTolerance ) || ( natX > ( 1.0 + mPBCTolerance ) ) ) continue; 
                    if ( ( natY < -mPBCTolerance ) || ( natY > ( 1.0 + mPBCTolerance ) ) ) continue; 
                    if ( ( natX + natY )           >           ( 1.0 + mPBCTolerance )   ) continue;

                    // This node is going to be assigned
                    AssignedRightNodes.push_back( *it_nRl );

                    // If the Right surface is tilted and the node is in the Z-axis then make Ex=Ey=0 and continue
                    if ( ( mIsRightPBCTilted == true ) && ( abs( prNode->X() ) < mPBCTolerance ) && ( abs( prNode->Y() ) < mPBCTolerance ) )
                    {
                        if ( mNormals.find( *it_nRl ) == mNormals.end() )
                        {
                            std::complex<double> cZero( 0.00, 0.00 );

                            if ( mPotentials_On )
                            {
                                (*mpModel)( cAx, *prNode ) = cZero; prNode->pDofcAx()->FixDof();
                                (*mpModel)( cAy, *prNode ) = cZero; prNode->pDofcAy()->FixDof();
                            }                            
                            else                         
                            {                            
                                (*mpModel)( cEx, *prNode ) = cZero; prNode->pDofcEx()->FixDof();
                                (*mpModel)( cEy, *prNode ) = cZero; prNode->pDofcEy()->FixDof();                  
                            }                                                        
                        }

                        continue;
                    }

                    // Fix all DOFs of the assigned Right node
                    Fix_All_DOFs_In_Node( prNode );

                    // Cheking if the Right node is projected onto a Left node
                    int lLocalId = PBC_Node_To_Node( natX, natY );

                    // Only nodes with the same DOFs are assigned
                    if ( lLocalId > -1 )
                    {
                        if ( ( ( rNodeInVertice[ prNode->Id() ] == true  ) && ( lLocalId < 3 ) ) ||
                             ( ( rNodeInVertice[ prNode->Id() ] == false ) && ( lLocalId > 2 ) )  )
                        {
                            mPBC_NodeNodePairs_Right[ prNode->Id() ] = (*it_eL)[ lLocalId ]; 

                            continue;                         
                        }
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

                    // If Right node is not projected onto a node/edge then it is projected onto a Left element
                    std::vector<double> NatCoord_InElement( 2 );

                    NatCoord_InElement[0] = natX;
                    NatCoord_InElement[1] = natY;

                    mPBC_NodeElementPairs_Right[ prNode->Id() ] = (*it_eL);
                    mPBC_NodeElementCoord_Right[ prNode->Id() ] = NatCoord_InElement;

                }// End loop over the list of nodes in the Right PBC

                // Removing assigned Right nodes from the RightNodesList 
                std::vector<int>::iterator ita;

                for( ita = AssignedRightNodes.begin(); ita != AssignedRightNodes.end(); ++ita ) 
                {
                    RightNodesList.remove( *ita );
                }

            }// End loop through Left elements

            // Assign the remaining nodes to the closest Left element edge 
            std::list<int>::iterator it_nRl;

            for( it_nRl = RightNodesList.begin(); it_nRl != RightNodesList.end(); ++it_nRl )
            {
                Assign_Right_Node_To_Closest_Edge( *it_nRl, elemSetLeft );
            }

        }// End loop through materials in Left PBC
            
        // Clearing PBC element maps
		mPBCRightElements.clear();
		std::map< unsigned int, std::vector< std::vector<int> > >().swap( mPBCRightElements );
        
		mPBCLeftElements.clear();
		std::map< unsigned int, std::vector< std::vector<int> > >().swap( mPBCLeftElements );
	}

	//*****************************************************************************************************************
	// - Rotates Right node around Z-axis from Right plane to Left plane (00-YZ-plane)
	//*****************************************************************************************************************
	double Modeler::PBC_RotateY( Node::Pointer prNode )
	{
		if ( mIsRightPBCTilted == false ) return ( prNode->Y() );

		double X = prNode->X();
		double Y = prNode->Y();

		double modXY2 = X*X + Y*Y;

		return ( modXY2 / sqrt( modXY2 ) );
	}

	//*****************************************************************************************************************
    //* - Calculates the volume integral of the fields
    //*****************************************************************************************************************
    void Modeler::GenerateVolIntElement( int* NodesId, unsigned int PropertiesId )
    {
		std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        if ( mElementOrder == 0 ) nodes.resize( 4 );

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double volumeID = cIBC[0];

		std::vector<double> cX;
		std::vector<double>	cY;
		std::vector<double> cZ;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

		Matrix<double> N;

        if      ( mElementOrder == 0 ) Lagrange3D_Ni_1st( N, cX, cY, cZ );
		else if ( mElementOrder == 1 ) Lagrange3D_Ni_1st( N, cX, cY, cZ );
        else if ( mElementOrder == 2 ) Lagrange3D_Ni_2nd( N, cX, cY, cZ );

		double elementVolume = 0.00;

		if ( mQuadraticGeometry )
		{
			std::vector<double> detJ;

			Calculate_detJ( detJ, cX, cY, cZ, nodes );

			for ( int gp=0; gp<nGaussPoints; ++gp ) W[gp]         *= detJ[gp];
			for ( int gp=0; gp<nGaussPoints; ++gp ) elementVolume +=    W[gp];
		}
		else
		{
			elementVolume = VolumeOfTetrahedra(nodes);

			double jacob  = 6.00 * elementVolume;

			for ( int gp=0; gp<nGaussPoints; ++gp ) W[gp] *= jacob;
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

            for ( int gp=0; gp<nGaussPoints; ++gp ) vIntNi += W[gp] * N[in][gp];

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
    void Modeler::NaturalDerivatives3D_2ndQ( Matrix<double>&      dNda, 
                                             Matrix<double>&      dNdb, 
                                             Matrix<double>&      dNdu, 
                                             std::vector<double>&   cX, 
                                             std::vector<double>&   cY, 
                                             std::vector<double>&   cZ )
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
    void Modeler::Calculate_detJ( std::vector<double>& detJ, 
                                  std::vector<double>&   cX, 
                                  std::vector<double>&   cY, 
                                  std::vector<double>&   cZ, 
                                  std::vector<Node::Pointer>& cnodes )
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
    double Modeler::VolumeOfTetrahedra( std::vector<Node::Pointer>& cnodes )
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
    double Modeler::AreaOfTriangle( std::vector<Node::Pointer>& cnodes )
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

    //***********************************************************************************
    //* - Calculation of the Jacobian determinant
    //***********************************************************************************
    void Modeler::Calculate_detJ( std::vector<double>& detJ, 
                                  std::vector<double>&   cX, 
                                  std::vector<double>&   cY, 
                                  std::vector<Node::Pointer>& cnodes )
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
    void Modeler::GenerateSurfIntElement( int* NodesId, unsigned int PropertiesId )
    {
		std::vector<Node::Pointer> nodes;

        Get_Surface_Element_Nodes( NodesId, nodes );

        if ( mElementOrder == 0 ) nodes.resize( 3 );

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*properties)(COMPLEX_IBC);

		double surfaceID = cIBC[0];

		int gp;

		std::vector<double> cX;
		std::vector<double>	cY;
		std::vector<double>  W;

		int nGaussPoints = GaussPoints2D_Order12( cX, cY, W );

		Matrix<double> N;

        if      ( mElementOrder == 0 ) Lagrange2D_Ni_1st( N, cX, cY );
		else if ( mElementOrder == 1 ) Lagrange2D_Ni_1st( N, cX, cY );
        else if ( mElementOrder == 2 ) Lagrange2D_Ni_2nd( N, cX, cY );

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

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the volumetric element nodes
    //**********************************************************************************************************************
    void Modeler::Get_Volume_Element_Nodes( int* NodesId, std::vector<Node::Pointer>& pNodes )
    {
        std::vector<int> HONodesId;

        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the volumetric element nodes and the high-order nodeIds vector
    //**********************************************************************************************************************
    void Modeler::Get_Volume_Element_Nodes( int* NodesId, std::vector<int>& HONodesId, std::vector<Node::Pointer>& pNodes )
    {
        PushHONodesOnVolume( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the surface element nodes 
    //**********************************************************************************************************************
    void Modeler::Get_Surface_Element_Nodes( int* NodesId, std::vector<Node::Pointer>& pNodes )
    {
        std::vector<int> HONodesId;

        PushHONodesOnSurface( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the surface element nodes and the high-order nodeIds vector
    //**********************************************************************************************************************
    void Modeler::Get_Surface_Element_Nodes( int* NodesId, std::vector<int>& HONodesId, std::vector<Node::Pointer>& pNodes )
    {
        PushHONodesOnSurface( NodesId, HONodesId );

        std::vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1) );
        }
    }

    //*****************************************************************************************************************
    //* - Calculates volumetric integrals of the fields in electrostatic mode.
    //*****************************************************************************************************************
    void Modeler::GenerateVolIntElement_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId = 0 then the element is used for surface integrals.
        if ( PropertiesId == 0 )
        {
            GetGrad_forSurfInt_SmoothingOff_Electrostatic( NodesId ); 
            return;
        }

        // Getting all element nodes
        std::vector<Node::Pointer> nodes; 
        
        Get_Volume_Element_Nodes( NodesId, nodes );

        // Getting Volume Id
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Vector<double> cIBC = (*properties)( COMPLEX_IBC );

        double volumeID = cIBC[0];

        // Integration Gauss points
        int gp; std::vector<double> cX, cY, cZ, W;

        int nGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

        Matrix<double> N;

        Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // Volume of the element
        double elementVolume = VolumeOfTetrahedra( nodes );
        double jacob         = 6.00 * elementVolume;

        // Volumetric integral weights
        for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob;

        // Volumetric integrals vector
        Vector<double> volIntegralE( 5, 0.00 );

        // Electrostatic electric field on nodes
        Vector< Vector<double> > EEfield_OnNodes( nodes.size() );

        if ( mSmoothing_On == true )
        { 
            for ( int i=0; i<nodes.size(); i++ ) EEfield_OnNodes[i] = (*properties)( ELECTROSTATIC_ELECTRIC_FIELD, *nodes[i] );
        }
        else
        {
            Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );
            
            pElement->Calculate_Gradient_OnNodes( EEfield_OnNodes );  

            for ( int i=0; i<nodes.size(); i++ ) EEfield_OnNodes[i] *= -1;
        }

        // Integration on volume over nodes
        for ( int i=0; i<nodes.size(); i++ )
        {
            // Volume integral
            double vIntNi = 0.00; 
            
            for ( gp=0; gp<nGaussPoints; gp++ ) vIntNi += W[gp] * N[i][gp];

             // Electrostatic electric field on nodes
            Vector<double> EEfield = EEfield_OnNodes[i];
            
            double modE2 = EEfield[0]*EEfield[0] + EEfield[1]*EEfield[1] + EEfield[2]*EEfield[2];
            double modE  = sqrt( modE2 );

            // Volume integral of Ex, Ey, Ez
            volIntegralE[0] += EEfield[0] * vIntNi;
            volIntegralE[1] += EEfield[1] * vIntNi;
            volIntegralE[2] += EEfield[2] * vIntNi;

            // Volume integral of |E|^2, |E|
            volIntegralE[3] += modE2 * vIntNi;
            volIntegralE[4] += modE  * vIntNi;
        }

        // Adding the volume integrals of this element to the total
        if ( mVolumeValue.find(volumeID) != mVolumeValue.end() )
        {
            mVolumeValue               [ volumeID ] += elementVolume;
            mVolIntegralE_Electrostatic[ volumeID ] += volIntegralE ;
        }                                         
        else                                      
        {                                         
            mVolumeValue               [ volumeID ] = elementVolume;
            mVolIntegralE_Electrostatic[ volumeID ] = volIntegralE ;
        }
    }

    //*****************************************************************************************************************
    //* - Calculates gradients for the non-smoothed fields surface integrals.
    //*   If PropertiesId = 0 then the PVIE element is used for getting field gradients.
    //*   These gradients are used in the surface integrals of non-smoothed fields.
    //*****************************************************************************************************************
    void Modeler::GetGrad_forSurfInt_SmoothingOff_Electrostatic( int* NodesId )
    {
        // Getting all element nodes
        std::vector<Node::Pointer> nodes; Get_Volume_Element_Nodes( NodesId, nodes );

        // Dummy property, no materials are used for the gradients
        Properties::Pointer properties = mpModel->GetProperties( 1 );

        // Setting element
        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );

        // Calculate gradients on nodes (non-smoothed gradients are constant in 1st order elements)
        Vector< Vector<double> > gradOnNodes;

        pElement->Calculate_Gradient_OnNodes( gradOnNodes );

        // E = -grad( V )
        Vector<double> Efield_OnNodes = gradOnNodes[0];

        Efield_OnNodes *= -1;

        // Adds gradient on Gauss point to the mSurfInt_NonSmooth_E_Electrostatic vector
        mSurfInt_NonSmooth_E_Electrostatic.push_back( Efield_OnNodes );
    }

    //*****************************************************************************************************************
    //* - Calculates surface integral of the fields for the electrostatic case.
    //*****************************************************************************************************************
    void Modeler::GenerateSurfIntElement_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // Getting all element nodes
        std::vector<Node::Pointer> nodes;

        Get_Surface_Element_Nodes( NodesId, nodes );

        // Getting Surface_Id
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Vector<double> cIBC = (*properties)(COMPLEX_IBC);

        double surfaceID = cIBC[0];

        // Integration Gauss points
        int gp; std::vector<double> cX, cY, W;

        int nGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        Matrix<double> N;

        Lagrange2D_Ni_1st( N, cX, cY );

        // Surface of the element
        double elementArea = AreaOfTriangle( nodes );
        double jacob       = 2.00 * elementArea;

        for ( gp=0; gp<nGaussPoints; gp++ ) W[gp] *= jacob;

        // Surface integrals vectors
        Vector<double> surfIntegralE( 5, 0.00 );

        // Electrostatic electric field on nodes 
        Vector< Vector<double> > EEfield_OnNodes( nodes.size() );

        if ( mSmoothing_On == true )
        { 
            for ( int i=0; i<nodes.size(); i++ ) EEfield_OnNodes[i] = (*properties)(ELECTROSTATIC_ELECTRIC_FIELD, *nodes[i]);
        }
        else
        {
            for ( int i=0; i<nodes.size(); i++ ) EEfield_OnNodes[i] = mSurfInt_NonSmooth_E_Electrostatic[ mNonSmooth_Surf_Counter ];
        }

        // Integration on volume over nodes
        for ( int i=0; i<nodes.size(); i++ )
        {
            // Nodal surface integral
            double sIntNi = 0.00;

            for ( gp=0; gp<nGaussPoints; gp++ ) sIntNi += W[gp]*N[i][gp];

            // Electrostatic electric field on node
            Vector<double> EEfield = EEfield_OnNodes[i];

            double modE2 = EEfield[0]*EEfield[0] + EEfield[1]*EEfield[1] + EEfield[2]*EEfield[2];
            double modE  = sqrt( modE2 );

            // Surface integral of Ex, Ey, Ez
            surfIntegralE[0] += EEfield[0] * sIntNi;
            surfIntegralE[1] += EEfield[1] * sIntNi;
            surfIntegralE[2] += EEfield[2] * sIntNi;

            // Surface integral of |E|^2, |E|
            surfIntegralE[3] += modE2 * sIntNi;
            surfIntegralE[4] += modE  * sIntNi;
        }

        // Adding surface integrals of this element to the total
        if ( mSurfaceValue.find( surfaceID ) != mSurfaceValue.end() )
        {
            mSurfaceValue               [ surfaceID ] += elementArea;
            mSurfIntegralE_Electrostatic[ surfaceID ] += surfIntegralE;
        }
        else
        {
            mSurfaceValue               [ surfaceID ] = elementArea;
            mSurfIntegralE_Electrostatic[ surfaceID ] = surfIntegralE;
        }

        // Increase global non-smoothing surface_id counter
        mNonSmooth_Surf_Counter++;
    }

	//*******************************************************************************************************
    //* - Projecting resulting field on plane
    //*******************************************************************************************************
    void Modeler::GenerateProjectionRWPortTE10( int*  NodesId, unsigned int PropertiesId )
    {
        std::vector<Node::Pointer> nodes;

        Get_Surface_Element_Nodes( NodesId, nodes );

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new RWPortTE10_3sb_FullWave ( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new RWPortTE10_1st_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new RWPortTE10_2ndQ_FullWave( nodes, properties ) );

        double portNumber = pElement->GetPortNumber();

		std::complex<double> TE10Proyec    = pElement->GetProjection   ();
		std::complex<double> Normalization = pElement->GetNormalization();

		// Elemental projection
		if ( mProjection.find(portNumber) != mProjection.end() ) mProjection[portNumber] += TE10Proyec;
        else                                                     mProjection[portNumber]  = TE10Proyec;

		// Elemental normalization
		if ( mNormalization.find(portNumber) != mNormalization.end() ) mNormalization[portNumber] += Normalization;
        else                                                           mNormalization[portNumber]  = Normalization;
    }

	//***********************************************************************************************************
    //* - Projecting resulting field on plane
    //***********************************************************************************************************
    void Modeler::GenerateProjectionCoaxPortTEM( int*  NodesId, unsigned int PropertiesId )
    {
		std::vector<Node::Pointer> nodes;

        Get_Surface_Element_Nodes( NodesId, nodes );

		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new CoaxialPortTEM_3sb_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st_FullWave ( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new CoaxialPortTEM_2ndQ_FullWave( nodes, properties ) );
		
        double portNumber = pElement->GetPortNumber();

		std::complex<double> TEMProyec     = pElement->GetProjection   ();
		std::complex<double> Normalization = pElement->GetNormalization();

		// Elemental projection
		if ( mProjection.find(portNumber) != mProjection.end() ) mProjection[portNumber] += TEMProyec;
        else                                                     mProjection[portNumber]  = TEMProyec;

		// Elemental normalization
		if ( mNormalization.find(portNumber) != mNormalization.end() ) mNormalization[portNumber] += Normalization;
        else                                                           mNormalization[portNumber]  = Normalization;
    }

    //*************************************************************************************************
    //* - Calculates average normal in PEC nodes
    //*************************************************************************************************
    void Modeler::GeneratePECElement( std::vector<int>& LONodesId )
    {
        if ( mElectrostaticMode ) return;
 
        std::vector<int> NodesId;

		PushHONodesOnSurface( LONodesId, NodesId );

		std::vector<int>::const_iterator id_it;

        std::vector<Node::Pointer> nodes;

        Node::Pointer pNode;
        
        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );
            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PEC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetPECnormals[*id_it].push_back( normal );
		}
	}

    //*************************************************************************************************
    //* - Calculates average normal in Plasma Sheath Boundary Condition (PSBC) nodes
    //*************************************************************************************************
    void Modeler::GeneratePSBCElement( std::vector<int>& LONodesId )
    {
        if ( mColdPlasmaMode == false ) return;

        std::vector<int> NodesId;

		PushHONodesOnSurface( LONodesId, NodesId );

		std::vector<int>::const_iterator id_it;

        std::vector<Node::Pointer> nodes;

        Node::Pointer pNode;

        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );
            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PSBC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetPSBCnormals[*id_it].push_back( normal );
		}
    }

    //*************************************************************************************************
    //* - Calculates normal in PMC nodes (symmetry planes)
    //*************************************************************************************************
    void Modeler::GeneratePMCElement( std::vector<int>& LONodesId )
    {
        if ( mElectrostaticMode ) return;

        if ( mColdPlasmaMode )
        {
            String sEparTol = mpColdPlasma->Get_Eparallel_Tolerance();

            if ( sEparTol != "Off" )
            {
                Send_Error_Msg( "PMC boundary condition incompatible with active E par tolerance",
                                "Please, eliminate PMC or change E par tolerance to Off.", 1 );
            }
        }

	    std::vector<int> NodesId;

		PushHONodesOnSurface( LONodesId, NodesId );

		std::vector<int>::const_iterator id_it;

        std::vector<Node::Pointer> nodes;

        Node::Pointer pNode;

        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );
            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the PMC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetPMCnormals[*id_it].push_back( normal );
		}
    }

	//*************************************************************************************************
    //* - Calculates average normal in TE PMC nodes (symmetry plane ortogonal to a PMC plane)
    //*************************************************************************************************
    void Modeler::GenerateTEPMCElement( std::vector<int>& LONodesId )
    {
        if ( mElectrostaticMode == true ) return;

        std::vector<int> NodesId;

		PushHONodesOnSurface( LONodesId, NodesId );

		std::vector<int>::const_iterator id_it;

        std::vector<Node::Pointer> nodes;

        Node::Pointer pNode;

        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );
            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z()) -
                    (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y());

        normal[1] = (nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X()) -
                    (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z());

        normal[2] = (nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y()) -
                    (nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X());

        // Pushing back the normal vector of the element to the TEPMC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetTEPMCnormals[*id_it].push_back( normal );
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

	        PrintSijParameters ( mProjectionNormalized, mProblemFrequency );
            WriteSijHeaderFiles( mProjectionNormalized );                  
            WriteSijInFiles    ( mProjectionNormalized, mProblemFrequency );

		    mProjectionNormalized.clear();
		    std::map<double, std::complex<double> >().swap( mProjectionNormalized );
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
	        PrintViParameters ( mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
				                mVolumeValue    , mProblemFrequency );
                                
            WriteViHeaderFiles( mVolumeIntegralE );
                                
            WriteViInFiles    ( mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
						        mProblemFrequency );

		    mVolumeIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralE );
                                                                                      
			mVolumeIntegralH.clear();                                                 
            std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralH );
                                                                                      
			mVolumeIntegralJ.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralJ );

			mVolumeValue.clear();
			std::map<double, double>().swap( mVolumeValue );
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters ( mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                   mSurfaceValue    , mProblemFrequency );
                                   
            WriteSurfiHeaderFiles( mSurfaceIntegralE );
                                   
            WriteSurfiInFiles    ( mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                   mProblemFrequency );

		    mSurfaceIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralE );
                                                                                       
			mSurfaceIntegralH.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralH );
                                                                                       
			mSurfaceIntegralJ.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralJ );

			mSurfaceValue.clear();
			std::map<double, double>().swap( mSurfaceValue );
        }

        b_vector.resize( mSystemSize, 0.00 );
        A_matrix.Resize( mSystemSize );

        if ( mA_matrix_aux_Required == true )
        {
            A_matrix_aux.Resize( mSystemSize );
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

	        PrintSijParameters( mProjectionNormalized, mProblemFrequency );
            WriteSijInFiles   ( mProjectionNormalized, mProblemFrequency );

		    mProjectionNormalized.clear();
		    std::map<double, std::complex<double> >().swap( mProjectionNormalized );
        }

		if ( mVolumeIntegralE.size() > 0 )
        {
            PrintViParameters( mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
				               mVolumeValue    , mProblemFrequency );
                               
            WriteViInFiles   ( mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
							   mProblemFrequency );

		    mVolumeIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralE );
                                                                                      
			mVolumeIntegralH.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralH );
                                                                                      
			mVolumeIntegralJ.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mVolumeIntegralJ );

			mVolumeValue.clear();
			std::map<double, double>().swap( mVolumeValue );
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters( mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                  mSurfaceValue    , mProblemFrequency );

            WriteSurfiInFiles( mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				               mProblemFrequency );

		    mSurfaceIntegralE.clear();
		    std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralE );
                                                                                       
			mSurfaceIntegralH.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralH );
                                                                                       
			mSurfaceIntegralJ.clear();                                                 
			std::map<double, Vector<std::complex<double> > >().swap( mSurfaceIntegralJ );

			mSurfaceValue.clear();
			std::map<double, double>().swap( mSurfaceValue );
        }

        b_vector.resize( mSystemSize, 0.00 );
        A_matrix.Resize( mSystemSize );

        if ( mA_matrix_aux_Required == true )
        {
            A_matrix_aux.Resize( mSystemSize );
        }
	}

    //********************************************************************************************************************************
    //* - Solving linear problem ( full wave mode )
    //********************************************************************************************************************************
    void Modeler::SolveProblem_FullWave()
    {
        // Solving linear system Ax=b
        SolveLinearSystem();

        // Saving solution vector on Variables 
		if ( mPotentials_On )
        {
            SaveSystemSolution_FullWave( cAx );
            SaveSystemSolution_FullWave( cAy );
            SaveSystemSolution_FullWave( cAz );    
            SaveSystemSolution_FullWave( cVs );    
        }
        else
        {
            SaveSystemSolution_FullWave( cEx );
            SaveSystemSolution_FullWave( cEy );
            SaveSystemSolution_FullWave( cEz );        
        }

        // Updating fixed DOFs
        if ( mColdPlasmaMode ) UpdateCoord_PlasmaRLP ();
        if ( mAxisymmetric   ) UpdateAxisToCartesian();

		UpdateCoord_EPEC   ();
		UpdateCoord_Contact();
		UpdateCoord_EPBC   ();
    }

    //********************************************************************************************************************************
    //* - Solving linear problem ( electrostatic mode )
    //********************************************************************************************************************************
    void Modeler::SolveProblem_Electrostatic()
    {
        // Solving linear system Ax=b
        SolveLinearSystem();

        // Saving solution vector on VOLTAGE variable
        SaveSystemSolution_Electrostatic( VOLTAGE );

        // Updating fixed DOFs
        UpdateFixedVOLTAGE();
    }

    //********************************************************************************************************************************
    //* - Updating fixed VOLTAGE DOFs
    //********************************************************************************************************************************
    void Modeler::UpdateFixedVOLTAGE()
    {
        // Updating fixed DOFs
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            if ( (*it)->pDofV()->IsFixed() )
            {
                (*mpModel)( VOLTAGE, **it ) = mFixVoltage[ (*it)->Id() ];
            }
        }
    }

    //********************************************************************************************************************************
    //* - Updates coordinate system in plasma mode
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_PlasmaRLP()
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

	//********************************************************************************************************************************
    //* - Updates coordinate system in axisymmetric problems
    //********************************************************************************************************************************
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
                Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( nodeId, n, t, b );

                std::complex<double> En = (*mpModel)(cEx, **nodes_it),
					                 Et = (*mpModel)(cEy, **nodes_it),
                                     Eb = (*mpModel)(cEz, **nodes_it);

                (*mpModel)(cEx, **nodes_it) = n[0]*En + t[0]*Et + b[0]*Eb;
                (*mpModel)(cEy, **nodes_it) = n[1]*En + t[1]*Et + b[1]*Eb;
                (*mpModel)(cEz, **nodes_it) = n[2]*En + t[2]*Et + b[2]*Eb;
            }
        }
    }

	//********************************************************************************************************************************
    //* - Change local coordinates in PEC for cartesian coordinates 
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_EPEC()
    {
        if ( ( mNormals.size() == 0 ) || ( mSetAllExyToZero ) ) 
        {
            return;
        }

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator nodes_it;

        for ( nodes_it = nodes.begin(); nodes_it != nodes.end(); ++nodes_it )
        {
            if ( mNormals.find( (*nodes_it)->Id() ) == mNormals.end() ) 
            {
                continue;
            }

            // From D_ntb to E_ntb / A_tnb ( plasma sheath BC )
            if ( ( mColdPlasmaMode == true ) && ( mType_Of_BC_Normal[ (*nodes_it)->Id() ] == 'D' ) )
            {
                Matrix< std::complex<double> > RiEpRt;

                Get_En_R_invEp_Rt_Dn( RiEpRt, (*nodes_it)->Id() );

                if ( mPotentials_On )
                {
                    std::complex<double> Dn = (*mpModel)( cAx, **nodes_it );
                    std::complex<double> Dt = (*mpModel)( cAy, **nodes_it );
                    std::complex<double> Db = (*mpModel)( cAz, **nodes_it );
                    
                    // {A_ntb} = [R * invEp * Rt] * {D_ntb}
                    (*mpModel)( cAx, **nodes_it ) = RiEpRt[0][0]*Dn + RiEpRt[0][1]*Dt + RiEpRt[0][2]*Db;
                    (*mpModel)( cAy, **nodes_it ) = RiEpRt[1][0]*Dn + RiEpRt[1][1]*Dt + RiEpRt[1][2]*Db;
                    (*mpModel)( cAz, **nodes_it ) = RiEpRt[2][0]*Dn + RiEpRt[2][1]*Dt + RiEpRt[2][2]*Db;
                }
                else
                {
                    std::complex<double> Dn = (*mpModel)( cEx, **nodes_it );
                    std::complex<double> Dt = (*mpModel)( cEy, **nodes_it );
                    std::complex<double> Db = (*mpModel)( cEz, **nodes_it );
                    
                    // {E_ntb} = [R * invEp * Rt] * {D_ntb}
                    (*mpModel)( cEx, **nodes_it ) = RiEpRt[0][0]*Dn + RiEpRt[0][1]*Dt + RiEpRt[0][2]*Db;
                    (*mpModel)( cEy, **nodes_it ) = RiEpRt[1][0]*Dn + RiEpRt[1][1]*Dt + RiEpRt[1][2]*Db;
                    (*mpModel)( cEz, **nodes_it ) = RiEpRt[2][0]*Dn + RiEpRt[2][1]*Dt + RiEpRt[2][2]*Db;                
                }
            }

            // From E_ntb / A_ntb to E_xyz / A_xyz ( PEC / PMC BC )
            Vector<double> n( mNormals[ (*nodes_it)->Id() ] );
			Vector<double> t( 3 ); 
            Vector<double> b( 3 );

            TangencialCoordinates( n, t, b );

            if ( mPotentials_On )
            {
                std::complex<double> An = (*mpModel)( cAx, **nodes_it );
			    std::complex<double> At = (*mpModel)( cAy, **nodes_it );
                std::complex<double> Ab = (*mpModel)( cAz, **nodes_it );
                
                // {A_xyz} = [Rt] * {A_ntb}
                (*mpModel)( cAx, **nodes_it ) = n[0]*An + t[0]*At + b[0]*Ab;
                (*mpModel)( cAy, **nodes_it ) = n[1]*An + t[1]*At + b[1]*Ab;
                (*mpModel)( cAz, **nodes_it ) = n[2]*An + t[2]*At + b[2]*Ab;            
            }
            else
            {
                std::complex<double> En = (*mpModel)( cEx, **nodes_it );
			    std::complex<double> Et = (*mpModel)( cEy, **nodes_it );
                std::complex<double> Eb = (*mpModel)( cEz, **nodes_it );
                
                // {E_xyz} = [Rt] * {E_ntb}
                (*mpModel)( cEx, **nodes_it ) = n[0]*En + t[0]*Et + b[0]*Eb;
                (*mpModel)( cEy, **nodes_it ) = n[1]*En + t[1]*Et + b[1]*Eb;
                (*mpModel)( cEz, **nodes_it ) = n[2]*En + t[2]*Et + b[2]*Eb;            
            }
        }
    }

	//********************************************************************************************************************************
    //* - Update electric field values in PBC nodes
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_EPBC()
	{
		std::vector<int> FrontLeftNodesId;

		UpdateCoord_RL_PBC( FrontLeftNodesId );
		UpdateCoord_FB_PBC(                  );
		UpdateCoord_LF_PBC( FrontLeftNodesId );
	}

	//********************************************************************************************************************************
    //* - Check if the left node is also a front node
    //********************************************************************************************************************************
	bool Modeler::PBC_LeftNodeIsFrontNode( int rightNodeId, std::vector<int>& leftNodesId, std::vector<int>& FrontLeftNodesId )
	{
		bool LeftNodeIsFrontNode = false;

		for ( int i=0; i<leftNodesId.size(); i++ )
		{
			if ( ( mPBC_NodeNodePairs_Front.find   ( leftNodesId[ i ] ) != mPBC_NodeNodePairs_Front.end   () ) ||
				 ( mPBC_NodeEdgePairs_Front.find   ( leftNodesId[ i ] ) != mPBC_NodeEdgePairs_Front.end   () ) ||
				 ( mPBC_NodeElementPairs_Front.find( leftNodesId[ i ] ) != mPBC_NodeElementPairs_Front.end() )  )
			{
				LeftNodeIsFrontNode = true;

				FrontLeftNodesId.push_back( rightNodeId );

				break;
		    }
		}

		return LeftNodeIsFrontNode;
	}

	//********************************************************************************************************************************
    //* - Update electric field values in PBC nodes ( Right-Left faces )
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_RL_PBC( std::vector<int>& FrontLeftNodesId )
	{
		if ( ( mPBC_NodeNodePairs_Right.size() + mPBC_NodeEdgePairs_Right.size() + mPBC_NodeElementPairs_Right.size() ) == 0 )
        {
            return;
        }

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator nit;

		for ( nit=nodes.begin(); nit!=nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();

            std::complex<double> lcFx( 0.00, 0.00 );
			std::complex<double> lcFy( 0.00, 0.00 );
			std::complex<double> lcFz( 0.00, 0.00 );
            std::complex<double> lcVs( 0.00, 0.00 );

			if ( mPBC_NodeNodePairs_Right.find( NodeId ) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNodeId = mPBC_NodeNodePairs_Right[ NodeId ];

				std::vector<int> leftNodesId( 1 ); 
                
                leftNodesId[0] = lNodeId;

				if ( PBC_LeftNodeIsFrontNode( NodeId, leftNodesId, FrontLeftNodesId ) == true ) 
                {
                    continue;
                }

		        Node::Pointer plNode = mpModel->GetNode( lNodeId-1 );

                if ( mPotentials_On )
                {
				    lcFx = (*mpModel)( cAx, *plNode );
				    lcFy = (*mpModel)( cAy, *plNode );
				    lcFz = (*mpModel)( cAz, *plNode ); 

                    if ( plNode->pDof( cVs ) != NULL ) 
                    {
                        lcVs = (*mpModel)( cVs, *plNode ); 
                    }
                }
                else
                {
				    lcFx = (*mpModel)( cEx, *plNode );
				    lcFy = (*mpModel)( cEy, *plNode );
				    lcFz = (*mpModel)( cEz, *plNode );                
                }
            }
			else if ( mPBC_NodeEdgePairs_Right.find( NodeId ) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<int> lNodes( mPBC_NodeEdgePairs_Right[ NodeId ] );

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) 
                {
                    continue;
                }

                double natCoord = mPBC_NodeEdgeCoord_Right[ NodeId ];

		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i=0; i<lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[ i ]-1 );

                    if ( mPotentials_On )
                    {
                        lcFx += N[i] * (*mpModel)( cAx, *plNode );
				        lcFy += N[i] * (*mpModel)( cAy, *plNode );
				        lcFz += N[i] * (*mpModel)( cAz, *plNode );

                        if ( ( (*nit)->pDof( cVs ) != NULL ) && ( i < 2 ) )
                        {
                            lcVs += N[i] * (*mpModel)( cVs, *plNode );
                        }
                    }
                    else
                    {
                        lcFx += N[i] * (*mpModel)( cEx, *plNode );
				        lcFy += N[i] * (*mpModel)( cEy, *plNode );
				        lcFz += N[i] * (*mpModel)( cEz, *plNode );                    
                    }
                }
			}
			else if ( mPBC_NodeElementPairs_Right.find( NodeId ) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<int> lNodes( mPBC_NodeElementPairs_Right[ NodeId ] );

				if ( PBC_LeftNodeIsFrontNode( NodeId, lNodes, FrontLeftNodesId ) == true ) 
                {
                    continue;
                }

		        std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );

                std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i=0; i<lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    if ( mPotentials_On )
                    {
                        lcFx += N[i] * (*mpModel)( cAx, *plNode );
				        lcFy += N[i] * (*mpModel)( cAy, *plNode );
				        lcFz += N[i] * (*mpModel)( cAz, *plNode );

                        if ( ( (*nit)->pDof( cVs ) != NULL ) && ( i < 3 ) )
                        {
                            lcVs += N[i] * (*mpModel)( cVs, *plNode );
                        }
                    }
                    else
                    {
                        lcFx += N[i] * (*mpModel)( cEx, *plNode );
				        lcFy += N[i] * (*mpModel)( cEy, *plNode );
				        lcFz += N[i] * (*mpModel)( cEz, *plNode );                    
                    }
                }
			}
            else
            {
                continue;
            }

            double cos_A = 1.0;
			double sin_A = 0.0;

			if ( mIsRightPBCTilted == true ) 
            {
                PBC_InvRot( NodeId, cos_A, sin_A );
            }

            if ( mPotentials_On )
            {
			    (*mpModel)( cAx, **nit ) =  cos_A * lcFx + sin_A * lcFy;
			    (*mpModel)( cAy, **nit ) = -sin_A * lcFx + cos_A * lcFy;
			    (*mpModel)( cAz, **nit ) =  lcFz; 

                if ( (*nit)->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, **nit ) = lcVs; 
                }
            }
            else
            {
			    (*mpModel)( cEx, **nit ) =  cos_A * lcFx + sin_A * lcFy;
			    (*mpModel)( cEy, **nit ) = -sin_A * lcFx + cos_A * lcFy;
			    (*mpModel)( cEz, **nit ) =  lcFz;                
            }
		}
	}

	//********************************************************************************************************************************
    //* - Computes inverse rotation of a point
    //********************************************************************************************************************************
    void Modeler::PBC_InvRot( int NodeId, double& cos_A, double& sin_A )
    {
        double X = ( mpModel->GetNode( NodeId - 1 ) )->X();
		double Y = ( mpModel->GetNode( NodeId - 1 ) )->Y();

		double mod_r = sqrt(X*X + Y*Y);

		cos_A = Y / mod_r;
		sin_A = X / mod_r;
	}

	//********************************************************************************************************************************
    //* - Updates field values in PBC nodes ( Front-Back surfaces )
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_FB_PBC()
	{
		if ( ( mPBC_NodeNodePairs_Front.size() + mPBC_NodeEdgePairs_Front.size() + mPBC_NodeElementPairs_Front.size() ) == 0 )
        {
            return;
        }

		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		std::vector<Node::Pointer>::iterator nit;

		for ( nit=nodes.begin(); nit!=nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();

            std::complex<double> bcFx( 0.00, 0.00 );
			std::complex<double> bcFy( 0.00, 0.00 );
			std::complex<double> bcFz( 0.00, 0.00 );
            std::complex<double> bcVs( 0.00, 0.00 );

			if ( mPBC_NodeNodePairs_Front.find( NodeId ) != mPBC_NodeNodePairs_Front.end() )
            {
				int bNodeId = mPBC_NodeNodePairs_Front[ NodeId ];

		        Node::Pointer pbNode = mpModel->GetNode( bNodeId-1 );

                if ( mPotentials_On )
                {
				    bcFx = (*mpModel)( cAx, *pbNode );
				    bcFy = (*mpModel)( cAy, *pbNode );
				    bcFz = (*mpModel)( cAz, *pbNode ); 

                    if ( pbNode->pDof( cVs ) != NULL ) 
                    {
                        bcVs = (*mpModel)( cVs, *pbNode ); 
                    }
                }
                else
                {
				    bcFx = (*mpModel)( cEx, *pbNode );
				    bcFy = (*mpModel)( cEy, *pbNode );
				    bcFz = (*mpModel)( cEz, *pbNode );                
                }
            }
			else if ( mPBC_NodeEdgePairs_Front.find( NodeId ) != mPBC_NodeEdgePairs_Front.end() )
			{
				std::vector<int>    bNodes   ( mPBC_NodeEdgePairs_Front[ NodeId ] );
                double              natCoord = mPBC_NodeEdgeCoord_Front[ NodeId ]  ;
		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i=0; i<bNodes.size(); i++ )
                {
					Node::Pointer pbNode = mpModel->GetNode( bNodes[ i ]-1 );

                    if ( mPotentials_On )
                    {
                        bcFx += N[i] * (*mpModel)( cAx, *pbNode );
				        bcFy += N[i] * (*mpModel)( cAy, *pbNode );
				        bcFz += N[i] * (*mpModel)( cAz, *pbNode );

                        if ( ( (*nit)->pDof( cVs ) != NULL ) && ( i < 2 ) )
                        {
                            bcVs += N[i] * (*mpModel)( cVs, *pbNode );
                        }
                    }
                    else
                    {
                        bcFx += N[i] * (*mpModel)( cEx, *pbNode );
				        bcFy += N[i] * (*mpModel)( cEy, *pbNode );
				        bcFz += N[i] * (*mpModel)( cEz, *pbNode );                    
                    }
                }
			}
			else if ( mPBC_NodeElementPairs_Front.find( NodeId ) != mPBC_NodeElementPairs_Front.end() )
			{
				std::vector<int>    bNodes  ( mPBC_NodeElementPairs_Front[ NodeId ] );
                std::vector<double> natCoord( mPBC_NodeElementCoord_Front[ NodeId ] );
		        std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i=0; i<bNodes.size(); i++ )
                {
					Node::Pointer pbNode = mpModel->GetNode( bNodes[i] - 1 );

                    if ( mPotentials_On )
                    {
                        bcFx += N[i] * (*mpModel)( cAx, *pbNode );
				        bcFy += N[i] * (*mpModel)( cAy, *pbNode );
				        bcFz += N[i] * (*mpModel)( cAz, *pbNode );

                        if ( ( (*nit)->pDof( cVs ) != NULL ) && ( i < 3 ) )
                        {
                            bcVs += N[i] * (*mpModel)( cVs, *pbNode );
                        }
                    }
                    else
                    {
                        bcFx += N[i] * (*mpModel)( cEx, *pbNode );
				        bcFy += N[i] * (*mpModel)( cEy, *pbNode );
				        bcFz += N[i] * (*mpModel)( cEz, *pbNode );                    
                    }
                }
			}
            else
            {
                 continue;
            }

            if ( mPotentials_On )
            {
			    (*mpModel)( cAx, **nit ) = bcFx;
			    (*mpModel)( cAy, **nit ) = bcFy;
			    (*mpModel)( cAz, **nit ) = bcFz; 

                if ( (*nit)->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, **nit ) = bcVs; 
                }
            }
            else
            {
			    (*mpModel)( cEx, **nit ) = bcFx;
			    (*mpModel)( cEy, **nit ) = bcFy;
			    (*mpModel)( cEz, **nit ) = bcFz;                
            }
        }
	}

	//********************************************************************************************************************************
    //* - Update electric field values in PBC nodes ( Left-Front nodes )
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_LF_PBC( std::vector<int>& FrontLeftNodes )
	{
		if ( FrontLeftNodes.size() == 0 ) 
        { 
            return; 
        }

		std::vector<int>::iterator it_fln;

		for ( it_fln=FrontLeftNodes.begin(); it_fln!=FrontLeftNodes.end(); ++it_fln )
        {
			int NodeId = (*it_fln);

            Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

            std::complex<double> lcFx( 0.00, 0.00 );
			std::complex<double> lcFy( 0.00, 0.00 );
			std::complex<double> lcFz( 0.00, 0.00 );
            std::complex<double> lcVs( 0.00, 0.00 );
                  
			if ( mPBC_NodeNodePairs_Right.find( NodeId ) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNode = mPBC_NodeNodePairs_Right[ NodeId ];

				Node::Pointer plNode = mpModel->GetNode( lNode-1 );

                if ( mPotentials_On )
                {
				    lcFx = (*mpModel)( cAx, *plNode );
				    lcFy = (*mpModel)( cAy, *plNode );
				    lcFz = (*mpModel)( cAz, *plNode ); 

                    if ( plNode->pDof( cVs ) != NULL ) 
                    {
                        lcVs = (*mpModel)( cVs, *plNode ); 
                    }
                }
                else
                {
				    lcFx = (*mpModel)( cEx, *plNode );
				    lcFy = (*mpModel)( cEy, *plNode );
				    lcFz = (*mpModel)( cEz, *plNode );                
                }
            }
			else if ( mPBC_NodeEdgePairs_Right.find( NodeId ) != mPBC_NodeEdgePairs_Right.end() )
			{
				std::vector<int>    lNodes   ( mPBC_NodeEdgePairs_Right[ NodeId ] );
                double              natCoord = mPBC_NodeEdgeCoord_Right[ NodeId ]  ;
		        std::vector<double> N;

                PBC_NaturalBase_N_Line( N, natCoord );

                for ( int i=0; i<lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[ i ]-1 );

                    if ( mPotentials_On )
                    {
                        lcFx += N[i] * (*mpModel)( cAx, *plNode );
				        lcFy += N[i] * (*mpModel)( cAy, *plNode );
				        lcFz += N[i] * (*mpModel)( cAz, *plNode );

                        if ( ( pNode->pDof( cVs ) != NULL ) && ( i < 2 ) )
                        {
                            lcVs += N[i] * (*mpModel)( cVs, *plNode );
                        }
                    }
                    else
                    {
                        lcFx += N[i] * (*mpModel)( cEx, *plNode );
				        lcFy += N[i] * (*mpModel)( cEy, *plNode );
				        lcFz += N[i] * (*mpModel)( cEz, *plNode );                    
                    }
                }
			}
			else if ( mPBC_NodeElementPairs_Right.find( NodeId ) != mPBC_NodeElementPairs_Right.end() )
			{
				std::vector<int>    lNodes  ( mPBC_NodeElementPairs_Right[NodeId] );
		        std::vector<double> natCoord( mPBC_NodeElementCoord_Right[NodeId] );
                std::vector<double> N;

                PBC_NaturalBase_N_Surface( N, natCoord );

                for ( int i=0; i<lNodes.size(); i++ )
                {
					Node::Pointer plNode = mpModel->GetNode( lNodes[i] - 1 );

                    if ( mPotentials_On )
                    {
                        lcFx += N[i] * (*mpModel)( cAx, *plNode );
				        lcFy += N[i] * (*mpModel)( cAy, *plNode );
				        lcFz += N[i] * (*mpModel)( cAz, *plNode );

                        if ( ( pNode->pDof( cVs ) != NULL ) && ( i < 3 ) )
                        {
                            lcVs += N[i] * (*mpModel)( cVs, *plNode );
                        }
                    }
                    else
                    {
                        lcFx += N[i] * (*mpModel)( cEx, *plNode );
				        lcFy += N[i] * (*mpModel)( cEy, *plNode );
				        lcFz += N[i] * (*mpModel)( cEz, *plNode );                    
                    }
                }
     		}
            else
            {
                continue;
            }

			double cos_A = 1.0;
			double sin_A = 0.0;

			if ( mIsRightPBCTilted == true ) 
            {
                PBC_InvRot( NodeId, cos_A, sin_A );
            }

            if ( mPotentials_On )
            {
			    (*mpModel)( cAx, *pNode ) =  cos_A * lcFx + sin_A * lcFy;
			    (*mpModel)( cAy, *pNode ) = -sin_A * lcFx + cos_A * lcFy;
			    (*mpModel)( cAz, *pNode ) =  lcFz; 

                if ( pNode->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, *pNode ) = lcVs; 
                }
            }
            else
            {
			    (*mpModel)( cEx, *pNode ) =  cos_A * lcFx + sin_A * lcFy;
			    (*mpModel)( cEy, *pNode ) = -sin_A * lcFx + cos_A * lcFy;
			    (*mpModel)( cEz, *pNode ) =  lcFz;                
            }
		}
	}

	//********************************************************************************************************************************
    //* - Updating contact coordinates
    //********************************************************************************************************************************
    void Modeler::UpdateCoord_Contact()
    {
        if ( mContactPairs.size() == 0 ) 
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for ( nit=nodes.begin(); nit!= nodes.end(); ++nit )
        {
			if ( mContactPairs.find( (*nit)->Id() ) != mContactPairs.end() )
            {
				Node::Pointer pLNode = mpModel->GetNode( mContactPairs[(*nit)->Id()].ctcNode - 1 );

                if ( mPotentials_On )
                {
                    (*mpModel)( cAx, **nit ) = (*mpModel)( cAx, *pLNode );
                    (*mpModel)( cAy, **nit ) = (*mpModel)( cAy, *pLNode );
				    (*mpModel)( cAz, **nit ) = (*mpModel)( cAz, *pLNode );   

                    if ( pLNode->pDof( cVs ) != NULL ) 
                    {
                        (*mpModel)( cVs, **nit ) = (*mpModel)( cVs, *pLNode );  
                    }
                }
                else
                {
                    Matrix< std::complex<double> > CtCMatrix( 3, 3 );
                    
                    if ( mColdPlasmaMode )
                    {
                        Get_Contact_Matrix_ColdPlasma( CtCMatrix, (*nit)->Id() );
                    }
                    else
                    {
                        Get_Contact_Matrix_FullWave  ( CtCMatrix, (*nit)->Id() );
                    }
                    
                    std::complex<double> ExL = (*mpModel)( cEx, *pLNode );
                    std::complex<double> EyL = (*mpModel)( cEy, *pLNode );
                    std::complex<double> EzL = (*mpModel)( cEz, *pLNode );
                    
                    (*mpModel)( cEx, **nit ) = CtCMatrix[0][0]*ExL + CtCMatrix[0][1]*EyL + CtCMatrix[0][2]*EzL;
                    (*mpModel)( cEy, **nit ) = CtCMatrix[1][0]*ExL + CtCMatrix[1][1]*EyL + CtCMatrix[1][2]*EzL;
				    (*mpModel)( cEz, **nit ) = CtCMatrix[2][0]*ExL + CtCMatrix[2][1]*EyL + CtCMatrix[2][2]*EzL;                
                }
            }
        }
    }

	//********************************************************************************************************************************
    //* - Returns the area weighted normal of a set of normals.
    //* - The module of the normals must be proportional to the element area.
    //********************************************************************************************************************************
    Vector<double> Modeler::Calculate_Area_Weighted_Normal( std::vector< std::vector<double> >& vectorSet )
    {
        // Area weighted normal
        Vector<double> aw_normal( 3, 0.00 );

        // Iterator over the set of normals
        std::vector< std::vector<double> >::iterator it_vs;

        // Adding all the normals in a node
        for( it_vs = vectorSet.begin(); it_vs != vectorSet.end(); ++it_vs )
        {
            Vector<double> nInSet( *it_vs );

            aw_normal += nInSet;
        }

        // Normalizing the total normal
        aw_normal *= ( 1.00 / sqrt( std::inner_product( aw_normal.begin(), aw_normal.end(), aw_normal.begin(), 0.00 ) ) );

        return aw_normal;
    }

    //********************************************************************************************************************************
    //* - Returns the geometric averaged normal of a set of normals.
    //********************************************************************************************************************************
    Vector<double> Modeler::Calculate_Geom_Averaged_Normal( std::vector< std::vector<double> >& vectorSet )
    {
        // Geometric average of a set of vectors
        Vector<double> ga_normal( 3, 0.00 );

        // Number of vectors in the set
        int setSize = vectorSet.size();

        // Group of vectors with a separation angle between them lower than 45 degrees
        Vector< Vector<double> > vector_group;

        // True is the vector has been assigned to a group
        Vector<bool> used_vector( setSize, 0 );

        // Vectors to compare
        Vector<double> v1, v2;

        // Dot product of v1 and v2
        double v1_dot_v2;

        // Iteration over vector set
        for ( int ivs = 0; ivs < setSize; ivs++ )
        {
            // If the vector is assigned -> next vector
            if ( used_vector[ivs] == true ) continue;

            // Take vector from set and normalize
            v1  = vectorSet[ivs];
            v1 *= ( 1.00 / sqrt( std::inner_product( v1.begin(), v1.end(), v1.begin(), 0.00 ) ) );

            // Initiate group of vectors
            Vector<double> v_group(v1);

            // Look for the vectors with a separation angle with v1 lower than 45 degrees
            for ( int jvs = ivs+1; jvs < setSize; jvs++ )
            {
                // If the vector has been assigned -> next vector
                if ( used_vector[jvs] == true ) continue;

                // Take vector from set and normalize
                v2  = vectorSet[jvs];
                v2 *= ( 1.00 / sqrt( std::inner_product( v2.begin(), v2.end(), v2.begin(), 0.00 ) ) );

                // Dot product v1*v2
                v1_dot_v2 = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

                // If v1 and v2 angle separation is lower than 45 -> add to n1 group
                if ( v1_dot_v2 > 0.5 )
                {
                    v_group += v2;

                    used_vector[jvs] = true;
                }
            }

            // Normalize the group resultant
            v_group *= ( 1.00 / sqrt( std::inner_product( v_group.begin(), v_group.end(), v_group.begin(), 0.00 ) ) );

            // Add group resultant to the geometric average normal
            ga_normal += v_group;
        }

        // Normalize geometric average normal
        ga_normal *= ( 1.00 / sqrt( std::inner_product( ga_normal.begin(), ga_normal.end(), ga_normal.begin(), 0.00 ) ) );

        return ga_normal;
    }

	//********************************************************************************************************************************
    //* - Setting normals in a PEC node ( nxE = 0 )
    //********************************************************************************************************************************
	void Modeler::SettingPECNormals()
    {
		std::map< unsigned int, std::vector< std::vector<double> > >::iterator it_Set;

        for( it_Set = mSetPECnormals.begin(); it_Set != mSetPECnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd ) 
            {
                normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            }
            else    
            {
                normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );
            }

			// Including the calculated normal in the mNormal vector
            mNormals[ it_Set->first ] = normal_atNode;

            // Setting degrees of freedom for a PEC node ( nxE = 0 )
			Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

            // Setting DOFs in a PEC node
            std::complex<double> cZero( 0.00, 0.00 );

            if ( mPotentials_On )
            {
                (*mpModel)( cAy, *pNode ) = cZero;
			    (*mpModel)( cAz, *pNode ) = cZero;

                pNode->pDofcAy()->FixDof();
                pNode->pDofcAz()->FixDof();  

                if ( pNode->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, *pNode ) = cZero;

                    pNode->pDofcVs()->FixDof(); 
                }
            }
            else
            {
                (*mpModel)( cEy, *pNode ) = cZero;
			    (*mpModel)( cEz, *pNode ) = cZero;

                pNode->pDofcEy()->FixDof();
                pNode->pDofcEz()->FixDof();            
            }
            
            mType_Of_BC_Normal[ pNode->Id() ] = 'E';
        }

        // Cleaning mSetPECnormals 
		mSetPECnormals.clear();
		std::map<unsigned int, std::vector<std::vector<double> > >().swap( mSetPECnormals );
	}

    //**********************************************************************************************************************
    //* - Setting normals in a PSBC node ( n*D = 0 )
    //**********************************************************************************************************************
	void Modeler::SettingPSBCNormals()
    {
		std::map< unsigned int, std::vector< std::vector<double> > >::iterator it_Set;

        for ( it_Set = mSetPSBCnormals.begin(); it_Set != mSetPSBCnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            if ( mNormalsAreaWtd == true ) normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );
            else                           normal_atNode = Calculate_Geom_Averaged_Normal( vectorSet );

            // If the node is a PEC -> calculate the average normal
			if ( mNormals.find( it_Set->first ) != mNormals.end() )
			{
				Vector<double> vPEC ( mNormals[it_Set->first] );
				Vector<double> vPSBC( normal_atNode           );

				vPSBC += vPEC;
                vPSBC *= ( 1.00 / sqrt( std::inner_product( vPSBC.begin(), vPSBC.end(), vPSBC.begin(), 0.00 ) ) );

				mNormals[it_Set->first] = vPSBC;
			}
            // if the node is not PEC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[it_Set->first] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

                std::complex<double> cZero( 0.00, 0.00 );

                if ( mPotentials_On ) 
                { 
                    (*mpModel)( cAx, *pNode ) = cZero; 
                    
                    pNode->pDofcAx()->FixDof(); 
                }
                else                  
                { 
                    (*mpModel)( cEx, *pNode ) = cZero; 
                    
                    pNode->pDofcEx()->FixDof(); 
                }
			        
                mType_Of_BC_Normal[pNode->Id()] = 'D';
			}
        }

        // Cleaning mSetPSBCnormals 
		mSetPSBCnormals.clear();
		std::map<unsigned int, std::vector<std::vector<double> > >().swap( mSetPSBCnormals );
	}

	//**********************************************************************************************************************
    //* - Setting normals in a PMC node ( n*E = 0 )
    //**********************************************************************************************************************
	void Modeler::SettingPMCNormals()
    {
		std::map<unsigned int, std::vector<std::vector<double> > >::iterator it_Set;

		for ( it_Set = mSetPMCnormals.begin(); it_Set != mSetPMCnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC -> make the normals ortogonal
			if ( mNormals.find( it_Set->first ) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[ it_Set->first ] );

				Vector<double> vOrto( 3, 0.00 );

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[ it_Set->first ] = vOrto;
			}
            // if the node is not PEC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[ it_Set->first ] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

                std::complex<double> cZero( 0.00, 0.00 );

                if ( mPotentials_On ) 
                { 
                    (*mpModel)( cAx, *pNode ) = cZero; 
                    
                    pNode->pDofcAx()->FixDof(); 
                }
                else                  
                { 
                    (*mpModel)( cEx, *pNode ) = cZero; 
                    
                    pNode->pDofcEx()->FixDof(); 
                }

                mType_Of_BC_Normal[ pNode->Id() ] = 'H';
			}
        }

        // Cleaning mSetPMCnormals 
        mSetPMCnormals.clear();
        std::map<unsigned int, std::vector<std::vector<double> > >().swap( mSetPMCnormals );
	}

	//**********************************************************************************************************************
    //* - Setting normals in a TE port ( ortogonal to a PMC surface )
    //**********************************************************************************************************************
	void Modeler::SettingTEPMCNormals()
    {
		std::map<unsigned int, std::vector<std::vector<double> > >::iterator it_Set;

		for ( it_Set = mSetTEPMCnormals.begin(); it_Set != mSetTEPMCnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            std::vector< std::vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // TE_PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC or PMC -> make the normals ortogonal
			if ( mNormals.find( it_Set->first ) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[ it_Set->first ] );

				Vector<double> vOrto( 3, 0.00 );

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[ it_Set->first ] = vOrto;
			}
            // If the node is not PEC/PMC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[ it_Set->first ] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

                std::complex<double> cZero( 0.00, 0.00 );

                if ( mPotentials_On ) 
                { 
                    (*mpModel)( cAx, *pNode ) = cZero; 
                    
                    pNode->pDofcAx()->FixDof(); 
                }
                else                  
                { 
                    (*mpModel)( cEx, *pNode ) = cZero; 
                    
                    pNode->pDofcEx()->FixDof(); 
                }

				mType_Of_BC_Normal[ pNode->Id() ] = 'T';
			}
        }

        // Cleaning mSetTEPMCnormals 
        mSetTEPMCnormals.clear();
        std::map<unsigned int, std::vector<std::vector<double> > >().swap( mSetTEPMCnormals );
	}

	//**********************************************************************************************************************
    //* - Making the vector vToBeOrto ortogonal to the vector vAxis
    //**********************************************************************************************************************
	void Modeler::Ortogonalization( Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto )
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

	//**********************************************************************************************************************
    //* - Ortogonalization of the contact normals respect to PEC and PMC normals
    //**********************************************************************************************************************
	void Modeler::OrtogonalizeContactNormals()
    {
		std::map<unsigned int, ContactPairData>::iterator itCtc;

        for( itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc )
        {
			unsigned int fCtcNode  = itCtc->first;

			if ( ( mNormals.find( fCtcNode ) != mNormals.end() ) )
            {
				Vector<double> sNormal  ( mNormals     [ fCtcNode ]           );
				Vector<double> vToBeOrto( mContactPairs[ fCtcNode ].ctcNormal );
				Vector<double> vOrto    ( 3, 0.00 );

				Ortogonalization( vToBeOrto, sNormal, vOrto );

				mContactPairs[ fCtcNode ].ctcNormal = vOrto;
			}
		}
	}

    //**********************************************************************************************************************
    //* - Setting normals in PEC and PMC surfaces
    //**********************************************************************************************************************
    void Modeler::SetNormals()
    {
        // We join to the Dirichlet sets the contacts nodes normals
		JoinContactNormalsInDirichletS();

		// Setting Dirichlet normals
		SettingPECNormals ();
        SettingPSBCNormals();

		// Dirichlet normals are forced in dielectric-Dirichlet intersections.
		ForceDirichletNormalsInContact();

		// Setting simmetry normals
		SettingPMCNormals  ();
		SettingTEPMCNormals();

		// Contact normals must be ortogonal to simmetry normals
		OrtogonalizeContactNormals();
    }

    //*******************************************************************************
    //* - Initialize Building process
    //*******************************************************************************
    void Modeler::Ini_Building()
    {
        std::cout << "Nodes read."    << std::endl << std::endl;

        std::cout << "Setting normals and DOFs..." << std::endl;

		SetContacts();
		SetNormals ();
		SetPBC     ();

		if ( mSetAllEzToZero  ) SetAllcEzToCero ();
		if ( mSetAllExyToZero ) SetAllcExyToCero();
		if ( mAxisymmetric    ) SetAllcEyzToCero();

        SetDofSet();

		std::cout << "Done." << std::endl << std::endl;
    }

	//*******************************************************************************
    //* - Set all cEz to cero ( 3D-Exy problems )
    //*******************************************************************************
    void Modeler::SetAllcEzToCero()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		
        std::vector<Node::Pointer>::iterator it;
		
        std::complex<double> cZero( 0.00, 0.00 );

		for ( it = nodes.begin(); it != nodes.end(); ++it )
		{
			Node::Pointer pNode = (*it);

			(*mpModel)( cEz, *pNode ) = cZero;

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

    //**********************************************************************************************************************
    //* - Solve linear system
    //**********************************************************************************************************************
    void Modeler::SolveLinearSystem()
    {
        if ( mReadSolutionMode )
		{
			ReadSolutionVector( x_vector );
			A_matrix.FreeData();
			b_vector.FreeData();
            A_matrix_aux.FreeData();
            return;
		}

		if ( mSolveWithExternal )
		{
			WriteMatrixInFile( A_matrix );
			WriteVectorInFile( b_vector );

			A_matrix.FreeData();
			b_vector.FreeData();

            if ( mA_matrix_aux_Required == true )
            {
                WriteAuxMatrixInFile( A_matrix_aux );
            }

            A_matrix_aux.FreeData();

			if ( mExternalSolverPath == "" )
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

            system( MSDOSOrder );

			std::cout << std::endl << "External solver has finished." << std::endl << std::endl;

			ReadSolutionVector( x_vector );
		}
		else
		{
			std::cout << "Solving linear system..." << std::endl << std::endl;

			if ( mReadFileInitialGuess ) ReadInitialGuess( x_vector );

			mpLinearSolver->SetWriteStepResult( mWriteResultEveryStep );

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

            if ( mWriteResultFinalStep ) WriteResultsInFile( x_vector );

            std::cout << "Solver finished." << std::endl << std::endl;
		}

        A_matrix    .FreeData();
        b_vector    .FreeData();
        A_matrix_aux.FreeData();
    }

    //*************************************************
    //* - Set singular nodes
    //*************************************************
    void Modeler::FixDof( unsigned int NodeId, const Variable<double>& rVariable, const double& Value )
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
    void Modeler::FixCDof( unsigned int NodeId, const Variable<double>& rVariable, const std::vector<double>& vValue )
    {
    }

    //*******************************************************************************
    //* - Set Linear solver in plasma mode
    //*******************************************************************************
    void Modeler::GenerateLinearSolver_ColdPlasma( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
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

    //*********************************************************************************************
    //* - Set Dofs
    //*********************************************************************************************
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

        if ( mColdPlasmaMode == true )
        {
            if ( !mpColdPlasma->Is_Full_Matrix() )
            {
                mA_matrix_aux_Required = true;
                A_matrix_aux.Resize(mSystemSize);
            }
        }
    }

	//*********************************************************************************************
    //* - Print projections
    //*********************************************************************************************
	void Modeler::PrintProjection()
    {
		if ( mProjection.size() > 0 )
		{
		    CalculateSijParameters();
		    PrintSijParameters( mProjectionNormalized, mProblemFrequency );
		}

		if ( mVolumeIntegralE.size() > 0 )
        {
            PrintViParameters( mVolumeIntegralE, mVolumeIntegralH, mVolumeIntegralJ,
				               mVolumeValue    , mProblemFrequency );
        }

		if ( mSurfaceIntegralE.size() > 0 )
        {
	        PrintSurfiParameters( mSurfaceIntegralE, mSurfaceIntegralH, mSurfaceIntegralJ,
				                  mSurfaceValue    , mProblemFrequency );
		}
    }

    //*********************************************************************************************
    //* - Print the results of the field integrals in electrostatic mode.
    //*********************************************************************************************
    void Modeler::PrintProjection_Electrostatic()
    {
        if ( mVolumeValue.size() > 0 )
        {
            PrintViParameters_Electrostatic( mVolIntegralE_Electrostatic, mVolumeValue );
        }

        if ( mSurfaceValue.size() > 0 )
        {
            PrintSurfiParameters_Electrostatic( mSurfIntegralE_Electrostatic, mSurfaceValue );
        }
    }

    //*********************************************************************************************
    //* - Setting problem type
    //*********************************************************************************************
    void Modeler::SetProblemType( const String& ProblemType )
    {
		if ( ProblemType == "Full_wave" )
		{
			mProblemType       = E3D;
			mSetAllEzToZero    = false;
			mSetAllExyToZero   = false;
			mAxisymmetric      = false;
			mColdPlasmaMode    = false;
            mElectrostaticMode = false;
            mFullWaveMode      = true ;
		}
        else if ( ProblemType == "E3D" )
		{
			mSetAllEzToZero    = false;
			mSetAllExyToZero   = false;
			mAxisymmetric      = false;
		}
		else if ( ProblemType == "Ez3D" )
		{
			mSetAllEzToZero    = false;
			mSetAllExyToZero   = true ;
			mAxisymmetric      = false;
		}
		else if ( ProblemType == "Exy3D" )
		{
			mSetAllEzToZero    = true ;
			mSetAllExyToZero   = false;
			mAxisymmetric      = false;
		}
		else if ( ProblemType == "Ea3D" )
		{
			mSetAllEzToZero    = false;
			mSetAllExyToZero   = false;
			mAxisymmetric      = true ;
		}
		else if ( ProblemType == "Cold_plasma" )
		{
			mProblemType       = E3D;
			mSetAllEzToZero    = false;
			mSetAllExyToZero   = false;
			mAxisymmetric      = false;
			mColdPlasmaMode    = true ;
            mElectrostaticMode = false;
            mFullWaveMode      = false;
		}
        else if ( ProblemType == "Electrostatic" )
        {
            mProblemType       = E3D;
            mSetAllEzToZero    = false;
            mSetAllExyToZero   = false;
            mAxisymmetric      = false;
            mColdPlasmaMode    = false;
            mElectrostaticMode = true ;
            mFullWaveMode      = false;
        }
		else if ( ProblemType == "RELSSOL" )
		{
			mReleaseSolutionMode = true ;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = false;
		}
		else if ( ProblemType == "DEBGSOL" )
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = true ;
			mReadSolutionMode    = false;
		}
		else if ( ProblemType == "READSOL" )
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = true ;
		}
		else if ( ProblemType == "Results_On_Nodes" )
		{
			mSmoothing_On = true;
            mResultsOnGPs = 0;
		}
		else if ( ProblemType == "Results_On_1GP" )
		{
			mSmoothing_On = false;
            mResultsOnGPs = 1;
		}
        else if ( ProblemType == "Results_On_4GP" )
		{
			mSmoothing_On = false;
            mResultsOnGPs = 4;
		}
        else if ( ProblemType == "Results_On_10GP" )
		{
			mSmoothing_On = false;
            mResultsOnGPs = 10;
		}
        else if ( ProblemType == "Potentials_On" )
		{
			mPotentials_On = true;
		}
		else if ( ProblemType == "Potentials_Off" )
		{
			mPotentials_On = false;
		}
        else if ( ProblemType == "AV_continuity_On" )
		{
			mAVContinuity_On = true;
		}
		else if ( ProblemType == "AV_continuity_Off" )
		{
			mAVContinuity_On = false;
		}
        else if ( ProblemType == "Ini_Surf_Counter" )
		{
			mNonSmooth_Surf_Counter = 0;
		}
	    else if ( ProblemType == "Project" )
		{
			PrintProjection();
		}
        else if ( ProblemType == "Project_Electrostatic_Fields" )
        {
            PrintProjection_Electrostatic();
        }
		else if ( ProblemType == "GAv" )
		{
			mNormalsGeomAvg = true ;
			mNormalsAreaWtd = false;
		}
		else if ( ProblemType == "AWg" )
		{
			mNormalsGeomAvg = false;
			mNormalsAreaWtd = true ;
		}
		else if ( ProblemType == "RFIG" )
		{
			mReadFileInitialGuess = true;
		}
		else if ( ProblemType == "NRFIG" )
		{
			mReadFileInitialGuess = false;
		}
		else if ( ProblemType == "NWRIF" )
		{
			mWriteResultEveryStep = false;
			mWriteResultFinalStep = false;
		}
		else if ( ProblemType == "ESWRIF" )
		{
			mWriteResultEveryStep = true;
			mWriteResultFinalStep = true;
		}
		else if ( ProblemType == "FSWRIF" )
		{
			mWriteResultEveryStep = false;
			mWriteResultFinalStep = true ;
		}
		else if ( ProblemType == "OFFBIN" )
		{
			mOutputFileFormatIsBIN   = true ;
			mOutputFileFormatIsASCII = false;
		}
		else if ( ProblemType == "OFFASCII" )
		{
			mOutputFileFormatIsBIN   = false;
			mOutputFileFormatIsASCII = true ;
		}
        else if ( ProblemType == "IMPJON" )
		{
			mImportCurrents = true;
		}
        else if ( ProblemType == "IMPJOFF" )
		{
			mImportCurrents = false;
		}
		else if ( ProblemType == "QE" )
		{
			mQuadraticGeometry = true;
		}
		else if ( ProblemType == "LE" )
		{
			mQuadraticGeometry = false;
		}
        else if ( ProblemType == "CheckConsistency" )
        {
            CheckConsistency();
        }

        else if ( ProblemType == "LL2P_3sb" ) { mElementOrder = 0; }
        else if ( ProblemType ==  "RME_1st" ) { mElementOrder = 1; }
        else if ( ProblemType ==  "RME_2nd" ) { mElementOrder = 2; }

        else if ( ProblemType == "GiDTol3"  ) { mGiDTolerance = 1e-3 ; }
        else if ( ProblemType == "GiDTol6"  ) { mGiDTolerance = 1e-6 ; }
        else if ( ProblemType == "GiDTol9"  ) { mGiDTolerance = 1e-9 ; }
        else if ( ProblemType == "GiDTol12" ) { mGiDTolerance = 1e-12; }

		else if ( ProblemType ==   "1pr" ) { mItSolverNumThreads =   1; }
        else if ( ProblemType ==   "2pr" ) { mItSolverNumThreads =   2; }
        else if ( ProblemType ==   "4pr" ) { mItSolverNumThreads =   4; }
        else if ( ProblemType ==   "8pr" ) { mItSolverNumThreads =   8; }
		else if ( ProblemType ==  "16pr" ) { mItSolverNumThreads =  16; }
        else if ( ProblemType ==  "32pr" ) { mItSolverNumThreads =  32; }
        else if ( ProblemType ==  "64pr" ) { mItSolverNumThreads =  64; }
        else if ( ProblemType == "128pr" ) { mItSolverNumThreads = 128; }

        else if ( ProblemType == "LLP2_hk_0"  ) { mLL2P_hk_factor = 0.0;  }
        else if ( ProblemType == "LLP2_hk_1"  ) { mLL2P_hk_factor = 1.0;  }
        else if ( ProblemType == "LLP2_hk_e1" ) { mLL2P_hk_factor = 1e-1; }
        else if ( ProblemType == "LLP2_hk_e2" ) { mLL2P_hk_factor = 1e-2; }
        else if ( ProblemType == "LLP2_hk_e3" ) { mLL2P_hk_factor = 1e-3; }
        else if ( ProblemType == "LLP2_hk_e4" ) { mLL2P_hk_factor = 1e-4; }
		else if ( ProblemType == "LLP2_hk_e5" ) { mLL2P_hk_factor = 1e-5; }
        else if ( ProblemType == "LLP2_hk_e6" ) { mLL2P_hk_factor = 1e-6; }
        else if ( ProblemType == "LLP2_hk_e7" ) { mLL2P_hk_factor = 1e-7; }
        else if ( ProblemType == "LLP2_hk_e8" ) { mLL2P_hk_factor = 1e-8; }
        else if ( ProblemType == "LLP2_hk_e9" ) { mLL2P_hk_factor = 1e-9; }
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
    //* - Returns true if the problem mode is Electrostatic.
    //*******************************************************
    bool Modeler::Is_ElectrostaticMode()
    {
        return mElectrostaticMode;
    }

    //*******************************************************
    //* - Returns true if the problem mode is Cold_plasma.
    //*******************************************************
    bool Modeler::Is_ColdPlasmaMode()
    {
        return mColdPlasmaMode;
    }

    //*******************************************************
    //* - Returns true if the problem mode is Full_wave.
    //*******************************************************
    bool Modeler::Is_FullWaveMode()
    {
        return mFullWaveMode;
    }

    //*******************************************************
    //* - Returns true if we are using a quadratic mesh.
    //*******************************************************
    bool Modeler::Is_QuadraticMesh()
    {
        return mQuadraticGeometry;
    }

    //*******************************************************
    //* - Returns true when smoothing is activated 
    //*******************************************************
    bool Modeler::Is_SmoothingOn()
    {
        return mSmoothing_On;
    }

    //*******************************************************
    //* - Returns true when potentials are in use
    //*******************************************************
    bool Modeler::Is_PotentialsOn()
    {
        return mPotentials_On;
    }

    //*******************************************************
    //* - Returns element order.
    //*******************************************************
    int Modeler::Get_Element_Order()
    {
        return mElementOrder;
    }

    //*******************************************************
    //* - Returns problem frequency.
    //*******************************************************
    double Modeler::Get_Problem_Frequency()
    {
        return mProblemFrequency;
    }

    //*******************************************************
    //* - Setting frequency mode
    //*******************************************************
    void Modeler::SetFrequencyMode( bool SweepingFreq )
    {
        mFrequencySweep = SweepingFreq;
    }

    //*******************************************************
    //* - Setting frequency
    //*******************************************************
    void Modeler::SetFrequency( double ProblemFrequency )
    {
        mProblemFrequency = ProblemFrequency;
    }

    //******************************************************************************
    //* - Gives boundary normal at NodeId
    //******************************************************************************
    void Modeler::Get_Boundary_Normal( Vector<double>& BoundaryNormal, int NodeId )
    {
        BoundaryNormal.resize( 3 );

        if ( mNormals.find( NodeId ) != mNormals.end() )
        {
            BoundaryNormal[0] = mNormals[ NodeId ][ 0 ];
            BoundaryNormal[1] = mNormals[ NodeId ][ 1 ];
            BoundaryNormal[2] = mNormals[ NodeId ][ 2 ];
        }
        else
        {
            BoundaryNormal[0] = 0.00;
            BoundaryNormal[1] = 0.00;
            BoundaryNormal[2] = 0.00;        
        }
    }

    //******************************************************************************
    //* - Gives contact pair normal at NodeId
    //******************************************************************************
    void Modeler::Get_ContactP_Normal( Vector<double>& ContactPNormal, int NodeId )
    {
        ContactPNormal.resize( 3 );

        if ( mContactPairs.find( NodeId ) != mContactPairs.end() )
        {
            ContactPNormal[0] = mContactPairs[ NodeId ].ctcNormal[ 0 ];
            ContactPNormal[1] = mContactPairs[ NodeId ].ctcNormal[ 1 ];
            ContactPNormal[2] = mContactPairs[ NodeId ].ctcNormal[ 2 ];
        }
        else
        {
            ContactPNormal[0] = 0.00;
            ContactPNormal[1] = 0.00;
            ContactPNormal[2] = 0.00;        
        }
    }

	//*************************************************************************************************************************
	//* - Creates plasma object and reads plasma data from files
	//*************************************************************************************************************************
	void Modeler::LoadPlasmaParameters( String plasma_file_name )
	{
		mpColdPlasma = ColdPlasma::Pointer(new ColdPlasma());
		mpColdPlasma->Load_Data(plasma_file_name);
	}

	//*************************************************************************************************************************
	//* - Obtains electron density at every node
	//*************************************************************************************************************************
	void Modeler::Extract_ColdPlasma_ElectronDensity()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			(*mpModel)(ELECTRON_DENSITY, **it) = mpColdPlasma->Get_ElectronDensity_InNode( *it );
		}
	}

	//*************************************************************************************************************************
	//* - Obtains applied external B field at every node
	//*************************************************************************************************************************
	void Modeler::Extract_ColdPlasma_Bexternal()
	{
		Model::NodesArrayType::iterator it;

		for (it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it)
		{
			(*mpModel)(B_EXT, **it) = mpColdPlasma->Get_Bexternal_InNode( *it );
		}
	}

	//*************************************************************************************************************************
	//* - Extracts S,D,P,R,L components of the plasma permittivity tensor at every node
	//*************************************************************************************************************************
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

    //*************************************************************************************************************************
    //* - Obtains E field component parallel to B_ext
    //*************************************************************************************************************************
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

    //*************************************************************************************************************************
    //* - Obtains E field component perpendicular to B_ext
    //*************************************************************************************************************************
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

	//*************************************************************************************************************************
	//* - Extracts all the plasma parameters at every node
	//*************************************************************************************************************************
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

    //*************************************************************************************************************************
    //* - Set to zero de parallel component of E when tolerance criteria is met
    //*************************************************************************************************************************
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

    //*************************************************************************************************************************
    //* - Applies PEC and PMC conditions to H
    //*************************************************************************************************************************
    void Modeler::Apply_PEC_and_PMC_To_H()
    {
        if ( mNormals.size() == 0 ) 
        {
            return;
        }

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        for ( it = nodes.begin(); it != nodes.end(); ++it )
        {
            int NodeId = (*it)->Id();

            if ( mNormals.find( NodeId ) != mNormals.end() )
            {
                Vector<double> n( mNormals[ NodeId ] );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

                TangencialCoordinates( n, t, b );

                Vector<double> free_H_real = (*mpModel)( REAL_H, **it );
                Vector<double> free_H_imag = (*mpModel)( IMAG_H, **it );

                Vector<double> fixd_H_real( 3, 0.00 );
                Vector<double> fixd_H_imag( 3, 0.00 );

                // If PEC make n*H = 0
                if ( mType_Of_BC_Normal[ NodeId ] == 'E' )
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
                else if ( mType_Of_BC_Normal[ NodeId ] == 'H' )
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

                (*mpModel)( REAL_H, **it ) = fixd_H_real;
                (*mpModel)( IMAG_H, **it ) = fixd_H_imag;
            }
        }
    }

    //*************************************************************************************************************************
    //* - Applies PEC and PMC conditions to E
    //*************************************************************************************************************************
    void Modeler::Apply_PEC_and_PMC_To_E()
    {
        if ( mNormals.size() == 0 ) 
        {
            return;
        }

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        for ( it = nodes.begin(); it != nodes.end(); ++it )
        {
            int NodeId = (*it)->Id();

            if ( mNormals.find( NodeId ) != mNormals.end() )
            {
                Vector<double> n( mNormals[ NodeId ] );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

                TangencialCoordinates( n, t, b );

                Vector<double> free_E_real = (*mpModel)( REAL_E, **it );
                Vector<double> free_E_imag = (*mpModel)( IMAG_E, **it );

                Vector<double> fixd_E_real( 3, 0.00 );
                Vector<double> fixd_E_imag( 3, 0.00 );

                // If PEC make nxE = 0
                if ( mType_Of_BC_Normal[ NodeId ] == 'E' )                
                {
                    double En_real = n[0] * free_E_real[0] + n[1] * free_E_real[1] + n[2] * free_E_real[2];
                    double En_imag = n[0] * free_E_imag[0] + n[1] * free_E_imag[1] + n[2] * free_E_imag[2];

                    fixd_E_real[0] = n[0] * En_real;
                    fixd_E_real[1] = n[1] * En_real;
                    fixd_E_real[2] = n[2] * En_real;

                    fixd_E_imag[0] = n[0] * En_imag;
                    fixd_E_imag[1] = n[1] * En_imag;
                    fixd_E_imag[2] = n[2] * En_imag;
                }
                // If PMC make n*E = 0
                else if ( mType_Of_BC_Normal[ NodeId ] == 'H' )
                {
                    double Et_real = t[0] * free_E_real[0] + t[1] * free_E_real[1] + t[2] * free_E_real[2];
                    double Eb_real = b[0] * free_E_real[0] + b[1] * free_E_real[1] + b[2] * free_E_real[2];

                    double Et_imag = t[0] * free_E_imag[0] + t[1] * free_E_imag[1] + t[2] * free_E_imag[2];
                    double Eb_imag = b[0] * free_E_imag[0] + b[1] * free_E_imag[1] + b[2] * free_E_imag[2];

                    fixd_E_real[0] = t[0] * Et_real + b[0] * Eb_real;
                    fixd_E_real[1] = t[1] * Et_real + b[1] * Eb_real;
                    fixd_E_real[2] = t[2] * Et_real + b[2] * Eb_real;

                    fixd_E_imag[0] = t[0] * Et_imag + b[0] * Eb_imag;
                    fixd_E_imag[1] = t[1] * Et_imag + b[1] * Eb_imag;
                    fixd_E_imag[2] = t[2] * Et_imag + b[2] * Eb_imag;
                }
                // If not PEC or PMC let E as calculated
                else
                {
                    fixd_E_real[0] = free_E_real[0];
                    fixd_E_real[1] = free_E_real[1];
                    fixd_E_real[2] = free_E_real[2];

                    fixd_E_imag[0] = free_E_imag[0];
                    fixd_E_imag[1] = free_E_imag[1];
                    fixd_E_imag[2] = free_E_imag[2];                
                }

                (*mpModel)( REAL_E, **it ) = fixd_E_real;
                (*mpModel)( IMAG_E, **it ) = fixd_E_imag;
            }
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to calculate electrostatic E field.
    //**********************************************************************************************************************
    void Modeler::Ini_Electrostatic_E_Derivation()
    {
        if ( mSmoothing_On == false ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) mNeighbourElements[i] = 0;

        Vector<double> vZero( 3, 0.00 );

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it ) = vZero;
        }
    }

    //**********************************************************************************************************************
    //* - Ending calculation of electrostatic E field.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_E_Derivation()
    {
        if ( mSmoothing_On == false ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        Vector<double> Electrostatic_E;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_E = (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it );

            if ( NumNeighbours > 0 ) Electrostatic_E /= NumNeighbours;
  
            (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it ) = Electrostatic_E;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to calculate electrostatic J current density.
    //**********************************************************************************************************************
    void Modeler::Ini_Electrostatic_J_Derivation()
    {
        if ( mSmoothing_On == false ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) mNeighbourElements[i] = 0;

        Vector<double> vZero( 3, 0.00 );

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it ) = vZero;
        }
    }

    //**********************************************************************************************************************
    //* - Ending calculation of electrostatic J current density.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_J_Derivation()
    {
        if ( mSmoothing_On == false ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        Vector<double> Electrostatic_J;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_J = (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it );

            if ( NumNeighbours > 0 ) Electrostatic_J /= NumNeighbours;

            (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it ) = Electrostatic_J;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to calculate electrostatic Joule heating.
    //**********************************************************************************************************************
    void Modeler::Ini_Electrostatic_JouleH_Derivation()
    {
        if ( mSmoothing_On == false ) return;

        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) mNeighbourElements[i] = 0;

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it ) = 0.0;
        }
    }

    //**********************************************************************************************************************
    //* - Ending calculation of electrostatic Joule heating.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_JouleH_Derivation()
    {
        if ( mSmoothing_On == false ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        double Electrostatic_JouleH;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_JouleH = (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it );

            if ( NumNeighbours > 0 ) Electrostatic_JouleH /= NumNeighbours;

            (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it ) = Electrostatic_JouleH;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to derive E field
    //**********************************************************************************************************************
	void Modeler::Ini_E_Derivation()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        if ( mPotentials_On )
        {
            Vector<double> vZero( 3, 0.00 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
		        (*mpModel)( REAL_E, **it ) = vZero;
                (*mpModel)( IMAG_E, **it ) = vZero;
            }                    
            
            mNeighbourElements.resize( nodes.size() ); 
        
            for ( int i=0; i<nodes.size(); i++ ) 
            {
                mNeighbourElements[i] = 0;
            }
        }
        else
        {
		    Vector<double> E_real( 3 );
            Vector<double> E_imag( 3 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                E_real[0] = std::real( (*mpModel)( cEx, **it ) );
                E_real[1] = std::real( (*mpModel)( cEy, **it ) );
                E_real[2] = std::real( (*mpModel)( cEz, **it ) );
                                                             
                E_imag[0] = std::imag( (*mpModel)( cEx, **it ) );
                E_imag[1] = std::imag( (*mpModel)( cEy, **it ) );
                E_imag[2] = std::imag( (*mpModel)( cEz, **it ) );

                (*mpModel)( REAL_E, **it ) = E_real;
                (*mpModel)( IMAG_E, **it ) = E_imag;
            }        
        }
	}

    //**********************************************************************************************************************
    //* - End calculate E field
    //**********************************************************************************************************************
	void Modeler::End_E_Derivation()
	{
		if ( !mPotentials_On ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;
        
        Vector<double> E_real;
        Vector<double> E_imag;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            E_real = (*mpModel)( REAL_E, **it );
            E_imag = (*mpModel)( IMAG_E, **it );

            if ( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
            {
                E_real /= mNeighbourElements[ (*it)->Id()-1 ];
                E_imag /= mNeighbourElements[ (*it)->Id()-1 ];
            }

            (*mpModel)( REAL_E, **it ) = E_real;
            (*mpModel)( IMAG_E, **it ) = E_imag;
        }
        
        Apply_PEC_and_PMC_To_E();
         
        std::complex<double> cUnit( 0.00, 1.00 );
        
        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( cEx, **it ) = (*mpModel)( REAL_E, **it )[0] + cUnit * (*mpModel)( IMAG_E, **it )[0];
            (*mpModel)( cEy, **it ) = (*mpModel)( REAL_E, **it )[1] + cUnit * (*mpModel)( IMAG_E, **it )[1];
            (*mpModel)( cEz, **it ) = (*mpModel)( REAL_E, **it )[2] + cUnit * (*mpModel)( IMAG_E, **it )[2];
        }
	}

    //**********************************************************************************************************************
    //* - Calculates E field on elements
    //**********************************************************************************************************************
    void Modeler::Calculate_E_Element( int* NodesId, unsigned int PropertiesId, 
                                       ResultsOnGPsType& REAL_E_OnGP, 
                                       ResultsOnGPsType& IMAG_E_OnGP )
    {
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

        Set_Volume_Element_Parameters( NodesId, pElement );

        if ( mPotentials_On ) 
        {
            cVector2Type cE_OnNodes;
            
            pElement->Calculate_E_field_OnNodes( cE_OnNodes );
            
            Vector2Type rE_OnNodes, iE_OnNodes;
            
            Split_Complex_VectorVector( rE_OnNodes, iE_OnNodes, cE_OnNodes );
            
		    for ( int n=0; n<nodes.size(); n++ )
            {
                (*mpModel)( REAL_E, *nodes[n] ) += rE_OnNodes[n];
                (*mpModel)( IMAG_E, *nodes[n] ) += iE_OnNodes[n];
            
                mNeighbourElements[ nodes[n]->Id() - 1 ]++;
            }
        }

        if ( !mSmoothing_On ) 
        {
            cVector2Type cE_OnGP;
            
            pElement->Calculate_E_field_OnGaussPoints( cE_OnGP, mResultsOnGPs );

            Vector2Type rE_OnGP, iE_OnGP;
            
            Split_Complex_VectorVector( rE_OnGP, iE_OnGP, cE_OnGP );
            
            REAL_E_OnGP.push_back( rE_OnGP );
            IMAG_E_OnGP.push_back( iE_OnGP );
        }
    }

	//**********************************************************************************************************************
    //* - Gets the real and imaginary parts of a Vector< Vector<complex> > type
    //**********************************************************************************************************************
    void Modeler::Split_Complex_VectorVector( Vector2Type& real_VV, Vector2Type& imag_VV, cVector2Type& complex_VV )
    {
        Vector<double> real_V( 3 );
        Vector<double> imag_V( 3 );
 
        cVector2Type::iterator it;

        for ( it = complex_VV.begin(); it != complex_VV.end(); it++ )
        {
            real_V[0] = std::real( (*it)[0] );  
            real_V[1] = std::real( (*it)[1] );  
            real_V[2] = std::real( (*it)[2] ); 

            imag_V[0] = std::imag( (*it)[0] ); 
            imag_V[1] = std::imag( (*it)[1] ); 
            imag_V[2] = std::imag( (*it)[2] );

            real_VV.push_back( real_V );
            imag_VV.push_back( imag_V );
        }
    }

	//**********************************************************************************************************************
    //* - Initial set up to calculate H field
    //**********************************************************************************************************************
	void Modeler::Ini_H_Derivation()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		std::vector<Node::Pointer>::iterator it;
        
        Vector<double> vZero( 3, 0.00 );

		for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
		    (*mpModel)( REAL_H, **it ) = vZero;
            (*mpModel)( IMAG_H, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
	}

    //**********************************************************************************************************************
    //* - End calculate H field.
    //**********************************************************************************************************************
	void Modeler::End_H_Derivation()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;
        
        Vector<double> H_real;
        Vector<double> H_imag;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            H_real = (*mpModel)( REAL_H, **it );
            H_imag = (*mpModel)( IMAG_H, **it );

            if ( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
            {
                H_real /= mNeighbourElements[ (*it)->Id()-1 ];
                H_imag /= mNeighbourElements[ (*it)->Id()-1 ];
            }

            (*mpModel)( REAL_H, **it ) = H_real;
            (*mpModel)( IMAG_H, **it ) = H_imag;
        }
        
        Apply_PEC_and_PMC_To_H();
	}

    //**********************************************************************************************************************
    //* - Initial set up to calculate B field.
    //**********************************************************************************************************************
	void Modeler::Ini_B_Derivation()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		std::vector<Node::Pointer>::iterator it;

        Vector<double> vZero( 3, 0.00 );

		for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
		    (*mpModel)( REAL_B, **it ) = vZero;
            (*mpModel)( IMAG_B, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
	}

    //**********************************************************************************************************************
    //* - End calculate B field.
    //**********************************************************************************************************************
	void Modeler::End_B_Derivation()
	{
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;
        
        Vector<double> B_real;
        Vector<double> B_imag;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            B_real = (*mpModel)( REAL_B, **it );
            B_imag = (*mpModel)( IMAG_B, **it );

            if ( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
            {
                B_real /= mNeighbourElements[ (*it)->Id()-1 ];
                B_imag /= mNeighbourElements[ (*it)->Id()-1 ];
            }

            (*mpModel)( REAL_B, **it ) = B_real;
            (*mpModel)( IMAG_B, **it ) = B_imag;
        }
	}

    //**********************************************************************************************************************
    //* - Initial set up to calculate J.
    //**********************************************************************************************************************
	void Modeler::Ini_J_Calculation()
	{
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		std::vector<Node::Pointer>::iterator it;

        Vector<double> vZero( 3, 0.00 );

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
	        (*mpModel)( REAL_J, **it ) = vZero;
            (*mpModel)( IMAG_J, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
	}

	//**********************************************************************************************************************
    //* - Calculate H field on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_H_Element( int* NodesId, unsigned int PropertiesId, 
                                       ResultsOnGPsType& REAL_H_OnGP, 
                                       ResultsOnGPsType& IMAG_H_OnGP )
	{
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

        Set_Volume_Element_Parameters( NodesId, pElement );

        // Element material properties
        double mo      = ( 4.00e-7 )  *  3.141592653589793238462643383279; 
        double mu_real = (*properties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		double mu_imag = (*properties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;
        
        std::complex<double> cMu  ( mu_real, mu_imag );
        std::complex<double> cUnit(    0.00,    1.00 );
        std::complex<double> cCteMat;

        if ( mPotentials_On ) cCteMat = 1.00 / ( cMu );
        else                  cCteMat = 1.00 / ( cMu * cUnit * mProblemFrequency );
        
        // Calculate H on nodes
        cVector2Type cH_OnNodes;

        pElement->Calculate_Rotational_OnNodes( cH_OnNodes );

        for ( int n=0; n<cH_OnNodes.size(); n++ ) 
        {
            cH_OnNodes[n] *= cCteMat;
        }

        Vector2Type rH_OnNodes, iH_OnNodes;

        Split_Complex_VectorVector( rH_OnNodes, iH_OnNodes, cH_OnNodes );

		for ( int n=0; n<nodes.size(); n++ )
        {
            (*mpModel)( REAL_H, *nodes[n] ) += rH_OnNodes[n];
            (*mpModel)( IMAG_H, *nodes[n] ) += iH_OnNodes[n];

            mNeighbourElements[ nodes[n]->Id() - 1 ]++;
        }

        // If smoothing is off then calculate H on GPs 
        if ( !mSmoothing_On ) 
        {
            cVector2Type cH_OnGPs;
            
            pElement->Calculate_Rotational_OnGaussPoints( cH_OnGPs, mResultsOnGPs );
            
            for ( int gp=0; gp<cH_OnGPs.size(); gp++ ) 
            {
                cH_OnGPs[gp] *= cCteMat;
            }
            
            Vector2Type rH_OnGPs, iH_OnGPs;
            
            Split_Complex_VectorVector( rH_OnGPs, iH_OnGPs, cH_OnGPs );
            
            REAL_H_OnGP.push_back( rH_OnGPs ); 
            IMAG_H_OnGP.push_back( iH_OnGPs ); 
        }
	}

	//**********************************************************************************************************************
    //* - Calculate B field on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_B_Element( int* NodesId, unsigned int PropertiesId, 
                                       ResultsOnGPsType& REAL_B_OnGP, 
                                       ResultsOnGPsType& IMAG_B_OnGP )
	{
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        if ( mSmoothing_On ) 
        {   
            double mo      = ( 4.00e-7 )  *  3.141592653589793238462643383279; 
            double mu_real = (*properties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		    double mu_imag = (*properties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

            std::complex<double> cMu  ( mu_real, mu_imag );
            std::complex<double> cUnit(    0.00,    1.00 );
            std::complex<double> cB_OnNodes;

            Vector<double> rH_OnNodes( 3 );
            Vector<double> iH_OnNodes( 3 );

            Vector<double> rB_OnNodes( 3 );
            Vector<double> iB_OnNodes( 3 );

            for ( int n=0; n<nodes.size(); n++ )
            {
                rH_OnNodes = (*mpModel)( REAL_H, *nodes[n] );
                iH_OnNodes = (*mpModel)( IMAG_H, *nodes[n] );

                for ( int i=0; i<3; i++ )
                {
                    cB_OnNodes    = cMu * ( rH_OnNodes[i] + cUnit * iH_OnNodes[i] );
                    rB_OnNodes[i] = std::real( cB_OnNodes );
                    iB_OnNodes[i] = std::imag( cB_OnNodes );
                }

                (*mpModel)( REAL_B, *nodes[n] ) += rB_OnNodes;
                (*mpModel)( IMAG_B, *nodes[n] ) += iB_OnNodes;
                
                mNeighbourElements[ nodes[n]->Id() - 1 ]++;
            }
        }
        else
        {
            Element::Pointer pElement;
            
            if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            std::complex<double> cUnit( 0.00, 1.00 );
            std::complex<double> cCteMat;
 
            if ( mPotentials_On ) cCteMat = 1.00;
            else                  cCteMat = 1.00 / ( cUnit * mProblemFrequency );

            cVector2Type cB_OnGPs;
            
            pElement->Calculate_Rotational_OnGaussPoints( cB_OnGPs, mResultsOnGPs );
               
            for ( int gp=0; gp<cB_OnGPs.size(); gp++ ) 
            {
                cB_OnGPs[gp] *= cCteMat;
            }
            
            Vector2Type rB_OnGPs, iB_OnGPs;
            
            Split_Complex_VectorVector( rB_OnGPs, iB_OnGPs, cB_OnGPs );

            REAL_B_OnGP.push_back( rB_OnGPs ); 
            IMAG_B_OnGP.push_back( iB_OnGPs ); 
        }
	}

    //********************************************************************************************************************************
    //* - Calculates electrostatic E field E = -grad(V)
    //********************************************************************************************************************************
    void Modeler::Calculate_E_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );

        Vector2Type gradV;

        if ( mSmoothing_On )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator node_it;
            
            int n = 0;
            
            for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
            {
                (*mpModel)(ELECTROSTATIC_ELECTRIC_FIELD, **node_it) -= gradV[n];
            
                mNeighbourElements[ (*node_it)->Id() - 1 ]++;

                n++;
            }
        }
        else
        {
            pElement->Calculate_Gradient_OnGaussPoints( gradV, mResultsOnGPs );
            
            Vector2Type Electrostatic_E_OnCentralGP = gradV;

            Electrostatic_E_OnCentralGP[0] *= -1.00;
            
            ResultsOnGPs.push_back( Electrostatic_E_OnCentralGP );
        }
    }

    //********************************************************************************************************************************
    //* - Calculates electrostatic current J = -sigma*grad(V)
    //********************************************************************************************************************************
    void Modeler::Calculate_J_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );

        Vector2Type gradV;
         
        if ( mSmoothing_On )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator node_it;
            
            int n = 0;
            
            for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
            {
                Vector<double> temp_J = gradV[n];
            
                temp_J *= sigma;
            
                (*mpModel)(ELECTROSTATIC_CURRENT_DENSITY, **node_it) -= temp_J;
            
                mNeighbourElements[ (*node_it)->Id() - 1 ]++;
            
                n++;
            }
        }
        else
        {
            pElement->Calculate_Gradient_OnGaussPoints( gradV, mResultsOnGPs );
            
            Vector2Type Electrostatic_J_OnCentralGP = gradV;
            
            Electrostatic_J_OnCentralGP[0] *= -sigma;
            
            ResultsOnGPs.push_back( Electrostatic_J_OnCentralGP );
        }
    }

    //********************************************************************************************************************************
    //* - It writes element nodes Ids and the calculated electrostatic current on the export currents file.
    //********************************************************************************************************************************
    void Modeler::Export_Electrostatic_J_Element_GP( int* NodesId, unsigned int PropertiesId )
    {
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

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
            nodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );

        Vector2Type gradOnGaussPoints;

        pElement->Calculate_Gradient_OnGaussPoints( gradOnGaussPoints, mResultsOnGPs );

        // J = -sigma*gradV in the central Gauss point
        Vector<double> Electrostatic_J_OnCentralGP = gradOnGaussPoints[0];

        Electrostatic_J_OnCentralGP *= -sigma;

        double modJ2 = Electrostatic_J_OnCentralGP[0] * Electrostatic_J_OnCentralGP[0] + 
                       Electrostatic_J_OnCentralGP[1] * Electrostatic_J_OnCentralGP[1] +
                       Electrostatic_J_OnCentralGP[2] * Electrostatic_J_OnCentralGP[2] ;

        // If |J| == 0.0 do not write the result in the export file
        if ( modJ2 <= 0.0 ) return;

        // Adding element nodes to the export file
        mExportCurrentFile << NodesId[0] << "  " << NodesId[1] << "  " 
                           << NodesId[2] << "  " << NodesId[3] << "  "; 
          
        // Adding current source components to the export file
        mExportCurrentFile << Electrostatic_J_OnCentralGP[0] << "  " 
                           << Electrostatic_J_OnCentralGP[1] << "  " 
                           << Electrostatic_J_OnCentralGP[2] ;

        // End of line
        mExportCurrentFile << std::endl; 
    }

    //********************************************************************************************************************************
    //* - It opens the export current file and writes the phase of the current.
    //********************************************************************************************************************************
    void Modeler::Ini_Exporting_Electrostatic_Current()
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

    //********************************************************************************************************************************
    //* - It closes the export current file.
    //********************************************************************************************************************************
    void Modeler::End_Exporting_Electrostatic_Current()
    {
        mExportCurrentFile.close();
    }

    //********************************************************************************************************************************
    //* - Add to the b_vector the J currents from the folder "Exp_J_Sources/"
    //********************************************************************************************************************************
    void Modeler::Build_Imported_Electrostatic_Currents()
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
            double J_phase    ;

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

	//********************************************************************************************************************************
    //* - Assembling imported current element
    //********************************************************************************************************************************
    void Modeler::Assemble_Imported_Current_Element( int* NodesId, double* Jvec, double Jpha )
    {
        // Get volume element nodes
		Vector<Node::Pointer> pNodes;

        Vector<int> HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Cartesian J current source
        Vector<double> Cartesian_J_Source( 6 );
        
        Cartesian_J_Source[0] = Jvec[0]; Cartesian_J_Source[1] = Jpha;
        Cartesian_J_Source[2] = Jvec[1]; Cartesian_J_Source[3] = Jpha;
        Cartesian_J_Source[4] = Jvec[2]; Cartesian_J_Source[5] = Jpha;

        // Axysimmetric J current source
        Vector<double> Axysimmetric_J_Source(2);

        Axysimmetric_J_Source[0] = 0.0; 
        Axysimmetric_J_Source[1] = 0.0;

        // Element properties 
        Properties::Pointer Properties( new Properties() );

        // Source properties
 	    PropertyFunction< Vector<double> >::Pointer    Cartesian_J( new ConstantProperty< Vector<double> >(    Cartesian_J_Source ) );
        PropertyFunction< Vector<double> >::Pointer Axysimmetric_J( new ConstantProperty< Vector<double> >( Axysimmetric_J_Source ) );

	    Properties->SetProperty( SINUSOIDAL_SURFACE_CURRENT,    Cartesian_J );       
        Properties->SetProperty( COMPLEX_IBC               , Axysimmetric_J );  

        // Problem frequency
        PropertyFunction< double >::Pointer ProblemFrequency( new ConstantProperty< double >( mProblemFrequency ) );

        Properties->SetProperty( FREQUENCY, ProblemFrequency );   

        // Volumetric source element definition
        Element::Pointer pElement;

        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, Properties ) );
	    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, Properties ) );
	    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, Properties ) );
	    
        // Activates potentials on element
        pElement->SetPotentials( mPotentials_On );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental residual vector
        Vector<std::complex<double> > EleResVector; 

        pElement->GetResidualVector( EleResVector );

        // Apply boundary conditions to force vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Cleaning vectors 
        EleIdsVector.FreeData();
        EleResVector.FreeData();
    }

    //********************************************************************************************************************************
    //* - Calculates Joule heating Q = sigma * ||grad(V)||^2 
    //********************************************************************************************************************************
    void Modeler::Calculate_JouleH_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        double sigma = (*properties)(IHL_ELECTRIC_CONDUCTIVITY);

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( nodes, properties ) );

        Vector2Type gradV;

        double gradV2, JouleHeating;

        if ( mSmoothing_On == true )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator node_it;
            
            int n = 0;
            
            for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
            {
                gradV2       = gradV[n][0]*gradV[n][0] + gradV[n][1]*gradV[n][1] + gradV[n][2]*gradV[n][2];
                JouleHeating = sigma * gradV2;

                (*mpModel)(ELECTROSTATIC_JOULE_HEATING, **node_it) += JouleHeating;
            
                mNeighbourElements[ (*node_it)->Id() - 1 ]++;
                n++;
            }        
        }
        else
        {
            pElement->Calculate_Gradient_OnGaussPoints( gradV, mResultsOnGPs );
            
            gradV2       = gradV[0][0]*gradV[0][0] + gradV[0][1]*gradV[0][1] + gradV[0][2]*gradV[0][2];
            JouleHeating = sigma * gradV2;

            Vector<double> JouleHeating_OnCentralGP; 
            JouleHeating_OnCentralGP.push_back( JouleHeating );

            Vector2Type JouleHeating_OnGPs; 
            JouleHeating_OnGPs.push_back( JouleHeating_OnCentralGP );
            
            ResultsOnGPs.push_back( JouleHeating_OnGPs );
        }
    }

    //********************************************************************************************************************************
    //* - Calculate Joule heating.
    //********************************************************************************************************************************
    void Modeler::Calculate_Joule_Heating_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP )
    {
        if ( mSmoothing_On )
        {
            Calculate_Joule_Heating_OnNodes( NodesId, PropertiesId );
        }
        else
        {
            Calculate_Joule_Heating_OnGaussPoints( NodesId, PropertiesId, JOULE_HEATING_OnGP );
        }
    }

    //********************************************************************************************************************************
    //* - Initiates Joule heating calculation
    //********************************************************************************************************************************
    void Modeler::Ini_Joule_Heating_Calculation()
    {
        if ( !mSmoothing_On ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

	    std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
	        (*mpModel)( JOULE_HEATING, **it ) = 0.00;
        }    

        mNeighbourElements.resize( nodes.size() ); 
        
        for ( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
    }

    //********************************************************************************************************************************
    //* - Ends Joule heating calculation
    //********************************************************************************************************************************
    void Modeler::End_Joule_Heating_Calculation()
    {
		if ( !mSmoothing_On ) return;
        
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        double Joule_Heating;
        
        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            Joule_Heating = (*mpModel)( JOULE_HEATING, **it );
           
            if ( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
            {
                Joule_Heating /= mNeighbourElements[ (*it)->Id()-1 ];
            }

            (*mpModel)( JOULE_HEATING, **it ) = Joule_Heating;           
        }    
    }

	//********************************************************************************************************************************
    //* - Calculate Joule heating on nodes
    //********************************************************************************************************************************
	void Modeler::Calculate_Joule_Heating_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
		// Joule heating in IHL material with losses
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		double eo   = 8.8541878176e-12;
		double freq = mProblemFrequency;

		double sigma    = (*properties)(  IHL_ELECTRIC_CONDUCTIVITY );
		double eps_imag = (*properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		double eff_cond = freq * eps_imag + sigma;

        if ( eff_cond > 0 )
        {
		    std::vector<Node::Pointer> nodes;
            
            Get_Volume_Element_Nodes( NodesId, nodes );

            this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

            Element::Pointer pElement;
            
            if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );
            
            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cE;
            
            pElement->Calculate_E_field_OnNodes( cE );

            Vector2Type rE, iE;
            
            Split_Complex_VectorVector( rE, iE, cE );

            double modE2;

            for ( int n=0; n<nodes.size(); n++ )
            {
                modE2 = rE[n][0]*rE[n][0] + iE[n][0]*iE[n][0] +
		    			rE[n][1]*rE[n][1] + iE[n][1]*iE[n][1] +
		    			rE[n][2]*rE[n][2] + iE[n][2]*iE[n][2] ;

                (*mpModel)( JOULE_HEATING, *nodes[n] ) += 0.5 * eff_cond * modE2;
                
                mNeighbourElements[ nodes[n]->Id() - 1 ]++;
            }
        }
	}

    //********************************************************************************************************************************
    //* - Calculate Joule heating on Gauss points
    //********************************************************************************************************************************
	void Modeler::Calculate_Joule_Heating_OnGaussPoints( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP )
	{
        int numGaussPoints = mResultsOnGPs;

        if ( ( mElementOrder < 2 ) && ( mResultsOnGPs > 4 ) ) numGaussPoints = 4;

        Vector2Type Joule_Heating_OnElement( numGaussPoints );

        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            Joule_Heating_OnElement[gp].resize(1);
            Joule_Heating_OnElement[gp][0] = 0.00;
        }

        // Joule heating in IHL material with losses
		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		double eo   = 8.8541878176e-12;
		double freq = mProblemFrequency;

		double sigma    = (*properties)( IHL_ELECTRIC_CONDUCTIVITY  );
		double eps_imag = (*properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		double eff_cond = freq * eps_imag + sigma;

        if ( eff_cond > 0 )
        {
            std::vector<Node::Pointer> nodes;

            Get_Volume_Element_Nodes( NodesId, nodes );

            Element::Pointer pElement;

            if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cE;

            pElement->Calculate_E_field_OnGaussPoints( cE, mResultsOnGPs );

            VectorType rE( 3 ), iE( 3 );

            double modE2;

            for ( int gp=0; gp<numGaussPoints; gp++ ) 
            {
                rE[0] = std::real( cE[gp][0] );
                rE[1] = std::real( cE[gp][1] );
                rE[2] = std::real( cE[gp][2] );

                iE[0] = std::imag( cE[gp][0] );
                iE[1] = std::imag( cE[gp][1] );
                iE[2] = std::imag( cE[gp][2] );

		    	modE2 = rE[0]*rE[0] + rE[1]*rE[1] + rE[2]*rE[2] + 
                        iE[0]*iE[0] + iE[1]*iE[1] + iE[2]*iE[2] ;
            
		    	Joule_Heating_OnElement[gp][0] = 0.5 * eff_cond * modE2;
            }
        }

        JOULE_HEATING_OnGP.push_back( Joule_Heating_OnElement );
    }

    //**********************************************************************************************************************
    //* - Calculate J induced on nodes (IHL materials)
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Induced_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		double sigma = (*properties)( IHL_ELECTRIC_CONDUCTIVITY );

        if ( sigma > 0.00 ) 
        {
		    std::vector<Node::Pointer> nodes;

            Get_Volume_Element_Nodes( NodesId, nodes );

            this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

            Element::Pointer pElement;
            
            if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );
            
            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cJ_OnNodes;
            
            pElement->Calculate_E_field_OnNodes( cJ_OnNodes );

            for ( int n=0; n<nodes.size(); n++ ) cJ_OnNodes[n] *= sigma;
            
            Vector2Type rJ_OnNodes, iJ_OnNodes;
            
            Split_Complex_VectorVector( rJ_OnNodes, iJ_OnNodes, cJ_OnNodes );
            
		    for ( int n=0; n<nodes.size(); n++ )
            {
                (*mpModel)( REAL_J, *nodes[n] ) += rJ_OnNodes[n];
                (*mpModel)( IMAG_J, *nodes[n] ) += iJ_OnNodes[n];

                mNeighbourElements[ nodes[n]->Id() - 1 ]++;
            }
        }
    }

    //**********************************************************************************************************************
    //* - Calculate J induced on Gauss points (IHL materials)
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Induced_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_eddy_OnGP )
    {
 		int numGaussPoints = mResultsOnGPs;

        if ( ( mElementOrder < 2 ) && ( mResultsOnGPs > 4 ) ) numGaussPoints = 4;

        cJ_eddy_OnGP.resize( numGaussPoints );

        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            cJ_eddy_OnGP[gp].resize( 3 );
            cJ_eddy_OnGP[gp][0] = std::complex<double>( 0.00, 0.00 );
            cJ_eddy_OnGP[gp][1] = std::complex<double>( 0.00, 0.00 );
            cJ_eddy_OnGP[gp][2] = std::complex<double>( 0.00, 0.00 );
        }

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		double sigma = (*properties)( IHL_ELECTRIC_CONDUCTIVITY );

        if ( sigma > 0.00 ) 
        {
		    std::vector<Node::Pointer> nodes;

            Get_Volume_Element_Nodes( NodesId, nodes );

            Element::Pointer pElement;

            if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
		    else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cE;

            pElement->Calculate_E_field_OnGaussPoints( cE, mResultsOnGPs );

            for ( int gp=0; gp<numGaussPoints; gp++ ) 
            {
                cJ_eddy_OnGP[gp][0] = sigma * cE[gp][0];
                cJ_eddy_OnGP[gp][1] = sigma * cE[gp][1];
                cJ_eddy_OnGP[gp][2] = sigma * cE[gp][2];
            }
        } 
    }

	//**********************************************************************************************************************
    //* - Calculate J imposed on nodes
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Imposed_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
		Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

        std::vector<double> vJ = (*properties)( SINUSOIDAL_SURFACE_CURRENT );

		std::complex<double> cJx( vJ[0] * cos( vJ[1] ), vJ[0] * sin( vJ[1] ) );
		std::complex<double> cJy( vJ[2] * cos( vJ[3] ), vJ[2] * sin( vJ[3] ) );
		std::complex<double> cJz( vJ[4] * cos( vJ[5] ), vJ[4] * sin( vJ[5] ) );

		std::vector<double> vJa   = (*properties)( COMPLEX_IBC          );
		std::vector<double> Jaxis = (*properties)( COMPLEX_NEUMANN_FLOW );

		std::complex<double> cJa( vJa[0] * cos( vJa[1] ), vJa[0] * sin( vJa[1] ) );

		std::vector<Node::Pointer> nodes;

        Get_Volume_Element_Nodes( NodesId, nodes );

		Vector<Node::Pointer>::iterator node_it;

		for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
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

            Vector<double> J_vector( 3 );

			int nodeId = (*node_it)->Id();

			J_vector[0] = std::real( cJx ) + std::real( cJa )*rx;
            J_vector[1] = std::real( cJy ) + std::real( cJa )*ry;
            J_vector[2] = std::real( cJz ) + std::real( cJa )*rz;

			tmp_rJ[nodeId] = J_vector;

			J_vector[0] = std::imag( cJx ) + std::imag( cJa )*rx;
            J_vector[1] = std::imag( cJy ) + std::imag( cJa )*ry;
            J_vector[2] = std::imag( cJz ) + std::imag( cJa )*rz;

			tmp_iJ[nodeId] = J_vector;
        }
	}

	//**********************************************************************************************************************
    //* - Add imported currents to element nodes
    //**********************************************************************************************************************
	void Modeler::Calculate_Imported_Electrostatic_Currents()
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
            int            NodesId[4] ;
            double         J_phase    ;
            Vector<double> J_vector(3);

            // Dummy vectors
            Vector<double> rtmp(3);
            Vector<double> itmp(3);

            // Counter of neighbords elements for smoothing
            std::map<unsigned int, int> numNeighbors;

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

                // Get high order nodes
                std::vector<int> HONodesId;

                PushHONodesOnVolume( NodesId, HONodesId );

                // Projecting imported J on element nodes
                for ( int i=0; i<HONodesId.size(); i++ )
                {
                    int nodeId = HONodesId[i];

                    // Check for the first occurrence of the node
                    if ( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() )
                    {
                        for ( int j=0; j<3; j++ )
                        {
                            tmp_ist_rJ[nodeId][j] += J_vector[j] * cos( J_phase ); 
                            tmp_ist_iJ[nodeId][j] += J_vector[j] * sin( J_phase ); 
                        }

                        numNeighbors[nodeId]++;
                    }
                    else
                    {
                        for ( int j=0; j<3; j++ )
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

                for ( int j=0; j<3; j++ )
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

	//**********************************************************************************************************************
    //* - Calculate J imposed/imported on Gauss points
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Imposed_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_imps_OnGP )
    {
 		// Complex current density imposed on Gauss points
        int numGaussPoints = mResultsOnGPs;

        if ( ( mElementOrder < 2 ) && ( mResultsOnGPs > 4 ) ) numGaussPoints = 4;

        // cJ vector on Gauss points
        cJ_imps_OnGP.resize( numGaussPoints );

        // Initializes cJ vector on Gauss points
        for ( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            cJ_imps_OnGP[gp].resize( 3 );
            cJ_imps_OnGP[gp][0] = std::complex<double>( 0.00, 0.00 );
            cJ_imps_OnGP[gp][1] = std::complex<double>( 0.00, 0.00 );
            cJ_imps_OnGP[gp][2] = std::complex<double>( 0.00, 0.00 );
        }
        
        // Number of vertices of a tetrahedral element
        int numVerticeNodes = 4;

        // cJ vector on nodes
        cVector2Type cJ_OnNodes( numVerticeNodes );

        // Initializes cJ vector on nodes
        for ( int n=0; n<numVerticeNodes; n++ ) 
        {
            cJ_OnNodes[n].resize( 3 );
            cJ_OnNodes[n][0] = std::complex<double>( 0.00, 0.00 );
            cJ_OnNodes[n][1] = std::complex<double>( 0.00, 0.00 );
            cJ_OnNodes[n][2] = std::complex<double>( 0.00, 0.00 );
        }

        // Check if the element is inside an imposed J source volume
        bool Element_Is_Inside_ImposedJ_Source;

        if ( ( tmp_rJ.find( NodesId[0] ) != tmp_rJ.end() ) &&
             ( tmp_rJ.find( NodesId[1] ) != tmp_rJ.end() ) &&
             ( tmp_rJ.find( NodesId[2] ) != tmp_rJ.end() ) &&
             ( tmp_rJ.find( NodesId[3] ) != tmp_rJ.end() )  )
        {
            Element_Is_Inside_ImposedJ_Source = true;
        }
        else 
        {
            Element_Is_Inside_ImposedJ_Source = false;
        }

        // If the element is inside an imposed J source volume then add cJ values on vertices
        if ( Element_Is_Inside_ImposedJ_Source )
        {
            for ( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_OnNodes[n][0] += std::complex<double>( tmp_rJ[NodesId[n]][0], tmp_iJ[NodesId[n]][0] );
                cJ_OnNodes[n][1] += std::complex<double>( tmp_rJ[NodesId[n]][1], tmp_iJ[NodesId[n]][1] );
                cJ_OnNodes[n][2] += std::complex<double>( tmp_rJ[NodesId[n]][2], tmp_iJ[NodesId[n]][2] );
            }
        }

        // Check if the element is inside an imported J source volume and outside an imposed J source
        bool Element_Is_Inside_ImportedJ_Source;

        if ( ( Element_Is_Inside_ImposedJ_Source == false        ) && 
             ( tmp_ist_rJ.find( NodesId[0] ) != tmp_ist_rJ.end() ) &&
             ( tmp_ist_rJ.find( NodesId[1] ) != tmp_ist_rJ.end() ) &&
             ( tmp_ist_rJ.find( NodesId[2] ) != tmp_ist_rJ.end() ) &&
             ( tmp_ist_rJ.find( NodesId[3] ) != tmp_ist_rJ.end() )  )
        {
            Element_Is_Inside_ImportedJ_Source = true;
        }
        else
        {
            Element_Is_Inside_ImportedJ_Source = false;
        }
         
        // If the element is inside an imported J source volume then add cJ values on vertices
        if ( Element_Is_Inside_ImportedJ_Source )
        {
            for( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_OnNodes[n][0] += std::complex<double>( tmp_ist_rJ[NodesId[n]][0], tmp_ist_iJ[NodesId[n]][0] );
                cJ_OnNodes[n][1] += std::complex<double>( tmp_ist_rJ[NodesId[n]][1], tmp_ist_iJ[NodesId[n]][1] );
                cJ_OnNodes[n][2] += std::complex<double>( tmp_ist_rJ[NodesId[n]][2], tmp_ist_iJ[NodesId[n]][2] );
            }
        }

        // If the element is not inside an imposed or an imported J source then finish 
        if ( ( Element_Is_Inside_ImposedJ_Source == false ) && ( Element_Is_Inside_ImportedJ_Source == false ) )
        {
            return;
        }

        // If the element is inside an imported/imposed J source then interpolate nodal values to Gauss points
        std::vector<Node::Pointer> nodes;
        
        Get_Volume_Element_Nodes( NodesId, nodes );

        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );
        
        Element::Pointer pElement;
        
        if      ( mElementOrder == 0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( nodes, properties ) );
        else if ( mElementOrder == 1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( nodes, properties ) );
        else if ( mElementOrder == 2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( nodes, properties ) );
        
        std::vector<double> cX, cY, cZ; 
        
        pElement->GetInnerGiDGaussPoints( cX, cY, cZ, mResultsOnGPs );
        
        Matrix<double> N; 
        
        Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // Interpolating nodal values to Gauss points
        for( int gp=0; gp<numGaussPoints; gp++ ) 
        {
            for( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_imps_OnGP[gp][0] += N[n][gp] * cJ_OnNodes[n][0];
                cJ_imps_OnGP[gp][1] += N[n][gp] * cJ_OnNodes[n][1];
                cJ_imps_OnGP[gp][2] += N[n][gp] * cJ_OnNodes[n][2];
            }
        }
    }

	//**********************************************************************************************************************
    //* - Calculate J on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Element( int* NodesId, unsigned int PropertiesId, 
                                       ResultsOnGPsType& REAL_J_OnGP, 
                                       ResultsOnGPsType& IMAG_J_OnGP )
	{
        Calculate_J_Induced_OnNodes( NodesId, PropertiesId );
        
        if ( mSmoothing_On ) 
        {
            return;
        }

        cVector2Type cJ_eddy_OnGP;

        Calculate_J_Induced_OnGaussPoints( NodesId, PropertiesId, cJ_eddy_OnGP );

        cVector2Type cJ_imps_OnGP;

        Calculate_J_Imposed_OnGaussPoints( NodesId, PropertiesId, cJ_imps_OnGP );

        cVector2Type cJ_tot_OnGP( cJ_eddy_OnGP.size() );

        for( int gp=0; gp<cJ_tot_OnGP.size(); gp++ ) 
        {
            cJ_tot_OnGP[gp].resize( 3 );
            cJ_tot_OnGP[gp][0] = cJ_eddy_OnGP[gp][0] + cJ_imps_OnGP[gp][0];
            cJ_tot_OnGP[gp][1] = cJ_eddy_OnGP[gp][1] + cJ_imps_OnGP[gp][1];
            cJ_tot_OnGP[gp][2] = cJ_eddy_OnGP[gp][2] + cJ_imps_OnGP[gp][2];
        }

        Vector2Type rJ_tot_OnGP, iJ_tot_OnGP;
        
        Split_Complex_VectorVector( rJ_tot_OnGP, iJ_tot_OnGP, cJ_tot_OnGP );
        
        REAL_J_OnGP.push_back( rJ_tot_OnGP ); 
        IMAG_J_OnGP.push_back( iJ_tot_OnGP ); 
    }

	//**********************************************************************************************************************
    //* - Adding all the current contributions ( eddy, imposed and imported )
    //**********************************************************************************************************************
	void Modeler::End_J_Calculation()
	{
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		std::vector<Node::Pointer>::iterator node_it;

		Vector<double> rEddy( 3, 0.00 ); 
        Vector<double> iEddy( 3, 0.00 );

		Vector<double> rJtot( 3, 0.00 ); 
        Vector<double> iJtot( 3, 0.00 );

        for( node_it = nodes.begin(); node_it != nodes.end(); ++node_it )
        {
            int nodeId = (*node_it)->Id();

            // Eddy currents contribution
			rEddy = (*mpModel)( REAL_J, **node_it );
            iEddy = (*mpModel)( IMAG_J, **node_it );

            if ( mNeighbourElements[ nodeId-1 ] > 0 ) 
            {
                rEddy /= mNeighbourElements[ nodeId-1 ];
                iEddy /= mNeighbourElements[ nodeId-1 ];
            }
        
            for( int i=0; i<3; i++ )
            {
                rJtot[i] = rEddy[i];
                iJtot[i] = iEddy[i]; 
            }

            // Adding imposed currents contribution
			if ( tmp_rJ.find( nodeId ) != tmp_rJ.end() )
			{
                for( int i=0; i<3; i++ )
                {
                    rJtot[i] += tmp_rJ[nodeId][i];
                    iJtot[i] += tmp_iJ[nodeId][i];
                }
			}

            // Adding imported electrostatic currents contribution
            if ( ( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() ) && 
                 (     tmp_rJ.find( nodeId ) ==     tmp_rJ.end() )  )
			{
                for( int i=0; i<3; i++ )
                {
                    rJtot[i] += tmp_ist_rJ[nodeId][i];
                    iJtot[i] += tmp_ist_iJ[nodeId][i];
                }
			}

            // Total current distribution
            (*mpModel)( REAL_J, **node_it ) = rJtot;
			(*mpModel)( IMAG_J, **node_it ) = iJtot;
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

    //*************************************************************************************************************
    //* - Calculate Sij parameters
    //*************************************************************************************************************
    void Modeler::CalculateSijParameters()
	{
		std::map<double,  std::complex<double> >::iterator proyec_it;

		for( proyec_it = mProjection.begin(); proyec_it != mProjection.end(); proyec_it++ )
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

    //*****************************************************************************************************************
    //* - Calculates of nodal magnitudes
    //*****************************************************************************************************************
    void Modeler::CalculateNodal( int variableKey )
    {
		std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();
        std::vector<Node::Pointer>::iterator it;

        // E(t)
        if ( variableKey == E.getKey() )
        {
            Vector<double> E_total( 3 );
            Vector<double> E_real;
            Vector<double> E_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );

                E_total[0] = E_real[0] * cos( freq * time ) + E_imag[0] * sin( freq * time );                                                  
                E_total[1] = E_real[1] * cos( freq * time ) + E_imag[1] * sin( freq * time );                                                   
                E_total[2] = E_real[2] * cos( freq * time ) + E_imag[2] * sin( freq * time );

                (*mpModel)(E, **it) = E_total;
            }
        }
		// B(t)
        else if ( variableKey == B.getKey() )
        {
            Vector<double> B_total( 3 );
            Vector<double> B_real;
            Vector<double> B_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                B_real = (*mpModel)( REAL_B, **it );
                B_imag = (*mpModel)( IMAG_B, **it );

                B_total[0] = B_real[0] * cos( freq * time ) + B_imag[0] * sin( freq * time );
                B_total[1] = B_real[1] * cos( freq * time ) + B_imag[1] * sin( freq * time );
                B_total[2] = B_real[2] * cos( freq * time ) + B_imag[2] * sin( freq * time );

                (*mpModel)(B, **it) = B_total;
            }
        }
		// H(t)
        else if ( variableKey == H.getKey() )
        {
            Vector<double> H_total( 3 );
            Vector<double> H_real;
            Vector<double> H_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                H_real = (*mpModel)( REAL_H, **it );
                H_imag = (*mpModel)( IMAG_H, **it );

                H_total[0] = H_real[0] * cos( freq * time ) + H_imag[0] * sin( freq * time );
                H_total[1] = H_real[1] * cos( freq * time ) + H_imag[1] * sin( freq * time );
                H_total[2] = H_real[2] * cos( freq * time ) + H_imag[2] * sin( freq * time );

                (*mpModel)(H, **it) = H_total;
            }
        }
		// J(t)
        else if ( variableKey == J.getKey() )
        {
            Vector<double> J_total( 3 );
            Vector<double> J_real;
            Vector<double> J_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                J_real = (*mpModel)( REAL_J, **it );
                J_imag = (*mpModel)( IMAG_J, **it );

                J_total[0] = J_real[0] * cos( freq * time ) + J_imag[0] * sin( freq * time );
                J_total[1] = J_real[1] * cos( freq * time ) + J_imag[1] * sin( freq * time );
                J_total[2] = J_real[2] * cos( freq * time ) + J_imag[2] * sin( freq * time );

                (*mpModel)(J, **it) = J_total;
            }
        }
		// MOD_E
        else if ( variableKey == MOD_E.getKey() )
        {
		    Vector<double> E_real( 3 );
            Vector<double> E_imag( 3 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );
               
                (*mpModel)(MOD_E, **it) = sqrt( E_real[0]*E_real[0] + E_imag[0]*E_imag[0] +
                                                E_real[1]*E_real[1] + E_imag[1]*E_imag[1] +
                                                E_real[2]*E_real[2] + E_imag[2]*E_imag[2] );
            }
        }
		// MOD_H
        else if ( variableKey == MOD_H.getKey() )
        {
		    Vector<double> H_real( 3 );
            Vector<double> H_imag( 3 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                H_real = (*mpModel)( REAL_H, **it );
                H_imag = (*mpModel)( IMAG_H, **it );
               
                (*mpModel)(MOD_H, **it) = sqrt( H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                                                H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                                                H_real[2]*H_real[2] + H_imag[2]*H_imag[2] );
            }
        }
		// MOD_J
        else if ( variableKey == MOD_J.getKey() )
        {
		    Vector<double> J_real( 3 );
            Vector<double> J_imag( 3 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                J_real = (*mpModel)( REAL_J, **it );
                J_imag = (*mpModel)( IMAG_J, **it );
               
                (*mpModel)(MOD_J, **it) = sqrt( J_real[0]*J_real[0] + J_imag[0]*J_imag[0] +
                                                J_real[1]*J_real[1] + J_imag[1]*J_imag[1] +
                                                J_real[2]*J_real[2] + J_imag[2]*J_imag[2] );
            }
        }
		// MOD_B
        else if ( variableKey == MOD_B.getKey() )
        {
		    Vector<double> B_real( 3 );
            Vector<double> B_imag( 3 );

		    for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                B_real = (*mpModel)( REAL_B, **it );
                B_imag = (*mpModel)( IMAG_B, **it );
               
                (*mpModel)(MOD_B, **it) = sqrt( B_real[0]*B_real[0] + B_imag[0]*B_imag[0] +
                                                B_real[1]*B_real[1] + B_imag[1]*B_imag[1] +
                                                B_real[2]*B_real[2] + B_imag[2]*B_imag[2] );
            }
        }
		// BOUNDARY_NORMALS
        else if ( variableKey == BOUNDARY_NORMALS.getKey() )
        {
            for ( it=nodes.begin(); it!= nodes.end(); ++it )
            {
			    Vector<double> BoundaryNormal;

                Get_Boundary_Normal( BoundaryNormal, (*it)->Id() );

                (*mpModel)( BOUNDARY_NORMALS, **it ) = BoundaryNormal;
            }
        }
        // CONTACT_NORMALS
        else if ( variableKey == CONTACT_NORMALS.getKey() )
        {
            for ( it=nodes.begin(); it!= nodes.end(); ++it )
            {
			    Vector<double> ContactPNormal;

                Get_ContactP_Normal( ContactPNormal, (*it)->Id() );

                (*mpModel)( CONTACT_NORMALS, **it ) = ContactPNormal;
            }
        }
        // E field parallel to B_ext
        else if ( variableKey == E_PARALLEL_T.getKey() )
        {
            Vector<double> E_total( 3 );
            Vector<double> E_real;
            Vector<double> E_imag;

            Vector<double> E_parallel_t( 3 );

            Vector<double> b( 3 );
            Vector<double> Bext;

            double Bnorm;
            double EparDotb;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                // External magnetic flux density
                Bext = mpColdPlasma->Get_Bexternal_InNode( *it );

                // b = B / |B|
                Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

                b[0] = Bext[0] / Bnorm;
                b[1] = Bext[1] / Bnorm;
                b[2] = Bext[2] / Bnorm;

                // E field
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );

                E_total[0] = E_real[0] * cos( freq * time ) + E_imag[0] * sin( freq * time );
                E_total[1] = E_real[1] * cos( freq * time ) + E_imag[1] * sin( freq * time );
                E_total[2] = E_real[2] * cos( freq * time ) + E_imag[2] * sin( freq * time );

                // E field proyection on B_ext
                EparDotb = E_total[0]*b[0] + E_total[1]*b[1] + E_total[2]*b[2];

                E_parallel_t[0] = EparDotb * b[0];
                E_parallel_t[1] = EparDotb * b[1];
                E_parallel_t[2] = EparDotb * b[2];

                (*mpModel)( E_PARALLEL_T, **it ) = E_parallel_t;
            }
        }
        // E field perpendicular to B_ext
        else if ( variableKey == E_PERPENDICULAR_T.getKey() )
        {
            Vector<double> E_total( 3 );
            Vector<double> E_real;
            Vector<double> E_imag;

            Vector<double> E_perpendicular_t( 3 );

            Vector<double> b( 3 );
            Vector<double> Bext;

            double Bnorm;  
            double EparDotb;     

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                // External magnetic flux density
                Bext = mpColdPlasma->Get_Bexternal_InNode( *it );

                // b = B / |B|
                Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

                b[0] = Bext[0] / Bnorm;
                b[1] = Bext[1] / Bnorm;
                b[2] = Bext[2] / Bnorm;

                // E field
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );

                E_total[0] = E_real[0] * cos( freq * time ) + E_imag[0] * sin( freq * time );
                E_total[1] = E_real[1] * cos( freq * time ) + E_imag[1] * sin( freq * time );
                E_total[2] = E_real[2] * cos( freq * time ) + E_imag[2] * sin( freq * time );

                // E field proyection on the perpendicular direction to B_ext (Eper = E - Epar)
                EparDotb = E_total[0]*b[0] + E_total[1]*b[1] + E_total[2]*b[2];

                E_perpendicular_t[0] = E_total[0] - EparDotb * b[0];
                E_perpendicular_t[1] = E_total[1] - EparDotb * b[1];
                E_perpendicular_t[2] = E_total[2] - EparDotb * b[2];

                (*mpModel)( E_PERPENDICULAR_T, **it ) = E_perpendicular_t;
            }
        }
        // Poynting vector
        else if ( variableKey == POYNTING_VECTOR.getKey() )
        {
            std::complex<double> cUnit( 0.00, 1.00 );

            Vector<double> E_real, E_imag;
            Vector<double> H_real, H_imag;
            Vector<double> avgS( 3 );

            Vector< std::complex<double> > cEv( 3 );
            Vector< std::complex<double> > cHv( 3 );

            Vector< std::complex<double> > E_X_conj_H(3);

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                // E field
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );

                cEv[0] = E_real[0] + cUnit * E_imag[0];
                cEv[1] = E_real[1] + cUnit * E_imag[1];
                cEv[2] = E_real[2] + cUnit * E_imag[2];

                // H field
                H_real = (*mpModel)( REAL_H, **it );
                H_imag = (*mpModel)( IMAG_H, **it );

                cHv[0] = H_real[0] + cUnit * H_imag[0];
                cHv[1] = H_real[1] + cUnit * H_imag[1];
                cHv[2] = H_real[2] + cUnit * H_imag[2];

                // E X conj(H)
                E_X_conj_H[0] = cEv[1] * std::conj( cHv[2] ) - cEv[2] * std::conj( cHv[1] );
                E_X_conj_H[1] = cEv[2] * std::conj( cHv[0] ) - cEv[0] * std::conj( cHv[2] );
                E_X_conj_H[2] = cEv[0] * std::conj( cHv[1] ) - cEv[1] * std::conj( cHv[0] );

                // < S > = 0.5 * real( E X conj(H) )
                avgS[0] = 0.5 * std::real( E_X_conj_H[0] );
                avgS[1] = 0.5 * std::real( E_X_conj_H[1] );
                avgS[2] = 0.5 * std::real( E_X_conj_H[2] );

                (*mpModel)( POYNTING_VECTOR, **it ) = avgS;
            }
        }
        // Lorentz force
        else if ( variableKey == LORENTZ_FORCE.getKey() )
        {
            std::complex<double> cUnit( 0.00, 1.00 );

            Vector<double> J_real, J_imag;
            Vector<double> B_real, B_imag;
            Vector<double> avgF( 3 );

            Vector< std::complex<double> > cJv( 3 );
            Vector< std::complex<double> > cBv( 3 );

            Vector< std::complex<double> > J_X_conj_B( 3 );

            for ( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                // J total
                J_real = (*mpModel)( REAL_J, **it );
                J_imag = (*mpModel)( IMAG_J, **it );

                cJv[0] = J_real[0] + cUnit * J_imag[0];
                cJv[1] = J_real[1] + cUnit * J_imag[1];
                cJv[2] = J_real[2] + cUnit * J_imag[2];

                // H field
                B_real = (*mpModel)( REAL_B, **it );
                B_imag = (*mpModel)( IMAG_B, **it );

                cBv[0] = B_real[0] + cUnit * B_imag[0];
                cBv[1] = B_real[1] + cUnit * B_imag[1];
                cBv[2] = B_real[2] + cUnit * B_imag[2];

                // J X conj(B)
                J_X_conj_B[0] = cJv[1] * std::conj( cBv[2] ) - cJv[2] * std::conj( cBv[1] );
                J_X_conj_B[1] = cJv[2] * std::conj( cBv[0] ) - cJv[0] * std::conj( cBv[2] );
                J_X_conj_B[2] = cJv[0] * std::conj( cBv[1] ) - cJv[1] * std::conj( cBv[0] );

                // < F > = 0.5 * real( J X conj(B) )
                avgF[0] = 0.5 * std::real( J_X_conj_B[0] );
                avgF[1] = 0.5 * std::real( J_X_conj_B[1] );
                avgF[2] = 0.5 * std::real( J_X_conj_B[2] );

                (*mpModel)( LORENTZ_FORCE, **it ) = avgF;
            }
        }
    }
} 
