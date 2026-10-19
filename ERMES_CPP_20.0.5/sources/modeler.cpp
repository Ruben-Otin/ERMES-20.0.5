
#include "../ERMES/VolumeElement_1eg_FullWave.h"
#include "../ERMES/VolumeElement_2eg_FullWave.h"
#include "../ERMES/VolumeElement_1st_FullWave.h"
#include "../ERMES/VolumeElement_2nd_FullWave.h"
#include "../ERMES/VolumeElement_3sb_FullWave.h"

#include "../ERMES/VolumeElement_1eg_ColdPlasma.h"
#include "../ERMES/VolumeElement_2eg_ColdPlasma.h"
#include "../ERMES/VolumeElement_1st_ColdPlasma.h"
#include "../ERMES/VolumeElement_2nd_ColdPlasma.h"
#include "../ERMES/VolumeElement_3sb_ColdPlasma.h"
#include "../ERMES/VolumeElement_1st_Electrostatic.h"

#include "../ERMES/JSource_1eg_FullWave.h"
#include "../ERMES/JSource_2eg_FullWave.h"
#include "../ERMES/JSource_1st_FullWave.h"
#include "../ERMES/JSource_2nd_FullWave.h"
#include "../ERMES/JSource_3sb_FullWave.h"

#include "../ERMES/GenericRobin_1eg_FullWave.h"
#include "../ERMES/GenericRobin_2eg_FullWave.h"
#include "../ERMES/GenericRobin_1st_FullWave.h"
#include "../ERMES/GenericRobin_2nd_FullWave.h"
#include "../ERMES/GenericRobin_3sb_FullWave.h"
#include "../ERMES/GenericRobin_1st_Electrostatic.h"

#include "../ERMES/FarField_1eg_FullWave.h"
#include "../ERMES/FarField_2eg_FullWave.h"
#include "../ERMES/FarField_1st_FullWave.h"
#include "../ERMES/FarField_2nd_FullWave.h"
#include "../ERMES/FarField_3sb_FullWave.h"

#include "../ERMES/FarField_1eg_ColdPlasma.h"
#include "../ERMES/FarField_2eg_ColdPlasma.h"
#include "../ERMES/FarField_1st_ColdPlasma.h"
#include "../ERMES/FarField_2nd_ColdPlasma.h"
#include "../ERMES/FarField_3sb_ColdPlasma.h"

#include "../ERMES/CoaxialPortTEM_1eg_FullWave.h"
#include "../ERMES/CoaxialPortTEM_2eg_FullWave.h"
#include "../ERMES/CoaxialPortTEM_1st_FullWave.h"
#include "../ERMES/CoaxialPortTEM_2nd_FullWave.h"
#include "../ERMES/CoaxialPortTEM_3sb_FullWave.h"

#include "../ERMES/RWPortTE10_1eg_FullWave.h"
#include "../ERMES/RWPortTE10_2eg_FullWave.h"
#include "../ERMES/RWPortTE10_1st_FullWave.h"
#include "../ERMES/RWPortTE10_2nd_FullWave.h"
#include "../ERMES/RWPortTE10_3sb_FullWave.h"

#include "../ERMES/Utils/MatrixFileWritting.h"
#include "../ERMES/Utils/LagrangeBases.h"
#include "../ERMES/Utils/GaussIntegrationTables.h"

#include "../includes/modeler.h"
#include "../includes/model.h"
#include "../linear_solvers/kratos_linear_solvers.h"

#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <stdlib.h>
#include <omp.h>

#ifdef _WIN64
    #include <direct.h>
    #include "dirent/include/dirent.h"
#else
    #include <dirent.h>
    #include <sys/stat.h>
#endif

namespace Kratos
{
    //*************************************************************************************************
    //* - Generates nodes.
    //*************************************************************************************************
    void Modeler::Generate_Node( IndexType NodeId, double X, double Y, double Z )
    {
        Node::Pointer new_node( new Node( NodeId, X, Y, Z ) );
    
        mpModel->AddNode( new_node );
    
        if ( mElectrostaticMode )
        {
            new_node->Add_V_Dof( mDofSet );
        }
        else if ( mPotentials ) 
        {
            new_node->Add_cAV_Dofs( mDofSet );
        }
        else   
        {
            new_node->Add_cE_Dofs( mDofSet );
        }
    }
    
    //*************************************************************************************************
    //* - Prints error messages and exit ERMES.
    //*************************************************************************************************
    void Modeler::Send_Error_Msg( String ErrorMsg, String CheckMsg, int VerticalSpaces )
    {
         for( int i=0; i<VerticalSpaces; i++ )
         {
             std::cout << std::endl;
         }
         std::cout << "!!!!!!!!!! " << ErrorMsg << " !!!!!!!!!!" << std::endl;
         std::cout << std::endl;
         std::cout << CheckMsg << std::endl;
         std::cout << std::endl;
         std::cout << "ERMES analysis finished." << std::endl;

         // Also report on stderr and return a failure status, so batch systems (SLURM afterok, scripts) see the error
         std::cerr << "ERMES error: " << ErrorMsg << std::endl;

         exit( EXIT_FAILURE );
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////
    // 
	// - Contact discontinuity in dielectrics surface.
    // 
	///////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************
    //* - Fills nodes properties vector map.
    //*************************************************************************************************
	void Modeler::Contact_Properties( int* NodesId, unsigned int PropertiesId )
	{
		 for( int i=0; i<4; i++ ) 
         {
             mTempNodeProperties[ NodesId[ i ] ] = PropertiesId;
         }
	}

	//*************************************************************************************************
    //* - Generates pairs of nodes in a contact surface.
    //*************************************************************************************************
	void Modeler::Generate_Contact_Pairs( Vector<int>& CtCNodesId )
	{
        if( mElectrostaticMode == true ) 
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
		double eo   = 8.8541878176e-12 ;

		double real_sgm_First3 = (*First3_prop)( REAL_ELECTRIC_CONDUCTIVITY );
        double imag_sgm_First3 = (*First3_prop)( IMAG_ELECTRIC_CONDUCTIVITY );

		double real_eps_First3 = (*First3_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double imag_eps_First3 = (*First3_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> cEps_First3( real_eps_First3 - ( imag_sgm_First3 / freq ), imag_eps_First3 + ( real_sgm_First3 / freq ) );

		double mod_epc_First3 = std::abs(cEps_First3);

		double real_sgm_Last3 = (*Last3_prop)( REAL_ELECTRIC_CONDUCTIVITY );
        double imag_sgm_Last3 = (*Last3_prop)( IMAG_ELECTRIC_CONDUCTIVITY );

		double real_eps_Last3 = (*Last3_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
		double imag_eps_Last3 = (*Last3_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		std::complex<double> cEps_Last3( real_eps_Last3 - ( imag_sgm_Last3 / freq ), imag_eps_Last3 + ( real_sgm_Last3 / freq ) );

		double mod_epc_Last3 = std::abs( cEps_Last3 );

		double mod_epc_R;
        double mod_epc_L;

		if( mod_epc_First3 > mod_epc_Last3 )
		{
			Vector<int> OneSideNodesId( 3 );

			for( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i ];
            }

            Push_HONodes_OnSurface( OneSideNodesId, sR_NodesId );

			for( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i + 3 ];
            }

            Push_HONodes_OnSurface( OneSideNodesId, sL_NodesId );

			Calculate_Contact_Normals( sR_NodesId );
			Calculate_Contact_Normals( sL_NodesId );

			R_PropertyId = mTempNodeProperties[ CtCNodesId[0] ];
			L_PropertyId = mTempNodeProperties[ CtCNodesId[3] ];

			mod_epc_R = mod_epc_First3;
			mod_epc_L = mod_epc_Last3 ;
		}
		else
		{
			Vector<int> OneSideNodesId( 3 );

			for( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i ];
            }
			
            Push_HONodes_OnSurface( OneSideNodesId, sL_NodesId );

			for( int i=0; i<3; i++ ) 
            {
                OneSideNodesId[ i ] = CtCNodesId[ i + 3 ];
            }
			
            Push_HONodes_OnSurface( OneSideNodesId, sR_NodesId );

			Calculate_Contact_Normals( sL_NodesId );
			Calculate_Contact_Normals( sR_NodesId );

			L_PropertyId = mTempNodeProperties[ CtCNodesId[0] ];
			R_PropertyId = mTempNodeProperties[ CtCNodesId[3] ];

			mod_epc_L = mod_epc_First3;
			mod_epc_R = mod_epc_Last3;
		}

        // Setting temporal contact pairs
        if( mPotentials )
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

		for( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			mTempNodeProperties[ (*R_Idit) ] = R_PropertyId;
		    mTempNodeProperties[ (*L_Idit) ] = L_PropertyId;
        }
	}

	//**********************************************************************************************************************
    //* - Sets temporal contact pairs in E field formulation.
    //**********************************************************************************************************************
	void Modeler::Set_TempContactPairs_E( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L )
    {
		Vector<int>::const_iterator R_Idit = sR_NodesId.begin();
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			int RNodeId = (*R_Idit);
			int LNodeId = (*L_Idit);

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

			if( pRNode->pDofcEx()->IsFixed() )
			{
				if( mTempContactPairs[ RNodeId ] != LNodeId )
				{
					int currLNodeId = mTempContactPairs[ RNodeId ];

					Properties::Pointer currL_prop = mpModel->GetProperties( mTempNodeProperties[ currLNodeId ] );

                    double freq = mProblemFrequency;
		            double eo   = 8.8541878176e-12 ;

					double real_sgm_currL = (*currL_prop)( REAL_ELECTRIC_CONDUCTIVITY );
                    double imag_sgm_currL = (*currL_prop)( IMAG_ELECTRIC_CONDUCTIVITY );

					double real_eps_currL = (*currL_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
					double imag_eps_currL = (*currL_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

					std::complex<double> cEps_currL( real_eps_currL - ( imag_sgm_currL / freq ), imag_eps_currL + ( real_sgm_currL / freq ) );

					double mod_epc_currL = std::abs( cEps_currL );

					if( mod_epc_currL > mod_epc_L )
					{
						Node::Pointer pcurrLNode = mpModel->GetNode( currLNodeId-1 );

						pcurrLNode->pDofcEx()->FixDof();

						pcurrLNode->pDofcEy()->FreeDof();

						pLNode->pDofcEy()->FixDof();

						mTempContactPairs[ RNodeId     ] = LNodeId;
						mTempContactPairs[ currLNodeId ] = LNodeId;

						std::map<unsigned int, unsigned int>::iterator cPair_it;

						for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
						{
							if( cPair_it->second == currLNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
						}
					}
					else
					{
						pLNode->pDofcEx()->FixDof();

						mTempContactPairs[ LNodeId ] = currLNodeId;
					}
				}
			}
			else if( pLNode->pDofcEx()->IsFixed() )
			{
				if( !pRNode->pDofcEy()->IsFixed() )
				{
					pRNode->pDofcEx()->FixDof();

					mTempContactPairs[ RNodeId ] = mTempContactPairs[ LNodeId ];
				}
			}
			else if( pRNode->pDofcEy()->IsFixed() )
			{
				pRNode->pDofcEx()->FixDof();

			    pRNode->pDofcEy()->FreeDof();

				pLNode->pDofcEy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;

				std::map<unsigned int, unsigned int>::iterator cPair_it;

		        for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
				{
					if( cPair_it->second == RNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
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
    //* - Sets temporal contact pairs in potential A formulation.
    //**********************************************************************************************************************
	void Modeler::Set_TempContactPairs_A( Vector<int>& sR_NodesId, Vector<int>& sL_NodesId, double mod_epc_L )
    {
		Vector<int>::const_iterator R_Idit = sR_NodesId.begin();
		Vector<int>::const_iterator L_Idit = sL_NodesId.begin();

		for( R_Idit = sR_NodesId.begin(); R_Idit != sR_NodesId.end(); ++R_Idit, ++L_Idit )
		{
			int RNodeId = (*R_Idit);
			int LNodeId = (*L_Idit);

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

			if( pRNode->pDofcAx()->IsFixed() )
			{
				if( mTempContactPairs[ RNodeId ] != LNodeId )
				{
					int currLNodeId = mTempContactPairs[ RNodeId ];

					Properties::Pointer currL_prop = mpModel->GetProperties( mTempNodeProperties[ currLNodeId ] );

                    double freq = mProblemFrequency;
		            double eo   = 8.8541878176e-12 ;

					double real_sgm_currL = (*currL_prop)( REAL_ELECTRIC_CONDUCTIVITY );
                    double imag_sgm_currL = (*currL_prop)( IMAG_ELECTRIC_CONDUCTIVITY );

					double real_eps_currL = (*currL_prop)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
					double imag_eps_currL = (*currL_prop)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

					std::complex<double> cEps_currL( real_eps_currL - ( imag_sgm_currL / freq ), imag_eps_currL + ( real_sgm_currL / freq ) );

					double mod_epc_currL = std::abs( cEps_currL );

					if( mod_epc_currL > mod_epc_L )
					{
						Node::Pointer pcurrLNode = mpModel->GetNode( currLNodeId-1 );

						pcurrLNode->pDofcAx()->FixDof();

						pcurrLNode->pDofcAy()->FreeDof();

						pLNode->pDofcAy()->FixDof();

						mTempContactPairs[ RNodeId     ] = LNodeId;
						mTempContactPairs[ currLNodeId ] = LNodeId;

						std::map<unsigned int, unsigned int>::iterator cPair_it;

						for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
						{
							if( cPair_it->second == currLNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
						}
					}
					else
					{
						pLNode->pDofcAx()->FixDof();

						mTempContactPairs[ LNodeId ] = currLNodeId;
					}
				}
			}
			else if( pLNode->pDofcAx()->IsFixed() )
			{
				if( !pRNode->pDofcAy()->IsFixed() )
				{
					pRNode->pDofcAx()->FixDof();

					mTempContactPairs[ RNodeId ] = mTempContactPairs[ LNodeId ];
				}
			}
			else if( pRNode->pDofcAy()->IsFixed() )
			{
				pRNode->pDofcAx()->FixDof();

			    pRNode->pDofcAy()->FreeDof();

				pLNode->pDofcAy()->FixDof();

				mTempContactPairs[ RNodeId ] = LNodeId;

				std::map<unsigned int, unsigned int>::iterator cPair_it;

		        for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
				{
					if( cPair_it->second == RNodeId ) mTempContactPairs[ cPair_it->first ] = LNodeId;
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
    //* - Calculates geometric averaged normal in contact surfaces.
    //**********************************************************************************************************************
	void Modeler::Calculate_Contact_Normals( Vector<int>& NodesId )
	{
		Node::Pointer pNode;

		Vector<Node::Pointer> nodes;

		Vector<int>::const_iterator id_it;

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
    //* - Sets average normals in contact surfaces.
    //**********************************************************************************************************************
	void Modeler::Set_AvgNormals_Contact()
	{
		// Setting average normals in surface dilelectric nodes
		std::map< unsigned int, Vector< Vector<double> > >::iterator it_Set;

        for ( it_Set = mSetContactNormals.begin(); it_Set != mSetContactNormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            Vector< Vector<double> > vectorSet( it_Set->second );

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
        std::map< unsigned int, Vector< Vector<double> > >().swap( mSetContactNormals );
	}

	//**********************************************************************************************************************
    //* - Sets discontinuities on dielectric surfaces.
    //**********************************************************************************************************************
	void Modeler::Set_Contacts()
	{
        Set_AvgNormals_Contact();

		std::map<unsigned int, unsigned int>::iterator cPair_it;

		for( cPair_it = mTempContactPairs.begin(); cPair_it != mTempContactPairs.end(); ++cPair_it )
        {
			unsigned int RNodeId = cPair_it->first;
			unsigned int LNodeId = mTempContactPairs[ RNodeId ];

			Node::Pointer pRNode = mpModel->GetNode( RNodeId-1 );
			Node::Pointer pLNode = mpModel->GetNode( LNodeId-1 );

            std::complex<double> cZero( 0.00, 0.00 );

            if( mPotentials )
            {
                (*mpModel)( cAx, *pRNode ) = cZero;
			    (*mpModel)( cAy, *pRNode ) = cZero;
			    (*mpModel)( cAz, *pRNode ) = cZero;
                (*mpModel)( cVs, *pRNode ) = cZero;

                pRNode->pDofcAx()->FixDof();
                pRNode->pDofcAy()->FixDof();
			    pRNode->pDofcAz()->FixDof();

                if( pRNode->pDof( cVs ) != NULL ) 
                {
                    pRNode->pDofcVs()->FixDof(); 
                }

			    pLNode->pDofcAx()->FreeDof();
			    pLNode->pDofcAy()->FreeDof();
			    pLNode->pDofcAz()->FreeDof();

                if( pLNode->pDof( cVs ) != NULL ) 
                {
                    pLNode->pDofcVs()->FreeDof(); 
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
    //* - Deletes temporal vectors generated along contact discontinuity making.
    //**********************************************************************************************************************
	void Modeler::DelTempContactVectors()
	{
		mTempNodeProperties.clear();
		std::map<unsigned int, unsigned short int>().swap( mTempNodeProperties );

		mContactNormals.clear();
		std::map<unsigned int, Vector<double> >().swap( mContactNormals );

		mTempContactPairs.clear();
		std::map<unsigned int, unsigned int>().swap( mTempContactPairs );
	}

	//**********************************************************************************************************************
    //* - Joins normals to second node from the first node.
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

				Vector< Vector<double> > FirstVectorSet ( mSetPECnormals[ FirstNode ] );

				Vector< Vector<double> >::iterator it_Fvs;

				for ( it_Fvs = FirstVectorSet.begin(); it_Fvs != FirstVectorSet.end(); ++it_Fvs )
				{
					Vector<double> vectorF( *it_Fvs );

					mSetPECnormals[ SecondNode ].push_back( vectorF );
				}
			}
		}
	}

	//**********************************************************************************************************************
    //* - Forces Dirichlet normals in all contact nodes that are in a Dirichlet surface.
    //**********************************************************************************************************************
	void Modeler::ForceDirichletNormalsInContact()
	{
		std::map<unsigned int, ContactPairData>::iterator itCtc;

		// New dielectric normal. All nodes in a Dirichlet surface have the same normal.
		// PEC / PMC / TEC dominates normal direction.
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
    // 
	// - Fill singularity maps.
    // 
	////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //**************************************************************************************************
    //* - Fills mNSingular2L map.
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
    //* - Fills mNSingular3L map.
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
    //* - Fills mNSingular4L map.
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
    //* - Fills mNSingular5L map.
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
    //* - Fills mNSingular6L map.
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

	//**************************************************************************************************
    //* - Collapses bubbles for the LL2P formulation. 
    //**************************************************************************************************
	void Modeler::Collapse_LL2P_Bubbles( int* NodesId )
	{
        if( mElementOrder != 0 ) 
        {
            return;
        }
        
        Vector<Node::Pointer> pNodes;
            
        Get_Volume_Element_Nodes( NodesId, pNodes );
            
        bool IsUgLayer = false;

        for( int i=0; i<4; i++ )
		{
            int NodeId = pNodes[ i ]->Id();

			if( ( mSingular   .find( NodeId ) != mSingular   .end() ) || ( mNSingular2L.find( NodeId ) != mNSingular2L.end() ) ||
			    ( mNSingular3L.find( NodeId ) != mNSingular3L.end() ) || ( mNSingular4L.find( NodeId ) != mNSingular4L.end() ) ||
			    ( mNSingular5L.find( NodeId ) != mNSingular5L.end() ) || ( mNSingular6L.find( NodeId ) != mNSingular6L.end() )  )
			{
			    IsUgLayer = true;
			}
		}

        std::complex<double> cZero( 0.0, 0.0 );

        if( IsUgLayer )
        {
            if( !mPotentials )
            {
                for( int i=4; i<17; i++ )
                {
                    (*mpModel)( cEx, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEx()->FixDof(); 
                    (*mpModel)( cEy, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEy()->FixDof(); 
                    (*mpModel)( cEz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEz()->FixDof(); 
                }
            }
            else
            {
                for( int i=4; i<17; i++ )
                {
                    (*mpModel)( cAx, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAx()->FixDof(); 
                    (*mpModel)( cAy, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAy()->FixDof(); 
                    (*mpModel)( cAz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAz()->FixDof(); 
                }
            }
        }
    }

    //**************************************************************************************************
    //* - Clears mSingular maps.
    //**************************************************************************************************
	void Modeler::Clear_Singularities( )
	{
        mSingular   .clear(); std::map<unsigned int, double>().swap( mSingular    );
        mNSingular2L.clear(); std::map<unsigned int, double>().swap( mNSingular2L );
        mNSingular3L.clear(); std::map<unsigned int, double>().swap( mNSingular3L );
        mNSingular4L.clear(); std::map<unsigned int, double>().swap( mNSingular4L );
        mNSingular5L.clear(); std::map<unsigned int, double>().swap( mNSingular5L );
        mNSingular6L.clear(); std::map<unsigned int, double>().swap( mNSingular6L );
    }

	////////////////////////////////////////////////////////////////////////////////////////////////////
    // 
	// - High Order Elements.
    // 
	////////////////////////////////////////////////////////////////////////////////////////////////////

	//**************************************************************************************************
    //* - Adds high order elements to the HONodesId vector in volume elements.
    //**************************************************************************************************
	void Modeler::Push_HONodes_OnVolume( int* NodesId, Vector<int>& HONodesId )
	{
		// Add corner nodes to any element order 
		HONodesId.push_back( NodesId[ 0 ] );
        HONodesId.push_back( NodesId[ 1 ] );
        HONodesId.push_back( NodesId[ 2 ] );
		HONodesId.push_back( NodesId[ 3 ] );

		// LL2P 3sb elements
        if( mElementOrder == 0 ) 
        {
            HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 1 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );

            HONodes_OnVolume( NodesId, HONodesId );
        }
		// RME 1st order 
		else if( mElementOrder == 1 ) 
		{
            return;
		}
		// RME 2nd order 
		else if( mElementOrder == 2 ) 
		{
            HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		    HONodes_OnEdge( NodesId[ 0 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		}
		// RME 3rd order 
		else if( mElementOrder == 3 ) 
		{
            HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		    HONodes_OnEdge( NodesId[ 0 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 3 ], HONodesId );
                                     
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 1 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		}
		// RME 4th order 
		else if( mElementOrder == 4 ) 
		{
            HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		    HONodes_OnEdge( NodesId[ 0 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 3 ], HONodesId );
						 
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 1 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );

		    HONodes_OnVolume( NodesId, HONodesId );
		}
		// Edge element 1st order 
		else if( mElementOrder == 11 ) 
		{
            HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		    HONodes_OnEdge( NodesId[ 0 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		}
		// Edge element 2nd order 
		else if( mElementOrder == 12 ) 
		{
            HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		    HONodes_OnEdge( NodesId[ 0 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 3 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 3 ], HONodesId );

		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 1 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		    HONodes_OnFace( NodesId[ 0 ], NodesId[ 2 ], NodesId[ 3 ], HONodesId );
		}
	}

	//**************************************************************************************************
    //* - Adds high order elements to the HONodesId vector in volume elements.
    //**************************************************************************************************
	void Modeler::Push_HONodes_OnVolume( Vector<int>& NodesId, Vector<int>& HONodesId )
	{
		int aNodesId[] = { NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], NodesId[ 3 ] };

		Push_HONodes_OnVolume( aNodesId, HONodesId );
	}

	//*************************************************************************************************
    //* - Adds high order elements to the HONodesId vector in surface elements.
    //*************************************************************************************************
	void Modeler::Push_HONodes_OnSurface( int* NodesId, Vector<int>& HONodesId )
	{
		// Add corner nodes to any element order 
		HONodesId.push_back( NodesId[ 0 ] );
        HONodesId.push_back( NodesId[ 1 ] );
        HONodesId.push_back( NodesId[ 2 ] );

		// LL2P 3sb elements
        if( mElementOrder == 0 ) 
        {
            HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
        }
		// RME 1st order 
		else if( mElementOrder == 1 ) 
		{
            return;
		}
		// RME 2nd order 
		else if( mElementOrder == 2 ) 
		{
			HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		}
		// RME 3rd order 
		else if( mElementOrder == 3 ) 
		{
			HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );

			HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		}
		// RME 4th order 
		else if( mElementOrder == 4 ) 
		{
			HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );

			HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		}
		// Edge element 1st order 
		else if( mElementOrder == 11 ) 
		{
			HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );
		}
		// Edge element 2nd order 
		else if( mElementOrder == 12 ) 
		{
			HONodes_OnEdge( NodesId[ 0 ], NodesId[ 1 ], HONodesId );
            HONodes_OnEdge( NodesId[ 1 ], NodesId[ 2 ], HONodesId );
            HONodes_OnEdge( NodesId[ 2 ], NodesId[ 0 ], HONodesId );

			HONodes_OnFace( NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ], HONodesId );
		}
	}

	//*****************************************************************************************************************
    //* - Adds high order elements to the HONodesId vector in surface elements.
    //*****************************************************************************************************************
	void Modeler::Push_HONodes_OnSurface( Vector<int>& NodesId, Vector<int>& HONodesId )
	{
		int aNodesId[] = { NodesId[ 0 ], NodesId[ 1 ], NodesId[ 2 ] };

		Push_HONodes_OnSurface( aNodesId, HONodesId );
	}

	//*****************************************************************************************************************
    //* - Looks for high order nodes on edges.
    //*****************************************************************************************************************
    void Modeler::HONodes_OnEdge( int IdNode1, int IdNode2, Vector<int>& HONodesId )
    {
        bool inverseOrder;

        unsigned int minId, maxId;

        if( IdNode1 < IdNode2 )
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

        Vector< Vector<unsigned int> > SameKeyEdges( mEdgeHONodes[ minId*maxId ] );

        Vector< Vector<unsigned int> >::iterator it_edges;

        for( it_edges = SameKeyEdges.begin(); it_edges != SameKeyEdges.end(); ++it_edges )
        {
            Vector<unsigned int> edge( *it_edges );

            if( ( edge[0] == minId ) && ( edge[1] == maxId ) )
            {
                if( !inverseOrder )
                {
                    Vector<unsigned int>::iterator it_hon = edge.begin()+2;

                    while( it_hon != edge.end() )
                    {
                         HONodesId.push_back( *it_hon );

                         it_hon++;
                    }
                }
                else
                {
                    Vector<unsigned int>::reverse_iterator rit_hon = edge.rbegin();

                    while( rit_hon != edge.rend()-2 )
                    {
                        HONodesId.push_back( *rit_hon );

                        rit_hon++;
                    }
                }
            }
        }
    }

	//*****************************************************************************************************************
    //* - Looks for high order nodes on faces.
    //*****************************************************************************************************************
	void Modeler::HONodes_OnFace( int IdNode1, int IdNode2, int IdNode3, Vector<int>& HONodesId )
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

		Vector< Vector<unsigned int> > SameKeyFaces( mFaceHONodes[ minId*medId*maxId ] );

        Vector< Vector<unsigned int> >::iterator it_faces;

		for( it_faces = SameKeyFaces.begin(); it_faces != SameKeyFaces.end(); ++it_faces )
        {
            Vector<unsigned int> face( *it_faces );

            if( ( face[0] == minId ) && ( face[1] == medId ) && ( face[2] == maxId ) )
            {
			    if     ( mElementOrder ==  0 ) Reorder_HONodes_OnFace_4th( orderingCase, face, HONodesId );
                else if( mElementOrder ==  3 ) HONodesId.push_back       ( face[3]                       );
			    else if( mElementOrder ==  4 ) Reorder_HONodes_OnFace_4th( orderingCase, face, HONodesId );
				else if( mElementOrder == 12 ) HONodesId.push_back       ( face[3]                       );
            }
        }
	}

	//*****************************************************************************************************************
    //* - Reorders nodes for a giving surface order.
    //*****************************************************************************************************************
	void Modeler::Reorder_HONodes_OnFace_4th( int orderingCase, Vector<unsigned int>& face, Vector<int>& HONodesId )
	{
		unsigned int minf = face[ 3 ];
		unsigned int medf = face[ 4 ];
		unsigned int maxf = face[ 5 ];
        
		//[n1,n2,n3] = [minId, medId, maxId]
		if( orderingCase == 1 )
		{
		    HONodesId.push_back( minf );
		    HONodesId.push_back( medf );
		    HONodesId.push_back( maxf );
		}
		//[n1,n3,n2] = [minId, maxId, medId]
		else if( orderingCase == 2 )
		{
		    HONodesId.push_back( minf );
		    HONodesId.push_back( maxf );
		    HONodesId.push_back( medf );
		}
		//[n3,n1,n2] = [medId, maxId, minId]
		else if( orderingCase == 3 )
		{
		    HONodesId.push_back( medf );
		    HONodesId.push_back( maxf );
		    HONodesId.push_back( minf );
		}
		//[n2,n1,n3] = [medId, minId, maxId]
		else if( orderingCase == 4 )
		{
		    HONodesId.push_back( medf );
		    HONodesId.push_back( minf );
		    HONodesId.push_back( maxf );
		}
		//[n2,n3,n1] = [maxId, minId, medId]
		else if( orderingCase == 5 )
		{
		    HONodesId.push_back( maxf );
		    HONodesId.push_back( minf );
		    HONodesId.push_back( medf );
		}
		//[n3,n2,n1] = [maxId, medId, minId]
		else if( orderingCase == 6 )
		{
		    HONodesId.push_back( maxf );
		    HONodesId.push_back( medf );
		    HONodesId.push_back( minf );
	    }
	}

	//*************************************************************************************************
    //* - Looks for high order nodes on volumes.
    //*************************************************************************************************
	void Modeler::HONodes_OnVolume( int* NodesId, Vector<int>& HONodesId )
	{
	    unsigned int key = NodesId[0]*NodesId[1]*NodesId[2]*NodesId[3];
        
        Vector< Vector<unsigned int> > SameKeyElement( mVolmHONodes[key] );

        Vector< Vector<unsigned int> >::iterator it_elements;
        
        for( it_elements = SameKeyElement.begin(); it_elements != SameKeyElement.end(); ++it_elements )
        {
            Vector<unsigned int> nInElement( *it_elements );
        
            if( ( nInElement[0] == NodesId[0] ) && ( nInElement[1] == NodesId[1] ) && 
                ( nInElement[2] == NodesId[2] ) && ( nInElement[3] == NodesId[3] ) )
            {
                HONodesId.push_back( nInElement[4] );
            }
        }
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
    // 
	// - High order nodes creation.
    // 
	///////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************
    //* - Creates high order nodes.
    //*************************************************************************************************
    void Modeler::Create_HONodes( int* NodesId )
    {
        // LL2P 3sb element
        if( mElementOrder == 0 )
        {
            Create_HONodes_OnFaces( NodesId );
            Create_HONodes_OnVolum( NodesId );
        }
        // RME 1st order element
        else if( mElementOrder == 1 )
        {
            return;
        }
        // RME 2nd order element
        else if( mElementOrder == 2 )
        {
            Create_HONodes_OnEdges( NodesId );
        }
        // RME 3rd order element
        else if( mElementOrder == 3 )
        {
            Create_HONodes_OnEdges( NodesId );
            Create_HONodes_OnFaces( NodesId );
        }
        // RME 4th order element
        else if( mElementOrder == 4 )
        {
            Create_HONodes_OnEdges( NodesId );
            Create_HONodes_OnFaces( NodesId );
            Create_HONodes_OnVolum( NodesId );
        }
        // Edge element 1st order
        else if( mElementOrder == 11 )
        {
            Create_HONodes_OnEdges( NodesId );
        }
		// Edge element 2nd order
        else if( mElementOrder == 12 )
        {
            Create_HONodes_OnEdges( NodesId );
			Create_HONodes_OnFaces( NodesId );
        }
	}

	//*************************************************************************************************
    //* - Creates high order nodes on edges.
    //*************************************************************************************************
	void Modeler::Create_HONodes_OnEdges( int* NodesId )
	{
		unsigned int IdNode1, IdNode2;

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

				if( mEdgeHONodes.find( key ) != mEdgeHONodes.end() )
				{
					Vector< Vector<unsigned int> > ExistantNodesOnEdges( mEdgeHONodes[ key ] );

					Vector< Vector<unsigned int> >::iterator it_ee;

					int repe = 0;

					for( it_ee = ExistantNodesOnEdges.begin(); it_ee != ExistantNodesOnEdges.end(); ++it_ee )
					{
						Vector<unsigned int> NodesOnEdge(*it_ee);

						if( ( IdNode1 == NodesOnEdge[0] ) && ( IdNode2 == NodesOnEdge[1] ) ) repe++;
					}

					if( repe == 0 )
					{
						if     ( mElementOrder ==  2 ) Create_HONodes_OnEdge_2nd( IdNode1, IdNode2 );
						else if( mElementOrder ==  3 ) Create_HONodes_OnEdge_3rd( IdNode1, IdNode2 );
						else if( mElementOrder ==  4 ) Create_HONodes_OnEdge_4th( IdNode1, IdNode2 );
                        else if( mElementOrder == 11 ) Create_HONodes_OnEdge_2nd( IdNode1, IdNode2 );
						else if( mElementOrder == 12 ) Create_HONodes_OnEdge_3rd( IdNode1, IdNode2 );
					}
				}
				else
				{
					if     ( mElementOrder ==  2 ) Create_HONodes_OnEdge_2nd( IdNode1, IdNode2 );
					else if( mElementOrder ==  3 ) Create_HONodes_OnEdge_3rd( IdNode1, IdNode2 );
					else if( mElementOrder ==  4 ) Create_HONodes_OnEdge_4th( IdNode1, IdNode2 );
                    else if( mElementOrder == 11 ) Create_HONodes_OnEdge_2nd( IdNode1, IdNode2 );
					else if( mElementOrder == 12 ) Create_HONodes_OnEdge_3rd( IdNode1, IdNode2 );
				}
			}
		}
	}

	//*************************************************************************************************
    //* - Creates high order nodes on faces.
    //*************************************************************************************************
	void Modeler::Create_HONodes_OnFaces( int* NodesId )
	{
		unsigned int IdNode1, IdNode2, IdNode3;
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

					if ( mFaceHONodes.find( key ) != mFaceHONodes.end() )
					{
						Vector< Vector<unsigned int> > ExistantNodesOnFaces( mFaceHONodes[ key ] );

						Vector< Vector<unsigned int> >::iterator it_ef;

						int repe = 0;

						for( it_ef = ExistantNodesOnFaces.begin(); it_ef != ExistantNodesOnFaces.end(); ++it_ef )
						{
							Vector<unsigned int> NodesOnFace(*it_ef);

							if( ( IdNode1==NodesOnFace[0] ) && ( IdNode2==NodesOnFace[1] ) && ( IdNode3==NodesOnFace[2] ) ) repe++;
						}

						if( repe == 0 )
						{
							if     ( mElementOrder ==  0 ) Create_HONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
                            else if( mElementOrder ==  3 ) Create_HONodes_OnFace_3rd( IdNode1, IdNode2, IdNode3 );
							else if( mElementOrder ==  4 ) Create_HONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
							else if( mElementOrder == 12 ) Create_HONodes_OnFace_3rd( IdNode1, IdNode2, IdNode3 );
						}
					}
					else
					{
					    if     ( mElementOrder ==  0 ) Create_HONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
                        else if( mElementOrder ==  3 ) Create_HONodes_OnFace_3rd( IdNode1, IdNode2, IdNode3 );
					    else if( mElementOrder ==  4 ) Create_HONodes_OnFace_4th( IdNode1, IdNode2, IdNode3 );
						else if( mElementOrder == 12 ) Create_HONodes_OnFace_3rd( IdNode1, IdNode2, IdNode3 );
					}
				}
			}
		}
	}

	//*************************************************************************************************
    //* - Creates new higher order nodes on edges ( 2nd order ).
    //*************************************************************************************************
    void Modeler::Create_HONodes_OnEdge_2nd( int IdNode1, int IdNode2 )
    {
		Vector<unsigned int> edgeNodes( 3 );

        edgeNodes[ 0 ] = IdNode1;
        edgeNodes[ 1 ] = IdNode2;

        edgeNodes[ 2 ] = ( mpModel->GetNodesArray() ).size() + 1;

        mEdgeHONodes[ IdNode1*IdNode2 ].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1 - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2 - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        // Creating new high order node on edge number 1
        alpha = 0.5;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[ 2 ], hoX, hoY, hoZ ) );

        // Adding node to Model
        mpModel->AddNode( new_ho_node1 );

        // Adding DOFs
        if( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
        }
        else        
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
        }
	}

	//*************************************************************************************************
    //* - Creates new higher order nodes on edges ( 3th order ).
    //*************************************************************************************************
    void Modeler::Create_HONodes_OnEdge_3rd( int IdNode1, int IdNode2 )
    {
		Vector<unsigned int> edgeNodes( 4 );

        edgeNodes[ 0 ] = IdNode1;
        edgeNodes[ 1 ] = IdNode2;
                   
        edgeNodes[ 2 ] = ( mpModel->GetNodesArray() ).size() + 1;
        edgeNodes[ 3 ] = edgeNodes[ 2 ] + 1;

        mEdgeHONodes[ IdNode1*IdNode2 ].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1 - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2 - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        // Creating new high order node on edge number 1
        alpha = 1.0 / 3.0;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[ 2 ], hoX, hoY, hoZ ) );

		// Creating new high order node on edge number 2
        alpha = 2.0 / 3.0;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( edgeNodes[ 3 ], hoX, hoY, hoZ ) );

        // Adding nodes to Model
        mpModel->AddNode( new_ho_node1 );
        mpModel->AddNode( new_ho_node2 );

        // Adding DOFs
        if( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
            new_ho_node2->Add_cA_Dofs( mDofSet );
        }
        else     
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
            new_ho_node2->Add_cE_Dofs( mDofSet );
        }
	}

	//*************************************************************************************************
    //* - Creates a new higher order node on faces. 
    //*************************************************************************************************
    void Modeler::Create_HONodes_OnFace_3rd( int IdNode1, int IdNode2, int IdNode3 )
    {
		Vector<unsigned int> faceNodes( 4 );

        faceNodes[ 0 ] = IdNode1;
        faceNodes[ 1 ] = IdNode2;
        faceNodes[ 2 ] = IdNode3;

        faceNodes[ 3 ] = ( mpModel->GetNodesArray() ).size() + 1;

        mFaceHONodes[ IdNode1*IdNode2*IdNode3 ].push_back( faceNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNode1 - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2 - 1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNode3 - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();

        // Creating new high order node on face number 1
        double alpha = 1.0 / 3.0;
        double beta  = 1.0 / 3.0;

        double hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        double hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		double hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( faceNodes[ 3 ], hoX, hoY, hoZ ) );

        // Adding node to Model
        mpModel->AddNode( new_ho_node1 );

        // Adding DOFs
        if( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
        }
        else     
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
        }
	}

	//*************************************************************************************************
    //* - Creates new higher order nodes on edges ( 4th order ).
    //*************************************************************************************************
    void Modeler::Create_HONodes_OnEdge_4th( int IdNode1, int IdNode2 )
    {
        Vector<unsigned int> edgeNodes( 5 );

        edgeNodes[ 0 ] = IdNode1;
        edgeNodes[ 1 ] = IdNode2;

        edgeNodes[ 2 ] = ( mpModel->GetNodesArray() ).size() + 1;
        edgeNodes[ 3 ] = edgeNodes[ 2 ] + 1;
        edgeNodes[ 4 ] = edgeNodes[ 3 ] + 1;

        mEdgeHONodes[ IdNode1*IdNode2 ].push_back( edgeNodes );

        Node::Pointer pNode1 = mpModel->GetNode( IdNode1 - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2 - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();

		double hoX, hoY, hoZ;
        double alpha;

        // Creating new high order node on edge number 1
        alpha = 0.25;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( edgeNodes[ 2 ], hoX, hoY, hoZ ) );

		// Creating new high order node on edge number 2
        alpha = 0.50;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( edgeNodes[ 3 ], hoX, hoY, hoZ ) );

		// Creating new high order node on edge number 3
        alpha = 0.75;

        hoX = alpha * ( X2 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + Z1;

        Node::Pointer new_ho_node3( new Node( edgeNodes[ 4 ], hoX, hoY, hoZ ) );

        // Adding nodes to Model
        mpModel->AddNode( new_ho_node1 );
        mpModel->AddNode( new_ho_node2 );
        mpModel->AddNode( new_ho_node3 );

        // Adding DOFs
        if( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
            new_ho_node2->Add_cA_Dofs( mDofSet );
            new_ho_node3->Add_cA_Dofs( mDofSet );
        }
        else       
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
            new_ho_node2->Add_cE_Dofs( mDofSet );
            new_ho_node3->Add_cE_Dofs( mDofSet );
        }
    }

	//*************************************************************************************************
    //* - Creates new higher order nodes on faces ( 4th order ).
    //*************************************************************************************************
    void Modeler::Create_HONodes_OnFace_4th( int IdNode1, int IdNode2, int IdNode3 )
    {
        Vector<unsigned int> faceNodes( 6 );

        faceNodes[ 0 ] = IdNode1;
        faceNodes[ 1 ] = IdNode2;
        faceNodes[ 2 ] = IdNode3;

        faceNodes[ 3 ] = ( mpModel->GetNodesArray() ).size() + 1;
        faceNodes[ 4 ] = faceNodes[ 3 ] + 1;
        faceNodes[ 5 ] = faceNodes[ 4 ] + 1;

        mFaceHONodes[ IdNode1*IdNode2*IdNode3 ].push_back( faceNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNode1 - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNode2 - 1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNode3 - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();

		double hoX, hoY, hoZ;
        double alpha, beta;

        // Creating new high order node on face number 1
        alpha = 0.25;
        beta  = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( faceNodes[ 3 ], hoX, hoY, hoZ ) );

		// Creating new high order node on face number 2
        alpha = 0.50;
        beta  = 0.25;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node2( new Node( faceNodes[ 4 ], hoX, hoY, hoZ ) );

		// Creating new high order node on face number 3
        alpha = 0.25;
        beta  = 0.50;

        hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + X1;
        hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + Y1;
		hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + Z1;

        Node::Pointer new_ho_node3( new Node( faceNodes[ 5 ], hoX, hoY, hoZ ) );

        // Adding nodes to Model
        mpModel->AddNode( new_ho_node1 );
        mpModel->AddNode( new_ho_node2 );
        mpModel->AddNode( new_ho_node3 );

        // Adding DOFs
        if ( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
            new_ho_node2->Add_cA_Dofs( mDofSet );
            new_ho_node3->Add_cA_Dofs( mDofSet );
        }
        else   
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
            new_ho_node2->Add_cE_Dofs( mDofSet );
            new_ho_node3->Add_cE_Dofs( mDofSet );
        }
    }

	//*************************************************************************************************
    //* - Creates new higher order nodes inside volumes. 
    //*************************************************************************************************
	void Modeler::Create_HONodes_OnVolum( int* IdNodes )
    {
		Vector<unsigned int> volumNodes( 5 );

        volumNodes[ 0 ] = IdNodes[ 0 ];
        volumNodes[ 1 ] = IdNodes[ 1 ];
        volumNodes[ 2 ] = IdNodes[ 2 ];
		volumNodes[ 3 ] = IdNodes[ 3 ];

        volumNodes[ 4 ] = ( mpModel->GetNodesArray() ).size() + 1;

        mVolmHONodes[ IdNodes[ 0 ]*IdNodes[ 1 ]*IdNodes[ 2 ]*IdNodes[ 3 ] ].push_back( volumNodes );

		Node::Pointer pNode1 = mpModel->GetNode( IdNodes[ 0 ] - 1 );
        Node::Pointer pNode2 = mpModel->GetNode( IdNodes[ 1 ] - 1 );
        Node::Pointer pNode3 = mpModel->GetNode( IdNodes[ 2 ] - 1 );
		Node::Pointer pNode4 = mpModel->GetNode( IdNodes[ 3 ] - 1 );

		double X1 = pNode1->X(), Y1 = pNode1->Y(), Z1 = pNode1->Z();
		double X2 = pNode2->X(), Y2 = pNode2->Y(), Z2 = pNode2->Z();
		double X3 = pNode3->X(), Y3 = pNode3->Y(), Z3 = pNode3->Z();
		double X4 = pNode4->X(), Y4 = pNode4->Y(), Z4 = pNode4->Z();

		// Creating new high order node on volume number 1
        double alpha = 0.25;
        double beta  = 0.25;
		double gamma = 0.25;

        double hoX = alpha * ( X2 - X1 ) + beta * ( X3 - X1 ) + gamma * ( X4 - X1 ) + X1;
        double hoY = alpha * ( Y2 - Y1 ) + beta * ( Y3 - Y1 ) + gamma * ( Y4 - Y1 ) + Y1;
		double hoZ = alpha * ( Z2 - Z1 ) + beta * ( Z3 - Z1 ) + gamma * ( Z4 - Z1 ) + Z1;

        Node::Pointer new_ho_node1( new Node( volumNodes[ 4 ], hoX, hoY, hoZ ) );

        // Adding node to Model
        mpModel->AddNode( new_ho_node1 );

        // Adding DOFs
        if( mPotentials ) 
        {
            new_ho_node1->Add_cA_Dofs( mDofSet );
        }
        else  
        {
            new_ho_node1->Add_cE_Dofs( mDofSet );
        }
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 
	// - Building.
    // 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*****************************************************************************************************************
    //* - Checks if the problem set-up is consistent.
    //*****************************************************************************************************************
    void Modeler::Check_Consistency()
    {
        if( mElectrostaticMode )  
        {
            mElementOrder     = 1  ;
            mProblemFrequency = 0.0;

            mPotentials     = true ;
            mEdgeElementsOn = false;

            mComplxFreqMode = false;
            mFrequencySweep = false;

            mImportRobinBCs = false;
            mImportCurrents = false;
            mImportEleMatrx = false;

            if( mResultsOnGPs > 1 ) 
            {
                mResultsOnGPs = 1;
            }
        }
        else
        {
            if( ( mProblemFrequency == 0.0 ) && ( mComplxFreqMode == false ) ) 
            {
                Send_Error_Msg( "Problem frequency equal to 0.0", "Please, change problem frequency.", 0 );
            }        
        }

        if( mEdgeElementsOn )
        {
            if( mAxisymmetric ) 
            {
                Send_Error_Msg( "AY-symmetry mode not defined for edge elements", "Please, change element type or symmetry mode.", 0 );
            }

            if( mSetAllFzoToZero ) 
            {
                Send_Error_Msg( "XY-symmetry mode not defined for edge elements", "Please, change element type or symmetry mode.", 0 );
            }

            if( mSetAllFxyToZero ) 
            {
                Send_Error_Msg( "ZZ-symmetry mode not defined for edge elements", "Please, change element type or symmetry mode.", 0 );
            }

            mAVContinuity = true; 
        }

        if( ( mElementOrder != 0 ) || ( mResultsOnNodes == false ) ) 
        {
            mLL2PSmoothing = false;
        }  

        if( mComplxFreqMode )
        {
            if( std::abs( mComplexFrequency ) == 0.0 ) 
            {
                Send_Error_Msg( "Complex frequency equal to 0.0", "Please, change complex frequency.", 0 );
            }        

            if( mColdPlasmaMode ) 
            {
                Send_Error_Msg( "Complex frequency not implemented for cold plasma mode", "Please, change problem or frequency mode.", 0 );
            }   

            if( std::imag( mComplexFrequency ) != 0.0 )
            {
                double pi2 = 6.283185307179586476925286766559;

                mProblemFrequency = pi2 * std::imag( mComplexFrequency );
            }
            else
            {
                mProblemFrequency = std::real( mComplexFrequency );
            }

            mPotentials     = true;
            mFrequencySweep = false;
        }

        Print_Set_Up_Info();

        std::cout << "Reading nodes and conditions..." << std::endl; 
    }

    //*****************************************************************************************************************
    //* - Prints on screen set-up info.
    //*****************************************************************************************************************
    void Modeler::Print_Set_Up_Info()
    {
        std::cout << "Problem  type: ";

        if     ( mFullWaveMode      ) std::cout << "Full wave"     << std::endl;
        else if( mElectrostaticMode ) std::cout << "Electrostatic" << std::endl; 
        else if( mColdPlasmaMode    ) std::cout << "Cold plasma"   << std::endl; 

        std::cout << "Geometry type: ";

        if     ( mSetAllFzoToZero ) std::cout << "XY" << std::endl;
        else if( mSetAllFxyToZero ) std::cout << "ZZ" << std::endl; 
        else if( mAxisymmetric    ) std::cout << "AY" << std::endl;
        else                        std::cout << "3D" << std::endl;    

        std::cout << "Element  type: ";

        if     ( mElementOrder ==  0 ) std::cout << "LL2P 3sb" << std::endl;
        else if( mElementOrder ==  1 ) std::cout << "RME 1st"  << std::endl; 
        else if( mElementOrder ==  2 ) std::cout << "RME 2nd"  << std::endl; 
        else if( mElementOrder == 11 ) std::cout << "EDGE 1st" << std::endl; 
        else if( mElementOrder == 12 ) std::cout << "EDGE 2nd" << std::endl; 

        std::cout << "AV potentials: ";

        if  ( mPotentials ) std::cout << "On"  << std::endl; 
        else                std::cout << "Off" << std::endl;  

        std::cout << "Stabilization: ";

        if     ( ( mElementOrder ==  0 ) && ( mLL2P_hk_factor !=  0.0 ) ) std::cout << "On"  << std::endl; 
        else if( ( mElementOrder ==  1 ) && ( mRMED_stabilize == true ) ) std::cout << "On"  << std::endl;  
        else if( ( mElementOrder ==  2 ) && ( mRMED_stabilize == true ) ) std::cout << "On"  << std::endl;  
        else if( ( mElementOrder == 11 ) && ( mEDGE_stabilize == true ) ) std::cout << "On"  << std::endl;  
        else if( ( mElementOrder == 12 ) && ( mEDGE_stabilize == true ) ) std::cout << "On"  << std::endl;  
        else                                                              std::cout << "Off" << std::endl; 

        double pi2 = 6.283185307179586476925286766559;

        if( mFrequencySweep ) 
        {
            std::cout << "Ini angl freq: " << mProblemFrequency / pi2 << " Hz" << std::endl;   
        }
        else if( mComplxFreqMode ) 
        {
            std::cout << "Complex  freq: " << std::real( mComplexFrequency ) << " + j * " << std::imag( mComplexFrequency ) << std::endl; 
        }
        else
        {
            std::cout << "Angular  freq: " << mProblemFrequency / pi2 << " Hz" << std::endl; ;
        }

        std::cout << std::endl; 
    }

    //*****************************************************************************************************************
    //* - Assembles volume elements in electro-static mode.
    //*****************************************************************************************************************
    void Modeler::Generate_Volume_Element_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId == 0 throw error, material not assigned
        if( PropertiesId == 0 ) throw 0; 

        // Get local nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get material properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        // Create element
        Element::Pointer pElement;

        pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector
        Vector<double> EleResVector;

        pElement->GetResidualVector_Dirichlet( mFix_Static_Voltage, EleStiffMatrix, EleResVector );

        // Assemble elemental stiffness matrix to symmetric global matrix
        Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );

        // Assemble elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Clean objects
        EleIdsVector.FreeData();
        EleResVector.FreeData();
        EleStiffMatrix.FreeData();
    }

    //*****************************************************************************************************************
    //* - Sets element parameters ( singularity weight, normal type, LL2P factor ).
    //*****************************************************************************************************************
    void Modeler::Set_Volume_Element_Parameters( int* NodesId, Element::Pointer pElement )
    {
		// Checking if the element is inside a sigularity layer
		bool IsUgLayer = false;

		for( int i=0; i<4; i++ )
		{
			if( ( mSingular   .find( NodesId[i] ) != mSingular   .end() ) || ( mNSingular2L.find( NodesId[i] ) != mNSingular2L.end() ) ||
			    ( mNSingular3L.find( NodesId[i] ) != mNSingular3L.end() ) || ( mNSingular4L.find( NodesId[i] ) != mNSingular4L.end() ) ||
			    ( mNSingular5L.find( NodesId[i] ) != mNSingular5L.end() ) || ( mNSingular6L.find( NodesId[i] ) != mNSingular6L.end() )  )
			{
				IsUgLayer = true;
			}
		}

        // Setting regularization weight 
		if( IsUgLayer || !mRMED_stabilize ) 
        {
            pElement->SetPeso( 0.0 );
        }
		else       
        {
            pElement->SetPeso( 1.0 );
        }

        // Setting parameters for LL2P bubble elements
        if( mElementOrder == 0 ) 
        {
            pElement->SetNormalType( mType_Of_BC_Normal );

            pElement->SetFactor( mLL2P_hk_factor );
        }

        // Setting plasma model
        if( mColdPlasmaMode )        
        {
            pElement->SetPlasmaModel( mpColdPlasma );
        }

        // Activates stabilizer in edge elements
        if( mEdgeElementsOn )
        {
            pElement->SetStabilizer( mEDGE_stabilize );
        }

        // Activates potentials on element
        pElement->SetPotentials( mPotentials );

        // Set complex frequency
        if( mComplxFreqMode )
        {
            pElement->SetFrequency( mComplexFrequency );
        }
    }

    //**********************************************************************************************************************************************************
    //* - Applies boundary conditions to element stiffness matrix.
    //**********************************************************************************************************************************************************
    void Modeler::Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        // PBC conditions ( Trans(T)*A*T )
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_RightLeft_TKT( HONodesId, EleIdVector, EleStiffMatrix );
        PBC_FrontBack_TKT( HONodesId, EleIdVector, EleStiffMatrix );

        // Contact discontinuities ( Trans(T)*A*T )
        Apply_Contact_Matrix_TKT( HONodesId, EleIdVector, EleStiffMatrix );
       
        // PEC, PMC conditions ( Trans(T)*A*T )
        Apply_Rotation_To_Matrix( HONodesId, EleStiffMatrix );

        // Axisymmetric mode ( Trans(T)*A*T ) 
		Apply_Axisymmetric_TKT( HONodesId, EleStiffMatrix );

        // Cold plasma RLP base ( Inv(T)*A*T )
        if( mColdPlasmaMode ) 
        {
            ColdPlasma_TKT( HONodesId, EleStiffMatrix );
        }
    }

    //**********************************************************************************************************************************************************
    //* - Applies boundary conditions to element force vector.
    //**********************************************************************************************************************************************************
    void Modeler::Apply_Element_Boundary_Conditions( Vector<int>& HONodesId, Vector<int>& EleIdVector, Vector< std::complex<double> >& EleResVector )
    {
	    // PBC conditions ( Inv(T)*B ) 
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );
	    PBC_RightLeft_TR( HONodesId, EleIdVector, EleResVector );
		PBC_FrontBack_TR( HONodesId, EleIdVector, EleResVector );

		// Contact discontinuities ( Inv(T)*B )
        Apply_Contact_Vector_TR( HONodesId, EleIdVector, EleResVector );

        // PEC, PMC conditions ( Inv(T)*B ) 
        Apply_Rotation_To_Vector( HONodesId, EleResVector );

		// Axisymmetric mode ( Inv(T)*B )
        Apply_Axisymmetric_TR( HONodesId, EleResVector );

        // Cold plasma RLP base ( Inv(T)*B )
        if( mColdPlasmaMode ) 
        {
            ColdPlasma_TR( HONodesId, EleResVector );
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental stiffness matrix (real) to global symmetric matrix.
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Symmetric( Vector<int>& EleIdVector, Matrix<double>& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        std::complex<double> cKij;

		for ( int i_local = 0; i_local < eleSize; i_local++ )
		{
			int i_global = EleIdVector[ i_local ];

			if ( i_global < mSystemSize )
			{
				for ( int j_local = 0; j_local < eleSize; j_local++ )
				{
					int j_global = EleIdVector[ j_local ];

					if ( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[ i_local ][ j_local ];

                        if ( abs( cKij ) > 0.0 ) A_matrix( i_global, j_global ) += cKij;
                    }
				}
			}
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental stiffness matrix (complex) to global symmetric matrix.
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Symmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;

		for( int i_local = 0; i_local < eleSize; i_local++ )
		{
			int i_global = EleIdVector[ i_local ];

			if( i_global < mSystemSize )
			{
				for( int j_local = 0; j_local < eleSize; j_local++ )
				{
					int j_global = EleIdVector[ j_local ];

					if( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[ i_local ][ j_local ];

                        if( abs( cKij ) > Kij_Tol ) A_matrix( i_global, j_global ) += cKij;
                    }
				}
			}
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental stiffness matrix (complex) to global non-symmetric matrix.
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_NonSymmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;

        for( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdVector[i_local];
        
            if( i_global < mSystemSize )
            {
                for( int j_local = 0; j_local < eleSize; j_local++ )
                {
                    int j_global = EleIdVector[j_local];
        
                    if( ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
        
                        if( abs( cKij ) > Kij_Tol ) A_matrix( i_global, j_global ) += cKij;
                    }
                }
            }
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental stiffness matrix (complex) to auxiliar symmetric global matrix.
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Aux_Symmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;        
        
        for( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdVector[i_local];
        
            if( i_global < mSystemSize )
            {
                for( int j_local = 0; j_local < eleSize; j_local++ )
                {
                    int j_global = EleIdVector[j_local];
        
                    if( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];
        
                        if( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global, j_global ) += cKij;
                    }
                }
            }
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental stiffness matrix (complex) to auxiliar non-symmetric global matrix.
    //*****************************************************************************************************************************************
    void Modeler::Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( Vector<int>& EleIdVector, Matrix< std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();

        double Kij_Tol = 0.00;

        if( mColdPlasmaMode ) Kij_Tol = mpColdPlasma->Get_Kij_Tol();

        std::complex<double> cKij;        

        for( int i_local = 0; i_local < eleSize; i_local++ )
        {
            int i_global = EleIdVector[i_local];

            if( i_global < mSystemSize )
            {
                for( int j_local = 0; j_local < eleSize; j_local++ )
                {
                    int j_global = EleIdVector[j_local];

                    if( ( j_global < mSystemSize ) )
                    {
                        cKij = EleStiffMatrix[i_local][j_local];

                        if( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global, j_global ) += cKij;
                    }
                }
            }
        }
    }

    //*****************************************************************************************************************************************
    //* - Assembles elemental residual vector (complex) to global residual vector.
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
    //* - Assembles elemental residual vector (real) to global residual vector.
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
    //* - Assembles volume elements in cold plasma mode. 
    //*************************************************************************************************************************
    void Modeler::Generate_Volume_Element_ColdPlasma( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId == 0 throw error, material not assigned
        if( PropertiesId == 0 ) throw 0; 

        // Get volume element nodes
        Vector<int> HONodesId;

		Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Plasma element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_ColdPlasma( pNodes, pProperties ) );
		else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_ColdPlasma( pNodes, pProperties ) );

        // Set element parameters
        Set_Volume_Element_Parameters( NodesId, pElement );

        // Get DOFs Id vector
        Vector<int> EleIdsVector; 

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix; 

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector 
        Vector< std::complex<double> > EleResVector;
            
        pElement->GetResidualVector_Dirichlet( mFix_FullWv_Voltage, EleStiffMatrix, EleResVector );

        // Apply boundary conditions to stiffness matrix
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Apply boundary conditions to residual vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assemble elemental stiffness matrix to global matrix 
        if     ( mpColdPlasma->Is_Full_Matrix    () ) Assemble_ElementMatrix_to_Global_NonSymmetric    ( EleIdsVector, EleStiffMatrix );
        else if(    !pElement->GetIsIHL          () ) Assemble_ElementMatrix_to_Global_Symmetric       ( EleIdsVector, EleStiffMatrix );
        else if( mpColdPlasma->Is_HermSymm_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_Symmetric   ( EleIdsVector, EleStiffMatrix );           
        else if( mpColdPlasma->Is_HermFull_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( EleIdsVector, EleStiffMatrix );

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );
     
        // Clean objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*****************************************************************************************************************
    //* - Assembles volume elements in Full-Wave mode. 
    //*****************************************************************************************************************
    void Modeler::Generate_Volume_Element_FullWave( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId == 0 throw error, material not assigned
        if( PropertiesId == 0 ) throw 0; 

        // Get volume element nodes
        Vector<int> HONodesId;

		Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId ); 
        
        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Volume element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
		else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

        // Set element parameters
        Set_Volume_Element_Parameters( NodesId, pElement );

        // Get DOFs Id vector
        Vector<int> EleIdsVector; 

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix; 

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector 
        Vector< std::complex<double> > EleResVector;
            
        pElement->GetResidualVector_Dirichlet( mFix_FullWv_Voltage, EleStiffMatrix, EleResVector );

        // Apply boundary conditions to stiffness matrix
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Apply boundary conditions to residual vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assemble elemental stiffness matrix to symmetric global matrix
        Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Clean objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembles high order source elements.
    //*************************************************************************************************************************************
    void Modeler::Generate_Source_Element( int* NodesId, unsigned int PropertiesId )
    {
        // Get volume element nodes
        Vector<int> HONodesId;

		Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Volumetric source element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new JSource_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new JSource_2eg_FullWave( pNodes, pProperties ) );
	    
        // Activates potentials on element
        pElement->SetPotentials( mPotentials );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental residual vector
        Vector< std::complex<double> > EleResVector; 

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
    //* - Assembles boundary plasma element matrix on auxiliar global matrix.  
    //*************************************************************************************************************************************
    void Modeler::Assemble_BC_Plasma_Element_On_AuxMatrix( Vector<int>& EleIdVector, Matrix<std::complex<double> >& EleStiffMatrix )
    {
        int eleSize = EleIdVector.size();
         
        std::complex<double> cKij;
        
        double Kij_Tol = mpColdPlasma->Get_Kij_Tol();
        
        if( mpColdPlasma->Is_HermSymm_Matrix() )
        {
            for( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[ i_local ];
        
                if( i_global < mSystemSize )
                {
                    for( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[ j_local ];
        
                        if( ( i_global <= j_global ) && ( j_global < mSystemSize ) )
                        {
                            cKij = EleStiffMatrix[ i_local ][ j_local ];

                            if( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else if( mpColdPlasma->Is_HermFull_Matrix() )
        {
            for( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[ i_local ];
        
                if( i_global < mSystemSize )
                {
                    for( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[ j_local ];
        
                        if( j_global < mSystemSize )
                        {
                            cKij = EleStiffMatrix[ i_local ][ j_local ];

                            if( abs( cKij ) > Kij_Tol ) A_matrix_aux( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
        else
        {
            for( int i_local = 0; i_local < eleSize; i_local++ )
            {
                int i_global = EleIdVector[ i_local ];
        
                if( i_global < mSystemSize )
                {
                    for( int j_local = 0; j_local < eleSize; j_local++ )
                    {
                        int j_global = EleIdVector[ j_local ];
        
                        if( j_global < mSystemSize )
                        {
                            cKij = EleStiffMatrix[ i_local ][ j_local ];

                            if( abs( cKij ) > Kij_Tol ) A_matrix( i_global,j_global ) += cKij;
                        }
                    }
                }
            }
        }
    }

    //*************************************************************************************************************************************
    //* - Assembles cold plasma far field elements.
    //*************************************************************************************************************************************
    void Modeler::Generate_FarField_Element_ColdPlasma( int* NodesId, unsigned int PropertiesId )
    {
        // Get surface element nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

		Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Element properties
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        // Far field (cold plasma) element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new FarField_3sb_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new FarField_1st_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new FarField_2nd_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new FarField_1eg_ColdPlasma( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new FarField_2eg_ColdPlasma( pNodes, pProperties ) );

        // Set element parameters
        pElement->SetPlasmaModel( mpColdPlasma      );
        pElement->SetPotentials ( mPotentials       );
        pElement->SetFrequency  ( mProblemFrequency );

        if( !mRMED_stabilize ) pElement->SetPeso( 0.0 );

        // Get DOFs Id vector
		Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get element stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Apply boundary conditions
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Assembling to global matrix 
        if( mColdPlasmaMode ) 
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix );    
        }
        else    
        {
            Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
        }
 
        // Cleaning elemental matrix
        EleIdsVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembles full wave far field elements. 
    //*************************************************************************************************************************************
    void Modeler::Generate_FarField_Element_FullWave( int* NodesId )
    {
        // Get surface element nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

		Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Far field element definition
        Element::Pointer pElement;
 
        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new FarField_3sb_FullWave( pNodes ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new FarField_1st_FullWave( pNodes ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new FarField_2nd_FullWave( pNodes ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new FarField_1eg_FullWave( pNodes ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new FarField_2eg_FullWave( pNodes ) );

        // Set element parameters
        pElement->SetPotentials( mPotentials );
		pElement->SetFrequency( mProblemFrequency );

        // Get DOFs Id vector
		Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get element stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Apply boundary conditions
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Assembling to global matrix 
        if( mColdPlasmaMode ) 
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix ); 
        }
        else  
        {
            Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
        }

        // Cleaning objects
        EleIdsVector  .FreeData();
        EleStiffMatrix.FreeData();
	}

	//*************************************************************************************************************************************
    //* - Assembles generic Robin elements.
    //*************************************************************************************************************************************
    void Modeler::Generate_GenericRobin_Element_FullWave( int* NodesId, unsigned int PropertiesId )
    {
        // Get surface element nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

		Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Element properties
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        // Generic Robin element definition
        Element::Pointer pElement;
		
        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new GenericRobin_3sb_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new GenericRobin_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new GenericRobin_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new GenericRobin_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new GenericRobin_2eg_FullWave( pNodes, pProperties ) );

        // Set element parameters
        pElement->SetPotentials( mPotentials ); 
        pElement->SetFrequency( mProblemFrequency );

        // Get DOFs Id vector
		Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get element stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector 
        Vector< std::complex<double> > EleResVector;

        // Calculate residual vector for a quasi-static flux
        pElement->GetResidualVector_Flux( mFix_FullWv_Current, EleResVector );

        // Calculate residual vector for plane waves and Gaussian beams
        pElement->GetResidualVector_Flux( mPlaneWavesParameters, EleResVector );

        // Apply boundary conditions to stiffness matrix
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        // Apply boundary conditions to residual vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assembling to global matrix 
        if( mColdPlasmaMode ) 
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix );    
        }
        else    
        {
            Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
        }

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Clean objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
	}

    //*************************************************************************************************************************************
    //* - Assembles generic Robin elements (Electrostatic mode).
    //*************************************************************************************************************************************
    void Modeler::Generate_GenericRobin_Element_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // Getting local nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

        Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Getting material properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        // Creating element
        Element::Pointer pElement;

        pElement = Element::Pointer( new GenericRobin_1st_Electrostatic( pNodes, pProperties ) );

        // Getting DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Getting elemental stiffness matrix
        Matrix<double> EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Including Dirichlet BC in the residual vector
        Vector<double> EleResVector;

        pElement->GetResidualVector( EleResVector );

        // Assembling to global matrix
        Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );

        // Assemble elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Cleaning objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembles rectangular waveguide port elements ( TE10 mode ).
    //*************************************************************************************************************************************
    void Modeler::Generate_RectPort_TE10_Element( int* NodesId, unsigned int PropertiesId )
    {
		// Get surface element nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

		Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Element properties
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Rectangular waveguide TE10 mode element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new RWPortTE10_3sb_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  1 ) pElement = Element::Pointer( new RWPortTE10_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new RWPortTE10_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new RWPortTE10_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new RWPortTE10_2eg_FullWave( pNodes, pProperties ) );

        // Set element parameters
        pElement->SetPotentials( mPotentials );

        // Get DOFs Id vector
		Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get element stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector
		Vector< std::complex<double> > EleResVector  ;

        pElement->GetResidualVector( EleResVector );

        // Apply boundary conditions to matrix and force vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assembling to global matrix 
        if( mColdPlasmaMode ) 
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix ); 
        }
        else  
        {
            Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
        }

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Cleaning objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	//*************************************************************************************************************************************
    //* - Assembles coaxial waveguide port elements ( TEM mode ).
    //*************************************************************************************************************************************
    void Modeler::Generate_CoaxPort_TEM_Element( int* NodesId, unsigned int PropertiesId )
    {
		// Get surface element nodes
        Vector<int> HONodesId;

        Vector<Node::Pointer> pNodes;

		Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

        // Element properties
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Coaxial waveguide TEM mode element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new CoaxialPortTEM_3sb_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new CoaxialPortTEM_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new CoaxialPortTEM_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new CoaxialPortTEM_2eg_FullWave( pNodes, pProperties ) );

        // Set element parameters
        pElement->SetPotentials( mPotentials );

        // Get DOFs Id vector
		Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get element stiffness matrix
        Matrix< std::complex<double> > EleStiffMatrix;

        pElement->GetStiffnessMatrix( EleStiffMatrix );

        // Get residual vector
		Vector< std::complex<double> > EleResVector  ;

        pElement->GetResidualVector( EleResVector );

        // Apply boundary conditions to matrix and force vector
        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

        Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

        // Assembling to global matrix 
        if( mColdPlasmaMode ) 
        {
            Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix );    
        }
        else    
        {
            Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
        }

        // Assembling elemental residual vector to global residual vector
        Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

        // Cleaning objects
        EleIdsVector  .FreeData();
        EleResVector  .FreeData();
        EleStiffMatrix.FreeData();
    }

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 
	// - Matrix rotations.
    // 
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*****************************************************************************************************************
    // - Trans(T) * K * T 
    //*****************************************************************************************************************
    void Modeler::TransT_K_T( Matrix< std::complex<double> >& K, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex )
    {
        int NumDofs = 3 * NodesIds.size();

        if( mPotentials ) 
        {
            for( int i=0; i<NodesIds.size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( NodesIds[ i ]-1 );
                
                if( pNode->pDof( cVs ) != NULL ) 
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
        for( int i=0; i<NumDofs; i++ )
        {
            Dof_X = K[ i ][ NodeIndex          ];
            Dof_Y = K[ i ][ NodeIndex +   Step ];
            Dof_Z = K[ i ][ NodeIndex + 2*Step ];
        
            K[ i ][ NodeIndex          ] = T[ 0 ][ 0 ] * Dof_X + T[ 0 ][ 1 ] * Dof_Y + T[ 0 ][ 2 ] * Dof_Z;
            K[ i ][ NodeIndex +   Step ] = T[ 1 ][ 0 ] * Dof_X + T[ 1 ][ 1 ] * Dof_Y + T[ 1 ][ 2 ] * Dof_Z;
            K[ i ][ NodeIndex + 2*Step ] = T[ 2 ][ 0 ] * Dof_X + T[ 2 ][ 1 ] * Dof_Y + T[ 2 ][ 2 ] * Dof_Z;
        }
        
        // Transpose(T) * K
        for( int j=0; j<NumDofs; j++ )
        {
            Dof_X = K[ NodeIndex          ][ j ];
            Dof_Y = K[ NodeIndex +   Step ][ j ];
            Dof_Z = K[ NodeIndex + 2*Step ][ j ];
        
        	K[ NodeIndex          ][ j ] = T[ 0 ][ 0 ] * Dof_X + T[ 0 ][ 1 ] * Dof_Y + T[ 0 ][ 2 ] * Dof_Z;
            K[ NodeIndex +   Step ][ j ] = T[ 1 ][ 0 ] * Dof_X + T[ 1 ][ 1 ] * Dof_Y + T[ 1 ][ 2 ] * Dof_Z;
            K[ NodeIndex + 2*Step ][ j ] = T[ 2 ][ 0 ] * Dof_X + T[ 2 ][ 1 ] * Dof_Y + T[ 2 ][ 2 ] * Dof_Z;
        }
    }

    //*****************************************************************************************************************
    // - Trans(T) * B 
    //*****************************************************************************************************************
    void Modeler::TransT_B( Vector< std::complex<double> >& B, Matrix< std::complex<double> >& T, Vector<int>& NodesIds, int NodeIndex )
    {
        int Step = NodesIds.size();
   
        std::complex<double> Dof_X;
        std::complex<double> Dof_Y;
        std::complex<double> Dof_Z;
        
	    Dof_X = B[ NodeIndex          ];
        Dof_Y = B[ NodeIndex +   Step ];
        Dof_Z = B[ NodeIndex + 2*Step ];
        
        B[ NodeIndex          ] = T[ 0 ][ 0 ] * Dof_X + T[ 0 ][ 1 ] * Dof_Y + T[ 0 ][ 2 ] * Dof_Z;
        B[ NodeIndex +   Step ] = T[ 1 ][ 0 ] * Dof_X + T[ 1 ][ 1 ] * Dof_Y + T[ 1 ][ 2 ] * Dof_Z;
        B[ NodeIndex + 2*Step ] = T[ 2 ][ 0 ] * Dof_X + T[ 2 ][ 1 ] * Dof_Y + T[ 2 ][ 2 ] * Dof_Z;
    }

	//*****************************************************************************************************************
    // - Rotates the stiffness matrix to the n, t, b coordinate system ( for PEC/PMC boundary conditions ).
    //*****************************************************************************************************************
    void Modeler::Apply_Rotation_To_Matrix( Vector<int>& HONodesId, Matrix<std::complex<double> >& StiffMatrix )
    {
		if( ( mNormals.size() == 0 ) || ( mSetAllFxyToZero ) || ( mEdgeElementsOn ) )
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            int NodeId = HONodesId[ local_i ];
            
            if( ( !Is_DisctAVNode( NodeId ) ) && ( mNormals.find( NodeId ) != mNormals.end() ) )
            {
                Matrix< std::complex<double> > RotMatrix;

                Vector<double> n( mNormals[ NodeId ] );
        	    
                Get_Rotation_Matrix( RotMatrix, n );
                
                TransT_K_T( StiffMatrix, RotMatrix, HONodesId, local_i );
            }
        }
    }

	//*****************************************************************************************************************
    // - Calculates rotation matrix to n, t, b coordinate system.
    //*****************************************************************************************************************
    void Modeler::Get_Rotation_Matrix( Matrix< std::complex<double> >& RotMatrix, Vector<double>& n )
    {
        RotMatrix.Resize( 3, 3 );
        
        Vector<double> t( 3 );
        Vector<double> b( 3 );

        TangencialCoordinates( n, t, b );

        RotMatrix[ 0 ][ 0 ] = n[ 0 ]; RotMatrix[ 0 ][ 1 ] = n[ 1 ]; RotMatrix[ 0 ][ 2 ] = n[ 2 ];
        RotMatrix[ 1 ][ 0 ] = t[ 0 ]; RotMatrix[ 1 ][ 1 ] = t[ 1 ]; RotMatrix[ 1 ][ 2 ] = t[ 2 ];
        RotMatrix[ 2 ][ 0 ] = b[ 0 ]; RotMatrix[ 2 ][ 1 ] = b[ 1 ]; RotMatrix[ 2 ][ 2 ] = b[ 2 ];
    }

	//**********************************************************************************************************************
    //* - Inv(T) * K * T to rotate stiffness matrix to the axisymmetric coordinate system.
    //**********************************************************************************************************************
    void Modeler::Apply_Axisymmetric_TKT( Vector<int>& HONodesId, Matrix< std::complex<double> >& StiffMatrix )
    {
        if( mAxisymmetric == false ) 
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); ++local_i )
        {
			if( ( mNormals.find ( HONodesId[ local_i ] ) == mNormals.end() ) && ( DistanceToAxis( HONodesId[ local_i ] ) > mGeoTolerance ) )
            {
				Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( HONodesId[ local_i ], n, t, b );

                Matrix< std::complex<double> > RotMatrix( 3, 3 );
        
                RotMatrix[ 0 ][ 0 ] = n[ 0 ]; RotMatrix[ 0 ][ 1 ] = n[ 1 ]; RotMatrix[ 0 ][ 2 ] = n[ 2 ];
                RotMatrix[ 1 ][ 0 ] = t[ 0 ]; RotMatrix[ 1 ][ 1 ] = t[ 1 ]; RotMatrix[ 1 ][ 2 ] = t[ 2 ];
                RotMatrix[ 2 ][ 0 ] = b[ 0 ]; RotMatrix[ 2 ][ 1 ] = b[ 1 ]; RotMatrix[ 2 ][ 2 ] = b[ 2 ];

                TransT_K_T( StiffMatrix, RotMatrix, HONodesId, local_i );
            }
        }
    }

    //**********************************************************************************************************************
    //* - Inv(T)*StiffMatrix*T to rotate matrix to the diagonal permittivity tensor base (RLP)
    //* - Being:
    //* - K_RPL = [U][T] * K_CC * [Trans(T)][Herm(U)]
    //* - E_SDP = [T] * E_CC  //  E_RLP = [U] * E_SDP
    //**********************************************************************************************************************
    void Modeler::ColdPlasma_TKT( Vector<int>& HONodesId, Matrix< std::complex<double> >& StiffMatrix )
    {
        if( mEdgeElementsOn || mpColdPlasma->Is_Epar_Tol_Off() ) 
        {
            return;
        }

        int local_i = 0;

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            if( mNormals.find( *hoit ) == mNormals.end() )
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
                Matrix< std::complex<double> > T(3,3);

                Rotation_Matrix_SDP_T_CC( T, *hoit );

                // inv[R] = [U][T]
                Matrix< std::complex<double> > invR(3,3);

                invR = U * T;

                // [R] = [Trans(T)][Herm(U)]
                Matrix< std::complex<double> > R(3,3);

                for( int i=0; i<3; i++ )
                {
                    for( int j=0; j<3; j++ )
                    {
                        U[i][j] = std::conj(U[i][j]);
                    }
                }

                R = T.Transpose() * U.Transpose();

                // Inv(R) * StiffMatrix * R  
                int step    = HONodesId.size();
                int step2   = 2*step;
                int dofSize = 3*step;

                std::complex<double> dof_x, dof_y, dof_z;

                // StiffMatrix * R
                for( int i=0; i<dofSize; i++ )
                {
                    dof_x = StiffMatrix[i][local_i        ];
                    dof_y = StiffMatrix[i][local_i + step ];
                    dof_z = StiffMatrix[i][local_i + step2];

                    StiffMatrix[i][local_i        ] = R[0][0]*dof_x + R[1][0]*dof_y + R[2][0]*dof_z;
                    StiffMatrix[i][local_i + step ] = R[0][1]*dof_x + R[1][1]*dof_y + R[2][1]*dof_z;
                    StiffMatrix[i][local_i + step2] = R[0][2]*dof_x + R[1][2]*dof_y + R[2][2]*dof_z;
                }

                // Inv(R) * StiffMatrix
                for( int i=0; i<dofSize; i++ )
                {
                    dof_x = StiffMatrix[local_i        ][i];
                    dof_y = StiffMatrix[local_i + step ][i];
                    dof_z = StiffMatrix[local_i + step2][i];

                    StiffMatrix[local_i        ][i] = invR[0][0]*dof_x + invR[0][1]*dof_y + invR[0][2]*dof_z;
                    StiffMatrix[local_i + step ][i] = invR[1][0]*dof_x + invR[1][1]*dof_y + invR[1][2]*dof_z;
                    StiffMatrix[local_i + step2][i] = invR[2][0]*dof_x + invR[2][1]*dof_y + invR[2][2]*dof_z;
                }
            }

            local_i++;
        }
    }

	//**********************************************************************************************************************
    // - Rotates the force vector to the n, t, b coordinate system ( for PEC / PMC boundary conditions ).
    //**********************************************************************************************************************
    void Modeler::Apply_Rotation_To_Vector( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector )
    {
        if( ( mNormals.size() == 0 ) || ( mSetAllFxyToZero ) || ( mEdgeElementsOn ) ) 
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            int NodeId = HONodesId[ local_i ];
            
            if( ( !Is_DisctAVNode( NodeId ) ) && ( mNormals.find( NodeId ) != mNormals.end() ) )
            {
                Matrix< std::complex<double> > RotMatrix;

                Vector<double> n( mNormals[ NodeId ] );
        	    
                Get_Rotation_Matrix( RotMatrix, n );

                TransT_B( ResVector, RotMatrix, HONodesId, local_i );
            }
        }
    }

	//**********************************************************************************************************************
    //* - Inv(T) * B for axisymmetric problems.
    //**********************************************************************************************************************
    void Modeler::Apply_Axisymmetric_TR( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector )
    {
        if( mAxisymmetric == false ) 
        {
            return;
        }
        
        for( int local_i = 0; local_i < HONodesId.size(); ++local_i )
        {
			if( ( mNormals.find ( HONodesId[ local_i ] ) == mNormals.end() ) && 
                ( DistanceToAxis( HONodesId[ local_i ] )  > mGeoTolerance  )  )
            {
				Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( HONodesId[ local_i ], n, t, b );

                int Step = HONodesId.size();

                std::complex<double> Dof_X = ResVector[ local_i          ];
                std::complex<double> Dof_Y = ResVector[ local_i +   Step ];
                std::complex<double> Dof_Z = ResVector[ local_i + 2*Step ];

                ResVector[ local_i          ] = n[0]*Dof_X + n[1]*Dof_Y + n[2]*Dof_Z;
                ResVector[ local_i +   Step ] = t[0]*Dof_X + t[1]*Dof_Y + t[2]*Dof_Z;
                ResVector[ local_i + 2*Step ] = b[0]*Dof_X + b[1]*Dof_Y + b[2]*Dof_Z;
            }
        }
    }

    //**********************************************************************************************************************
    //* - Inv(T)*b for cold plasma mode.
    //**********************************************************************************************************************
    void Modeler::ColdPlasma_TR( Vector<int>& HONodesId, Vector< std::complex<double> >& ResVector )
    {
        if( mEdgeElementsOn || mpColdPlasma->Is_Epar_Tol_Off() ) 
        {
            return;
        }

        int local_i = 0;

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            if( mNormals.find( *hoit ) == mNormals.end() )
            {
                // Complex constants
                std::complex<double> cZero    ( 0.0        , 0.0         );
                std::complex<double> cUnit    ( 0.0        , 1.0         );
                std::complex<double> rUnit    ( 1.0        , 0.0         );
                std::complex<double> rInvSqrt2( 1.0/sqrt(2), 0.0         );
                std::complex<double> cInvSqrt2( 0.0        , 1.0/sqrt(2) );

                // E_RLP = [U] * E_SDP
                Matrix< std::complex<double> > U(3,3);

                U[0][0] = rInvSqrt2; U[0][1] =  cInvSqrt2; U[0][2] = cZero;
                U[1][0] = rInvSqrt2; U[1][1] = -cInvSqrt2; U[1][2] = cZero;
                U[2][0] = cZero    ; U[2][1] =  cZero    ; U[2][2] = rUnit;

                // E_SDP = [T] * E_CC
                Matrix< std::complex<double> > T(3,3);

                Rotation_Matrix_SDP_T_CC( T, *hoit );

                // inv[R] = [U][T]
                Matrix< std::complex<double> > invR(3,3);

                invR = U * T;

                // Inv(T)*b  
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

	//**********************************************************************************************************************
    //* - Local system of coordinates in PEC surface for n.
    //**********************************************************************************************************************
    void Modeler::TangencialCoordinates( Vector<double>& n, Vector<double>& t, Vector<double>& b )
    {
        double tnorm = sqrt( n[0]*n[0] + n[1]*n[1] );

        if( tnorm > 0.00 )
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

    //**********************************************************************************************************************
    //* -  Calculates the rotation matrix from CC to SDP ( E_SDP = [T] * E_CC ).
    //**********************************************************************************************************************
    void Modeler::Rotation_Matrix_SDP_T_CC( Matrix<std::complex<double> >& T, unsigned int NodeId )
    {
        // Unit vectors
        Vector<double> n(3), t(3), b(3);

        // Get node
        Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

        // Plasma parameters on node
        Vector<double> Bext = mpColdPlasma->Get_Bext_OnNode( pNode );

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

        if( nnorm > 0.0 )
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
    //* - Axisymmetric unitary vectors.
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
    //* - Trans(T) * StiffMatrix * T for PBC nodes.
    //**************************************************************************************************************************************************************
    void Modeler::PBC_FrontBack_TKT( Vector<int>& HONodesId, Vector<int>& IdVector, Matrix< std::complex<double> >& StiffMatrix )
    {
        if( mPBC_NodeNodePairs_Front.size() == 0 )
        {
            return;
        }

        for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if( mPBC_NodeNodePairs_Front.find( HONodesId[ local_i ] ) != mPBC_NodeNodePairs_Front.end() )
            {
                int bNodeId = mPBC_NodeNodePairs_Front[ HONodesId[ local_i ] ];
                int Step    = HONodesId.size();
                int DofSize = IdVector .size();

		        Node::Pointer pbNode = mpModel->GetNode( bNodeId-1 );

                // StiffnesMatrix * T
                for( int i=0; i<DofSize; i++ )
                {
                    StiffMatrix[ i ][ local_i          ] = mExp_jPha_FB * StiffMatrix[ i ][ local_i          ];
                    StiffMatrix[ i ][ local_i +   Step ] = mExp_jPha_FB * StiffMatrix[ i ][ local_i +   Step ];
                    StiffMatrix[ i ][ local_i + 2*Step ] = mExp_jPha_FB * StiffMatrix[ i ][ local_i + 2*Step ];
                    
                    if( ( local_i + 3*Step ) < DofSize )
                    {
                        StiffMatrix[ i ][ local_i + 3*Step ] = mExp_jPha_FB * StiffMatrix[ i ][ local_i + 3*Step ];
                    }
                }
            
                // Transpose(T) * StiffnesMatrix
                for( int i=0; i<DofSize; i++ )
                {
            	    StiffMatrix[ local_i          ][ i ] = mExp_jPha_FB * StiffMatrix[ local_i          ][ i ];
                    StiffMatrix[ local_i +   Step ][ i ] = mExp_jPha_FB * StiffMatrix[ local_i +   Step ][ i ];
                    StiffMatrix[ local_i + 2*Step ][ i ] = mExp_jPha_FB * StiffMatrix[ local_i + 2*Step ][ i ];

                    if( ( local_i + 3*Step ) < DofSize )
                    {
                        StiffMatrix[ local_i + 3*Step ][ i ] = mExp_jPha_FB * StiffMatrix[ local_i + 3*Step ][ i ];
                    }
                }   

                // Change DOFs IDs: Front (fixed) <-> Back (solve)
                if( mPotentials )
                {
                    IdVector[ local_i          ] = pbNode->pDofcAx()->EquationId();
		            IdVector[ local_i +   Step ] = pbNode->pDofcAy()->EquationId();
		            IdVector[ local_i + 2*Step ] = pbNode->pDofcAz()->EquationId();   
                
                    if( ( local_i + 3*Step ) < DofSize ) 
                    {
                        IdVector[ local_i + 3*Step ] = pbNode->pDofcVs()->EquationId();  
                    }
                }
                else
                {
                    IdVector[ local_i          ] = pbNode->pDofcEx()->EquationId();
		            IdVector[ local_i +   Step ] = pbNode->pDofcEy()->EquationId();
		            IdVector[ local_i + 2*Step ] = pbNode->pDofcEz()->EquationId();        
                }

                // Change node IDs: Front (fixed) <-> Back (solve)
                HONodesId[ local_i ] = pbNode->Id();  
            }
        }
    }

	//**************************************************************************************************************************************************************
    //* - Trans(T) * b for PBC nodes.
    //**************************************************************************************************************************************************************
	void Modeler::PBC_FrontBack_TR( Vector<int>& HONodesId, Vector<int>& IdVector, Vector< std::complex<double> >& ResVector )
	{
        if( mPBC_NodeNodePairs_Front.size() == 0 )
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if( mPBC_NodeNodePairs_Front.find( HONodesId[ local_i ] ) != mPBC_NodeNodePairs_Front.end() )
            {
                int bNodeId = mPBC_NodeNodePairs_Front[ HONodesId[ local_i ] ];
		        int Step    = HONodesId.size();

		        Node::Pointer pbNode = mpModel->GetNode( bNodeId-1 );

                // Transpose(T) * ResVector
                ResVector[ local_i          ] = mExp_jPha_FB * ResVector[ local_i          ];
                ResVector[ local_i +   Step ] = mExp_jPha_FB * ResVector[ local_i +   Step ];
                ResVector[ local_i + 2*Step ] = mExp_jPha_FB * ResVector[ local_i + 2*Step ];

                // Change DOFs IDs: Front (fixed) <-> Back (solve)
                if( mPotentials )
                {
                    IdVector[ local_i          ] = pbNode->pDofcAx()->EquationId();
		            IdVector[ local_i +   Step ] = pbNode->pDofcAy()->EquationId();
		            IdVector[ local_i + 2*Step ] = pbNode->pDofcAz()->EquationId();   
                }
                else
                {
                    IdVector[ local_i          ] = pbNode->pDofcEx()->EquationId();
		            IdVector[ local_i +   Step ] = pbNode->pDofcEy()->EquationId();
		            IdVector[ local_i + 2*Step ] = pbNode->pDofcEz()->EquationId();        
                }
                
                // Change node IDs: Front (fixed) <-> Back (solve)
                HONodesId[ local_i ] = pbNode->Id();  
			}
        }
	}

	//**************************************************************************************************************************************************************
    //* - Trans(T) * StiffnesMatrix * T for PBC nodes.
    //**************************************************************************************************************************************************************
    void Modeler::PBC_RightLeft_TKT( Vector<int>& HONodesId, Vector<int>& IdVector, Matrix< std::complex<double> >& StiffMatrix )
    {
        if( mPBC_NodeNodePairs_Right.size() == 0 ) 
        {
            return;
        }

        for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if( mPBC_NodeNodePairs_Right.find( HONodesId[ local_i ] ) != mPBC_NodeNodePairs_Right.end() )
            {
                int lNodeId = mPBC_NodeNodePairs_Right[ HONodesId[ local_i ] ];
                int Step    = HONodesId.size();
                int DofSize = IdVector .size();

		        Node::Pointer plNode = mpModel->GetNode( lNodeId-1 );

                // Rotation matrix
                if( ( mIsRightPBCTilted == true ) && ( mEdgeElementsOn == false ) ) 
                {
                    double cos_A, sin_A;
                    
                    PBC_InvRot( HONodesId[ local_i ], cos_A, sin_A );

                    Matrix< std::complex<double> > RotMatrix( 3, 3 );
                    
                    RotMatrix[0][0] = cos_A; RotMatrix[0][1] =-sin_A; RotMatrix[0][2] = 0.0;
                    RotMatrix[1][0] = sin_A; RotMatrix[1][1] = cos_A; RotMatrix[1][2] = 0.0;
                    RotMatrix[2][0] =   0.0; RotMatrix[2][1] =   0.0; RotMatrix[2][2] = 1.0;
                    
                    TransT_K_T( StiffMatrix, RotMatrix, HONodesId, local_i );
                }

                // StiffnesMatrix * T
                for( int i=0; i<DofSize; i++ )
                {
                    StiffMatrix[ i ][ local_i          ] = mExp_jPha_RL * StiffMatrix[ i ][ local_i          ];
                    StiffMatrix[ i ][ local_i +   Step ] = mExp_jPha_RL * StiffMatrix[ i ][ local_i +   Step ];
                    StiffMatrix[ i ][ local_i + 2*Step ] = mExp_jPha_RL * StiffMatrix[ i ][ local_i + 2*Step ];
                    
                    if( ( local_i + 3*Step ) < DofSize )
                    {
                        StiffMatrix[ i ][ local_i + 3*Step ] = mExp_jPha_RL * StiffMatrix[ i ][ local_i + 3*Step ];
                    }
                }
            
                // Transpose(T) * StiffnesMatrix
                for( int i=0; i<DofSize; i++ )
                {
            	    StiffMatrix[ local_i          ][ i ] = mExp_jPha_RL * StiffMatrix[ local_i          ][ i ];
                    StiffMatrix[ local_i +   Step ][ i ] = mExp_jPha_RL * StiffMatrix[ local_i +   Step ][ i ];
                    StiffMatrix[ local_i + 2*Step ][ i ] = mExp_jPha_RL * StiffMatrix[ local_i + 2*Step ][ i ];

                    if( ( local_i + 3*Step ) < DofSize )
                    {
                        StiffMatrix[ local_i + 3*Step ][ i ] = mExp_jPha_RL * StiffMatrix[ local_i + 3*Step ][ i ];
                    }
                } 

                // Change DOFs IDs: Right (fixed) <-> Left (solve)
                if( mPotentials )
                {
                    IdVector[ local_i          ] = plNode->pDofcAx()->EquationId();
		            IdVector[ local_i +   Step ] = plNode->pDofcAy()->EquationId();
		            IdVector[ local_i + 2*Step ] = plNode->pDofcAz()->EquationId();   
                
                    if( ( local_i + 3*Step ) < DofSize ) 
                    {
                        IdVector[ local_i + 3*Step ] = plNode->pDofcVs()->EquationId();  
                    }
                }
                else
                {
                    IdVector[ local_i          ] = plNode->pDofcEx()->EquationId();
		            IdVector[ local_i +   Step ] = plNode->pDofcEy()->EquationId();
		            IdVector[ local_i + 2*Step ] = plNode->pDofcEz()->EquationId();        
                }

                // Change node IDs: Right (fixed) <-> Left (solve)
                HONodesId[ local_i ] = plNode->Id();  
            }
        }
    }

	//**************************************************************************************************************************************************************
    //* - Trans(T) * b for PBC nodes.
    //**************************************************************************************************************************************************************
	void Modeler::PBC_RightLeft_TR( Vector<int>& HONodesId, Vector<int>& IdVector, Vector< std::complex<double> >& ResVector )
	{
	    if( mPBC_NodeNodePairs_Right.size() == 0 ) 
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if( mPBC_NodeNodePairs_Right.find( HONodesId[ local_i ] ) != mPBC_NodeNodePairs_Right.end() )
            {
                int lNodeId = mPBC_NodeNodePairs_Right[ HONodesId[ local_i ] ];
                int Step    = HONodesId.size();
		        
		        Node::Pointer plNode = mpModel->GetNode( lNodeId-1 );

                // Rotation matrix
                if( ( mIsRightPBCTilted == true ) && ( mEdgeElementsOn == false ) ) 
                {
                    double cos_A, sin_A;
                    
		            PBC_InvRot( HONodesId[ local_i ], cos_A, sin_A );

                    Matrix< std::complex<double> > RotMatrix( 3, 3 );
                    
                    RotMatrix[0][0] = cos_A; RotMatrix[0][1] =-sin_A; RotMatrix[0][2] = 0.0;
                    RotMatrix[1][0] = sin_A; RotMatrix[1][1] = cos_A; RotMatrix[1][2] = 0.0;
                    RotMatrix[2][0] =   0.0; RotMatrix[2][1] =   0.0; RotMatrix[2][2] = 1.0;  

                    TransT_B( ResVector, RotMatrix, HONodesId, local_i );
                }

                // Transpose(T) * ResVector
                ResVector[ local_i          ] = mExp_jPha_RL * ResVector[ local_i          ];
                ResVector[ local_i +   Step ] = mExp_jPha_RL * ResVector[ local_i +   Step ];
                ResVector[ local_i + 2*Step ] = mExp_jPha_RL * ResVector[ local_i + 2*Step ];

                // Change DOFs IDs: Right (fixed) <-> Left (solve)
                if( mPotentials )
                {
                    IdVector[ local_i          ] = plNode->pDofcAx()->EquationId();
		            IdVector[ local_i +   Step ] = plNode->pDofcAy()->EquationId();
		            IdVector[ local_i + 2*Step ] = plNode->pDofcAz()->EquationId();   
                }
                else
                {
                    IdVector[ local_i          ] = plNode->pDofcEx()->EquationId();
		            IdVector[ local_i +   Step ] = plNode->pDofcEy()->EquationId();
		            IdVector[ local_i + 2*Step ] = plNode->pDofcEz()->EquationId();        
                }

                // Change node IDs: Right (fixed) <-> Left (solve)
                HONodesId[ local_i ] = plNode->Id();  
            }
        }
	}

    //****************************************************************************************************************************************************
    // - Applies Trans(T) * StiffMatrix * T to eliminate one of the nodes in a contact pair.
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Matrix_TKT( Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix )
    {
		if( mContactPairs.size() == 0 ) 
        {
            return;
        }

        for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            if( Is_DisctAVNode( HONodesId[ local_i ] ) ) 
            {
                Apply_DisctAV_TKT_AV( local_i, HONodesId, DIdVector, StiffMatrix );
            }            
            else if( mContactPairs.find( HONodesId[ local_i ] ) != mContactPairs.end() )
            {
                if( mPotentials )
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
    // - Applies TKT to discontinuous A nodes in AV potentials fromulation.
    //****************************************************************************************************************************************************
    void Modeler::Apply_DisctAV_TKT_AV( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix )
    {
        // Get node Id
        int NodeId = HONodesId[ local_i ];

        // Get rotation matrix
        Matrix< std::complex<double> > RotMatrix;

        Get_Rotation_Matrix_DisctAV( RotMatrix, NodeId );
        
        // Apply rotation
        TransT_K_T( StiffMatrix, RotMatrix, HONodesId, local_i );
        
        // If it is a CtC node then change tangential DOFs Ids
        if( mContactPairs.find( NodeId ) != mContactPairs.end() )
        {
            Node::Pointer pNode = mpModel->GetNode( mContactPairs[ NodeId ].ctcNode - 1 );
            
            int Step    = HONodesId.size();
            int DofSize = DIdVector.size();

            if( mNormals.find( NodeId ) != mNormals.end() )
            {
                DIdVector[ local_i          ] = pNode->pDofcAx()->EquationId();
                DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();            
            }
            else
            {
                DIdVector[ local_i +   Step ] = pNode->pDofcAy()->EquationId();
                DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();            
            }

            if( ( local_i + 3*Step ) < DofSize ) 
            {
                DIdVector[ local_i + 3*Step ] = pNode->pDofcVs()->EquationId();
            }   
        }
    }

    //****************************************************************************************************************************************************
    // - Calculates rotation matrix for discontinuous nodes in AV potential formulation.
    //****************************************************************************************************************************************************
    void Modeler::Get_Rotation_Matrix_DisctAV( Matrix< std::complex<double> >& RotMatrix, int NodeId )
    {
        // Set rotation basis
        Vector<double> n( 3 );
        Vector<double> t( 3 );
        Vector<double> b( 3 );
        
        if( mNormals.find( NodeId ) != mNormals.end() )
        {
            n = mNormals[ NodeId ];

            if( mContactPairs.find( NodeId ) != mContactPairs.end() ) 
            {
                t = mContactPairs[ NodeId ].ctcNormal;
            }
            else
            {
                t = mCtCNormals[ NodeId ];
            }

            b[ 0 ] = n[ 1 ]*t[ 2 ] - n[ 2 ]*t[ 1 ];
            b[ 1 ] = n[ 2 ]*t[ 0 ] - n[ 0 ]*t[ 2 ];
            b[ 2 ] = n[ 0 ]*t[ 1 ] - n[ 1 ]*t[ 0 ];
        }
        else
        {
            if( mContactPairs.find( NodeId ) != mContactPairs.end() ) 
            {
                n = mContactPairs[ NodeId ].ctcNormal;
            }
            else          
            {
                n = mCtCNormals[ NodeId ];
            }

            TangencialCoordinates( n, t, b );
        }

        // Set rotation matrix
        RotMatrix.Resize( 3, 3 );
                
        RotMatrix[ 0 ][ 0 ] = n[ 0 ]; RotMatrix[ 0 ][ 1 ] = n[ 1 ]; RotMatrix[ 0 ][ 2 ] = n[ 2 ];
        RotMatrix[ 1 ][ 0 ] = t[ 0 ]; RotMatrix[ 1 ][ 1 ] = t[ 1 ]; RotMatrix[ 1 ][ 2 ] = t[ 2 ];
        RotMatrix[ 2 ][ 0 ] = b[ 0 ]; RotMatrix[ 2 ][ 1 ] = b[ 1 ]; RotMatrix[ 2 ][ 2 ] = b[ 2 ];    
    }

    //****************************************************************************************************************************************************
    // - Applies Trans(T) * StiffMatrix * T to eliminate one of the nodes in a contact pair ( E field fromulation ).
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_TKT_Ef( int local_i, Vector<int>& HONodesId, Vector<int>& DIdVector, Matrix< std::complex<double> >& StiffMatrix )
    {   
        // Step size 
        int Step = HONodesId.size();
 
        // If a nodal element formulation is active then rotate stiffness matrix
        if( mEdgeElementsOn == false )
        {
            // Get contact matrix
            Matrix< std::complex<double> > CtCMatrix;

            if( mColdPlasmaMode )
            {
                Get_Contact_Matrix_ColdPlasma( CtCMatrix, HONodesId[ local_i ] );
            }
            else
            {
                Get_Contact_Matrix_FullWave( CtCMatrix, HONodesId[ local_i ] );
            }
            
            // DOF size
            int DofSize = 3*HONodesId.size();
            
            // Inv(T) * StiffnesMatrix * T     
            std::complex<double> dof_x;
            std::complex<double> dof_y;
            std::complex<double> dof_z;
            
            // StiffnesMatrix * T
            for( int i=0; i<DofSize; i++ )
            {
                dof_x = StiffMatrix[ i ][ local_i          ];
                dof_y = StiffMatrix[ i ][ local_i +   Step ];
                dof_z = StiffMatrix[ i ][ local_i + 2*Step ];
            
                StiffMatrix[ i ][ local_i          ] = CtCMatrix[ 0 ][ 0 ]*dof_x + CtCMatrix[ 1 ][ 0 ]*dof_y + CtCMatrix[ 2 ][ 0 ]*dof_z;
                StiffMatrix[ i ][ local_i +   Step ] = CtCMatrix[ 0 ][ 1 ]*dof_x + CtCMatrix[ 1 ][ 1 ]*dof_y + CtCMatrix[ 2 ][ 1 ]*dof_z;
                StiffMatrix[ i ][ local_i + 2*Step ] = CtCMatrix[ 0 ][ 2 ]*dof_x + CtCMatrix[ 1 ][ 2 ]*dof_y + CtCMatrix[ 2 ][ 2 ]*dof_z;
            }
            
            // Inv(T) * StiffnesMatrix
            for( int i=0; i<DofSize; i++ )
            {
                dof_x = StiffMatrix[ local_i          ][ i ];
                dof_y = StiffMatrix[ local_i +   Step ][ i ];
                dof_z = StiffMatrix[ local_i + 2*Step ][ i ];
            
            	StiffMatrix[ local_i          ][ i ] = CtCMatrix[ 0 ][ 0 ]*dof_x + CtCMatrix[ 1 ][ 0 ]*dof_y + CtCMatrix[ 2 ][ 0 ]*dof_z;
                StiffMatrix[ local_i +   Step ][ i ] = CtCMatrix[ 0 ][ 1 ]*dof_x + CtCMatrix[ 1 ][ 1 ]*dof_y + CtCMatrix[ 2 ][ 1 ]*dof_z;
                StiffMatrix[ local_i + 2*Step ][ i ] = CtCMatrix[ 0 ][ 2 ]*dof_x + CtCMatrix[ 1 ][ 2 ]*dof_y + CtCMatrix[ 2 ][ 2 ]*dof_z;
            }           
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
    // - Calculates contact matrix for dielectric interfaces.
    //****************************************************************************************************************************************************
    void Modeler::Get_Contact_Matrix_FullWave( Matrix< std::complex<double> >& CtCMatrix, int RNodesId )
    {
        ContactPairData& pContactPair = mContactPairs[ RNodesId ];
 
        Properties::Pointer Properties;
        
        double eo   = 8.8541878176e-12 ;
        double freq = mProblemFrequency;
        
        double sgm_real;
        double sgm_imag;

        double eps_real;
        double eps_imag;
        
        // Material R
        Properties = mpModel->GetProperties( pContactPair.matR );
        
        sgm_real = (*Properties)( REAL_ELECTRIC_CONDUCTIVITY );
        sgm_imag = (*Properties)( IMAG_ELECTRIC_CONDUCTIVITY );

        eps_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
        eps_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
        
        std::complex<double> cEpsR( eps_real - ( sgm_imag / freq ), eps_imag + ( sgm_real / freq ) );
        
        // Material L
        Properties = mpModel->GetProperties( pContactPair.matL );
        
        sgm_real = (*Properties)( REAL_ELECTRIC_CONDUCTIVITY );
        sgm_imag = (*Properties)( IMAG_ELECTRIC_CONDUCTIVITY );

        eps_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY ) * eo;
        eps_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;
        
        std::complex<double> cEpsL( eps_real - ( sgm_imag / freq ), eps_imag + ( sgm_real / freq ) );
                
        std::complex<double> rmat = ( cEpsL / cEpsR ) - 1.0;
        
        Vector<double> n = pContactPair.ctcNormal;

        // Set contact matrix
        CtCMatrix.Resize( 3, 3 );
                                                                                                                                    
        CtCMatrix[0][0] = ( n[0]*n[0] ) * rmat + 1.0; CtCMatrix[0][1] = ( n[0]*n[1] ) * rmat      ; CtCMatrix[0][2] = ( n[0]*n[2] ) * rmat      ;
        CtCMatrix[1][0] = ( n[1]*n[0] ) * rmat      ; CtCMatrix[1][1] = ( n[1]*n[1] ) * rmat + 1.0; CtCMatrix[1][2] = ( n[1]*n[2] ) * rmat      ;
        CtCMatrix[2][0] = ( n[2]*n[0] ) * rmat      ; CtCMatrix[2][1] = ( n[2]*n[1] ) * rmat      ; CtCMatrix[2][2] = ( n[2]*n[2] ) * rmat + 1.0;
    }

    //****************************************************************************************************************************************************
    // - Eliminates one of the nodes in a contact pair ( AV potentials fromulation ).
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

        if( ( local_i + 3*Step ) < DofSize ) 
        {
            DIdVector[ local_i + 3*Step ] = pNode->pDofcVs()->EquationId();
        }   

        // Change node Id
        HONodesId[ local_i ] = pNode->Id();
    }

	//****************************************************************************************************************************************************
    // - Inv(T) * b in contact nodes.
    //****************************************************************************************************************************************************
    void Modeler::Apply_Contact_Vector_TR( Vector<int>& HONodesId, Vector<int>& DIdVector, Vector<std::complex<double> >& ResVector )
    {
        if( mContactPairs.size() == 0 ) 
        {
            return;
        }

		for( int local_i = 0; local_i < HONodesId.size(); local_i++ )
        {
            int NodeId = HONodesId[ local_i ];
            
            if( Is_DisctAVNode( NodeId ) ) 
            {
                Matrix< std::complex<double> > RotMatrix;

                Get_Rotation_Matrix_DisctAV( RotMatrix, NodeId );

                TransT_B( ResVector, RotMatrix, HONodesId, local_i );

                if( mContactPairs.find( NodeId ) != mContactPairs.end() )
                {
                    Node::Pointer pNode = mpModel->GetNode( mContactPairs[ NodeId ].ctcNode - 1 );
                    
                    int Step = HONodesId.size();
                    
                    if ( mNormals.find( NodeId ) != mNormals.end() )
                    {
                        DIdVector[ local_i          ] = pNode->pDofcAx()->EquationId();
                        DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();            
                    }
                    else
                    {
                        DIdVector[ local_i +   Step ] = pNode->pDofcAy()->EquationId();
                        DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();            
                    }
                }
            }                        
            else if( mContactPairs.find( NodeId ) != mContactPairs.end() )
            {
				Node::Pointer pNode = mpModel->GetNode( mContactPairs[ NodeId ].ctcNode - 1 );
                    
                int Step = HONodesId.size();
 
                if( mPotentials )
                {
				    DIdVector[ local_i          ] = pNode->pDofcAx()->EquationId();
				    DIdVector[ local_i +   Step ] = pNode->pDofcAy()->EquationId();
				    DIdVector[ local_i + 2*Step ] = pNode->pDofcAz()->EquationId();                
                }
                else
                {
                    if( mEdgeElementsOn == false )
                    {
                        Matrix< std::complex<double> > CtCMatrix; 

                        if( mColdPlasmaMode )
                        {
                            Get_Contact_Matrix_ColdPlasma( CtCMatrix, NodeId );
                        }
                        else
                        {
                            Get_Contact_Matrix_FullWave( CtCMatrix, NodeId );
                        }

                        TransT_B( ResVector, CtCMatrix.Transpose(), HONodesId, local_i );
                    }

				    DIdVector[ local_i          ] = pNode->pDofcEx()->EquationId();
				    DIdVector[ local_i +   Step ] = pNode->pDofcEy()->EquationId();
				    DIdVector[ local_i + 2*Step ] = pNode->pDofcEz()->EquationId();
                }

				HONodesId[ local_i ] = pNode->Id();
            }
        }
    }

    //****************************************************************************************************************************************************
    // - Calculates the contact matrix for Dielectric-Plasma interfaces.
    //****************************************************************************************************************************************************
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

        double sgm_real, sgm_imag;
        double epr_real, epr_imag;
        double mur_real, mur_imag;

        Matrix< std::complex<double> > TEn( 3, 3 );

        // Material R
        Properties = mpModel->GetProperties( pContactPair.matR );

        sgm_real = (*Properties)( REAL_ELECTRIC_CONDUCTIVITY );
        sgm_imag = (*Properties)( IMAG_ELECTRIC_CONDUCTIVITY );

        epr_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
        epr_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );

        mur_real = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
        mur_imag = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );

        // Plasma permittivity tensor
        if( ( sgm_real == 0.0 ) && ( sgm_imag == 0.0 ) &&
            ( epr_real == 1.0 ) && ( epr_imag == 0.0 ) && 
            ( mur_real == 1.0 ) && ( mur_imag == 0.0 )  )
        {
            mpColdPlasma->Get_PermittivityTensor_OnNode( TEn, prNode, freq );
        }
        // IHL permittivity tensor
        else
        {
            std::complex<double> cEps( eo*epr_real - ( sgm_imag / freq ), eo*epr_imag + ( sgm_real / freq ) );

            TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero;
            TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero;
            TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
        }

        // [TR] * {Er} = [TL] * {El}
        Matrix< std::complex<double> > TR( 3, 3 );

        TR[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0];
        TR[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1];
        TR[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2];

        TR[1][0] = t[0]; TR[1][1] = t[1]; TR[1][2] = t[2];
        TR[2][0] = b[0]; TR[2][1] = b[1]; TR[2][2] = b[2];

        // Material L
        Properties = mpModel->GetProperties( pContactPair.matL );

        sgm_real = (*Properties)( REAL_ELECTRIC_CONDUCTIVITY );
        sgm_imag = (*Properties)( IMAG_ELECTRIC_CONDUCTIVITY );

        epr_real = (*Properties)( REAL_ELECTRIC_PERMITTIVITY );
        epr_imag = (*Properties)( IMAG_ELECTRIC_PERMITTIVITY );

        mur_real = (*Properties)( REAL_MAGNETIC_PERMEABILITY );
        mur_imag = (*Properties)( IMAG_MAGNETIC_PERMEABILITY );

        // Plasma permittivity tensor
        if( ( sgm_real == 0.0 ) && ( sgm_imag == 0.0 ) &&
            ( epr_real == 1.0 ) && ( epr_imag == 0.0 ) && 
            ( mur_real == 1.0 ) && ( mur_imag == 0.0 )  )
        {
            mpColdPlasma->Get_PermittivityTensor_OnNode( TEn, plNode, freq );
        }
        // IHL permittivity tensor
        else
        {
            std::complex<double> cEps( eo*epr_real - ( sgm_imag / freq ), eo*epr_imag + ( sgm_real / freq ) );

            TEn[0][0] = cEps ; TEn[0][1] = cZero; TEn[0][2] = cZero;
            TEn[1][0] = cZero; TEn[1][1] = cEps ; TEn[1][2] = cZero;
            TEn[2][0] = cZero; TEn[2][1] = cZero; TEn[2][2] = cEps ;
        }

        // [TR] * {Er} = [TL] * {El}
        Matrix< std::complex<double> > TL( 3, 3 );

        TL[0][0] = n[0]*TEn[0][0] + n[1]*TEn[1][0] + n[2]*TEn[2][0];
        TL[0][1] = n[0]*TEn[0][1] + n[1]*TEn[1][1] + n[2]*TEn[2][1];
        TL[0][2] = n[0]*TEn[0][2] + n[1]*TEn[1][2] + n[2]*TEn[2][2];

        TL[1][0] = t[0]; TL[1][1] = t[1]; TL[1][2] = t[2];
        TL[2][0] = b[0]; TL[2][1] = b[1]; TL[2][2] = b[2];

        // {Er} = ( [invTR] * [TL] ) * {El}
        Matrix< std::complex<double> > invTR( 3, 3 );

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
        CtCMatrix.Resize( 3, 3 );

        CtCMatrix = invTR * TL;
    }

    //*************************************************************************************************************************************
    // - Returns the determinant of a square matrix M of size n (for complex matrices).
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
    // - Fixes all DOFs in a node to zero.
    //**********************************************************************************************************************
    void Modeler::Fix_All_DOFs_In_Node( Node::Pointer pNode )
    {
        std::complex<double> cZero( 0.00, 0.00 );
        
        if ( mPotentials )
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
    // 
	// - PEC, PMC, PBC and projection elements.
    // 
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//**********************************************************************************************************************
    // - Assigns each PBC node to its corresponding std::map.
    //**********************************************************************************************************************
    void Modeler::Generate_PBC_Element( int* NodesId, unsigned int PropertiesId )
    {
        // Properties of the PBC element
        Properties::Pointer properties = mpModel->GetProperties( PropertiesId );

		int    FaceTypeId = (*properties)( COMPLEX_IBC )[0];
        double Phase_Diff = (*properties)( COMPLEX_IBC )[1];

        // Phase difference between faces
        if( FaceTypeId == 11 || FaceTypeId == 12 ) 
        {
            mExp_jPha_FB = std::complex<double>( cos( Phase_Diff ), sin( Phase_Diff ) );
        }
        
        if( FaceTypeId == 21 || FaceTypeId == 22 ) 
        {
            mExp_jPha_RL = std::complex<double>( cos( Phase_Diff ), sin( Phase_Diff ) );
        }
        
        // Each material is assigned to a different group to speed the PBC set-up.
        int MaterialId = NodesId[3];

		Vector<int> HONodesId;

        Push_HONodes_OnSurface( NodesId, HONodesId );

		if     ( FaceTypeId == 11 ) mPBCFrontElements[ MaterialId ].push_back( HONodesId );
		else if( FaceTypeId == 12 ) mPBCBackElements [ MaterialId ].push_back( HONodesId );
		else if( FaceTypeId == 21 ) mPBCRightElements[ MaterialId ].push_back( HONodesId );
		else if( FaceTypeId == 22 ) mPBCLeftElements [ MaterialId ].push_back( HONodesId );
	}

	//**********************************************************************************************************************
    // - Assigns nodes in a PBC surface to nodes in its corresponding PBC surface (Front-Back, Right-Left).
	// - Geometry must be placed carefully:
	//   1) Front-Back surfaces in the XY-plane.
	//   2) Left surface in the 00-YZ-Plane.
	//	 3) For cylindrical symmetry the central axis must be placed along the Z axis.
    //**********************************************************************************************************************
	void Modeler::Set_PBC()
	{
        // Setting PBC for Front-Back surfaces
		if( mPBCFrontElements.size() > 0 ) 
        {
            Set_FrontBack_PBC();
        }

        // Setting PBC for Right-Left surfaces
        if( mPBCRightElements.size() > 0 ) 
        {
            Set_RightLeft_PBC();
        }
	}

	//**********************************************************************************************************************
    // - Assigns nodes in the Front PBC surface to nodes in the Back PBC surface.
	// - Front-Back surfaces must be placed in XY-planes.
    //**********************************************************************************************************************
    void Modeler::Set_FrontBack_PBC()
    {
	    time_t bStart; time( &bStart ); 
        
        // Loop through materials in Back PBC
        std::map<unsigned int, Vector< Vector<int> > >::iterator it_Back;

        for( it_Back = mPBCBackElements.begin(); it_Back != mPBCBackElements.end(); ++it_Back )
        {
            // Material ID
			int pbcMaterial = it_Back->first;

            // Group of Back elements made of the same material (pbcMaterial)
			Vector< Vector<int> > elemSetBack( it_Back->second );

            // Group of Frontal elements made of the same material (pbcMaterial)
			Vector< Vector<int> > elemSetFront( mPBCFrontElements[ pbcMaterial ] );

            // Filling an ordered unique list with the Front nodes
            std::list<int> FrontNodesList;

            Vector< Vector<int> >::iterator it_eF;

            for( it_eF = elemSetFront.begin(); it_eF != elemSetFront.end(); ++it_eF )
            {
                Vector<int>::iterator it_nF;

                for( it_nF = (*it_eF).begin(); it_nF != (*it_eF).end(); ++it_nF )
				{
                    FrontNodesList.push_back( *it_nF );
                }
            }

            FrontNodesList.sort  ();
            FrontNodesList.unique();

            // Filling an ordered unique list with the Back nodes 
            std::list<int> BackNodesList;

            Vector< Vector<int> >::iterator it_eB;

            for( it_eB = elemSetBack.begin(); it_eB != elemSetBack.end(); ++it_eB )
            {
                Vector<int>::iterator it_nB;

                for( it_nB = (*it_eB).begin(); it_nB != (*it_eB).end(); ++it_nB )
				{
                    BackNodesList.push_back( *it_nB );
                }
            }

            BackNodesList.sort  ();
            BackNodesList.unique();

            // Get the closest Back node to the Front node proyection on Back plane
            Node::Pointer pFNode;
            Node::Pointer pBNode;

            double MinDistance_Old; 
            double MinDistance_New; 

            double FX, FY;
            double BX, BY;

            int Closest_BNodeId; 
            int FNodeId;

            std::list<int> BackList_Updated;

            std::list<int>::iterator lit;

            Vector<int> FrontNodesVector; 

            for ( lit = FrontNodesList.begin(); lit != FrontNodesList.end(); ++lit )
            {
                FrontNodesVector.push_back( *lit );
            }

            omp_set_num_threads( mItSolverNumThreads );

		    int nthreads = omp_get_num_threads();

            #pragma omp parallel private( MinDistance_Old, MinDistance_New, Closest_BNodeId, FNodeId, FX, FY, BX, BY, pFNode, pBNode, lit, BackList_Updated ) shared( FrontNodesVector, BackNodesList )   
            {
                #pragma omp for 
                for( int i=0; i<FrontNodesVector.size(); i++ )
                {
                    FNodeId = FrontNodesVector[i];
                
                    pFNode = mpModel->GetNode( FNodeId-1 );

                    FX = pFNode->X();
                    FY = pFNode->Y();

                    Fix_All_DOFs_In_Node( pFNode );
                   
                    MinDistance_Old = 1e12;

                    Closest_BNodeId = -1;

                    #pragma omp critical
                    {
                       BackList_Updated = BackNodesList;
                    }

                    for( lit = BackList_Updated.begin(); lit != BackList_Updated.end(); lit++ )
                    {
                        pBNode = mpModel->GetNode( (*lit)-1 );
                    
                        BX = pBNode->X();
                        BY = pBNode->Y();
                    
                        MinDistance_New = sqrt( ( FX - BX )*( FX - BX ) + ( FY - BY )*( FY - BY ) );
                    
                        if ( MinDistance_New <= MinDistance_Old )
                        {
                            Closest_BNodeId = (*lit);
                            
                            MinDistance_Old = MinDistance_New;
                        }
                    }
                    
                    #pragma omp critical
                    {
                        // No partner node left ( different number of nodes on the periodic faces )
                        if( Closest_BNodeId < 0 ) Send_Error_Msg( "Periodic boundary condition: no matching node on the opposite face", "Please, use identical ( copied ) meshes on both periodic faces.", 1 );

                        mPBC_NodeNodePairs_Front[ FNodeId ] = Closest_BNodeId;

                        BackNodesList.remove( Closest_BNodeId );
                    }
                }
            }
        }

        time_t bEnd; time( &bEnd );
        
        std::cout << "Front-Back PBC set in " << difftime( bEnd, bStart ) << " seconds." << std::endl;

        // Clearing PBC element maps
	    mPBCFrontElements.clear(); 
        std::map< unsigned int, Vector< Vector<int> > >().swap( mPBCFrontElements );

		mPBCBackElements.clear(); 
        std::map< unsigned int, Vector< Vector<int> > >().swap( mPBCBackElements );
    }

	//***********************************************************************************************************************************************
    // - Assigns nodes in the Right PBC surface to nodes in the Left PBC surface.
	// - Left surface must be placed in the 00-YZ-Plane.
	// - For cylindrical symmetry (cyclic periodic) the central axis must be placed along the Z axis.
    //***********************************************************************************************************************************************
	void Modeler::Set_RightLeft_PBC()
	{
        time_t bStart; time( &bStart ); 

        // Loop through materials in Left PBC
		std::map<unsigned int, Vector< Vector<int> > >::iterator it_Left;

        for( it_Left = mPBCLeftElements.begin(); it_Left != mPBCLeftElements.end(); ++it_Left )
        {
            // Material ID
			int pbcMaterial = it_Left->first;

            // Group of Left elements made of the same material (pbcMaterial)
			Vector< Vector<int> > elemSetLeft( it_Left->second );

            // Group of Right elements made of the same material (pbcMaterial)
			Vector< Vector<int> > elemSetRight( mPBCRightElements[ pbcMaterial ] );

            // List of nodes (sorted and unique) in the Right PBC
            std::list<int> RightNodesList;
            
            // Filling an ordered unique list with the Right nodes 
            Vector< Vector<int> >::iterator it_eR;

            for( it_eR = elemSetRight.begin(); it_eR != elemSetRight.end(); ++it_eR )
            {
                Vector<int>::iterator it_nR;

                for( it_nR = (*it_eR).begin(); it_nR != (*it_eR).end(); ++it_nR )
				{
                    RightNodesList.push_back( *it_nR );
                }
            }

            RightNodesList.sort  ();
            RightNodesList.unique();

            // Filling an ordered unique list with the Left nodes 
            std::list<int> LeftNodesList;
            
            Vector< Vector<int> >::iterator it_eL;

            for( it_eL = elemSetLeft.begin(); it_eL != elemSetLeft.end(); ++it_eL )
            {
                Vector<int>::iterator it_nL;

                for( it_nL = (*it_eL).begin(); it_nL != (*it_eL).end(); ++it_nL )
				{
                    LeftNodesList.push_back( *it_nL );
                }
            }

            LeftNodesList.sort  ();
            LeftNodesList.unique();

            // Get the closest Left node to the Right node proyection on Left plane
            Node::Pointer pRNode;
            Node::Pointer pLNode;

            double MinDistance_Old; 
            double MinDistance_New; 

            double RX, RY;
            double LX, LY;

            int Closest_LNodeId; 
            int RNodeId;

            std::complex<double> cZero( 0.00, 0.00 );

            std::list<int> LeftList_Updated;

            std::list<int>::iterator lit;

            Vector<int> RightNodesVector; 

            for ( lit = RightNodesList.begin(); lit != RightNodesList.end(); ++lit )
            {
                RightNodesVector.push_back( *lit );
            }

            omp_set_num_threads( mItSolverNumThreads );

		    int nthreads = omp_get_num_threads();

            #pragma omp parallel private( MinDistance_Old, MinDistance_New, Closest_LNodeId, RNodeId, RX, RY, LX, LY, pRNode, pLNode, lit, LeftList_Updated ) shared( RightNodesVector, LeftNodesList )   
            {
                #pragma omp for 
                for( int i=0; i<RightNodesVector.size(); i++ )
                {
                    RNodeId = RightNodesVector[i];
                
                    pRNode = mpModel->GetNode( RNodeId-1 );
                    
                    RX = pRNode->Z();

                    if ( mIsRightPBCTilted ) 
                    {
                        RY = sqrt( pRNode->X() * pRNode->X() + pRNode->Y() * pRNode->Y() );
                    }
                    else
                    {
                        RY = pRNode->Y();
                    }

                    if ( ( mIsRightPBCTilted ) && ( abs( pRNode->X() ) < mGeoTolerance ) && ( abs( pRNode->Y() ) < mGeoTolerance ) )
                    {
                        if ( mNormals.find( RNodeId ) == mNormals.end() )
                        {                
                            if ( mPotentials )
                            {
                                (*mpModel)( cAx, *pRNode ) = cZero; pRNode->pDofcAx()->FixDof();
                                (*mpModel)( cAy, *pRNode ) = cZero; pRNode->pDofcAy()->FixDof();
                            }                            
                            else                         
                            {                            
                                (*mpModel)( cEx, *pRNode ) = cZero; pRNode->pDofcEx()->FixDof();
                                (*mpModel)( cEy, *pRNode ) = cZero; pRNode->pDofcEy()->FixDof();                  
                            }                                                        
                        }
                    }
                    else
                    {
                        Fix_All_DOFs_In_Node( pRNode );
                        
                        MinDistance_Old = 1e12;

                    Closest_LNodeId = -1;

                        #pragma omp critical
                        {
                           LeftList_Updated = LeftNodesList;
                        }
                        
                        for( lit = LeftList_Updated.begin(); lit != LeftList_Updated.end(); lit++ )
                        {
                            pLNode = mpModel->GetNode( (*lit)-1 );
                        
                            LX = pLNode->Z();
                            LY = pLNode->Y();
                        
                            MinDistance_New = sqrt( ( RX - LX )*( RX - LX ) + ( RY - LY )*( RY - LY ) );
                        
                            if ( MinDistance_New <= MinDistance_Old )
                            {
                                Closest_LNodeId = (*lit);
                                
                                MinDistance_Old = MinDistance_New;
                            }
                        }
                        
                        #pragma omp critical
                        {
                            // No partner node left ( different number of nodes on the periodic faces )
                        if( Closest_LNodeId < 0 ) Send_Error_Msg( "Periodic boundary condition: no matching node on the opposite face", "Please, use identical ( copied ) meshes on both periodic faces.", 1 );

                        mPBC_NodeNodePairs_Right[ RNodeId ] = Closest_LNodeId;
                            
                            LeftNodesList.remove( Closest_LNodeId );
                        }
                    }
                }
            }
        }

        time_t bEnd; time( &bEnd );
        
        std::cout << "Left-Right PBC set in " << difftime( bEnd, bStart ) << " seconds." << std::endl;
 
        // Clearing PBC element maps
		mPBCRightElements.clear();
		std::map< unsigned int, Vector< Vector<int> > >().swap( mPBCRightElements );
        
		mPBCLeftElements.clear();
		std::map< unsigned int, Vector< Vector<int> > >().swap( mPBCLeftElements );
	}

	//*****************************************************************************************************************
    //* - Calculates the volume integral of fields.
    //*****************************************************************************************************************
    void Modeler::Generate_Vol_Intgrl_Element_FullWave( int* NodesId, unsigned int PropertiesId )
    {
		// Element nodes
        Vector<Node::Pointer> pNodes; 
        
        Get_Volume_Element_Nodes( NodesId, pNodes );

        pNodes.resize( 4 );

        // Volume ID
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*pProperties)( COMPLEX_IBC );

		double VolumeID = cIBC[ 0 ];

        // Add nodes IDs to volume elements list 
        if( ( mExportFields ) && ( mVolIntg_Nodes.size() == 0 ) )
        {
            Vector<int> HONodesId; 
            
            Push_HONodes_OnVolume( NodesId, HONodesId );

            HONodesId.resize( 4 );
            
            mVolIntg_Elements[ VolumeID ].push_back( HONodesId );           
        }

        // Volume integral
		double EleVolume  = VolumeOfTetrahedra( pNodes );
        double VolIntg_Ni = EleVolume / 4.0;

        // Field volume integrals vectors
        std::complex<double> cZero( 0.0, 0.0 );

		Vector< std::complex<double> > VolIntg_E( 5, cZero );
        Vector< std::complex<double> > VolIntg_H( 5, cZero );
        Vector< std::complex<double> > VolIntg_B( 5, cZero );
        Vector< std::complex<double> > VolIntg_J( 5, cZero );
        Vector< std::complex<double> > VolIntg_F( 5, cZero );

        // Volume integral of fields
        for( int i = 0; i< pNodes.size(); i++ )
        {
			// Electric field
            Vector<double> E_real = (*mpModel)( REAL_E, *pNodes[ i ] );
            Vector<double> E_imag = (*mpModel)( IMAG_E, *pNodes[ i ] );

            Add_NodalFieldIntegral( E_real, E_imag, VolIntg_Ni, VolIntg_E );

            // Magnetic field
			Vector<double> H_real = (*mpModel)( REAL_H, *pNodes[ i ] );
            Vector<double> H_imag = (*mpModel)( IMAG_H, *pNodes[ i ] );

            Add_NodalFieldIntegral( H_real, H_imag, VolIntg_Ni, VolIntg_H );
                                                        
            // Magnetic flux density
			Vector<double> B_real = (*mpModel)( REAL_B, *pNodes[ i ] );
            Vector<double> B_imag = (*mpModel)( IMAG_B, *pNodes[ i ] );
                                                        
            Add_NodalFieldIntegral( B_real, B_imag, VolIntg_Ni, VolIntg_B );

			// Current density
			Vector<double> J_real = (*mpModel)( REAL_J, *pNodes[ i ] );
            Vector<double> J_imag = (*mpModel)( IMAG_J, *pNodes[ i ] );

            Add_NodalFieldIntegral( J_real, J_imag, VolIntg_Ni, VolIntg_J );

            // Lorentz force ( J x B )
            Vector< std::complex<double> > cE( 3 ), cB( 3 ), cJ( 3 );

            for( int i=0; i<3; i++ )
            {
                cE[ i ] = std::complex<double>( E_real[ i ], E_imag[ i ] );
                cB[ i ] = std::complex<double>( B_real[ i ], B_imag[ i ] );
                cJ[ i ] = std::complex<double>( J_real[ i ], J_imag[ i ] );
            }

            // F = J x B
            Vector< std::complex<double> > cF( 3 );

            cF[ 0 ] = cJ[ 1 ] * cB[ 2 ] - cJ[ 2 ] * cB[ 1 ];
            cF[ 1 ] = cJ[ 2 ] * cB[ 0 ] - cJ[ 0 ] * cB[ 2 ];
            cF[ 2 ] = cJ[ 0 ] * cB[ 1 ] - cJ[ 1 ] * cB[ 0 ];

            // mod( J x B )
            double modF; 
            
            modF = sqrt( std::real( cF[ 0 ] )*std::real( cF[ 0 ] ) + std::imag( cF[ 0 ] )*std::imag( cF[ 0 ] ) +
                         std::real( cF[ 1 ] )*std::real( cF[ 1 ] ) + std::imag( cF[ 1 ] )*std::imag( cF[ 1 ] ) +
                         std::real( cF[ 2 ] )*std::real( cF[ 2 ] ) + std::imag( cF[ 2 ] )*std::imag( cF[ 2 ] ) );

            // J* . E
            std::complex<double> conjJpE;
            
            conjJpE = std::conj( cJ[ 0 ] )*cE[ 0 ] + std::conj( cJ[ 1 ] )*cE[ 1 ] + std::conj( cJ[ 2 ] )*cE[ 2 ] ;

            VolIntg_F[ 0 ] += cF[ 0 ] * VolIntg_Ni;
            VolIntg_F[ 1 ] += cF[ 1 ] * VolIntg_Ni;
            VolIntg_F[ 2 ] += cF[ 2 ] * VolIntg_Ni;
            VolIntg_F[ 3 ] += modF    * VolIntg_Ni;
            VolIntg_F[ 4 ] += conjJpE * VolIntg_Ni;
        }

		if( mVolIntg_E.find( VolumeID ) != mVolIntg_E.end() )
		{
			mVolume_m3[ VolumeID ] += EleVolume;
			mVolIntg_E[ VolumeID ] += VolIntg_E;
			mVolIntg_H[ VolumeID ] += VolIntg_H;
            mVolIntg_B[ VolumeID ] += VolIntg_B;
			mVolIntg_J[ VolumeID ] += VolIntg_J;
            mVolIntg_F[ VolumeID ] += VolIntg_F;
		}                              
        else                           
		{                              
			mVolume_m3[ VolumeID ] = EleVolume;
			mVolIntg_E[ VolumeID ] = VolIntg_E;
			mVolIntg_H[ VolumeID ] = VolIntg_H;
            mVolIntg_B[ VolumeID ] = VolIntg_B;
			mVolIntg_J[ VolumeID ] = VolIntg_J;
            mVolIntg_F[ VolumeID ] = VolIntg_F;
		}                     
	}

    //***************************************************************************************************************************
    //* - Adds nodal field integrals.
    //***************************************************************************************************************************
    void Modeler::Add_NodalFieldIntegral( Vector<double>& F_real, Vector<double>& F_imag, double Intg_Ni, Vector< std::complex<double> >& FieldIntg_F )
    {
	    // Field components
        std::complex<double> cFx( F_real[ 0 ], F_imag[ 0 ] );
	    std::complex<double> cFy( F_real[ 1 ], F_imag[ 1 ] );
	    std::complex<double> cFz( F_real[ 2 ], F_imag[ 2 ] );
        
        // mod( F ), mod( F )^2
	    double modF, modF2;
        
        modF2 = F_real[ 0 ]*F_real[ 0 ] + F_imag[ 0 ]*F_imag[ 0 ] +
                F_real[ 1 ]*F_real[ 1 ] + F_imag[ 1 ]*F_imag[ 1 ] +
                F_real[ 2 ]*F_real[ 2 ] + F_imag[ 2 ]*F_imag[ 2 ] ;
        
        modF  = sqrt( modF2 );
        
        // Field integral
	    FieldIntg_F[ 0 ] += cFx   * Intg_Ni;
	    FieldIntg_F[ 1 ] += cFy   * Intg_Ni;
	    FieldIntg_F[ 2 ] += cFz   * Intg_Ni;
        FieldIntg_F[ 3 ] += modF  * Intg_Ni;
	    FieldIntg_F[ 4 ] += modF2 * Intg_Ni;
    }

    //***************************************************************************************************************************
    //* - Calculates the volume of a linear tetrahedra.
    //***************************************************************************************************************************
    double Modeler::VolumeOfTetrahedra( Vector< Node::Pointer >& pNodes )
    {
        double X1 = pNodes[ 0 ]->X(), Y1 = pNodes[ 0 ]->Y(), Z1 = pNodes[ 0 ]->Z();
		double X2 = pNodes[ 1 ]->X(), Y2 = pNodes[ 1 ]->Y(), Z2 = pNodes[ 1 ]->Z();
		double X3 = pNodes[ 2 ]->X(), Y3 = pNodes[ 2 ]->Y(), Z3 = pNodes[ 2 ]->Z();
		double X4 = pNodes[ 3 ]->X(), Y4 = pNodes[ 3 ]->Y(), Z4 = pNodes[ 3 ]->Z();

		double Det;

        Det = + X2*Y3*Z4 + X4*Y2*Z3 + X3*Y4*Z2 - X4*Y3*Z2 - X2*Y4*Z3 - X3*Y2*Z4
              - X1*Y3*Z4 - X4*Y1*Z3 - X3*Y4*Z1 + X4*Y3*Z1 + X1*Y4*Z3 + X3*Y1*Z4
              + X1*Y2*Z4 + X4*Y1*Z2 + X2*Y4*Z1 - X4*Y2*Z1 - X1*Y4*Z2 - X2*Y1*Z4
              - X1*Y2*Z3 - X3*Y1*Z2 - X2*Y3*Z1 + X3*Y2*Z1 + X1*Y3*Z2 + X2*Y1*Z3;

        return fabs( Det / 6.0 );
    }

	//***************************************************************************************************************************
    //* - Calculates the surface of a linear tetrahedra.
    //***************************************************************************************************************************
    double Modeler::AreaOfTriangle( Vector< Node::Pointer >& pNodes )
    {
		double na[ 3 ], v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = pNodes[ 2 ]->X() - pNodes[ 0 ]->X();
        v2[ 1 ] = pNodes[ 2 ]->Y() - pNodes[ 0 ]->Y();
        v2[ 2 ] = pNodes[ 2 ]->Z() - pNodes[ 0 ]->Z();
                                 
        v1[ 0 ] = pNodes[ 1 ]->X() - pNodes[ 0 ]->X();
        v1[ 1 ] = pNodes[ 1 ]->Y() - pNodes[ 0 ]->Y();
        v1[ 2 ] = pNodes[ 1 ]->Z() - pNodes[ 0 ]->Z();

        // Area = 0.5 * ( v2 x v1 )
        na[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        na[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        na[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

		return ( 0.5 * sqrt( na[ 0 ]*na[ 0 ] + na[ 1 ]*na[ 1 ] + na[ 2 ]*na[ 2 ] ) );
    }

	//***************************************************************************************************************************
    //* - Calculates the surface integral of fields.
    //***************************************************************************************************************************
    void Modeler::Generate_Srf_Intgrl_Element_FullWave( int* NodesId, unsigned int PropertiesId )
    {
	    // Element nodes
        Vector<Node::Pointer> pNodes;

        Get_Surface_Element_Nodes( NodesId, pNodes );

        pNodes.resize( 3 );

        // Surface ID
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		Vector<double> cIBC = (*pProperties)( COMPLEX_IBC );

		double SurfaceID = cIBC[ 0 ];

        // Add nodes to surface element list 
        if( ( mExportFields ) && ( mSrfIntg_Nodes.size() == 0 ) )
        {
            Vector<int> HONodesId; Push_HONodes_OnSurface( NodesId, HONodesId );

            HONodesId.resize( 3 );
            
            mSrfIntg_Elements[ SurfaceID ].push_back( HONodesId );           
        }

        // Surface integral
		double SurfcArea  = AreaOfTriangle( pNodes );
		double SrfIntg_Ni = SurfcArea / 3.0;

        // Field surface integrals vectors
        std::complex<double> cZero( 0.0, 0.0 );

		Vector< std::complex<double> > SrfIntg_E( 5, cZero );
        Vector< std::complex<double> > SrfIntg_H( 5, cZero );
        Vector< std::complex<double> > SrfIntg_B( 5, cZero );
        Vector< std::complex<double> > SrfIntg_J( 5, cZero );
        Vector< std::complex<double> > SrfIntg_S( 5, cZero );

        // External normal
        Vector< double > extN( 3 );

        double v1[ 3 ], v2[ 3 ];

        v2[ 0 ] = pNodes[ 2 ]->X() - pNodes[ 0 ]->X();
        v2[ 1 ] = pNodes[ 2 ]->Y() - pNodes[ 0 ]->Y();
        v2[ 2 ] = pNodes[ 2 ]->Z() - pNodes[ 0 ]->Z();

        v1[ 1 ] = pNodes[ 1 ]->Y() - pNodes[ 0 ]->Y();
        v1[ 2 ] = pNodes[ 1 ]->Z() - pNodes[ 0 ]->Z();
        v1[ 0 ] = pNodes[ 1 ]->X() - pNodes[ 0 ]->X();

        extN[ 0 ] = v2[ 1 ]*v1[ 2 ] - v2[ 2 ]*v1[ 1 ];
        extN[ 1 ] = v2[ 2 ]*v1[ 0 ] - v2[ 0 ]*v1[ 2 ];
        extN[ 2 ] = v2[ 0 ]*v1[ 1 ] - v2[ 1 ]*v1[ 0 ];

        extN *= ( 1.0 / sqrt( extN[ 0 ]*extN[ 0 ] + extN[ 1 ]*extN[ 1 ] + extN[ 2 ]*extN[ 2 ] ) );

        // Surface integral of fields
        for( int i=0; i<pNodes.size(); i++ )
        {
			// Electric field
            Vector<double> E_real = (*mpModel)( REAL_E, *pNodes[ i ] );
            Vector<double> E_imag = (*mpModel)( IMAG_E, *pNodes[ i ] );

            Add_NodalFieldIntegral( E_real, E_imag, SrfIntg_Ni, SrfIntg_E );

            // Magnetic field
			Vector<double> H_real = (*mpModel)( REAL_H, *pNodes[ i ] );
            Vector<double> H_imag = (*mpModel)( IMAG_H, *pNodes[ i ] );

            Add_NodalFieldIntegral( H_real, H_imag, SrfIntg_Ni, SrfIntg_H );
                                                        
            // Magnetic flux density
			Vector<double> B_real = (*mpModel)( REAL_B, *pNodes[ i ] );
            Vector<double> B_imag = (*mpModel)( IMAG_B, *pNodes[ i ] );
                                                        
            Add_NodalFieldIntegral( B_real, B_imag, SrfIntg_Ni, SrfIntg_B );

			// Current density
			Vector<double> J_real = (*mpModel)( REAL_J, *pNodes[ i ] );
            Vector<double> J_imag = (*mpModel)( IMAG_J, *pNodes[ i ] );

            Add_NodalFieldIntegral( J_real, J_imag, SrfIntg_Ni, SrfIntg_J );

            // Poynting vector S = E x conj(H)
            Vector< std::complex<double> > cE( 3 ), cH( 3 ), cJ( 3 );

            for( int i=0; i<3; i++ )
            {
                cE[ i ] = std::complex<double>( E_real[ i ], E_imag[ i ] );
                cH[ i ] = std::complex<double>( H_real[ i ], H_imag[ i ] );
                cJ[ i ] = std::complex<double>( J_real[ i ], J_imag[ i ] );
            }

            // S = E x conj(H)
            Vector< std::complex<double> > cS( 3 );

            cS[ 0 ] = cE[ 1 ]*std::conj( cH[ 2 ] ) - cE[ 2 ]*std::conj( cH[ 1 ] );
            cS[ 1 ] = cE[ 2 ]*std::conj( cH[ 0 ] ) - cE[ 0 ]*std::conj( cH[ 2 ] );
            cS[ 2 ] = cE[ 0 ]*std::conj( cH[ 1 ] ) - cE[ 1 ]*std::conj( cH[ 0 ] );

            std::complex<double> cSn = cS[ 0 ]*extN[ 0 ] + cS[ 1 ]*extN[ 1 ] + cS[ 2 ]*extN[ 2 ];
            std::complex<double> cJn = cJ[ 0 ]*extN[ 0 ] + cJ[ 1 ]*extN[ 1 ] + cJ[ 2 ]*extN[ 2 ];
            
            SrfIntg_S[ 0 ] += cS[ 0 ] * SrfIntg_Ni;
            SrfIntg_S[ 1 ] += cS[ 1 ] * SrfIntg_Ni;
            SrfIntg_S[ 2 ] += cS[ 2 ] * SrfIntg_Ni;
            SrfIntg_S[ 3 ] += cSn     * SrfIntg_Ni;
            SrfIntg_S[ 4 ] += cJn     * SrfIntg_Ni;
        }

		if( mSrfIntg_E.find( SurfaceID ) != mSrfIntg_E.end() )
		{
			mSurfce_m2[ SurfaceID ] += SurfcArea;
			mSrfIntg_E[ SurfaceID ] += SrfIntg_E;
			mSrfIntg_H[ SurfaceID ] += SrfIntg_H;
            mSrfIntg_B[ SurfaceID ] += SrfIntg_B;
			mSrfIntg_J[ SurfaceID ] += SrfIntg_J;
            mSrfIntg_S[ SurfaceID ] += SrfIntg_S;
		}
        else
		{
			mSurfce_m2[ SurfaceID ] = SurfcArea;
			mSrfIntg_E[ SurfaceID ] = SrfIntg_E;
			mSrfIntg_H[ SurfaceID ] = SrfIntg_H;
            mSrfIntg_B[ SurfaceID ] = SrfIntg_B;
			mSrfIntg_J[ SurfaceID ] = SrfIntg_J;
            mSrfIntg_S[ SurfaceID ] = SrfIntg_S;
		}
	}

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the volumetric element nodes.
    //**********************************************************************************************************************
    void Modeler::Get_Volume_Element_Nodes( int* NodesId, Vector< Node::Pointer >& pNodes )
    {
        Vector<int> HONodesId;

        Push_HONodes_OnVolume( NodesId, HONodesId );

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the volumetric element nodes and the high-order nodeIds vector.
    //**********************************************************************************************************************
    void Modeler::Get_Volume_Element_Nodes( int* NodesId, Vector<int>& HONodesId, Vector<Node::Pointer>& pNodes )
    {
        Push_HONodes_OnVolume( NodesId, HONodesId );

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the surface element nodes.
    //**********************************************************************************************************************
    void Modeler::Get_Surface_Element_Nodes( int* NodesId, Vector<Node::Pointer>& pNodes )
    {
        Vector<int> HONodesId;

        Push_HONodes_OnSurface( NodesId, HONodesId );

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }
    }

    //**********************************************************************************************************************
    //* - Returns a vector of pointers to the surface element nodes and the high-order nodeIds vector.
    //**********************************************************************************************************************
    void Modeler::Get_Surface_Element_Nodes( int* NodesId, Vector<int>& HONodesId, Vector<Node::Pointer>& pNodes )
    {
        Push_HONodes_OnSurface( NodesId, HONodesId );

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }
    }

    //**********************************************************************************************************************
    //* - Calculates volumetric integrals of the fields in electrostatic mode.
    //**********************************************************************************************************************
    void Modeler::Generate_Vol_Intgrl_Element_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // If PropertiesId = 0 then the element is used for surface integrals.
        if( PropertiesId == 0 )
        {
            GetGrad_forSrfIntg_SmoothingOff_Electrostatic( NodesId ); 

            return;
        }

        // Getting element nodes
        Vector<Node::Pointer> pNodes; 
        
        Get_Volume_Element_Nodes( NodesId, pNodes );

        // Getting Volume ID
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        Vector<double> cIBC = (*pProperties)( COMPLEX_IBC );

        double VolumeID = cIBC[ 0 ];

        // Integration Gauss points
        Vector<double> cX, cY, cZ, W;

        int numGaussPoints = GaussPoints3D_Order5( cX, cY, cZ, W );

        Matrix<double> N;

        Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // Volume of the element
        double EleVolume = VolumeOfTetrahedra( pNodes );
        double Jacob     = 6.0 * EleVolume;

        // Volumetric integrals vector
        Vector<double> VolIntg_E( 5, 0.0 );

        // Electrostatic electric field on nodes
        Vector< Vector<double> > EEfield_OnNodes( pNodes.size() );

        if( mResultsOnNodes )
        { 
            for( int i=0; i<pNodes.size(); i++ ) 
            {
                EEfield_OnNodes[ i ] = (*pProperties)( ELECTROSTATIC_ELECTRIC_FIELD, *pNodes[ i ] );
            }
        }
        else
        {
            Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );
            
            pElement->Calculate_Gradient_OnNodes( EEfield_OnNodes );  

            for( int i=0; i<pNodes.size(); i++ ) 
            {
                EEfield_OnNodes[ i ] *= -1;
            }
        }

        // Integration on volume over nodes
        for( int i=0; i<pNodes.size(); i++ )
        {
            double VolIntg_Ni = 0.0; 
            
            for( int gp=0; gp<numGaussPoints; gp++ )
            {
                VolIntg_Ni += Jacob * W[ gp ] * N[ i ][ gp ];
            }

            Vector<double> EEfield = EEfield_OnNodes[i];
            
            double modE2 = EEfield[ 0 ]*EEfield[ 0 ] + EEfield[ 1 ]*EEfield[ 1 ] + EEfield[ 2 ]*EEfield[ 2 ];
            double modE  = sqrt( modE2 );

            VolIntg_E[ 0 ] += EEfield[ 0 ] * VolIntg_Ni;
            VolIntg_E[ 1 ] += EEfield[ 1 ] * VolIntg_Ni;
            VolIntg_E[ 2 ] += EEfield[ 2 ] * VolIntg_Ni;
            VolIntg_E[ 3 ] += modE         * VolIntg_Ni;
            VolIntg_E[ 4 ] += modE2        * VolIntg_Ni;
        }

        // Adding volume integrals of this element to the total
        if( mVolume_m3.find( VolumeID ) != mVolume_m3.end() )
        {
            mVolume_m3              [ VolumeID ] += EleVolume;
            mVolIntg_E_Electrostatic[ VolumeID ] += VolIntg_E;
        }                                         
        else                                      
        {                                         
            mVolume_m3              [ VolumeID ] = EleVolume;
            mVolIntg_E_Electrostatic[ VolumeID ] = VolIntg_E;
        }
    }

    //**********************************************************************************************************************
    //* - Calculates gradients for the non-smoothed fields surface integrals.
    //*   If PropertiesId = 0 then the PVIE element is used for getting field gradients.
    //*   These gradients are used in the surface integrals of non-smoothed fields.
    //**********************************************************************************************************************
    void Modeler::GetGrad_forSrfIntg_SmoothingOff_Electrostatic( int* NodesId )
    {
        // Getting all element nodes
        Vector<Node::Pointer> nodes; 
        
        Get_Volume_Element_Nodes( NodesId, nodes );

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

    //**********************************************************************************************************************
    //* - Calculates surface integral of the fields for the electrostatic case.
    //**********************************************************************************************************************
    void Modeler::Generate_Srf_Intgrl_Element_Electrostatic( int* NodesId, unsigned int PropertiesId )
    {
        // Getting element nodes
        Vector< Node::Pointer > pNodes;

        Get_Surface_Element_Nodes( NodesId, pNodes );

        // Getting Surface ID
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        Vector<double> cIBC = (*pProperties)( COMPLEX_IBC );

        double SurfaceID = cIBC[ 0 ];

        // Gauss points
        Vector<double> cX, cY, W;

        int numGaussPoints = GaussPoints2D_Order12( cX, cY, W );

        // Linear base
        Matrix<double> N;

        Lagrange2D_Ni_1st( N, cX, cY );

        // Surface of the element
        double ElementArea = AreaOfTriangle( pNodes );
        double Jacob       = 2.0 * ElementArea;

        // Surface integrals vector
        Vector<double> SurfIntg_E( 5, 0.0 );

        // Electrostatic electric field on nodes 
        Vector< Vector<double> > EEfield_OnNodes( pNodes.size() );

        if( mResultsOnNodes )
        { 
            for( int i=0; i<pNodes.size(); i++ ) 
            {
                EEfield_OnNodes[ i ] = (*pProperties)( ELECTROSTATIC_ELECTRIC_FIELD, *pNodes[ i ] );
            }
        }
        else
        {
            for( int i=0; i<pNodes.size(); i++ ) 
            {
                EEfield_OnNodes[ i ] = mSurfInt_NonSmooth_E_Electrostatic[ mNonSmooth_Surf_Counter ];
            }
        }

        // Integration on volume over nodes
        for( int i=0; i<pNodes.size(); i++ )
        {
            double SurfIntg_Ni = 0.0;

            for( int gp=0; gp<numGaussPoints; gp++ ) 
            {
                SurfIntg_Ni += Jacob * W[ gp ] * N[ i ][ gp ];
            }

            Vector<double> EEfield = EEfield_OnNodes[ i ];

            double modE2 = EEfield[ 0 ]*EEfield[ 0 ] + EEfield[ 1 ]*EEfield[ 1 ] + EEfield[ 2 ]*EEfield[ 2 ];

            double modE = sqrt( modE2 );

            SurfIntg_E[ 0 ] += EEfield[ 0 ] * SurfIntg_Ni;
            SurfIntg_E[ 1 ] += EEfield[ 1 ] * SurfIntg_Ni;
            SurfIntg_E[ 2 ] += EEfield[ 2 ] * SurfIntg_Ni;
            SurfIntg_E[ 3 ] += modE         * SurfIntg_Ni;
            SurfIntg_E[ 4 ] += modE2        * SurfIntg_Ni;
        }                            

        // Adding surface integrals of this element to the total
        if( mSurfce_m2.find( SurfaceID ) != mSurfce_m2.end() )
        {
            mSurfce_m2              [ SurfaceID ] += ElementArea;
            mSrfIntg_E_Electrostatic[ SurfaceID ] += SurfIntg_E;
        }
        else
        {
            mSurfce_m2              [ SurfaceID ]  = ElementArea;
            mSrfIntg_E_Electrostatic[ SurfaceID ]  = SurfIntg_E;
        }

        // Increase global non-smoothing surface_id counter
        mNonSmooth_Surf_Counter++;
    }

	//**********************************************************************************************************************
    //* - Projects the field onto a rectangular waveguide port TE10 mode.  
    //**********************************************************************************************************************
    void Modeler::Generate_Projection_RectPort_TE10( int* NodesId, unsigned int PropertiesId )
    {
        Vector<Node::Pointer> pNodes;

        Get_Surface_Element_Nodes( NodesId, pNodes );

		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new RWPortTE10_3sb_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  1 ) pElement = Element::Pointer( new RWPortTE10_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new RWPortTE10_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new RWPortTE10_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new RWPortTE10_2eg_FullWave( pNodes, pProperties ) );

        pElement->SetPotentials( mPotentials );

        double PortNumber = pElement->GetPortNumber();

        // Elemental projection
        std::complex<double> TE10Proyec = pElement->GetProjection();

		if( mProjection.find( PortNumber ) != mProjection.end() ) 
        {
            mProjection[ PortNumber ] += TE10Proyec;
        }
        else     
        {
            mProjection[ PortNumber ] = TE10Proyec;
        }

        // Elemental normalization
        std::complex<double> Normalization = pElement->GetNormalization();

		if( mNormalization.find( PortNumber ) != mNormalization.end() )
        {
            mNormalization[ PortNumber ] += Normalization;
        }
        else   
        {
            mNormalization[ PortNumber ] = Normalization;
        }

        // Input ports 
		double Mod_Ein = (*pProperties)( COMPLEX_IBC_2o )[ 0 ];
	    double Pha_Ein = (*pProperties)( COMPLEX_IBC_2o )[ 1 ];
		
        if( ( Mod_Ein != 0.0 ) && ( mInputPorts.find( PortNumber ) == mInputPorts.end() ) )
        {
            mInputPorts[ PortNumber ] = std::complex<double>( Mod_Ein * cos( Pha_Ein ), Mod_Ein * sin( Pha_Ein ) );
        }
    }

	//**********************************************************************************************************************
    //* - Projects the field onto a coaxial waveguide port TEM mode. 
    //**********************************************************************************************************************
    void Modeler::Generate_Projection_CoaxPort_TEM( int* NodesId, unsigned int PropertiesId )
    {
		Vector<Node::Pointer> pNodes;

        Get_Surface_Element_Nodes( NodesId, pNodes );

		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new CoaxialPortTEM_3sb_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new CoaxialPortTEM_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new CoaxialPortTEM_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new CoaxialPortTEM_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new CoaxialPortTEM_2eg_FullWave( pNodes, pProperties ) );

        pElement->SetPotentials( mPotentials );
		
        double PortNumber = pElement->GetPortNumber();  

        // Elemental projection
        std::complex<double> TEMProyec = pElement->GetProjection();

		if( mProjection.find( PortNumber ) != mProjection.end() ) 
        {
            mProjection[ PortNumber ] += TEMProyec;
        }
        else      
        {
            mProjection[ PortNumber ] = TEMProyec;
        }

		// Elemental normalization
        std::complex<double> Normalization = pElement->GetNormalization();

		if( mNormalization.find( PortNumber ) != mNormalization.end() )
        {
            mNormalization[ PortNumber ] += Normalization;
        }
        else  
        {
            mNormalization[ PortNumber ] = Normalization;
        }

        // Input ports 
		double Mod_Ein = (*pProperties)( COMPLEX_IBC_2o )[ 0 ];
	    double Pha_Ein = (*pProperties)( COMPLEX_IBC_2o )[ 1 ];
		
        if( ( Mod_Ein != 0.0 ) && ( mInputPorts.find( PortNumber ) == mInputPorts.end() ) )
        {
            mInputPorts[ PortNumber ] = std::complex<double>( Mod_Ein * cos( Pha_Ein ), Mod_Ein * sin( Pha_Ein ) );
        }
    }

    //**********************************************************************************************************************
    //* - Calculates average normal in PEC nodes.
    //**********************************************************************************************************************
    void Modeler::Generate_PEC_Element( Vector<int>& LONodesId )
    {
        if( mElectrostaticMode )
        {
            return;
        }
 
        Vector<int> NodesId;

		Push_HONodes_OnSurface( LONodesId, NodesId );

		Vector<int>::const_iterator id_it;

        Vector<Node::Pointer> nodes;

        Node::Pointer pNode;
        
        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );

            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = ( nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->Z() - nodes[0]->Z() ) -
                    ( nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->Y() - nodes[0]->Y() );
                                                                                      
        normal[1] = ( nodes[2]->Z() - nodes[0]->Z()) * (nodes[1]->X() - nodes[0]->X() ) -
                    ( nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Z() - nodes[0]->Z() );
                                                                                      
        normal[2] = ( nodes[2]->X() - nodes[0]->X()) * (nodes[1]->Y() - nodes[0]->Y() ) -
                    ( nodes[2]->Y() - nodes[0]->Y()) * (nodes[1]->X() - nodes[0]->X() );

        // Pushing back the normal vector of the element to the PEC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetPECnormals[*id_it].push_back( normal );
		}
	}

    //**********************************************************************************************************************
    //* - Calculates normal in PMC nodes ( symmetry planes ).
    //**********************************************************************************************************************
    void Modeler::Generate_PMC_Element( Vector<int>& LONodesId )
    {
        if( mElectrostaticMode ) 
        {
            return;
        }

        if( mColdPlasmaMode )
        {
            if( !mEdgeElementsOn && !mpColdPlasma->Is_Epar_Tol_Off() )
            {
                Send_Error_Msg( "PMC boundary condition incompatible with active E par tolerance",
                                "Please, eliminate PMC or change E par tolerance to Off.", 1    );
            }
        }

	    Vector<int> NodesId;

		Push_HONodes_OnSurface( LONodesId, NodesId );

		Vector<int>::const_iterator id_it;

        Vector<Node::Pointer> nodes;

        Node::Pointer pNode;

        for( id_it=NodesId.begin(); id_it!=NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );

            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = ( nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z() ) -
                    ( nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y() );
                                                                                    
        normal[1] = ( nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X() ) -
                    ( nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z() );
                                                                                    
        normal[2] = ( nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y() ) -
                    ( nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X() );

        // Pushing back the normal vector of the element to the PMC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetPMCnormals[*id_it].push_back( normal );
		}
    }

	//**********************************************************************************************************************
    //* - Calculates average normal in TEC nodes ( symmetry plane ortogonal to a PMC plane ).
    //**********************************************************************************************************************
    void Modeler::Generate_TEC_Element( Vector<int>& LONodesId )
    {
        if( mElectrostaticMode ) 
        {
            return;
        }

        Vector<int> NodesId;

		Push_HONodes_OnSurface( LONodesId, NodesId );

		Vector<int>::const_iterator id_it;

        Vector<Node::Pointer> nodes;

        Node::Pointer pNode;

        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
        {
            pNode = mpModel->GetNode( (*id_it)-1 );

            nodes.push_back( pNode );
        }

        // Exterior normal of the element
        Vector<double> normal( 3 );

        // Normal = V2 x V1 = 2 * element_area * n
        normal[0] = ( nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->Z() - nodes[0]->Z() ) -
                    ( nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->Y() - nodes[0]->Y() );
                                                                                    
        normal[1] = ( nodes[2]->Z() - nodes[0]->Z())*(nodes[1]->X() - nodes[0]->X() ) -
                    ( nodes[2]->X() - nodes[0]->X())*(nodes[1]->Z() - nodes[0]->Z() );
                                                                                    
        normal[2] = ( nodes[2]->X() - nodes[0]->X())*(nodes[1]->Y() - nodes[0]->Y() ) -
                    ( nodes[2]->Y() - nodes[0]->Y())*(nodes[1]->X() - nodes[0]->X() );

        // Pushing back the normal vector of the element to the TEC normal set at each node
        for( id_it = NodesId.begin(); id_it != NodesId.end(); ++id_it )
		{
			mSetTECnormals[*id_it].push_back( normal );
		}
	}

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	// - Setting, fixing, solving, updating, steps.
    //
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //********************************************************************************************************************************
    //* - Finishes the first frequency sweep step.
    //********************************************************************************************************************************
	void Modeler::Finish_First_Step()
	{
		if( mProjection.size() > 0 )
        {
            Create_Directory( "Integrals/" );
            
            Calculate_Sij_Parameters();

	        Print_Sij_OnScreen( mProjectionNormalized, mInputPorts, mProblemFrequency );
            Write_Sij_Headers ( mProjectionNormalized, mInputPorts                    );    
            Write_Sij_OnFiles ( mProjectionNormalized, mInputPorts, mProblemFrequency );

		    Clear_Projection_Maps();
        }

		if( mVolIntg_E.size() > 0 )
        {            
            Create_Directory( "Integrals/" );  
            
            Print_VolIntgs_OnScreen( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F, mVolume_m3, mProblemFrequency );
            Write_VolIntgs_Headers ( mVolIntg_E );
            Write_VolIntgs_OnFiles ( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F,             mProblemFrequency );

            if ( mExportFields ) 
            {
                Initz_VolIntgs_Fields_OnFiles();
                Write_VolIntgs_Fields_OnFiles();
            }

            Clear_VolumeIntg_Maps();
        }

		if( mSrfIntg_E.size() > 0 )
        {
	        Create_Directory( "Integrals/" );
            
            Print_SrfIntgs_OnScreen( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S, mSurfce_m2, mProblemFrequency );  
            Write_SrfIntgs_Headers ( mSrfIntg_E );     
            Write_SrfIntgs_OnFiles ( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S,             mProblemFrequency );

            if( mExportFields ) 
            {
                Initz_SrfIntgs_Fields_OnFiles();
                Write_SrfIntgs_Fields_OnFiles();
            }

            Clear_SurfacIntg_Maps();
        }

        b_vector.resize( mSystemSize, 0.00 );

        A_matrix.Resize( mSystemSize );

        if( mA_matrix_aux_Required )
        {
            A_matrix_aux.Resize( mSystemSize );
        }
	}

    //********************************************************************************************************************************
    //* - Creates directory. 
    //********************************************************************************************************************************
	void Modeler::Create_Directory( String Directory_Name )
    {
        #ifdef _WIN64
            mkdir( Directory_Name );
        #else
            mkdir( Directory_Name, 0755 );
        #endif    
    }

    //********************************************************************************************************************************
    //* - Prepares the next step in frequency sweep mode.
    //********************************************************************************************************************************
	void Modeler::Finish_Step()
	{
		if( mProjection.size() > 0 )
        {
            Calculate_Sij_Parameters();

	        Print_Sij_OnScreen( mProjectionNormalized, mInputPorts, mProblemFrequency );
            Write_Sij_OnFiles ( mProjectionNormalized, mInputPorts, mProblemFrequency );

		    Clear_Projection_Maps();
        }

		if( mVolIntg_E.size() > 0 )
        {           
            Print_VolIntgs_OnScreen( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F, mVolume_m3, mProblemFrequency );
            Write_VolIntgs_OnFiles ( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F,             mProblemFrequency );

            if( mExportFields ) 
            {
                Write_VolIntgs_Fields_OnFiles();
            }

            Clear_VolumeIntg_Maps();
        }

		if( mSrfIntg_E.size() > 0 )
        {
	        Print_SrfIntgs_OnScreen( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S, mSurfce_m2, mProblemFrequency );   
            Write_SrfIntgs_OnFiles ( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S,             mProblemFrequency );

            if( mExportFields ) 
            {
                Write_SrfIntgs_Fields_OnFiles();
            }
            
            Clear_SurfacIntg_Maps();
        }

        b_vector.resize( mSystemSize, 0.00 );

        A_matrix.Resize( mSystemSize );

        if( mA_matrix_aux_Required )
        {
            A_matrix_aux.Resize( mSystemSize );
        }
	}

    //********************************************************************************************************************************
    //* - Cleans electrostatic maps. 
    //********************************************************************************************************************************
    void Modeler::Clear_ElecStatic_Maps()
    {
        mVolume_m3.clear(); std::map<double, double>().swap( mVolume_m3 );
        mSurfce_m2.clear(); std::map<double, double>().swap( mSurfce_m2 );

        mVolIntg_E_Electrostatic.clear(); std::map<double, Vector<double> >().swap( mVolIntg_E_Electrostatic );
        mSrfIntg_E_Electrostatic.clear(); std::map<double, Vector<double> >().swap( mSrfIntg_E_Electrostatic );
    }

    //********************************************************************************************************************************
    //* - Cleans projection maps. 
    //********************************************************************************************************************************
    void Modeler::Clear_Projection_Maps()
    {
        mProjection          .clear(); std::map<double, std::complex<double> >().swap( mProjection           );
        mInputPorts          .clear(); std::map<double, std::complex<double> >().swap( mInputPorts           );
        mNormalization       .clear(); std::map<double, std::complex<double> >().swap( mNormalization        );
        mProjectionNormalized.clear(); std::map<double, std::complex<double> >().swap( mProjectionNormalized );
    }

    //********************************************************************************************************************************
    //* - Cleans volume integral maps. 
    //********************************************************************************************************************************
    void Modeler::Clear_VolumeIntg_Maps()
    {
        mVolume_m3.clear(); std::map<double, double>                         ().swap( mVolume_m3 );
        mVolIntg_E.clear(); std::map<double, Vector< std::complex<double> > >().swap( mVolIntg_E );                                                                    
        mVolIntg_H.clear(); std::map<double, Vector< std::complex<double> > >().swap( mVolIntg_H );
        mVolIntg_B.clear(); std::map<double, Vector< std::complex<double> > >().swap( mVolIntg_B );                                                                   
        mVolIntg_J.clear(); std::map<double, Vector< std::complex<double> > >().swap( mVolIntg_J );
        mVolIntg_F.clear(); std::map<double, Vector< std::complex<double> > >().swap( mVolIntg_F );
    }

    //********************************************************************************************************************************
    //* - Cleans surface integral maps. 
    //********************************************************************************************************************************
    void Modeler::Clear_SurfacIntg_Maps()
    {	    
        mSurfce_m2.clear(); std::map<double, double>                         ().swap( mSurfce_m2 );
	    mSrfIntg_E.clear(); std::map<double, Vector< std::complex<double> > >().swap( mSrfIntg_E );                                                              
	    mSrfIntg_H.clear(); std::map<double, Vector< std::complex<double> > >().swap( mSrfIntg_H );
        mSrfIntg_B.clear(); std::map<double, Vector< std::complex<double> > >().swap( mSrfIntg_B );                                                                   
	    mSrfIntg_J.clear(); std::map<double, Vector< std::complex<double> > >().swap( mSrfIntg_J );
        mSrfIntg_S.clear(); std::map<double, Vector< std::complex<double> > >().swap( mSrfIntg_S );
    }

    //********************************************************************************************************************************
    //* - Initializes files to write volume integral fields. 
    //********************************************************************************************************************************
    void Modeler::Initz_VolIntgs_Fields_OnFiles()
    {	
        Create_Directory( "Fields/" );

        std::cout << "Writing volumes meshes on files..." << std::endl;

        Write_VolIntgs_Elements_OnFiles( mVolIntg_Elements );

        std::map<int, Vector< Vector<int> > >::iterator vit;

		for ( vit = mVolIntg_Elements.begin(); vit != mVolIntg_Elements.end(); vit++ )
        {
            std::list<int> NodesListInVol;
            
            Vector< Vector<int> > ElementsInVol = vit->second;

            Vector< Vector<int> >::iterator eit;

            for ( eit = ElementsInVol.begin(); eit != ElementsInVol.end(); eit++ )
            {
                Vector<int> Element = *eit;

                for ( int i=0; i<Element.size(); i++ ) 
                {
                    NodesListInVol.push_back( Element[ i ] ); 
                }
            }

            NodesListInVol.sort();

            NodesListInVol.unique();

            Vector<int> NodesIdsInVol;

            Vector< Vector<double> > NodesCoordsInVol;

            std::list<int>::iterator nit;

            for ( nit = NodesListInVol.begin(); nit != NodesListInVol.end(); nit++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit) - 1 );

                Vector<double> NodeCoords( 3 );

                NodeCoords[ 0 ] = pNode->X();
                NodeCoords[ 1 ] = pNode->Y();
                NodeCoords[ 2 ] = pNode->Z();

                NodesCoordsInVol.push_back( NodeCoords );
                NodesIdsInVol   .push_back( *nit       );
            }

            Write_VolIntgs_Nodes_OnFile( NodesIdsInVol, NodesCoordsInVol, vit->first );

            // Add volume nodes list to global volumes nodes lists 
            mVolIntg_Nodes.push_back( NodesIdsInVol );

            // Cleaning previous fields files
            String FieldsFileBaseName; 
            
            FieldsFileBaseName << "Fields/Volume-" << vit->first;
        
            std::fstream FieldsFile_E( FieldsFileBaseName + "-E.dat", std::ios::out );
            std::fstream FieldsFile_H( FieldsFileBaseName + "-H.dat", std::ios::out );
            std::fstream FieldsFile_B( FieldsFileBaseName + "-B.dat", std::ios::out );
            std::fstream FieldsFile_J( FieldsFileBaseName + "-J.dat", std::ios::out );

            FieldsFile_E.close();
            FieldsFile_H.close();
            FieldsFile_B.close();
            FieldsFile_J.close();
        }

        std::cout << "Done." << std::endl << std::endl;

        mVolIntg_Elements.clear(); std::map<int, Vector< Vector<int> > >().swap( mVolIntg_Elements );
    }

    //********************************************************************************************************************************
    //* - Writes volume integral fields on files.
    //********************************************************************************************************************************
    void Modeler::Write_VolIntgs_Fields_OnFiles()
    {	
        std::cout << "Writing volumes fields on files..." << std::endl;

        std::map<double, double>::iterator vit = mVolume_m3.begin();

        Vector< Vector<int> >::iterator nit;

        for ( nit = mVolIntg_Nodes.begin(); nit != mVolIntg_Nodes.end(); nit++, vit++ )
        {
            // Vector with the fields on the volume nodes
            Vector< Vector< std::complex<double> > > cFvn( (*nit).size() );

            // Field on node
            Vector< std::complex<double> > cFn( 3 );

            // Field name
            String FieldName;

            // Electric field
            FieldName = "E";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> E_real = (*mpModel)( REAL_E, *pNode );
                Vector<double> E_imag = (*mpModel)( IMAG_E, *pNode );   

                cFn[ 0 ] = std::complex<double>( E_real[ 0 ], E_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( E_real[ 1 ], E_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( E_real[ 2 ], E_imag[ 2 ] );

                cFvn[ i ] = cFn;
            } 

            if( mComplxFreqMode ) 
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, vit->first );
            }
            else
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, vit->first );
            }

            // Magnetic field
            FieldName = "H";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> H_real = (*mpModel)( REAL_H, *pNode );
                Vector<double> H_imag = (*mpModel)( IMAG_H, *pNode );   

                cFn[ 0 ] = std::complex<double>( H_real[ 0 ], H_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( H_real[ 1 ], H_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( H_real[ 2 ], H_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, vit->first );
            }
            else
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, vit->first );
            }

            // Magnetic flux density
            FieldName = "B";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> B_real = (*mpModel)( REAL_B, *pNode );
                Vector<double> B_imag = (*mpModel)( IMAG_B, *pNode );   

                cFn[ 0 ] = std::complex<double>( B_real[ 0 ], B_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( B_real[ 1 ], B_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( B_real[ 2 ], B_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, vit->first );
            }
            else
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, vit->first );
            }

            // Current density
            FieldName = "J";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> J_real = (*mpModel)( REAL_J, *pNode );
                Vector<double> J_imag = (*mpModel)( IMAG_J, *pNode );   

                cFn[ 0 ] = std::complex<double>( J_real[ 0 ], J_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( J_real[ 1 ], J_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( J_real[ 2 ], J_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, vit->first );
            }
            else
            {
                Write_VolIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, vit->first );
            }
        }

        std::cout << "Done." << std::endl << std::endl;
    }

    //********************************************************************************************************************************
    //* - Initializes files to write surface integral fields.  
    //********************************************************************************************************************************
    void Modeler::Initz_SrfIntgs_Fields_OnFiles()
    {	
        Create_Directory( "Fields/" );

        std::cout << "Writing surfaces meshes on files..." << std::endl;

        Write_SrfIntgs_Elements_OnFiles( mSrfIntg_Elements );

        std::map<int, Vector< Vector<int> > >::iterator sit;

		for ( sit = mSrfIntg_Elements.begin(); sit != mSrfIntg_Elements.end(); sit++ )
        {
            std::list<int> NodesListInSrf;
            
            Vector< Vector<int> > ElementsInSrf = sit->second;

            Vector< Vector<int> >::iterator eit;

            for ( eit = ElementsInSrf.begin(); eit != ElementsInSrf.end(); eit++ )
            {
                Vector<int> Element = *eit;

                for ( int i=0; i<Element.size(); i++ ) 
                {
                    NodesListInSrf.push_back( Element[ i ] ); 
                }
            }

            NodesListInSrf.sort();

            NodesListInSrf.unique();

            Vector<int> NodesIdsInSrf;

            Vector< Vector<double> > NodesCoordsInSrf;

            std::list<int>::iterator nit;

            for ( nit = NodesListInSrf.begin(); nit != NodesListInSrf.end(); nit++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit) - 1 );

                Vector<double> NodeCoords( 3 );

                NodeCoords[ 0 ] = pNode->X();
                NodeCoords[ 1 ] = pNode->Y();
                NodeCoords[ 2 ] = pNode->Z();

                NodesCoordsInSrf.push_back( NodeCoords );
                NodesIdsInSrf   .push_back( *nit       );
            }

            Write_SrfIntgs_Nodes_OnFile( NodesIdsInSrf, NodesCoordsInSrf, sit->first );

            // Add surface nodes list to global surface nodes lists 
            mSrfIntg_Nodes.push_back( NodesIdsInSrf );

            // Cleaning previous fields files
            String FieldsFileBaseName; 
            
            FieldsFileBaseName << "Fields/Surface-" << sit->first;
        
            std::fstream FieldsFile_E( FieldsFileBaseName + "-E.dat", std::ios::out ); 
            std::fstream FieldsFile_H( FieldsFileBaseName + "-H.dat", std::ios::out ); 
            std::fstream FieldsFile_B( FieldsFileBaseName + "-B.dat", std::ios::out ); 
            std::fstream FieldsFile_J( FieldsFileBaseName + "-J.dat", std::ios::out ); 

            FieldsFile_E.close();
            FieldsFile_H.close();
            FieldsFile_B.close();
            FieldsFile_J.close();
        }

        std::cout << "Done." << std::endl << std::endl;

        mSrfIntg_Elements.clear(); std::map<int, Vector< Vector<int> > >().swap( mSrfIntg_Elements );
    }

    //********************************************************************************************************************************
    //* - Writes surface integral fields on files. 
    //********************************************************************************************************************************
    void Modeler::Write_SrfIntgs_Fields_OnFiles()
    {	
        std::cout << "Writing surfaces fields on files..." << std::endl;

        std::map<double, double>::iterator sit = mSurfce_m2.begin();

        Vector< Vector<int> >::iterator nit;

        for ( nit = mSrfIntg_Nodes.begin(); nit != mSrfIntg_Nodes.end(); nit++, sit++ )
        {
            // Vector with the fields on the volume nodes
            Vector< Vector< std::complex<double> > > cFvn( (*nit).size() );

            // Field on node
            Vector< std::complex<double> > cFn( 3 );

            // Field name
            String FieldName;

            // Electric field
            FieldName = "E";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> E_real = (*mpModel)( REAL_E, *pNode );
                Vector<double> E_imag = (*mpModel)( IMAG_E, *pNode );   

                cFn[ 0 ] = std::complex<double>( E_real[ 0 ], E_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( E_real[ 1 ], E_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( E_real[ 2 ], E_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, sit->first );
            }
            else
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, sit->first );
            }
        
            // Magnetic field
            FieldName = "H";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> H_real = (*mpModel)( REAL_H, *pNode );
                Vector<double> H_imag = (*mpModel)( IMAG_H, *pNode );   

                cFn[ 0 ] = std::complex<double>( H_real[ 0 ], H_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( H_real[ 1 ], H_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( H_real[ 2 ], H_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, sit->first );
            }
            else
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, sit->first );
            }

            // Magnetic flux density
            FieldName = "B";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> B_real = (*mpModel)( REAL_B, *pNode );
                Vector<double> B_imag = (*mpModel)( IMAG_B, *pNode );   

                cFn[ 0 ] = std::complex<double>( B_real[ 0 ], B_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( B_real[ 1 ], B_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( B_real[ 2 ], B_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, sit->first );
            }
            else
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, sit->first );
            }

            // Current density
            FieldName = "J";

            for ( int i=0; i<(*nit).size(); i++ )
            {
                Node::Pointer pNode = mpModel->GetNode( (*nit)[ i ] - 1 );

                Vector<double> J_real = (*mpModel)( REAL_J, *pNode );
                Vector<double> J_imag = (*mpModel)( IMAG_J, *pNode );   

                cFn[ 0 ] = std::complex<double>( J_real[ 0 ], J_imag[ 0 ] );
                cFn[ 1 ] = std::complex<double>( J_real[ 1 ], J_imag[ 1 ] );
                cFn[ 2 ] = std::complex<double>( J_real[ 2 ], J_imag[ 2 ] );

                cFvn[ i ] = cFn;
            }

            if( mComplxFreqMode ) 
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mComplexFrequency, sit->first );
            }
            else
            {
                Write_SrfIntgs_Fields_OnFile( FieldName, cFvn, mProblemFrequency, sit->first );
            }
        }

        std::cout << "Done." << std::endl << std::endl;
    }

    //***********************************************************************************************************************************************
    //* - Solves the linear problem in full wave mode.
    //***********************************************************************************************************************************************
    void Modeler::Solve_Problem_FullWave()
    {
        // Solving linear system Ax=b
        Solve_Linear_System();

        // Saving solution vector on Variables 
		if( mPotentials )
        {
            Save_System_Solution_FullWave( cAx );
            Save_System_Solution_FullWave( cAy );
            Save_System_Solution_FullWave( cAz );    
            Save_System_Solution_FullWave( cVs );    
        }
        else
        {
            Save_System_Solution_FullWave( cEx );
            Save_System_Solution_FullWave( cEy );
            Save_System_Solution_FullWave( cEz );        
        }

        // Updating fixed DOFs
        if( mPotentials     ) Update_FullWv_Voltage();
        if( mColdPlasmaMode ) UpdateCoord_PlasmaRLP();
        if( mAxisymmetric   ) UpdateAxisToCartesian();

		UpdateCoord_EPEC   ();
        UpdateCoord_DisctAV();
		UpdateCoord_Contact();
		UpdateCoord_EPBC   ();
    }

    //***********************************************************************************************************************************************
    //* - Solves the linear problem in electrostatic mode.
    //***********************************************************************************************************************************************
    void Modeler::Solve_Problem_Electrostatic()
    {
        // Solving linear system Ax=b
        Solve_Linear_System();

        // Saving solution vector on VOLTAGE variable
        Save_System_Solution_Electrostatic( VOLTAGE );

        // Updating fixed DOFs
        Update_Static_Voltage();
    }

    //***********************************************************************************************************************************************
    //* - Updates fixed electrostatic voltage DOFs.
    //***********************************************************************************************************************************************
    void Modeler::Update_Static_Voltage()
    {
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for( nit=nodes.begin(); nit!=nodes.end(); nit++ )
        {
            if( (*nit)->pDofV()->IsFixed() )
            {
                (*mpModel)( VOLTAGE, **nit ) = mFix_Static_Voltage[ (*nit)->Id() ];
            }
        }
    }

    //***********************************************************************************************************************************************
    //* - Updates fixed full wave voltage DOFs.
    //***********************************************************************************************************************************************
    void Modeler::Update_FullWv_Voltage()
    {
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for( nit=nodes.begin(); nit!=nodes.end(); nit++ )
        {
            if( mFix_FullWv_Voltage.find( (*nit)->Id() ) != mFix_FullWv_Voltage.end() )
            {
                (*mpModel)( cVs, **nit ) = mFix_FullWv_Voltage[ (*nit)->Id() ];
            }
        }
    }

    //***********************************************************************************************************************************************
    //* - Updates coordinate system in plasma mode.
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_PlasmaRLP()
    {
        if( mEdgeElementsOn || mpColdPlasma->Is_Epar_Tol_Off() )
        {
            return;
        }

        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nodes_it;

        for( nodes_it=nodes.begin(); nodes_it!= nodes.end(); ++nodes_it )
        {
            unsigned int nodeId = (*nodes_it)->Id();

            if( mNormals.find( nodeId ) == mNormals.end() )
            {
                // Complex constants
                std::complex<double> cZero    ( 0.0        , 0.0         );
                std::complex<double> cUnit    ( 0.0        , 1.0         );
                std::complex<double> rUnit    ( 1.0        , 0.0         );
                std::complex<double> rInvSqrt2( 1.0/sqrt(2), 0.0         );
                std::complex<double> cInvSqrt2( 0.0        , 1.0/sqrt(2) );

                // E_RLP = [U] * E_SDP
                Matrix<std::complex<double> > U_Herm( 3, 3 );

                U_Herm[0][0] = rInvSqrt2; U_Herm[0][1] = rInvSqrt2; U_Herm[0][2] = cZero;
                U_Herm[1][0] =-cInvSqrt2; U_Herm[1][1] = cInvSqrt2; U_Herm[1][2] = cZero;
                U_Herm[2][0] = cZero    ; U_Herm[2][1] = cZero    ; U_Herm[2][2] = rUnit;

                // E_SDP = [T] * E_CC
                Matrix<std::complex<double> > T( 3, 3 );

                Rotation_Matrix_SDP_T_CC( T, nodeId );

                // [R] = [Trans(T)][Herm(U)]
                Matrix<std::complex<double> > R( 3, 3 );

                R = T.Transpose() * U_Herm;

                std::complex<double> Er = (*mpModel)( cEx, **nodes_it );
                std::complex<double> El = (*mpModel)( cEy, **nodes_it );
                std::complex<double> Ep = (*mpModel)( cEz, **nodes_it );

                (*mpModel)( cEx, **nodes_it ) = R[0][0]*Er + R[0][1]*El + R[0][2]*Ep;
                (*mpModel)( cEy, **nodes_it ) = R[1][0]*Er + R[1][1]*El + R[1][2]*Ep;
                (*mpModel)( cEz, **nodes_it ) = R[2][0]*Er + R[2][1]*El + R[2][2]*Ep;
            }
        }
    }

	//***********************************************************************************************************************************************
    //* - Updates coordinate system in axisymmetric problems.
    //***********************************************************************************************************************************************
    void Modeler::UpdateAxisToCartesian()
    {
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;
		
        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			if( ( mNormals.find ( (*nit)->Id() ) == mNormals.end() )  && 
                ( DistanceToAxis( (*nit)->Id() )  > mGeoTolerance  )   )
            {
                Vector<double> n( 3 );
                Vector<double> t( 3 );
                Vector<double> b( 3 );

				AxisymmetricCoordinates( (*nit)->Id(), n, t, b );

                if( mPotentials )
                {
                    std::complex<double> An = (*mpModel)( cAx, **nit );
				    std::complex<double> At = (*mpModel)( cAy, **nit );
                    std::complex<double> Ab = (*mpModel)( cAz, **nit );
                    
                    (*mpModel)( cAx, **nit ) = n[0]*An + t[0]*At + b[0]*Ab;
                    (*mpModel)( cAy, **nit ) = n[1]*An + t[1]*At + b[1]*Ab;
                    (*mpModel)( cAz, **nit ) = n[2]*An + t[2]*At + b[2]*Ab;                
                }
                else
                {
                    std::complex<double> En = (*mpModel)( cEx, **nit );
				    std::complex<double> Et = (*mpModel)( cEy, **nit );
                    std::complex<double> Eb = (*mpModel)( cEz, **nit );
                    
                    (*mpModel)( cEx, **nit ) = n[0]*En + t[0]*Et + b[0]*Eb;
                    (*mpModel)( cEy, **nit ) = n[1]*En + t[1]*Et + b[1]*Eb;
                    (*mpModel)( cEz, **nit ) = n[2]*En + t[2]*Et + b[2]*Eb;
                }
            }
        }
    }

	//***********************************************************************************************************************************************
    //* - Changes local coordinates in PEC for cartesian coordinates. 
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_EPEC()
    {
        if( ( mNormals.size() == 0 ) || ( mSetAllFxyToZero ) || ( mEdgeElementsOn ) ) 
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
            int NodeId = (*nit)->Id();
            
            if( ( Is_DisctAVNode( NodeId ) ) || ( mNormals.find( NodeId ) == mNormals.end() ) ) 
            {
                continue;
            }

            // From E_ntb | A_ntb to E_xyz | A_xyz ( PEC, PMC BC )
            Vector<double> n( mNormals[ NodeId ] );
			Vector<double> t( 3 ); 
            Vector<double> b( 3 );

            TangencialCoordinates( n, t, b );

            if ( mPotentials )
            {
                std::complex<double> An = (*mpModel)( cAx, **nit );
			    std::complex<double> At = (*mpModel)( cAy, **nit );
                std::complex<double> Ab = (*mpModel)( cAz, **nit );
                
                // {A_xyz} = [Rt] * {A_ntb}
                (*mpModel)( cAx, **nit ) = n[0]*An + t[0]*At + b[0]*Ab;
                (*mpModel)( cAy, **nit ) = n[1]*An + t[1]*At + b[1]*Ab;
                (*mpModel)( cAz, **nit ) = n[2]*An + t[2]*At + b[2]*Ab;            
            }
            else
            {
                std::complex<double> En = (*mpModel)( cEx, **nit );
			    std::complex<double> Et = (*mpModel)( cEy, **nit );
                std::complex<double> Eb = (*mpModel)( cEz, **nit );
                
                // {E_xyz} = [Rt] * {E_ntb}
                (*mpModel)( cEx, **nit ) = n[0]*En + t[0]*Et + b[0]*Eb;
                (*mpModel)( cEy, **nit ) = n[1]*En + t[1]*Et + b[1]*Eb;
                (*mpModel)( cEz, **nit ) = n[2]*En + t[2]*Et + b[2]*Eb;            
            }
        }
    }

	//***********************************************************************************************************************************************
    //* - Updates electric field values in PBC nodes.
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_EPBC()
	{
		Vector<int> FrontLeftNodesId;

		UpdateCoord_RL_PBC( FrontLeftNodesId );
		UpdateCoord_FB_PBC(                  );
		UpdateCoord_LF_PBC( FrontLeftNodesId );
	}

	//***********************************************************************************************************************************************
    //* - Checks if the left node is also a front node.
    //***********************************************************************************************************************************************
	bool Modeler::PBC_LeftNodeIsFrontNode( int rightNodeId, Vector<int>& leftNodesId, Vector<int>& FrontLeftNodesId )
	{
		bool LeftNodeIsFrontNode = false;

		for( int i=0; i<leftNodesId.size(); i++ )
		{
			if( mPBC_NodeNodePairs_Front.find( leftNodesId[ i ] ) != mPBC_NodeNodePairs_Front.end() ) 
			{
				LeftNodeIsFrontNode = true;

				FrontLeftNodesId.push_back( rightNodeId );

				break;
		    }
		}

		return LeftNodeIsFrontNode;
	}

	//***********************************************************************************************************************************************
    //* - Updates electric field values in PBC nodes ( Right-Left faces ).
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_RL_PBC( Vector<int>& FrontLeftNodesId )
	{
		if( mPBC_NodeNodePairs_Right.size() == 0 )
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

		for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			if( mPBC_NodeNodePairs_Right.find( (*nit)->Id() ) != mPBC_NodeNodePairs_Right.end() )
            {
				int lNodeId = mPBC_NodeNodePairs_Right[ (*nit)->Id() ];

				Vector<int> leftNodesId( 1 ); leftNodesId[0] = lNodeId;

				if( PBC_LeftNodeIsFrontNode( (*nit)->Id(), leftNodesId, FrontLeftNodesId ) ) 
                {
                    continue;
                }

		        Node::Pointer plNode = mpModel->GetNode( lNodeId-1 );

                std::complex<double> lcFx( 0.00, 0.00 );
			    std::complex<double> lcFy( 0.00, 0.00 );
			    std::complex<double> lcFz( 0.00, 0.00 );

                if( mPotentials )
                {
				    lcFx = mExp_jPha_RL * (*mpModel)( cAx, *plNode );
				    lcFy = mExp_jPha_RL * (*mpModel)( cAy, *plNode );
				    lcFz = mExp_jPha_RL * (*mpModel)( cAz, *plNode ); 

                    if( plNode->pDof( cVs ) != NULL ) 
                    {
                        (*mpModel)( cVs, **nit ) = mExp_jPha_RL * (*mpModel)( cVs, *plNode ); 
                    }
                }
                else
                {
				    lcFx = mExp_jPha_RL * (*mpModel)( cEx, *plNode );
				    lcFy = mExp_jPha_RL * (*mpModel)( cEy, *plNode );
				    lcFz = mExp_jPha_RL * (*mpModel)( cEz, *plNode );                
                }

                double cos_A = 1.0;
			    double sin_A = 0.0;
                
			    if( ( mIsRightPBCTilted == true ) && ( mEdgeElementsOn == false ) ) 
                {
                    PBC_InvRot( (*nit)->Id(), cos_A, sin_A );
                }
                
                if( mPotentials )
                {
			        (*mpModel)( cAx, **nit ) =  cos_A * lcFx + sin_A * lcFy;
			        (*mpModel)( cAy, **nit ) = -sin_A * lcFx + cos_A * lcFy;
			        (*mpModel)( cAz, **nit ) =  lcFz; 
                }
                else
                {
			        (*mpModel)( cEx, **nit ) =  cos_A * lcFx + sin_A * lcFy;
			        (*mpModel)( cEy, **nit ) = -sin_A * lcFx + cos_A * lcFy;
			        (*mpModel)( cEz, **nit ) =  lcFz;                
                }
            }
		}
	}

	//***********************************************************************************************************************************************
    //* - Computes inverse rotation of a point.
    //***********************************************************************************************************************************************
    void Modeler::PBC_InvRot( int NodeId, double& cos_A, double& sin_A )
    {
        double X = ( mpModel->GetNode( NodeId - 1 ) )->X();
		double Y = ( mpModel->GetNode( NodeId - 1 ) )->Y();

		double mod_r = sqrt(X*X + Y*Y);

		cos_A = Y / mod_r;
		sin_A = X / mod_r;
	}

	//***********************************************************************************************************************************************
    //* - Updates field values in PBC nodes ( Front-Back surfaces ).
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_FB_PBC()
	{
		if( mPBC_NodeNodePairs_Front.size() == 0 )
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator nit;

		for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			if( mPBC_NodeNodePairs_Front.find( (*nit)->Id() ) != mPBC_NodeNodePairs_Front.end() )
            {
				int bNodeId = mPBC_NodeNodePairs_Front[ (*nit)->Id() ];

		        Node::Pointer pbNode = mpModel->GetNode( bNodeId-1 );

                if( mPotentials )
                {
				    (*mpModel)( cAx, **nit ) = mExp_jPha_FB * (*mpModel)( cAx, *pbNode );
				    (*mpModel)( cAy, **nit ) = mExp_jPha_FB * (*mpModel)( cAy, *pbNode );
				    (*mpModel)( cAz, **nit ) = mExp_jPha_FB * (*mpModel)( cAz, *pbNode ); 

                    if( pbNode->pDof( cVs ) != NULL ) 
                    {
                        (*mpModel)( cVs, **nit ) = mExp_jPha_FB * (*mpModel)( cVs, *pbNode ); 
                    }
                }
                else
                {
				    (*mpModel)( cEx, **nit ) = mExp_jPha_FB * (*mpModel)( cEx, *pbNode );
				    (*mpModel)( cEy, **nit ) = mExp_jPha_FB * (*mpModel)( cEy, *pbNode );
				    (*mpModel)( cEz, **nit ) = mExp_jPha_FB * (*mpModel)( cEz, *pbNode );                
                }
            }
        }
	}

	//***********************************************************************************************************************************************
    //* - Updates electric field values in PBC nodes ( Left-Front nodes ).
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_LF_PBC( Vector<int>& FrontLeftNodes )
	{
		if( FrontLeftNodes.size() == 0 ) 
        { 
            return; 
        }

		Vector<int>::iterator it_fln;

		for( it_fln = FrontLeftNodes.begin(); it_fln != FrontLeftNodes.end(); ++it_fln )
        {
			if( mPBC_NodeNodePairs_Right.find( (*it_fln) ) != mPBC_NodeNodePairs_Right.end() )
            {
                int lNode = mPBC_NodeNodePairs_Right[ (*it_fln) ];

				Node::Pointer plNode = mpModel->GetNode(   lNode  -1 );
                Node::Pointer prNode = mpModel->GetNode( (*it_fln)-1 );

                std::complex<double> lcFx( 0.00, 0.00 );
			    std::complex<double> lcFy( 0.00, 0.00 );
			    std::complex<double> lcFz( 0.00, 0.00 );

                if( mPotentials )
                {
				    lcFx = mExp_jPha_RL * (*mpModel)( cAx, *plNode );
				    lcFy = mExp_jPha_RL * (*mpModel)( cAy, *plNode );
				    lcFz = mExp_jPha_RL * (*mpModel)( cAz, *plNode ); 

                    if( plNode->pDof( cVs ) != NULL ) 
                    {
                        (*mpModel)( cVs, *prNode ) = mExp_jPha_RL * (*mpModel)( cVs, *plNode ); 
                    }
                }
                else
                {
				    lcFx = mExp_jPha_RL * (*mpModel)( cEx, *plNode );
				    lcFy = mExp_jPha_RL * (*mpModel)( cEy, *plNode );
				    lcFz = mExp_jPha_RL * (*mpModel)( cEz, *plNode );                
                }

			    double cos_A = 1.0;
			    double sin_A = 0.0;
                
			    if( ( mIsRightPBCTilted == true ) && ( mEdgeElementsOn == false ) ) 
                {
                    PBC_InvRot( (*it_fln), cos_A, sin_A );
                }

                if( mPotentials )
                {
			        (*mpModel)( cAx, *prNode ) =  cos_A * lcFx + sin_A * lcFy;
			        (*mpModel)( cAy, *prNode ) = -sin_A * lcFx + cos_A * lcFy;
			        (*mpModel)( cAz, *prNode ) =  lcFz; 
                }
                else
                {
			        (*mpModel)( cEx, *prNode ) =  cos_A * lcFx + sin_A * lcFy;
			        (*mpModel)( cEy, *prNode ) = -sin_A * lcFx + cos_A * lcFy;
			        (*mpModel)( cEz, *prNode ) =  lcFz;                
                }
            }
		}
	}

	//***********************************************************************************************************************************************
    //* - Updates discontinuity nodes for the AV formulation.
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_DisctAV()
    {
        if( ( mContactPairs.size() == 0 ) || ( !mPotentials ) || ( mAVContinuity ) )
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        // Update discontinuity 
        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();
            
            if( ( Is_DisctAVNode( NodeId ) ) && ( mContactPairs.find( NodeId ) != mContactPairs.end() ) )
            {
                Node::Pointer pLNode = mpModel->GetNode( mContactPairs[ NodeId ].ctcNode - 1 );

                if( mNormals.find( NodeId ) != mNormals.end() )
                {
                    (*mpModel)( cAx, **nit ) = (*mpModel)( cAx, *pLNode );
			        (*mpModel)( cAz, **nit ) = (*mpModel)( cAz, *pLNode );                   
                }
                else
                {
                    (*mpModel)( cAy, **nit ) = (*mpModel)( cAy, *pLNode );
			        (*mpModel)( cAz, **nit ) = (*mpModel)( cAz, *pLNode );                   
                }
                                
                if( pLNode->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, **nit ) = (*mpModel)( cVs, *pLNode );  
                }                            
            }
        }

        // Rotate nodes
        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();
            
            if( Is_DisctAVNode( NodeId ) ) 
            {
                Matrix< std::complex<double> > RotMatrix;

                Get_Rotation_Matrix_DisctAV( RotMatrix, NodeId );          

                std::complex<double> An = (*mpModel)( cAx, **nit );
			    std::complex<double> At = (*mpModel)( cAy, **nit );
                std::complex<double> Ab = (*mpModel)( cAz, **nit );
                
                // {A_xyz} = [Rt] * {A_ntb}
                (*mpModel)( cAx, **nit ) = RotMatrix[ 0 ][ 0 ]*An + RotMatrix[ 1 ][ 0 ]*At + RotMatrix[ 2 ][ 0 ]*Ab;
                (*mpModel)( cAy, **nit ) = RotMatrix[ 0 ][ 1 ]*An + RotMatrix[ 1 ][ 1 ]*At + RotMatrix[ 2 ][ 1 ]*Ab;
                (*mpModel)( cAz, **nit ) = RotMatrix[ 0 ][ 2 ]*An + RotMatrix[ 1 ][ 2 ]*At + RotMatrix[ 2 ][ 2 ]*Ab;         
            }
        }
    }

	//***********************************************************************************************************************************************
    //* - Updates contact nodes.
    //***********************************************************************************************************************************************
    void Modeler::UpdateCoord_Contact()
    {
        if( mContactPairs.size() == 0 ) 
        {
            return;
        }

		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
			int NodeId = (*nit)->Id();
            
            if( ( !Is_DisctAVNode( NodeId ) ) && ( mContactPairs.find( NodeId ) != mContactPairs.end() ) )
            {
                Node::Pointer pLNode = mpModel->GetNode( mContactPairs[ NodeId ].ctcNode - 1 );

                if( mPotentials )
                {
                    (*mpModel)( cAx, **nit ) = (*mpModel)( cAx, *pLNode );
                    (*mpModel)( cAy, **nit ) = (*mpModel)( cAy, *pLNode );
				    (*mpModel)( cAz, **nit ) = (*mpModel)( cAz, *pLNode );   

                    if( pLNode->pDof( cVs ) != NULL ) 
                    {
                        (*mpModel)( cVs, **nit ) = (*mpModel)( cVs, *pLNode );  
                    }
                }
                else if( mEdgeElementsOn )
                {
                    (*mpModel)( cEx, **nit ) = (*mpModel)( cEx, *pLNode );
                    (*mpModel)( cEy, **nit ) = (*mpModel)( cEy, *pLNode );
				    (*mpModel)( cEz, **nit ) = (*mpModel)( cEz, *pLNode );   
                }
                else
                {
                    Matrix< std::complex<double> > CtCMatrix;
                    
                    if( mColdPlasmaMode )
                    {
                        Get_Contact_Matrix_ColdPlasma( CtCMatrix, NodeId );
                    }
                    else
                    {
                        Get_Contact_Matrix_FullWave( CtCMatrix, NodeId );
                    }
                    
                    std::complex<double> ExL = (*mpModel)( cEx, *pLNode );
                    std::complex<double> EyL = (*mpModel)( cEy, *pLNode );
                    std::complex<double> EzL = (*mpModel)( cEz, *pLNode );
                    
                    (*mpModel)( cEx, **nit ) = CtCMatrix[ 0 ][ 0 ]*ExL + CtCMatrix[ 0 ][ 1 ]*EyL + CtCMatrix[ 0 ][ 2 ]*EzL;
                    (*mpModel)( cEy, **nit ) = CtCMatrix[ 1 ][ 0 ]*ExL + CtCMatrix[ 1 ][ 1 ]*EyL + CtCMatrix[ 1 ][ 2 ]*EzL;
				    (*mpModel)( cEz, **nit ) = CtCMatrix[ 2 ][ 0 ]*ExL + CtCMatrix[ 2 ][ 1 ]*EyL + CtCMatrix[ 2 ][ 2 ]*EzL;                
                }
            }
        }
    }

	//***********************************************************************************************************************************************
    //* - Returns the area weighted normal of a set of normals.
    //* - The module of the normals must be proportional to the element area.
    //***********************************************************************************************************************************************
    Vector<double> Modeler::Calculate_Area_Weighted_Normal( Vector< Vector<double> >& vectorSet )
    {
        // Area weighted normal
        Vector<double> aw_normal( 3, 0.00 );

        // Iterator over the set of normals
        Vector< Vector<double> >::iterator it_vs;

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

    //***********************************************************************************************************************************************
    //* - Returns the geometric averaged normal of a set of normals.
    //***********************************************************************************************************************************************
    Vector<double> Modeler::Calculate_Geom_Averaged_Normal( Vector< Vector<double> >& vectorSet )
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
            v1 *= ( 1.0 / sqrt( std::inner_product( v1.begin(), v1.end(), v1.begin(), 0.0 ) ) );

            // Initiate group of vectors
            Vector<double> v_group(v1);

            // Look for the vectors with a separation angle with v1 lower than 45 degrees
            for ( int jvs = ivs+1; jvs < setSize; jvs++ )
            {
                // If the vector has been assigned -> next vector
                if ( used_vector[jvs] == true ) continue;

                // Take vector from set and normalize
                v2  = vectorSet[jvs];
                v2 *= ( 1.0 / sqrt( std::inner_product( v2.begin(), v2.end(), v2.begin(), 0.0 ) ) );

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
            v_group *= ( 1.0 / sqrt( std::inner_product( v_group.begin(), v_group.end(), v_group.begin(), 0.0 ) ) );

            // Add group resultant to the geometric average normal
            ga_normal += v_group;
        }

        // Normalize geometric average normal
        ga_normal *= ( 1.0 / sqrt( std::inner_product( ga_normal.begin(), ga_normal.end(), ga_normal.begin(), 0.0 ) ) );

        return ga_normal;
    }

	//***********************************************************************************************************************************************
    //* - Sets normals in a PEC node ( nxE = 0 ).
    //***********************************************************************************************************************************************
	void Modeler::Set_PEC_Normals()
    {
		std::map< unsigned int, Vector< Vector<double> > >::iterator it_Set;

        for( it_Set = mSetPECnormals.begin(); it_Set != mSetPECnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            Vector< Vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            if( mNormalsAreaWtd ) 
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
            std::complex<double> cZero( 0.0, 0.0 );

            if( mPotentials )
            {
                (*mpModel)( cAy, *pNode ) = cZero; pNode->pDofcAy()->FixDof();
			    (*mpModel)( cAz, *pNode ) = cZero; pNode->pDofcAz()->FixDof();

                if( pNode->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, *pNode ) = cZero; pNode->pDofcVs()->FixDof(); 
                }
            }
            else
            {
                (*mpModel)( cEy, *pNode ) = cZero; pNode->pDofcEy()->FixDof();
			    (*mpModel)( cEz, *pNode ) = cZero; pNode->pDofcEz()->FixDof();        
            }
            
            mType_Of_BC_Normal[ pNode->Id() ] = 'E';
        }

        // Cleaning mSetPECnormals 
		mSetPECnormals.clear(); std::map<unsigned int, Vector< Vector<double> > >().swap( mSetPECnormals );
	}

	//***********************************************************************************************************************************************
    //* - Sets normals in a PMC node ( n*E = 0 ).
    //***********************************************************************************************************************************************
	void Modeler::Set_PMC_Normals()
    {
		std::map<unsigned int, Vector< Vector<double> > >::iterator it_Set;

		for( it_Set = mSetPMCnormals.begin(); it_Set != mSetPMCnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            Vector< Vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC -> make the normals ortogonal
			if( mNormals.find( it_Set->first ) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[ it_Set->first ] );

				Vector<double> vOrto( 3, 0.0 );

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[ it_Set->first ] = vOrto;
			}
            // if the node is not PEC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[ it_Set->first ] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

                std::complex<double> cZero( 0.0, 0.0 );

                if( mPotentials ) 
                { 
                    (*mpModel)( cAx, *pNode ) = cZero; pNode->pDofcAx()->FixDof(); 
                }
                else                  
                { 
                    (*mpModel)( cEx, *pNode ) = cZero; pNode->pDofcEx()->FixDof(); 
                }

                mType_Of_BC_Normal[ pNode->Id() ] = 'H';
			}
        }

        // Cleaning mSetPMCnormals 
        mSetPMCnormals.clear(); std::map<unsigned int, Vector< Vector<double> > >().swap( mSetPMCnormals );
	}

	//***********************************************************************************************************************************************
    //* - Sets normals in a TEC surfaces ( ortogonal to a PMC surfaces ).
    //***********************************************************************************************************************************************
	void Modeler::Set_TEC_Normals()
    {
		std::map<unsigned int, Vector< Vector<double> > >::iterator it_Set;

		for ( it_Set = mSetTECnormals.begin(); it_Set != mSetTECnormals.end(); ++it_Set )
        {
            // Set of normals per node, being the module of each normal proportional to the area of the element.
            Vector< Vector<double> > vectorSet( it_Set->second );

            // Resultant normal at each node
            Vector<double> normal_atNode;

            // TE_PMC is a symmetry plane -> no need to calculate geometric average
            normal_atNode = Calculate_Area_Weighted_Normal( vectorSet );

            // If the node is a PEC or PMC -> make the normals ortogonal
			if ( mNormals.find( it_Set->first ) != mNormals.end() )
			{
				Vector<double> vToBeOrto( mNormals[ it_Set->first ] );

				Vector<double> vOrto( 3, 0.0 );

				Ortogonalization( vToBeOrto, normal_atNode, vOrto );

				mNormals[ it_Set->first ] = vOrto;
			}
            // If the node is not PEC/PMC -> add normal_atNode to mNodes and set DOFs
			else
			{
				mNormals[ it_Set->first ] = normal_atNode;

				Node::Pointer pNode = mpModel->GetNode( (it_Set->first) - 1 );

                std::complex<double> cZero( 0.0, 0.0 );

                if ( mPotentials ) 
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

        // Cleaning mSetTECnormals 
        mSetTECnormals.clear(); std::map<unsigned int, Vector< Vector<double> > >().swap( mSetTECnormals );
	}

	//***********************************************************************************************************************************************
    //* - Makes the vector vToBeOrto ortogonal to the vector vAxis.
    //***********************************************************************************************************************************************
	void Modeler::Ortogonalization( Vector<double>& vToBeOrto, Vector<double>& vAxis, Vector<double>& vOrto )
    {
		// vOrto = ( vAxis x vToBeOrto ) x vAxis,
		// vOrto is ortogonal to vAxis and pointing in the same direction as vToBeOrto
		Vector<double> nDummy( 3, 0.0 );

		// nDummy = vAxis x vToBeOrto
		nDummy[ 0 ] = vAxis[ 1 ]*vToBeOrto[ 2 ] - vAxis[ 2 ]*vToBeOrto[ 1 ];
		nDummy[ 1 ] = vAxis[ 2 ]*vToBeOrto[ 0 ] - vAxis[ 0 ]*vToBeOrto[ 2 ];
		nDummy[ 2 ] = vAxis[ 0 ]*vToBeOrto[ 1 ] - vAxis[ 1 ]*vToBeOrto[ 0 ];

		// vOrto = nDummy x vAxis
		vOrto[ 0 ] = nDummy[ 1 ]*vAxis[ 2 ] - nDummy[ 2 ]*vAxis[ 1 ];
		vOrto[ 1 ] = nDummy[ 2 ]*vAxis[ 0 ] - nDummy[ 0 ]*vAxis[ 2 ];
		vOrto[ 2 ] = nDummy[ 0 ]*vAxis[ 1 ] - nDummy[ 1 ]*vAxis[ 0 ];

		vOrto *= ( 1.0 / sqrt( std::inner_product( vOrto.begin(), vOrto.end(), vOrto.begin(), 0.0 ) ) );
	}

	//***********************************************************************************************************************************************
    //* - Ortogonalization of the contact normals respect to PEC and PMC normals.
    //***********************************************************************************************************************************************
	void Modeler::OrtogonalizeContactNormals()
    {
		std::map<unsigned int, ContactPairData>::iterator itCtc;

        for( itCtc = mContactPairs.begin(); itCtc != mContactPairs.end(); ++itCtc )
        {
			unsigned int fCtcNode  = itCtc->first;

			if( ( mNormals.find( fCtcNode ) != mNormals.end() ) )
            {
				Vector<double> sNormal  ( mNormals     [ fCtcNode ]           );
				Vector<double> vToBeOrto( mContactPairs[ fCtcNode ].ctcNormal );
				Vector<double> vOrto    ( 3, 0.0 );

				Ortogonalization( vToBeOrto, sNormal, vOrto );

				mContactPairs[ fCtcNode ].ctcNormal = vOrto;
			}
		}
	}

    //***********************************************************************************************************************************************
    //* - Sets normals in PEC and PMC surfaces.
    //***********************************************************************************************************************************************
    void Modeler::Set_Normals()
    {
        // Joining the contacts nodes normals to the Dirichlet set
		JoinContactNormalsInDirichletS();

		// Setting Dirichlet normals
		Set_PEC_Normals();

		// Dirichlet normals are forced in dielectric-Dirichlet intersections.
		ForceDirichletNormalsInContact();

		// Setting simmetry normals
		Set_PMC_Normals();
		Set_TEC_Normals();

		// Contact normals must be ortogonal to simmetry normals
		OrtogonalizeContactNormals();

        // Setting discontinous AV nodes
        Set_DisctAV_Nodes();
    }

    //***********************************************************************************************************************************************
    //* - Sets discontinous AV nodes.
    //***********************************************************************************************************************************************
    void Modeler::Set_DisctAV_Nodes()
    {
        if( mAVContinuity || !mPotentials )
        {
            return;
        }
        
        std::list<int> CtCNodesList;

        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator nit;

        for( nit = nodes.begin(); nit != nodes.end(); nit++ )
        {
            if( mContactPairs.find( (*nit)->Id() ) != mContactPairs.end() )
            {
                CtCNodesList.push_back( (*nit)->Id() );

                CtCNodesList.push_back( mContactPairs[ (*nit)->Id() ].ctcNode );

                mCtCNormals[ mContactPairs[ (*nit)->Id() ].ctcNode ] = mContactPairs[ (*nit)->Id() ].ctcNormal;
            }        
        }

        CtCNodesList.sort();

        CtCNodesList.unique();

        std::list<int>::iterator lit;

        for( lit = CtCNodesList.begin();  lit != CtCNodesList.end(); lit++ ) 
        {
            if( mNormals.find( *lit ) != mNormals.end() )
            {
                if( mType_Of_BC_Normal[ *lit ] == 'E' )
                {
                    continue;
                }

                Node::Pointer pNode = mpModel->GetNode( (*lit)-1 );

                pNode->pDofcAy()->FreeDof();
            }
            else
            {
                Node::Pointer pNode = mpModel->GetNode( (*lit)-1 );

                pNode->pDofcAx()->FreeDof();   
            }
            
            mDisctAVNodes.push_back( *lit );
        }
    }

    //***********************************************************************************************************************************************
    //* - Checks if NodeId is on the discontinous AV nodes list.
    //***********************************************************************************************************************************************
    bool Modeler::Is_DisctAVNode( int NodeId )
    {
        if( mAVContinuity || !mPotentials )
        {
            return false;
        }

        Vector<int>::iterator dit = std::find( mDisctAVNodes.begin(), mDisctAVNodes.end(), NodeId );

        if( dit != mDisctAVNodes.end() )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //***********************************************************************************************************************************************
    //* - Sets nodal elements Degrees Of Freedom (DOFs).
    //***********************************************************************************************************************************************
    void Modeler::Set_Node_Elements_DOFs()
    {
        Set_Contacts();
		Set_Normals ();
		Set_PBC     ();
		
		if     ( mSetAllFzoToZero ) SetAll_cFzo_ToZero();
		else if( mSetAllFxyToZero ) SetAll_cFxy_ToZero();
		else if( mAxisymmetric    ) SetAll_cFyz_ToZero();
    }

    //***********************************************************************************************************************************************
    //* - Sets edge elements Degrees Of Freedom (DOFs).
    //***********************************************************************************************************************************************
    void Modeler::Set_Edge_Elements_DOFs( int* NodesId )
	{
		// Get element nodes pointers
		Vector<Node::Pointer> pNodes;

		Get_Volume_Element_Nodes( NodesId, pNodes );

		// Complex zero
		std::complex<double> cZero( 0.0, 0.0 );		

		// Fx = 0.0 for all the element nodes
		for( int i=0; i<pNodes.size(); i++ )
		{
		    if( mPotentials ) 
			{
				(*mpModel)( cAx, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAx()->FixDof();
			}
            else    
			{
				(*mpModel)( cEx, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEx()->FixDof();
			}
		}

		// Fy = 0.0 for the corner nodes
		for( int i=0; i<4; i++ )
		{
		    if( mPotentials ) 
			{
				(*mpModel)( cAy, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAy()->FixDof();
			}
            else    
			{
				(*mpModel)( cEy, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEy()->FixDof();
			}
		}

		// If stabilization is OFF then Fz = 0.0 for the corner nodes 
		if( mEDGE_stabilize == false )
		{
			for( int i=0; i<4; i++ )
		    {
		        if( mPotentials ) 
		    	{
		    		(*mpModel)( cAz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAz()->FixDof();
		    	}
                else    
		    	{
		    		(*mpModel)( cEz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEz()->FixDof();
		    	}
		    }
		}

		// Fz = 0.0 for all the edge nodes
		int LastEdgeNode = 4;   // no edge nodes to fix for other element orders

		if     ( mElementOrder == 11 ) LastEdgeNode = pNodes.size();
		else if( mElementOrder == 12 ) LastEdgeNode = pNodes.size() - 4;  
		
        for( int i=4; i<LastEdgeNode; i++ )
		{
		    if( mPotentials ) 
			{
				(*mpModel)( cAz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcAz()->FixDof();
			}
            else    
			{
				(*mpModel)( cEz, *pNodes[ i ] ) = cZero; pNodes[ i ]->pDofcEz()->FixDof();
			}
		}
	}

	//***********************************************************************************************************************************************
    //* - Sets all cFz to zero ( 3D-Fxy problems ).
    //***********************************************************************************************************************************************
    void Modeler::SetAll_cFzo_ToZero()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();
		
        Vector<Node::Pointer>::iterator it;
		
        std::complex<double> cZero( 0.0, 0.0 );

		for( it = nodes.begin(); it != nodes.end(); ++it )
		{
			Node::Pointer pNode = (*it);

            if( mPotentials )
            {
                (*mpModel)( cAz, *pNode ) = cZero;

			    pNode->pDofcAz()->FixDof();
            }
            else
            {
                (*mpModel)( cEz, *pNode ) = cZero;

			    pNode->pDofcEz()->FixDof();
            }
		}
	}

	//***********************************************************************************************************************************************
    //* - Sets all cFx cFy to zero ( 3D-Fz problems ).
    //***********************************************************************************************************************************************
    void Modeler::SetAll_cFxy_ToZero()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator it;

		std::complex<double> cZero( 0.00, 0.00 );

		for( it = nodes.begin(); it != nodes.end(); ++it )
		{
			Node::Pointer pNode = (*it);

            if( mPotentials )
            {
			    (*mpModel)( cAx, *pNode ) = cZero;
			    (*mpModel)( cAy, *pNode ) = cZero;

			    pNode->pDofcAx()->FixDof();
			    pNode->pDofcAy()->FixDof();

                if( pNode->pDof( cVs ) != NULL ) 
                {
                    (*mpModel)( cVs, *pNode ) = cZero; 
                    
                    pNode->pDofcVs()->FixDof();
                }
            }
            else
            {
			    (*mpModel)( cEx, *pNode ) = cZero;
			    (*mpModel)( cEy, *pNode ) = cZero;

			    pNode->pDofcEx()->FixDof();
			    pNode->pDofcEy()->FixDof();
            }
		}
	}

	//***********************************************************************************************************************************************
    //* - Sets all cFy cFz to cero ( Axisymmetric problems, actually Fr and Fphs ).
    //***********************************************************************************************************************************************
    void Modeler::SetAll_cFyz_ToZero()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator it;

		std::complex<double> cZero( 0.00, 0.00 );

		for( it = nodes.begin(); it != nodes.end(); ++it )
		{
			Node::Pointer pNode = (*it);
			
            if( mPotentials )
            {
			    (*mpModel)( cAy, *pNode ) = cZero;
			    (*mpModel)( cAz, *pNode ) = cZero;

			    pNode->pDofcAy()->FixDof();
			    pNode->pDofcAz()->FixDof();

                if( ( mNormals.find ( pNode->Id() ) != mNormals.end() ) ||  
                    ( DistanceToAxis( pNode->Id() ) <= mGeoTolerance  )  )
                {
				    (*mpModel)( cAx, *pNode ) = cZero;

				    pNode->pDofcAx()->FixDof();
			    }

                if( pNode->pDof( cVs ) != NULL ) 
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

                if( ( mNormals.find ( pNode->Id() ) != mNormals.end() ) ||  
                    ( DistanceToAxis( pNode->Id() ) <= mGeoTolerance  )  )
                {
				    (*mpModel)( cEx, *pNode ) = cZero;

				    pNode->pDofcEx()->FixDof();
			    }
            }
		}
	}

    //***********************************************************************************************************************************************
    //* - Shows global matrix.
    //***********************************************************************************************************************************************
    void Modeler::Show_Global_Matrix()
    {
        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        for( int i = 0; i < A_matrix.RowsNumber(); i++ )
        {
            data_iterator     = A_matrix[i].begin();
            end_data_iterator = A_matrix[i].end();

            while( data_iterator != end_data_iterator )
            {
                int         i_col = data_iterator->first;
                ComplexType aij   = data_iterator->second;

                std::cout << i << " , " << i_col << " : " << aij << std::endl;

                data_iterator++;
            }
        }

        if( mA_matrix_aux_Required == false ) return;

        for( int i = 0; i < A_matrix_aux.RowsNumber(); i++ )
        {
            data_iterator     = A_matrix_aux[i].begin();
            end_data_iterator = A_matrix_aux[i].end();

            while( data_iterator != end_data_iterator )
            {
                int         i_col = data_iterator->first;
                ComplexType aij   = data_iterator->second;

                std::cout << i << " , " << i_col << " (aux): " << aij << std::endl;

                data_iterator++;
            }
        }
    }

    //***********************************************************************************************************************************************
    //* - Solves linear system.
    //***********************************************************************************************************************************************
    void Modeler::Solve_Linear_System()
    {
        if( mReadSolutionMode )
		{
			if( !ReadSolutionVector( x_vector ) )
			{
				Send_Error_Msg( "Solution file Vector_Xo.bin not found or incomplete", "Please, check that the external solver finished correctly.", 1 );
			}

			A_matrix    .FreeData();
			b_vector    .FreeData();
            A_matrix_aux.FreeData();

            return;
		}

		if( mSolveWithExternal )
		{
			WriteMatrixInFile( A_matrix );
			WriteVectorInFile( b_vector );

			// Remove auxiliary matrix files left by a previous run in a Hermitian storage: external solvers
			// ( e.g. ERMESPETScSolver ) detect the storage format from the presence of these files.
			if( !mA_matrix_aux_Required )
			{
				std::remove( "Matrix_A_aux_int.bin"   );
				std::remove( "Matrix_A_aux_cmplx.bin" );
			}

			A_matrix.FreeData();
			b_vector.FreeData();

            if( mA_matrix_aux_Required )
            {
                WriteAuxMatrixInFile( A_matrix_aux );
            }

            A_matrix_aux.FreeData();

			if( mExternalSolverPath == "" )
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

			// A solution left by a previous run must never be taken as the result of this one
			std::remove( "Vector_Xo.bin" );

            int SolverStatus = system( MSDOSOrder );

			std::cout << std::endl << "External solver has finished ( exit status " << SolverStatus << " )." << std::endl << std::endl;

			if( !ReadSolutionVector( x_vector ) )
			{
				Send_Error_Msg( "External solver did not write a complete solution ( Vector_Xo.bin )", "Please, check the external solver output above and its settings.", 1 );
			}
		}
		else
		{
			std::cout << "Solving linear system..." << std::endl << std::endl;

			if( mReadFileInitialGuess ) 
			{
				ReadInitialGuess( x_vector );
			}

			mpLinearSolver->SetWriteStepResult( mWriteSolutionEveryStep );

            if( !mA_matrix_aux_Required )
            {
                mpLinearSolver->Solve( A_matrix, b_vector, x_vector );
            }
            else if( mpColdPlasma->Is_HermSymm_Matrix() )
            {
                mpLinearSolver->Solve_HermSymm( A_matrix, A_matrix_aux, b_vector, x_vector );
            }
            else
            {
                mpLinearSolver->Solve_HermFull( A_matrix, A_matrix_aux, b_vector, x_vector );
            }

            if( mWriteSolutionFinalStep ) 
			{
				WriteResultsInFile( x_vector );
			}

            std::cout << "Solver finished." << std::endl << std::endl;
		}

        A_matrix    .FreeData();
        b_vector    .FreeData();
        A_matrix_aux.FreeData();
    }

    //***********************************************************************************************************************************************
    //* - Fixes singular nodes and static voltage.
    //***********************************************************************************************************************************************
    void Modeler::FixDof( unsigned int NodeId, const Variable<double>& rVariable, const double& Value )
    {
		int varKey = rVariable.getKey();

		if( varKey == Sg.getKey() )
        {
            mSingular[ NodeId ] = Value;
        }
        else if( varKey == VOLTAGE.getKey() )
        {
            Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

            pNode->pDofV()->FixDof();

            mFix_Static_Voltage[ NodeId ] = Value;
        }
    }

    //***********************************************************************************************************************************************
    //* - Fixes full wave voltage.
    //***********************************************************************************************************************************************
    void Modeler::FixCDof( unsigned int NodeId, const Variable<double>& rVariable, const Vector<double>& vValue )
    {
        if( mPotentials ) 
        {
            Node::Pointer pNode = mpModel->GetNode( NodeId-1 );

            pNode->pDofcVs()->FixDof();

            std::complex<double> mjw( 0.0, -mProblemFrequency );

            std::complex<double> cValue( vValue[ 0 ] * cos( vValue[ 1 ] ), vValue[ 0 ] * sin( vValue[ 1 ] ) );

            mFix_FullWv_Voltage[ NodeId ] = cValue / mjw;
        }
    }

    //***********************************************************************************************************************************************
    //* - True if the internal direct solver is requested: SolverType "Direct_solver" ( from the GiD problemtype ) or the environment variable
    //*   ERMES_SOLVER=Direct, which overrides the solver chosen in GiD without modifying the problemtype.
    //***********************************************************************************************************************************************
    bool Modeler::Use_Direct_Solver( const String& SolverType )
    {
        if( SolverType == "Direct_solver" ) return true;

        const char* pEnv = std::getenv( "ERMES_SOLVER" );

        if( pEnv != NULL )
        {
            std::string Env( pEnv );

            if( ( Env == "Direct" ) || ( Env == "direct" ) || ( Env == "Direct_solver" ) )
            {
                std::cout << "ERMES_SOLVER=" << Env << ": using the internal direct solver." << std::endl << std::endl;

                return true;
            }
        }

        return false;
    }

    //***********************************************************************************************************************************************
    //* - Sets Linear solver in plasma mode.
    //***********************************************************************************************************************************************
    void Modeler::Generate_LinearSolver_ColdPlasma( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
    {
        mSolveWithExternal  = false;

        mExternalSolverPath = "";

        int NumTreads = mItSolverNumThreads;

        if( Use_Direct_Solver( SolverType ) )
        {
            // Only the Full_matrix storage reaches Solve() in plasma mode ( Herm storages use Solve_Herm* )
            mpLinearSolver = ComplexSolver::Pointer( new DirectComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name, true ) );
        }
        else if( SolverType == "Bi_Conjugate_Gradient"  )
        {
            mpLinearSolver = ComplexSolver::Pointer( new BiCG_CP_ComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name ) );
        }
        else if( SolverType == "Quasi_Minimal_Residual" )
        {
            mpLinearSolver = ComplexSolver::Pointer( new QMR_CP_ComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name ) );
        }
        else if( SolverType == "External_solver" )
        {
            mSolveWithExternal = true;

            ReadExternalPathFromFile( mBaseFileName, mExternalSolverPath );
        }
        else
        {
            mpLinearSolver = ComplexSolver::Pointer( new BiCG_CP_ComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name) );
        }
    }

    //***********************************************************************************************************************************************
    //* - Sets linear solver.
    //***********************************************************************************************************************************************
    void Modeler::Generate_LinearSolver_FullWave( String Name, String SolverType, int MaxNumberOfIterations, int StepIterations, double Tolerance )
    {
        mSolveWithExternal  = false;

		mExternalSolverPath = "";

		int NumTreads = mItSolverNumThreads;

		if( Use_Direct_Solver( SolverType ) )
		{
			// Full-wave matrices are complex symmetric ( upper triangle stored )
			mpLinearSolver = ComplexSolver::Pointer( new DirectComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name, false ) );
		}
		else if( SolverType == "Bi_Conjugate_Gradient" )
		{
			mpLinearSolver = ComplexSolver::Pointer( new BiCGComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name ) );
		}
		else if( SolverType == "Quasi_Minimal_Residual" )
		{
			mpLinearSolver = ComplexSolver::Pointer( new QMRComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name ) );
		}
		else if( SolverType == "External_solver" )
		{
			mSolveWithExternal = true;

			ReadExternalPathFromFile( mBaseFileName, mExternalSolverPath );
		}
		else
		{
			mpLinearSolver = ComplexSolver::Pointer( new QMRComplexSolver( NumTreads, MaxNumberOfIterations, StepIterations, Tolerance, Name ) );
		}
    }

    //***********************************************************************************************************************************************
    //* - Sets Degrees Of Freedom (DOF) for the system sparse matrix.
    //***********************************************************************************************************************************************
    void Modeler::Set_System_Matrix_Size()
    {
        int i = 0;

        int j = mDofSet.size();

        std::set<Dof::Pointer,ComparePDof>::iterator it3;

        for( it3 = mDofSet.begin(); it3 != mDofSet.end(); ++it3 )
        {
            if( (*it3)->IsFixed() )
            {
	            (*it3)->SetEquationId( --j );
            }
            else
            {
	            (*it3)->SetEquationId( i++ );
            }
        }

        mSystemSize = j;

        x_vector.resize( mSystemSize, 0.00 );
        b_vector.resize( mSystemSize, 0.00 );
        A_matrix.Resize( mSystemSize       );

        mA_matrix_aux_Required = false;

        if( mColdPlasmaMode )
        {
            if( !mpColdPlasma->Is_Full_Matrix() )
            {
                mA_matrix_aux_Required = true;

                A_matrix_aux.Resize( mSystemSize );
            }
        }
    }

    //***********************************************************************************************************************************************
    //* - Resizes plasma profile vectors if necessary.
    //***********************************************************************************************************************************************
    bool Modeler::Is_Plasma_Profile_Resized()
    {
        if( ( mColdPlasmaMode ) && ( !mPotentials ) && ( mContactPairs.size() > 0 ) && ( mElementOrder == 0 || mElementOrder == 2 ) )
        {    
            mpColdPlasma->Resize_Profile_3D_Vectors( mpModel->GetNodesArraySize() ); 
               
            return true;
        }

        return false;    
    }

    //***********************************************************************************************************************************************
    //* - Interpolates plasma profile on high-order nodes.
    //***********************************************************************************************************************************************
    void Modeler::Plasma_Profile_Interpolate( int* NodesId )
    {
        Vector<int> HONodesId;

        Push_HONodes_OnVolume( NodesId, HONodesId );

        mpColdPlasma->Linear_Profile_Interpolate( HONodesId );
    }

	//***********************************************************************************************************************************************
    //* - Prints projections.
    //***********************************************************************************************************************************************
	void Modeler::Print_Projection()
    {
		if( mProjection.size() > 0 )
		{
		    Create_Directory( "Integrals/" );
            
            Calculate_Sij_Parameters();

		    Print_Sij_OnScreen( mProjectionNormalized, mInputPorts, mProblemFrequency );
            Write_Sij_Headers ( mProjectionNormalized, mInputPorts                    );    
            Write_Sij_OnFiles ( mProjectionNormalized, mInputPorts, mProblemFrequency );    

		    Clear_Projection_Maps();
		}

		if( mVolIntg_E.size() > 0 )
        {
            Create_Directory( "Integrals/" );  
            
            if( mComplxFreqMode )
            {
                Print_VolIntgs_OnScreen( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F, mVolume_m3, mComplexFrequency );  
                Write_VolIntgs_HeadersC( mVolIntg_E );
                Write_VolIntgs_OnFiles ( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F,             mComplexFrequency );                       
            }
            else
            {
                Print_VolIntgs_OnScreen( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F, mVolume_m3, mProblemFrequency );  
                Write_VolIntgs_Headers ( mVolIntg_E );
                Write_VolIntgs_OnFiles ( mVolIntg_E, mVolIntg_H, mVolIntg_B, mVolIntg_J, mVolIntg_F,             mProblemFrequency );            
            }
            
            if( mExportFields )
            {
                Initz_VolIntgs_Fields_OnFiles();
                Write_VolIntgs_Fields_OnFiles();
            }

            Clear_VolumeIntg_Maps();
        }

		if( mSrfIntg_E.size() > 0 )
        {
	        Create_Directory( "Integrals/" );

            if( mComplxFreqMode )
            {
                Print_SrfIntgs_OnScreen( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S, mSurfce_m2, mComplexFrequency );  
                Write_SrfIntgs_HeadersC( mSrfIntg_E );    
                Write_SrfIntgs_OnFiles ( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S,             mComplexFrequency );                      
            }
            else
            {
                Print_SrfIntgs_OnScreen( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S, mSurfce_m2, mProblemFrequency );  
                Write_SrfIntgs_Headers ( mSrfIntg_E );    
                Write_SrfIntgs_OnFiles ( mSrfIntg_E, mSrfIntg_H, mSrfIntg_B, mSrfIntg_J, mSrfIntg_S,             mProblemFrequency );     
            }
           
            if( mExportFields )
            {
                Initz_SrfIntgs_Fields_OnFiles();
                Write_SrfIntgs_Fields_OnFiles();
            }

            Clear_SurfacIntg_Maps();
		}
    }

    //***********************************************************************************************************************************************
    //* - Prints the results of the field integrals in electrostatic mode.
    //***********************************************************************************************************************************************
    void Modeler::Print_Projection_Electrostatic()
    {
        if( mVolume_m3.size() > 0 )
        {
            Print_VolIntg_Electrostatic_OnScreen( mVolIntg_E_Electrostatic, mVolume_m3 );
        }

        if( mSurfce_m2.size() > 0 )
        {
            Print_SrfIntg_Electrostatic_OnScreen( mSrfIntg_E_Electrostatic, mSurfce_m2 );
        }

        Clear_ElecStatic_Maps();
    }

    //***********************************************************************************************************************************************
    //* - Sets problem type.
    //***********************************************************************************************************************************************
    void Modeler::Set_Problem_Type( const String& ProblemType )
    {
		if( ProblemType == "Full_Wave" )
		{
			mProblemType       = E3D  ;
			mSetAllFzoToZero   = false;
			mSetAllFxyToZero   = false;
			mAxisymmetric      = false;
			mColdPlasmaMode    = false;
            mElectrostaticMode = false;
            mFullWaveMode      = true ;
		}
        else if( ProblemType == "Exyz_3D" )
		{
			mSetAllFzoToZero   = false;
			mSetAllFxyToZero   = false;
			mAxisymmetric      = false;
		}
		else if( ProblemType == "Ez_3D" )
		{
			mSetAllFzoToZero   = false;
			mSetAllFxyToZero   = true ;
			mAxisymmetric      = false;
		}
		else if( ProblemType == "Exy_3D" )
		{
			mSetAllFzoToZero   = true ;
			mSetAllFxyToZero   = false;
			mAxisymmetric      = false;
		}
		else if( ProblemType == "Ea_3D" )
		{
			mSetAllFzoToZero   = false;
			mSetAllFxyToZero   = false;
			mAxisymmetric      = true ;
		}
		else if( ProblemType == "Cold_Plasma" )
		{
			mProblemType       = E3D  ;
			mSetAllFzoToZero   = false;
			mSetAllFxyToZero   = false;
			mAxisymmetric      = false;
			mColdPlasmaMode    = true ;
            mElectrostaticMode = false;
            mFullWaveMode      = false;
		}
        else if( ProblemType == "Electrostatic" )
        {
            mProblemType       = E3D  ;
            mSetAllFzoToZero   = false;
            mSetAllFxyToZero   = false;
            mAxisymmetric      = false;
            mColdPlasmaMode    = false;
            mElectrostaticMode = true ;
            mFullWaveMode      = false;
        }
		else if( ProblemType == "Release_Mode" )
		{
			mReleaseSolutionMode = true ;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = false;
		}
		else if( ProblemType == "Debug_Mode" )
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = true ;
			mReadSolutionMode    = false;
		}
		else if( ProblemType == "Read_Mode" )
		{
			mReleaseSolutionMode = false;
			mDebugSolutionMode   = false;
			mReadSolutionMode    = true ;
		}
		else if( ProblemType == "Results_On_Nodes" )
		{
			mResultsOnNodes = true;
            mResultsOnGPs   = 0   ;
		}
		else if( ProblemType == "Results_On_1GP" )
		{
			mResultsOnNodes = false;
            mResultsOnGPs   = 1    ;
		}
        else if( ProblemType == "Results_On_4GP" )
		{
			mResultsOnNodes = false;
            mResultsOnGPs   = 4    ;
		}
        else if( ProblemType == "LL2P_Smoothing_On" )
		{
			mLL2PSmoothing = true;
		}
        else if( ProblemType == "LL2P_Smoothing_Off" )
		{
			mLL2PSmoothing = false;
		}
        else if( ProblemType == "Potentials_On" )
		{
			mPotentials = true;
		}
		else if( ProblemType == "Potentials_Off" )
		{
			mPotentials = false;
		}
        else if( ProblemType == "AV_Continuity_On" )
		{
            mAVContinuity = true;
		}
		else if( ProblemType == "AV_Continuity_Off" )
		{
            mAVContinuity = false;
		}
        else if( ProblemType == "Ini_Surf_Counter" )
		{
			mNonSmooth_Surf_Counter = 0;
		}
	    else if( ProblemType == "Print_Field_Projections" )
		{
			Print_Projection();
		}
        else if( ProblemType == "Print_Static_Projections" )
        {
            Print_Projection_Electrostatic();
        }
		else if( ProblemType == "Geometric_Average_Normals" )
		{
			mNormalsGeomAvg = true ;
			mNormalsAreaWtd = false;
		}
		else if( ProblemType == "Area_Weighted_Normals" )
		{
			mNormalsGeomAvg = false;
			mNormalsAreaWtd = true ;
		}
		else if( ProblemType == "PBC_Cyclic" )
		{
			mIsRightPBCTilted = true;
		}
		else if( ProblemType == "PBC_Periodic" )
		{
            mIsRightPBCTilted = false;
		}
		else if( ProblemType == "Read_Guess_On" )
		{
			mReadFileInitialGuess = true;
		}
		else if( ProblemType == "Read_Guess_Off" )
		{
			mReadFileInitialGuess = false;
		}
		else if( ProblemType == "Write_Solution_Off" )
		{
			mWriteSolutionEveryStep = false;
			mWriteSolutionFinalStep = false;
		}
		else if( ProblemType == "Write_Solution_Steps" )
		{
			mWriteSolutionEveryStep = true;
			mWriteSolutionFinalStep = true;
		}
		else if( ProblemType == "Write_Solution_Final" )
		{
			mWriteSolutionEveryStep = false;
			mWriteSolutionFinalStep = true ;
		}
        else if( ProblemType == "Import_Robin_On" )
		{
			mImportRobinBCs = true;
		}
        else if( ProblemType == "Import_Robin_Off" )
		{
			mImportRobinBCs = false;
		}
        else if( ProblemType == "Import_J_On" )
		{
			mImportCurrents = true;
		}
        else if( ProblemType == "Import_J_Off" )
		{
			mImportCurrents = false;
		}
        else if( ProblemType == "Import_VEM_On" )
		{
			mImportEleMatrx = true;
		}
        else if( ProblemType == "Import_VEM_Off" )
		{
			mImportEleMatrx = false;
		}
        else if( ProblemType == "Export_Fields_On" )
		{
			mExportFields = true;
		}
        else if( ProblemType == "Export_Fields_Off" )
		{
			mExportFields = false;
		}
		else if( ProblemType == "Results_Format_Binary" )
		{
			mOutputFileFormatIsBIN   = true ;
			mOutputFileFormatIsASCII = false;
		}
		else if( ProblemType == "Results_Format_Ascii" )
		{
			mOutputFileFormatIsBIN   = false;
			mOutputFileFormatIsASCII = true ;
		}
        else if( ProblemType == "Check_Consistency" )
        {
            Check_Consistency();
        }
        else if( ProblemType == "Clear_Singularities" )
        {
            Clear_Singularities();
        }

        else if( ProblemType == "LL2P_3sb" ) { mElementOrder =  0; mEdgeElementsOn = false; }
        else if( ProblemType ==  "RME_1st" ) { mElementOrder =  1; mEdgeElementsOn = false; }
        else if( ProblemType ==  "RME_2nd" ) { mElementOrder =  2; mEdgeElementsOn = false; }
        else if( ProblemType ==  "EDG_1st" ) { mElementOrder = 11; mEdgeElementsOn = true ; }
        else if( ProblemType ==  "EDG_2nd" ) { mElementOrder = 12; mEdgeElementsOn = true ; }

        else if( ProblemType == "GeoTol_0"     ) { mGeoTolerance = 00.00; }
        else if( ProblemType == "GeoTol_1"     ) { mGeoTolerance = 1e-01; }
        else if( ProblemType == "GeoTol_1em2"  ) { mGeoTolerance = 1e-02; }
        else if( ProblemType == "GeoTol_1em3"  ) { mGeoTolerance = 1e-03; }
        else if( ProblemType == "GeoTol_1em4"  ) { mGeoTolerance = 1e-04; }
        else if( ProblemType == "GeoTol_1em5"  ) { mGeoTolerance = 1e-05; }
        else if( ProblemType == "GeoTol_1em6"  ) { mGeoTolerance = 1e-06; }
        else if( ProblemType == "GeoTol_1em7"  ) { mGeoTolerance = 1e-07; }
        else if( ProblemType == "GeoTol_1em8"  ) { mGeoTolerance = 1e-08; }
        else if( ProblemType == "GeoTol_1em9"  ) { mGeoTolerance = 1e-09; }
        else if( ProblemType == "GeoTol_1em10" ) { mGeoTolerance = 1e-10; }
        else if( ProblemType == "GeoTol_1em11" ) { mGeoTolerance = 1e-11; }
        else if( ProblemType == "GeoTol_1em12" ) { mGeoTolerance = 1e-12; }

		else if( ProblemType ==  "1pr" ) { mItSolverNumThreads =  1; }
        else if( ProblemType ==  "2pr" ) { mItSolverNumThreads =  2; }
        else if( ProblemType ==  "4pr" ) { mItSolverNumThreads =  4; }
        else if( ProblemType ==  "6pr" ) { mItSolverNumThreads =  6; }
        else if( ProblemType ==  "8pr" ) { mItSolverNumThreads =  8; }
        else if( ProblemType == "10pr" ) { mItSolverNumThreads = 10; }
        else if( ProblemType == "12pr" ) { mItSolverNumThreads = 12; }
        else if( ProblemType == "14pr" ) { mItSolverNumThreads = 14; }
		else if( ProblemType == "16pr" ) { mItSolverNumThreads = 16; }
        else if( ProblemType == "32pr" ) { mItSolverNumThreads = 32; }
        else if( ProblemType == "64pr" ) { mItSolverNumThreads = 64; }

        else if( ProblemType == "RMED_Stab_On"  ) { mRMED_stabilize = true ; }
        else if( ProblemType == "RMED_Stab_Off" ) { mRMED_stabilize = false; }

        else if( ProblemType == "EDGE_Stab_On"  ) { mEDGE_stabilize = true ; }
        else if( ProblemType == "EDGE_Stab_Off" ) { mEDGE_stabilize = false; }

        else if( ProblemType == "LL2P_Stab_On"  ) { mLL2P_hk_factor = 1.0; }
        else if( ProblemType == "LL2P_Stab_Off" ) { mLL2P_hk_factor = 0.0; }
    }

	//***********************************************************************************
	//* - Returns true if the output file format is binary.
	//***********************************************************************************
	bool Modeler::Is_OutputFileFormat_BIN()
	{
		return mOutputFileFormatIsBIN;
	}

	//***********************************************************************************
	//* - Returns true if the output file format is ASCII.
	//***********************************************************************************
	bool Modeler::Is_OutputFileFormat_ASCII()
	{
		return mOutputFileFormatIsASCII;
	}

    //***********************************************************************************
    //* - Returns true if the problem mode is Electrostatic.
    //***********************************************************************************
    bool Modeler::Is_ElectrostaticMode()
    {
        return mElectrostaticMode;
    }

    //***********************************************************************************
    //* - Returns true if the problem mode is Cold_plasma.
    //***********************************************************************************
    bool Modeler::Is_ColdPlasmaMode()
    {
        return mColdPlasmaMode;
    }

    //***********************************************************************************
    //* - Returns true if the problem mode is Full_wave.
    //***********************************************************************************
    bool Modeler::Is_FullWaveMode()
    {
        return mFullWaveMode;
    }

    //***********************************************************************************
    //* - Returns true when smoothing is activated. 
    //***********************************************************************************
    bool Modeler::Is_ResultsOnNodesOn()
    {
        return mResultsOnNodes;
    }

    //***********************************************************************************
    //* - Returns true when LL2P smoothing is activated. 
    //***********************************************************************************
    bool Modeler::Is_LL2PSmoothingOn()
    {
        return mLL2PSmoothing;
    }

    //***********************************************************************************
    //* - Returns true when potentials are in use.
    //***********************************************************************************
    bool Modeler::Is_PotentialsOn()
    {
        return mPotentials;
    }

    //***********************************************************************************
    //* - Returns true when edge elements are in use.
    //***********************************************************************************
    bool Modeler::Is_EdgeElementsOn()
    {
        return mEdgeElementsOn;
    }

    //***********************************************************************************
    //* - Returns element order.
    //***********************************************************************************
    int Modeler::Get_Element_Order()
    {
        return mElementOrder;
    }

    //***********************************************************************************
    //* - Returns number of Gauss points.
    //***********************************************************************************
    int Modeler::Get_Number_Of_GPs()
    {
        return mResultsOnGPs;
    }

    //***********************************************************************************
    //* - Returns problem frequency.
    //***********************************************************************************
    double Modeler::Get_Problem_Frequency()
    {
        return mProblemFrequency;
    }

    //***********************************************************************************
    //* - Sets frequency mode.
    //***********************************************************************************
    void Modeler::Set_Frequency_Mode( bool SweepingFreq )
    {
        mFrequencySweep = SweepingFreq;
        mComplxFreqMode = false;
    }

    //***********************************************************************************
    //* - Sets frequency.
    //***********************************************************************************
    void Modeler::Set_Frequency( double ProblemFrequency )
    {
        mProblemFrequency = ProblemFrequency;
    }

    //**********************************************************************************************************************
    //* - Gives boundary normal at NodeId.
    //**********************************************************************************************************************
    void Modeler::Get_Boundary_Normal( Vector<double>& BoundaryNormal, int NodeId )
    {
        BoundaryNormal.resize( 3 );

        if( mNormals.find( NodeId ) != mNormals.end() )
        {
            BoundaryNormal[ 0 ] = mNormals[ NodeId ][ 0 ];
            BoundaryNormal[ 1 ] = mNormals[ NodeId ][ 1 ];
            BoundaryNormal[ 2 ] = mNormals[ NodeId ][ 2 ];
        }
        else
        {
            BoundaryNormal[ 0 ] = 0.0;
            BoundaryNormal[ 1 ] = 0.0;
            BoundaryNormal[ 2 ] = 0.0;        
        }
    }

    //**********************************************************************************************************************
    //* - Gives contact pair normal at NodeId.
    //**********************************************************************************************************************
    void Modeler::Get_ContactP_Normal( Vector<double>& ContactPNormal, int NodeId )
    {
        ContactPNormal.resize( 3 );

        if( mContactPairs.find( NodeId ) != mContactPairs.end() )
        {
            ContactPNormal[ 0 ] = mContactPairs[ NodeId ].ctcNormal[ 0 ];
            ContactPNormal[ 1 ] = mContactPairs[ NodeId ].ctcNormal[ 1 ];
            ContactPNormal[ 2 ] = mContactPairs[ NodeId ].ctcNormal[ 2 ];
        }
        else
        {
            ContactPNormal[ 0 ] = 0.0;
            ContactPNormal[ 1 ] = 0.0;
            ContactPNormal[ 2 ] = 0.0;        
        }
    }

	//**********************************************************************************************************************
	//* - Reads complex frequency file.
	//**********************************************************************************************************************
	void Modeler::LoadComplexFrequency( String cpxfrq_file_name )
	{
		mComplxFreqMode = true;
        mFrequencySweep = false;
        
        std::fstream CpxFrqFile( cpxfrq_file_name, std::fstream::in );

        Vector< double > CpxFrq( 2 );

		while( true )
		{
			CpxFrqFile >> CpxFrq[ 0 ];
			
            if( CpxFrqFile.eof() ) 
			{
				break;
			}

            CpxFrqFile >> CpxFrq[ 1 ];
		}

        std::complex<double> CpxFreq( CpxFrq[ 0 ], CpxFrq[ 1 ] );

        mComplexFrequency = CpxFreq;
        
		CpxFrqFile.close();
	}

	//**********************************************************************************************************************
	//* - Reads plane waves parameters file.
	//**********************************************************************************************************************
	void Modeler::LoadPWavesParameters( String pwaves_file_name )
	{
        // Open plane waves parameters file
		std::fstream PWavesFile( pwaves_file_name, std::fstream::in );

        // Plane wave parameter data vector
        Vector< double > PWavesData( 8 );

        // Read file until end
		while( true )
		{
            // Mod_F
			PWavesFile >> PWavesData[ 0 ];
			
            if( PWavesFile.eof() ) 
			{
				break;
			}

            // Phase_F
            PWavesFile >> PWavesData[ 1 ];

            // Polarization X, Y, Z
            PWavesFile >> PWavesData[ 2 ];
            PWavesFile >> PWavesData[ 3 ];
            PWavesFile >> PWavesData[ 4 ];

            // Wave vector X, Y, Z
            PWavesFile >> PWavesData[ 5 ];
            PWavesFile >> PWavesData[ 6 ];
            PWavesFile >> PWavesData[ 7 ];

            mPlaneWavesParameters.push_back( PWavesData );
		}

        // Close file
		PWavesFile.close();
	}

	//**********************************************************************************************************************
	//* - Creates plasma object and reads plasma data from files.
	//**********************************************************************************************************************
	void Modeler::LoadPlasmaParameters( String plasma_file_name )
	{
		mpColdPlasma = ColdPlasma::Pointer( new ColdPlasma( mBaseFileName ) );

        mpColdPlasma->Load_Data( plasma_file_name );
	}

	//**********************************************************************************************************************
	//* - Obtains electron density at every node.
	//**********************************************************************************************************************
	void Modeler::Extract_ColdPlasma_ElectronDensity()
	{
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )		
        {
			(*mpModel)( ELECTRON_DENSITY, **it ) = mpColdPlasma->Get_eDensity_OnNode( *it );
		}
	}

	//**********************************************************************************************************************
	//* - Obtains applied external B field at every node.
	//**********************************************************************************************************************
	void Modeler::Extract_ColdPlasma_Bexternal()
	{
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )		
        {
			(*mpModel)( B_EXT, **it ) = mpColdPlasma->Get_Bext_OnNode( *it );
		}
	}

	//**********************************************************************************************************************
	//* - Extracts S,D,P,R,L components of the plasma permittivity tensor at every node.
	//**********************************************************************************************************************
	void Modeler::Extract_ColdPlasma_PermittivityTensor()
	{
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )	
        {
			Vector< std::complex<double> > SDPRL = mpColdPlasma->Get_SDPRL_OnNode( *it, mProblemFrequency );

			(*mpModel)( PLASMA_PERMITTIVITY_S, **it ) = std::real( SDPRL[0] );
			(*mpModel)( PLASMA_PERMITTIVITY_D, **it ) = std::real( SDPRL[1] );
			(*mpModel)( PLASMA_PERMITTIVITY_P, **it ) = std::real( SDPRL[2] );
			(*mpModel)( PLASMA_PERMITTIVITY_R, **it ) = std::real( SDPRL[3] );
			(*mpModel)( PLASMA_PERMITTIVITY_L, **it ) = std::real( SDPRL[4] );
		}
	}

    //**********************************************************************************************************************
    //* - Obtains E field component parallel to B_ext.
    //**********************************************************************************************************************
    void Modeler::Extract_ColdPlasma_Eparallel()
    {
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )
        {
            // External magnetic flux density
            Vector<double> Bext = mpColdPlasma->Get_Bext_OnNode( *it );

            // b = B / |B|
            double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

            Vector<double> b(3);

            b[0] = Bext[0] / Bnorm;
            b[1] = Bext[1] / Bnorm;
            b[2] = Bext[2] / Bnorm;

            // E field
            Vector<double> E_real = (*mpModel)( REAL_E, **it );
            Vector<double> E_imag = (*mpModel)( IMAG_E, **it );

            // E field proyection onto B_ext
            double E_par_real = E_real[0]*b[0] + E_real[1]*b[1] + E_real[2]*b[2];
            double E_par_imag = E_imag[0]*b[0] + E_imag[1]*b[1] + E_imag[2]*b[2];

            double mod_E_par  = std::sqrt( E_par_real*E_par_real + E_par_imag*E_par_imag );

            (*mpModel)( MOD_E_PARALLEL, **it ) = mod_E_par;
        }
    }

    //**********************************************************************************************************************
    //* - Obtains E field component perpendicular to B_ext.
    //**********************************************************************************************************************
    void Modeler::Extract_ColdPlasma_Eperpendicular()
    {
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )
        {
            // External magnetic flux density
            Vector<double> Bext = mpColdPlasma->Get_Bext_OnNode( *it );

            // b = B / |B|
            double Bnorm = std::sqrt( Bext[0]*Bext[0] + Bext[1]*Bext[1] + Bext[2]*Bext[2] );

            Vector<double> b( 3 );

            b[0] = Bext[0] / Bnorm;
            b[1] = Bext[1] / Bnorm;
            b[2] = Bext[2] / Bnorm;

            // E field
            Vector<double> E_real = (*mpModel)( REAL_E, **it );
            Vector<double> E_imag = (*mpModel)( IMAG_E, **it );

            // E field proyection onto B_ext
            double E_par_real = E_real[0]*b[0] + E_real[1]*b[1] + E_real[2]*b[2];
            double E_par_imag = E_imag[0]*b[0] + E_imag[1]*b[1] + E_imag[2]*b[2];

            // Real part of the perpendicular component (Eper = E - Epar)
            Vector<double> E_per_real( 3 );

            E_per_real[0] = E_real[0] - E_par_real*b[0];
            E_per_real[1] = E_real[1] - E_par_real*b[1];
            E_per_real[2] = E_real[2] - E_par_real*b[2];

            // Imaginary part of the perpendicular component (Eper = E - Epar)
            Vector<double> E_per_imag( 3 );

            E_per_imag[0] = E_imag[0] - E_par_imag*b[0];
            E_per_imag[1] = E_imag[1] - E_par_imag*b[1];
            E_per_imag[2] = E_imag[2] - E_par_imag*b[2];

            (*mpModel)( MOD_E_PERPENDICULAR, **it ) = std::sqrt( E_per_real[0]*E_per_real[0] + E_per_imag[0]*E_per_imag[0] +
                                                                 E_per_real[1]*E_per_real[1] + E_per_imag[1]*E_per_imag[1] +
                                                                 E_per_real[2]*E_per_real[2] + E_per_imag[2]*E_per_imag[2] );
        }
    }

	//**********************************************************************************************************************
	//* - Extracts all the plasma parameters at every node.
	//**********************************************************************************************************************
	void Modeler::Extract_ColdPlasma_AllParameters()
	{
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;		
        
        for( it = nodes.begin(); it != nodes.end(); ++it )
		{
			double eDensity;

            Vector<double> Bext;

            Vector< std::complex<double> > SDPRL;

            mpColdPlasma->Get_AllPlasmaParameters_OnNode( eDensity, Bext, SDPRL, *it, mProblemFrequency );

			(*mpModel)( ELECTRON_DENSITY,      **it ) = eDensity;
			(*mpModel)( B_EXT,                 **it ) = Bext    ;
            (*mpModel)( PLASMA_PERMITTIVITY_S, **it ) = std::real( SDPRL[0] );
            (*mpModel)( PLASMA_PERMITTIVITY_D, **it ) = std::real( SDPRL[1] );
            (*mpModel)( PLASMA_PERMITTIVITY_P, **it ) = std::real( SDPRL[2] );
            (*mpModel)( PLASMA_PERMITTIVITY_R, **it ) = std::real( SDPRL[3] );
            (*mpModel)( PLASMA_PERMITTIVITY_L, **it ) = std::real( SDPRL[4] );
		}
	}

    //**********************************************************************************************************************
    //* - Sets to zero all E parallel components when tolerance criteria is met.
    //**********************************************************************************************************************
    void Modeler::SetAll_Epar_ToZero()
    {
        if( mEdgeElementsOn || mpColdPlasma->Is_Epar_Tol_Off() )
        {
            return;
        }

        // E parallel equal to zero when P < 0 and |P| > Tol
        double Tol = mpColdPlasma->Get_Epar_Tol();

        Model::NodesArrayType::iterator it;

        for( it = mpModel->GetNodesArray().begin(); it != mpModel->GetNodesArray().end(); ++it )
        {
            Vector< std::complex<double> > SDPRL = mpColdPlasma->Get_SDPRL_OnNode( *it, mProblemFrequency );

            double P = std::real( SDPRL[2] );

            if( ( P < 0 ) && ( abs( P ) > Tol ) )
            {
                std::complex<double> cZero( 0.00, 0.00 );

                (*mpModel)( cEz, **it ) = cZero;

                (*it)->pDofcEz()->FixDof();
            }
        }
    }

    //**********************************************************************************************************************
    //* - Applies PEC and PMC conditions to H.
    //**********************************************************************************************************************
    void Modeler::Apply_PEC_and_PMC_To_H()
    {
        if ( mNormals.size() == 0 ) 
        {
            return;
        }

        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

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

    //**********************************************************************************************************************
    //* - Applies PEC and PMC conditions to E.
    //**********************************************************************************************************************
    void Modeler::Apply_PEC_and_PMC_To_E()
    {
        if( mNormals.size() == 0 ) 
        {
            return;
        }

        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

        for( it = nodes.begin(); it != nodes.end(); ++it )
        {
            int NodeId = (*it)->Id();

            if( mNormals.find( NodeId ) != mNormals.end() )
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
                if( mType_Of_BC_Normal[ NodeId ] == 'E' )                
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
                else if( mType_Of_BC_Normal[ NodeId ] == 'H' )
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
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }

        Vector<double> vZero( 3, 0.0 );

        Vector<Node::Pointer>::iterator it;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it ) = vZero;
        }
    }

    //**********************************************************************************************************************
    //* - Ends calculation of electrostatic E field.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_E_Derivation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

        Vector<double> Electrostatic_E;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_E = (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it );

            if( NumNeighbours > 0 ) 
            {
                Electrostatic_E /= NumNeighbours;
            }
  
            (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it ) = Electrostatic_E;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to calculate electrostatic J current density.
    //**********************************************************************************************************************
    void Modeler::Ini_Electrostatic_J_Derivation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }

        Vector<double> vZero( 3, 0.0 );

        Vector<Node::Pointer>::iterator it;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it ) = vZero;
        }
    }

    //**********************************************************************************************************************
    //* - Ends calculation of electrostatic J current density.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_J_Derivation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

        Vector<double> Electrostatic_J;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_J = (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it );

            if( NumNeighbours > 0 ) 
            {
                Electrostatic_J /= NumNeighbours;
            }

            (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it ) = Electrostatic_J;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to calculate electrostatic Joule heating.
    //**********************************************************************************************************************
    void Modeler::Ini_Electrostatic_JouleH_Derivation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }

        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }

        Vector<Node::Pointer>::iterator it;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it ) = 0.0;
        }
    }

    //**********************************************************************************************************************
    //* - Ends calculation of electrostatic Joule heating.
    //**********************************************************************************************************************
    void Modeler::End_Electrostatic_JouleH_Derivation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

        double Electrostatic_JouleH;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int NumNeighbours = mNeighbourElements[ (*it)->Id()-1 ];

            Electrostatic_JouleH = (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it );

            if( NumNeighbours > 0 ) 
            {
                Electrostatic_JouleH /= NumNeighbours;
            }

            (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it ) = Electrostatic_JouleH;
        }
    }

    //**********************************************************************************************************************
    //* - Initial set up to derive E field.
    //**********************************************************************************************************************
	void Modeler::Ini_E_Derivation()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;

        if( mPotentials || mLL2PSmoothing || mEdgeElementsOn )
        {
            Vector<double> vZero( 3, 0.0 );

		    for( it=nodes.begin(); it!=nodes.end(); ++it )
            {
		        (*mpModel)( REAL_E, **it ) = vZero;
                (*mpModel)( IMAG_E, **it ) = vZero;
            }                    
            
            mNeighbourElements.resize( nodes.size() ); 
        
            for( int i=0; i<nodes.size(); i++ ) 
            {
                mNeighbourElements[i] = 0;
            }
        }
        else
        {
		    Vector<double> E_real( 3 );
            Vector<double> E_imag( 3 );

		    for( it=nodes.begin(); it!=nodes.end(); ++it )
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
    //* - Ends calculate E field.
    //**********************************************************************************************************************
	void Modeler::End_E_Derivation()
	{
		if( mPotentials || mLL2PSmoothing || mEdgeElementsOn )
        {
            Vector<Node::Pointer> nodes = mpModel->GetNodesArray();
            
            Vector<Node::Pointer>::iterator it;
            
            Vector<double> E_real;
            Vector<double> E_imag;
            
            for( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );
            
                if( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
                {
                    E_real /= mNeighbourElements[ (*it)->Id()-1 ];
                    E_imag /= mNeighbourElements[ (*it)->Id()-1 ];
                }
            
                (*mpModel)( REAL_E, **it ) = E_real;
                (*mpModel)( IMAG_E, **it ) = E_imag;
            }
            
            Apply_PEC_and_PMC_To_E();

            if( mPotentials )
            { 
                std::complex<double> cUnit( 0.0, 1.0 );
                
                for( it=nodes.begin(); it!=nodes.end(); ++it )
                {
                    (*mpModel)( cEx, **it ) = (*mpModel)( REAL_E, **it )[0] + cUnit * (*mpModel)( IMAG_E, **it )[0];
                    (*mpModel)( cEy, **it ) = (*mpModel)( REAL_E, **it )[1] + cUnit * (*mpModel)( IMAG_E, **it )[1];
                    (*mpModel)( cEz, **it ) = (*mpModel)( REAL_E, **it )[2] + cUnit * (*mpModel)( IMAG_E, **it )[2];
                }
            }
        }
	}

    //**********************************************************************************************************************
    //* - Calculates E field on elements.
    //**********************************************************************************************************************
    void Modeler::Calculate_E_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_E_OnGP, ResultsOnGPsType& IMAG_E_OnGP )
    {
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

        Set_Volume_Element_Parameters( NodesId, pElement );

        if( mLL2PSmoothing )
        {
            cVector2Type cE_OnGP;
            
            pElement->Calculate_E_field_OnGaussPoints( cE_OnGP, 1 );

            Vector2Type rE_OnGP, iE_OnGP;
            
            Split_Complex_VectorVector( rE_OnGP, iE_OnGP, cE_OnGP );   

		    for( int n=0; n<pNodes.size(); n++ )
            {
                (*mpModel)( REAL_E, *pNodes[ n ] ) += rE_OnGP[ 0 ];
                (*mpModel)( IMAG_E, *pNodes[ n ] ) += iE_OnGP[ 0 ];
            
                mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
            }    
        }
        else if( mPotentials || mEdgeElementsOn ) 
        {
            cVector2Type cE_OnNodes;
            
            pElement->Calculate_E_field_OnNodes( cE_OnNodes );
            
            Vector2Type rE_OnNodes, iE_OnNodes;
            
            Split_Complex_VectorVector( rE_OnNodes, iE_OnNodes, cE_OnNodes );
            
		    for( int n=0; n<pNodes.size(); n++ )
            {
                (*mpModel)( REAL_E, *pNodes[ n ] ) += rE_OnNodes[ n ];
                (*mpModel)( IMAG_E, *pNodes[ n ] ) += iE_OnNodes[ n ];
            
                mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
            }
        }

        if( mResultsOnNodes == false ) 
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
    //* - Gets the real and imaginary parts of a Vector< Vector<complex> > type.
    //**********************************************************************************************************************
    void Modeler::Split_Complex_VectorVector( Vector2Type& real_VV, Vector2Type& imag_VV, cVector2Type& complex_VV )
    {
        Vector<double> real_V( 3 );
        Vector<double> imag_V( 3 );
 
        cVector2Type::iterator it;

        for( it=complex_VV.begin(); it!=complex_VV.end(); it++ )
        {
            real_V[ 0 ] = std::real( (*it)[ 0 ] );  
            real_V[ 1 ] = std::real( (*it)[ 1 ] );  
            real_V[ 2 ] = std::real( (*it)[ 2 ] ); 

            imag_V[ 0 ] = std::imag( (*it)[ 0 ] ); 
            imag_V[ 1 ] = std::imag( (*it)[ 1 ] ); 
            imag_V[ 2 ] = std::imag( (*it)[ 2 ] );

            real_VV.push_back( real_V );
            imag_VV.push_back( imag_V );
        }
    }

	//**********************************************************************************************************************
    //* - Initial set up to calculate H field.
    //**********************************************************************************************************************
	void Modeler::Ini_H_Derivation()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator it;
        
        Vector<double> vZero( 3, 0.0 );

		for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
		    (*mpModel)( REAL_H, **it ) = vZero;
            (*mpModel)( IMAG_H, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
	}

    //**********************************************************************************************************************
    //* - Ends calculate H field.
    //**********************************************************************************************************************
	void Modeler::End_H_Derivation()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;
        
        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            Vector<double> H_real = (*mpModel)( REAL_H, **it );
            Vector<double> H_imag = (*mpModel)( IMAG_H, **it );

            if( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
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
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator it;

        Vector<double> vZero( 3, 0.0 );

		for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
		    (*mpModel)( REAL_B, **it ) = vZero;
            (*mpModel)( IMAG_B, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
	}

    //**********************************************************************************************************************
    //* - Ends calculate B field.
    //**********************************************************************************************************************
	void Modeler::End_B_Derivation()
	{
		Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        Vector<Node::Pointer>::iterator it;
        
        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            Vector<double> B_real = (*mpModel)( REAL_B, **it );
            Vector<double> B_imag = (*mpModel)( IMAG_B, **it );

            if( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
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
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator it;

        Vector<double> vZero( 3, 0.0 );

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
	        (*mpModel)( REAL_J, **it ) = vZero;
            (*mpModel)( IMAG_J, **it ) = vZero;
        }

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[ i ] = 0;
        }
	}

	//**********************************************************************************************************************
    //* - Calculates H field on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_H_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_H_OnGP, ResultsOnGPsType& IMAG_H_OnGP )
	{
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

        Set_Volume_Element_Parameters( NodesId, pElement );

        // Element material properties
        double mo      = ( 4.0e-7 ) * 3.141592653589793238462643383279; 
        double mu_real = (*pProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		double mu_imag = (*pProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;
        
        std::complex<double> cMu  ( mu_real, mu_imag );
        std::complex<double> cUnit(     0.0,     1.0 );
        std::complex<double> cCteMat;

        if( mPotentials ) 
        {
            cCteMat = 1.0 / ( cMu );
        }
        else
        {
            cCteMat = 1.0 / ( cMu * cUnit * mProblemFrequency );
        }
        
        // Calculate H on nodes
        cVector2Type cH_OnNodes;

        pElement->Calculate_Rotational_OnNodes( cH_OnNodes );

        for( int n=0; n<cH_OnNodes.size(); n++ ) 
        {
            cH_OnNodes[ n ] *= cCteMat;
        }

        Vector2Type rH_OnNodes, iH_OnNodes;

        Split_Complex_VectorVector( rH_OnNodes, iH_OnNodes, cH_OnNodes );

		for( int n=0; n<pNodes.size(); n++ )
        {
            (*mpModel)( REAL_H, *pNodes[ n ] ) += rH_OnNodes[ n ];
            (*mpModel)( IMAG_H, *pNodes[ n ] ) += iH_OnNodes[ n ];

            mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
        }

        // If results on nodes are off then calculate H on GPs 
        if( mResultsOnNodes == false ) 
        {
            cVector2Type cH_OnGPs;
            
            pElement->Calculate_Rotational_OnGaussPoints( cH_OnGPs, mResultsOnGPs );
            
            for( int gp=0; gp<cH_OnGPs.size(); gp++ ) 
            {
                cH_OnGPs[ gp ] *= cCteMat;
            }
            
            Vector2Type rH_OnGPs, iH_OnGPs;
            
            Split_Complex_VectorVector( rH_OnGPs, iH_OnGPs, cH_OnGPs );
            
            REAL_H_OnGP.push_back( rH_OnGPs ); 
            IMAG_H_OnGP.push_back( iH_OnGPs ); 
        }
	}

	//**********************************************************************************************************************
    //* - Calculates B field on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_B_Element( int* NodesId, unsigned int PropertiesId, 
                                       ResultsOnGPsType& REAL_B_OnGP, 
                                       ResultsOnGPsType& IMAG_B_OnGP )
	{
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        double mo      = ( 4.0e-7 ) * 3.141592653589793238462643383279; 
        double mu_real = (*pProperties)( REAL_MAGNETIC_PERMEABILITY ) * mo;
		double mu_imag = (*pProperties)( IMAG_MAGNETIC_PERMEABILITY ) * mo;

        std::complex<double> cMu  ( mu_real, mu_imag );
        std::complex<double> cUnit(     0.0,     1.0 );

        for( int n=0; n<pNodes.size(); n++ )
        {
            Vector<double> rH_OnNodes = (*mpModel)( REAL_H, *pNodes[ n ] );
            Vector<double> iH_OnNodes = (*mpModel)( IMAG_H, *pNodes[ n ] );

            Vector<double> rB_OnNodes( 3 );
            Vector<double> iB_OnNodes( 3 );

            for( int i=0; i<3; i++ )
            {
                std::complex<double> cB_OnNodes = cMu * ( rH_OnNodes[ i ] + cUnit * iH_OnNodes[ i ] );

                rB_OnNodes[ i ] = std::real( cB_OnNodes );
                iB_OnNodes[ i ] = std::imag( cB_OnNodes );
            }

            (*mpModel)( REAL_B, *pNodes[ n ] ) += rB_OnNodes;
            (*mpModel)( IMAG_B, *pNodes[ n ] ) += iB_OnNodes;
            
            mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
        }
       
        if( mResultsOnNodes == false ) 
        {
            Element::Pointer pElement;
            
            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            std::complex<double> cCteMat;
 
            if( mPotentials ) 
            {
                cCteMat = 1.0;
            }
            else  
            {
                cCteMat = 1.0 / ( cUnit * mProblemFrequency );
            }

            cVector2Type cB_OnGPs;
            
            pElement->Calculate_Rotational_OnGaussPoints( cB_OnGPs, mResultsOnGPs );
               
            for( int gp=0; gp<cB_OnGPs.size(); gp++ ) 
            {
                cB_OnGPs[ gp ] *= cCteMat;
            }
            
            Vector2Type rB_OnGPs, iB_OnGPs;
            
            Split_Complex_VectorVector( rB_OnGPs, iB_OnGPs, cB_OnGPs );

            REAL_B_OnGP.push_back( rB_OnGPs ); 
            IMAG_B_OnGP.push_back( iB_OnGPs ); 
        }
	}

    //********************************************************************************************************************************
    //* - Calculates electrostatic E field E = -grad(V).
    //********************************************************************************************************************************
    void Modeler::Calculate_E_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );

        Vector2Type gradV;

        if( mResultsOnNodes )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator it;
            
            int n = 0;
            
            for( it = pNodes.begin(); it != pNodes.end(); ++it )
            {
                (*mpModel)( ELECTROSTATIC_ELECTRIC_FIELD, **it ) -= gradV[n];
            
                mNeighbourElements[ (*it)->Id() - 1 ]++;

                n++;
            }
        }
        else
        {
            pElement->Calculate_Gradient_OnGaussPoints( gradV, mResultsOnGPs );
            
            Vector2Type Electrostatic_E_OnCentralGP = gradV;

            Electrostatic_E_OnCentralGP[0] *= -1.0;
            
            ResultsOnGPs.push_back( Electrostatic_E_OnCentralGP );
        }
    }

    //********************************************************************************************************************************
    //* - Calculates electrostatic current J = -sigma * grad(V).
    //********************************************************************************************************************************
    void Modeler::Calculate_J_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        double sigma = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );

        Vector2Type gradV;
         
        if( mResultsOnNodes )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator it;
            
            int n = 0;
            
            for( it = pNodes.begin(); it != pNodes.end(); ++it )
            {
                Vector<double> temp_J = gradV[n];
            
                temp_J *= sigma;
            
                (*mpModel)( ELECTROSTATIC_CURRENT_DENSITY, **it ) -= temp_J;
            
                if( sigma != 0.0 ) 
                {
                    mNeighbourElements[ (*it)->Id() - 1 ]++;
                }
            
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
    //* - Writes element nodes Ids and the calculated electrostatic current on the export currents file.
    //********************************************************************************************************************************
    void Modeler::Export_Electrostatic_J_Element_GP( int* NodesId, unsigned int PropertiesId )
    {
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        // Checking if the element has an electrical conductivity > 0.0
        double sigma = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        
        if( sigma == 0.0 ) 
        {
            return;
        }
        
        // Continue the computation of J only in elements of conductive materials
        Vector<Node::Pointer> pNodes;

        Vector<int> HONodesId;

        Push_HONodes_OnVolume( NodesId, HONodesId );

        Vector<int>::iterator hoit;

        for( hoit = HONodesId.begin(); hoit != HONodesId.end(); ++hoit )
        {
            pNodes.push_back( mpModel->GetNode( (*hoit)-1 ) );
        }

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );

        Vector2Type gradOnGaussPoints;

        pElement->Calculate_Gradient_OnGaussPoints( gradOnGaussPoints, mResultsOnGPs );

        // J = -sigma*gradV in the central Gauss point
        Vector<double> Electrostatic_J_OnCentralGP = gradOnGaussPoints[0];

        Electrostatic_J_OnCentralGP *= -sigma;

        double modJ2 = Electrostatic_J_OnCentralGP[0] * Electrostatic_J_OnCentralGP[0] + 
                       Electrostatic_J_OnCentralGP[1] * Electrostatic_J_OnCentralGP[1] +
                       Electrostatic_J_OnCentralGP[2] * Electrostatic_J_OnCentralGP[2] ;

        // If |J| == 0.0 do not write the result in the export file
        if( modJ2 <= 0.0 ) 
        {
            return;
        }

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
    //* - Opens the export current file and writes the phase of the current.
    //********************************************************************************************************************************
    void Modeler::Ini_Exporting_Electrostatic_Current()
    {
        String ExpCurrent_FileName; 
        double ExpCurrent_Phase;

        // Reading info from *-9.dat to export currents
        ReadExpCurrentInfoFromFile( mBaseFileName, ExpCurrent_FileName, ExpCurrent_Phase );

        std::cout << "Writing current density on file " << ExpCurrent_FileName << "..." << std::endl;

        // Creating directory to store exported currents
        Create_Directory( "Export_J_Sources/" );
        
        // Open file to write exported currents
        mExportCurrentFile.open( "Export_J_Sources/" + ExpCurrent_FileName, std::fstream::out );

        // Writing phase
        mExportCurrentFile << ExpCurrent_Phase << std::endl;        
    }

    //********************************************************************************************************************************
    //* - Closes the export current file.
    //********************************************************************************************************************************
    void Modeler::End_Exporting_Electrostatic_Current()
    {
        mExportCurrentFile.close();
    }

    //********************************************************************************************************************************
    //* - Adds volumetric element matrices and vectors imported from the binary files "Matrix_K_IVEM.bin" and "Vector_R_IVEM.bin".
    //********************************************************************************************************************************
    void Modeler::Build_Imported_VolElmts()
    {  
        // Check if importing volumetric element matrices mode is activated
        if( !mImportEleMatrx ) return; 

        std::cout << "Importing volume element matrices..." << std::endl;

        // Open volume elements file
        std::fstream VolElements_file( mBaseFileName + "-4.dat", std::fstream::in );  

        // Open R vector binary file
        std::ifstream R_file( "Vector_R_IVEM.bin", std::ios::binary );

        // Open K matrix binary file
        std::ifstream K_file( "Matrix_K_IVEM.bin", std::ios::binary );

        // Line from VolElements_file
        std::string line;

        // Element nodes Ids and material Id from VolElements_file
        int NodesId[ 4 ], MatId;

        // Read files and assemble imported volumetric element matrices and vectors to global system
        while( std::getline( VolElements_file, line ) ) 
        {
            // Ignore comments
            if( line.substr( 0, 2 ) == "//" ) continue;

            // Read formated line
            std::sscanf( line.c_str(), "VE(%d,%d,%d,%d,%d);", &NodesId[0], &NodesId[1], &NodesId[2], &NodesId[3], &MatId );

            // Get high-order volume element nodes
            Vector<int> HONodesId;

		    Vector<Node::Pointer> pNodes;

            Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

            // Volume element definition
            Element::Pointer pElement;

		    Properties::Pointer pProperties = mpModel->GetProperties( MatId ); 

            this->SetProperties( MatId, FREQUENCY, mProblemFrequency );

            if( mColdPlasmaMode ) 
            {
                if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_ColdPlasma( pNodes, pProperties ) );
                else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_ColdPlasma( pNodes, pProperties ) );
		        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_ColdPlasma( pNodes, pProperties ) );
                else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_ColdPlasma( pNodes, pProperties ) );
                else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_ColdPlasma( pNodes, pProperties ) );
            }
            else
            {
                if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
		        else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
                else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );            
            }

            pElement->SetPotentials( mPotentials );

            // Get DOFs Ids vector
		    Vector<int> EleIdsVector;

            pElement->GetEquationIdVector( EleIdsVector );

            // Calculate vector and matrix sizes
            int VectorSize = EleIdsVector.size(); 
            int MatrixSize = VectorSize * VectorSize; 

            int Chunck_Size_Vector = VectorSize * sizeof( std::complex<double> );
            int Chunck_Size_Matrix = MatrixSize * sizeof( std::complex<double> );

            // Read element stiffness matrix from file
            Matrix< std::complex<double> > EleStiffMatrix;

            EleStiffMatrix.Resize( VectorSize, VectorSize, std::complex<double>(0.0, 0.0) );

            Vector< std::complex<double> > K( MatrixSize, std::complex<double>(0.0, 0.0) );

            K_file.read( reinterpret_cast<char*>( K.data() ), Chunck_Size_Matrix );

            if( K_file.gcount() == Chunck_Size_Matrix )
            {
                for( int i=0; i<VectorSize; i++ )
                {
                    for( int j=0; j<VectorSize; j++ )
                    {
                        EleStiffMatrix[ i ][ j ] = K[ i * VectorSize + j ];
                    }
                }
            }

            // Read element residual vector from file
            Vector< std::complex<double> > EleResVector;

            pElement->GetResidualVector_Dirichlet( mFix_FullWv_Voltage, EleStiffMatrix, EleResVector );

            Vector< std::complex<double> > R( VectorSize, std::complex<double>(0.0, 0.0) );

            R_file.read( reinterpret_cast<char*>( R.data() ), Chunck_Size_Vector );

            if( R_file.gcount() == Chunck_Size_Vector ) 
            {
                EleResVector += R;            
            }

            // Apply boundary conditions to stiffness matrix
            Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

            // Apply boundary conditions to residual vector
            Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

            // Assembling to global matrix 
            if( mColdPlasmaMode ) 
            {
                if     ( mpColdPlasma->Is_Full_Matrix    () ) Assemble_ElementMatrix_to_Global_NonSymmetric    ( EleIdsVector, EleStiffMatrix );
                else if(    !pElement->GetIsIHL          () ) Assemble_ElementMatrix_to_Global_Symmetric       ( EleIdsVector, EleStiffMatrix );
                else if( mpColdPlasma->Is_HermSymm_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_Symmetric   ( EleIdsVector, EleStiffMatrix );           
                else if( mpColdPlasma->Is_HermFull_Matrix() ) Assemble_ElementMatrix_to_Global_Aux_NonSymmetric( EleIdsVector, EleStiffMatrix ); 
            }
            else    
            {
                Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
            }

            // Assembling elemental residual vector to global residual vector
            Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

            // Clean objects
            K             .FreeData();
            R             .FreeData();
            EleIdsVector  .FreeData();
            EleResVector  .FreeData();
            EleStiffMatrix.FreeData();
        }

        K_file          .close();
        R_file          .close();
        VolElements_file.close();

        std::cout << "Importing volume element matrices finished." << std::endl;
    }

    //********************************************************************************************************************************
    //* - Adds Robin boundary conditions imported from the binary files "Matrix_P_IRBC.bin" and "Vector_U_IRBC.bin".
    //********************************************************************************************************************************
    void Modeler::Build_Imported_RobinBCs()
    {  
        // Check if importing Robin boundary condition mode is activated
        if( !mImportRobinBCs ) return; 

        std::cout << "Importing Robin boundary conditions..." << std::endl;

        // Open Robin boundary conditions elements file
        std::fstream RBCElements_file( mBaseFileName + "-19.dat", std::fstream::in );  

        // Open U vector binary file
        std::ifstream U_file( "Vector_U_IRBC.bin", std::ios::binary );

        // Open P matrix binary file
        std::ifstream P_file( "Matrix_P_IRBC.bin", std::ios::binary );

        // Line from RBCElements_file
        std::string line;

        // Surface element nodes Ids and Surface Id from RBCElements_file
        int NodesId[ 3 ], SurfId;

        // Read files and assemble Robin condition to global system
        while( std::getline( RBCElements_file, line ) ) 
        {
            std::stringstream LineStrm( line );

            LineStrm >> NodesId[ 0 ]; 
            LineStrm >> NodesId[ 1 ];
            LineStrm >> NodesId[ 2 ];
            LineStrm >> SurfId;

            Vector<int> HONodesId;

            Vector<Node::Pointer> pNodes;

		    Get_Surface_Element_Nodes( NodesId, HONodesId, pNodes );

            Element::Pointer pElement;

		    Properties::Pointer pProperties; 

            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new GenericRobin_3sb_FullWave( pNodes, pProperties ) );
            else if( mElementOrder ==  1 ) pElement = Element::Pointer( new GenericRobin_1st_FullWave( pNodes, pProperties ) );
            else if( mElementOrder ==  2 ) pElement = Element::Pointer( new GenericRobin_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new GenericRobin_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new GenericRobin_2eg_FullWave( pNodes, pProperties ) );

            pElement->SetPotentials( mPotentials ); 
            
            // Get DOFs Ids vector
		    Vector<int> EleIdsVector;

            pElement->GetEquationIdVector( EleIdsVector );

            // Calculate vector and matrix sizes
            int VectorSize = EleIdsVector.size(); 
            int MatrixSize = VectorSize * VectorSize; 

            int Chunck_Size_Vector = VectorSize * sizeof( std::complex<double> );
            int Chunck_Size_Matrix = MatrixSize * sizeof( std::complex<double> );

            // Read element residual vector from file
            Vector< std::complex<double> > EleResVector( VectorSize, std::complex<double>(0.0, 0.0) );
            Vector< std::complex<double> > U           ( VectorSize, std::complex<double>(0.0, 0.0) );

            U_file.read( reinterpret_cast<char*>( U.data() ), Chunck_Size_Vector );

            if( U_file.gcount() == Chunck_Size_Vector ) 
            {
                EleResVector = U;            
            }

            // Read element stiffness matrix from file
            Matrix< std::complex<double> > EleStiffMatrix;

            EleStiffMatrix.Resize( VectorSize, VectorSize, std::complex<double>(0.0, 0.0) );

            Vector< std::complex<double> > P( MatrixSize, std::complex<double>(0.0, 0.0) );

            P_file.read( reinterpret_cast<char*>( P.data() ), Chunck_Size_Matrix );

            if( P_file.gcount() == Chunck_Size_Matrix )
            {
                for( int i=0; i<VectorSize; i++ )
                {
                    for( int j=0; j<VectorSize; j++ )
                    {
                        EleStiffMatrix[ i ][ j ] = P[ i * VectorSize + j ];
                    }
                }
            }

            // Apply boundary conditions to stiffness matrix
            Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleStiffMatrix );

            // Apply boundary conditions to residual vector
            Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );

            // Assembling to global matrix 
            if( mColdPlasmaMode ) 
            {
                Assemble_BC_Plasma_Element_On_AuxMatrix( EleIdsVector, EleStiffMatrix );    
            }
            else    
            {
                Assemble_ElementMatrix_to_Global_Symmetric( EleIdsVector, EleStiffMatrix );  
            }

            // Assembling elemental residual vector to global residual vector
            Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );

            // Clean objects
            P             .FreeData();
            U             .FreeData();
            EleIdsVector  .FreeData();
            EleResVector  .FreeData();
            EleStiffMatrix.FreeData();
        }

        P_file          .close();
        U_file          .close();
        RBCElements_file.close();

        std::cout << "Importing Robin boundary conditions finished." << std::endl;
    }

    //********************************************************************************************************************************
    //* - Adds to the b_vector the J currents from the folder "Export_J_Sources/".
    //********************************************************************************************************************************
    void Modeler::Build_Imported_Currents()
    {        
        // Check if import currents mode is activated
        if( !mImportCurrents ) return; 

        // Objectes from "dirent.h" 
        DIR           *ExpJs_Dir;
        struct dirent *Dir_Ent  ;

        // Name of the folder where the exported J currents are stored
        String ExpJs_Dir_Name = "Export_J_Sources/";

        // If the directory "Export_J_Sources/" does not exist then exit function
        if( ( ExpJs_Dir = opendir( ExpJs_Dir_Name ) ) == NULL ) return;

        std::cout << "Importing J sources..." << std::endl;

        // Vector with cJ on nodes
        Vector< Vector< std::complex<double> > > Import_cJ_Nodal; 
       
        // Read all the files inside the directoy "Export_J_Sources/"
        while( ( Dir_Ent = readdir( ExpJs_Dir ) ) != NULL ) 
        {
            // If the entity is not a file go to the next entity
            if( Dir_Ent->d_type != DT_REG ) continue; 
                      
            // Saving file name as a string
            std::string File_Name = std::string( ExpJs_Dir_Name ) + std::string( Dir_Ent->d_name );

            // Open J source file
            std::fstream ImpCurrentFile( File_Name, std::fstream::in );  

            // Read first line
            std::string first_line; 
            
            ImpCurrentFile >> first_line;

            // If first character == 'N': file is read in nodal mode.
            // If first character == 'N': read J nodal file, store it in
            // mImported_cJ_Nodal, close file, and continue to next file.
            if( first_line[0] == 'N' )
            { 
                Read_Imported_Nodal_Current( ImpCurrentFile, Import_cJ_Nodal );
                continue;
            }

            // If first character != 'N': file is read in volumetric element mode.
            // In volumetric element mode, two options are available:
            // - is_a_global_phase_file == true : all the Jxyz in the file have the same phase.
            // - is_a_global_phase_file == false: each Jxyz in each element has its own phase.
            bool is_a_global_phase_file;

            // Global phase value
            double global_phase_value;
            
            // If first character == 'V', the file has a local phase for each element and J component.
            // If first character != 'V', it will by a number, which is the global phase of the file.
            if( first_line[0] == 'V' )
            { 
                is_a_global_phase_file = false; 
            }
            else                      
            { 
                is_a_global_phase_file = true; 
                global_phase_value     = std::stod( first_line );
            }

            // Element nodes Ids
            int NodesId[ 4 ];

            // Element nodes Jxyz components and Jxyz phases
            Vector<double> Element_J( 6 ); 

            // Reading the rest of the file
            while( true ) 
            {
                // Reading Node Ids from file
                ImpCurrentFile >> NodesId[0]; 
                ImpCurrentFile >> NodesId[1];
                ImpCurrentFile >> NodesId[2];
                ImpCurrentFile >> NodesId[3];
                
                // Reading J source from file
                if( is_a_global_phase_file )
                {
                    ImpCurrentFile >> Element_J[0]; Element_J[1] = global_phase_value;
                    ImpCurrentFile >> Element_J[2]; Element_J[3] = global_phase_value;
                    ImpCurrentFile >> Element_J[4]; Element_J[5] = global_phase_value;
                }
                else
                {
                    ImpCurrentFile >> Element_J[0]; ImpCurrentFile >> Element_J[1];
                    ImpCurrentFile >> Element_J[2]; ImpCurrentFile >> Element_J[3];
                    ImpCurrentFile >> Element_J[4]; ImpCurrentFile >> Element_J[5];                
                }

                // Check if End Of File
                if( ImpCurrentFile.eof() ) break;

                // Assembling imported source element to b_vector
                Assemble_Imported_Current_Element( NodesId, Element_J );
            }

            // Closing J source file
            ImpCurrentFile.close(); 
        }

        // Close directory Export_J_Sources/
        closedir( ExpJs_Dir );

        // Assemble all nodal J files
        Assemble_Imported_Current_OnNodes( Import_cJ_Nodal );
       
        std::cout << "Importing J sources finished." << std::endl;
    }

	//********************************************************************************************************************************
    //* - It reads J nodal file, store it in mImported_cJ_Nodal vector, and close file.
    //********************************************************************************************************************************
    void Modeler::Read_Imported_Nodal_Current( std::fstream& Import_J_File, Vector< Vector< std::complex<double> > >& Import_cJ_Nodal )
    {
        // Initiate Import_cJ_Nodal vector
        if( Import_cJ_Nodal.size() == 0 )
        {
            std::fstream Nodes_File( mBaseFileName + "-1.dat", std::fstream::in );  

            std::string line;

            int NodeId; double X, Y, Z;

            int MaxNodeId = 0;
                
            while( std::getline( Nodes_File, line ) )
            {
                if( line.substr( 0, 2 ) == "//" ) continue;

                std::sscanf( line.c_str(), "No[%d] = p(%lf,%lf,%lf);", &NodeId, &X, &Y, &Z );

                if( NodeId >= MaxNodeId ) MaxNodeId = NodeId;
            }

            Nodes_File.close();
            
            Import_cJ_Nodal.resize( MaxNodeId );

            for( int i=0; i<MaxNodeId; i++)
            { 
                Import_cJ_Nodal[ i ].resize( 3, std::complex<double>( 0.0, 0.0 ) );
            }
        }

        // Read imported nodal current file
        int    nId;
        double mod_Jx, pha_Jx;
        double mod_Jy, pha_Jy;
        double mod_Jz, pha_Jz;

        while( true ) 
        {
            Import_J_File >> nId; 
            Import_J_File >> mod_Jx; Import_J_File >> pha_Jx;
            Import_J_File >> mod_Jy; Import_J_File >> pha_Jy;
            Import_J_File >> mod_Jz; Import_J_File >> pha_Jz;

            if( Import_J_File.eof() ) break;

            if( nId > Import_cJ_Nodal.size() )
            {
                Send_Error_Msg( "Index out of range in imported J file", "Please, check imported J files.", 1 );
            }
            
            std::complex<double> cJx( mod_Jx * cos( pha_Jx ),  mod_Jx * sin( pha_Jx ) );
            std::complex<double> cJy( mod_Jy * cos( pha_Jy ),  mod_Jy * sin( pha_Jy ) );
            std::complex<double> cJz( mod_Jz * cos( pha_Jz ),  mod_Jz * sin( pha_Jz ) );

            Import_cJ_Nodal[ nId - 1 ][ 0 ] += cJx;
            Import_cJ_Nodal[ nId - 1 ][ 1 ] += cJy;
            Import_cJ_Nodal[ nId - 1 ][ 2 ] += cJz;
        }

        Import_J_File.close();
    }

	//********************************************************************************************************************************
    //* - Assembles imported current in nodal format.
    //********************************************************************************************************************************
    void Modeler::Assemble_Imported_Current_OnNodes( Vector< Vector< std::complex<double> > >& Import_cJ_Nodal )
    {
        // If Import_cJ_Nodal empty, return
        if( Import_cJ_Nodal.size() == 0 ) return;

        // Open volume elements file
        std::fstream VolElements_file( mBaseFileName + "-4.dat", std::fstream::in );  

        std::string line;

        int NodesId[ 4 ], MatId;

        // Read volume elements files and assemble imported nodal J to global system
        while( std::getline( VolElements_file, line ) ) 
        {
            if( line.substr( 0, 2 ) == "//" ) continue;

            std::sscanf( line.c_str(), "VE(%d,%d,%d,%d,%d);", &NodesId[0], &NodesId[1], &NodesId[2], &NodesId[3], &MatId );

		    Vector<Node::Pointer> pNodes; Vector<int> HONodesId;

            Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

            Properties::Pointer pProperties( new Kratos::Properties() );

 	        PropertyFunction< Vector<double> >::Pointer Cartesian_J( new ConstantProperty< Vector<double> >( Vector<double>( 6, 0.0 ) ) );
            PropertyFunction< Vector<double> >::Pointer Null_Prop_2( new ConstantProperty< Vector<double> >( Vector<double>( 2, 0.0 ) ) );
            PropertyFunction< Vector<double> >::Pointer Null_Prop_6( new ConstantProperty< Vector<double> >( Vector<double>( 6, 0.0 ) ) );

            // Properties definition as in 11-Materials_Properties.bas input file
            pProperties->SetProperty( COMPLEX_IBC               , Null_Prop_2 );
            pProperties->SetProperty( COMPLEX_IBC_2o            , Null_Prop_2 );
            pProperties->SetProperty( COMPLEX_NEUMANN_FLOW      , Null_Prop_6 );
            pProperties->SetProperty( SINUSOIDAL_SURFACE_CURRENT, Cartesian_J ); 

            // Set frequency in element
            PropertyFunction< double >::Pointer ProblemFrequency( new ConstantProperty< double >( mProblemFrequency ) );

            pProperties->SetProperty( FREQUENCY, ProblemFrequency );   

            // Define volumetric source element
            Element::Pointer pElement;

            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, pProperties ) );
	        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, pProperties ) );
	        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new JSource_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new JSource_2eg_FullWave( pNodes, pProperties ) );

            // Activates potentials on element
            pElement->SetPotentials( mPotentials );

            // Setting nodal J for this element
            Vector< Vector< std::complex<double> > > Nodal_J( 4 );

            for( int i=0; i<Nodal_J.size(); i++ )
            {
                Nodal_J[ i ] = Import_cJ_Nodal[ NodesId[ i ] - 1 ];
            }

            pElement->Set_Cartesian_J_Nodal( Nodal_J );

            // Get DOFs Id vector
            Vector<int> EleIdsVector;

            pElement->GetEquationIdVector( EleIdsVector );

            // Get elemental residual vector
            Vector< std::complex<double> > EleResVector; 

            pElement->GetResidualVector( EleResVector );
  
            // Apply boundary conditions to force vector
            Apply_Element_Boundary_Conditions( HONodesId, EleIdsVector, EleResVector );
	  
            // Assembling elemental residual vector to global residual vector
            Assemble_ResidualVector_to_Global( EleIdsVector, EleResVector );
    
            // Cleaning vectors 
            EleIdsVector.FreeData();
            EleResVector.FreeData();
        }

        // Clear Import_cJ_Nodal
        Import_cJ_Nodal.FreeData();
    }

	//********************************************************************************************************************************
    //* - Assembles imported current element.
    //********************************************************************************************************************************
    void Modeler::Assemble_Imported_Current_Element( int* NodesId, Vector<double>& Element_J )
    {
        // Get volume element nodes
		Vector<Node::Pointer> pNodes;

        Vector<int> HONodesId;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Element properties 
        Properties::Pointer pProperties( new Kratos::Properties() );

        // Source properties
 	    PropertyFunction< Vector<double> >::Pointer Cartesian_J( new ConstantProperty< Vector<double> >( Element_J                ) );
        PropertyFunction< Vector<double> >::Pointer Null_Prop_2( new ConstantProperty< Vector<double> >( Vector<double>( 2, 0.0 ) ) );
        PropertyFunction< Vector<double> >::Pointer Null_Prop_6( new ConstantProperty< Vector<double> >( Vector<double>( 6, 0.0 ) ) );

        // Properties definition as in 11-Materials_Properties.bas input file
        pProperties->SetProperty( COMPLEX_IBC               , Null_Prop_2 );
        pProperties->SetProperty( COMPLEX_IBC_2o            , Null_Prop_2 );
        pProperties->SetProperty( COMPLEX_NEUMANN_FLOW      , Null_Prop_6 );
        pProperties->SetProperty( SINUSOIDAL_SURFACE_CURRENT, Cartesian_J ); 

        // Problem frequency
        PropertyFunction< double >::Pointer ProblemFrequency( new ConstantProperty< double >( mProblemFrequency ) );

        pProperties->SetProperty( FREQUENCY, ProblemFrequency );   

        // Volumetric source element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new JSource_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new JSource_2eg_FullWave( pNodes, pProperties ) );

        // Activates potentials on element
        pElement->SetPotentials( mPotentials );

        // Get DOFs Id vector
        Vector<int> EleIdsVector;

        pElement->GetEquationIdVector( EleIdsVector );

        // Get elemental residual vector
        Vector< std::complex<double> > EleResVector; 

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
    //* - Calculates Joule heating Q = sigma * ||grad(V)||^2. 
    //********************************************************************************************************************************
    void Modeler::Calculate_JouleH_Element_Electrostatic( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& ResultsOnGPs )
    {
        Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        double sigma = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );

        Element::Pointer pElement = Element::Pointer( new VolumeElement_1st_Electrostatic( pNodes, pProperties ) );

        Vector2Type gradV;

        double gradV2, JouleHeating;

        if( mResultsOnNodes )
        {
            pElement->Calculate_Gradient_OnNodes( gradV );
            
            Vector<Node::Pointer>::iterator it;
            
            int n = 0;
            
            for( it = pNodes.begin(); it != pNodes.end(); ++it )
            {
                gradV2 = gradV[n][0]*gradV[n][0] + gradV[n][1]*gradV[n][1] + gradV[n][2]*gradV[n][2];

                JouleHeating = sigma * gradV2;

                (*mpModel)( ELECTROSTATIC_JOULE_HEATING, **it ) += JouleHeating;

                if( sigma != 0.0 ) 
                {
                    mNeighbourElements[ (*it)->Id() - 1 ]++;
                }

                n++;
            }        
        }
        else
        {
            pElement->Calculate_Gradient_OnGaussPoints( gradV, mResultsOnGPs );
            
            gradV2 = gradV[0][0]*gradV[0][0] + gradV[0][1]*gradV[0][1] + gradV[0][2]*gradV[0][2];

            JouleHeating = sigma * gradV2;

            Vector<double> JouleHeating_OnCentralGP; 

            JouleHeating_OnCentralGP.push_back( JouleHeating );

            Vector2Type JouleHeating_OnGPs; 

            JouleHeating_OnGPs.push_back( JouleHeating_OnCentralGP );
            
            ResultsOnGPs.push_back( JouleHeating_OnGPs );
        }
    }

    //********************************************************************************************************************************
    //* - Calculates Joule heating.
    //********************************************************************************************************************************
    void Modeler::Calculate_JouleH_Element_FullWave( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP )
    {
        if( mResultsOnNodes )
        {
            Calculate_Joule_Heating_OnNodes( NodesId, PropertiesId );
        }
        else
        {
            Calculate_Joule_Heating_OnGaussPoints( NodesId, PropertiesId, JOULE_HEATING_OnGP );
        }
    }

    //********************************************************************************************************************************
    //* - Initiates Joule heating calculation.
    //********************************************************************************************************************************
    void Modeler::Ini_Joule_Heating_Calculation()
    {
        if( mResultsOnNodes == false ) 
        {
            return;
        }
        
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

	    Vector<Node::Pointer>::iterator it;

        for( it=nodes.begin(); it!=nodes.end(); ++it )
        {
	        (*mpModel)( JOULE_HEATING, **it ) = 0.0;
        }    

        mNeighbourElements.resize( nodes.size() ); 
        
        for( int i=0; i<nodes.size(); i++ ) 
        {
            mNeighbourElements[i] = 0;
        }
    }

    //********************************************************************************************************************************
    //* - Ends Joule heating calculation.
    //********************************************************************************************************************************
    void Modeler::End_Joule_Heating_Calculation()
    {
        if( mResultsOnNodes ) 
        {
            Vector<Node::Pointer> nodes = mpModel->GetNodesArray();
            
            Vector<Node::Pointer>::iterator it;
            
            double Joule_Heating;
            
            for( it=nodes.begin(); it!=nodes.end(); ++it )
            {
                Joule_Heating = (*mpModel)( JOULE_HEATING, **it );
               
                if( mNeighbourElements[ (*it)->Id()-1 ] > 0 ) 
                {
                    Joule_Heating /= mNeighbourElements[ (*it)->Id()-1 ];
                }
            
                (*mpModel)( JOULE_HEATING, **it ) = Joule_Heating;           
            }   
        }
    }

	//********************************************************************************************************************************
    //* - Calculates Joule heating on nodes.
    //********************************************************************************************************************************
	void Modeler::Calculate_Joule_Heating_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
		// Joule heating in IHL material with losses
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		double eo   = 8.8541878176e-12;
		double freq = mProblemFrequency;

		double sgm_real = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double eps_imag = (*pProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		double eff_cond = freq * eps_imag + sgm_real;

        if( eff_cond > 0 )
        {
		    Vector<Node::Pointer> pNodes;
            
            Get_Volume_Element_Nodes( NodesId, pNodes );

            this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

            Element::Pointer pElement;
            
            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );
            
            Set_Volume_Element_Parameters( NodesId, pElement );

            if( mLL2PSmoothing )
            {
                cVector2Type cE_OnGP;
            
                pElement->Calculate_E_field_OnGaussPoints( cE_OnGP, 1 );

                Vector2Type rE_OnGP, iE_OnGP;
            
                Split_Complex_VectorVector( rE_OnGP, iE_OnGP, cE_OnGP );   

                double modE2 = rE_OnGP[ 0 ][ 0 ]*rE_OnGP[ 0 ][ 0 ] + iE_OnGP[ 0 ][ 0 ]*iE_OnGP[ 0 ][ 0 ] +
		        			   rE_OnGP[ 0 ][ 1 ]*rE_OnGP[ 0 ][ 1 ] + iE_OnGP[ 0 ][ 1 ]*iE_OnGP[ 0 ][ 1 ] +
		        			   rE_OnGP[ 0 ][ 2 ]*rE_OnGP[ 0 ][ 2 ] + iE_OnGP[ 0 ][ 2 ]*iE_OnGP[ 0 ][ 2 ] ;
                
		        for( int n=0; n<pNodes.size(); n++ )
                {
                    (*mpModel)( JOULE_HEATING, *pNodes[ n ] ) += 0.5 * eff_cond * modE2;
                
                    mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
                }    
            }
            else
            {
                cVector2Type cE;
                
                pElement->Calculate_E_field_OnNodes( cE );
                
                Vector2Type rE, iE;
                
                Split_Complex_VectorVector( rE, iE, cE );
                
                double modE2;
                
                for( int n=0; n<pNodes.size(); n++ )
                {
                    modE2 = rE[ n ][ 0 ]*rE[ n ][ 0 ] + iE[ n ][ 0 ]*iE[ n ][ 0 ] +
		        			rE[ n ][ 1 ]*rE[ n ][ 1 ] + iE[ n ][ 1 ]*iE[ n ][ 1 ] +
		        			rE[ n ][ 2 ]*rE[ n ][ 2 ] + iE[ n ][ 2 ]*iE[ n ][ 2 ] ;
                
                    (*mpModel)( JOULE_HEATING, *pNodes[ n ] ) += 0.5 * eff_cond * modE2;
                    
                    mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
                }
            }
        }
	}

    //********************************************************************************************************************************
    //* - Calculates Joule heating on Gauss points.
    //********************************************************************************************************************************
	void Modeler::Calculate_Joule_Heating_OnGaussPoints( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& JOULE_HEATING_OnGP )
	{
        Vector2Type Joule_Heating_OnElement( mResultsOnGPs );

        for( int gp=0; gp<mResultsOnGPs; gp++ ) 
        {
            Joule_Heating_OnElement[ gp ].resize( 1 );
            Joule_Heating_OnElement[ gp ][ 0 ] = 0.0;
        }

        // Joule heating in IHL material with losses
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		double eo   = 8.8541878176e-12;
		double freq = mProblemFrequency;

		double sgm_real = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );
		double eps_imag = (*pProperties)( IMAG_ELECTRIC_PERMITTIVITY ) * eo;

		double eff_cond = freq * eps_imag + sgm_real;

        if( eff_cond > 0 )
        {
            Vector<Node::Pointer> pNodes;

            Get_Volume_Element_Nodes( NodesId, pNodes );

            Element::Pointer pElement;

            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cE;

            pElement->Calculate_E_field_OnGaussPoints( cE, mResultsOnGPs );

            VectorType rE( 3 ), iE( 3 );

            double modE2;

            for( int gp=0; gp<mResultsOnGPs; gp++ ) 
            {
                rE[ 0 ] = std::real( cE[ gp ][ 0 ] );
                rE[ 1 ] = std::real( cE[ gp ][ 1 ] );
                rE[ 2 ] = std::real( cE[ gp ][ 2 ] );

                iE[ 0 ] = std::imag( cE[ gp ][ 0 ] );
                iE[ 1 ] = std::imag( cE[ gp ][ 1 ] );
                iE[ 2 ] = std::imag( cE[ gp ][ 2 ] );

		    	modE2 = rE[ 0 ]*rE[ 0 ] + rE[ 1 ]*rE[ 1 ] + rE[ 2 ]*rE[ 2 ] + 
                        iE[ 0 ]*iE[ 0 ] + iE[ 1 ]*iE[ 1 ] + iE[ 2 ]*iE[ 2 ] ;
            
		    	Joule_Heating_OnElement[ gp ][ 0 ] = 0.5 * eff_cond * modE2;
            }
        }

        JOULE_HEATING_OnGP.push_back( Joule_Heating_OnElement );
    }

    //**********************************************************************************************************************
    //* - Calculates J induced on nodes ( IHL materials ).
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Induced_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
		Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		double sgm_real = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        double sgm_imag = (*pProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        std::complex<double> sigma( sgm_real, sgm_imag );

        if( std::abs( sigma ) > 0.0 ) 
        {
		    Vector<Node::Pointer> pNodes;

            Get_Volume_Element_Nodes( NodesId, pNodes );

            this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

            Element::Pointer pElement;
            
            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );
            
            Set_Volume_Element_Parameters( NodesId, pElement );

            if( mLL2PSmoothing )
            {
                cVector2Type cJ_OnGP;
            
                pElement->Calculate_E_field_OnGaussPoints( cJ_OnGP, 1 );

                cJ_OnGP[ 0 ] *= sigma;

                Vector2Type rJ_OnGP, iJ_OnGP;
            
                Split_Complex_VectorVector( rJ_OnGP, iJ_OnGP, cJ_OnGP );   

		        for( int n=0; n<pNodes.size(); n++ )
                {
                    (*mpModel)( REAL_J, *pNodes[ n ] ) += rJ_OnGP[ 0 ];
                    (*mpModel)( IMAG_J, *pNodes[ n ] ) += iJ_OnGP[ 0 ];
                
                    mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
                }    
            }
            else
            {
                cVector2Type cJ_OnNodes;
                
                pElement->Calculate_E_field_OnNodes( cJ_OnNodes );
                
                for( int n=0; n<pNodes.size(); n++ ) 
                {
                    cJ_OnNodes[ n ] *= sigma;
                }
                
                Vector2Type rJ_OnNodes, iJ_OnNodes;
                
                Split_Complex_VectorVector( rJ_OnNodes, iJ_OnNodes, cJ_OnNodes );
                
		        for( int n=0; n<pNodes.size(); n++ )
                {
                    (*mpModel)( REAL_J, *pNodes[ n ] ) += rJ_OnNodes[ n ];
                    (*mpModel)( IMAG_J, *pNodes[ n ] ) += iJ_OnNodes[ n ];
                
                    mNeighbourElements[ pNodes[ n ]->Id() - 1 ]++;
                }
            }
        }
    }

    //**********************************************************************************************************************
    //* - Calculates J induced on Gauss points ( IHL materials ).
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Induced_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_eddy_OnGP )
    {
 		cJ_eddy_OnGP.resize( mResultsOnGPs );

        for( int gp=0; gp<mResultsOnGPs; gp++ ) 
        {
            cJ_eddy_OnGP[ gp ].resize( 3 );
            cJ_eddy_OnGP[ gp ][ 0 ] = std::complex<double>( 0.0, 0.0 );
            cJ_eddy_OnGP[ gp ][ 1 ] = std::complex<double>( 0.0, 0.0 );
            cJ_eddy_OnGP[ gp ][ 2 ] = std::complex<double>( 0.0, 0.0 );
        }

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

		double sgm_real = (*pProperties)( REAL_ELECTRIC_CONDUCTIVITY );
        double sgm_imag = (*pProperties)( IMAG_ELECTRIC_CONDUCTIVITY );

        std::complex<double> sigma( sgm_real, sgm_imag );

        if( std::abs( sigma ) > 0.0 ) 
        {
		    Vector<Node::Pointer> pNodes;

            Get_Volume_Element_Nodes( NodesId, pNodes );

            Element::Pointer pElement;

            if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
		    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
            else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );

            Set_Volume_Element_Parameters( NodesId, pElement );

            cVector2Type cE;

            pElement->Calculate_E_field_OnGaussPoints( cE, mResultsOnGPs );

            for( int gp=0; gp<mResultsOnGPs; gp++ ) 
            {
                cJ_eddy_OnGP[ gp ][ 0 ] = sigma * cE[ gp ][ 0 ];
                cJ_eddy_OnGP[ gp ][ 1 ] = sigma * cE[ gp ][ 1 ];
                cJ_eddy_OnGP[ gp ][ 2 ] = sigma * cE[ gp ][ 2 ];
            }
        } 
    }

	//**********************************************************************************************************************
    //* - Calculates J imposed on nodes.
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Imposed_OnNodes( int* NodesId, unsigned int PropertiesId )
	{
        // Get volume element nodes
        Vector<int> HONodesId;

		Vector<Node::Pointer> pNodes;

        Get_Volume_Element_Nodes( NodesId, HONodesId, pNodes );

        // Get element properties
        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );

        this->SetProperties( PropertiesId, FREQUENCY, mProblemFrequency );

        // Volumetric source element definition
        Element::Pointer pElement;

        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new JSource_3sb_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  1 ) pElement = Element::Pointer( new JSource_1st_FullWave( pNodes, pProperties ) );
	    else if( mElementOrder ==  2 ) pElement = Element::Pointer( new JSource_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new JSource_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new JSource_2eg_FullWave( pNodes, pProperties ) );
	    
        // Activates potentials on element
        pElement->SetPotentials( mPotentials );

        // Get imposed J on element nodes
        pElement->Get_J_Imposed_OnNodes( tmp_rJ, tmp_iJ ); 
	}

	//**********************************************************************************************************************
    //* - Adds imported currents to element nodes.
    //**********************************************************************************************************************
	void Modeler::Calculate_Imported_Currents()
    {
        // Check if import currents mode is activated
        if( !mImportCurrents ) return;

        // Objectes from "dirent.h" 
        DIR *ExpJs_Dir;
        struct dirent *Dir_Ent;

        // Name of the folder where the exported J currents are stored
        String ExpJs_Dir_Name = "Export_J_Sources/";

        // If the directory "Export_J_Sources/" does not exist then exit function
        if( ( ExpJs_Dir = opendir( ExpJs_Dir_Name ) ) == NULL ) return;

        // Vector with imported cJ on nodes
        Vector< Vector< std::complex<double> > > Import_cJ_Nodal; 

        // Read all the files inside the directoy "Export_J_Sources/"
        while( ( Dir_Ent = readdir( ExpJs_Dir ) ) != NULL ) 
        {
            // If the entity is not a file go to the next entity
            if( Dir_Ent->d_type != DT_REG ) continue; 
                      
            // Saving file name as a string
            String File_Name( Dir_Ent->d_name );

            // Open J source file
            std::fstream ImpCurrentFile( ExpJs_Dir_Name + File_Name, std::fstream::in );  

            // Read first line
            std::string first_line; 
            
            ImpCurrentFile >> first_line;

            // If first character == 'N': file is read in nodal mode.
            // If first character == 'N': read J nodal file, store it in
            // mImported_cJ_Nodal, close file, and continue to next file.
            if( first_line[0] == 'N' )
            { 
                Read_Imported_Nodal_Current( ImpCurrentFile, Import_cJ_Nodal );
                continue;
            }

            // If first character != 'N': file is read in volumetric element mode.
            // In volumetric element mode, two options are available:
            // - is_a_global_phase_file == true : all the Jxyz in the file have the same phase.
            // - is_a_global_phase_file == false: each Jxyz in each element has its own phase.
            bool is_a_global_phase_file;

            // Global phase value
            double global_phase_value;
            
            // If first character == 'V', the file has a local phase for each element and J component.
            // If first character != 'V', it will by a number, which is the global phase of the file.
            if( first_line[0] == 'V' )
            { 
                is_a_global_phase_file = false; 
            }
            else                      
            { 
                is_a_global_phase_file = true; 
                global_phase_value     = std::stod( first_line );
            }

            // Element nodes Ids, Jxyz components, Jxyz phases
            Vector<int>    NodesId ( 4 );
            Vector<double> J_vector( 3 ); 
            Vector<double> J_phases( 3 );

            // Maps for volumetric file
            std::map<unsigned int, Vector<double> > rJ_file;
	        std::map<unsigned int, Vector<double> > iJ_file;
            std::map<unsigned int, int> numNeighbors;

            // Reading the rest of the file
            while( true )
            {
                // Reading Node Ids from file
                ImpCurrentFile >> NodesId[0]; 
                ImpCurrentFile >> NodesId[1];
                ImpCurrentFile >> NodesId[2];
                ImpCurrentFile >> NodesId[3];
                
                // Reading J source from file
                if( is_a_global_phase_file )
                {
                    ImpCurrentFile >> J_vector[0]; J_phases[0] = global_phase_value;
                    ImpCurrentFile >> J_vector[1]; J_phases[1] = global_phase_value;
                    ImpCurrentFile >> J_vector[2]; J_phases[2] = global_phase_value;
                }
                else
                {
                    ImpCurrentFile >> J_vector[0]; ImpCurrentFile >> J_phases[0];
                    ImpCurrentFile >> J_vector[1]; ImpCurrentFile >> J_phases[1];
                    ImpCurrentFile >> J_vector[2]; ImpCurrentFile >> J_phases[2];                
                }

                // Check if End Of File
                if( ImpCurrentFile.eof() ) break;

                // Projecting imported J on element nodes
                for( int i=0; i<NodesId.size(); i++ )
                {
                    int nodeId = NodesId[ i ];

                    // Check for the first occurrence of the node
                    if( rJ_file.find( nodeId ) != rJ_file.end() )
                    {
                        for( int j=0; j<3; j++ )
                        {
                            rJ_file[ nodeId ][ j ] += J_vector[ j ] * cos( J_phases[ j ] ); 
                            iJ_file[ nodeId ][ j ] += J_vector[ j ] * sin( J_phases[ j ] ); 
                        }

                        numNeighbors[ nodeId ]++;
                    }
                    else
                    {
                        Vector<double> rtmp( 3 );
                        Vector<double> itmp( 3 );
                        
                        for( int j=0; j<3; j++ )
                        {
                            rtmp[ j ] = J_vector[ j ] * cos( J_phases[ j ] ); 
                            itmp[ j ] = J_vector[ j ] * sin( J_phases[ j ] ); 
                        }

                        rJ_file[ nodeId ] = rtmp;
                        iJ_file[ nodeId ] = itmp;

                        numNeighbors[ nodeId ] = 1;
                    }
                }
            }

            // Close file
            ImpCurrentFile.close();

            // Assign J in file to global tmp_ist_rJ and tmp_ist_iJ
            std::map<unsigned int, int>::iterator nit;

            for( nit = numNeighbors.begin(); nit != numNeighbors.end(); nit++ )
            {
                int nodeId = nit->first;
                
                Vector<double> rtmp = rJ_file[ nodeId ]; 
                Vector<double> itmp = iJ_file[ nodeId ];

                for( int j=0; j<3; j++ )
                {
                    rtmp[ j ] /= nit->second; 
                    itmp[ j ] /= nit->second; 
                }

                if( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() )
                {
                    tmp_ist_rJ[ nodeId ] += rtmp; 
                    tmp_ist_iJ[ nodeId ] += itmp;                 
                }
                else
                {
                    tmp_ist_rJ[ nodeId ] = rtmp; 
                    tmp_ist_iJ[ nodeId ] = itmp;                 
                }
            }

            // Clear objects
            rJ_file     .clear(); std::map<unsigned int, Vector<double> >().swap( rJ_file      );
            iJ_file     .clear(); std::map<unsigned int, Vector<double> >().swap( iJ_file      );
            numNeighbors.clear(); std::map<unsigned int, int            >().swap( numNeighbors );
        }

        // Close directory Export_J_Sources/
        closedir( ExpJs_Dir );

        // Adding Import_cJ_Nodal to tmp_ist_rJ, tmp_ist_iJ
        if( Import_cJ_Nodal.size() > 0 )
        {
            for( int i=0; i<Import_cJ_Nodal.size(); i++ )
            {
                int nodeId = i + 1;
                
                Vector<double> rtmp( 3 );
                Vector<double> itmp( 3 );

                bool IsCero = true;

                for( int j=0; j<3; j++ )
                {
                    rtmp[ j ] = std::real( Import_cJ_Nodal[ i ][ j ] );
                    itmp[ j ] = std::imag( Import_cJ_Nodal[ i ][ j ] ); 

                    if( ( rtmp[ j ] != 0.0 ) || ( itmp[ j ] != 0.0 ) ) 
                    {
                        IsCero = false;
                    }
                }    

                if( IsCero ) continue;

                if( tmp_ist_rJ.find( nodeId ) != tmp_ist_rJ.end() )
                {
                    tmp_ist_rJ[ nodeId ] += rtmp;
                    tmp_ist_iJ[ nodeId ] += itmp;
                }
                else
                {
                    tmp_ist_rJ[ nodeId ] = rtmp;
                    tmp_ist_iJ[ nodeId ] = itmp;
                }
            }

            // Clear Import_cJ_Nodal
            Import_cJ_Nodal.clear(); 
            Vector< Vector< std::complex<double> > >().swap( Import_cJ_Nodal ); 
        }
	}

	//**********************************************************************************************************************
    //* - Calculates J imposed/imported on Gauss points.
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Imposed_OnGaussPoints( int* NodesId, unsigned int PropertiesId, cVector2Type& cJ_imps_OnGP )
    {
 		// cJ vector on Gauss points
        cJ_imps_OnGP.resize( mResultsOnGPs );

        // Initializes cJ vector on Gauss points
        for( int gp=0; gp<mResultsOnGPs; gp++ ) 
        {
            cJ_imps_OnGP[ gp ].resize( 3 );
            cJ_imps_OnGP[ gp ][ 0 ] = std::complex<double>( 0.0, 0.0 );
            cJ_imps_OnGP[ gp ][ 1 ] = std::complex<double>( 0.0, 0.0 );
            cJ_imps_OnGP[ gp ][ 2 ] = std::complex<double>( 0.0, 0.0 );
        }
        
        // Number of vertices of a tetrahedral element
        int numVerticeNodes = 4;

        // cJ vector on nodes
        cVector2Type cJ_OnNodes( numVerticeNodes );

        // Initializes cJ vector on nodes
        for( int n=0; n<numVerticeNodes; n++ ) 
        {
            cJ_OnNodes[ n ].resize( 3 );
            cJ_OnNodes[ n ][ 0 ] = std::complex<double>( 0.0, 0.0 );
            cJ_OnNodes[ n ][ 1 ] = std::complex<double>( 0.0, 0.0 );
            cJ_OnNodes[ n ][ 2 ] = std::complex<double>( 0.0, 0.0 );
        }

        // Check if the element is inside an imposed J source volume
        bool Element_Is_Inside_ImposedJ_Source;

        if( ( tmp_rJ.find( NodesId[ 0 ] ) != tmp_rJ.end() ) &&
            ( tmp_rJ.find( NodesId[ 1 ] ) != tmp_rJ.end() ) &&
            ( tmp_rJ.find( NodesId[ 2 ] ) != tmp_rJ.end() ) &&
            ( tmp_rJ.find( NodesId[ 3 ] ) != tmp_rJ.end() )  )
        {
            Element_Is_Inside_ImposedJ_Source = true;
        }
        else 
        {
            Element_Is_Inside_ImposedJ_Source = false;
        }

        // If the element is inside an imposed J source volume then add cJ values on vertices
        if( Element_Is_Inside_ImposedJ_Source )
        {
            for( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_OnNodes[ n ][ 0 ] += std::complex<double>( tmp_rJ[ NodesId[ n ] ][ 0 ], tmp_iJ[ NodesId[ n ] ][ 0 ] );
                cJ_OnNodes[ n ][ 1 ] += std::complex<double>( tmp_rJ[ NodesId[ n ] ][ 1 ], tmp_iJ[ NodesId[ n ] ][ 1 ] );
                cJ_OnNodes[ n ][ 2 ] += std::complex<double>( tmp_rJ[ NodesId[ n ] ][ 2 ], tmp_iJ[ NodesId[ n ] ][ 2 ] );
            }
        }

        // Check if the element is inside an imported J source volume and outside an imposed J source
        bool Element_Is_Inside_ImportedJ_Source;

        if( ( Element_Is_Inside_ImposedJ_Source == false          ) && 
            ( tmp_ist_rJ.find( NodesId[ 0 ] ) != tmp_ist_rJ.end() ) &&
            ( tmp_ist_rJ.find( NodesId[ 1 ] ) != tmp_ist_rJ.end() ) &&
            ( tmp_ist_rJ.find( NodesId[ 2 ] ) != tmp_ist_rJ.end() ) &&
            ( tmp_ist_rJ.find( NodesId[ 3 ] ) != tmp_ist_rJ.end() )  )
        {
            Element_Is_Inside_ImportedJ_Source = true;
        }
        else
        {
            Element_Is_Inside_ImportedJ_Source = false;
        }
         
        // If the element is inside an imported J source volume then add cJ values on vertices
        if( Element_Is_Inside_ImportedJ_Source )
        {
            for( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_OnNodes[ n ][ 0 ] += std::complex<double>( tmp_ist_rJ[ NodesId[ n ] ][ 0 ], tmp_ist_iJ[ NodesId[ n ] ][ 0 ] );
                cJ_OnNodes[ n ][ 1 ] += std::complex<double>( tmp_ist_rJ[ NodesId[ n ] ][ 1 ], tmp_ist_iJ[ NodesId[ n ] ][ 1 ] );
                cJ_OnNodes[ n ][ 2 ] += std::complex<double>( tmp_ist_rJ[ NodesId[ n ] ][ 2 ], tmp_ist_iJ[ NodesId[ n ] ][ 2 ] );
            }
        }

        // If the element is not inside an imposed or an imported J source then finish 
        if( ( Element_Is_Inside_ImposedJ_Source == false ) && ( Element_Is_Inside_ImportedJ_Source == false ) )
        {
            return;
        }

        // If the element is inside an imported/imposed J source then interpolate nodal values to Gauss points
        Vector<Node::Pointer> pNodes;
        
        Get_Volume_Element_Nodes( NodesId, pNodes );

        Properties::Pointer pProperties = mpModel->GetProperties( PropertiesId );
        
        Element::Pointer pElement;
        
        if     ( mElementOrder ==  0 ) pElement = Element::Pointer( new VolumeElement_3sb_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  1 ) pElement = Element::Pointer( new VolumeElement_1st_FullWave( pNodes, pProperties ) );
        else if( mElementOrder ==  2 ) pElement = Element::Pointer( new VolumeElement_2nd_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 11 ) pElement = Element::Pointer( new VolumeElement_1eg_FullWave( pNodes, pProperties ) );
        else if( mElementOrder == 12 ) pElement = Element::Pointer( new VolumeElement_2eg_FullWave( pNodes, pProperties ) );
        
        Vector<double> cX, cY, cZ; 
        
        pElement->GetInnerGiDGaussPoints( cX, cY, cZ, mResultsOnGPs );
        
        Matrix<double> N; 
        
        Lagrange3D_Ni_1st( N, cX, cY, cZ );

        // Interpolating nodal values to Gauss points
        for( int gp=0; gp<mResultsOnGPs; gp++ ) 
        {
            for( int n=0; n<numVerticeNodes; n++ )
            {
                cJ_imps_OnGP[ gp ][ 0 ] += N[ n ][ gp ] * cJ_OnNodes[ n ][ 0 ];
                cJ_imps_OnGP[ gp ][ 1 ] += N[ n ][ gp ] * cJ_OnNodes[ n ][ 1 ];
                cJ_imps_OnGP[ gp ][ 2 ] += N[ n ][ gp ] * cJ_OnNodes[ n ][ 2 ];
            }
        }
    }

	//**********************************************************************************************************************
    //* - Calculates J on elements.
    //**********************************************************************************************************************
	void Modeler::Calculate_J_Element( int* NodesId, unsigned int PropertiesId, ResultsOnGPsType& REAL_J_OnGP, ResultsOnGPsType& IMAG_J_OnGP )
	{
        Calculate_J_Induced_OnNodes( NodesId, PropertiesId );
        
        if( mResultsOnNodes ) 
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
            cJ_tot_OnGP[ gp ].resize( 3 );
            cJ_tot_OnGP[ gp ][ 0 ] = cJ_eddy_OnGP[ gp ][ 0 ] + cJ_imps_OnGP[ gp ][ 0 ];
            cJ_tot_OnGP[ gp ][ 1 ] = cJ_eddy_OnGP[ gp ][ 1 ] + cJ_imps_OnGP[ gp ][ 1 ];
            cJ_tot_OnGP[ gp ][ 2 ] = cJ_eddy_OnGP[ gp ][ 2 ] + cJ_imps_OnGP[ gp ][ 2 ];
        }

        Vector2Type rJ_tot_OnGP, iJ_tot_OnGP;
        
        Split_Complex_VectorVector( rJ_tot_OnGP, iJ_tot_OnGP, cJ_tot_OnGP );
        
        REAL_J_OnGP.push_back( rJ_tot_OnGP ); 
        IMAG_J_OnGP.push_back( iJ_tot_OnGP ); 
    }

	//**********************************************************************************************************************
    //* - Adds all the current contributions ( eddy, imposed and imported ).
    //**********************************************************************************************************************
	void Modeler::End_J_Calculation()
	{
        Vector<Node::Pointer> nodes = mpModel->GetNodesArray();

		Vector<Node::Pointer>::iterator nit;

		Vector<double> rEddy( 3, 0.0 ); 
        Vector<double> iEddy( 3, 0.0 );

		Vector<double> rJtot( 3, 0.0 ); 
        Vector<double> iJtot( 3, 0.0 );

        for( nit = nodes.begin(); nit != nodes.end(); ++nit )
        {
            int NodeId = (*nit)->Id();

            // Eddy currents contribution
			rEddy = (*mpModel)( REAL_J, **nit );
            iEddy = (*mpModel)( IMAG_J, **nit );

            if( mNeighbourElements[ NodeId-1 ] > 0 ) 
            {
                rEddy /= mNeighbourElements[ NodeId-1 ];
                iEddy /= mNeighbourElements[ NodeId-1 ];
            }
        
            for( int i=0; i<3; i++ )
            {
                rJtot[ i ] = rEddy[ i ];
                iJtot[ i ] = iEddy[ i ]; 
            }

            // Adding imposed currents contribution
			if( (              tmp_rJ.find( NodeId ) !=              tmp_rJ.end() ) && 
                ( mFix_FullWv_Voltage.find( NodeId ) == mFix_FullWv_Voltage.end() ) && 
                ( mFix_FullWv_Current.find( NodeId ) == mFix_FullWv_Current.end() )  )
			{
                for( int i=0; i<3; i++ )
                {
                    rJtot[ i ] += tmp_rJ[ NodeId ][ i ];
                    iJtot[ i ] += tmp_iJ[ NodeId ][ i ];
                }
			}

            // Adding imported electrostatic currents contribution
            if( (          tmp_ist_rJ.find( NodeId ) !=          tmp_ist_rJ.end() ) && 
                (              tmp_rJ.find( NodeId ) ==              tmp_rJ.end() ) && 
                ( mFix_FullWv_Voltage.find( NodeId ) == mFix_FullWv_Voltage.end() ) &&
                ( mFix_FullWv_Current.find( NodeId ) == mFix_FullWv_Current.end() )  )
			{
                for( int i=0; i<3; i++ )
                {
                    rJtot[ i ] += tmp_ist_rJ[ NodeId ][ i ];
                    iJtot[ i ] += tmp_ist_iJ[ NodeId ][ i ];
                }
			}

            // Total current distribution
            (*mpModel)( REAL_J, **nit ) = rJtot;
			(*mpModel)( IMAG_J, **nit ) = iJtot;
        }

		tmp_rJ.clear(); std::map< unsigned int, Vector<double> >().swap( tmp_rJ );
		tmp_iJ.clear(); std::map< unsigned int, Vector<double> >().swap( tmp_iJ );

        tmp_ist_rJ.clear(); std::map< unsigned int, Vector<double> >().swap( tmp_ist_rJ );
		tmp_ist_iJ.clear(); std::map< unsigned int, Vector<double> >().swap( tmp_ist_iJ );
	}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
	// - Projections
    //
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //*************************************************************************************************************
    //* - Calculates Sij parameters.
    //*************************************************************************************************************
    void Modeler::Calculate_Sij_Parameters()
	{
        // Calculate Sij only if there is 1 input port
        if( mInputPorts.size() != 1 ) 
        {
            return;
        }

        // Input port ID
        double InputPortID = mInputPorts.begin()->first;

        // Input port coefficient
        std::complex<double> InputPortCoeff = mInputPorts.begin()->second;
        
        // Iteration over all ports to calculate Sij
		std::map<double, std::complex<double> >::iterator pit;

		for( pit = mProjection.begin(); pit != mProjection.end(); pit++ )
		{
		    std::complex<double> NormCoeff = InputPortCoeff * mNormalization[ pit->first ];
        
            std::complex<double> NormProjection = pit->second / NormCoeff;

		    if( pit->first != InputPortID )
            {
			    mProjectionNormalized[ pit->first ] = NormProjection;
            }
		    else
            {
			    mProjectionNormalized[ pit->first ] = NormProjection - 1.0;
            }
        }
	}

    //*****************************************************************************************************************
    //* - Calculates of nodal magnitudes.
    //*****************************************************************************************************************
    void Modeler::Calculate_Nodal( int variableKey )
    {
		Vector<Node::Pointer> nodes = mpModel->GetPrintableNodesArray();

        Vector<Node::Pointer>::iterator it;

        // E(t)
        if( variableKey == E.getKey() )
        {
            Vector<double> E_total( 3 );
            Vector<double> E_real;
            Vector<double> E_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );

                E_total[0] = E_real[0] * cos( freq * time ) + E_imag[0] * sin( freq * time );                                                  
                E_total[1] = E_real[1] * cos( freq * time ) + E_imag[1] * sin( freq * time );                                                   
                E_total[2] = E_real[2] * cos( freq * time ) + E_imag[2] * sin( freq * time );

                (*mpModel)( E, **it ) = E_total;
            }
        }
		// B(t)
        else if( variableKey == B.getKey() )
        {
            Vector<double> B_total( 3 );
            Vector<double> B_real;
            Vector<double> B_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                B_real = (*mpModel)( REAL_B, **it );
                B_imag = (*mpModel)( IMAG_B, **it );

                B_total[0] = B_real[0] * cos( freq * time ) + B_imag[0] * sin( freq * time );
                B_total[1] = B_real[1] * cos( freq * time ) + B_imag[1] * sin( freq * time );
                B_total[2] = B_real[2] * cos( freq * time ) + B_imag[2] * sin( freq * time );

                (*mpModel)( B, **it ) = B_total;
            }
        }
		// H(t)
        else if( variableKey == H.getKey() )
        {
            Vector<double> H_total( 3 );
            Vector<double> H_real;
            Vector<double> H_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                H_real = (*mpModel)( REAL_H, **it );
                H_imag = (*mpModel)( IMAG_H, **it );

                H_total[0] = H_real[0] * cos( freq * time ) + H_imag[0] * sin( freq * time );
                H_total[1] = H_real[1] * cos( freq * time ) + H_imag[1] * sin( freq * time );
                H_total[2] = H_real[2] * cos( freq * time ) + H_imag[2] * sin( freq * time );

                (*mpModel)( H, **it ) = H_total;
            }
        }
		// J(t)
        else if( variableKey == J.getKey() )
        {
            Vector<double> J_total( 3 );
            Vector<double> J_real;
            Vector<double> J_imag;

            double time = mpModel->CurrentTime();
            double freq = mProblemFrequency;

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                J_real = (*mpModel)( REAL_J, **it );
                J_imag = (*mpModel)( IMAG_J, **it );

                J_total[0] = J_real[0] * cos( freq * time ) + J_imag[0] * sin( freq * time );
                J_total[1] = J_real[1] * cos( freq * time ) + J_imag[1] * sin( freq * time );
                J_total[2] = J_real[2] * cos( freq * time ) + J_imag[2] * sin( freq * time );

                (*mpModel)( J, **it ) = J_total;
            }
        }
		// MOD_E
        else if( variableKey == MOD_E.getKey() )
        {
		    Vector<double> E_real( 3 );
            Vector<double> E_imag( 3 );

		    for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                E_real = (*mpModel)( REAL_E, **it );
                E_imag = (*mpModel)( IMAG_E, **it );
               
                (*mpModel)( MOD_E, **it ) = sqrt( E_real[0]*E_real[0] + E_imag[0]*E_imag[0] +
                                                  E_real[1]*E_real[1] + E_imag[1]*E_imag[1] +
                                                  E_real[2]*E_real[2] + E_imag[2]*E_imag[2] );
            }
        }
		// MOD_H
        else if( variableKey == MOD_H.getKey() )
        {
		    Vector<double> H_real( 3 );
            Vector<double> H_imag( 3 );

		    for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                H_real = (*mpModel)( REAL_H, **it );
                H_imag = (*mpModel)( IMAG_H, **it );
               
                (*mpModel)( MOD_H, **it ) = sqrt( H_real[0]*H_real[0] + H_imag[0]*H_imag[0] +
                                                  H_real[1]*H_real[1] + H_imag[1]*H_imag[1] +
                                                  H_real[2]*H_real[2] + H_imag[2]*H_imag[2] );
            }
        }
		// MOD_J
        else if( variableKey == MOD_J.getKey() )
        {
		    Vector<double> J_real( 3 );
            Vector<double> J_imag( 3 );

		    for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                J_real = (*mpModel)( REAL_J, **it );
                J_imag = (*mpModel)( IMAG_J, **it );
               
                (*mpModel)( MOD_J, **it ) = sqrt( J_real[0]*J_real[0] + J_imag[0]*J_imag[0] +
                                                  J_real[1]*J_real[1] + J_imag[1]*J_imag[1] +
                                                  J_real[2]*J_real[2] + J_imag[2]*J_imag[2] );
            }
        }
		// MOD_B
        else if( variableKey == MOD_B.getKey() )
        {
		    Vector<double> B_real( 3 );
            Vector<double> B_imag( 3 );

		    for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                B_real = (*mpModel)( REAL_B, **it );
                B_imag = (*mpModel)( IMAG_B, **it );
               
                (*mpModel)( MOD_B, **it ) = sqrt( B_real[0]*B_real[0] + B_imag[0]*B_imag[0] +
                                                  B_real[1]*B_real[1] + B_imag[1]*B_imag[1] +
                                                  B_real[2]*B_real[2] + B_imag[2]*B_imag[2] );
            }
        }
		// BOUNDARY_NORMALS
        else if( variableKey == BOUNDARY_NORMALS.getKey() )
        {
            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
			    Vector<double> BoundaryNormal;

                Get_Boundary_Normal( BoundaryNormal, (*it)->Id() );

                (*mpModel)( BOUNDARY_NORMALS, **it ) = BoundaryNormal;
            }
        }
        // CONTACT_NORMALS
        else if( variableKey == CONTACT_NORMALS.getKey() )
        {
            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
			    Vector<double> ContactPNormal;

                Get_ContactP_Normal( ContactPNormal, (*it)->Id() );

                (*mpModel)( CONTACT_NORMALS, **it ) = ContactPNormal;
            }
        }
        // E field parallel to B_ext
        else if( variableKey == E_PARALLEL_T.getKey() )
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

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                // External magnetic flux density
                Bext = mpColdPlasma->Get_Bext_OnNode( *it );

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
        else if( variableKey == E_PERPENDICULAR_T.getKey() )
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

            for( it = nodes.begin(); it != nodes.end(); ++it )
            {
                // External magnetic flux density
                Bext = mpColdPlasma->Get_Bext_OnNode( *it );

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
        else if( variableKey == POYNTING_VECTOR.getKey() )
        {
            std::complex<double> cUnit( 0.00, 1.00 );

            Vector<double> E_real, E_imag;
            Vector<double> H_real, H_imag;
            Vector<double> avgS( 3 );

            Vector< std::complex<double> > cEv( 3 );
            Vector< std::complex<double> > cHv( 3 );

            Vector< std::complex<double> > E_X_conj_H(3);

            for( it = nodes.begin(); it != nodes.end(); ++it )
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
        else if( variableKey == LORENTZ_FORCE.getKey() )
        {
            std::complex<double> cUnit( 0.00, 1.00 );

            Vector<double> J_real, J_imag;
            Vector<double> B_real, B_imag;
            Vector<double> avgF( 3 );

            Vector< std::complex<double> > cJv( 3 );
            Vector< std::complex<double> > cBv( 3 );

            Vector< std::complex<double> > J_X_conj_B( 3 );

            for( it = nodes.begin(); it != nodes.end(); ++it )
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