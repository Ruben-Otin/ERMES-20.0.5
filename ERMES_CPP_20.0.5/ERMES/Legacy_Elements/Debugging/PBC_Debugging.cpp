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
        
		if ( mPBCFrontElements.size() > 0 ) Set_FrontBack_PBC();
		if ( mPBCRightElements.size() > 0 ) Set_RightLeft_PBC();

        // Debugging 
        std::vector<Node::Pointer> nodes = mpModel->GetNodesArray();

        std::vector<Node::Pointer>::iterator it;

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int nId = (*it)->Id();
            
            if ( mPBC_NodeNodePairs_Front.find( nId ) != mPBC_NodeNodePairs_Front.end() )
            {
			    std::cout<<"nId : "<<nId<<" -> "<<mPBC_NodeNodePairs_Front[nId]<<std::endl;   
            }
			else if ( mPBC_NodeEdgePairs_Front.find( nId ) != mPBC_NodeEdgePairs_Front.end() )
			{
			    std::cout<<"nId : "<<nId<<" -> ";  

                std::vector<int>::iterator n2eit;

                for ( n2eit=mPBC_NodeEdgePairs_Front[nId].begin(); n2eit!=mPBC_NodeEdgePairs_Front[nId].end(); ++n2eit )
                {
                    std::cout<< *n2eit << " - ";  
                }
                std::cout<< " ||| "<< mPBC_NodeEdgeCoord_Front[nId]<<std::endl; 
			}
			else if ( mPBC_NodeElementPairs_Front.find( nId ) != mPBC_NodeElementPairs_Front.end() )
			{
			    std::cout<<"nId : "<<nId<<" -> ";  

                std::vector<int>::iterator n2eit;

                for ( n2eit=mPBC_NodeElementPairs_Front[nId].begin(); n2eit!=mPBC_NodeElementPairs_Front[nId].end(); ++n2eit )
                {
                    std::cout<< *n2eit << " + ";  
                }
                std::cout<< " ||| "<< mPBC_NodeElementCoord_Front[nId][0] <<" , "<< mPBC_NodeElementCoord_Front[nId][1] <<std::endl; 
			}
        }

        //std::cout << "---Debugging---" << " Front ****** Back " << std::endl;
        Send_Error_Msg( "---Debugging---", " Front ****** Back ", 1 );

        for ( it=nodes.begin(); it!=nodes.end(); ++it )
        {
            int nId = (*it)->Id();
            
            if ( mPBC_NodeNodePairs_Right.find( nId ) != mPBC_NodeNodePairs_Right.end() )
            {
			    std::cout<<"nId : "<<nId<<" -> "<<mPBC_NodeNodePairs_Right[nId]<<std::endl;   
            }
			else if ( mPBC_NodeEdgePairs_Right.find( nId ) != mPBC_NodeEdgePairs_Right.end() )
			{
			    std::cout<<"nId : "<<nId<<" -> ";  

                std::vector<int>::iterator n2eit;

                for ( n2eit=mPBC_NodeEdgePairs_Right[nId].begin(); n2eit!=mPBC_NodeEdgePairs_Right[nId].end(); ++n2eit )
                {
                    std::cout<< *n2eit << " - ";  
                }
                std::cout<< " ||| "<< mPBC_NodeEdgeCoord_Right[nId]<<std::endl; 
			}
			else if ( mPBC_NodeElementPairs_Right.find( nId ) != mPBC_NodeElementPairs_Right.end() )
			{
			    std::cout<<"nId : "<<nId<<" -> ";  

                std::vector<int>::iterator n2eit;

                for ( n2eit=mPBC_NodeElementPairs_Right[nId].begin(); n2eit!=mPBC_NodeElementPairs_Right[nId].end(); ++n2eit )
                {
                    std::cout<< *n2eit << " + ";  
                }
                std::cout<< " ||| "<< mPBC_NodeElementCoord_Right[nId][0] <<" , "<< mPBC_NodeElementCoord_Right[nId][1] <<std::endl; 
			}
        }

        //std::cout << "---Debugging---" << " Right ****** Left " << std::endl;
        Send_Error_Msg( "---Debugging---", " Right ****** Left ", 1 );
		
		////////////////////////////////////////////////////////////////////////////////////////
		std::list<int>::iterator it_nFl;
		
        std::cout<< "+++After remaining nodes, material: " << pbcMaterial << std::endl;
		
        for( it_nFl = FrontNodesList.begin(); it_nFl != FrontNodesList.end(); ++it_nFl ) 
        {
            std::cout<< (*it_nFl) << " , " << mPBC_NodeEdgePairs_Front[*it_nFl][0] 
                                  << " , " << mPBC_NodeEdgePairs_Front[*it_nFl][1] 
                                  << " , " << mPBC_NodeEdgeCoord_Front[*it_nFl] << std::endl;            
        }    
        std::cout << std::endl;

        std::cout<< "+++After remaining nodes, material: " << pbcMaterial << std::endl;
		
        for( it_nFl = FrontNodesList.begin(); it_nFl != FrontNodesList.end(); ++it_nFl ) 
        {
            std::cout<< (*it_nFl) << " , " ;            
        }   		
        std::cout << std::endl;	
	}