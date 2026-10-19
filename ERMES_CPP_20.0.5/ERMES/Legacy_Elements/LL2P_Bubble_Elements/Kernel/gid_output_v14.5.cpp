
#include "../includes/gid_output.h"
#include "../external_libraries/gidpost/source/gidpost.h"

namespace Kratos
{	
	GidOutput::GidOutput( Model::Pointer pModel )
	{
		mFirstWrite = false;
		mpModel     = pModel;
	}

	GidOutput::GidOutput( Model::Pointer pModel, String FileName ) : mpModel( pModel ), mFileName( FileName )
	{
	}

	GidOutput::~GidOutput()
	{
		GiD_ClosePostResultFile();
	}

    //*****************************************************************    
	//* - Print nodal results on flavia.res file
	//*****************************************************************
	void GidOutput::Print( const Variable<double>& rVariable )
	{
		double time = mpModel->CurrentTime();

		if (time != 0.00)
		{
			if (time == mpModel->CurrentDeltaTime())
			{
				std::cout << "Printing double " << rVariable.Name() << "..." << std::endl;
			}

			GiD_BeginResult( (char*)(const char*)(rVariable.Name()),
							 "ERMES", 
							 time,
							 GiD_Scalar, 
							 GiD_OnNodes, NULL, NULL, 0, NULL );
		}
		else
		{
            std::cout << "Printing double " << rVariable.Name() << "..." << std::endl;

			GiD_BeginResult( (char*)(const char*)(rVariable.Name()),
							 "ERMES", 
							 0,
							 GiD_Scalar, 
							 GiD_OnNodes, NULL, NULL, 0, NULL );
		}

		// gets nodes on model;
		NodesArrayType nodes = mpModel->GetNodesArray();
		NodesArrayType::iterator it;

		for (it=nodes.begin();it!= nodes.end(); ++it)
		{
			 GiD_WriteScalar( (*it)->Id(), (*mpModel)(rVariable, **it) );
		}

		GiD_EndResult();
	}

    //*****************************************************************    
	//* - Print nodal results on flavia.res file
	//*****************************************************************
	void GidOutput::Print( const Variable< Vector<double> >& rVariable )
	{
		double gid_vec[3];
		int i;

		double time = mpModel->CurrentTime();

        if(time != 0.00)
		{
			if (time == mpModel->CurrentDeltaTime()) 
			{
                std::cout << "Printing vector " << rVariable.Name() << "..." << std::endl;
			}

			GiD_BeginResult( (char*)(const char*)(rVariable.Name()), 
							 "ERMES", 
							  time,
							  GiD_Vector, 
							  GiD_OnNodes, NULL, NULL, 0,NULL );
		}
		else
		{
            std::cout << "Printing vector " << rVariable.Name() << "..." << std::endl;

			GiD_BeginResult( (char*)(const char*)(rVariable.Name()), 
							 "ERMES", 
							 0,
							 GiD_Vector, 
							 GiD_OnNodes, NULL, NULL, 0,NULL );
		}

		// gets nodes on model;
		NodesArrayType nodes = mpModel->GetNodesArray();
		NodesArrayType::iterator it;

		for (it=nodes.begin();it!= nodes.end();++it)
		{
			Vector<double>& vec = (*mpModel)(rVariable, **it);

			if ( vec.size() > 3 )
			{
			  for ( i = 0; i < 3; i++ ) gid_vec[i] = vec[i];
			}
			else
			{
			  for (i=0; i < vec.size(); i++ ) gid_vec[i] = vec[i];
			  for (   ; i < 3         ; i++ ) gid_vec[i] = 0.0;
			}

			GiD_WriteVector( (*it)->Id(), gid_vec[0], gid_vec[1], gid_vec[2] );
		}    

		GiD_EndResult();
	}
  
    //*****************************************************************    
    //* - Print nodal results on flavia.res file
    //*****************************************************************
    void GidOutput::Print( const Variable<Matrix<double> >& rVariable )
    {
        int i;
        double time = mpModel->CurrentTime();

        if(time != 0.00)
        {
            if (time == mpModel->CurrentDeltaTime()) 
            {
                std::cout << "Printing matrix " << rVariable.Name() << "..." << std::endl;
            }

            GiD_BeginResult( (char *)(const char*)(rVariable.Name()), 
                             "ERMES", 
                             time,
                             GiD_Matrix, 
                             GiD_OnNodes, 
                             NULL, NULL, 0, NULL );
        }
        else
        {
            std::cout << "Printing matrix " << rVariable.Name() << "..." << std::endl;

            GiD_BeginResult( (char *)(const char*)(rVariable.Name()), 
                             "ERMES", 
                             0,
                             GiD_Matrix, 
                             GiD_OnNodes, 
                             NULL, NULL, 0, NULL );
        }

        // gets nodes on model;
        NodesArrayType nodes = mpModel->GetNodesArray();

        for (NodesArrayType::iterator it=nodes.begin(); it!= nodes.end();++it)
        {
            Matrix<double>  & mat = (*mpModel)(rVariable, **it);

            //for 2D problems (matrix size 1*3)
            if (mat.ColumnsNumber() == 3) GiD_Write2DMatrix( (*it)->Id(), mat(0,0), mat(0,1), mat(0,2) );
            //for 3D problems (matrix size 1*6)
            if (mat.ColumnsNumber() == 6) GiD_Write3DMatrix( (*it)->Id(), mat(0,0), mat(0,1), mat(0,2), mat(0,3), mat(0,4), mat(0,5) );
            //for plane deformation problems (matrix size 1*4)
            if (mat.ColumnsNumber() == 4) GiD_WritePlainDefMatrix( (*it)->Id(), mat(0,0), mat(0,1), mat(0,2), mat(0,3) );
        } 
        
        GiD_EndResult();
    }

    //*******************************************************************************************************************************    
    //* - Print vector results on Gauss points 
    //*******************************************************************************************************************************
    void GidOutput::PrintOnGaussPoints( const Variable< Vector<double> >& rVariable, Vector< Vector<double> >& ResultsOnGaussPoints )
    {
        std::cout << "Printing vector " << rVariable.Name() << "..." << std::endl;

        GiD_BeginResult( (char *)(const char*)(rVariable.Name()), "ERMES", 0, GiD_Vector, GiD_OnGaussPoints, "GPsEM", NULL, 0, NULL );

        int elementId = 1;
        
        Vector< Vector<double> >::iterator vit;
       
        for (vit = ResultsOnGaussPoints.begin(); vit != ResultsOnGaussPoints.end(); vit++)
        {
            GiD_WriteVector( elementId, (*vit)[0], (*vit)[1], (*vit)[2] );

            elementId++;
        }    

        GiD_EndResult();    
    }

    //********************************************************************************************************************    
    //* - Print scalar results on Gauss points 
    //********************************************************************************************************************
    void GidOutput::PrintOnGaussPoints( const Variable<double>& rVariable, Vector<double>& ResultsOnGaussPoints )
    {
        std::cout << "Printing double " << rVariable.Name() << "..." << std::endl;

        GiD_BeginResult( (char *)(const char*)(rVariable.Name()), "ERMES", 0, GiD_Scalar, GiD_OnGaussPoints, "GPsEM", NULL, 0, NULL );

        int elementId = 1;

        Vector<double>::iterator sit;

        for (sit = ResultsOnGaussPoints.begin(); sit != ResultsOnGaussPoints.end(); sit++)
        {
            GiD_WriteScalar( elementId, (*sit) );

            elementId++;
        }    

        GiD_EndResult();    
    }

	//********************************************************************************************************************    
    //* - Write new mesh with high order nodes
    //********************************************************************************************************************
    void GidOutput::PrintMeshNodes( int ElementOrder )
    {
         std::cout<<"Printing result mesh..."<<std::endl;
        
		 if   ( ElementOrder == 2 ) GiD_BeginMesh("HOMesh", GiD_3D, GiD_Tetrahedra, 10); 
		 else                       GiD_BeginMesh("HOMesh", GiD_3D, GiD_Tetrahedra, 4 ); 

         GiD_BeginCoordinates();

         Model::NodesArrayType::iterator node_it;

         for(node_it = mpModel->GetNodesArray().begin(); node_it != mpModel->GetNodesArray().end(); ++node_it)
         {
             GiD_WriteCoordinates( (*node_it)->Id(), (*node_it)->X(), (*node_it)->Y(), (*node_it)->Z() );
         }

         GiD_EndCoordinates();

         GiD_BeginElements();
    }

    //********************************************************************************************************************    
    //* - Print an element
    //********************************************************************************************************************
    void GidOutput::PrintHOElement( int elementId, int* NodesId )
    {  
		GiD_WriteElementMat( elementId, NodesId );
    }

    //********************************************************************************************************************    
    //* - Finish printing mesh
    //********************************************************************************************************************
    void GidOutput::CloseHOMeshFile()
    {
         GiD_EndElements();

         GiD_EndMesh();

         std::cout<<"Result mesh printed."<<std::endl;

         std::cout<<std::endl;
    }

    //********************************************************************************************************************    
    //* - Print header for Gauss points representation
    //********************************************************************************************************************
    void GidOutput::PrintGaussPointsHeader()
    {
        GiD_BeginGaussPoint( "GPsEM", GiD_Tetrahedra, NULL, 1, 0, 1 ); 

        GiD_EndGaussPoint();
    }

	/*
	//********************************************************************************************************************
	//* - Matrix output on gauss points. ATTENTION!!! -  to be called on a group of elements of the same type!!!
	//********************************************************************************************************************
	void GidOutput::PrintOnGaussPoints( const Variable<Matrix<double> >& rVariable, const String& rElementGroupName )
	{
		ProcessInfo& pCurrentProcessInfo = mpModel->GetCurrentProcessInfo();

		//Get Elements
		Model::ElementsArrayType& rElements = mpModel->GetElements(rElementGroupName);

		//defining the set of gauss points
		PElementArray::iterator first_element=rElements.begin();
		
		unsigned int NumberOfNodes       = (*first_element)->NumberOfNodes();
		unsigned int NumberOfGaussPoints = (*first_element)->GetGaussPointsVariables(rVariable,pCurrentProcessInfo).size();

		if (NumberOfNodes == 3 || NumberOfNodes ==6)
		{
		  GiD_BeginGaussPoint("element_gp", GiD_Triangle, NULL, NumberOfGaussPoints, 0, 1); //TO BE CHANGED!!
		}
		else if (NumberOfNodes == 4 || NumberOfNodes ==8 || NumberOfNodes ==9)
		{
		  GiD_BeginGaussPoint("element_gp", GiD_Quadrilateral, NULL, NumberOfGaussPoints, 0, 1); //TO BE CHANGED!!
		}

		GiD_EndGaussPoint();

		GiD_BeginResult( (char *)(const char*)(rVariable.Name()), 
		                 "Kratos", 
						 mpModel->GetSolutionStep(),
		                 GiD_Matrix, 
						 GiD_OnGaussPoints, 
						 "element_gp", NULL, 0, NULL );

        PElementArray::iterator it;

		for (it = rElements.begin(); it != rElements.end(); ++it)
		{
			//getting a Vector containing the values on all the gauss points
			Vector< Vector<double> > GaussDataVector = (*it)->GetGaussPointsVariables(rVariable,pCurrentProcessInfo);

			for(int i=0; i<GaussDataVector.size(); i++)
			{
				Vector<double> mat = GaussDataVector[i];

				//for 2D problems (matrix size 1*3)
				if (mat.size() == 3) GiD_Write2DMatrix( (*it)->Id(), mat[0], mat[1], mat[2] );

				//for 3D problems (matrix size 1*6)
				if (mat.size() == 6) GiD_Write3DMatrix( (*it)->Id(), mat[0], mat[1], mat[2], mat[3], mat[4], mat[5] );

				//for plane deformation problems (matrix size 1*4)
				if (mat.size() == 4) GiD_WritePlainDefMatrix( (*it)->Id(), mat[0], mat[1], mat[3] ,mat[2]);
			}    
		}

		GiD_EndResult();	  
	} 
	*/  

    /*
	//********************************************************************************************************************    
	//*
	//********************************************************************************************************************
	void GidOutput::PrintMesh( void )
	{
		KratosString filename = mFileName;
		filename <<  ".flavia.msh";

		GiD_OpenPostMeshFile((char *)(const char*)filename, GiD_PostAscii);

		GiD_BeginMesh("TestMsh",GiD_2D,GiD_Triangle,3);

		GiD_BeginCoordinates();

		Model::NodesArrayType::iterator node_iterator;

		for(node_iterator = mpModel->GetNodesArray().begin(); 
		    node_iterator!= mpModel->GetNodesArray().end() ; 
		    ++node_iterator)
		{
		    GiD_WriteCoordinates((*node_iterator)->Id(),
							     (*node_iterator)->X(),
							     (*node_iterator)->Y(),
							     (*node_iterator)->Z());
		}

		GiD_EndCoordinates();

		int nodes_id[3];

		GiD_BeginElements();

		Model::ElementsArrayType elements_array = mpModel->GetAllElements();

		Model::ElementsArrayType::iterator element_iterator:

		for(element_iterator = elements_array.begin(); 
		    element_iterator!= elements_array.end(); 
		    ++element_iterator)
		{
		    if( (*element_iterator)->GetNodes().size() == 3 )
		    {
				nodes_id[0] = (*element_iterator)->GetNodes()[0]->Id();
				nodes_id[1] = (*element_iterator)->GetNodes()[1]->Id();
				nodes_id[2] = (*element_iterator)->GetNodes()[2]->Id();

				GiD_WriteElement((*element_iterator)->Id(), nodes_id);
		    }
		}

		GiD_EndElements();

		GiD_EndMesh();

		GiD_ClosePostMeshFile();
	}
    */

	//********************************************************************************************************************    
    //* - Setting the output file properties
	//********************************************************************************************************************
	void  GidOutput::SetFileFormat( bool IsASCII )
	{
		mFirstWrite = true;
		
		// ASCII output file format  
		if ( IsASCII )
		{   
			String RestFileName = mFileName;
			String MeshFileName = mFileName;

			RestFileName = RestFileName << ".post.res";
			if ( GiD_OpenPostResultFile( (char*)(const char*)RestFileName, GiD_PostAscii ) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << RestFileName << "/";
				throw Exception(" GidOutput::GidOutput(String Filename)", std::runtime_error(buffer));
			}			
			
			MeshFileName = MeshFileName << ".post.msh";
			if ( GiD_OpenPostMeshFile( (char*)(const char*)MeshFileName, GiD_PostAscii ) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << MeshFileName << "/";
				throw Exception(" GidOutput::GidOutput(String Filename)", std::runtime_error(buffer));
			}
		}
		// BINARY output file format (default)
		else
		{
			String RestFileName = mFileName;

			RestFileName = RestFileName << ".flavia.res";
			if ( GiD_OpenPostResultFile( (char*)(const char*)RestFileName, GiD_PostBinary) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << RestFileName << "/";
				throw Exception(" GidOutput::GidOutput(String Filename)", std::runtime_error(buffer));
			}
	    }
	}
      
} // Namespace Kratos



