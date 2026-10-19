
#include "../includes/gid_output.h"
#include "../external_libraries/gidpost/source/gidpost.h"

namespace Kratos
{	
	GidOutput::GidOutput( Model::Pointer pModel )
	{
		mFirstWrite = false;
		mpModel     = pModel;
	}

	GidOutput::GidOutput( Model::Pointer pModel, String FileName ): mpModel( pModel ), mFileName( FileName )
	{
	}

	GidOutput::~GidOutput()
	{
		GiD_ClosePostResultFile();
	}

    //********************************************************************************************************************************    
	//* - Prints scalar results on nodes 
	//********************************************************************************************************************************
	void GidOutput::PrintOnNodes( const Variable<double>& rVariable )
	{
		double time = mpModel->CurrentTime();

		if( time <= mpModel->CurrentDeltaTime() ) 
		{
			std::cout << "Printing double " << rVariable.Name() << "..." << std::endl;
		}
			
		GiD_BeginResult( (char*)(const char*)( rVariable.Name() ), "ERMES", time, GiD_Scalar, GiD_OnNodes, NULL, NULL, 0, NULL );

		NodesArrayType nodes = mpModel->GetPrintableNodesArray();
		
		NodesArrayType::iterator it;

		for( it = nodes.begin(); it != nodes.end(); ++it )
		{
		    GiD_WriteScalar( (*it)->Id(), (*mpModel)( rVariable, **it ) );
		}

		GiD_EndResult();
	}

    //********************************************************************************************************************************    
	//* - Prints vector results on nodes 
	//********************************************************************************************************************************
	void GidOutput::PrintOnNodes( const Variable< Vector<double> >& rVariable )
	{
		double time = mpModel->CurrentTime();

		if( time <= mpModel->CurrentDeltaTime() ) 
		{
			std::cout << "Printing vector " << rVariable.Name() << "..." << std::endl;
		}

		GiD_BeginResult( (char*)(const char*)( rVariable.Name() ), "ERMES", time, GiD_Vector, GiD_OnNodes, NULL, NULL, 0, NULL );

		NodesArrayType nodes = mpModel->GetPrintableNodesArray();
		
		NodesArrayType::iterator it;

		int i;

        double gid_vec[ 3 ];

		for( it = nodes.begin(); it != nodes.end(); ++it )
		{
			Vector<double>& vec = (*mpModel)( rVariable, **it );

			if( vec.size() > 3 )
			{
			    for( i=0; i<3; i++ ) gid_vec[ i ] = vec[ i ];
			}
			else
			{
			    for( i=0; i<vec.size(); i++ ) gid_vec[ i ] = vec[ i ];
			    for(    ; i<3         ; i++ ) gid_vec[ i ] = 0.0;
			}

			GiD_WriteVector( (*it)->Id(), gid_vec[ 0 ], gid_vec[ 1 ], gid_vec[ 2 ] );
		}    

		GiD_EndResult();
	}

	//***************************************************************************************************************************    
    //* - Print scalar results on Gauss points 
    //***************************************************************************************************************************
	void GidOutput::PrintOnGaussPoints( const Variable< double >& rVariable, ResultsOnGPsType& ResultsOnGPsVector )
	{
        if( ResultsOnGPsVector.size() == 0 ) 
		{
			return;
		}
		
		double time = mpModel->CurrentTime();

		if( time <= mpModel->CurrentDeltaTime() ) 
		{
			std::cout << "Printing double " << rVariable.Name() << "..." << std::endl;
		}

		char* GaussPointsName;

		int numGaussPoints = ResultsOnGPsVector[ 0 ].size();

		if( numGaussPoints == 1 ) 
		{
			GaussPointsName = "GPsEM_01";
		}
		else 
		{
			GaussPointsName = "GPsEM_04";
		}

		GiD_BeginResult( (char*)(const char*)( rVariable.Name() ), "ERMES", time, GiD_Scalar, GiD_OnGaussPoints, GaussPointsName, NULL, 0, NULL );

        ResultsOnGPsType::iterator eit;

		Vector< Vector<double> >::iterator vit;

        int elementId = 1;

        for( eit = ResultsOnGPsVector.begin(); eit != ResultsOnGPsVector.end(); eit++ )
        {
			Vector< Vector<double> > valOnGPs = (*eit);

			for( vit = valOnGPs.begin(); vit != valOnGPs.end(); vit++ ) 
			{
				GiD_WriteScalar( elementId, (*vit)[ 0 ] );
			}
			
            elementId++;
        }    

        GiD_EndResult();    	
	}

	//***************************************************************************************************************************
    //* - Print vector results on Gauss points 
    //***************************************************************************************************************************
    void GidOutput::PrintOnGaussPoints( const Variable< Vector<double> >& rVariable, ResultsOnGPsType& ResultsOnGPsVector )
	{
        if( ResultsOnGPsVector.size() == 0 ) 
		{
			return;
		}
		
		double time = mpModel->CurrentTime();

		if( time <= mpModel->CurrentDeltaTime() ) 
		{
			std::cout << "Printing vector " << rVariable.Name() << "..." << std::endl;
		}

		char* GaussPointsName;

		int numGaussPoints = ResultsOnGPsVector[ 0 ].size();

		if( numGaussPoints == 1 ) 
		{
			GaussPointsName = "GPsEM_01";
		}
		else 
		{
			GaussPointsName = "GPsEM_04";
		}

		GiD_BeginResult( (char*)(const char*)( rVariable.Name() ), "ERMES", time, GiD_Vector, GiD_OnGaussPoints, GaussPointsName, NULL, 0, NULL );

		ResultsOnGPsType::iterator eit;

        Vector< Vector<double> >::iterator vit;

        int elementId = 1;
       
        for( eit = ResultsOnGPsVector.begin(); eit != ResultsOnGPsVector.end(); eit++ )
        {
            Vector< Vector<double> > valOnGPs = (*eit);

			for( vit = valOnGPs.begin(); vit != valOnGPs.end(); vit++ ) 
			{
				GiD_WriteVector( elementId, (*vit)[ 0 ], (*vit)[ 1 ], (*vit)[ 2 ] );
			}
			
			elementId++;
        }    

        GiD_EndResult();    	
	}

	//***************************************************************************************************************************    
    //* - Write new mesh with high order nodes
    //***************************************************************************************************************************
    void GidOutput::PrintMeshNodes( int ElementOrder )
    {
         std::cout << "Printing result mesh..." << std::endl;

		 int numNodes = 4;
        
		 GiD_BeginMesh( "HOMesh", GiD_3D, GiD_Tetrahedra, numNodes ); 

         GiD_BeginCoordinates();

		 NodesArrayType nodes = mpModel->GetPrintableNodesArray();

		 Model::NodesArrayType::iterator it;

         for( it = nodes.begin(); it != nodes.end(); ++it )
         {
             GiD_WriteCoordinates( (*it)->Id(), (*it)->X(), (*it)->Y(), (*it)->Z() );
         }

         GiD_EndCoordinates();

         GiD_BeginElements();
    }

    //***************************************************************************************************************************    
    //* - Print an element
    //***************************************************************************************************************************
    void GidOutput::PrintHOElement( int elementId, int* NodesIdMat )
    {  
		GiD_WriteElementMat( elementId, NodesIdMat );
    }

    //***************************************************************************************************************************    
    //* - Finish printing mesh
    //***************************************************************************************************************************
    void GidOutput::CloseHOMeshFile()
    {
         GiD_EndElements();
         GiD_EndMesh    ();

         std::cout << "Result mesh printed." << std::endl;
         std::cout << std::endl;
    }

	//***************************************************************************************************************************    
    //* - Prints header for Gauss points representation
    //***************************************************************************************************************************
    void GidOutput::PrintGaussPointsHeader()
    {
        GiD_BeginGaussPoint( "GPsEM_01", GiD_Tetrahedra, NULL, 1, 0, 1 ); GiD_EndGaussPoint();
		GiD_BeginGaussPoint( "GPsEM_04", GiD_Tetrahedra, NULL, 4, 0, 1 ); GiD_EndGaussPoint();
    }

	//***************************************************************************************************************************    
    //* - Setting the output file properties
	//***************************************************************************************************************************
	void  GidOutput::SetFileFormat( bool IsASCII )
	{
		mFirstWrite = true;
		
		// ASCII output file format  
		if( IsASCII )
		{   
			String RestFileName = mFileName;
			String MeshFileName = mFileName;

			RestFileName = RestFileName << ".post.res";

			if( GiD_OpenPostResultFile( (char*)(const char*)RestFileName, GiD_PostAscii ) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << RestFileName << "/";
				throw Exception( " GidOutput::GidOutput(String Filename)", std::runtime_error( buffer ) );
			}			
			
			MeshFileName = MeshFileName << ".post.msh";

			if( GiD_OpenPostMeshFile( (char*)(const char*)MeshFileName, GiD_PostAscii ) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << MeshFileName << "/";
				throw Exception( " GidOutput::GidOutput(String Filename)", std::runtime_error( buffer ) );
			}
		}
		// BINARY output file format (default)
		else
		{
			String RestFileName = mFileName;

			RestFileName = RestFileName << ".flavia.res";

			if( GiD_OpenPostResultFile( (char*)(const char*)RestFileName, GiD_PostBinary) )
			{
				String buffer;
				buffer << "error opening results file:" << "/" << RestFileName << "/";
				throw Exception( " GidOutput::GidOutput(String Filename)", std::runtime_error( buffer ) );
			}
	    }
	}   
} 



