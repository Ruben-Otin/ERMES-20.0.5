
#if !defined(KRATOS_GID_OUTPUT)
#define KRATOS_GID_OUTPUT

#include <fstream>

#include <memory>

#include "../includes/output.h"
#include "../includes/model.h"

namespace Kratos
{
    class GidOutput : public Output
	{
	  public:
      
	    typedef std::shared_ptr<GidOutput> Pointer; 

	    typedef std::vector<Node::Pointer> NodesArrayType;

		typedef Vector< Vector < Vector<double> > > ResultsOnGPsType;
	    
	    GidOutput( Model::Pointer pModel                  );
	    GidOutput( Model::Pointer pModel, String FileName );
	    
	    virtual ~GidOutput();
	    
	    void PrintOnNodes( const Variable<        double  >& rVariable );
	    void PrintOnNodes( const Variable< Vector<double> >& rVariable );
	    
        void PrintGaussPointsHeader();

		void PrintOnGaussPoints( const Variable<        double  >& rVariable, ResultsOnGPsType& ResultsOnGPsVector );
        void PrintOnGaussPoints( const Variable< Vector<double> >& rVariable, ResultsOnGPsType& ResultsOnGPsVector );
	  
	    void SetFileFormat ( bool IsASCII                   );      
	    void PrintHOElement( int elementId, int* NodesIdMat );
	    void PrintMeshNodes( int ElementOrder               );
	    
        void CloseHOMeshFile();

	  private:
  
		std::ofstream mResultsFile;
		  
		Model::Pointer mpModel;

		String mAnalisysName;
		String mFileName;

        bool mFirstWrite;
	}; 
} 

#endif 

