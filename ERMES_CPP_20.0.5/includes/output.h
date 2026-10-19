
#if !defined(KRATOS_OUTPUT)
#define KRATOS_OUTPUT

#include <memory>
#include "../includes/variables.h"

namespace Kratos
{
    class Output
    {
 	  public:
                
        typedef std::shared_ptr<Output> Pointer;

        typedef Vector< Vector < Vector<double> > > ResultsOnGPsType;
 	    
 	    Output(){}
 	    
 	    virtual ~Output(){}
 	    
 	    virtual void SetFileFormat( bool bFormat ){}
 	    
 	    virtual void PrintOnNodes( const Variable<        double  >& rVariable ){}
 	    virtual void PrintOnNodes( const Variable< Vector<double> >& rVariable ){}
        virtual void PrintOnNodes( const Variable< Matrix<double> >& rVariable ){}

        virtual void PrintGaussPointsHeader(){}

        virtual void PrintOnGaussPoints( const Variable<        double  >& rVariable, ResultsOnGPsType& ResultsOnGPsVector ){}
        virtual void PrintOnGaussPoints( const Variable< Vector<double> >& rVariable, ResultsOnGPsType& ResultsOnGPsVector ){}
        virtual void PrintOnGaussPoints( const Variable< Matrix<double> >& rVariable, ResultsOnGPsType& ResultsOnGPsVector ){}

 	    virtual void PrintHOElement( int elementId, int* NodesIdMat ){}
 	    virtual void PrintMeshNodes( int ElementOrder               ){}
 	    
 	    virtual void CloseHOMeshFile(){}
    }; 
} 

#endif 

