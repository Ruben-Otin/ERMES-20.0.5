
#if !defined(KRATOS_SOURCES)
#define KRATOS_SOURCES

#include "../includes/matrix.h"
#include "../includes/sources_base.h"

namespace Kratos
{
    class Sources : public SourcesBase<double>, SourcesBase<Vector<double> >, SourcesBase<Matrix<double> >
    {
      public:
     
        typedef std::shared_ptr<Sources> Pointer; 
      
        Sources()
        {
        }

        virtual ~Sources()
        {
        }
    
        template<class TDataType> 
        TDataType operator()( const Variable<TDataType>& rV )
        {
            return ( SourcesBase<TDataType>::GetSource( rV ) )->getValue();
        }

        template<class TDataType> 
        TDataType operator()( const Variable<TDataType>& rV, const Point<double>& rThisPoint, double ThisTime )
        {
            return ( *SourcesBase<TDataType>::GetSource( rV ) )( rThisPoint );
        }
      
	    template<class TDataType> 
	    void AddSource( const Variable<TDataType>& rV, typename SourceFunction<TDataType>::Pointer pNewFunction )
        {
            SourcesBase<TDataType>::AddSource( rV, pNewFunction );
        }

	    template<class TDataType> 
	    bool HasSourceFor( const Variable<TDataType>& rVariable )
        {
            return SourcesBase<TDataType>::HasSourceFor( rVariable );
        }
    }; 
}  

#endif 

