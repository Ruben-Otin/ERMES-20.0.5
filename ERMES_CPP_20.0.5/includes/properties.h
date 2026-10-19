
#if !defined(KRATOS_PROPERTIES)
#define KRATOS_PROPERTIES

#include "../includes/properties_base.h"
#include "../includes/vector.h"

namespace Kratos
{
    class Properties : public PropertiesBase< double >              , 
                              PropertiesBase< Matrix<double> >      , 
                              PropertiesBase< std::complex<double> >, 
                              PropertiesBase< std::vector<double> > , 
                              PropertiesBase< Vector<double> >
    {
      public:
     
        typedef std::shared_ptr<Properties> Pointer; 

        Properties()
        {
        }
          
	    template<class TDefaultFunction> 
        Properties( TDefaultFunction& rDefaultFunction ) : PropertiesBase< double >              ( rDefaultFunction ),
                                                           PropertiesBase< Matrix<double> >      ( rDefaultFunction ),
                                                           PropertiesBase< std::complex<double> >( rDefaultFunction ), 
		                                                   PropertiesBase< std::vector<double> > ( rDefaultFunction ),
		                                                   PropertiesBase< Vector<double> >      ( rDefaultFunction )
        {                                                                                        
        }

	    virtual ~Properties()
        {
        }
      
        template<class TDataType> 
        TDataType operator()( const Variable<TDataType>& rV )
        {
            return PropertiesBase<TDataType>::GetProperty( rV );
        }

        template<class TDataType> 
        TDataType operator()( const Variable<TDataType>& rV, const Node& rThisNode, double ThisTime )
        {
            return PropertiesBase<TDataType>::GetProperty( rV, rThisNode, ThisTime );
        }
        
        template<class TDataType> 
	    TDataType operator()( const Variable<TDataType>& rV, const Node::Pointer pThisNode, double ThisTime )
        {
            return PropertiesBase<TDataType>::GetProperty( rV, *pThisNode, ThisTime );
        }
        
        template<class TDataType> 
        TDataType operator()( const Variable<TDataType>& rV, const Node& rThisNode )
        {
            return PropertiesBase<TDataType>::GetProperty( rV, rThisNode );
        }
        
        template<class TDataType> 
	    TDataType operator()( const Variable<TDataType>& rV, const Node::Pointer pThisNode )
        {
            return PropertiesBase<TDataType>::GetProperty( rV, *pThisNode );
        }
      
        template<class TDataType> 
        void SetProperty( const Variable<TDataType>& rV, typename PropertiesBase<TDataType>::FunctionType::Pointer pNewFunction )
        {
            PropertiesBase<TDataType>::SetProperty( rV, pNewFunction );
        }
    }; 
}  

#endif 

