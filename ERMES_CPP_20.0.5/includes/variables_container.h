
#if !defined(KRATOS_VARIABLES_CONTAINER)
#define KRATOS_VARIABLES_CONTAINER

#include <memory>

#include "../includes/variables_container_base.h"
#include "../includes/matrix.h"

namespace Kratos
{
    class VariablesContainer : public VariablesContainerBase<int>                  , 
                                      VariablesContainerBase<double>               ,
                                      VariablesContainerBase<std::complex<double> >, 
                                      VariablesContainerBase<std::vector<double> > , 
                                      VariablesContainerBase<Vector<double> >      , 
                                      VariablesContainerBase<Matrix<double> >
    {
      public:
      
	    typedef std::shared_ptr<VariablesContainer> Pointer;
        
        VariablesContainer()
        {
        }
        
	    VariablesContainer( const VariablesContainer& rOther) : VariablesContainerBase<int>                  ( rOther ),
		                                                        VariablesContainerBase<double>               ( rOther ),
                                                                VariablesContainerBase<std::complex<double> >( rOther ),
		                                                        VariablesContainerBase<std::vector<double> > ( rOther ),
		                                                        VariablesContainerBase<Vector<double> >      ( rOther ),
		                                                        VariablesContainerBase<Matrix<double> >      ( rOther )
	    {
        }
        
        virtual ~VariablesContainer()
        {
        }
        
        template<class TDataType> 
        TDataType& operator()( const Variable<TDataType>& rV )
        {
            return VariablesContainerBase<TDataType>::GetValue( rV );
        }
        
        template<class TDataType> 
        const TDataType& operator()( const Variable<TDataType>& rV ) const
        {
            return VariablesContainerBase<TDataType>::GetValue( rV );
        }
        
        template<class TDataType> 
        TDataType& GetValue( const Variable<TDataType>& rV )
        {
            return VariablesContainerBase<TDataType>::GetValue( rV );
        }
        
        template<class TDataType> 
        const TDataType& GetValue( const Variable<TDataType>& rV ) const
        {
            return VariablesContainerBase<TDataType>::GetValue( rV );
        }              
   }; 
}  

#endif 

