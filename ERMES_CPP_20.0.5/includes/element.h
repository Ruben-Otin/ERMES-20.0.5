
#if !defined(DEFINE_KRATOS_ELEMENT)
#define DEFINE_KRATOS_ELEMENT

#include <fstream>
#include <vector>
#include <map>

#include <memory>

#include "../includes/define.h"
#include "../includes/dof.h"
#include "../includes/node.h"
#include "../includes/properties.h"
#include "../includes/element_sources.h"
#include "../includes/process_info.h"

#include "../ERMES/ColdPlasma.h"

namespace Kratos
{
    class Element
    {
      public:
             
        typedef std::shared_ptr<Element> Pointer;
        typedef std::vector<Node::Pointer> NodesArrayType;

        typedef Vector< double >                VectorType;
        typedef Vector< std::complex<double> > cVectorType;

        typedef Vector< Vector< double > >                Vector2Type;
        typedef Vector< Vector< std::complex<double> > > cVector2Type;

        typedef Vector< Vector< Vector< double > > > ResultsOnGPsType;

        Element( NodesArrayType& rNodes, Properties::Pointer rProperties ) : mNodes( rNodes ), mProperties( rProperties )
        {
        }

        Element( NodesArrayType& rNodes ) : mNodes( rNodes )
        {
        }

        // Virtual destructor: elements are polymorphic and are deleted through shared pointers
        virtual ~Element(){}

        virtual std::vector<Node::Pointer>& GetNodes( void )
        {
            return mNodes;
        }
        
        virtual int NumberOfNodes()
        {
	        return mNodes.size();
        }

        virtual int GetInnerGiDGaussPoints( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, int numResultsOnGPs )
        { 
            return int(); 
        }

	    virtual bool GetIsIHL()
	    {
	        return bool();
	    }

	    virtual double GetPortNumber()
	    {
	        return double();
	    }
        
	    virtual std::complex<double> GetProjection()
	    {
	        return std::complex<double>();
	    }
        
	    virtual std::complex<double> GetNormalization()
	    {
	        return std::complex<double>();
	    }

        virtual void GetEquationIdVector( int*                EqIdVector ){}
        virtual void GetEquationIdVector( Vector     < int >& EquationId ){}
        virtual void GetEquationIdVector( std::vector< int >& EquationId ){}
        
        virtual void GetStiffnessMatrix( Matrix< double               >& StiffMatrix ){}        
        virtual void GetStiffnessMatrix( Matrix< std::complex<double> >& StiffMatrix ){}
                                        
        virtual void GetResidualVector( std::complex<double>*           ResidualVector ){}        
        virtual void GetResidualVector( Vector< double               >& ResidualVector ){}
        virtual void GetResidualVector( Vector< std::complex<double> >& ResidualVector ){}

        virtual void GetResidualVector_Dirichlet( std::map< unsigned int, double               >& FixedValue, Matrix< double               >& EleStiffMatrix, Vector< double               >& ResidualVector ){}
        virtual void GetResidualVector_Dirichlet( std::map< unsigned int, std::complex<double> >& FixedValue, Matrix< std::complex<double> >& EleStiffMatrix, Vector< std::complex<double> >& ResidualVector ){}

        virtual void GetResidualVector_Flux( std::map< unsigned int, std::complex<double> >& FluxValue, Vector< std::complex<double> >& ResidualVector ){}
        virtual void GetResidualVector_Flux( Vector  < Vector<double>                     >& FluxValue, Vector< std::complex<double> >& ResidualVector ){}

        virtual void Get_J_Imposed_OnNodes( std::map< unsigned int, Vector<double> >& rJimp, std::map< unsigned int, Vector<double> >& iJimp ){}

        virtual void Calculate_Rotational_OnNodes( cVector2Type&  rot_OnNodes ){}
        virtual void Calculate_Gradient_OnNodes  (  Vector2Type& grad_OnNodes ){}
        virtual void Calculate_E_field_OnNodes   ( cVector2Type&   Ef_OnNodes ){}
 
        virtual void Calculate_Rotational_OnGaussPoints( cVector2Type&  rot_OnGPs, int numResultsOnGPs ){}
        virtual void Calculate_Gradient_OnGaussPoints  (  Vector2Type& grad_OnGPs, int numResultsOnGPs ){}   
        virtual void Calculate_E_field_OnGaussPoints   ( cVector2Type&   Ef_OnGPs, int numResultsOnGPs ){}

        virtual void Set_Cartesian_J_Nodal( Vector< Vector< std::complex<double> > >& Nodal_J ){}

        virtual void SetNormalType( std::map<unsigned int, char>& Type_Of_BCNormal ){}

        virtual void SetPlasmaModel( ColdPlasma::Pointer pColdPlasma ){}

        virtual void SetPotentials( bool PotentialsOn ){}
        virtual void SetStabilizer( bool StabilizerOn ){}

        virtual void SetFrequency( double NewFrequency ){}
        virtual void SetFrequency( std::complex<double> NewFrequency ){}

        virtual void SetPeso  ( double Peso   ){}
        virtual void SetFactor( double Factor ){}
      
      protected:
        
        Properties::Pointer mProperties;
        
        std::vector<Node::Pointer>& mNodes;
    }; 
    
    typedef std::vector<Element::Pointer> PElementArray;

    typedef std::map<const String, PElementArray> ElementsMap;  
}  

#endif 

