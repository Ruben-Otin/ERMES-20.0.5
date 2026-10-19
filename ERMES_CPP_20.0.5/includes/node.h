
#if !defined(DEFINE_KRATOS_NODE)
#define DEFINE_KRATOS_NODE

#include <vector>
#include <set>

#include <memory>

#include "../includes/dof.h"
#include "../includes/compare_pdof.h"
#include "../includes/point.h"
#include "../includes/exception.h"
#include "../includes/point_data.h"

namespace Kratos
{
    class Element;
    
    class Node : public Point<double>
    {
      public:
      
        typedef std::shared_ptr<Node> Pointer;
        
        Node( int NodeId, double Xval, double Yval, double Zval ) : Point<double>( Xval, Yval, Zval )
        {
            mId = NodeId; 
        }
        
        void SetId( int NodeId )
        {
            mId = NodeId;
        }
        
        virtual ~Node()
        {
        }
        
        int Id() const 
        {
            return  mId;
        }
        
        unsigned int Index() const 
        {
            return mpData->Index();
        }
        
        void Add_cE_Dofs( std::set<Dof::Pointer,ComparePDof>& DofSet )
        {
            Dof::Pointer new_dofcEx( new Dof( cEx, Index() ) ); 
            Dof::Pointer new_dofcEy( new Dof( cEy, Index() ) ); 
            Dof::Pointer new_dofcEz( new Dof( cEz, Index() ) ); 
            
            mDofVector[0] = new_dofcEx;
            mDofVector[1] = new_dofcEy;
            mDofVector[2] = new_dofcEz;
        
            DofSet.insert( new_dofcEx );
            DofSet.insert( new_dofcEy );
            DofSet.insert( new_dofcEz );
        }

        void Add_cA_Dofs( std::set<Dof::Pointer,ComparePDof>& DofSet )
        {
            Dof::Pointer new_dofcAx( new Dof( cAx, Index() ) ); 
            Dof::Pointer new_dofcAy( new Dof( cAy, Index() ) ); 
            Dof::Pointer new_dofcAz( new Dof( cAz, Index() ) ); 
                                     
            mDofVector[0] = new_dofcAx;
            mDofVector[1] = new_dofcAy;
            mDofVector[2] = new_dofcAz;
            mDofVector[3] = NULL;
        
            DofSet.insert( new_dofcAx );
            DofSet.insert( new_dofcAy );
            DofSet.insert( new_dofcAz );
        }

        void Add_cAV_Dofs( std::set<Dof::Pointer,ComparePDof>& DofSet )
        {
            Dof::Pointer new_dofcAx( new Dof( cAx, Index() ) ); 
            Dof::Pointer new_dofcAy( new Dof( cAy, Index() ) ); 
            Dof::Pointer new_dofcAz( new Dof( cAz, Index() ) ); 
            Dof::Pointer new_dofcVs( new Dof( cVs, Index() ) ); 
                                     
            mDofVector[0] = new_dofcAx;
            mDofVector[1] = new_dofcAy;
            mDofVector[2] = new_dofcAz;
            mDofVector[3] = new_dofcVs;
        
            DofSet.insert( new_dofcAx );
            DofSet.insert( new_dofcAy );
            DofSet.insert( new_dofcAz );
            DofSet.insert( new_dofcVs );
        }

        void Add_V_Dof( std::set<Dof::Pointer,ComparePDof>& DofSet )
        {
            Dof::Pointer new_dof( new Dof( VOLTAGE, Index() ) ); 
        
            mDofVector[0] = new_dof;
        
            DofSet.insert( new_dof );
        }
        
        template<class TDataType>
        Dof::Pointer pDof( const Variable<TDataType>& rDofVariable )
        {
            int varKey = rDofVariable.getKey();
        
            if      ( varKey == cEx.getKey() ) { return mDofVector[0]; }
            else if ( varKey == cEy.getKey() ) { return mDofVector[1]; }
            else if ( varKey == cEz.getKey() ) { return mDofVector[2]; }
                                                                       
            else if ( varKey == cAx.getKey() ) { return mDofVector[0]; }
            else if ( varKey == cAy.getKey() ) { return mDofVector[1]; }
            else if ( varKey == cAz.getKey() ) { return mDofVector[2]; }
            else if ( varKey == cVs.getKey() ) { return mDofVector[3]; }
                                                                     
            else if ( varKey == VOLTAGE.getKey() ) { return mDofVector[0]; }

            // Unknown DOF variable ( previously the function fell off the end: undefined behaviour )
            return Dof::Pointer();
        }
        
        Dof::Pointer pDofV()
        {
            return mDofVector[0];
        }
        
        Dof::Pointer pDofcEx()
        {
            return mDofVector[0];
        }
        
        Dof::Pointer pDofcEy()
        {
            return mDofVector[1];
        }
        
        Dof::Pointer pDofcEz()
        {
            return mDofVector[2];
        }
        
        Dof::Pointer pDofcAx()
        {
            return mDofVector[0];
        }
        
        Dof::Pointer pDofcAy()
        {
            return mDofVector[1];
        }
        
        Dof::Pointer pDofcAz()
        {
            return mDofVector[2];
        }
        
        Dof::Pointer pDofcVs()
        {
            return mDofVector[3];
        }

        PointData::Pointer pData()
        {
            return mpData;
        }
        
        void SetData( PointData::Pointer pNewData )
        {
            mpData = pNewData;
        }
           
      private:
     
        int mId;
        
        PointData::Pointer mpData;
        
        Dof::Pointer mDofVector[4];
        
        // Default constructor as private prevent us from defining a node without node id
        Node();
        
        // Copy constructor as private help us to prevents bugs before define it explicitly
        Node( const Node& SourceNode );
        
        // Assignment operator as private help us to prevents bugs before define it explicitly
        const Node& operator=( const Node& SourceNode );   
    }; 
}  

#endif 

