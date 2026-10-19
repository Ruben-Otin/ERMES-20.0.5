
#if !defined(KRATOS_LAGRANGE_BASES)
#define KRATOS_LAGRANGE_BASES
     
#include <vector> 
#include "../../includes/matrix.h"      
   
namespace Kratos
{
	// Lagrange elements
	void Lagrange3D_Ni_1st( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_1st( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );
	    
	void Lagrange3D_Ni_2nd( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_2nd( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );
						    				   					   
	void Lagrange3D_Ni_3rd( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_3rd( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );
						    				   					   
	void Lagrange3D_Ni_4th( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_4th( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );	
											   					   
	// Lagrange elements + bubble elements	   					   
	void Lagrange3D_Ni_0sb( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
											   					   
	void Lagrange3D_Ni_1sb( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_1sb( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );
											   					   
	void Lagrange3D_Ni_3sb( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	void Lagrange2D_Ni_3sb( Matrix<double>& N, Vector<double>& cX, Vector<double>& cY                     );
}

#endif