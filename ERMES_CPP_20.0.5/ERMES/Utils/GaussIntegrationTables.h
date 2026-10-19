
#if !defined(KRATOS_GAUSS_INTEGRATION_TABLES)
#define KRATOS_GAUSS_INTEGRATION_TABLES

#include <stdio.h>    
#include "../../includes/vector.h"
   
namespace Kratos
{
    // GiD inner Gauss points - 0rd order
	int GaussPoints3D_InnerGiD_0rd( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	// GiD inner Gauss points - 1st order
	int GaussPoints3D_InnerGiD_1st( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );
	// GiD inner Gauss points - 2nd order
    int GaussPoints3D_InnerGiD_2nd( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ );

	// 1 pto de Gauss
	int GaussPoints3D_Order1 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
	// 8 ptos de Gauss		   																									 
	int GaussPoints3D_Order2 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 27 ptos de Gauss		   																									 
	int GaussPoints3D_Order3 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
	// 64 ptos de Gauss		   																									 
	int GaussPoints3D_Order4 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 125 ptos de Gauss		   																									 
	int GaussPoints3D_Order5 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 216 ptos de Gauss		   																									 
	int GaussPoints3D_Order6 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 343 ptos de Gauss		   																									 
	int GaussPoints3D_Order7 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 512 ptos de Gauss		   																									 
	int GaussPoints3D_Order8 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
    // 729 ptos de Gauss		   																									 
	int GaussPoints3D_Order9 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& cZ, Vector<double>& W );
							   
	// 1 pto de Gauss		   
	int GaussPoints2D_Order1 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 4 ptos de Gauss		   																		    
	int GaussPoints2D_Order2 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 9 ptos de Gauss		   																		    
	int GaussPoints2D_Order3 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 16 ptos de Gauss		   																		    
	int GaussPoints2D_Order4 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 25 ptos de Gauss		   																		    
	int GaussPoints2D_Order5 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 36 ptos de Gauss		   																		    
	int GaussPoints2D_Order6 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 49 ptos de Gauss		   																		    
	int GaussPoints2D_Order7 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 64 ptos de Gauss		   																		    
	int GaussPoints2D_Order8 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 81 ptos de Gauss		   																		    
	int GaussPoints2D_Order9 ( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 100 ptos de Gauss		   																		    
	int GaussPoints2D_Order10( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 121 ptos de Gauss		   																		    
	int GaussPoints2D_Order11( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
	// 144 ptos de Gauss		   																		    
	int GaussPoints2D_Order12( Vector<double>& cX, Vector<double>& cY, Vector<double>& W );
}

#endif