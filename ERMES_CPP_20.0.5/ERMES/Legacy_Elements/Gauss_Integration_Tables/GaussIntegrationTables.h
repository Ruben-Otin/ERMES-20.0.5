
#if !defined(KRATOS_GAUSS_INTEGRATION_TABLES)
#define  KRATOS_GAUSS_INTEGRATION_TABLES

#include <stdio.h>      
#include <vector>       
   
namespace Kratos
{
	//1 pto de Gauss
	int GaussPoints3D_Order1 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//8 ptos de Gauss
	int GaussPoints3D_Order2 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //27 ptos de Gauss
	int GaussPoints3D_Order3 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//64 ptos de Gauss
	int GaussPoints3D_Order4 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //125 ptos de Gauss
	int GaussPoints3D_Order5 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //216 ptos de Gauss
	int GaussPoints3D_Order6 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //343 ptos de Gauss
	int GaussPoints3D_Order7 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //512 ptos de Gauss
	int GaussPoints3D_Order8 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
    //729 ptos de Gauss
	int GaussPoints3D_Order9 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//1000 ptos de Gauss
	int GaussPoints3D_Order10(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//1331 ptos de Gauss
	int GaussPoints3D_Order11(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//1728 ptos de Gauss
	int GaussPoints3D_Order12(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	//2197 ptos de Gauss
	int GaussPoints3D_Order13(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& cZ, std::vector<double>& W);
	

	//1 pto de Gauss
	int GaussPoints2D_Order1 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//4 ptos de Gauss
	int GaussPoints2D_Order2 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//9 ptos de Gauss
	int GaussPoints2D_Order3 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//16 ptos de Gauss
	int GaussPoints2D_Order4 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//25 ptos de Gauss
	int GaussPoints2D_Order5 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//36 ptos de Gauss
	int GaussPoints2D_Order6 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//49 ptos de Gauss
	int GaussPoints2D_Order7 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//64 ptos de Gauss
	int GaussPoints2D_Order8 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//81 ptos de Gauss
	int GaussPoints2D_Order9 (std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//100 ptos de Gauss
	int GaussPoints2D_Order10(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//121 ptos de Gauss
	int GaussPoints2D_Order11(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//144 ptos de Gauss
	int GaussPoints2D_Order12(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
	//169 ptos de Gauss
	int GaussPoints2D_Order13(std::vector<double>& cX, std::vector<double>& cY, std::vector<double>& W);
}

#endif