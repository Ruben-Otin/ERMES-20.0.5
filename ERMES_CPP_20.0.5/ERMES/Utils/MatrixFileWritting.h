
#if !defined(KRATOS_MATRIX_FILE_WRITTING)
#define KRATOS_MATRIX_FILE_WRITTING

#include <complex>
#include <map>

#include "../../includes/vector.h"
#include "../../includes/csr_matrix.h"   

namespace Kratos
{
    typedef std::complex<double> ComplexType;

    typedef CSRMatrix<ComplexType> ComplexMatrixType;

    typedef std::vector<ComplexType> ComplexVectorType;

	int ReadInSolvedFile();

	void WriteInSolvedFile( int bIsSolved );

    void ReadInitialGuess  ( ComplexVectorType& Xo );
	bool ReadSolutionVector( ComplexVectorType& Xo );   // false if Vector_Xo.bin is missing or too short
	
    void WriteMatrixInFile   ( ComplexMatrixType& A );
    void WriteAuxMatrixInFile( ComplexMatrixType& A );
    void WriteVectorInFile   ( ComplexVectorType& b );

	void WriteResultsInFile    ( ComplexVectorType& Xr );
	void WriteStepResultsInFile( ComplexVectorType& Xr );

	void ReadInfoSolverFromFile    ();
	void ReadExternalPathFromFile  ( String baseFileName, String& ExternalPath );
    void ReadExpCurrentInfoFromFile( String baseFileName, String& ExpCurrent_FileName, double& ExpCurrent_Phase );

    void Write_Sij_OnFiles ( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort, double Freq );
    void Write_Sij_Headers ( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort              );
    void Print_Sij_OnScreen( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort, double Freq );

	void Write_VolIntgs_Headers ( std::map<double, Vector< std::complex<double> > >& VolIntg_E );
	void Write_VolIntgs_HeadersC( std::map<double, Vector< std::complex<double> > >& VolIntg_E );

	void Write_SrfIntgs_Headers ( std::map<double, Vector< std::complex<double> > >& SrfIntg_E );
	void Write_SrfIntgs_HeadersC( std::map<double, Vector< std::complex<double> > >& SrfIntg_E );

    void Write_VolIntgs_OnFiles ( std::map<double, Vector< std::complex<double> > >& VolIntg_E,
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
						          std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J,
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F,
                                  double Freq );

    void Write_VolIntgs_OnFiles ( std::map<double, Vector< std::complex<double> > >& VolIntg_E,
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
						          std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J,
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F,
                                  std::complex<double> cFreq );

	void Write_SrfIntgs_OnFiles ( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
		                          std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
		                          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          double Freq );

	void Write_SrfIntgs_OnFiles ( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
		                          std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
		                          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          std::complex<double> cFreq );

    void Print_VolIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& VolIntg_E, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F, 
						          std::map<double, double>& Volume, 
		                          double Freq );

    void Print_VolIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& VolIntg_E, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F, 
						          std::map<double, double>& Volume, 
		                          std::complex<double> cFreq );

	void Print_SrfIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          std::map<double, double>& SfArea, 
		                          double Freq );

	void Print_SrfIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          std::map<double, double>& SfArea, 
		                          std::complex<double> cFreq );

    void Print_VolIntg_Electrostatic_OnScreen( std::map<double, Vector<double> >& VolIntg_E, std::map<double, double>& Volume );
    void Print_SrfIntg_Electrostatic_OnScreen( std::map<double, Vector<double> >& SrfIntg_E, std::map<double, double>& SfArea );

	void Write_VolIntgs_Elements_OnFiles( std::map<int, Vector< Vector<int> > >& VolIntg_Elements );
	void Write_SrfIntgs_Elements_OnFiles( std::map<int, Vector< Vector<int> > >& SrfIntg_Elements );

	void Write_VolIntgs_Nodes_OnFile( Vector<int>& NodesIDs, Vector< Vector<double> >& NodesCoords, int VolumeID );
	void Write_SrfIntgs_Nodes_OnFile( Vector<int>& NodesIDs, Vector< Vector<double> >& NodesCoords, int SrfaceID );

	void Write_VolIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnVol, double Freq, int VolumeID );
	void Write_SrfIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnSrf, double Freq, int SrfaceID );

    void Write_VolIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnVol, std::complex<double> cFreq, int VolumeID );
	void Write_SrfIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnSrf, std::complex<double> cFreq, int SrfaceID );
}

#endif