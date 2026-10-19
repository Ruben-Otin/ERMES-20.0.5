
#include <iostream>
#include <fstream>

#include "../../includes/kratos_string.h"
#include "../../includes/kernel.h"
#include <cstdlib>
#include <cstdio>

#include "MatrixFileWritting.h"

namespace Kratos
{
    //****************************************************************************************************************************************************
    //* - Write the system sparse matrix in a binary file 
    //****************************************************************************************************************************************************
    void WriteMatrixInFile( ComplexMatrixType& A )
    {
		std::cout << "Writing sparse matrix on file..." << std::endl;

		ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

		int row_s; 
        int  it_s;

		double val_real; 
        double val_imag;	

		std::fstream fouti( "Matrix_A_int.bin",std::ios::out | std::ios::binary );	
		
        if( !fouti.is_open() ) { std::cout << "ERMES error: cannot create system matrix file ( disk full or no write permission? )" << std::endl; std::cerr << "ERMES error: cannot create system matrix file" << std::endl; exit( EXIT_FAILURE ); }			

		std::fstream foutc( "Matrix_A_cmplx.bin",std::ios::out | std::ios::binary );	
		
        if( !foutc.is_open() ) { std::cout << "ERMES error: cannot create system matrix file ( disk full or no write permission? )" << std::endl; std::cerr << "ERMES error: cannot create system matrix file" << std::endl; exit( EXIT_FAILURE ); }			

		int row = 0;

        while( row_iterator != end_iterator )
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            while( data_iterator != end_data_iterator )
            {
				row_s = row+1;
				it_s  = data_iterator->first+1;

				val_real = std::real( data_iterator->second );
				val_imag = std::imag( data_iterator->second );

				fouti.write( (char*)&row_s, sizeof( int ) );			
				fouti.write( (char*)&it_s , sizeof( int ) );	

				foutc.write( (char*)&val_real, sizeof( double ) );	
				foutc.write( (char*)&val_imag, sizeof( double ) );						

                data_iterator++;
            }

            row_iterator++;
			row++;
        }

		fouti.close(); 
		foutc.close(); 

        std::cout << "Done." << std::endl << std::endl; 
    }

    //****************************************************************************************************************************************************
    //* - Write the system auxiliar sparse matrix in a binary file 
    //****************************************************************************************************************************************************
    void WriteAuxMatrixInFile( ComplexMatrixType& A )
    {
        std::cout << "Writing auxiliar sparse matrix on file..." << std::endl;

        ComplexMatrixType::DataArrayType::iterator row_iterator = A.Data().begin();
        ComplexMatrixType::DataArrayType::iterator end_iterator = A.Data().end();

        ComplexMatrixType::DataArrayType::value_type::iterator data_iterator;
        ComplexMatrixType::DataArrayType::value_type::iterator end_data_iterator;

        int row_s;	
        int  it_s;	

        double val_real;
        double val_imag;	

        std::fstream fouti( "Matrix_A_aux_int.bin",std::ios::out | std::ios::binary );	

        if( !fouti.is_open() ) { std::cout << "ERMES error: cannot create system matrix file ( disk full or no write permission? )" << std::endl; std::cerr << "ERMES error: cannot create system matrix file" << std::endl; exit( EXIT_FAILURE ); }			

        std::fstream foutc( "Matrix_A_aux_cmplx.bin",std::ios::out | std::ios::binary );	

        if( !foutc.is_open() ) { std::cout << "ERMES error: cannot create system matrix file ( disk full or no write permission? )" << std::endl; std::cerr << "ERMES error: cannot create system matrix file" << std::endl; exit( EXIT_FAILURE ); }			

        int row = 0;

        while( row_iterator != end_iterator )
        {
            data_iterator     = row_iterator->begin();
            end_data_iterator = row_iterator->end();

            while( data_iterator != end_data_iterator )
            {
                row_s = row+1;
                it_s  = data_iterator->first+1;

                val_real = std::real( data_iterator->second );
                val_imag = std::imag( data_iterator->second );

                fouti.write( (char*)&row_s, sizeof( int ) );			
                fouti.write( (char*)&it_s , sizeof( int ) );	

                foutc.write( (char*)&val_real, sizeof( double ) );	
                foutc.write( (char*)&val_imag, sizeof( double ) );						

                data_iterator++;
            }

            row_iterator++;
            row++;
        }

        fouti.close(); 
        foutc.close(); 

        std::cout << "Done." << std::endl << std::endl; 
    }

	//****************************************************************************************************************************************************
    //* - Write a vector on file "Vector_B.bin"
    //****************************************************************************************************************************************************
    void WriteVectorInFile( ComplexVectorType& b )
    {
		std::cout << "Writing dense vector on file..." << std::endl;
		
		std::fstream fout( "Vector_B.bin",std::ios::out | std::ios::binary );	

		if( !fout.is_open() ) return;		
        
        int size = b.size();

        for( int i=0; i<size; i++ ) 
        {
            fout.write( (char*)&b[ i ], sizeof( ComplexType ) );	
        }

		fout.close(); 

        std::cout << "Done." << std::endl << std::endl;
    }

	//****************************************************************************************************************************************************
	//* - Read solution vector from file "Vector_Xo.bin"
	//****************************************************************************************************************************************************
	bool ReadSolutionVector( ComplexVectorType& Xo )
	{
		// Returns false if "Vector_Xo.bin" does not exist or is shorter than the system size ( external solver failed ).
		// Previous versions returned silently and ERMES wrote a zero solution as if it were valid.
		std::fstream finput( "Vector_Xo.bin", std::ios::in | std::ios::binary );

		if( !finput.is_open() ) return false;
        
		std::cout << "Reading solution vector from file..." << std::endl;

		int size = Xo.size();

		for( int i = 0; i<size; i++ ) 
        {
            finput.read( (char*)&Xo[ i ], sizeof( ComplexType ) );

            if( finput.gcount() != (std::streamsize)sizeof( ComplexType ) )
            {
                finput.close();

                return false;
            }
        }

		finput.close();

		std::cout << "Done." << std::endl << std::endl;

		return true;
	}

	//****************************************************************************************************************************************************
    //* - Read initial guess from file "Vector_Xo.bin"
    //****************************************************************************************************************************************************
    void ReadInitialGuess( ComplexVectorType& Xo )
    {
		std::fstream finput( "Vector_Xo.bin", std::ios::in | std::ios::binary );
		
		if( !finput.is_open() ) return;	
        
		std::cout << "Reading initial guess from file..." << std::endl;

        int size = Xo.size();

        for( int i=0; i<size; i++ ) 
        {
            finput.read( (char*)&Xo[ i ], sizeof( ComplexType ) );	
        }

		finput.close(); 

        std::cout << "Done." << std::endl << std::endl;
    }

	//****************************************************************************************************************************************************
    //* - Write the results vector on file "Vector_Xo.bin"
    //****************************************************************************************************************************************************
    void WriteResultsInFile( ComplexVectorType& Xr )
    {
        std::fstream fout( "Vector_Xo.bin", std::ios::out | std::ios::binary );
		
		if( !fout.is_open() ) return;	
        
		std::cout << "Writing results vector on file..." << std::endl;	

        int size = Xr.size();

        for( int i=0; i<size; i++ ) 
        {
            fout.write( (char*)&Xr[ i ], sizeof( ComplexType ) );	
        }

		fout.close(); 

        std::cout << "Done." << std::endl << std::endl;
    }

	//****************************************************************************************************************************************************
    //* - Write the step results vector on file "Vector_Xo.bin"
    //****************************************************************************************************************************************************
    void WriteStepResultsInFile( ComplexVectorType& Xr )
    {
        std::fstream fout( "Vector_Xo.bin", std::ios::out | std::ios::binary );
		
		if( !fout.is_open() ) return;	
       
        int size = Xr.size();

        for( int i=0; i<size; i++ ) 
        {
            fout.write( (char*)&Xr[ i ], sizeof( ComplexType ) );	
        }

		fout.close(); 
    }

	//****************************************************************************************************************************************************
    //* - Write if the system is solved or not with MATLAB
    //****************************************************************************************************************************************************
    void WriteInSolvedFile( int bIsSolved )
    {
        std::fstream SolvedFile( "IsSolved.dat", std::ios::out | std::ios::binary );
    
        SolvedFile << bIsSolved << std::endl;

        SolvedFile.close();
    }

	//****************************************************************************************************************************************************
    //* - Read if the system is solved or not with MATLAB
    //****************************************************************************************************************************************************
    int ReadInSolvedFile()
    {
		int bIsSolved;

        std::fstream SolvedFile( "IsSolved.dat", std::fstream::in | std::ios::binary );  

        if( !SolvedFile )
		{
			bIsSolved = 0;
		}
		else
		{
			SolvedFile >> bIsSolved;

			SolvedFile.close();
		}
    
		return bIsSolved;
    }

	//****************************************************************************************************************************************************
    //* - Read info solver from file "info_solver.dat"
    //****************************************************************************************************************************************************
	void ReadInfoSolverFromFile()
    {
		std::fstream infoSolverFile( "info_solver.dat", std::fstream::in | std::ios::binary ); 

		String info;

		for( infoSolverFile >> info; !infoSolverFile.eof(); infoSolverFile >> info )
        {
			if( info != "ENDL" ) 
            {
                std::cout << info <<" ";
            }
			else   
            {
                std::cout << std::endl;
            }
		}

        infoSolverFile.close();

        std::cout << std::endl;
	}

	//****************************************************************************************************************************************************
    //* - Read from a file the external solver path.
    //****************************************************************************************************************************************************
	void ReadExternalPathFromFile( String baseFileName, String& ExternalPath )
    {   
		baseFileName << "-14.dat";
		
		std::fstream pathFile( baseFileName, std::fstream::in );

		String path;

		for( pathFile >> path; !pathFile.eof(); pathFile >> path )
		{
			ExternalPath << path << " ";
		}
		
        pathFile.close();
	}

    //****************************************************************************************************************************************************
    //* - Read from a file the export current filename Id and the phase of the exported current. 
    //****************************************************************************************************************************************************
	void ReadExpCurrentInfoFromFile( String baseFileName, String& ExpCurrent_FileName, double& ExpCurrent_Phase )
    {   
		baseFileName << "-15.dat";

        String ExpCurrent_Id;

		std::fstream ExpCurrentInfoFile( baseFileName, std::fstream::in );

        ExpCurrentInfoFile >> ExpCurrent_Id;
        ExpCurrentInfoFile >> ExpCurrent_Phase;

        ExpCurrentInfoFile.close();

        ExpCurrent_FileName = "ExportJSource-";
        ExpCurrent_FileName << ExpCurrent_Id << ".dat"; 
	}

    //****************************************************************************************************************************************************
    //* - Write Sij parameters files headers 
    //****************************************************************************************************************************************************
    void Write_Sij_Headers( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort )
    {
        // Print Sij files only if there is 1 input port
        if( InputPort.size() != 1 ) 
        {
            return;
        }

        // Input port ID
        int InputPort_ID = int( InputPort.begin()->first );
                
        // Iterate over ports
        std::map<double, ComplexType>::iterator pit;

        for( pit = TEProjection.begin(); pit != TEProjection.end(); pit++ )
		{
            String filename;

            int Port_ID = int( pit->first );

            filename << "Integrals/SParameter-S_" << Port_ID << "_" << InputPort_ID << ".dat";

            std::fstream save_Sij( filename, std::ios::out );

            save_Sij << "# Graf: " << "\"" << "SParameter-S" << Port_ID << InputPort_ID << "\""                                      << std::endl; 
            save_Sij << "#"                                                                                                          << std::endl; 
            save_Sij << "# X: \"Freq (Hz)\" Y: " << "\"" << "S" << Port_ID << InputPort_ID << " (dB)\"" << "    |Sij|    phase(rad)" << std::endl;
            save_Sij << "#"                                                                                                          << std::endl;

            save_Sij.close();
        }
    }

    //****************************************************************************************************************************************************
    //* - Write Sij parameters on file
    //****************************************************************************************************************************************************
    void Write_Sij_OnFiles( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort, double Freq )
    {
        // Print Sij files only if there is 1 input port
        if( InputPort.size() != 1 ) 
        {
            return;
        }

        // Input port ID
        int InputPort_ID = int( InputPort.begin()->first );
        
        // Problem frequency
        double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        // Write Sij parameters in file
        std::cout << "Writing scattering parameters on files..." << std::endl;

		std::map<double, ComplexType>::iterator pit;

        for( pit = TEProjection.begin(); pit!=TEProjection.end(); pit++ )
		{
            String filename;

            int Port_ID = int( pit->first );

            filename << "Integrals/SParameter-S_" << Port_ID << "_" << InputPort_ID << ".dat";

            std::fstream save_Sij( filename, std::ios::out | std::ios::app );

            save_Sij.setf( std::ios::scientific, std::ios::floatfield );

            double Mod_Proj = std::abs( pit->second );
            double Pha_Proj = std::arg( pit->second );

		    save_Sij << std::noshowpos << vf                   << "    " 
                     << std::showpos   << 20*log10( Mod_Proj ) << "    " 
                     << std::noshowpos << Mod_Proj             << "    " 
                     << std::showpos   << Pha_Proj             << std::endl; 

            save_Sij.close();
        }

        std::cout << "Done." << std::endl << std::endl;
    }

    //****************************************************************************************************************************************************
    //* - Prints Sij parameters on screen 
    //****************************************************************************************************************************************************
    void Print_Sij_OnScreen( std::map<double, ComplexType>& TEProjection, std::map<double, ComplexType>& InputPort, double Freq )
	{
        // Show Sij only if there is 1 input port
        if( InputPort.size() != 1 ) 
        {
            return;
        }

        // Input port ID
        int InputPort_ID = int( InputPort.begin()->first );

        // Problem frequency
        double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        std::cout << "------------------------------------------------------------" << std::endl;

		std::cout << "Frequency: " << vf << std::endl;

        std::cout << std::endl;
		
        // Show Sij parameters on screen
        std::map<double, ComplexType>::iterator pit;

		for( pit = TEProjection.begin(); pit != TEProjection.end(); pit++ )
		{
            int    Port_ID  = int     ( pit->first  );
            double Mod_Proj = std::abs( pit->second );
            double Pha_Proj = std::arg( pit->second );

            printf( "|S%i%i| = %+.6e dB [ %.6e, %+.6e rad ]\n", Port_ID, InputPort_ID, 20*log10( Mod_Proj ), Mod_Proj, Pha_Proj );
		}

        std::cout << "------------------------------------------------------------" << std::endl;

		std::cout << std::endl;
	}

    //****************************************************************************************************************************************************
    //* - Write headers for volume integral files 
    //****************************************************************************************************************************************************
    void Write_VolIntgs_Headers( std::map<double, Vector<std::complex<double> > >& VolIntg_E )
    {
        std::map<double, Vector<std::complex<double> > >::iterator vit;

        for( vit = VolIntg_E.begin(); vit!= VolIntg_E.end(); vit++ )
		{
            String filename;

            filename << "Integrals/Volume-" << vit->first << "-Integrals.dat";

            std::fstream save_Vi( filename, std::ios::out );

            save_Vi << "# Graf: " << "\"" << "Volume-" << vit->first << "-Integrals" << "\""                                         << std::endl; 
            save_Vi << "#"                                                                                                           << std::endl; 
            save_Vi << "# Freq(Hz)   |viEx|   phviEx(rad)   |viEy|   phviEy(rad)   |viEz|   phviEz(rad)   vi( |E| )   vi(  |E|^2  )" << std::endl;
			save_Vi << "# Freq(Hz)   |viHx|   phviHx(rad)   |viHy|   phviHy(rad)   |viHz|   phviHz(rad)   vi( |H| )   vi(  |H|^2  )" << std::endl;
            save_Vi << "# Freq(Hz)   |viBx|   phviBx(rad)   |viBy|   phviBy(rad)   |viBz|   phviBz(rad)   vi( |B| )   vi(  |B|^2  )" << std::endl;
			save_Vi << "# Freq(Hz)   |viJx|   phviJx(rad)   |viJy|   phviJy(rad)   |viJz|   phviJz(rad)   vi( |J| )   vi(  |J|^2  )" << std::endl;
            save_Vi << "# Freq(Hz)   |viFx|   phviFx(rad)   |viFy|   phvIFy(rad)   |viFz|   phviFz(rad)   vi( |F| )   vi( r(J*.E) )" << std::endl;
            save_Vi << "#"                                                                                                           << std::endl;

            save_Vi.close();
        }
    }

    //****************************************************************************************************************************************************
    //* - Write headers for volume integral files ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Write_VolIntgs_HeadersC( std::map<double, Vector<std::complex<double> > >& VolIntg_E )
    {
        std::map<double, Vector<std::complex<double> > >::iterator vit;

        for( vit = VolIntg_E.begin(); vit!= VolIntg_E.end(); vit++ )
		{
            String filename;

            filename << "Integrals/Volume-" << vit->first << "-Integrals.dat";

            std::fstream save_Vi( filename, std::ios::out );

            save_Vi << "# Graf: " << "\"" << "Volume-" << vit->first << "-Integrals" << "\""                                                        << std::endl; 
            save_Vi << "#"                                                                                                                          << std::endl; 
            save_Vi << "# Real(Freq)   Imag(Freq)   |viEx|   phviEx(rad)   |viEy|   phviEy(rad)   |viEz|   phviEz(rad)   vi( |E| )   vi(  |E|^2  )" << std::endl;
			save_Vi << "# Real(Freq)   Imag(Freq)   |viHx|   phviHx(rad)   |viHy|   phviHy(rad)   |viHz|   phviHz(rad)   vi( |H| )   vi(  |H|^2  )" << std::endl;
            save_Vi << "# Real(Freq)   Imag(Freq)   |viBx|   phviBx(rad)   |viBy|   phviBy(rad)   |viBz|   phviBz(rad)   vi( |B| )   vi(  |B|^2  )" << std::endl;
			save_Vi << "# Real(Freq)   Imag(Freq)   |viJx|   phviJx(rad)   |viJy|   phviJy(rad)   |viJz|   phviJz(rad)   vi( |J| )   vi(  |J|^2  )" << std::endl;
            save_Vi << "# Real(Freq)   Imag(Freq)   |viFx|   phviFx(rad)   |viFy|   phvIFy(rad)   |viFz|   phviFz(rad)   vi( |F| )   vi( r(J*.E) )" << std::endl;
            save_Vi << "#"                                                                                                                          << std::endl;

            save_Vi.close();
        }
    }

	//****************************************************************************************************************************************************
    //* - Write volume integrals in files 
    //****************************************************************************************************************************************************
    void Write_VolIntgs_OnFiles( std::map<double, Vector< std::complex<double> > >& VolIntg_E,
		                         std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
						         std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                         std::map<double, Vector< std::complex<double> > >& VolIntg_J,
                                 std::map<double, Vector< std::complex<double> > >& VolIntg_F,
                                 double Freq )
    {
        std::cout << "Writing volumes integrals on files..." << std::endl;

		double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        std::map<double, Vector< std::complex<double> > >::iterator vit;

		for( vit = VolIntg_E.begin(); vit != VolIntg_E.end(); vit++ )
		{
            String filename;
            
            filename << "Integrals/Volume-" << vit->first << "-Integrals.dat";

            std::fstream save_Vi( filename, std::ios::out | std::ios::app );

            save_Vi.setf( std::ios::scientific, std::ios::floatfield );

			Vector< std::complex<double> > viE = vit->second;
			Vector< std::complex<double> > viH = VolIntg_H[ vit->first ];
            Vector< std::complex<double> > viB = VolIntg_B[ vit->first ];
			Vector< std::complex<double> > viJ = VolIntg_J[ vit->first ];
            Vector< std::complex<double> > viF = VolIntg_F[ vit->first ];

		    save_Vi << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( viE[0] ) << "    " << std::showpos << std::arg( viE[0] ) << "    " 
                    << std::noshowpos << std::abs( viE[1] ) << "    " << std::showpos << std::arg( viE[1] ) << "    " 
                    << std::noshowpos << std::abs( viE[2] ) << "    " << std::showpos << std::arg( viE[2] ) << "    " 
                    << std::noshowpos << std::abs( viE[3] ) << "    " 
                    << std::noshowpos << std::abs( viE[4] ) << std::endl;     

		    save_Vi << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( viH[0] ) << "    " << std::showpos << std::arg( viH[0] ) << "    " 
                    << std::noshowpos << std::abs( viH[1] ) << "    " << std::showpos << std::arg( viH[1] ) << "    " 
                    << std::noshowpos << std::abs( viH[2] ) << "    " << std::showpos << std::arg( viH[2] ) << "    " 
                    << std::noshowpos << std::abs( viH[3] ) << "    " 
                    << std::noshowpos << std::abs( viH[4] ) << std::endl;   

		    save_Vi << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( viB[0] ) << "    " << std::showpos << std::arg( viB[0] ) << "    " 
                    << std::noshowpos << std::abs( viB[1] ) << "    " << std::showpos << std::arg( viB[1] ) << "    " 
                    << std::noshowpos << std::abs( viB[2] ) << "    " << std::showpos << std::arg( viB[2] ) << "    " 
                    << std::noshowpos << std::abs( viB[3] ) << "    " 
                    << std::noshowpos << std::abs( viB[4] ) << std::endl;   

		    save_Vi << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( viJ[0] ) << "    " << std::showpos << std::arg( viJ[0] ) << "    " 
                    << std::noshowpos << std::abs( viJ[1] ) << "    " << std::showpos << std::arg( viJ[1] ) << "    " 
                    << std::noshowpos << std::abs( viJ[2] ) << "    " << std::showpos << std::arg( viJ[2] ) << "    " 
                    << std::noshowpos << std::abs( viJ[3] ) << "    " 
                    << std::noshowpos << std::abs( viJ[4] ) << std::endl;   

		    save_Vi << std::noshowpos << vf                  << "    " 
                    << std::noshowpos << std::abs ( viF[0] ) << "    " << std::showpos << std::arg( viF[0] ) << "    " 
                    << std::noshowpos << std::abs ( viF[1] ) << "    " << std::showpos << std::arg( viF[1] ) << "    " 
                    << std::noshowpos << std::abs ( viF[2] ) << "    " << std::showpos << std::arg( viF[2] ) << "    " 
                    << std::noshowpos << std::abs ( viF[3] ) << "    " 
                    << std::noshowpos << std::real( viF[4] ) << std::endl;   
                                                                                  
            save_Vi.close();
        }

        std::cout << "Done." << std::endl << std::endl;
    }

	//****************************************************************************************************************************************************
    //* - Write volume integrals in files ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Write_VolIntgs_OnFiles( std::map<double, Vector< std::complex<double> > >& VolIntg_E,
		                         std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
						         std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                         std::map<double, Vector< std::complex<double> > >& VolIntg_J,
                                 std::map<double, Vector< std::complex<double> > >& VolIntg_F,
                                 std::complex<double> cFreq )
    {
        std::cout << "Writing volumes integrals on files..." << std::endl;

        double rFreq = std::real( cFreq );
        double iFreq = std::imag( cFreq );
        
        std::map<double, Vector< std::complex<double> > >::iterator vit;

		for( vit = VolIntg_E.begin(); vit != VolIntg_E.end(); vit++ )
		{
            String filename;
            
            filename << "Integrals/Volume-" << vit->first << "-Integrals.dat";

            std::fstream save_Vi( filename, std::ios::out | std::ios::app );

            save_Vi.setf( std::ios::scientific, std::ios::floatfield );

			Vector< std::complex<double> > viE = vit->second;
			Vector< std::complex<double> > viH = VolIntg_H[ vit->first ];
            Vector< std::complex<double> > viB = VolIntg_B[ vit->first ];
			Vector< std::complex<double> > viJ = VolIntg_J[ vit->first ];
            Vector< std::complex<double> > viF = VolIntg_F[ vit->first ];

		    save_Vi << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( viE[0] ) << "    " << std::showpos   << std::arg( viE[0] ) << "    " 
                    << std::noshowpos << std::abs( viE[1] ) << "    " << std::showpos   << std::arg( viE[1] ) << "    " 
                    << std::noshowpos << std::abs( viE[2] ) << "    " << std::showpos   << std::arg( viE[2] ) << "    " 
                    << std::noshowpos << std::abs( viE[3] ) << "    " 
                    << std::noshowpos << std::abs( viE[4] ) << std::endl;     

		    save_Vi << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( viH[0] ) << "    " << std::showpos   << std::arg( viH[0] ) << "    " 
                    << std::noshowpos << std::abs( viH[1] ) << "    " << std::showpos   << std::arg( viH[1] ) << "    " 
                    << std::noshowpos << std::abs( viH[2] ) << "    " << std::showpos   << std::arg( viH[2] ) << "    " 
                    << std::noshowpos << std::abs( viH[3] ) << "    " 
                    << std::noshowpos << std::abs( viH[4] ) << std::endl;   

		    save_Vi << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( viB[0] ) << "    " << std::showpos   << std::arg( viB[0] ) << "    " 
                    << std::noshowpos << std::abs( viB[1] ) << "    " << std::showpos   << std::arg( viB[1] ) << "    " 
                    << std::noshowpos << std::abs( viB[2] ) << "    " << std::showpos   << std::arg( viB[2] ) << "    " 
                    << std::noshowpos << std::abs( viB[3] ) << "    " 
                    << std::noshowpos << std::abs( viB[4] ) << std::endl;   

		    save_Vi << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( viJ[0] ) << "    " << std::showpos   << std::arg( viJ[0] ) << "    " 
                    << std::noshowpos << std::abs( viJ[1] ) << "    " << std::showpos   << std::arg( viJ[1] ) << "    " 
                    << std::noshowpos << std::abs( viJ[2] ) << "    " << std::showpos   << std::arg( viJ[2] ) << "    " 
                    << std::noshowpos << std::abs( viJ[3] ) << "    " 
                    << std::noshowpos << std::abs( viJ[4] ) << std::endl;   

		    save_Vi << std::noshowpos << rFreq               << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs ( viF[0] ) << "    " << std::showpos   << std::arg( viF[0] ) << "    " 
                    << std::noshowpos << std::abs ( viF[1] ) << "    " << std::showpos   << std::arg( viF[1] ) << "    " 
                    << std::noshowpos << std::abs ( viF[2] ) << "    " << std::showpos   << std::arg( viF[2] ) << "    " 
                    << std::noshowpos << std::abs ( viF[3] ) << "    " 
                    << std::noshowpos << std::real( viF[4] ) << std::endl;   

            save_Vi.close();
        }

        std::cout << "Done." << std::endl << std::endl;
    }

	//****************************************************************************************************************************************************
    //* - Prints volume integrals on screen 
    //****************************************************************************************************************************************************
    void Print_VolIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& VolIntg_E, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F, 
						          std::map<double, double>& Volume, 
                                  double Freq )
	{
        double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );
        
        std::map<double, Vector< std::complex<double> > >::iterator vit;

		for( vit = VolIntg_E.begin(); vit != VolIntg_E.end(); vit++ )
		{
			Vector< std::complex<double> > viE = vit->second;
			Vector< std::complex<double> > viH = VolIntg_H[ vit->first ];
            Vector< std::complex<double> > viB = VolIntg_B[ vit->first ];
			Vector< std::complex<double> > viJ = VolIntg_J[ vit->first ];
            Vector< std::complex<double> > viF = VolIntg_F[ vit->first ];

		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << "Volume[" << vit->first << "]: " << Volume[ vit->first ] << " m^3" << std::endl;

			std::cout << std::endl;

			std::cout << "Frequency: " << vf << " Hz" << std::endl;

            std::cout << std::endl;

			printf( "VolIntg( Ex ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[0] ), std::arg( viE[0] ) );
            printf( "VolIntg( Ey ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[1] ), std::arg( viE[1] ) );
            printf( "VolIntg( Ez ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[2] ), std::arg( viE[2] ) );
                                                             
			std::cout << std::endl;                          

            printf( "VolIntg( Hx ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[0] ), std::arg( viH[0] ) );
            printf( "VolIntg( Hy ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[1] ), std::arg( viH[1] ) );
            printf( "VolIntg( Hz ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[2] ), std::arg( viH[2] ) );
                                                                                                        
			std::cout << std::endl;                          
                                                             
            printf( "VolIntg( Bx ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[0] ), std::arg( viB[0] ) );
            printf( "VolIntg( By ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[1] ), std::arg( viB[1] ) );
            printf( "VolIntg( Bz ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[2] ), std::arg( viB[2] ) );
                                                                                                          
			std::cout << std::endl;                                                                                                                        
                                                                                                          
            printf( "VolIntg( Jx ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[0] ), std::arg( viJ[0] ) );
            printf( "VolIntg( Jy ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[1] ), std::arg( viJ[1] ) );
            printf( "VolIntg( Jz ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[2] ), std::arg( viJ[2] ) );
                                                            
			std::cout << std::endl;                         
                                                            
            printf( "VolIntg( |E| ) = %.6e\n", std::abs( viE[3] ) );
            printf( "VolIntg( |H| ) = %.6e\n", std::abs( viH[3] ) );
            printf( "VolIntg( |B| ) = %.6e\n", std::abs( viB[3] ) );
            printf( "VolIntg( |J| ) = %.6e\n", std::abs( viJ[3] ) );
			                                                
			std::cout << std::endl;                         
                                                            
            printf( "VolIntg( |E|^2 ) = %.6e\n", std::abs( viE[4] ) );
            printf( "VolIntg( |H|^2 ) = %.6e\n", std::abs( viH[4] ) );
            printf( "VolIntg( |B|^2 ) = %.6e\n", std::abs( viB[4] ) );
            printf( "VolIntg( |J|^2 ) = %.6e\n", std::abs( viJ[4] ) );

			std::cout << std::endl;                                                                                                                     
                                                                                                       
            printf( "VolIntg( (JxB)x ) = %+.6e + j * %+.6e\n", std::abs( viF[0] ), std::arg( viF[0] ) );
            printf( "VolIntg( (JxB)y ) = %+.6e + j * %+.6e\n", std::abs( viF[1] ), std::arg( viF[1] ) );
            printf( "VolIntg( (JxB)z ) = %+.6e + j * %+.6e\n", std::abs( viF[2] ), std::arg( viF[2] ) );
            printf( "VolIntg( |JxB|  ) = %+.6e            \n", std::abs( viF[3] )                     );

            std::cout << std::endl; 

            printf( "VolIntg(  E.J*  ) = %+.6e + j * %+.6e\n", std::real( viF[4] ), std::imag( viF[4] ) );
                                               
		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << std::endl;
		}
	}

	//****************************************************************************************************************************************************
    //* - Prints volume integrals on screen ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Print_VolIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& VolIntg_E, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_H, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_B, 
		                          std::map<double, Vector< std::complex<double> > >& VolIntg_J, 
                                  std::map<double, Vector< std::complex<double> > >& VolIntg_F, 
						          std::map<double, double>& Volume, 
                                  std::complex<double> cFreq )
	{
        double rFreq = std::real( cFreq );
        double iFreq = std::imag( cFreq );
        
        std::map<double, Vector< std::complex<double> > >::iterator vit;

		for( vit = VolIntg_E.begin(); vit != VolIntg_E.end(); vit++ )
		{
			Vector< std::complex<double> > viE = vit->second;
			Vector< std::complex<double> > viH = VolIntg_H[ vit->first ];
            Vector< std::complex<double> > viB = VolIntg_B[ vit->first ];
			Vector< std::complex<double> > viJ = VolIntg_J[ vit->first ];
            Vector< std::complex<double> > viF = VolIntg_F[ vit->first ];

		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << "Volume[" << vit->first << "]: " << Volume[ vit->first ] << " m^3" << std::endl;

			std::cout << std::endl;

			std::cout << "Frequency: " << rFreq << " + j * " << iFreq << std::endl;

            std::cout << std::endl;

			printf( "VolIntg( Ex ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[0] ), std::arg( viE[0] ) );
            printf( "VolIntg( Ey ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[1] ), std::arg( viE[1] ) );
            printf( "VolIntg( Ez ) = %.6e * exp( j * %+.6e )\n", std::abs( viE[2] ), std::arg( viE[2] ) );
                                                             
			std::cout << std::endl;                          

            printf( "VolIntg( Hx ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[0] ), std::arg( viH[0] ) );
            printf( "VolIntg( Hy ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[1] ), std::arg( viH[1] ) );
            printf( "VolIntg( Hz ) = %.6e * exp( j * %+.6e )\n", std::abs( viH[2] ), std::arg( viH[2] ) );
                                                                                                        
			std::cout << std::endl;                          
                                                             
            printf( "VolIntg( Bx ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[0] ), std::arg( viB[0] ) );
            printf( "VolIntg( By ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[1] ), std::arg( viB[1] ) );
            printf( "VolIntg( Bz ) = %.6e * exp( j * %+.6e )\n", std::abs( viB[2] ), std::arg( viB[2] ) );
                                                                                                          
			std::cout << std::endl;                                                                                                                        
                                                                                                          
            printf( "VolIntg( Jx ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[0] ), std::arg( viJ[0] ) );
            printf( "VolIntg( Jy ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[1] ), std::arg( viJ[1] ) );
            printf( "VolIntg( Jz ) = %.6e * exp( j * %+.6e )\n", std::abs( viJ[2] ), std::arg( viJ[2] ) );
                                                            
			std::cout << std::endl;                         
                                                            
            printf( "VolIntg( |E| ) = %.6e\n", std::abs( viE[3] ) );
            printf( "VolIntg( |H| ) = %.6e\n", std::abs( viH[3] ) );
            printf( "VolIntg( |B| ) = %.6e\n", std::abs( viB[3] ) );
            printf( "VolIntg( |J| ) = %.6e\n", std::abs( viJ[3] ) );
			                                                
			std::cout << std::endl;                         
                                                            
            printf( "VolIntg( |E|^2 ) = %.6e\n", std::abs( viE[4] ) );
            printf( "VolIntg( |H|^2 ) = %.6e\n", std::abs( viH[4] ) );
            printf( "VolIntg( |B|^2 ) = %.6e\n", std::abs( viB[4] ) );
            printf( "VolIntg( |J|^2 ) = %.6e\n", std::abs( viJ[4] ) );

			std::cout << std::endl;                                                                                                                     
                                                                                                       
            printf( "VolIntg( (JxB)x ) = %+.6e + j * %+.6e\n", std::abs( viF[0] ), std::arg( viF[0] ) );
            printf( "VolIntg( (JxB)y ) = %+.6e + j * %+.6e\n", std::abs( viF[1] ), std::arg( viF[1] ) );
            printf( "VolIntg( (JxB)z ) = %+.6e + j * %+.6e\n", std::abs( viF[2] ), std::arg( viF[2] ) );
            printf( "VolIntg( |JxB|  ) = %+.6e            \n", std::abs( viF[3] )                     );

            std::cout << std::endl; 

            printf( "VolIntg(  E.J*  ) = %+.6e + j * %+.6e\n", std::real( viF[4] ), std::imag( viF[4] ) );
                                               
		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << std::endl;
		}
	}

    //****************************************************************************************************************************************************
    //* - Prints volume integrals on screen 
    //****************************************************************************************************************************************************
    void Print_VolIntg_Electrostatic_OnScreen( std::map<double, Vector<double> >& VolIntg_E, std::map<double, double>& Volume )
    {
        std::map<double, Vector<double> >::iterator vit;

        for ( vit = VolIntg_E.begin(); vit != VolIntg_E.end(); vit++ )
        {
            Vector<double> viE = vit->second;
            
            std::cout << std::endl;
 
            std::cout << "------------------------------------------------------------" << std::endl;

            std::cout << "Volume[" << vit->first << "]: " << Volume[ vit->first ] << " m^3" << std::endl;

            std::cout << std::endl;
            
            printf( "VolIntg( Ex ) = %+.6e\n", viE[0] );
            printf( "VolIntg( Ey ) = %+.6e\n", viE[1] );
            printf( "VolIntg( Ez ) = %+.6e\n", viE[2] );
            
            std::cout << std::endl;
            
            printf( "VolIntg( |E|   ) = %.6e\n", viE[3] );
            printf( "VolIntg( |E|^2 ) = %.6e\n", viE[4] );
                          
            std::cout << "------------------------------------------------------------" << std::endl;
        }     
    }

    //****************************************************************************************************************************************************
    //* - Prints surface integrals on screen 
    //****************************************************************************************************************************************************
    void Print_SrfIntg_Electrostatic_OnScreen( std::map<double, Vector<double> >& SrfIntg_E, std::map< double, double>& SurfArea )
    {
        std::map<double, Vector<double> >::iterator sit;

        for( sit = SrfIntg_E.begin(); sit != SrfIntg_E.end(); sit++ )
        {
            Vector<double> siE = sit->second;
           
            std::cout << std::endl;

            std::cout << "------------------------------------------------------------" << std::endl;

            std::cout << "Surface[" << sit->first << "]: " << SurfArea[ sit->first ] << " m^2" << std::endl;

            std::cout << std::endl;
            
            printf( "SurfIntg( Ex ) = %+.6e\n", siE[0] );
            printf( "SurfIntg( Ey ) = %+.6e\n", siE[1] );
            printf( "SurfIntg( Ez ) = %+.6e\n", siE[2] );
            
            std::cout << std::endl;
            
            printf( "SurfIntg( |E|   ) = %.6e\n", siE[3] );
            printf( "SurfIntg( |E|^2 ) = %.6e\n", siE[4] );

            std::cout << "------------------------------------------------------------" << std::endl;
        }     
    }

    //****************************************************************************************************************************************************
    //* - Prints surface integrals on screen 
    //****************************************************************************************************************************************************
    void Print_SrfIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          std::map<double, double>& SfArea, 
                                  double Freq )
	{
		double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );
        
        std::map<double, Vector< std::complex<double> > >::iterator sit;

		for( sit = SrfIntg_E.begin(); sit != SrfIntg_E.end(); sit++ )
		{
			Vector< std::complex<double> > siE = sit->second;
			Vector< std::complex<double> > siH = SrfIntg_H[ sit->first ];
            Vector< std::complex<double> > siB = SrfIntg_B[ sit->first ];
			Vector< std::complex<double> > siJ = SrfIntg_J[ sit->first ];
            Vector< std::complex<double> > siS = SrfIntg_S[ sit->first ];

		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << "Surface[" << sit->first << "]: " << SfArea[ sit->first ] << " m^2" << std::endl;

			std::cout << std::endl;

			std::cout << "Frequency: " << vf << " Hz" << std::endl;

			std::cout << std::endl;

			printf( "SrfIntg( Ex ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[0] ), std::arg( siE[0] ) );
            printf( "SrfIntg( Ey ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[1] ), std::arg( siE[1] ) );
            printf( "SrfIntg( Ez ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[2] ), std::arg( siE[2] ) );
                                                             
			std::cout << std::endl;                          

            printf( "SrfIntg( Hx ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[0] ), std::arg( siH[0] ) );
            printf( "SrfIntg( Hy ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[1] ), std::arg( siH[1] ) );
            printf( "SrfIntg( Hz ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[2] ), std::arg( siH[2] ) );
                                                                                                        
			std::cout << std::endl;                          
                                                             
            printf( "SrfIntg( Bx ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[0] ), std::arg( siB[0] ) );
            printf( "SrfIntg( By ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[1] ), std::arg( siB[1] ) );
            printf( "SrfIntg( Bz ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[2] ), std::arg( siB[2] ) );
                                                                                                          
			std::cout << std::endl;                                                                                                                        
                                                                                                          
            printf( "SrfIntg( Jx ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[0] ), std::arg( siJ[0] ) );
            printf( "SrfIntg( Jy ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[1] ), std::arg( siJ[1] ) );
            printf( "SrfIntg( Jz ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[2] ), std::arg( siJ[2] ) );
            printf( "SrfIntg( Jn ) = %.6e * exp( j * %+.6e )\n", std::abs( siS[4] ), std::arg( siS[4] ) );
                                                            
			std::cout << std::endl;                         
                                                            
            printf( "SrfIntg( |E| ) = %.6e\n", std::abs( siE[3] ) );
            printf( "SrfIntg( |H| ) = %.6e\n", std::abs( siH[3] ) );
            printf( "SrfIntg( |B| ) = %.6e\n", std::abs( siB[3] ) );
            printf( "SrfIntg( |J| ) = %.6e\n", std::abs( siJ[3] ) );
			                                                
			std::cout << std::endl;                         
                                                            
            printf( "SrfIntg( |E|^2 ) = %.6e\n", std::abs( siE[4] ) );
            printf( "SrfIntg( |H|^2 ) = %.6e\n", std::abs( siH[4] ) );
            printf( "SrfIntg( |B|^2 ) = %.6e\n", std::abs( siB[4] ) );
            printf( "SrfIntg( |J|^2 ) = %.6e\n", std::abs( siJ[4] ) );

			std::cout << std::endl;                                                                                                                     
                                                                                                       
            printf( "SrfIntg( (ExH*)x ) = %+.6e + j * %+.6e\n", std::real( siS[0] ), std::imag( siS[0] ) );
            printf( "SrfIntg( (ExH*)y ) = %+.6e + j * %+.6e\n", std::real( siS[1] ), std::imag( siS[1] ) );
            printf( "SrfIntg( (ExH*)z ) = %+.6e + j * %+.6e\n", std::real( siS[2] ), std::imag( siS[2] ) );
            printf( "SrfIntg( (ExH*)n ) = %+.6e + j * %+.6e\n", std::real( siS[3] ), std::imag( siS[3] ) );
                                               
		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << std::endl;
		}
	}

    //****************************************************************************************************************************************************
    //* - Prints surface integrals on screen ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Print_SrfIntgs_OnScreen( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
								  std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
						          std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                          std::map<double, double>& SfArea, 
                                  std::complex<double> cFreq )
	{
        double rFreq = std::real( cFreq );
        double iFreq = std::imag( cFreq );

        std::map<double, Vector< std::complex<double> > >::iterator sit;

		for( sit = SrfIntg_E.begin(); sit != SrfIntg_E.end(); sit++ )
		{
			Vector< std::complex<double> > siE = sit->second;
			Vector< std::complex<double> > siH = SrfIntg_H[ sit->first ];
            Vector< std::complex<double> > siB = SrfIntg_B[ sit->first ];
			Vector< std::complex<double> > siJ = SrfIntg_J[ sit->first ];
            Vector< std::complex<double> > siS = SrfIntg_S[ sit->first ];

		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << "Surface[" << sit->first << "]: " << SfArea[ sit->first ] << " m^2" << std::endl;

			std::cout << std::endl;

			std::cout << "Frequency: " << rFreq << " + j * " << iFreq << std::endl;

			std::cout << std::endl;

			printf( "SrfIntg( Ex ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[0] ), std::arg( siE[0] ) );
            printf( "SrfIntg( Ey ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[1] ), std::arg( siE[1] ) );
            printf( "SrfIntg( Ez ) = %.6e * exp( j * %+.6e )\n", std::abs( siE[2] ), std::arg( siE[2] ) );
                                                             
			std::cout << std::endl;                          

            printf( "SrfIntg( Hx ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[0] ), std::arg( siH[0] ) );
            printf( "SrfIntg( Hy ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[1] ), std::arg( siH[1] ) );
            printf( "SrfIntg( Hz ) = %.6e * exp( j * %+.6e )\n", std::abs( siH[2] ), std::arg( siH[2] ) );
                                                                                                        
			std::cout << std::endl;                          
                                                             
            printf( "SrfIntg( Bx ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[0] ), std::arg( siB[0] ) );
            printf( "SrfIntg( By ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[1] ), std::arg( siB[1] ) );
            printf( "SrfIntg( Bz ) = %.6e * exp( j * %+.6e )\n", std::abs( siB[2] ), std::arg( siB[2] ) );
                                                                                                          
			std::cout << std::endl;                                                                                                                        
                                                                                                          
            printf( "SrfIntg( Jx ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[0] ), std::arg( siJ[0] ) );
            printf( "SrfIntg( Jy ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[1] ), std::arg( siJ[1] ) );
            printf( "SrfIntg( Jz ) = %.6e * exp( j * %+.6e )\n", std::abs( siJ[2] ), std::arg( siJ[2] ) );
            printf( "SrfIntg( Jn ) = %.6e * exp( j * %+.6e )\n", std::abs( siS[4] ), std::arg( siS[4] ) );
                                                            
			std::cout << std::endl;                         
                                                            
            printf( "SrfIntg( |E| ) = %.6e\n", std::abs( siE[3] ) );
            printf( "SrfIntg( |H| ) = %.6e\n", std::abs( siH[3] ) );
            printf( "SrfIntg( |B| ) = %.6e\n", std::abs( siB[3] ) );
            printf( "SrfIntg( |J| ) = %.6e\n", std::abs( siJ[3] ) );
			                                                
			std::cout << std::endl;                         
                                                            
            printf( "SrfIntg( |E|^2 ) = %.6e\n", std::abs( siE[4] ) );
            printf( "SrfIntg( |H|^2 ) = %.6e\n", std::abs( siH[4] ) );
            printf( "SrfIntg( |B|^2 ) = %.6e\n", std::abs( siB[4] ) );
            printf( "SrfIntg( |J|^2 ) = %.6e\n", std::abs( siJ[4] ) );

			std::cout << std::endl;                                                                                                                     
                                                                                                       
            printf( "SrfIntg( (ExH*)x ) = %+.6e + j * %+.6e\n", std::real( siS[0] ), std::imag( siS[0] ) );
            printf( "SrfIntg( (ExH*)y ) = %+.6e + j * %+.6e\n", std::real( siS[1] ), std::imag( siS[1] ) );
            printf( "SrfIntg( (ExH*)z ) = %+.6e + j * %+.6e\n", std::real( siS[2] ), std::imag( siS[2] ) );
            printf( "SrfIntg( (ExH*)n ) = %+.6e + j * %+.6e\n", std::real( siS[3] ), std::imag( siS[3] ) );
                                               
		    std::cout << "------------------------------------------------------------" << std::endl;

			std::cout << std::endl;
		}
	}

    //****************************************************************************************************************************************************
    //* - Write headers for surface integral file
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_Headers( std::map<double, Vector< std::complex<double> > >& SrfIntg_E )
    {
        std::map<double, Vector< std::complex<double> > >::iterator sit;

        for( sit = SrfIntg_E.begin(); sit!= SrfIntg_E.end(); sit++ )
		{
            String filename;

            filename << "Integrals/Surface-" << sit->first << "-Integrals.dat";

            std::fstream save_Si( filename, std::ios::out );

            save_Si << "# Graf: " << "\"" << "Surface-" << sit->first << "-Integrals" << "\""                                        << std::endl; 
            save_Si << "#"                                                                                                           << std::endl; 
            save_Si << "# Freq(Hz)   |siEx|   phsiEx(rad)   |siEy|   phsiEy(rad)   |siEz|    phsiEz(rad)   si( |E|  )   si( |E|^2 )" << std::endl;
			save_Si << "# Freq(Hz)   |siHx|   phsiHx(rad)   |siHy|   phsiHy(rad)   |siHz|    phsiHz(rad)   si( |H|  )   si( |H|^2 )" << std::endl;
            save_Si << "# Freq(Hz)   |siBx|   phsiBx(rad)   |siBy|   phsiBy(rad)   |siBz|    phsiBz(rad)   si( |B|  )   si( |B|^2 )" << std::endl;
			save_Si << "# Freq(Hz)   |siJx|   phsiJx(rad)   |siJy|   phsiJy(rad)   |siJz|    phsiJz(rad)   si( |J|  )   si( |J|^2 )" << std::endl;
            save_Si << "# Freq(Hz)   |siSx|   phsiSx(rad)   |siSy|   phsiSy(rad)   |siSz|    phsiSz(rad)   si( |Jn| )   si( r(Sn) )" << std::endl;
            save_Si << "#"                                                                                                           << std::endl;

            save_Si.close();
        }
    }

    //****************************************************************************************************************************************************
    //* - Write headers for surface integral file ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_HeadersC( std::map<double, Vector< std::complex<double> > >& SrfIntg_E )
    {
        std::map<double, Vector< std::complex<double> > >::iterator sit;

        for( sit = SrfIntg_E.begin(); sit!= SrfIntg_E.end(); sit++ )
		{
            String filename;

            filename << "Integrals/Surface-" << sit->first << "-Integrals.dat";

            std::fstream save_Si( filename, std::ios::out );

            save_Si << "# Graf: " << "\"" << "Surface-" << sit->first << "-Integrals" << "\""                                                      << std::endl; 
            save_Si << "#"                                                                                                                         << std::endl; 
            save_Si << "# Real(Freq)   Imag(Freq)   |siEx|   phsiEx(rad)   |siEy|   phsiEy(rad)   |siEz|   phsiEz(rad)   si( |E|  )   si( |E|^2 )" << std::endl;
			save_Si << "# Real(Freq)   Imag(Freq)   |siHx|   phsiHx(rad)   |siHy|   phsiHy(rad)   |siHz|   phsiHz(rad)   si( |H|  )   si( |H|^2 )" << std::endl;
            save_Si << "# Real(Freq)   Imag(Freq)   |siBx|   phsiBx(rad)   |siBy|   phsiBy(rad)   |siBz|   phsiBz(rad)   si( |B|  )   si( |B|^2 )" << std::endl;
			save_Si << "# Real(Freq)   Imag(Freq)   |siJx|   phsiJx(rad)   |siJy|   phsiJy(rad)   |siJz|   phsiJz(rad)   si( |J|  )   si( |J|^2 )" << std::endl;
            save_Si << "# Real(Freq)   Imag(Freq)   |siSx|   phsiSx(rad)   |siSy|   phsiSy(rad)   |siSz|   phsiSz(rad)   si( |Jn| )   si( r(Sn) )" << std::endl;
            save_Si << "#"                                                                                                                         << std::endl;

            save_Si.close();
        }
    }

    //****************************************************************************************************************************************************
    //* - Write surface integrals on file
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_OnFiles( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
						         std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
						         std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
		                         std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
		                         std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                         double Freq )
    {
        std::cout << "Writing surfaces integrals on files..." << std::endl;
        
        double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        std::map<double, Vector< std::complex<double> > >::iterator sit;

		for( sit = SrfIntg_E.begin(); sit != SrfIntg_E.end(); sit++ )
		{
            String filename;
            
            filename << "Integrals/Surface-" << sit->first << "-Integrals.dat";

            std::fstream save_Si( filename, std::ios::out | std::ios::app );

            save_Si.setf( std::ios::scientific, std::ios::floatfield );

			Vector< std::complex<double> > siE = sit->second;
			Vector< std::complex<double> > siH = SrfIntg_H[ sit->first ];
            Vector< std::complex<double> > siB = SrfIntg_B[ sit->first ];
			Vector< std::complex<double> > siJ = SrfIntg_J[ sit->first ];
            Vector< std::complex<double> > siS = SrfIntg_S[ sit->first ];

		    save_Si << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( siE[0] ) << "    " << std::showpos << std::arg( siE[0] ) << "    " 
                    << std::noshowpos << std::abs( siE[1] ) << "    " << std::showpos << std::arg( siE[1] ) << "    " 
                    << std::noshowpos << std::abs( siE[2] ) << "    " << std::showpos << std::arg( siE[2] ) << "    " 
                    << std::noshowpos << std::abs( siE[3] ) << "    " 
                    << std::noshowpos << std::abs( siE[4] ) << std::endl;     

		    save_Si << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( siH[0] ) << "    " << std::showpos << std::arg( siH[0] ) << "    " 
                    << std::noshowpos << std::abs( siH[1] ) << "    " << std::showpos << std::arg( siH[1] ) << "    " 
                    << std::noshowpos << std::abs( siH[2] ) << "    " << std::showpos << std::arg( siH[2] ) << "    " 
                    << std::noshowpos << std::abs( siH[3] ) << "    " 
                    << std::noshowpos << std::abs( siH[4] ) << std::endl;   

		    save_Si << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( siB[0] ) << "    " << std::showpos << std::arg( siB[0] ) << "    " 
                    << std::noshowpos << std::abs( siB[1] ) << "    " << std::showpos << std::arg( siB[1] ) << "    " 
                    << std::noshowpos << std::abs( siB[2] ) << "    " << std::showpos << std::arg( siB[2] ) << "    " 
                    << std::noshowpos << std::abs( siB[3] ) << "    " 
                    << std::noshowpos << std::abs( siB[4] ) << std::endl;   

		    save_Si << std::noshowpos << vf                 << "    " 
                    << std::noshowpos << std::abs( siJ[0] ) << "    " << std::showpos << std::arg( siJ[0] ) << "    " 
                    << std::noshowpos << std::abs( siJ[1] ) << "    " << std::showpos << std::arg( siJ[1] ) << "    " 
                    << std::noshowpos << std::abs( siJ[2] ) << "    " << std::showpos << std::arg( siJ[2] ) << "    " 
                    << std::noshowpos << std::abs( siJ[3] ) << "    " 
                    << std::noshowpos << std::abs( siJ[4] ) << std::endl;   

		    save_Si << std::noshowpos << vf                  << "    " 
                    << std::noshowpos << std::abs ( siS[0] ) << "    " << std::showpos << std::arg( siS[0] ) << "    " 
                    << std::noshowpos << std::abs ( siS[1] ) << "    " << std::showpos << std::arg( siS[1] ) << "    " 
                    << std::noshowpos << std::abs ( siS[2] ) << "    " << std::showpos << std::arg( siS[2] ) << "    " 
                    << std::noshowpos << std::abs ( siS[4] ) << "    " 
                    << std::noshowpos << std::real( siS[3] ) << std::endl;   

            save_Si.close();
        }

        std::cout << "Done." << std::endl << std::endl;
    }

    //****************************************************************************************************************************************************
    //* - Write surface integrals on file ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_OnFiles( std::map<double, Vector< std::complex<double> > >& SrfIntg_E, 
						         std::map<double, Vector< std::complex<double> > >& SrfIntg_H,
						         std::map<double, Vector< std::complex<double> > >& SrfIntg_B,
		                         std::map<double, Vector< std::complex<double> > >& SrfIntg_J,
		                         std::map<double, Vector< std::complex<double> > >& SrfIntg_S,
		                         std::complex<double> cFreq )
    {
        std::cout << "Writing surfaces integrals on files..." << std::endl;
        
        double rFreq = std::real( cFreq );
        double iFreq = std::imag( cFreq );

        std::map<double, Vector< std::complex<double> > >::iterator sit;

		for( sit = SrfIntg_E.begin(); sit != SrfIntg_E.end(); sit++ )
		{
            String filename;
            
            filename << "Integrals/Surface-" << sit->first << "-Integrals.dat";

            std::fstream save_Si( filename, std::ios::out | std::ios::app );

            save_Si.setf( std::ios::scientific, std::ios::floatfield );

			Vector< std::complex<double> > siE = sit->second;
			Vector< std::complex<double> > siH = SrfIntg_H[ sit->first ];
            Vector< std::complex<double> > siB = SrfIntg_B[ sit->first ];
			Vector< std::complex<double> > siJ = SrfIntg_J[ sit->first ];
            Vector< std::complex<double> > siS = SrfIntg_S[ sit->first ];

		    save_Si << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( siE[0] ) << "    " << std::showpos   << std::arg( siE[0] ) << "    " 
                    << std::noshowpos << std::abs( siE[1] ) << "    " << std::showpos   << std::arg( siE[1] ) << "    " 
                    << std::noshowpos << std::abs( siE[2] ) << "    " << std::showpos   << std::arg( siE[2] ) << "    " 
                    << std::noshowpos << std::abs( siE[3] ) << "    " 
                    << std::noshowpos << std::abs( siE[4] ) << std::endl;     

            save_Si << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( siH[0] ) << "    " << std::showpos   << std::arg( siH[0] ) << "    " 
                    << std::noshowpos << std::abs( siH[1] ) << "    " << std::showpos   << std::arg( siH[1] ) << "    " 
                    << std::noshowpos << std::abs( siH[2] ) << "    " << std::showpos   << std::arg( siH[2] ) << "    " 
                    << std::noshowpos << std::abs( siH[3] ) << "    " 
                    << std::noshowpos << std::abs( siH[4] ) << std::endl;   

		    save_Si << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( siB[0] ) << "    " << std::showpos   << std::arg( siB[0] ) << "    " 
                    << std::noshowpos << std::abs( siB[1] ) << "    " << std::showpos   << std::arg( siB[1] ) << "    " 
                    << std::noshowpos << std::abs( siB[2] ) << "    " << std::showpos   << std::arg( siB[2] ) << "    " 
                    << std::noshowpos << std::abs( siB[3] ) << "    " 
                    << std::noshowpos << std::abs( siB[4] ) << std::endl;   

		    save_Si << std::noshowpos << rFreq              << "    " << std::noshowpos << iFreq              << "    " 
                    << std::noshowpos << std::abs( siJ[0] ) << "    " << std::showpos   << std::arg( siJ[0] ) << "    " 
                    << std::noshowpos << std::abs( siJ[1] ) << "    " << std::showpos   << std::arg( siJ[1] ) << "    " 
                    << std::noshowpos << std::abs( siJ[2] ) << "    " << std::showpos   << std::arg( siJ[2] ) << "    " 
                    << std::noshowpos << std::abs( siJ[3] ) << "    " 
                    << std::noshowpos << std::abs( siJ[4] ) << std::endl;   

		    save_Si << std::noshowpos << rFreq               << "    " << std::noshowpos << iFreq              << "    "  
                    << std::noshowpos << std::abs ( siS[0] ) << "    " << std::showpos   << std::arg( siS[0] ) << "    " 
                    << std::noshowpos << std::abs ( siS[1] ) << "    " << std::showpos   << std::arg( siS[1] ) << "    " 
                    << std::noshowpos << std::abs ( siS[2] ) << "    " << std::showpos   << std::arg( siS[2] ) << "    " 
                    << std::noshowpos << std::abs ( siS[4] ) << "    " 
                    << std::noshowpos << std::real( siS[3] ) << std::endl;   

            save_Si.close();
        }

        std::cout << "Done." << std::endl << std::endl;
    }

    //****************************************************************************************************************************************************
    //* - Writes volume integral elements on files
    //****************************************************************************************************************************************************
    void Write_VolIntgs_Elements_OnFiles( std::map<int, Vector< Vector<int> > >& VolIntg_Elements )
    {
        std::map<int, Vector< Vector<int> > >::iterator vit;

		for( vit = VolIntg_Elements.begin(); vit != VolIntg_Elements.end(); vit++ )
		{
            String ElementsFileName;
            
            ElementsFileName << "Fields/Volume-" << vit->first << "-Elements.dat";

            std::fstream ElementsFile( ElementsFileName, std::ios::out );

            ElementsFile.setf( std::ios::scientific, std::ios::floatfield );

			Vector< Vector<int> > ElementsInVol = vit->second;

            Vector< Vector<int> >::iterator eit;

            for( eit = ElementsInVol.begin(); eit != ElementsInVol.end(); eit++ )
            {
                Vector<int> Element = *eit;
                
                for( int i=0; i<Element.size(); i++ ) 
                {
                    ElementsFile << std::noshowpos << Element[i] << "    "; 
                }
               
                ElementsFile << std::endl; 
            }

            ElementsFile.close();
        }
    }

    //****************************************************************************************************************************************************
    //* - Writes nodes IDs and coordinates of a volume on a file
    //****************************************************************************************************************************************************
    void Write_VolIntgs_Nodes_OnFile( Vector<int>& NodesIDs, Vector< Vector<double> >& NodesCoords, int VolumeID )
    {
        String NodesFileName;
        
        NodesFileName << "Fields/Volume-" << VolumeID << "-Nodes.dat";

        std::fstream NodesFile( NodesFileName, std::ios::out );

        NodesFile.setf( std::ios::scientific, std::ios::floatfield );
        
        Vector< Vector<double> >::iterator cit = NodesCoords.begin();

        Vector<int>::iterator nit;

		for( nit = NodesIDs.begin(); nit != NodesIDs.end(); nit++, cit++ )
		{
            NodesFile << std::noshowpos << (*nit)    << "    " 
                      << std::showpos   << (*cit)[0] << "    " 
                      << std::showpos   << (*cit)[1] << "    " 
                      << std::showpos   << (*cit)[2] << std::endl;        
        }

        NodesFile.close();
    }

    //****************************************************************************************************************************************************
    //* - Writes surface integral elements on files
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_Elements_OnFiles( std::map<int, Vector< Vector<int> > >& SrfIntg_Elements )
    {
        std::map<int, Vector< Vector<int> > >::iterator sit;

		for( sit = SrfIntg_Elements.begin(); sit != SrfIntg_Elements.end(); sit++ )
		{
            String ElementsFileName;
            
            ElementsFileName << "Fields/Surface-" << sit->first << "-Elements.dat";

            std::fstream ElementsFile( ElementsFileName, std::ios::out );

            ElementsFile.setf( std::ios::scientific, std::ios::floatfield );

			Vector< Vector<int> > ElementsInSrf = sit->second;

            Vector< Vector<int> >::iterator eit;

            for( eit = ElementsInSrf.begin(); eit != ElementsInSrf.end(); eit++ )
            {
                Vector<int> Element = *eit;
                
                for( int i=0; i<Element.size(); i++ ) 
                {
                    ElementsFile << std::noshowpos << Element[i] << "    "; 
                }
               
                ElementsFile << std::endl; 
            }

            ElementsFile.close();
        } 
    }

    //****************************************************************************************************************************************************
    //* - Writes nodes IDs and coordinates of a surface on a file
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_Nodes_OnFile( Vector<int>& NodesIDs, Vector< Vector<double> >& NodesCoords, int SrfaceID )
    {
        String NodesFileName;
        
        NodesFileName << "Fields/Surface-" << SrfaceID << "-Nodes.dat";

        std::fstream NodesFile( NodesFileName, std::ios::out );

        NodesFile.setf( std::ios::scientific, std::ios::floatfield );
        
        Vector< Vector<double> >::iterator cit = NodesCoords.begin();

        Vector<int>::iterator nit;

		for( nit = NodesIDs.begin(); nit != NodesIDs.end(); nit++, cit++ )
		{
            NodesFile << std::noshowpos << (*nit)    << "    " 
                      << std::showpos   << (*cit)[0] << "    " 
                      << std::showpos   << (*cit)[1] << "    " 
                      << std::showpos   << (*cit)[2] << std::endl;  
        }

        NodesFile.close();
    }

    //****************************************************************************************************************************************************
    //* - Writes volume fields on a file
    //****************************************************************************************************************************************************
	void Write_VolIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnVol, double Freq, int VolumeID )
    {
        String FieldsFileName; 
            
        FieldsFileName << "Fields/Volume-" << VolumeID << "-" << FieldName <<".dat";

        std::fstream FieldsFile( FieldsFileName, std::ios::out | std::ios::app );

        FieldsFile.setf( std::ios::scientific, std::ios::floatfield );

        // Write frequency on file
		double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        FieldsFile << "# " << vf << std::endl;

        // Write fields on file
        Vector< Vector< std::complex<double> > >::iterator fit;

		for( fit = FieldsOnVol.begin(); fit != FieldsOnVol.end(); fit++ )
		{
            FieldsFile << std::showpos << std::real( (*fit)[0] ) << "    " << std::showpos << std::imag( (*fit)[0] ) << "    " 
                       << std::showpos << std::real( (*fit)[1] ) << "    " << std::showpos << std::imag( (*fit)[1] ) << "    " 
                       << std::showpos << std::real( (*fit)[2] ) << "    " << std::showpos << std::imag( (*fit)[2] ) << std::endl;
        }  

        FieldsFile.close();
    }

    //****************************************************************************************************************************************************
    //* - Writes surface fields on a file
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnSrf, double Freq, int SrfaceID )
    {
        String FieldsFileName; 
            
        FieldsFileName << "Fields/Surface-" << SrfaceID << "-" << FieldName <<".dat";

        std::fstream FieldsFile( FieldsFileName, std::ios::out | std::ios::app );

        FieldsFile.setf( std::ios::scientific, std::ios::floatfield );

        // Write frequency on file
		double pi = 3.141592653589793;
        double vf = Freq / ( 2 * pi );

        FieldsFile << "# " << vf << std::endl;

        // Write fields on file
        Vector< Vector< std::complex<double> > >::iterator fit;

		for( fit = FieldsOnSrf.begin(); fit != FieldsOnSrf.end(); fit++ )
		{
            FieldsFile << std::showpos << std::real( (*fit)[0] ) << "    " << std::showpos << std::imag( (*fit)[0] ) << "    " 
                       << std::showpos << std::real( (*fit)[1] ) << "    " << std::showpos << std::imag( (*fit)[1] ) << "    " 
                       << std::showpos << std::real( (*fit)[2] ) << "    " << std::showpos << std::imag( (*fit)[2] ) << std::endl;
        }  

        FieldsFile.close();    
    }

    //****************************************************************************************************************************************************
    //* - Writes volume fields on a file ( complex frequency mode )
    //****************************************************************************************************************************************************
	void Write_VolIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnVol, std::complex<double> cFreq, int VolumeID )
    {
        String FieldsFileName; 
            
        FieldsFileName << "Fields/Volume-" << VolumeID << "-" << FieldName <<".dat";

        std::fstream FieldsFile( FieldsFileName, std::ios::out | std::ios::app );

        FieldsFile.setf( std::ios::scientific, std::ios::floatfield );

        // Write frequency on file
		FieldsFile << "# " << std::real( cFreq ) << "   " << std::imag( cFreq ) << std::endl;

        // Write fields on file
        Vector< Vector< std::complex<double> > >::iterator fit;

		for( fit = FieldsOnVol.begin(); fit != FieldsOnVol.end(); fit++ )
		{
            FieldsFile << std::showpos << std::real( (*fit)[0] ) << "    " << std::showpos << std::imag( (*fit)[0] ) << "    " 
                       << std::showpos << std::real( (*fit)[1] ) << "    " << std::showpos << std::imag( (*fit)[1] ) << "    " 
                       << std::showpos << std::real( (*fit)[2] ) << "    " << std::showpos << std::imag( (*fit)[2] ) << std::endl;
        }  

        FieldsFile.close();
    }

    //****************************************************************************************************************************************************
    //* - Writes surface fields on a file ( complex frequency mode )
    //****************************************************************************************************************************************************
    void Write_SrfIntgs_Fields_OnFile( String& FieldName, Vector< Vector< std::complex<double> > >& FieldsOnSrf, std::complex<double> cFreq, int SrfaceID )
    {
        String FieldsFileName; 
            
        FieldsFileName << "Fields/Surface-" << SrfaceID << "-" << FieldName <<".dat";

        std::fstream FieldsFile( FieldsFileName, std::ios::out | std::ios::app );

        FieldsFile.setf( std::ios::scientific, std::ios::floatfield );

        // Write frequency on file
		FieldsFile << "# " << std::real( cFreq ) << "   " << std::imag( cFreq ) << std::endl;

        // Write fields on file
        Vector< Vector< std::complex<double> > >::iterator fit;

		for( fit = FieldsOnSrf.begin(); fit != FieldsOnSrf.end(); fit++ )
		{
            FieldsFile << std::showpos << std::real( (*fit)[0] ) << "    " << std::showpos << std::imag( (*fit)[0] ) << "    " 
                       << std::showpos << std::real( (*fit)[1] ) << "    " << std::showpos << std::imag( (*fit)[1] ) << "    " 
                       << std::showpos << std::real( (*fit)[2] ) << "    " << std::showpos << std::imag( (*fit)[2] ) << std::endl;
        }  

        FieldsFile.close();    
    }
}