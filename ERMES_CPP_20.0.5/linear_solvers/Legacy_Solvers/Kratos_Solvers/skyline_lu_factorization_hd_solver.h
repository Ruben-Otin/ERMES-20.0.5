//*******************************************************************************  
//*   Project Name:        Kratos       
//*   Last Modified by:    $Author: pooyan $
//*   Date:                $Date: 2004/03/26 16:38:45 $
//*   Revision:            $Revision: 1.6 $
//*
//********************************************************************************


#if !defined(KRATOS_SKYLINE_LU_FACTORIZATION_HD_SOLVER_H_INCLUDED )
#define  KRATOS_SKYLINE_LU_FACTORIZATION_HD_SOLVER_H_INCLUDED

// System includes 
#include <iostream>
#include <fstream>
#include <stdio.h>


// External includes 
#include "../../external_libraries/boost/boost_1_35_0/boost/smart_ptr.hpp"

// Project includes
#include "../includes/define.h"
#include "direct_solver.h"
#include "../vectorial_spaces/kratos_vectorial_spaces.h"



namespace Kratos
{

class LUSkylineDumpUHD 
{
	public:

		int dumpUHD;        // dump U to hard drive when U larger than dumpUHD Mb
		const char* fileName; // name of the file used to store data
		int size;           // size of matrix
		double** U;         // columns of U
		int* hU   ;         // height of columns of U (including diagonal entry)
		double* LU;         // current row of L/U
		int smnzLU;         // smallest index to nonzero entry in L/U
		int grnzLU;         // greatest index to nonzero entry in L/U
		int* endRowSkyline; // end of skyline format rows
		int* minColRequired;// used to decide which columns of U can be saved
		int lastColumnSaved;// last column saved
        std::fstream* saveFile;  // stream used for saving data
		bool fileOpened;    // whether the file to contain data has been opened

		// WARNING: the following members should be integer, but they are defined 
		// as double to make sure they can have values over 2^32. The type of 
		// U and LU might be changed to float, complex, etc, but these members
		// could only be changed to long integer.
		double memoryU;     // memory occupied by U. It is double to make sure we
		                    // can deal with huge numbers.
		double limitU;      // dumpUHD in bytes
		double datasize;    // size in bytes of an entry of U(8 for double, etc)
		double fileSize;    // size of columns of U already dumped to file
		// the following 5 members are provided only for the user's information:
		// a little bit of time could be saved by commenting out the lines marked
		// with //#
		double totNZU;      // size of skyline U                              //#
		double totNZL;      // size of skyline L                              //#
		double areaSkyline; // provided only for user information             //#
		double nURNZ;       // elements in U that were really different to 0. //#
		double nLRNZ;       // elements in L that were really different to 0. //#

        LUSkylineDumpUHD(int dumpSize = 50, const char* name="KratosU.dat")
        {
            dumpUHD    = dumpSize;
            fileName   = name;
            fileOpened = false;
        } 

        //******************************************************************************************
		//* by default, dump when 50 megabytes to binary file KratosU.dat
        //********************************************************************************************
		void Solve(CSRMatrix<double> A, Vector<double>& rightHandSide, Vector<double>& solution)
		{
		    int i,j,k,m,index,indexj,p_sol;
		    double Lik;

		    // A is suposed to be reordered already!
		    size = A.RowsNumber();

		    if(size == 0)
		      return;

		    // initializations
		    U   = new double* [size];
		    hU  = new int     [size];
		    LU  = new double  [size];
		    for (i=0; i<size; i++) 
            {
		      U[i]=NULL;
		      hU[i]=0;
		      LU[i]=0.00;
		      solution[i]=rightHandSide[i];
		    }
		    memoryU= 0.00;
		    limitU= (double)(dumpUHD * 1048576L);

		    // somewhere else, it should have been checked that A[0] exists
		    datasize= (double)sizeof(double); 
		    totNZU= 0.00;      //#
		    totNZL= 0.00;      //#
		    nURNZ = 0.00;      //#
		    nLRNZ = 0.00;      //#
		    areaSkyline=0.00;  //#
		    lastColumnSaved=-1;
		    fileOpened=false;
		    fileSize=0;

		    // We will need some knowledge about the skyline format resulting from 
		    // this matrix. endRowSkyline[i] is the largest j of a nonzero in row i of
		    // the skyline pattern. minColRequired[i] is the minimum column that will 
		    // be required along the skyline computation starting from (i,i); it is 
		    // used to decide which columns can be saved.
		    // Initialization:
		    endRowSkyline = new int[size];
		    minColRequired= new int[size];
		    for (i=0; i<size; i++) 
            {
		      endRowSkyline[i] =-1;
		      minColRequired[i]=i;
		    }
		    // traverse matrix looking at nonzeros
		    for (i=0; i<size; i++) 
            {
                for (CSRMatrix<double>::RowIteratorType row_iterator =A[i].begin(); row_iterator != A[i].end(); ++row_iterator)
                {
			        j= row_iterator->first;
			        if (j>endRowSkyline [i]) endRowSkyline [i]=j;
			        if (j<minColRequired[i]) minColRequired[i]=j;
		        }
		    }
		    // actual data
		    for (i=1; i<size; i++) 
            {
		      if (endRowSkyline[i]<endRowSkyline[i-1]) 
              {
			     endRowSkyline[i] = endRowSkyline[i-1];
		      }
		    }
            for (i=size-2; i>=0; i--) 
            {
              if ( minColRequired[i] > minColRequired[i+1] )
              {
                  minColRequired[i] = minColRequired[i+1];
              }
            }

		    
		    // MAIN LOOP
		    for (i=0; i<size; i++) 
            { 
		          // copy row i of A into LU
		          // The "sparse" row LU is here implemented in yet another different way,
		          // since it will actually be quite dense.
		          // LU[i] is the actual entry; if the entry is a zero, then LU[i]=0.
		          // The only information on its structure that is kept is that all of its
		          // nonzeros are included in the interval [smnzLU,grnzLU].
		          smnzLU=size;
		          grnzLU=-1;
         
		          for (CSRMatrix<double>::RowIteratorType row_iterator = A[i].begin() ; row_iterator != A[i].end() ; ++row_iterator) 
                  {
		      	      j=row_iterator->first;
		      	      LU[j] = row_iterator->second;
		      	      if (j<smnzLU) smnzLU=j;
		      	      if (j>grnzLU) grnzLU=j;
		          }

		          if (smnzLU==size) 
                  {
			        error();
			        KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::solve", "matrix had an empty row", "");
		          }
		          
		          // Actually, we should add grnzLU-smnzLU+1, but it will be faster to add
		          // size outside of the i-loop.
		          areaSkyline += (double)(grnzLU-smnzLU); //#
		          totNZL += (double)(i-smnzLU); //#
		          
		          // zero out entries in L
		          for (k=smnzLU; k<i; k++) 
                  {
			        if (LU[k] != 0.00) 
                    {
			          nLRNZ += 1.00; //#

			          if (U[k][hU[k]-1]==0.00) 
                      {
			            error();
			            KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::solve","zero on diagonal, factorization", ""); 
			          }
			          Lik= LU[k]/U[k][hU[k]-1];  // diagonal U(k,k)=U[k][hU[k]-1]

			          solution[i]=solution[i]-Lik*solution[k];
			          LU[k] = 0.00; // clean up LU
			          
			          // the normal loop would be :
			          //      for j=k+1 to size :    LU(j)=LU(j)-Lik*U(k,j)
			          // however, we know where the nonzeros end in the skyline 
			          // format, so use   for j=k+1 to endRowSkyline[k]
			          for (j=k+1; j<=endRowSkyline[k]; j++) 
                      {
			            // U(k,j) would be U[j][hU[j]-1+k-j] if it was stored...
			            index=hU[j]-1+k-j;
			            if (index>=0) 
                        {
			              LU[j]=LU[j]-Lik*U[j][index];
			              if (j>grnzLU) grnzLU=j;
			            }
			          }
                    } // end if (LU[k] != 0.00)
		          } // end for k
		          
		          // what is left in LU is part of U; store it.
		          for (j=i; j<=grnzLU; j++) 
                  {
			        if (LU[j] != 0.00) 
                    {
			          if (U[j] != NULL) 
                      {
			            U[j][hU[j]-1+i-j]=LU[j];
			          }
                      else
                      {
			            // we have to allocate a new column of U starting here
			            hU[j]=j-i+1;
			            U[j]=new double[hU[j]];
			            U[j][0]=LU[j];
			            for (m=1; m<hU[j]; m++) U[j][m]=0.00;
			            totNZU += (double)hU[j]; //#
			            memoryU += datasize*(double)hU[j];
			            // using too much memory already? Before calling save, check also
                        // that something can truly be saved
			            if ( (memoryU > limitU) && (minColRequired[i] >= lastColumnSaved+2)) save(i);
			          }
			          LU[j]=0.00;
                    }
                  }
		    } // end of for i
		    
		    // count number of nonzeros in the non-saved part of U
		    // (when a column is saved, its nonzeros are counted)
		    for (i=lastColumnSaved+1; i<size; i++) 
            {    //#
		      for (j=0; j<hU[i]; j++) 
              {               //#
			      if (U[i][j] != 0.00) nURNZ += 1.00; //#
		      }                                       //#
		    }                                           //#
		    areaSkyline += (double)size;                //#
		    
		    // back solve;
		    for (i=size-1; i>=0; i--) 
            {
		      if (U[i] == NULL) load(i);
		      if (U[i][hU[i]-1]==0.00) 
              {
			    error();
			    KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::solve","zero on diagonal, backsolve", ""); 
		      }
		      solution[i]=solution[i]/U[i][hU[i]-1];
		      p_sol=i-hU[i]+1;
		      for (j=0; j<hU[i]-1; j++) 
              {
			    solution[p_sol]=solution[p_sol]-U[i][j]*solution[i];
			    p_sol++;
		      }
		      delete[] U[i]; // free up memory
		      U[i]=NULL;
		      memoryU -= datasize*(double)hU[i];
		    }
		    
		    // delete everything
		    for (i=0; i<size; i++) if (U[i] != NULL) delete[] U[i];
		    delete[] U;
		    delete[] hU;
		    delete[] LU;
		    delete[] endRowSkyline;
		    delete[] minColRequired;		    

		    //if this block is commented out, comment out also the lines with //#   
		    std::cout << "Linear solver DirectSolverHD memory report:"<<std::endl; 
            std::cout << std::endl;
		    std::cout << "Maximum RAM: " << limitU/1048576   <<" MB."<<std::endl;     
		    std::cout << "File size  : " << fileSize/1048576 <<" MB."<<std::endl; 
            std::cout << std::endl;
         
        }



        //************************************************************************************************************
        //* save to hard drive all the columns of U that will not be required
		//* in the factorization of the rows below workingRow
		//* The file is not opened until we are sure it is necessary.
        //************************************************************************************************************
		void save(int workingRow)
        {
              std::cout << "Saving in file " << fileName << "..." << std::endl; 

		      if (!fileOpened) 
              {
		            fileOpened=true;
		            // std::cout << "    Opening file " << fileName << "\n";
		            // The file will be binary, used for reading and writing, new data 
		            // will be added only at the end, and if it exists when it is opened, 
		            // then its contents will be erased. Beware; no two instances of this 
		            // program can run simultaneously unless some precaution is taken for
		            // their files to have different names.
                    saveFile = new std::fstream(fileName , 
                                                std::ios::binary | 
                                                std::ios::in     | 
                                                std::ios::out    | 
                                                std::ios::app );

                    if (!(saveFile->is_open()))
                    {
		                 error();
		                 KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::save",
                                      "Error when opening file.", ""); 
		            }

		      }

		      for (int i=lastColumnSaved+1; i<=minColRequired[workingRow]-1; i++) 
              {
		            // count the number of entries of U that are really nonzero, so that the
		            // user can find out whether this method is good enough.
		            for (int j=0; j<hU[i]; j++) 
                    {                                               //#
		                 if (U[i][j] != 0.00) nURNZ = nURNZ + 1.00; //#
                    }                                               //#

		            saveFile->write( (char*) U[i] , (int)(datasize)*hU[i] );
		            memoryU -= datasize*(double)hU[i];
		            fileSize += datasize*(double)hU[i];

		            if (saveFile->fail()) 
                    {
		                 error();
		                 KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::save",
                                       "Error when writing to file.(out of space in disk?) File erased", ""); 
		            }

		            delete[] U[i];
		            U[i]=NULL; // to remember it has been deallocated
		      }

		      lastColumnSaved= minColRequired[workingRow]-1;

              std::cout << "U size   : " << memoryU/1048576  <<" MB."<<std::endl; 
		      std::cout << "File size: " << fileSize/1048576 <<" MB."<<std::endl;

              std::cout << "Saving finished." << std::endl;
              std::cout << std::endl;        
        }



        //************************************************************************************************************
        //*  This method is called when a column of U is required for the backsolve 
		//*  but it has been saved to the hard drive. We will always read that column,
		//*  and maybe more untill reaching a total of limitU bytes. 
        //************************************************************************************************************  
		void load(int columnRequired)
        {
              std::cout << "Loading from file " << fileName << "..." << std::endl; 

		      int firstToRead;      // the first column that we will read
		      double memToRead;     // total amount of memory to read
		      
		      firstToRead = columnRequired;
		      memToRead   = datasize*(double)hU[columnRequired];

		      while ((firstToRead>0) && ( memToRead+datasize*(double)hU[firstToRead-1] < limitU)) 
              {
		           firstToRead--;
		           memToRead += datasize*(double)hU[firstToRead];
		      }

		      // rewind file;
              saveFile->seekg( -(int)memToRead , std::ios::cur);

		      for (int i=firstToRead; i<=columnRequired; i++) 
              {
		           U[i]= new double[hU[i]];
		           saveFile->read( (char*) U[i] , (int)(datasize)*hU[i] );
		           if (saveFile->fail()) 
                   {
		              error();
		              KRATOS_ERROR(std::runtime_error,"LUSkylineDumpUHD::load","Error! File erased", ""); 
		           }
		      }

		      // rewind file again
              saveFile->seekg( -(int)memToRead , std::ios::cur);

              std::cout << "Loading finished."<< std::endl; 
              std::cout << std::endl;  
		}


        //************************************************************************************************************
        //* an error has happened; before throwing an exception or whatever,
		//* check whether the file had been opened, and, if so, close the stream, 
		//* delete it, and erase the file.
        //************************************************************************************************************  
		void error()
		{
		    if (fileOpened) 
            {
		         saveFile->close();
		         delete saveFile;
		         _unlink(fileName);	
		    }
		}


        //************************************************************************************************************
		//* Destructor, check whether the file had been opened, and, if so, close the stream, 
		//* delete it, and erase the file.
        //************************************************************************************************************  
        ~LUSkylineDumpUHD()
        {
            if (fileOpened) 
            {
		         saveFile->close();
		         delete saveFile;
		         _unlink(fileName);	
		    }
        }


};
 



 
template<class TSparseSpaceType,class TDenseSpaceType,class TReordererType>
class SkylineLUFactorizationHDSolver : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
    public:
     
        typedef boost::shared_ptr<SkylineLUFactorizationHDSolver> Pointer;

        typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 

        typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

        typedef typename TSparseSpaceType::VectorType VectorType;

        typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
      
        /// Default constructor.
        SkylineLUFactorizationHDSolver(){}


        //Constructor
        SkylineLUFactorizationHDSolver(String fileName, int MaxRAM)
        {
            //para generación automática de nombres de archivo
            //usar la función: char *tmpnam( char *string );
              mFileName= fileName + ".dat";
              mMaxRAM  = MaxRAM;
        }

        /// Copy constructor.
        SkylineLUFactorizationHDSolver(const SkylineLUFactorizationHDSolver& Other) : BaseType(Other){}

        /// Destructor.
        virtual ~SkylineLUFactorizationHDSolver(){}

        /// Assignment operator.
        SkylineLUFactorizationHDSolver& operator=(const SkylineLUFactorizationHDSolver& Other)
        {
            BaseType::operator=(Other);
            return *this;
        }
        
        /** Normal solve method.
        Solves the linear system Ax=b and puts the result on SystemVector& rX. 
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector.
        @param rB. Right hand side vector.
        */
        bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
        {
            if(IsNotConsistent(rA, rX, rB)) return false;

            typedef CuthillMcKeeReorderer<CSRSpace<double>, DenseSpace<double> > CMReorderType;
            CMReorderType::Pointer CMReorder;
            CMReorder = CMReorderType::Pointer(new CuthillMcKeeReorderer<CSRSpace<double>, DenseSpace<double> >());

            LUSkylineDumpUHD myFactorization(mMaxRAM,mFileName);

            std::cout << "Applying Cuthill-McKee reorder..." << std::endl;

            CMReorder->Initialize(rA, rX, rB);
            CMReorder->Reorder(rA, rX, rB);

            std::cout << "Reorder finish." << std::endl; 
            std::cout << std::endl;  

            myFactorization.Solve(rA, rB, rX);

            std::cout << "Applying inverse reorder..." << std::endl;

            CMReorder->InverseReorder(rA, rX, rB);

            std::cout << "Inverse reorder finish." << std::endl;
            std::cout << std::endl;

            return true;
        }
      

        /** Multi solve method for solving a set of linear systems with same coefficient matrix.
        Solves the linear system Ax=b and puts the result on SystemVector& rX. 
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector.
        @param rB. Right hand side vector.*/
        bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
        {   
            return false;
        }

     
      
        /// Print information about this object.
        void  PrintInfo(std::ostream& rOStream) const
        {
            rOStream << "LU_HD factorization solver finished.";
        }

        /// Print object's data.
        void  PrintData(std::ostream& rOStream) const 
        {
        }

    private:
            
            int mMaxRAM;

            String mFileName;
        
}; // Class SkylineLUFactorizationHDSolver 

 
    /// input stream function
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, 
				                      SkylineLUFactorizationHDSolver<TSparseSpaceType, 
                                                                     TDenseSpaceType, 
				                                                     TReordererType>& rThis)
    {
    }

    /// output stream function
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
				                      const SkylineLUFactorizationHDSolver<TSparseSpaceType, 
                                                                           TDenseSpaceType, 
				                                                           TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }

}  // namespace Kratos.

#endif // KRATOS_SKYLINE_LU_FACTORIZATION_HD_SOLVER_H_INCLUDED  defined 


