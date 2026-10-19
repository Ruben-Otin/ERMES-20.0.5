#define YY_InputCParser_h_included

/*  A Bison++ parser, made from input_c_parser.y  */

 /* with Bison++ version bison++ Version 1.21-8, adapted from GNU bison by coetmeur@icdc.fr
  */

/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Bob Corbett and Richard Stallman

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* HEADER SECTION */
#if defined( _MSDOS ) || defined(MSDOS) || defined(__MSDOS__) 
#define __MSDOS_AND_ALIKE
#endif
#if defined(_WINDOWS) && defined(_MSC_VER)
#define __HAVE_NO_ALLOCA
#define __MSDOS_AND_ALIKE
#endif

#ifndef alloca
#if defined( __GNUC__)
#define alloca __builtin_alloca

#elif (!defined (__STDC__) && defined (sparc)) || defined (__sparc__) || defined (__sparc)  || defined (__sgi)
#include <alloca.h>

#elif defined (__MSDOS_AND_ALIKE)
#include <malloc.h>
#ifndef __TURBOC__
/* MS C runtime lib */
#define alloca _alloca
#endif

#elif defined(_AIX)
#include <malloc.h>
#pragma alloca

#elif defined(__hpux)
#ifdef __cplusplus
extern "C" {
void *alloca (unsigned int);
};
#else /* not __cplusplus */
void *alloca ();
#endif /* not __cplusplus */

#endif /* not _AIX  not MSDOS, or __TURBOC__ or _AIX, not sparc.  */
#endif /* alloca not defined.  */
#ifdef c_plusplus
#ifndef __cplusplus
#define __cplusplus
#endif
#endif
#ifdef __cplusplus
#ifndef YY_USE_CLASS
#define YY_USE_CLASS
#endif
#else
#ifndef __STDC__
#define const
#endif
#endif
#include <stdio.h>
#define YYBISON 1  

#include <functional>   // std::plus ( BinaryVariableStatement )
#include "../includes/variables.h"
  class Statement;
  class BlockStatement;
  template<class TDataType> class ValueStatement;

typedef union {
  int integer_value;
  double double_value;
  char string_value[255];
  std::istream* istream_value;
  Statement* statement_handler;
  BlockStatement* block_statement_handler;
  ValueStatement<double>* double_value_statement_handler;
  Kratos::Vector<int>* vector_integer_value;
  Kratos::Vector<Kratos::Vector<int> >* vector_vector_integer_value;
  Kratos::Vector<double>* vector_double_value;
  Kratos::Vector<Kratos::Vector<double> >* vector_vector_double_value;
  Kratos::Matrix<double>* matrix_double_value;
  const Kratos::Variable<double>* double_variable;
  const Kratos::Variable<Kratos::Vector<double> >* vector_double_variable;
  const Kratos::Variable<Kratos::Matrix<double> >* matrix_variable;
} yy_InputCParser_stype;
#define YY_InputCParser_STYPE yy_InputCParser_stype

#include "../includes/kernel.h"
#include "../includes/input_c_scanner.h"
#include "../includes/element_sources.h"

using Kratos::Exception;

  Kratos::String block_name;

  class Statement
  {
    public:
	
      Statement(){}

      // Virtual destructor: statements are deleted through Statement* (shared_ptr<Statement>, delete $n)
      virtual ~Statement(){}
	  
      virtual void Execute(Kratos::Kernel* pKernel){}
	  
      virtual int Value(Kratos::Kernel* pKernel){return 0;}
  };

  class BlockStatement : public Statement
  {
    public:

      BlockStatement(){}

      void AddStatement(Statement* pNewStatement)
      {
          mStatements.push_back(std::shared_ptr<Statement>(pNewStatement));
      }

      void Execute(Kratos::Kernel* pKernel)
      {
          for(Kratos::Vector<std::shared_ptr<Statement> >::iterator i = mStatements.begin() ; i != mStatements.end() ; i++) (*i)->Execute(pKernel);
      }

      Kratos::Vector<std::shared_ptr<Statement> > mStatements;
  };

  //****************************************************************************************
  //* - Generate nodes
  //****************************************************************************************
  int node_Id;

  double node_X, node_Y, node_Z;

  class GenerateNodeStatement : public Statement
  {
    public:
  
      GenerateNodeStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateNode(node_Id, node_X, node_Y, node_Z);
      }
  };

  //****************************************************************************************
  //* - Generate properties
  //****************************************************************************************
  template<class TDataType> class GeneratePropertiesStatement : public Statement
  {
    public:

      GeneratePropertiesStatement( int Id, 
                                   const Kratos::Variable<TDataType>& rVariable, 
                                   const TDataType& Value ) : mId(Id), 
                                                              mVariable(rVariable), 
                                                              mValue(Value){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->SetProperties(mId, mVariable, mValue);
      }

      int mId;
      Kratos::Variable<TDataType> mVariable;
      TDataType mValue;
  };

  //****************************************************************************************
  //* - Set elements group
  //****************************************************************************************
  class SetElementsGroup : public Statement
  {
    public:

      SetElementsGroup(const Kratos::String& ElementsGroupName) : mElementsGroupName(ElementsGroupName){}

      void Execute(Kratos::Kernel* pKernel)
      {
          block_name = mElementsGroupName;
		  
          pKernel->Ini_Building();
      }

      Kratos::String mElementsGroupName;
  };

  //****************************************************************************************
  //* - Set problem type
  //****************************************************************************************
  class SetProblemType : public Statement
  {
    public:

      SetProblemType(const Kratos::String& ProblemType):mProblemType(ProblemType){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->SetProblemType(mProblemType);
      }

      Kratos::String mProblemType;
  };

  //****************************************************************************************
  //* - Set frequency
  //****************************************************************************************
  class SetFrequency : public Statement
  {
    public:

      SetFrequency(double Freq):mFrequency(Freq){}

      void Execute(Kratos::Kernel* pKernel)
      {
         pKernel->SetFrequency(mFrequency);
      }

      double mFrequency;
  };

  //****************************************************************************************
  //* - Set sweep frequency
  //****************************************************************************************
  class SweepFrequency : public Statement
  {
    public:

      SweepFrequency(double IniFreq,
                     double EndFreq,
                     double StepFreq) : mIniFrequency(IniFreq),
                                        mEndFrequency(EndFreq),
                                        mStepFrequency(StepFreq){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->SweepFrequency(mIniFrequency,mEndFrequency,mStepFrequency);
      }

      double mIniFrequency;
      double mEndFrequency;
      double mStepFrequency;
  };

  //****************************************************************************************
  //* - Dynamic vectors
  //****************************************************************************************
  Kratos::Vector<int> sNodesId;
  
  int sPropertiesId;

  //****************************************************************************************
  //* - Generate PEC_Elements
  //****************************************************************************************
  class GeneratePECElementStatement: public Statement
  {
    public:

      GeneratePECElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GeneratePECElement(sNodesId);
          
		  sNodesId.FreeData();
      }
  };
  
  //****************************************************************************************
  //* - Generate PMC_Elements
  //****************************************************************************************
  class GeneratePMCElementStatement: public Statement
  {
    public:

      GeneratePMCElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GeneratePMCElement(sNodesId);
          
		  sNodesId.FreeData();
      }
  };

  //****************************************************************************************
  //* - Generate TEC_Elements
  //****************************************************************************************
  class GenerateTECElementStatement: public Statement
  {
    public:

      GenerateTECElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateTECElement(sNodesId);
          
		  sNodesId.FreeData();
      }
  };

  //****************************************************************************************
  //* - Contact generation
  //****************************************************************************************
  class GenerateContactStatement: public Statement
  {
    public:

      GenerateContactStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateContactPairs(sNodesId);
          
		  sNodesId.FreeData();
      }
  };

  //****************************************************************************************
  //****************************************************************************************
  //** - FAST element generation
  //****************************************************************************************
  //****************************************************************************************

  // Static arrays for nodes
  int stNodesId[4];

  //****************************************************************************************
  //* - Generate volume elements
  //****************************************************************************************
  class GenerateVolumeElementStatement: public Statement
  {
    public:

      GenerateVolumeElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateVolumeElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate source elements
  //****************************************************************************************
  class GenerateSourceElementStatement: public Statement
  {
    public:

      GenerateSourceElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateSourceElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate volume integration elements
  //****************************************************************************************
  class GenerateVolIntgElementStatement: public Statement
  {
    public:

      GenerateVolIntgElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateVolIntgElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate PBC elements
  //****************************************************************************************
  class GeneratePBCElementStatement: public Statement
  {
    public:

      GeneratePBCElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GeneratePBCElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate surface integration elements
  //****************************************************************************************
  class GenerateSrfIntgElementStatement: public Statement
  {
    public:

      GenerateSrfIntgElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateSrfIntgElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate generic Robin elements
  //****************************************************************************************
  class GenerateGenRobinElementStatement: public Statement
  {
    public:

      GenerateGenRobinElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateGenericRobinElement(stNodesId, sPropertiesId);
      }
  };
  
  //****************************************************************************************
  //* - Generate plasma far field elements
  //****************************************************************************************
  class GeneratePlasmaFarFieldElementStatement: public Statement
  {
    public:

      GeneratePlasmaFarFieldElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GeneratePlasmaFarFieldElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate rectangular waveguide port TE10 elements
  //****************************************************************************************
  class GenerateRWPortTE10ElementStatement: public Statement
  {
    public:

      GenerateRWPortTE10ElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateRWPortTE10Element(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate coaxial port TEM elements
  //****************************************************************************************
  class GenerateCoaxPortTEMElementStatement: public Statement
  {
    public:

      GenerateCoaxPortTEMElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateCoaxPortTEMElement(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate far field elements
  //****************************************************************************************
  class GenerateFarFieldElementStatement: public Statement
  {
    public:

      GenerateFarFieldElementStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateFarFieldElement(stNodesId);
      }
  };

  //****************************************************************************************
  //* - Generate projection over rectangular waveguide mode TE10
  //****************************************************************************************
  class GenerateProjectionRWPortTE10Statement: public Statement
  {
    public:

      GenerateProjectionRWPortTE10Statement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateProjectionRWPortTE10(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //* - Generate projection over coaxial port TEM mode
  //****************************************************************************************
  class GenerateProjectionCoaxialTEMStatement: public Statement
  {
    public:

      GenerateProjectionCoaxialTEMStatement(){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateProjectionCoaxialTEM(stNodesId, sPropertiesId);
      }
  };

  //****************************************************************************************
  //****************************************************************************************
  //** - Fixing degrees of freedom
  //****************************************************************************************
  //****************************************************************************************

  //****************************************************************************************
  //* - Fix degrees of freedom
  //****************************************************************************************
  class FixDofStatement : public Statement
  {
    public:

      FixDofStatement(int NodeId, 
                      const Kratos::Variable<double>& rVariable, 
                      const double& Value) : mId(NodeId), 
                                             mVariable(rVariable), 
                                             mValue(Value){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->FixDof(mId, mVariable, mValue);
      }

      int mId;
      Kratos::Variable<double> mVariable;
      double mValue;
  };

  //****************************************************************************************
  //* - Fix vector degrees of freedom
  //****************************************************************************************
  class FixCDofStatement : public Statement
  {
    public:

      FixCDofStatement(int NodeId, 
                       const Kratos::Variable<double>& rVariable, 
                       const Kratos::Vector<double>& Value) : mId(NodeId), 
                                                              mVariable(rVariable), 
                                                              mValue(Value){}

      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->FixCDof(mId, mVariable, mValue);
      }

      int mId;
      Kratos::Variable<double> mVariable;
      Kratos::Vector<double> mValue;
  };

  //****************************************************************************************
  //* - Generate linear solver
  //****************************************************************************************
  class GenerateLinearSolverStatement : public Statement
  {
    public:

      GenerateLinearSolverStatement(const Kratos::String& Name, 
                                    const Kratos::String& SolverType, 
                                    int MaxNumberOfIterations, 
                                    int StepIterations,
                                    double Tolerance) : mName(Name), 
                                                        mSolverType(SolverType), 
                                                        mMaxNumberOfIterations(MaxNumberOfIterations), 
                                                        mStepIterations(StepIterations), 
                                                        mTolerance(Tolerance){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->GenerateLinearSolver(mName, mSolverType, mMaxNumberOfIterations, mStepIterations, mTolerance);
		  
          mName.erase(mName.begin(), mName.end());
          
		  mSolverType.erase(mSolverType.begin(), mSolverType.end());
      }
	  
      Kratos::String mName;
      Kratos::String mSolverType;
	  double mTolerance;
      int mMaxNumberOfIterations;
      int mStepIterations;
  };

  //****************************************************************************************
  //* - Solve problem statement
  //****************************************************************************************
  class SolveStatement : public Statement
  {
    public:

      SolveStatement(const Kratos::String& StrategyName, 
                     const Kratos::String& ElementsGroup) : mStrategyName(StrategyName), 
                                                            mElementsGroup(ElementsGroup){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->Solve(mStrategyName, mElementsGroup);
	  
          mStrategyName.erase(mStrategyName.begin(), mStrategyName.end());
          
		  mElementsGroup.erase(mElementsGroup.begin(), mElementsGroup.end()); 
      }
	  
      Kratos::String mStrategyName;
      Kratos::String mElementsGroup;
  };

  //****************************************************************************************
  //* - For loop statement
  //****************************************************************************************
  class ForStatement : public Statement
  {
    public:

      ForStatement(Statement* pFirst, Statement* pSecond, Statement* pThird, Statement* pForth) :
                   mpFirst(pFirst), mpSecond(pSecond), mpThird(pThird), mpForth(pForth) {}
	  
      ~ForStatement(){delete mpFirst; delete mpSecond; delete mpThird; delete mpForth;}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          for(mpFirst->Execute(pKernel); mpSecond->Value(pKernel); mpThird->Execute(pKernel)) 
		  {
		      mpForth->Execute(pKernel);
		  }
      }
	  
      Statement* mpFirst;
      Statement* mpSecond;
      Statement* mpThird;
      Statement* mpForth;
  };

  //****************************************************************************************
  //* - Print on nodes or Gauss points
  //****************************************************************************************
  template<class TDataType> class PrintStatement : public Statement
  {
    public:

      PrintStatement(const Kratos::Variable<TDataType>& rVariable) : mVariable(rVariable){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->Print(mVariable);
      }
	  
      Kratos::Variable<TDataType> mVariable;
  };
  
  //****************************************************************************************
  //* - Print on nodes
  //****************************************************************************************
  template<class TDataType> class PrintOnNodesStatement : public Statement
  {
    public:

      PrintOnNodesStatement(const Kratos::Variable<TDataType>& rVariable) : mVariable(rVariable){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->PrintOnNodes(mVariable);
      }
	  
      Kratos::Variable<TDataType> mVariable;
  };

  //****************************************************************************************
  //* - Calculate variable on nodes or Gauss points
  //****************************************************************************************
  template<class TDataType> class CalculateStatement : public Statement
  {
    public:

      CalculateStatement(const Kratos::Variable<TDataType>& rVariable) : mVariable(rVariable){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->Calculate(mVariable);
      }
	  
      Kratos::Variable<TDataType> mVariable;
  };

  //****************************************************************************************
  //* - Calculate variable on nodes
  //****************************************************************************************
  template<class TDataType> class CalculateNodalStatement : public Statement
  {
    public:

      CalculateNodalStatement(const Kratos::Variable<TDataType>& rVariable) : mVariable(rVariable){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->CalculateNodal(mVariable);
      }
	  
      Kratos::Variable<TDataType> mVariable;
  };
  
  //****************************************************************************************
  //* - Create time step statement
  //****************************************************************************************
  class CreateTimeStepStatement : public Statement
  {
    public:

      CreateTimeStepStatement(){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->CreateTimeStep();
      }
  };

  //****************************************************************************************
  //* - Create solution step statement
  //****************************************************************************************
  class CreateSolutionStepStatement : public Statement
  {
    public:

      CreateSolutionStepStatement(){}
	  
      void Execute(Kratos::Kernel* pKernel)
      {
          pKernel->CreateSolutionStep();
      }
  };

  //****************************************************************************************
  //****************************************************************************************
  //** - Assignment statements
  //****************************************************************************************
  //****************************************************************************************

  //****************************************************************************************
  //* - Logical operators statements
  //****************************************************************************************
  template<class TFunction, class TDataType> class LogicalStatement : public Statement
  {
    public:

      LogicalStatement(const Kratos::Variable<TDataType>& rVariable, ValueStatement<TDataType>* pSecond) :
                       mVariable(rVariable), mpSecond(pSecond){}
	  
      ~LogicalStatement(){delete mpSecond;}
	  
      int Value(Kratos::Kernel* pKernel)
      {
          TFunction function; 
		  
          return function(pKernel->GetValue(mVariable), mpSecond->GetValue(pKernel));
      }
	  
      Kratos::Variable<TDataType> mVariable;
      ValueStatement<TDataType>* mpSecond;
  };

  //****************************************************************************************
  //* - Get value statement
  //****************************************************************************************
  template<class TDataType> class ValueStatement
  {
    public:

      ValueStatement() {}

      // Virtual destructor: derived statements are deleted through ValueStatement<T>*
      virtual ~ValueStatement() {}
	  
      virtual TDataType GetValue(Kratos::Kernel* pKernel){return TDataType();}
  };

  //****************************************************************************************
  //* - Get constante value
  //****************************************************************************************
  template<class TDataType> class ConstantValueStatement : public ValueStatement<TDataType> 
  {
    public:

      ConstantValueStatement(const TDataType& rValue) : mValue(rValue){}
	  
      TDataType GetValue(Kratos::Kernel* pKernel){return mValue;}
	  
      TDataType mValue;
  };

  //****************************************************************************************
  //* - Get variable
  //****************************************************************************************
  template<class TDataType> class VariableStatement : public ValueStatement<TDataType> 
  {
    public:

      VariableStatement(const Kratos::Variable<TDataType>& rVariable) : mVariable(rVariable){}
	  
      TDataType GetValue(Kratos::Kernel* pKernel)
	  {
	      return pKernel->GetValue(mVariable);
	  }
	  
      Kratos::Variable<TDataType> mVariable;
  };

  //****************************************************************************************
  //* - Get binary variable
  //****************************************************************************************
  template<class TFunction, class TDataType> class BinaryVariableStatement : public ValueStatement<TDataType> 
  {
    public:

      BinaryVariableStatement(ValueStatement<TDataType>* pFirst, 
	                          ValueStatement<TDataType>* pSecond) : 
                              mpFirst(pFirst), 
							  mpSecond(pSecond){}
	  
      ~BinaryVariableStatement(){delete mpFirst; delete mpSecond;}
	  
      TDataType GetValue(Kratos::Kernel* pKernel)
      {
          TFunction function;
		  
          return function(mpFirst->GetValue(pKernel), mpSecond->GetValue(pKernel));
      }
	  
      ValueStatement<TDataType>* mpFirst;
      ValueStatement<TDataType>* mpSecond;
  };

  //****************************************************************************************
  //* - Assign variable statement
  //****************************************************************************************
  template<class TDataType> class AssigningVariableStatement : public Statement 
  {
    public:

      AssigningVariableStatement(const Kratos::Variable<TDataType>& rVariable, 
	                             ValueStatement<TDataType>* pSecond) : 
                                 mVariable(rVariable), 
								 mpSecond(pSecond){}
	  
      ~AssigningVariableStatement(){delete mpSecond;}
	  
      void Execute(Kratos::Kernel* pKernel)
	  {
	      pKernel->SetValue(mVariable, mpSecond->GetValue(pKernel));
	  }
	  
      Kratos::Variable<TDataType> mVariable;
      ValueStatement<TDataType>* mpSecond;
  };

  //****************************************************************************************
  //* - Assign nodal variables statement
  //****************************************************************************************
  template<class TDataType> class AssigningNodalVariableStatement : public Statement 
  {
    public:

      AssigningNodalVariableStatement(int NewNodeId, 
	                                  const Kratos::Variable<TDataType>& rVariable, 
									  ValueStatement<TDataType>* pSecond) : 
                                      mNodeId(NewNodeId), 
									  mVariable(rVariable), 
									  mpSecond(pSecond){}
	  
      ~AssigningNodalVariableStatement(){delete mpSecond;}
	  
      void Execute(Kratos::Kernel* pKernel)
	  {
		  pKernel->SetValue(mNodeId, mVariable, mpSecond->GetValue(pKernel));
	  }
	  
      int mNodeId;
      Kratos::Variable<TDataType> mVariable;
      ValueStatement<TDataType>* mpSecond;
  };

#define YY_InputCParser_CONSTRUCTOR_PARAM  Kratos::Kernel* pKernel, std::istream* pNewInput, std::ostream* pNewOutput
#define YY_InputCParser_CONSTRUCTOR_INIT  : mInputCScanner(pNewInput, pNewOutput), mpKernel(pKernel), mpInput(pNewInput)
#define YY_InputCParser_YY_InputGid_CONSTRUCTOR_CODE  
#define YY_InputCParser_LEX_BODY  { return mInputCScanner.yylex(yylval); }
#define YY_InputCParser_MEMBERS   Kratos::InputCScanner mInputCScanner; Kratos::Kernel* mpKernel; std::istream* mpInput; Kratos::String mBlockName;
/* %{ and %header{ and %union, during decl */
#define YY_InputCParser_BISON 1
#ifndef YY_InputCParser_COMPATIBILITY
#ifndef YY_USE_CLASS
#define  YY_InputCParser_COMPATIBILITY 1
#else
#define  YY_InputCParser_COMPATIBILITY 0
#endif
#endif

#if YY_InputCParser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YYLTYPE
#ifndef YY_InputCParser_LTYPE
#define YY_InputCParser_LTYPE YYLTYPE
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_InputCParser_STYPE 
#define YY_InputCParser_STYPE YYSTYPE
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_InputCParser_DEBUG
#define  YY_InputCParser_DEBUG YYDEBUG
#endif
#endif
#ifdef YY_InputCParser_STYPE
#ifndef yystype
#define yystype YY_InputCParser_STYPE
#endif
#endif
/* use goto to be compatible */
#ifndef YY_InputCParser_USE_GOTO
#define YY_InputCParser_USE_GOTO 1
#endif
#endif

/* use no goto to be clean in C++ */
#ifndef YY_InputCParser_USE_GOTO
#define YY_InputCParser_USE_GOTO 0
#endif

#ifndef YY_InputCParser_PURE
/*  YY_InputCParser_PURE */
#endif

/* section apres lecture def, avant lecture grammaire S2 */
/* prefix */
#ifndef YY_InputCParser_DEBUG
/* YY_InputCParser_DEBUG */
#endif


#ifndef YY_InputCParser_LSP_NEEDED
 /* YY_InputCParser_LSP_NEEDED*/
#endif



/* DEFAULT LTYPE*/
#ifdef YY_InputCParser_LSP_NEEDED
#ifndef YY_InputCParser_LTYPE
typedef
  struct yyltype
    {
      int timestamp;
      int first_line;
      int first_column;
      int last_line;
      int last_column;
      char *text;
   }
  yyltype;

#define YY_InputCParser_LTYPE yyltype
#endif
#endif
/* DEFAULT STYPE*/
      /* We used to use `unsigned long' as YY_InputCParser_STYPE on MSDOS,
	 but it seems better to be consistent.
	 Most programs should declare their own type anyway.  */

#ifndef YY_InputCParser_STYPE
#define YY_InputCParser_STYPE int
#endif
/* DEFAULT MISCELANEOUS */
#ifndef YY_InputCParser_PARSE
#define YY_InputCParser_PARSE yyparse
#endif
#ifndef YY_InputCParser_LEX
#define YY_InputCParser_LEX yylex
#endif
#ifndef YY_InputCParser_LVAL
#define YY_InputCParser_LVAL yylval
#endif
#ifndef YY_InputCParser_LLOC
#define YY_InputCParser_LLOC yylloc
#endif
#ifndef YY_InputCParser_CHAR
#define YY_InputCParser_CHAR yychar
#endif
#ifndef YY_InputCParser_NERRS
#define YY_InputCParser_NERRS yynerrs
#endif
#ifndef YY_InputCParser_DEBUG_FLAG
#define YY_InputCParser_DEBUG_FLAG yydebug
#endif
#ifndef YY_InputCParser_ERROR
#define YY_InputCParser_ERROR yyerror
#endif
#ifndef YY_InputCParser_PARSE_PARAM
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
#define YY_InputCParser_PARSE_PARAM
#ifndef YY_InputCParser_PARSE_PARAM_DEF
#define YY_InputCParser_PARSE_PARAM_DEF
#endif
#endif
#endif
#endif
#ifndef YY_InputCParser_PARSE_PARAM
#define YY_InputCParser_PARSE_PARAM void
#endif
#endif
#if YY_InputCParser_COMPATIBILITY != 0
/* backward compatibility */
#ifdef YY_InputCParser_LTYPE
#ifndef YYLTYPE
#define YYLTYPE YY_InputCParser_LTYPE
#else
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
#endif
#endif
#ifndef YYSTYPE
#define YYSTYPE YY_InputCParser_STYPE
#else
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
#endif
#ifdef YY_InputCParser_PURE
#ifndef YYPURE
#define YYPURE YY_InputCParser_PURE
#endif
#endif
#ifdef YY_InputCParser_DEBUG
#ifndef YYDEBUG
#define YYDEBUG YY_InputCParser_DEBUG 
#endif
#endif
#ifndef YY_InputCParser_ERROR_VERBOSE
#ifdef YYERROR_VERBOSE
#define YY_InputCParser_ERROR_VERBOSE YYERROR_VERBOSE
#endif
#endif
#ifndef YY_InputCParser_LSP_NEEDED
#ifdef YYLSP_NEEDED
#define YY_InputCParser_LSP_NEEDED YYLSP_NEEDED
#endif
#endif
#endif
#ifndef YY_USE_CLASS
/* TOKEN C */
#define	INTEGER_TOKEN	258
#define	DOUBLE_TOKEN	259
#define	WORD_TOKEN	260
#define	DOUBLE_VARIABLE_TOKEN	261
#define	VECTOR_DOUBLE_VARIABLE_TOKEN	262
#define	MATRIX_VARIABLE_TOKEN	263
#define	OPEN_BRACKET_TOKEN	264
#define	CLOSE_BRACKET_TOKEN	265
#define	NODES_TOKEN	266
#define	NODE_TOKEN	267
#define	NO_TOKEN	268
#define	PROPERTIES_TOKEN	269
#define	ELEMENTS_GROUP_TOKEN	270
#define	PROBLEM_TYPE_TOKEN	271
#define	PROBLEM_FREQUENCY_TOKEN	272
#define	SWEEP_FREQUENCY_TOKEN	273
#define	VE_TOKEN	274
#define	FF_TOKEN	275
#define	GRC_TOKEN	276
#define	PFF_TOKEN	277
#define	COP_TOKEN	278
#define	RWP_TOKEN	279
#define	JE_TOKEN	280
#define	PVIE_TOKEN	281
#define	PSIE_TOKEN	282
#define	PCOP_TOKEN	283
#define	PRWP_TOKEN	284
#define	PES_TOKEN	285
#define	PEP_TOKEN	286
#define	PEC_TOKEN	287
#define	PBC_TOKEN	288
#define	PMC_TOKEN	289
#define	TEC_TOKEN	290
#define	CE_TOKEN	291
#define	ELEMENT_TOKEN	292
#define	BEGIN_TOKEN	293
#define	END_TOKEN	294
#define	FIX_TOKEN	295
#define	FIXC_TOKEN	296
#define	SOURCES_TOKEN	297
#define	FOR_TOKEN	298
#define	SOLVE_TOKEN	299
#define	PRINT_TOKEN	300
#define	PRINT_ON_NODES_TOKEN	301
#define	EXECUTE_TOKEN	302
#define	TRANSIENT_TOKEN	303
#define	ALPHA_TOKEN	304
#define	CALCULATE_TOKEN	305
#define	CALCULATE_NODAL_TOKEN	306
#define	LINEAR_SOLVER_TOKEN	307
#define	SOLVING_STRATEGY_TOKEN	308
#define	EQUAL_TOKEN	309
#define	NOT_EQUAL_TOKEN	310
#define	LESS_EQUAL_TOKEN	311
#define	GREATER_EQUAL_TOKEN	312
#define	CREATE_SOLUTION_STEP_TOKEN	313
#define	CREATE_TIME_STEP_TOKEN	314

 /* #defines tokens */
#else
/* CLASS */
#ifndef YY_InputCParser_CLASS
#define YY_InputCParser_CLASS InputCParser
#endif
#ifndef YY_InputCParser_INHERIT
#define YY_InputCParser_INHERIT
#endif
#ifndef YY_InputCParser_MEMBERS
#define YY_InputCParser_MEMBERS 
#endif
#ifndef YY_InputCParser_LEX_BODY
#define YY_InputCParser_LEX_BODY  
#endif
#ifndef YY_InputCParser_ERROR_BODY
#define YY_InputCParser_ERROR_BODY  
#endif
#ifndef YY_InputCParser_CONSTRUCTOR_PARAM
#define YY_InputCParser_CONSTRUCTOR_PARAM
#endif
#ifndef YY_InputCParser_CONSTRUCTOR_CODE
#define YY_InputCParser_CONSTRUCTOR_CODE
#endif
#ifndef YY_InputCParser_CONSTRUCTOR_INIT
#define YY_InputCParser_CONSTRUCTOR_INIT
#endif
/* choose between enum and const */
#ifndef YY_InputCParser_USE_CONST_TOKEN
#define YY_InputCParser_USE_CONST_TOKEN 0
/* yes enum is more compatible with flex,  */
/* so by default we use it */ 
#endif
#if YY_InputCParser_USE_CONST_TOKEN != 0
#ifndef YY_InputCParser_ENUM_TOKEN
#define YY_InputCParser_ENUM_TOKEN yy_InputCParser_enum_token
#endif
#endif

class YY_InputCParser_CLASS YY_InputCParser_INHERIT
{
public: 
#if YY_InputCParser_USE_CONST_TOKEN != 0
/* static const int token ... */
static const int INTEGER_TOKEN;
static const int DOUBLE_TOKEN;
static const int WORD_TOKEN;
static const int DOUBLE_VARIABLE_TOKEN;
static const int VECTOR_DOUBLE_VARIABLE_TOKEN;
static const int MATRIX_VARIABLE_TOKEN;
static const int OPEN_BRACKET_TOKEN;
static const int CLOSE_BRACKET_TOKEN;
static const int NODES_TOKEN;
static const int NODE_TOKEN;
static const int NO_TOKEN;
static const int PROPERTIES_TOKEN;
static const int ELEMENTS_GROUP_TOKEN;
static const int PROBLEM_TYPE_TOKEN;
static const int PROBLEM_FREQUENCY_TOKEN;
static const int SWEEP_FREQUENCY_TOKEN;
static const int VE_TOKEN;
static const int FF_TOKEN;
static const int GRC_TOKEN;
static const int PFF_TOKEN;
static const int COP_TOKEN;
static const int RWP_TOKEN;
static const int JE_TOKEN;
static const int PVIE_TOKEN;
static const int PSIE_TOKEN;
static const int PCOP_TOKEN;
static const int PRWP_TOKEN;
static const int PES_TOKEN;
static const int PEP_TOKEN;
static const int PEC_TOKEN;
static const int PBC_TOKEN;
static const int PMC_TOKEN;
static const int TEC_TOKEN;
static const int CE_TOKEN;
static const int ELEMENT_TOKEN;
static const int BEGIN_TOKEN;
static const int END_TOKEN;
static const int FIX_TOKEN;
static const int FIXC_TOKEN;
static const int SOURCES_TOKEN;
static const int FOR_TOKEN;
static const int SOLVE_TOKEN;
static const int PRINT_TOKEN;
static const int PRINT_ON_NODES_TOKEN;
static const int EXECUTE_TOKEN;
static const int TRANSIENT_TOKEN;
static const int ALPHA_TOKEN;
static const int CALCULATE_TOKEN;
static const int CALCULATE_NODAL_TOKEN;
static const int LINEAR_SOLVER_TOKEN;
static const int SOLVING_STRATEGY_TOKEN;
static const int EQUAL_TOKEN;
static const int NOT_EQUAL_TOKEN;
static const int LESS_EQUAL_TOKEN;
static const int GREATER_EQUAL_TOKEN;
static const int CREATE_SOLUTION_STEP_TOKEN;
static const int CREATE_TIME_STEP_TOKEN;

 /* decl const */
#else
enum YY_InputCParser_ENUM_TOKEN { YY_InputCParser_NULL_TOKEN=0
	,INTEGER_TOKEN=258
	,DOUBLE_TOKEN=259
	,WORD_TOKEN=260
	,DOUBLE_VARIABLE_TOKEN=261
	,VECTOR_DOUBLE_VARIABLE_TOKEN=262
	,MATRIX_VARIABLE_TOKEN=263
	,OPEN_BRACKET_TOKEN=264
	,CLOSE_BRACKET_TOKEN=265
	,NODES_TOKEN=266
	,NODE_TOKEN=267
	,NO_TOKEN=268
	,PROPERTIES_TOKEN=269
	,ELEMENTS_GROUP_TOKEN=270
	,PROBLEM_TYPE_TOKEN=271
	,PROBLEM_FREQUENCY_TOKEN=272
	,SWEEP_FREQUENCY_TOKEN=273
	,VE_TOKEN=274
	,FF_TOKEN=275
	,GRC_TOKEN=276
	,PFF_TOKEN=277
	,COP_TOKEN=278
	,RWP_TOKEN=279
	,JE_TOKEN=280
	,PVIE_TOKEN=281
	,PSIE_TOKEN=282
	,PCOP_TOKEN=283
	,PRWP_TOKEN=284
	,PES_TOKEN=285
	,PEP_TOKEN=286
	,PEC_TOKEN=287
	,PBC_TOKEN=288
	,PMC_TOKEN=289
	,TEC_TOKEN=290
	,CE_TOKEN=291
	,ELEMENT_TOKEN=292
	,BEGIN_TOKEN=293
	,END_TOKEN=294
	,FIX_TOKEN=295
	,FIXC_TOKEN=296
	,SOURCES_TOKEN=297
	,FOR_TOKEN=298
	,SOLVE_TOKEN=299
	,PRINT_TOKEN=300
	,PRINT_ON_NODES_TOKEN=301
	,EXECUTE_TOKEN=302
	,TRANSIENT_TOKEN=303
	,ALPHA_TOKEN=304
	,CALCULATE_TOKEN=305
	,CALCULATE_NODAL_TOKEN=306
	,LINEAR_SOLVER_TOKEN=307
	,SOLVING_STRATEGY_TOKEN=308
	,EQUAL_TOKEN=309
	,NOT_EQUAL_TOKEN=310
	,LESS_EQUAL_TOKEN=311
	,GREATER_EQUAL_TOKEN=312
	,CREATE_SOLUTION_STEP_TOKEN=313
	,CREATE_TIME_STEP_TOKEN=314

 /* enum token */
     }; /* end of enum declaration */
#endif
public:
 int YY_InputCParser_PARSE (YY_InputCParser_PARSE_PARAM);
 virtual void YY_InputCParser_ERROR(char *msg) YY_InputCParser_ERROR_BODY;
#ifdef YY_InputCParser_PURE
#ifdef YY_InputCParser_LSP_NEEDED
 virtual int  YY_InputCParser_LEX (YY_InputCParser_STYPE *YY_InputCParser_LVAL,YY_InputCParser_LTYPE *YY_InputCParser_LLOC) YY_InputCParser_LEX_BODY;
#else
 virtual int  YY_InputCParser_LEX (YY_InputCParser_STYPE *YY_InputCParser_LVAL) YY_InputCParser_LEX_BODY;
#endif
#else
 virtual int YY_InputCParser_LEX() YY_InputCParser_LEX_BODY;
 YY_InputCParser_STYPE YY_InputCParser_LVAL;
#ifdef YY_InputCParser_LSP_NEEDED
 YY_InputCParser_LTYPE YY_InputCParser_LLOC;
#endif
 int   YY_InputCParser_NERRS;
 int    YY_InputCParser_CHAR;
#endif
#if YY_InputCParser_DEBUG != 0
 int YY_InputCParser_DEBUG_FLAG;   /*  nonzero means print parse trace     */
#endif
public:
 YY_InputCParser_CLASS(YY_InputCParser_CONSTRUCTOR_PARAM);
public:
 YY_InputCParser_MEMBERS 
};
/* other declare folow */
#if YY_InputCParser_USE_CONST_TOKEN != 0
const int YY_InputCParser_CLASS::INTEGER_TOKEN=258;
const int YY_InputCParser_CLASS::DOUBLE_TOKEN=259;
const int YY_InputCParser_CLASS::WORD_TOKEN=260;
const int YY_InputCParser_CLASS::DOUBLE_VARIABLE_TOKEN=261;
const int YY_InputCParser_CLASS::VECTOR_DOUBLE_VARIABLE_TOKEN=262;
const int YY_InputCParser_CLASS::MATRIX_VARIABLE_TOKEN=263;
const int YY_InputCParser_CLASS::OPEN_BRACKET_TOKEN=264;
const int YY_InputCParser_CLASS::CLOSE_BRACKET_TOKEN=265;
const int YY_InputCParser_CLASS::NODES_TOKEN=266;
const int YY_InputCParser_CLASS::NODE_TOKEN=267;
const int YY_InputCParser_CLASS::NO_TOKEN=268;
const int YY_InputCParser_CLASS::PROPERTIES_TOKEN=269;
const int YY_InputCParser_CLASS::ELEMENTS_GROUP_TOKEN=270;
const int YY_InputCParser_CLASS::PROBLEM_TYPE_TOKEN=271;
const int YY_InputCParser_CLASS::PROBLEM_FREQUENCY_TOKEN=272;
const int YY_InputCParser_CLASS::SWEEP_FREQUENCY_TOKEN=273;
const int YY_InputCParser_CLASS::VE_TOKEN=274;
const int YY_InputCParser_CLASS::FF_TOKEN=275;
const int YY_InputCParser_CLASS::GRC_TOKEN=276;
const int YY_InputCParser_CLASS::PFF_TOKEN=277;
const int YY_InputCParser_CLASS::COP_TOKEN=278;
const int YY_InputCParser_CLASS::RWP_TOKEN=279;
const int YY_InputCParser_CLASS::JE_TOKEN=280;
const int YY_InputCParser_CLASS::PVIE_TOKEN=281;
const int YY_InputCParser_CLASS::PSIE_TOKEN=282;
const int YY_InputCParser_CLASS::PCOP_TOKEN=283;
const int YY_InputCParser_CLASS::PRWP_TOKEN=284;
const int YY_InputCParser_CLASS::PES_TOKEN=285;
const int YY_InputCParser_CLASS::PEP_TOKEN=286;
const int YY_InputCParser_CLASS::PEC_TOKEN=287;
const int YY_InputCParser_CLASS::PBC_TOKEN=288;
const int YY_InputCParser_CLASS::PMC_TOKEN=289;
const int YY_InputCParser_CLASS::TEC_TOKEN=290;
const int YY_InputCParser_CLASS::CE_TOKEN=291;
const int YY_InputCParser_CLASS::ELEMENT_TOKEN=292;
const int YY_InputCParser_CLASS::BEGIN_TOKEN=293;
const int YY_InputCParser_CLASS::END_TOKEN=294;
const int YY_InputCParser_CLASS::FIX_TOKEN=295;
const int YY_InputCParser_CLASS::FIXC_TOKEN=296;
const int YY_InputCParser_CLASS::SOURCES_TOKEN=297;
const int YY_InputCParser_CLASS::FOR_TOKEN=298;
const int YY_InputCParser_CLASS::SOLVE_TOKEN=299;
const int YY_InputCParser_CLASS::PRINT_TOKEN=300;
const int YY_InputCParser_CLASS::PRINT_ON_NODES_TOKEN=301;
const int YY_InputCParser_CLASS::EXECUTE_TOKEN=302;
const int YY_InputCParser_CLASS::TRANSIENT_TOKEN=303;
const int YY_InputCParser_CLASS::ALPHA_TOKEN=304;
const int YY_InputCParser_CLASS::CALCULATE_TOKEN=305;
const int YY_InputCParser_CLASS::CALCULATE_NODAL_TOKEN=306;
const int YY_InputCParser_CLASS::LINEAR_SOLVER_TOKEN=307;
const int YY_InputCParser_CLASS::SOLVING_STRATEGY_TOKEN=308;
const int YY_InputCParser_CLASS::EQUAL_TOKEN=309;
const int YY_InputCParser_CLASS::NOT_EQUAL_TOKEN=310;
const int YY_InputCParser_CLASS::LESS_EQUAL_TOKEN=311;
const int YY_InputCParser_CLASS::GREATER_EQUAL_TOKEN=312;
const int YY_InputCParser_CLASS::CREATE_SOLUTION_STEP_TOKEN=313;
const int YY_InputCParser_CLASS::CREATE_TIME_STEP_TOKEN=314;

 /* const YY_InputCParser_CLASS::token */
#endif
/*apres const  */
YY_InputCParser_CLASS::YY_InputCParser_CLASS(YY_InputCParser_CONSTRUCTOR_PARAM) YY_InputCParser_CONSTRUCTOR_INIT
{
#if YY_InputCParser_DEBUG != 0
YY_InputCParser_DEBUG_FLAG=0;
#endif
YY_InputCParser_CONSTRUCTOR_CODE;
};
#endif


#define	YYFINAL		378
#define	YYFLAG		32768
#define	YYNTBASE	74

#define YYTRANSLATE(x) ((unsigned)(x) <= 314 ? yytranslate[x] : 127)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    64,
    66,     2,    72,    65,    73,    67,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,    60,    68,
    63,    69,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
    70,     2,    71,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,    62,     2,    61,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
    26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
    36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
    46,    47,    48,    49,    50,    51,    52,    53,    54,    55,
    56,    57,    58,    59
};

#if YY_InputCParser_DEBUG != 0
static const short yyprhs[] = {     0,
     0,     1,     4,     7,     9,    11,    13,    15,    17,    19,
    21,    23,    25,    27,    29,    31,    33,    35,    37,    39,
    41,    43,    45,    47,    49,    51,    53,    55,    57,    59,
    61,    63,    65,    67,    69,    71,    73,    75,    77,    79,
    81,    84,    87,    90,   101,   107,   113,   119,   123,   127,
   131,   142,   155,   168,   179,   192,   205,   214,   225,   236,
   247,   258,   269,   280,   287,   294,   301,   308,   317,   326,
   332,   344,   354,   361,   366,   371,   376,   381,   386,   391,
   396,   401,   406,   411,   416,   421,   425,   429,   432,   435,
   439,   443,   447,   451,   455,   459,   463,   467,   471,   475,
   479,   483,   487,   489,   493,   497,   501,   505,   509,   513,
   515,   517,   521,   525,   527,   530,   533,   535
};

static const short yyrhs[] = {    -1,
    74,    75,     0,    76,    60,     0,    77,     0,   105,     0,
   115,     0,   116,     0,    79,     0,    80,     0,    81,     0,
    82,     0,    83,     0,    84,     0,    85,     0,    86,     0,
    88,     0,    89,     0,    87,     0,    90,     0,    91,     0,
    92,     0,    94,     0,    93,     0,    95,     0,    96,     0,
    97,     0,    98,     0,    99,     0,   100,     0,   101,     0,
   102,     0,   103,     0,   106,     0,   107,     0,   108,     0,
   109,     0,   110,     0,   104,     0,   111,     0,   112,     0,
    78,    61,     0,    62,    75,     0,    78,    75,     0,   114,
    63,     5,    64,   125,    65,   125,    65,   125,    66,     0,
   113,    67,     6,    63,   125,     0,   113,    67,     8,    63,
   117,     0,   113,    67,     7,    63,   119,     0,    15,    63,
     5,     0,    16,    63,     5,     0,    17,    63,   125,     0,
    18,    63,     5,    64,   125,    65,   125,    65,   125,    66,
     0,    19,    64,   126,    65,   126,    65,   126,    65,   126,
    65,   126,    66,     0,    25,    64,   126,    65,   126,    65,
   126,    65,   126,    65,   126,    66,     0,    27,    64,   126,
    65,   126,    65,   126,    65,   126,    66,     0,    26,    64,
   126,    65,   126,    65,   126,    65,   126,    65,   126,    66,
     0,    33,    64,   126,    65,   126,    65,   126,    65,   126,
    65,   126,    66,     0,    20,    64,   126,    65,   126,    65,
   126,    66,     0,    21,    64,   126,    65,   126,    65,   126,
    65,   126,    66,     0,    22,    64,   126,    65,   126,    65,
   126,    65,   126,    66,     0,    24,    64,   126,    65,   126,
    65,   126,    65,   126,    66,     0,    23,    64,   126,    65,
   126,    65,   126,    65,   126,    66,     0,    28,    64,   126,
    65,   126,    65,   126,    65,   126,    66,     0,    29,    64,
   126,    65,   126,    65,   126,    65,   126,    66,     0,    32,
    63,     5,    64,   120,    66,     0,    34,    63,     5,    64,
   120,    66,     0,    35,    63,     5,    64,   120,    66,     0,
    36,    63,     5,    64,   120,    66,     0,   114,    67,     6,
    67,    40,    64,   125,    66,     0,   114,    67,     6,    67,
    41,    64,   119,    66,     0,   114,    67,     6,    63,   124,
     0,    52,     5,    63,     5,    64,   126,    65,   126,    65,
   125,    66,     0,    43,    64,   115,    60,   116,    60,    76,
    66,    75,     0,     5,    67,    44,    64,     5,    66,     0,
    45,    64,     6,    66,     0,    45,    64,     7,    66,     0,
    45,    64,     8,    66,     0,    46,    64,     6,    66,     0,
    46,    64,     7,    66,     0,    46,    64,     8,    66,     0,
    50,    64,     6,    66,     0,    50,    64,     7,    66,     0,
    50,    64,     8,    66,     0,    51,    64,     6,    66,     0,
    51,    64,     7,    66,     0,    51,    64,     8,    66,     0,
    58,    64,    66,     0,    59,    64,    66,     0,    14,   121,
     0,    13,   121,     0,     6,    63,   124,     0,     6,    68,
   124,     0,     6,    69,   124,     0,     6,    56,   124,     0,
     6,    57,   124,     0,     6,    54,   124,     0,     6,    55,
   124,     0,    70,   119,    71,     0,    70,   118,    71,     0,
   119,    65,   119,     0,   118,    65,   119,     0,    70,   125,
    71,     0,    70,   122,    71,     0,   121,     0,    70,   123,
    71,     0,    70,   126,    71,     0,   125,    65,   125,     0,
   122,    65,   125,     0,   126,    65,   126,     0,   123,    65,
   126,     0,   125,     0,     6,     0,   124,    72,     6,     0,
   124,    72,   125,     0,     4,     0,    72,   125,     0,    73,
   125,     0,   126,     0,     3,     0
};

#endif

#if YY_InputCParser_DEBUG != 0
static const short yyrline[] = { 0,
   950,   951,   956,   957,   958,   961,   962,   963,   964,   965,
   966,   967,   968,   969,   970,   971,   972,   973,   974,   975,
   976,   977,   978,   979,   980,   981,   982,   983,   984,   985,
   986,   987,   988,   989,   990,   991,   992,   993,   994,   995,
   998,  1003,  1008,  1014,  1027,  1032,  1038,  1045,  1050,  1055,
  1061,  1069,  1085,  1101,  1115,  1131,  1147,  1159,  1173,  1187,
  1201,  1215,  1229,  1243,  1251,  1259,  1267,  1275,  1281,  1288,
  1294,  1302,  1311,  1316,  1320,  1324,  1329,  1333,  1337,  1342,
  1346,  1350,  1355,  1359,  1363,  1368,  1373,  1378,  1383,  1385,
  1390,  1394,  1398,  1402,  1406,  1410,  1415,  1420,  1426,  1435,
  1442,  1447,  1449,  1454,  1456,  1458,  1464,  1470,  1476,  1482,
  1486,  1490,  1494,  1499,  1500,  1501,  1502,  1504
};

static const char * const yytname[] = {   "$","error","$illegal.","INTEGER_TOKEN",
"DOUBLE_TOKEN","WORD_TOKEN","DOUBLE_VARIABLE_TOKEN","VECTOR_DOUBLE_VARIABLE_TOKEN",
"MATRIX_VARIABLE_TOKEN","OPEN_BRACKET_TOKEN","CLOSE_BRACKET_TOKEN","NODES_TOKEN",
"NODE_TOKEN","NO_TOKEN","PROPERTIES_TOKEN","ELEMENTS_GROUP_TOKEN","PROBLEM_TYPE_TOKEN",
"PROBLEM_FREQUENCY_TOKEN","SWEEP_FREQUENCY_TOKEN","VE_TOKEN","FF_TOKEN","GRC_TOKEN",
"PFF_TOKEN","COP_TOKEN","RWP_TOKEN","JE_TOKEN","PVIE_TOKEN","PSIE_TOKEN","PCOP_TOKEN",
"PRWP_TOKEN","PES_TOKEN","PEP_TOKEN","PEC_TOKEN","PBC_TOKEN","PMC_TOKEN","TEC_TOKEN",
"CE_TOKEN","ELEMENT_TOKEN","BEGIN_TOKEN","END_TOKEN","FIX_TOKEN","FIXC_TOKEN",
"SOURCES_TOKEN","FOR_TOKEN","SOLVE_TOKEN","PRINT_TOKEN","PRINT_ON_NODES_TOKEN",
"EXECUTE_TOKEN","TRANSIENT_TOKEN","ALPHA_TOKEN","CALCULATE_TOKEN","CALCULATE_NODAL_TOKEN",
"LINEAR_SOLVER_TOKEN","SOLVING_STRATEGY_TOKEN","EQUAL_TOKEN","NOT_EQUAL_TOKEN",
"LESS_EQUAL_TOKEN","GREATER_EQUAL_TOKEN","CREATE_SOLUTION_STEP_TOKEN","CREATE_TIME_STEP_TOKEN",
"';'","'}'","'{'","'='","'('","','","')'","'.'","'<'","'>'","'['","']'","'+'",
"'-'","statements","statement","expresion","block","block_generating","node_generating",
"properties_adding","set_elements_group","set_problem_type","set_problem_frequency",
"set_sweep_frequency","volume_element_generating","source_element_generating",
"srfIntg_element_generating","volIntg_element_generating","PBC_element_generating",
"farfield_element_generating","genRobin_element_generating","plasmafarfield_element_generating",
"rwportTE10_element_generating","coaxportTEM_element_generating","projection_coaxTEM_generating",
"projection_rwTE10_generating","pec_element_generating","pmc_element_generating",
"tec_element_generating","contact_generation","nodes_variables_fixing","nodes_cvariables_fixing",
"nodes_variables_setting","linear_solver_generating","for_loop","solve_process",
"print_process","print_on_nodes_process","calculate_process","calculate_nodal_process",
"step_creating","time_step_creating","properties_array","nodes_array","assignment_expresion",
"logical_expresion","matrix","vector_double_sequence","vector_double","vector_integer",
"integer_index","double_sequence","integer_sequence","double_variable_expresion",
"double_expresion","integer_expresion",""
};
#endif

static const short yyr1[] = {     0,
    74,    74,    75,    75,    75,    76,    76,    76,    76,    76,
    76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
    76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
    76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
    77,    78,    78,    79,    80,    80,    80,    81,    82,    83,
    84,    85,    86,    87,    88,    89,    90,    91,    92,    93,
    94,    95,    96,    97,    98,    99,   100,   101,   102,   103,
   104,   105,   106,   107,   107,   107,   108,   108,   108,   109,
   109,   109,   110,   110,   110,   111,   112,   113,   114,   115,
   116,   116,   116,   116,   116,   116,   117,   117,   118,   118,
   119,   119,   120,   120,   121,   122,   122,   123,   123,   124,
   124,   124,   124,   125,   125,   125,   125,   126
};

static const short yyr2[] = {     0,
     0,     2,     2,     1,     1,     1,     1,     1,     1,     1,
     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
     2,     2,     2,    10,     5,     5,     5,     3,     3,     3,
    10,    12,    12,    10,    12,    12,     8,    10,    10,    10,
    10,    10,    10,     6,     6,     6,     6,     8,     8,     5,
    11,     9,     6,     4,     4,     4,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     3,     3,     2,     2,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     3,     3,     1,     3,     3,     3,     3,     3,     3,     1,
     1,     3,     3,     1,     2,     2,     1,     1
};

static const short yydefact[] = {     1,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     2,     0,     4,     0,     8,     9,
    10,    11,    12,    13,    14,    15,    18,    16,    17,    19,
    20,    21,    23,    22,    24,    25,    26,    27,    28,    29,
    30,    31,    32,    38,     5,    33,    34,    35,    36,    37,
    39,    40,     0,     0,     6,     7,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    89,    88,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,    42,     3,    41,    43,     0,
     0,     0,     0,   118,   114,   111,     0,     0,    95,   110,
   117,    96,    93,    94,    90,    91,    92,     0,    48,    49,
    50,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,    86,    87,     0,     0,     0,     0,     0,
     0,   115,   116,     0,   105,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,    74,    75,    76,    77,    78,    79,    80,
    81,    82,    83,    84,    85,     0,     0,     0,     0,     0,
     0,     0,     0,   112,   113,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,   103,
     0,     0,     0,     0,     0,     0,     0,    45,     0,    47,
     0,    46,     0,    70,     0,     0,    73,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    64,     0,    65,    66,    67,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,   104,     0,
     0,     0,     0,     0,   102,     0,   101,     0,    98,     0,
    97,     0,     0,     0,     0,     0,    57,     0,     0,     0,
     0,     0,     0,     0,     0,     0,   109,   108,     0,     0,
     0,   107,   106,   100,    99,     0,    68,    69,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
    72,     0,     0,    51,     0,    58,    59,    61,    60,     0,
     0,    54,    62,    63,     0,     0,    44,     0,     0,     0,
     0,    71,    52,    53,    55,    56,     0,     0
};

static const short yydefgoto[] = {     1,
    35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
    45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
    55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
    65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
    75,    76,   252,   281,   250,   239,   240,   279,   270,   129,
   130,   131
};

static const short yypact[] = {-32768,
   284,   -48,   -21,   -40,   -40,    11,    16,    36,    37,   -23,
    50,    52,    54,    55,    56,    57,    58,    59,    61,    63,
    39,    65,    67,    68,    69,    70,    71,    72,    73,    74,
    96,    76,    77,   384,-32768,     2,-32768,   334,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,    66,    -2,-32768,-32768,    80,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,   114,-32768,-32768,   134,   137,     3,
   138,   114,   114,   114,   114,   114,   114,   114,   114,   114,
   114,   114,   139,   114,   140,   141,   142,   143,    38,    60,
    75,    78,    85,    84,    86,-32768,-32768,-32768,-32768,    81,
   146,   147,    92,-32768,-32768,-32768,     3,     3,    83,-32768,
-32768,    83,    83,    83,    83,    83,    83,    87,-32768,-32768,
-32768,    93,    95,   100,   101,   102,   103,   104,   105,   106,
   107,   108,   109,    97,   121,   123,   124,   126,   129,   133,
   128,   130,   131,   132,   135,   136,   144,   145,   149,   151,
   152,   154,   190,-32768,-32768,   150,   158,   159,   148,    10,
   194,-32768,-32768,    25,-32768,     3,   114,   114,   114,   114,
   114,   114,   114,   114,   114,   114,   114,   153,   114,   153,
   153,   153,   197,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,   160,     3,   157,   161,     3,
    -1,    53,   162,-32768,-32768,   164,   165,   168,   178,   179,
   182,   184,   185,   186,   187,   188,   189,   114,   192,-32768,
   191,   193,   195,   196,   -17,   200,   114,-32768,     3,-32768,
   157,-32768,   198,    83,   201,   202,-32768,     3,   114,   114,
   114,   114,   114,   114,   114,   114,   114,   114,   114,   -22,
   -15,-32768,   114,-32768,-32768,-32768,   434,   199,   -12,   -11,
    -8,    -7,     3,     3,   157,   203,   204,   206,   205,   208,
   209,   210,   214,   215,   217,   218,   220,   114,-32768,   114,
   221,   222,   114,     3,-32768,     3,-32768,   157,-32768,   157,
-32768,   226,   227,   229,     3,   114,-32768,   114,   114,   114,
   114,   114,   114,   114,   114,   114,-32768,-32768,   114,   384,
   249,-32768,-32768,-32768,-32768,     3,-32768,-32768,   255,   250,
   256,   257,   258,   259,   261,   263,   265,   266,   267,   272,
-32768,     3,   275,-32768,   114,-32768,-32768,-32768,-32768,   114,
   114,-32768,-32768,-32768,   114,   278,-32768,   279,   298,   299,
   305,-32768,-32768,-32768,-32768,-32768,   207,-32768
};

static const short yypgoto[] = {-32768,
   -34,   -73,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,-32768,
    98,     5,-32768,-32768,  -230,  -110,    91,-32768,-32768,   -57,
   -58,   -84
};


#define	YYLAST		493


static const short yytable[] = {   116,
   138,   124,   125,   119,   126,   124,   125,   143,   144,   145,
   146,   147,   148,   149,   150,   151,   152,   153,    77,   155,
   282,   132,   133,   134,   135,   136,   137,   124,   125,    85,
   224,   141,    78,    79,    80,    81,    78,    79,    80,    81,
    92,    82,   298,   161,   162,   163,    83,    84,   299,   300,
    83,    84,   304,   306,   314,   185,   308,   310,   305,   307,
   121,   117,   309,   311,   122,   164,   165,   166,   182,   183,
   127,   128,   221,    88,   127,   128,   222,   334,    89,   335,
   167,   168,   169,   170,   171,   172,   176,   177,   178,   242,
   243,   244,   255,   256,    86,    87,   127,   128,    90,    91,
   113,   103,   227,   228,   229,   230,   231,   232,   233,   234,
   235,   236,   237,    93,   241,    94,   124,    95,    96,    97,
    98,    99,   100,   123,   101,   225,   102,   226,   104,   105,
   106,   107,   120,   108,   109,   110,   111,   112,   139,   114,
   115,   140,   142,   154,   156,   157,   158,   173,   159,   174,
   179,   175,   180,   271,   184,   181,   186,   185,   248,   187,
   198,   253,   278,   254,   188,   189,   190,   191,   192,   193,
   194,   195,   196,   197,   287,   288,   289,   290,   291,   292,
   293,   294,   295,   296,   297,   199,   200,   201,   301,   202,
   280,    82,   203,   204,   216,   205,   206,   207,   223,   286,
   208,   209,   245,   302,     0,   160,   378,   246,     0,   210,
   211,   220,   217,   327,   212,   328,   213,   214,   331,   215,
   218,   219,   238,   247,   312,   313,   249,   257,   258,   259,
   251,   340,   260,   341,   342,   343,   344,   345,   346,   347,
   348,   349,   261,   262,   350,   332,   263,   333,   264,   265,
   266,   267,   268,   269,     0,   273,   339,   272,   274,   277,
   275,   276,   283,   303,   284,   285,     0,   315,   316,   318,
   368,   317,   319,   320,   321,   369,   370,   353,   322,   323,
   371,   324,   325,   377,   326,   329,     0,   330,     2,     3,
   336,     0,   337,   366,   338,   351,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,   352,   355,    21,    22,    23,    24,    25,
   354,   356,   357,   358,   359,   360,    26,   361,    27,    28,
   362,   363,   364,    29,    30,    31,   365,     0,     2,     3,
   367,    32,    33,   372,   373,    34,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,   374,   375,    21,    22,    23,    24,    25,
   376,     0,     0,     0,     0,     0,    26,     0,    27,    28,
     0,     0,     0,    29,    30,    31,     0,     0,     2,     3,
     0,    32,    33,     0,   118,    34,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,     0,     0,    21,    22,    23,    24,    25,
     0,     0,     0,     0,     0,     0,    26,     0,    27,    28,
     0,     0,     0,    29,    30,    31,     0,     0,     2,     3,
     0,    32,    33,     0,     0,    34,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,     0,     0,    21,    22,    23,    24,    25,
     0,     0,     0,     0,     0,     0,     0,     0,    27,    28,
     0,     0,     0,    29,    30,    31,     0,     0,     0,     0,
     0,    32,    33
};

static const short yycheck[] = {    34,
    85,     3,     4,    38,     6,     3,     4,    92,    93,    94,
    95,    96,    97,    98,    99,   100,   101,   102,    67,   104,
   251,    79,    80,    81,    82,    83,    84,     3,     4,    70,
     6,    90,    54,    55,    56,    57,    54,    55,    56,    57,
    64,    63,    65,     6,     7,     8,    68,    69,    71,    65,
    68,    69,    65,    65,   285,    71,    65,    65,    71,    71,
    63,    60,    71,    71,    67,     6,     7,     8,   127,   128,
    72,    73,    63,    63,    72,    73,    67,   308,    63,   310,
     6,     7,     8,     6,     7,     8,     6,     7,     8,   200,
   201,   202,    40,    41,     4,     5,    72,    73,    63,    63,
     5,    63,   187,   188,   189,   190,   191,   192,   193,   194,
   195,   196,   197,    64,   199,    64,     3,    64,    64,    64,
    64,    64,    64,    44,    64,   184,    64,   186,    64,    63,
    63,    63,    67,    64,    64,    64,    64,    64,     5,    64,
    64,     5,     5,     5,     5,     5,     5,    63,     6,    66,
     5,    66,     6,   238,    72,    64,    64,    71,   217,    65,
    64,   220,   247,   221,    65,    65,    65,    65,    65,    65,
    65,    65,    65,    65,   259,   260,   261,   262,   263,   264,
   265,   266,   267,   268,   269,    65,    64,    64,   273,    64,
   249,    63,    60,    66,     5,    66,    66,    66,     5,   258,
    66,    66,     6,   277,    -1,   108,     0,   203,    -1,    66,
    66,    64,    63,   298,    66,   300,    66,    66,   303,    66,
    63,    63,    70,    64,   283,   284,    70,    66,    65,    65,
    70,   316,    65,   318,   319,   320,   321,   322,   323,   324,
   325,   326,    65,    65,   329,   304,    65,   306,    65,    65,
    65,    65,    65,    65,    -1,    65,   315,    66,    66,    60,
    66,    66,    65,    65,    64,    64,    -1,    65,    65,    65,
   355,    66,    65,    65,    65,   360,   361,   336,    65,    65,
   365,    65,    65,     0,    65,    65,    -1,    66,     5,     6,
    65,    -1,    66,   352,    66,   330,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    27,    28,    29,    65,    65,    32,    33,    34,    35,    36,
    66,    66,    66,    66,    66,    65,    43,    65,    45,    46,
    66,    66,    66,    50,    51,    52,    65,    -1,     5,     6,
    66,    58,    59,    66,    66,    62,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    27,    28,    29,    66,    66,    32,    33,    34,    35,    36,
    66,    -1,    -1,    -1,    -1,    -1,    43,    -1,    45,    46,
    -1,    -1,    -1,    50,    51,    52,    -1,    -1,     5,     6,
    -1,    58,    59,    -1,    61,    62,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    27,    28,    29,    -1,    -1,    32,    33,    34,    35,    36,
    -1,    -1,    -1,    -1,    -1,    -1,    43,    -1,    45,    46,
    -1,    -1,    -1,    50,    51,    52,    -1,    -1,     5,     6,
    -1,    58,    59,    -1,    -1,    62,    13,    14,    15,    16,
    17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
    27,    28,    29,    -1,    -1,    32,    33,    34,    35,    36,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    45,    46,
    -1,    -1,    -1,    50,    51,    52,    -1,    -1,    -1,    -1,
    -1,    58,    59
};
 /* fattrs + tables */

/* parser code folow  */


/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: dollar marks section change
   the next  is replaced by the list of actions, each action
   as one case of the switch.  */ 

#if YY_InputCParser_USE_GOTO != 0
/* 
 SUPRESSION OF GOTO : on some C++ compiler (sun c++)
  the goto is strictly forbidden if any constructor/destructor
  is used in the whole function (very stupid isn't it ?)
 so goto are to be replaced with a 'while/switch/case construct'
 here are the macro to keep some apparent compatibility
*/
#define YYGOTO(lb) {yy_gotostate=lb;continue;}
#define YYBEGINGOTO  enum yy_labels yy_gotostate=yygotostart; \
                     for(;;) switch(yy_gotostate) { case yygotostart: {
#define YYLABEL(lb) } case lb: {
#define YYENDGOTO } } 
#define YYBEGINDECLARELABEL enum yy_labels {yygotostart
#define YYDECLARELABEL(lb) ,lb
#define YYENDDECLARELABEL  };
#else
/* macro to keep goto */
#define YYGOTO(lb) goto lb
#define YYBEGINGOTO 
#define YYLABEL(lb) lb:
#define YYENDGOTO
#define YYBEGINDECLARELABEL 
#define YYDECLARELABEL(lb)
#define YYENDDECLARELABEL 
#endif
/* LABEL DECLARATION */
YYBEGINDECLARELABEL
  YYDECLARELABEL(yynewstate)
  YYDECLARELABEL(yybackup)
/* YYDECLARELABEL(yyresume) */
  YYDECLARELABEL(yydefault)
  YYDECLARELABEL(yyreduce)
  YYDECLARELABEL(yyerrlab)   /* here on detecting error */
  YYDECLARELABEL(yyerrlab1)   /* here on error raised explicitly by an action */
  YYDECLARELABEL(yyerrdefault)  /* current state does not do anything special for the error token. */
  YYDECLARELABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */
  YYDECLARELABEL(yyerrhandle)  
YYENDDECLARELABEL
/* ALLOCA SIMULATION */
/* __HAVE_NO_ALLOCA */
#ifdef __HAVE_NO_ALLOCA
int __alloca_free_ptr(char *ptr,char *ref)
{if(ptr!=ref) free(ptr);
 return 0;}

#define __ALLOCA_alloca(size) malloc(size)
#define __ALLOCA_free(ptr,ref) __alloca_free_ptr((char *)ptr,(char *)ref)

#ifdef YY_InputCParser_LSP_NEEDED
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		    __ALLOCA_free(yyls,yylsa)+\
		   (num))
#else
#define __ALLOCA_return(num) \
            return( __ALLOCA_free(yyss,yyssa)+\
		    __ALLOCA_free(yyvs,yyvsa)+\
		   (num))
#endif
#else
#define __ALLOCA_return(num) return(num)
#define __ALLOCA_alloca(size) alloca(size)
#define __ALLOCA_free(ptr,ref) 
#endif

/* ENDALLOCA SIMULATION */

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (YY_InputCParser_CHAR = YYEMPTY)
#define YYEMPTY         -2
#define YYEOF           0
#define YYACCEPT        __ALLOCA_return(0)
#define YYABORT         __ALLOCA_return(1)
#define YYERROR         YYGOTO(yyerrlab1)
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL          YYGOTO(yyerrlab)
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do                                                              \
  if (YY_InputCParser_CHAR == YYEMPTY && yylen == 1)                               \
    { YY_InputCParser_CHAR = (token), YY_InputCParser_LVAL = (value);                 \
      yychar1 = YYTRANSLATE (YY_InputCParser_CHAR);                                \
      YYPOPSTACK;                                               \
      YYGOTO(yybackup);                                            \
    }                                                           \
  else                                                          \
    { YY_InputCParser_ERROR ("syntax error: cannot back up"); YYERROR; }   \
while (0)

#define YYTERROR        1
#define YYERRCODE       256

#ifndef YY_InputCParser_PURE
/* UNPURE */
#define YYLEX           YY_InputCParser_LEX()
#ifndef YY_USE_CLASS
/* If nonreentrant, and not class , generate the variables here */
int     YY_InputCParser_CHAR;                      /*  the lookahead symbol        */
YY_InputCParser_STYPE      YY_InputCParser_LVAL;              /*  the semantic value of the */
				/*  lookahead symbol    */
int YY_InputCParser_NERRS;                 /*  number of parse errors so far */
#ifdef YY_InputCParser_LSP_NEEDED
YY_InputCParser_LTYPE YY_InputCParser_LLOC;   /*  location data for the lookahead     */
			/*  symbol                              */
#endif
#endif


#else
/* PURE */
#ifdef YY_InputCParser_LSP_NEEDED
#define YYLEX           YY_InputCParser_LEX(&YY_InputCParser_LVAL, &YY_InputCParser_LLOC)
#else
#define YYLEX           YY_InputCParser_LEX(&YY_InputCParser_LVAL)
#endif
#endif
#ifndef YY_USE_CLASS
#if YY_InputCParser_DEBUG != 0
int YY_InputCParser_DEBUG_FLAG;                    /*  nonzero means print parse trace     */
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif
#endif



/*  YYINITDEPTH indicates the initial size of the parser's stacks       */

#ifndef YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif


#if __GNUC__ > 1                /* GNU C and GNU C++ define this.  */
#define __yy_bcopy(FROM,TO,COUNT)       __builtin_memcpy(TO,FROM,COUNT)
#else                           /* not GNU C or C++ */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */

#ifdef __cplusplus
static void __yy_bcopy (char *from, char *to, int count)
#else
#ifdef __STDC__
static void __yy_bcopy (char *from, char *to, int count)
#else
static void __yy_bcopy (from, to, count)
     char *from;
     char *to;
     int count;
#endif
#endif
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}
#endif

int
#ifdef YY_USE_CLASS
 YY_InputCParser_CLASS::
#endif
     YY_InputCParser_PARSE(YY_InputCParser_PARSE_PARAM)
#ifndef __STDC__
#ifndef __cplusplus
#ifndef YY_USE_CLASS
/* parameter definition without protypes */
YY_InputCParser_PARSE_PARAM_DEF
#endif
#endif
#endif
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YY_InputCParser_STYPE *yyvsp;
  int yyerrstatus;      /*  number of tokens to shift before error messages enabled */
  int yychar1=0;          /*  lookahead token as an internal (translated) token number */

  short yyssa[YYINITDEPTH];     /*  the state stack                     */
  YY_InputCParser_STYPE yyvsa[YYINITDEPTH];        /*  the semantic value stack            */

  short *yyss = yyssa;          /*  refer to the stacks thru separate pointers */
  YY_InputCParser_STYPE *yyvs = yyvsa;     /*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YY_InputCParser_LSP_NEEDED
  YY_InputCParser_LTYPE yylsa[YYINITDEPTH];        /*  the location stack                  */
  YY_InputCParser_LTYPE *yyls = yylsa;
  YY_InputCParser_LTYPE *yylsp;

#define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YY_InputCParser_PURE
  int YY_InputCParser_CHAR;
  YY_InputCParser_STYPE YY_InputCParser_LVAL;
  int YY_InputCParser_NERRS;
#ifdef YY_InputCParser_LSP_NEEDED
  YY_InputCParser_LTYPE YY_InputCParser_LLOC;
#endif
#endif

  YY_InputCParser_STYPE yyval;             /*  the variable used to return         */
				/*  semantic values from the action     */
				/*  routines                            */

  int yylen;
/* start loop, in which YYGOTO may be used. */
YYBEGINGOTO

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    fprintf(stderr, "Starting parse\n");
#endif
  yystate = 0;
  yyerrstatus = 0;
  YY_InputCParser_NERRS = 0;
  YY_InputCParser_CHAR = YYEMPTY;          /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YY_InputCParser_LSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
YYLABEL(yynewstate)

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YY_InputCParser_STYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YY_InputCParser_LSP_NEEDED
      YY_InputCParser_LTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
#ifdef YY_InputCParser_LSP_NEEDED
      /* This used to be a conditional around just the two extra args,
	 but that might be undefined if yyoverflow is a macro.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yyls1, size * sizeof (*yylsp),
		 &yystacksize);
#else
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
		 &yystacksize);
#endif

      yyss = yyss1; yyvs = yyvs1;
#ifdef YY_InputCParser_LSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  YY_InputCParser_ERROR("parser stack overflow");
	  __ALLOCA_return(2);
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
      yyss = (short *) __ALLOCA_alloca (yystacksize * sizeof (*yyssp));
      __yy_bcopy ((char *)yyss1, (char *)yyss, size * sizeof (*yyssp));
      __ALLOCA_free(yyss1,yyssa);
      yyvs = (YY_InputCParser_STYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yyvsp));
      __yy_bcopy ((char *)yyvs1, (char *)yyvs, size * sizeof (*yyvsp));
      __ALLOCA_free(yyvs1,yyvsa);
#ifdef YY_InputCParser_LSP_NEEDED
      yyls = (YY_InputCParser_LTYPE *) __ALLOCA_alloca (yystacksize * sizeof (*yylsp));
      __yy_bcopy ((char *)yyls1, (char *)yyls, size * sizeof (*yylsp));
      __ALLOCA_free(yyls1,yylsa);
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YY_InputCParser_LSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YY_InputCParser_DEBUG != 0
      if (YY_InputCParser_DEBUG_FLAG)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

  YYGOTO(yybackup);
YYLABEL(yybackup)

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* YYLABEL(yyresume) */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yydefault);

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (YY_InputCParser_CHAR == YYEMPTY)
    {
#if YY_InputCParser_DEBUG != 0
      if (YY_InputCParser_DEBUG_FLAG)
	fprintf(stderr, "Reading a token: ");
#endif
      YY_InputCParser_CHAR = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (YY_InputCParser_CHAR <= 0)           /* This means end of input. */
    {
      yychar1 = 0;
      YY_InputCParser_CHAR = YYEOF;                /* Don't call YYLEX any more */

#if YY_InputCParser_DEBUG != 0
      if (YY_InputCParser_DEBUG_FLAG)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(YY_InputCParser_CHAR);

#if YY_InputCParser_DEBUG != 0
      if (YY_InputCParser_DEBUG_FLAG)
	{
	  fprintf (stderr, "Next token is %d (%s", YY_InputCParser_CHAR, yytname[yychar1]);
	  /* Give the individual parser a way to print the precise meaning
	     of a token, for further debugging info.  */
#ifdef YYPRINT
	  YYPRINT (stderr, YY_InputCParser_CHAR, YY_InputCParser_LVAL);
#endif
	  fprintf (stderr, ")\n");
	}
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    YYGOTO(yydefault);

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrlab);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrlab);

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    fprintf(stderr, "Shifting token %d (%s), ", YY_InputCParser_CHAR, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (YY_InputCParser_CHAR != YYEOF)
    YY_InputCParser_CHAR = YYEMPTY;

  *++yyvsp = YY_InputCParser_LVAL;
#ifdef YY_InputCParser_LSP_NEEDED
  *++yylsp = YY_InputCParser_LLOC;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  YYGOTO(yynewstate);

/* Do the default action for the current state.  */
YYLABEL(yydefault)

  yyn = yydefact[yystate];
  if (yyn == 0)
    YYGOTO(yyerrlab);

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
YYLABEL(yyreduce)
  yylen = yyr2[yyn];
  if (yylen > 0)
    yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symbols being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 2:
{   
                            yyvsp[0].statement_handler->Execute(mpKernel); 
                            delete yyvsp[0].statement_handler;
                        ;
    break;}
case 3:
{yyval.statement_handler = yyvsp[-1].statement_handler;;
    break;}
case 4:
{yyval.statement_handler = yyvsp[0].block_statement_handler;;
    break;}
case 5:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 6:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 7:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 8:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 9:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 10:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 11:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 12:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 13:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 14:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 15:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 16:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 17:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 18:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 19:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 20:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 21:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 22:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 23:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 24:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 25:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 26:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 27:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 28:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 29:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 30:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 31:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 32:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 33:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 34:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 35:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 36:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 37:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 38:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 39:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 40:
{yyval.statement_handler = yyvsp[0].statement_handler;;
    break;}
case 41:
{
		    yyval.block_statement_handler = yyvsp[-1].block_statement_handler;
		;
    break;}
case 42:
{
			            yyval.block_statement_handler = new BlockStatement();
			            yyval.block_statement_handler->AddStatement(yyvsp[0].statement_handler);
			        ;
    break;}
case 43:
{
			            yyvsp[-1].block_statement_handler->AddStatement(yyvsp[0].statement_handler);
			            yyval.block_statement_handler = yyvsp[-1].block_statement_handler;
			        ;
    break;}
case 44:
{
                      node_Id = yyvsp[-9].integer_value;
                      node_X  = yyvsp[-5].double_value;
                      node_Y  = yyvsp[-3].double_value;
                      node_Z  = yyvsp[-1].double_value;

                      yyval.statement_handler = new GenerateNodeStatement();
                  ;
    break;}
case 45:
{
			            yyval.statement_handler = new GeneratePropertiesStatement<double>(yyvsp[-4].integer_value, *yyvsp[-2].double_variable, yyvsp[0].double_value);
			        ;
    break;}
case 46:
{
			            yyval.statement_handler = new GeneratePropertiesStatement<Kratos::Matrix<double> >(yyvsp[-4].integer_value, *yyvsp[-2].matrix_variable, *yyvsp[0].matrix_double_value);
                        delete yyvsp[0].matrix_double_value;
                    ;
    break;}
case 47:
{
			            yyval.statement_handler = new GeneratePropertiesStatement<Kratos::Vector<double> >(yyvsp[-4].integer_value, *yyvsp[-2].vector_double_variable, *yyvsp[0].vector_double_value);
                        delete yyvsp[0].vector_double_value;
                    ;
    break;}
case 48:
{
                         yyval.statement_handler = new SetElementsGroup(yyvsp[0].string_value);
                     ;
    break;}
case 49:
{
                       yyval.statement_handler = new SetProblemType(yyvsp[0].string_value);
                   ;
    break;}
case 50:
{
                            yyval.statement_handler = new SetFrequency(yyvsp[0].double_value);
                        ;
    break;}
case 51:
{
                          yyval.statement_handler = new SweepFrequency(yyvsp[-5].double_value, yyvsp[-3].double_value, yyvsp[-1].double_value);
                      ;
    break;}
case 52:
{
                                stNodesId[0]  = yyvsp[-9].integer_value;
                                stNodesId[1]  = yyvsp[-7].integer_value;
                                stNodesId[2]  = yyvsp[-5].integer_value;
                                stNodesId[3]  = yyvsp[-3].integer_value;
                                sPropertiesId = yyvsp[-1].integer_value;
                              
                                yyval.statement_handler = new GenerateVolumeElementStatement();
                            ;
    break;}
case 53:
{
                                stNodesId[0]  = yyvsp[-9].integer_value;
                                stNodesId[1]  = yyvsp[-7].integer_value;
                                stNodesId[2]  = yyvsp[-5].integer_value;
                                stNodesId[3]  = yyvsp[-3].integer_value;
                                sPropertiesId = yyvsp[-1].integer_value;
                              
                                yyval.statement_handler = new GenerateSourceElementStatement();
                            ;
    break;}
case 54:
{
                                 stNodesId[0]  = yyvsp[-7].integer_value;
                                 stNodesId[1]  = yyvsp[-5].integer_value;
                                 stNodesId[2]  = yyvsp[-3].integer_value;
                                 sPropertiesId = yyvsp[-1].integer_value;
                               
                                 yyval.statement_handler = new GenerateSrfIntgElementStatement();
                             ;
    break;}
case 55:
{
                                 stNodesId[0]  = yyvsp[-9].integer_value;
                                 stNodesId[1]  = yyvsp[-7].integer_value;
                                 stNodesId[2]  = yyvsp[-5].integer_value;
                                 stNodesId[3]  = yyvsp[-3].integer_value;
                                 sPropertiesId = yyvsp[-1].integer_value;
                               
                                 yyval.statement_handler = new GenerateVolIntgElementStatement();
                             ;
    break;}
case 56:
{
                             stNodesId[0]  = yyvsp[-9].integer_value;
                             stNodesId[1]  = yyvsp[-7].integer_value;
                             stNodesId[2]  = yyvsp[-5].integer_value;
                             stNodesId[3]  = yyvsp[-3].integer_value;
                             sPropertiesId = yyvsp[-1].integer_value;
                           
                             yyval.statement_handler = new GeneratePBCElementStatement();
                         ;
    break;}
case 57:
{
                                  stNodesId[0] = yyvsp[-5].integer_value;
                                  stNodesId[1] = yyvsp[-3].integer_value;
                                  stNodesId[2] = yyvsp[-1].integer_value;
                                
                                  yyval.statement_handler = new GenerateFarFieldElementStatement();
                              ;
    break;}
case 58:
{
                                  stNodesId[0]  = yyvsp[-7].integer_value;
                                  stNodesId[1]  = yyvsp[-5].integer_value;
                                  stNodesId[2]  = yyvsp[-3].integer_value;
                                  sPropertiesId = yyvsp[-1].integer_value;
                                
                                  yyval.statement_handler = new GenerateGenRobinElementStatement();
                              ;
    break;}
case 59:
{
                                        stNodesId[0]  = yyvsp[-7].integer_value;
                                        stNodesId[1]  = yyvsp[-5].integer_value;
                                        stNodesId[2]  = yyvsp[-3].integer_value;
                                        sPropertiesId = yyvsp[-1].integer_value;
                                      
                                        yyval.statement_handler = new GeneratePlasmaFarFieldElementStatement();
                                    ;
    break;}
case 60:
{
                                    stNodesId[0]  = yyvsp[-7].integer_value;
                                    stNodesId[1]  = yyvsp[-5].integer_value;
                                    stNodesId[2]  = yyvsp[-3].integer_value;
                                    sPropertiesId = yyvsp[-1].integer_value;
                                
                                    yyval.statement_handler = new GenerateRWPortTE10ElementStatement();
                                ;
    break;}
case 61:
{
                                     stNodesId[0]  = yyvsp[-7].integer_value;
                                     stNodesId[1]  = yyvsp[-5].integer_value;
                                     stNodesId[2]  = yyvsp[-3].integer_value;
                                     sPropertiesId = yyvsp[-1].integer_value;
                                 
                                     yyval.statement_handler = new GenerateCoaxPortTEMElementStatement();
                                 ;
    break;}
case 62:
{
                                    stNodesId[0]  = yyvsp[-7].integer_value;
                                    stNodesId[1]  = yyvsp[-5].integer_value;
                                    stNodesId[2]  = yyvsp[-3].integer_value;
                                    sPropertiesId = yyvsp[-1].integer_value;
                                 
                                    yyval.statement_handler = new GenerateProjectionCoaxialTEMStatement();
                                ;
    break;}
case 63:
{
                                   stNodesId[0]  = yyvsp[-7].integer_value;
                                   stNodesId[1]  = yyvsp[-5].integer_value;
                                   stNodesId[2]  = yyvsp[-3].integer_value;
                                   sPropertiesId = yyvsp[-1].integer_value;
                                
                                   yyval.statement_handler = new GenerateProjectionRWPortTE10Statement();
                               ;
    break;}
case 64:
{
                                sNodesId = *yyvsp[-1].vector_integer_value;
                                yyval.statement_handler = new GeneratePECElementStatement();
			                    delete yyvsp[-1].vector_integer_value; 
                            ;
    break;}
case 65:
{
                                sNodesId = *yyvsp[-1].vector_integer_value;
                                yyval.statement_handler = new GeneratePMCElementStatement();
			                    delete yyvsp[-1].vector_integer_value; 
                            ;
    break;}
case 66:
{
                                sNodesId = *yyvsp[-1].vector_integer_value;
                                yyval.statement_handler = new GenerateTECElementStatement();
			                    delete yyvsp[-1].vector_integer_value; 
                            ;
    break;}
case 67:
{
                                sNodesId = *yyvsp[-1].vector_integer_value;
                                yyval.statement_handler = new GenerateContactStatement();
			                    delete yyvsp[-1].vector_integer_value; 
                            ;
    break;}
case 68:
{
			                    yyval.statement_handler = new FixDofStatement(yyvsp[-7].integer_value, *yyvsp[-5].double_variable, yyvsp[-1].double_value);
                            ;
    break;}
case 69:
{
			                    yyval.statement_handler = new FixCDofStatement(yyvsp[-7].integer_value, *yyvsp[-5].double_variable, *yyvsp[-1].vector_double_value);
                                delete yyvsp[-1].vector_double_value;
                            ;
    break;}
case 70:
{
			                    yyval.statement_handler = new AssigningNodalVariableStatement<double>(yyvsp[-4].integer_value, *yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
                            ;
    break;}
case 71:
{
			                    yyval.statement_handler = new GenerateLinearSolverStatement(yyvsp[-9].string_value, yyvsp[-7].string_value, yyvsp[-5].integer_value, yyvsp[-3].integer_value, yyvsp[-1].double_value);
			                ;
    break;}
case 72:
{
			                    yyval.statement_handler = new ForStatement(yyvsp[-6].statement_handler, yyvsp[-4].statement_handler, yyvsp[-2].statement_handler, yyvsp[0].statement_handler);
                            ;
    break;}
case 73:
{
			                    yyval.statement_handler = new SolveStatement(yyvsp[-5].string_value, yyvsp[-1].string_value);
			                ;
    break;}
case 74:
{
			                    yyval.statement_handler = new PrintStatement<double>(*yyvsp[-1].double_variable);
			                ;
    break;}
case 75:
{
			                    yyval.statement_handler = new PrintStatement<Kratos::Vector<double> >(*yyvsp[-1].vector_double_variable);
			                ;
    break;}
case 76:
{
			                    yyval.statement_handler = new PrintStatement<Kratos::Matrix<double> >(*yyvsp[-1].matrix_variable);
			                ;
    break;}
case 77:
{
			                    yyval.statement_handler = new PrintOnNodesStatement<double>(*yyvsp[-1].double_variable);
			                ;
    break;}
case 78:
{
			                    yyval.statement_handler = new PrintOnNodesStatement<Kratos::Vector<double> >(*yyvsp[-1].vector_double_variable);
			                ;
    break;}
case 79:
{
			                    yyval.statement_handler = new PrintOnNodesStatement<Kratos::Matrix<double> >(*yyvsp[-1].matrix_variable);
			                ;
    break;}
case 80:
{
			                    yyval.statement_handler = new CalculateStatement<double>(*yyvsp[-1].double_variable);
			                ;
    break;}
case 81:
{
			                    yyval.statement_handler = new CalculateStatement<Kratos::Vector<double> >(*yyvsp[-1].vector_double_variable);
			                ;
    break;}
case 82:
{
			                    yyval.statement_handler = new CalculateStatement<Kratos::Matrix<double> >(*yyvsp[-1].matrix_variable);
			                ;
    break;}
case 83:
{
			                    yyval.statement_handler = new CalculateNodalStatement<double>(*yyvsp[-1].double_variable);
			                ;
    break;}
case 84:
{
			                    yyval.statement_handler = new CalculateNodalStatement<Kratos::Vector<double> >(*yyvsp[-1].vector_double_variable);
			                ;
    break;}
case 85:
{
			                    yyval.statement_handler = new CalculateNodalStatement<Kratos::Matrix<double> >(*yyvsp[-1].matrix_variable);
			                ;
    break;}
case 86:
{
			                    yyval.statement_handler = new CreateSolutionStepStatement();
			                ;
    break;}
case 87:
{
			                    yyval.statement_handler = new CreateTimeStepStatement();
			                ;
    break;}
case 88:
{ 
                                yyval.integer_value = yyvsp[0].integer_value; 
                            ;
    break;}
case 89:
{ yyval.integer_value = yyvsp[0].integer_value; ;
    break;}
case 90:
{
			                    yyval.statement_handler = new AssigningVariableStatement<double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 91:
{
			                    yyval.statement_handler = new LogicalStatement<std::less<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 92:
{
			                    yyval.statement_handler = new LogicalStatement<std::greater<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 93:
{
			                    yyval.statement_handler = new LogicalStatement<std::less_equal<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 94:
{
			                    yyval.statement_handler = new LogicalStatement<std::greater_equal<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 95:
{
			                    yyval.statement_handler = new LogicalStatement<std::equal_to<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 96:
{
			                    yyval.statement_handler = new LogicalStatement<std::not_equal_to<double>, double>(*yyvsp[-2].double_variable, yyvsp[0].double_value_statement_handler);
			                ;
    break;}
case 97:
{
			                    yyval.matrix_double_value = new Kratos::Matrix<double>(*yyvsp[-1].vector_double_value);
                                delete yyvsp[-1].vector_double_value;
			                ;
    break;}
case 98:
{
			                    yyval.matrix_double_value = new Kratos::Matrix<double>(*yyvsp[-1].vector_vector_double_value);
                                delete yyvsp[-1].vector_vector_double_value;
			                ;
    break;}
case 99:
{
			                    yyval.vector_vector_double_value = new Kratos::Vector<Kratos::Vector<double> >();
                                yyval.vector_vector_double_value->push_back(*yyvsp[-2].vector_double_value);
                                yyval.vector_vector_double_value->push_back(*yyvsp[0].vector_double_value);

                                delete yyvsp[-2].vector_double_value;
                                delete yyvsp[0].vector_double_value;
			                ;
    break;}
case 100:
{
			                    yyval.vector_vector_double_value = yyvsp[-2].vector_vector_double_value;
                                yyval.vector_vector_double_value->push_back(*yyvsp[0].vector_double_value);
                                delete yyvsp[0].vector_double_value;
			                ;
    break;}
case 101:
{
			                    yyval.vector_double_value = new Kratos::Vector<double>();
                                yyval.vector_double_value->push_back(yyvsp[-1].double_value);
			                ;
    break;}
case 102:
{ yyval.vector_double_value = yyvsp[-1].vector_double_value; ;
    break;}
case 103:
{
			                    yyval.vector_integer_value = new Kratos::Vector<int>();
                                yyval.vector_integer_value->push_back(yyvsp[0].integer_value);
			                ;
    break;}
case 104:
{ yyval.vector_integer_value = yyvsp[-1].vector_integer_value; ;
    break;}
case 105:
{ yyval.integer_value = yyvsp[-1].integer_value; ;
    break;}
case 106:
{
			                    yyval.vector_double_value = new Kratos::Vector<double>();
                                yyval.vector_double_value->push_back(yyvsp[-2].double_value);
                                yyval.vector_double_value->push_back(yyvsp[0].double_value);
			                ;
    break;}
case 107:
{
			                    yyval.vector_double_value = yyvsp[-2].vector_double_value;
                                yyval.vector_double_value->push_back(yyvsp[0].double_value);
			                ;
    break;}
case 108:
{
			                    yyval.vector_integer_value = new Kratos::Vector<int>();
                                yyval.vector_integer_value->push_back(yyvsp[-2].integer_value);
                                yyval.vector_integer_value->push_back(yyvsp[0].integer_value);
			                ;
    break;}
case 109:
{
			                    yyval.vector_integer_value = yyvsp[-2].vector_integer_value;
                                yyval.vector_integer_value->push_back(yyvsp[0].integer_value);
			                ;
    break;}
case 110:
{
			                    yyval.double_value_statement_handler = new ConstantValueStatement<double>(yyvsp[0].double_value);
			                ;
    break;}
case 111:
{
			                    yyval.double_value_statement_handler = new VariableStatement<double>(*yyvsp[0].double_variable);
			                ;
    break;}
case 112:
{
			                    yyval.double_value_statement_handler = new BinaryVariableStatement<std::plus<double>, double>(yyvsp[-2].double_value_statement_handler, new VariableStatement<double>(*yyvsp[0].double_variable));
			                ;
    break;}
case 113:
{
			                    yyval.double_value_statement_handler = new BinaryVariableStatement<std::plus<double>, double>(yyvsp[-2].double_value_statement_handler, new ConstantValueStatement<double>(yyvsp[0].double_value));
			                ;
    break;}
case 114:
{yyval.double_value =  yyvsp[0].double_value;;
    break;}
case 115:
{yyval.double_value =  yyvsp[0].double_value;;
    break;}
case 116:
{yyval.double_value = -yyvsp[0].double_value;;
    break;}
case 117:
{yyval.double_value =  yyvsp[0].integer_value;;
    break;}
case 118:
{yyval.integer_value =  yyvsp[0].integer_value;;
    break;}
}
   /* the action file gets copied in in place of this dollarsign  */
  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YY_InputCParser_LSP_NEEDED
  yylsp -= yylen;
#endif

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YY_InputCParser_LSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = YY_InputCParser_LLOC.first_line;
      yylsp->first_column = YY_InputCParser_LLOC.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  YYGOTO(yynewstate);

YYLABEL(yyerrlab)   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++YY_InputCParser_NERRS;

#ifdef YY_InputCParser_ERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  /* Start X at -yyn if nec to avoid negative indexes in yycheck.  */
	  for (x = (yyn < 0 ? -yyn : 0);
	       x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) malloc(size + 15);
	  if (msg != 0)
	    {
	      strcpy(msg, "parse error");

	      if (count < 5)
		{
		  count = 0;
		  for (x = (yyn < 0 ? -yyn : 0);
		       x < (sizeof(yytname) / sizeof(char *)); x++)
		    if (yycheck[x + yyn] == x)
		      {
			strcat(msg, count == 0 ? ", expecting `" : " or `");
			strcat(msg, yytname[x]);
			strcat(msg, "'");
			count++;
		      }
		}
	      YY_InputCParser_ERROR(msg);
	      free(msg);
	    }
	  else
	    YY_InputCParser_ERROR ("parse error; also virtual memory exceeded");
	}
      else
#endif /* YY_InputCParser_ERROR_VERBOSE */
	YY_InputCParser_ERROR("parse error");
    }

  YYGOTO(yyerrlab1);
YYLABEL(yyerrlab1)   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (YY_InputCParser_CHAR == YYEOF)
	YYABORT;

#if YY_InputCParser_DEBUG != 0
      if (YY_InputCParser_DEBUG_FLAG)
	fprintf(stderr, "Discarding token %d (%s).\n", YY_InputCParser_CHAR, yytname[yychar1]);
#endif

      YY_InputCParser_CHAR = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;              /* Each real token shifted decrements this */

  YYGOTO(yyerrhandle);

YYLABEL(yyerrdefault)  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) YYGOTO(yydefault);
#endif

YYLABEL(yyerrpop)   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YY_InputCParser_LSP_NEEDED
  yylsp--;
#endif

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

YYLABEL(yyerrhandle)

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    YYGOTO(yyerrdefault);

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    YYGOTO(yyerrdefault);

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	YYGOTO(yyerrpop);
      yyn = -yyn;
      YYGOTO(yyreduce);
    }
  else if (yyn == 0)
    YYGOTO(yyerrpop);

  if (yyn == YYFINAL)
    YYACCEPT;

#if YY_InputCParser_DEBUG != 0
  if (YY_InputCParser_DEBUG_FLAG)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = YY_InputCParser_LVAL;
#ifdef YY_InputCParser_LSP_NEEDED
  *++yylsp = YY_InputCParser_LLOC;
#endif

  yystate = yyn;
  YYGOTO(yynewstate);
/* end loop, in which YYGOTO may be used. */
  YYENDGOTO
}

/* END */

using Kratos::Exception;
void InputCParser::yyerror( char *s) 
{
    Kratos::String buffer;
    buffer << "Pars error in line no " << mInputCScanner.rNumberOfLines() << ", lass token was : " << mInputCScanner.GetYYText();
    KRATOS_ERROR(std::runtime_error, "Reading Input file", buffer, "")
}
