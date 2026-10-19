%name InputCParser

%header{
#include "../includes/variables.h"
  class Statement;
  class BlockStatement;
  template<class TDataType> class ValueStatement;
%}

/* This union must be before include input_c_scanner.h, cause yylex needs yystype */
%union {
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
};

%header{
#include "../includes/kernel.h"
#include "../includes/input_c_scanner.h"
#include "../includes/element_sources.h"

using Kratos::Exception;
%}
%{
  Kratos::String block_name;

  class Statement
  {
    public:
	
      Statement(){}
	  
      virtual void Execute(Kratos::Kernel* pKernel){}
	  
      virtual int Value(Kratos::Kernel* pKernel){return 0;}
  };

  class BlockStatement : public Statement
  {
    public:

      BlockStatement(){}

      void AddStatement(Statement* pNewStatement)
      {
          mStatements.push_back(boost::shared_ptr<Statement>(pNewStatement));
      }

      void Execute(Kratos::Kernel* pKernel)
      {
          for(Kratos::Vector<boost::shared_ptr<Statement> >::iterator i = mStatements.begin() ; i != mStatements.end() ; i++) (*i)->Execute(pKernel);
      }

      Kratos::Vector<boost::shared_ptr<Statement> > mStatements;
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

%}

%define CONSTRUCTOR_PARAM Kratos::Kernel* pKernel, std::istream* pNewInput, std::ostream* pNewOutput
%define CONSTRUCTOR_INIT : mInputCScanner(pNewInput, pNewOutput), mpKernel(pKernel), mpInput(pNewInput)
%define YY_InputGid_CONSTRUCTOR_CODE 
%define LEX_BODY { return mInputCScanner.yylex(yylval); }
%define MEMBERS  Kratos::InputCScanner mInputCScanner; Kratos::Kernel* mpKernel; std::istream* mpInput; Kratos::String mBlockName;

%token <integer_value> INTEGER_TOKEN
%token <double_value> DOUBLE_TOKEN
%token <string_value> WORD_TOKEN
%token <double_variable> DOUBLE_VARIABLE_TOKEN
%token <vector_double_variable> VECTOR_DOUBLE_VARIABLE_TOKEN
%token <matrix_variable> MATRIX_VARIABLE_TOKEN
%token OPEN_BRACKET_TOKEN
%token CLOSE_BRACKET_TOKEN
%token NODES_TOKEN
%token NODE_TOKEN
%token NO_TOKEN
%token PROPERTIES_TOKEN
%token ELEMENTS_GROUP_TOKEN
%token PROBLEM_TYPE_TOKEN
%token PROBLEM_FREQUENCY_TOKEN
%token SWEEP_FREQUENCY_TOKEN
%token VE_TOKEN
%token FF_TOKEN
%token GRC_TOKEN
%token PFF_TOKEN
%token COP_TOKEN
%token RWP_TOKEN
%token JE_TOKEN
%token PVIE_TOKEN
%token PSIE_TOKEN
%token PCOP_TOKEN
%token PRWP_TOKEN
%token PES_TOKEN
%token PEP_TOKEN
%token PEC_TOKEN
%token PBC_TOKEN
%token PMC_TOKEN
%token TEC_TOKEN
%token CE_TOKEN
%token ELEMENT_TOKEN
%token BEGIN_TOKEN
%token END_TOKEN
%token FIX_TOKEN
%token FIXC_TOKEN
%token SOURCES_TOKEN
%token FOR_TOKEN
%token SOLVE_TOKEN
%token PRINT_TOKEN
%token PRINT_ON_NODES_TOKEN
%token EXECUTE_TOKEN
%token TRANSIENT_TOKEN
%token ALPHA_TOKEN
%token CALCULATE_TOKEN
%token CALCULATE_NODAL_TOKEN
%token LINEAR_SOLVER_TOKEN
%token SOLVING_STRATEGY_TOKEN
%token EQUAL_TOKEN
%token NOT_EQUAL_TOKEN
%token LESS_EQUAL_TOKEN
%token GREATER_EQUAL_TOKEN
%token CREATE_SOLUTION_STEP_TOKEN
%token CREATE_TIME_STEP_TOKEN

%type <statement_handler> statement
%type <statement_handler> expresion
%type <statement_handler> node_generating
%type <statement_handler> properties_adding
%type <statement_handler> set_elements_group
%type <statement_handler> set_problem_type
%type <statement_handler> set_problem_frequency
%type <statement_handler> set_sweep_frequency
%type <statement_handler> volume_element_generating
%type <statement_handler> source_element_generating
%type <statement_handler> volIntg_element_generating
%type <statement_handler> PBC_element_generating
%type <statement_handler> srfIntg_element_generating
%type <statement_handler> farfield_element_generating
%type <statement_handler> genRobin_element_generating
%type <statement_handler> plasmafarfield_element_generating
%type <statement_handler> coaxportTEM_element_generating
%type <statement_handler> rwportTE10_element_generating
%type <statement_handler> projection_coaxTEM_generating
%type <statement_handler> projection_rwTE10_generating
%type <statement_handler> pec_element_generating
%type <statement_handler> pmc_element_generating
%type <statement_handler> tec_element_generating
%type <statement_handler> contact_generation
%type <statement_handler> nodes_variables_fixing
%type <statement_handler> nodes_cvariables_fixing
%type <statement_handler> nodes_variables_setting
%type <statement_handler> linear_solver_generating
%type <statement_handler> solve_process
%type <statement_handler> print_process
%type <statement_handler> print_on_nodes_process
%type <statement_handler> calculate_process
%type <statement_handler> calculate_nodal_process
%type <statement_handler> step_creating
%type <statement_handler> time_step_creating
%type <statement_handler> for_loop
%type <statement_handler> assignment_expresion
%type <statement_handler> logical_expresion
%type <block_statement_handler> block_generating
%type <block_statement_handler> block
%type <double_value_statement_handler> double_variable_expresion
%type <integer_value> integer_expresion
%type <integer_value> integer_index
%type <integer_value> nodes_array
%type <integer_value> properties_array
%type <double_value> double_expresion
%type <vector_integer_value> integer_sequence
%type <vector_integer_value> vector_integer
%type <vector_double_value> vector_double
%type <matrix_double_value> matrix
%type <vector_double_value> double_sequence
%type <vector_vector_double_value> vector_double_sequence

%%
statements              :   
                        |   statements statement{   
                            $2->Execute(mpKernel); 
                            delete $2;
                        };

statement               :   expresion ';' {$$ = $1;}
                        |   block         {$$ = $1;}
                        |   for_loop      {$$ = $1;}
                        ;

expresion               :   assignment_expresion             {$$ = $1;}
                        |   logical_expresion                {$$ = $1;}
                        |   node_generating                  {$$ = $1;}
                        |   properties_adding                {$$ = $1;}
                        |   set_elements_group               {$$ = $1;}
                        |   set_problem_type                 {$$ = $1;}
                        |   set_problem_frequency            {$$ = $1;}
                        |   set_sweep_frequency              {$$ = $1;}
                        |   volume_element_generating        {$$ = $1;}
                        |   source_element_generating        {$$ = $1;}
                        |   volIntg_element_generating       {$$ = $1;}
                        |   PBC_element_generating           {$$ = $1;}
                        |   srfIntg_element_generating       {$$ = $1;}
                        |   farfield_element_generating      {$$ = $1;}
                        |   genRobin_element_generating      {$$ = $1;}
                        |   plasmafarfield_element_generating{$$ = $1;}
                        |   coaxportTEM_element_generating   {$$ = $1;}
                        |   rwportTE10_element_generating    {$$ = $1;}
                        |   projection_coaxTEM_generating    {$$ = $1;}
                        |   projection_rwTE10_generating     {$$ = $1;}
                        |   pec_element_generating           {$$ = $1;}
                        |   pmc_element_generating           {$$ = $1;}
                        |   tec_element_generating           {$$ = $1;}
                        |   contact_generation               {$$ = $1;}    
                        |   nodes_variables_fixing           {$$ = $1;}
                        |   nodes_cvariables_fixing          {$$ = $1;}
                        |   nodes_variables_setting          {$$ = $1;}
                        |   solve_process                    {$$ = $1;}
                        |   print_process                    {$$ = $1;}
						|   print_on_nodes_process           {$$ = $1;}
                        |   calculate_process                {$$ = $1;}
                        |   calculate_nodal_process          {$$ = $1;}
                        |   linear_solver_generating         {$$ = $1;}
                        |   step_creating                    {$$ = $1;} 
                        |   time_step_creating               {$$ = $1;} 
                        ;

block : block_generating '}' 
        {
		    $$ = $1;
		}

block_generating : '{' statement
                    {
			            $$ = new BlockStatement();
			            $$->AddStatement($2);
			        }
                    |   block_generating statement
                    {
			            $1->AddStatement($2);
			            $$ = $1;
			        }

node_generating : nodes_array '=' WORD_TOKEN '(' 
                  double_expresion ',' 
                  double_expresion ','
                  double_expresion ')' 
                  {
                      node_Id = $1;
                      node_X  = $5;
                      node_Y  = $7;
                      node_Z  = $9;

                      $$ = new GenerateNodeStatement();
                  };

properties_adding : properties_array '.' DOUBLE_VARIABLE_TOKEN '='
                    double_expresion 
                    {
			            $$ = new GeneratePropertiesStatement<double>($1, *$3, $5);
			        }
                    |   properties_array '.' MATRIX_VARIABLE_TOKEN '=' 
					    matrix  
                    {
			            $$ = new GeneratePropertiesStatement<Kratos::Matrix<double> >($1, *$3, *$5);
                        delete $5;
                    }
                    |   properties_array '.' VECTOR_DOUBLE_VARIABLE_TOKEN '=' 
					    vector_double
                    {
			            $$ = new GeneratePropertiesStatement<Kratos::Vector<double> >($1, *$3, *$5);
                        delete $5;
                    }

set_elements_group : ELEMENTS_GROUP_TOKEN '=' WORD_TOKEN
                     {
                         $$ = new SetElementsGroup($3);
                     }

set_problem_type : PROBLEM_TYPE_TOKEN '=' WORD_TOKEN
                   {
                       $$ = new SetProblemType($3);
                   }

set_problem_frequency : PROBLEM_FREQUENCY_TOKEN '='  
                        double_expresion 
                        {
                            $$ = new SetFrequency($3);
                        }

set_sweep_frequency : SWEEP_FREQUENCY_TOKEN '=' WORD_TOKEN '(' 
                      double_expresion ',' 
                      double_expresion ','
                      double_expresion ')'  
                      {
                          $$ = new SweepFrequency($5, $7, $9);
                      }

volume_element_generating : VE_TOKEN '('
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ')' 
                            {
                                stNodesId[0]  = $3;
                                stNodesId[1]  = $5;
                                stNodesId[2]  = $7;
                                stNodesId[3]  = $9;
                                sPropertiesId = $11;
                              
                                $$ = new GenerateVolumeElementStatement();
                            }

source_element_generating : JE_TOKEN '('
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ','
                            integer_expresion ')' 
                            {
                                stNodesId[0]  = $3;
                                stNodesId[1]  = $5;
                                stNodesId[2]  = $7;
                                stNodesId[3]  = $9;
                                sPropertiesId = $11;
                              
                                $$ = new GenerateSourceElementStatement();
                            }

srfIntg_element_generating : PSIE_TOKEN '('
                             integer_expresion ','
                             integer_expresion ','
                             integer_expresion ','
                             integer_expresion ')' 
                             {
                                 stNodesId[0]  = $3;
                                 stNodesId[1]  = $5;
                                 stNodesId[2]  = $7;
                                 sPropertiesId = $9;
                               
                                 $$ = new GenerateSrfIntgElementStatement();
                             }

volIntg_element_generating : PVIE_TOKEN '('
                             integer_expresion ','
                             integer_expresion ','
                             integer_expresion ','
                             integer_expresion ','
                             integer_expresion ')' 
                             {
                                 stNodesId[0]  = $3;
                                 stNodesId[1]  = $5;
                                 stNodesId[2]  = $7;
                                 stNodesId[3]  = $9;
                                 sPropertiesId = $11;
                               
                                 $$ = new GenerateVolIntgElementStatement();
                             }
							
PBC_element_generating : PBC_TOKEN '('
                         integer_expresion ','
                         integer_expresion ','
                         integer_expresion ','
                         integer_expresion ','
                         integer_expresion ')' 
                         {
                             stNodesId[0]  = $3;
                             stNodesId[1]  = $5;
                             stNodesId[2]  = $7;
                             stNodesId[3]  = $9;
                             sPropertiesId = $11;
                           
                             $$ = new GeneratePBCElementStatement();
                         }

farfield_element_generating : FF_TOKEN '('
                              integer_expresion ','
                              integer_expresion ','
                              integer_expresion ')' 
                              {
                                  stNodesId[0] = $3;
                                  stNodesId[1] = $5;
                                  stNodesId[2] = $7;
                                
                                  $$ = new GenerateFarFieldElementStatement();
                              }

genRobin_element_generating : GRC_TOKEN '('
                              integer_expresion ','
                              integer_expresion ','
                              integer_expresion ','
                              integer_expresion ')' 
                              {
                                  stNodesId[0]  = $3;
                                  stNodesId[1]  = $5;
                                  stNodesId[2]  = $7;
                                  sPropertiesId = $9;
                                
                                  $$ = new GenerateGenRobinElementStatement();
                              }
							
plasmafarfield_element_generating : PFF_TOKEN '('
                                    integer_expresion ','
                                    integer_expresion ','
                                    integer_expresion ','
                                    integer_expresion ')' 
                                    {
                                        stNodesId[0]  = $3;
                                        stNodesId[1]  = $5;
                                        stNodesId[2]  = $7;
                                        sPropertiesId = $9;
                                      
                                        $$ = new GeneratePlasmaFarFieldElementStatement();
                                    }

rwportTE10_element_generating : RWP_TOKEN '('
                                integer_expresion ','
                                integer_expresion ','
                                integer_expresion ','
                                integer_expresion ')' 
                                {
                                    stNodesId[0]  = $3;
                                    stNodesId[1]  = $5;
                                    stNodesId[2]  = $7;
                                    sPropertiesId = $9;
                                
                                    $$ = new GenerateRWPortTE10ElementStatement();
                                }

coaxportTEM_element_generating : COP_TOKEN '('
                                 integer_expresion ','
                                 integer_expresion ','
                                 integer_expresion ','
                                 integer_expresion ')' 
                                 {
                                     stNodesId[0]  = $3;
                                     stNodesId[1]  = $5;
                                     stNodesId[2]  = $7;
                                     sPropertiesId = $9;
                                 
                                     $$ = new GenerateCoaxPortTEMElementStatement();
                                 }

projection_coaxTEM_generating : PCOP_TOKEN '('
                                integer_expresion ','
                                integer_expresion ','
                                integer_expresion ','
                                integer_expresion ')' 
                                {
                                    stNodesId[0]  = $3;
                                    stNodesId[1]  = $5;
                                    stNodesId[2]  = $7;
                                    sPropertiesId = $9;
                                 
                                    $$ = new GenerateProjectionCoaxialTEMStatement();
                                }

projection_rwTE10_generating : PRWP_TOKEN '('
                               integer_expresion ','
                               integer_expresion ','
                               integer_expresion ','
                               integer_expresion ')' 
                               {
                                   stNodesId[0]  = $3;
                                   stNodesId[1]  = $5;
                                   stNodesId[2]  = $7;
                                   sPropertiesId = $9;
                                
                                   $$ = new GenerateProjectionRWPortTE10Statement();
                               }

pec_element_generating  :   PEC_TOKEN '=' WORD_TOKEN '('
                            vector_integer')' 
                            {
                                sNodesId = *$5;
                                $$ = new GeneratePECElementStatement();
			                    delete $5; 
                            }
							
pmc_element_generating  :   PMC_TOKEN '=' WORD_TOKEN '('
                            vector_integer')' 
                            {
                                sNodesId = *$5;
                                $$ = new GeneratePMCElementStatement();
			                    delete $5; 
                            }

tec_element_generating:     TEC_TOKEN '=' WORD_TOKEN '('
                            vector_integer')' 
                            {
                                sNodesId = *$5;
                                $$ = new GenerateTECElementStatement();
			                    delete $5; 
                            }

contact_generation      :   CE_TOKEN '=' WORD_TOKEN '(' 
                            vector_integer ')' 
                            {
                                sNodesId = *$5;
                                $$ = new GenerateContactStatement();
			                    delete $5; 
                            }

nodes_variables_fixing  :   nodes_array '.' DOUBLE_VARIABLE_TOKEN '.' FIX_TOKEN '(' 
                            double_expresion ')'
                            {
			                    $$ = new FixDofStatement($1, *$3, $7);
                            }

nodes_cvariables_fixing :   nodes_array '.' DOUBLE_VARIABLE_TOKEN '.' FIXC_TOKEN '(' 
                            vector_double ')'
                            {
			                    $$ = new FixCDofStatement($1, *$3, *$7);
                                delete $7;
                            }

nodes_variables_setting :   nodes_array '.' DOUBLE_VARIABLE_TOKEN '='
                            double_variable_expresion 
                            {
			                    $$ = new AssigningNodalVariableStatement<double>($1, *$3, $5);
                            }

linear_solver_generating:   LINEAR_SOLVER_TOKEN WORD_TOKEN '=' WORD_TOKEN '(' 
                            integer_expresion ',' 
                            integer_expresion ','
                            double_expresion  ')'
                            {
			                    $$ = new GenerateLinearSolverStatement($2, $4, $6, $8, $10);
			                }

for_loop                :   FOR_TOKEN '('
                            assignment_expresion ';' 
                            logical_expresion    ';'
                            expresion            ')'
                            statement
                            {
			                    $$ = new ForStatement($3, $5, $7, $9);
                            }

solve_process           :   WORD_TOKEN '.' SOLVE_TOKEN '(' WORD_TOKEN ')'
                            {
			                    $$ = new SolveStatement($1, $5);
			                }
                            
print_process           :   PRINT_TOKEN '(' DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintStatement<double>(*$3);
			                }
                        |   PRINT_TOKEN '(' VECTOR_DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintStatement<Kratos::Vector<double> >(*$3);
			                }
                        |   PRINT_TOKEN '(' MATRIX_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintStatement<Kratos::Matrix<double> >(*$3);
			                }
							
print_on_nodes_process  :   PRINT_ON_NODES_TOKEN '(' DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintOnNodesStatement<double>(*$3);
			                }
                        |   PRINT_ON_NODES_TOKEN '(' VECTOR_DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintOnNodesStatement<Kratos::Vector<double> >(*$3);
			                }
                        |   PRINT_ON_NODES_TOKEN '(' MATRIX_VARIABLE_TOKEN ')'
                            {
			                    $$ = new PrintOnNodesStatement<Kratos::Matrix<double> >(*$3);
			                }

calculate_process       :   CALCULATE_TOKEN '(' DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateStatement<double>(*$3);
			                }
                        |   CALCULATE_TOKEN '(' VECTOR_DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateStatement<Kratos::Vector<double> >(*$3);
			                }
                        |   CALCULATE_TOKEN '(' MATRIX_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateStatement<Kratos::Matrix<double> >(*$3);
			                }

calculate_nodal_process :   CALCULATE_NODAL_TOKEN '(' DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateNodalStatement<double>(*$3);
			                }
                        |   CALCULATE_NODAL_TOKEN '(' VECTOR_DOUBLE_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateNodalStatement<Kratos::Vector<double> >(*$3);
			                }
                        |   CALCULATE_NODAL_TOKEN '(' MATRIX_VARIABLE_TOKEN ')'
                            {
			                    $$ = new CalculateNodalStatement<Kratos::Matrix<double> >(*$3);
			                }

step_creating           :   CREATE_SOLUTION_STEP_TOKEN  '(' ')'
                            {
			                    $$ = new CreateSolutionStepStatement();
			                };

time_step_creating      :   CREATE_TIME_STEP_TOKEN  '(' ')'
                            {
			                    $$ = new CreateTimeStepStatement();
			                };

properties_array        :   PROPERTIES_TOKEN integer_index 
                            { 
                                $$ = $2; 
                            }

nodes_array             :   NO_TOKEN integer_index { $$ = $2; };

assignment_expresion    :   DOUBLE_VARIABLE_TOKEN '=' double_variable_expresion 
                            {
			                    $$ = new AssigningVariableStatement<double>(*$1, $3);
			                };

logical_expresion       :   DOUBLE_VARIABLE_TOKEN '<' double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::less<double>, double>(*$1, $3);
			                }
                        |   DOUBLE_VARIABLE_TOKEN '>' double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::greater<double>, double>(*$1, $3);
			                }
                        |   DOUBLE_VARIABLE_TOKEN LESS_EQUAL_TOKEN double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::less_equal<double>, double>(*$1, $3);
			                }
                        |   DOUBLE_VARIABLE_TOKEN GREATER_EQUAL_TOKEN double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::greater_equal<double>, double>(*$1, $3);
			                }
                        |   DOUBLE_VARIABLE_TOKEN EQUAL_TOKEN double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::equal_to<double>, double>(*$1, $3);
			                }
                        |   DOUBLE_VARIABLE_TOKEN NOT_EQUAL_TOKEN double_variable_expresion
                            {
			                    $$ = new LogicalStatement<std::not_equal_to<double>, double>(*$1, $3);
			                };

matrix                  :   '[' vector_double']'
                            {
			                    $$ = new Kratos::Matrix<double>(*$2);
                                delete $2;
			                }
                        |   '[' vector_double_sequence ']'
                            {
			                    $$ = new Kratos::Matrix<double>(*$2);
                                delete $2;
			                };

vector_double_sequence  :   vector_double ',' vector_double
                            {
			                    $$ = new Kratos::Vector<Kratos::Vector<double> >();
                                $$->push_back(*$1);
                                $$->push_back(*$3);

                                delete $1;
                                delete $3;
			                }
                        |   vector_double_sequence ',' vector_double
                            {
			                    $$ = $1;
                                $$->push_back(*$3);
                                delete $3;
			                };

vector_double           :   '[' double_expresion ']'
                            {
			                    $$ = new Kratos::Vector<double>();
                                $$->push_back($2);
			                }
                        |   '[' double_sequence ']'  { $$ = $2; };

vector_integer          :   integer_index
                            {
			                    $$ = new Kratos::Vector<int>();
                                $$->push_back($1);
			                }
                        |   '[' integer_sequence ']'  { $$ = $2; };

integer_index           :   '[' integer_expresion ']' { $$ = $2; };

double_sequence         :   double_expresion ',' double_expresion
                            {
			                    $$ = new Kratos::Vector<double>();
                                $$->push_back($1);
                                $$->push_back($3);
			                }
                        |   double_sequence ',' double_expresion
                            {
			                    $$ = $1;
                                $$->push_back($3);
			                }; 

integer_sequence        :   integer_expresion ',' integer_expresion
                            {
			                    $$ = new Kratos::Vector<int>();
                                $$->push_back($1);
                                $$->push_back($3);
			                }
                        |   integer_sequence ',' integer_expresion
                            {
			                    $$ = $1;
                                $$->push_back($3);
			                };

double_variable_expresion : double_expresion
                            {
			                    $$ = new ConstantValueStatement<double>($1);
			                }
                        |   DOUBLE_VARIABLE_TOKEN
                            {
			                    $$ = new VariableStatement<double>(*$1);
			                }
                        |   double_variable_expresion '+' DOUBLE_VARIABLE_TOKEN
                            {
			                    $$ = new BinaryVariableStatement<std::plus<double>, double>($1, new VariableStatement<double>(*$3));
			                }
                        |   double_variable_expresion '+'  double_expresion
                            {
			                    $$ = new BinaryVariableStatement<std::plus<double>, double>($1, new ConstantValueStatement<double>($3));
			                }

double_expresion        :   DOUBLE_TOKEN         {$$ =  $1;}
                        |   '+' double_expresion {$$ =  $2;}
                        |   '-' double_expresion {$$ = -$2;}
                        |   integer_expresion    {$$ =  $1;};

integer_expresion       :   INTEGER_TOKEN        {$$ =  $1;};

%%
using Kratos::Exception;
void InputCParser::yyerror( char *s) 
{
    Kratos::String buffer;
    buffer << "Pars error in line no " << mInputCScanner.rNumberOfLines() << ", lass token was : " << mInputCScanner.GetYYText();
    KRATOS_ERROR(std::runtime_error, "Reading Input file", buffer, "")
}
