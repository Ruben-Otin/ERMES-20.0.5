#ifndef YY_InputCParser_h_included
#define YY_InputCParser_h_included
/* before anything */
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
#endif
#include <stdio.h>

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
#define YY_InputCParser_CONSTRUCTOR_PARAM  Kratos::Kernel* pKernel, std::istream* pNewInput, std::ostream* pNewOutput
#define YY_InputCParser_CONSTRUCTOR_INIT  : mInputCScanner(pNewInput, pNewOutput), mpKernel(pKernel), mpInput(pNewInput)
#define YY_InputCParser_YY_InputGid_CONSTRUCTOR_CODE  
#define YY_InputCParser_LEX_BODY  { return mInputCScanner.yylex(yylval); }
#define YY_InputCParser_MEMBERS   Kratos::InputCScanner mInputCScanner; Kratos::Kernel* mpKernel; std::istream* mpInput; Kratos::String mBlockName;
 /* %{ and %header{ and %union, during decl */
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
/* WARNING obsolete !!! user defined YYLTYPE not reported into generated header */
/* use %define LTYPE */
#endif
#endif
#ifdef YYSTYPE
#ifndef YY_InputCParser_STYPE 
#define YY_InputCParser_STYPE YYSTYPE
/* WARNING obsolete !!! user defined YYSTYPE not reported into generated header */
/* use %define STYPE */
#endif
#endif
#ifdef YYDEBUG
#ifndef YY_InputCParser_DEBUG
#define  YY_InputCParser_DEBUG YYDEBUG
/* WARNING obsolete !!! user defined YYDEBUG not reported into generated header */
/* use %define DEBUG */
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
/* YY_InputCParser_PURE */
#endif
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

/* TOKEN C */
#ifndef YY_USE_CLASS

#ifndef YY_InputCParser_PURE
extern YY_InputCParser_STYPE YY_InputCParser_LVAL;
#endif

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

 /* #defines token */
/* after #define tokens, before const tokens S5*/
#else
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
 int YY_InputCParser_PARSE(YY_InputCParser_PARSE_PARAM);
 virtual void YY_InputCParser_ERROR(char *msg) YY_InputCParser_ERROR_BODY;
#ifdef YY_InputCParser_PURE
#ifdef YY_InputCParser_LSP_NEEDED
 virtual int  YY_InputCParser_LEX(YY_InputCParser_STYPE *YY_InputCParser_LVAL,YY_InputCParser_LTYPE *YY_InputCParser_LLOC) YY_InputCParser_LEX_BODY;
#else
 virtual int  YY_InputCParser_LEX(YY_InputCParser_STYPE *YY_InputCParser_LVAL) YY_InputCParser_LEX_BODY;
#endif
#else
 virtual int YY_InputCParser_LEX() YY_InputCParser_LEX_BODY;
 YY_InputCParser_STYPE YY_InputCParser_LVAL;
#ifdef YY_InputCParser_LSP_NEEDED
 YY_InputCParser_LTYPE YY_InputCParser_LLOC;
#endif
 int YY_InputCParser_NERRS;
 int YY_InputCParser_CHAR;
#endif
#if YY_InputCParser_DEBUG != 0
public:
 int YY_InputCParser_DEBUG_FLAG;	/*  nonzero means print parse trace	*/
#endif
public:
 YY_InputCParser_CLASS(YY_InputCParser_CONSTRUCTOR_PARAM);
public:
 YY_InputCParser_MEMBERS 
};
/* other declare folow */
#endif


#if YY_InputCParser_COMPATIBILITY != 0
/* backward compatibility */
#ifndef YYSTYPE
#define YYSTYPE YY_InputCParser_STYPE
#endif

#ifndef YYLTYPE
#define YYLTYPE YY_InputCParser_LTYPE
#endif
#ifndef YYDEBUG
#ifdef YY_InputCParser_DEBUG 
#define YYDEBUG YY_InputCParser_DEBUG
#endif
#endif

#endif
/* END */
#endif
