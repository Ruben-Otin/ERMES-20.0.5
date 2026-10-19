%option noyywrap

%{
#define YY_DECL int Kratos::InputCScanner::yylex(yy_InputCParser_stype &yylval)

#include <iostream>
#include <malloc.h>
#include "../includes/input_c_parser.h"
#include "../includes/input_c_scanner.h"

using namespace std;

int yyFlexLexer::yylex()    
{                           
    return 0;
}
  

%}

%x COMMENT

%%

\n                      { ++mNumberOfLines; }

"//".*$                 /* remove one line coments */

"/*"			   { BEGIN COMMENT;  /* Enter comment eating state */     }
<COMMENT>"/*"      { fprintf(stderr, "Warning! Comment within comment."); }
<COMMENT>\n		   { ++mNumberOfLines;                                    }
<COMMENT>[^*/\n]*	  ;
<COMMENT>"*"+[^*/\n]* ;
<COMMENT>"/"+[^*/\n]* ;
<COMMENT>"*"+"/"	{ BEGIN INITIAL; }

\[|\]|"("|")"|"{"|"}" { return yytext[0]; }

"="                   { return yytext[0]; }

","|";"|"."|":"       { return yytext[0]; }

"<"|">"               { return yytext[0]; }

"=="                  { return InputCParser::EQUAL_TOKEN;         }

"!="                  { return InputCParser::NOT_EQUAL_TOKEN;     }

"<="                  { return InputCParser::LESS_EQUAL_TOKEN;    }

">="                  { return InputCParser::GREATER_EQUAL_TOKEN; }

"+"|"-"               { return yytext[0]; }

[ \t]+ /* remove whitespaces */

[[:digit:]]+ {
          sscanf(yytext,"%d",&yylval.integer_value);
          return InputCParser::INTEGER_TOKEN;
         }

[[:digit:]]+"."[[:digit:]]* {
                     sscanf(yytext,"%lf",&yylval.double_value);
                     return InputCParser::DOUBLE_TOKEN;
                    }

NODES            { return InputCParser::NODES_TOKEN;             }

PROPERTIES       { return InputCParser::PROPERTIES_TOKEN;        }

Node             { return InputCParser::NODE_TOKEN;              }

No               { return InputCParser::NO_TOKEN;                }
 
Fix              { return InputCParser::FIX_TOKEN;               }

FixC             { return InputCParser::FIXC_TOKEN;              }

ElementsGroup    { return InputCParser::ELEMENTS_GROUP_TOKEN;    }

ProblemType      { return InputCParser::PROBLEM_TYPE_TOKEN;      }

ProblemFrequency { return InputCParser::PROBLEM_FREQUENCY_TOKEN; }

SweepFrequency   { return InputCParser::SWEEP_FREQUENCY_TOKEN;   }

VE               { return InputCParser::VE_TOKEN;   }    /* Volume elements */
																								 
RWP              { return InputCParser::RWP_TOKEN;  }    /* RW port elements */
													
COP              { return InputCParser::COP_TOKEN;  }    /* Coaxial port elements */
													
GRC              { return InputCParser::GRC_TOKEN;  }    /* Generic Robin elements */
													
PFF              { return InputCParser::PFF_TOKEN;  }    /* Plasma far field elements */
													
FF               { return InputCParser::FF_TOKEN;   }    /* Far field elements */
													
CE               { return InputCParser::CE_TOKEN;   }    /* Contact elements */
													
JE               { return InputCParser::JE_TOKEN;   }    /* Source elements */
																									 													 
PCOP             { return InputCParser::PCOP_TOKEN; }    /* Projection over coaxial port elements */
													
PRWP             { return InputCParser::PRWP_TOKEN; }    /* Projection over rw port elements */
													
PES              { return InputCParser::PES_TOKEN;  }    /* Projection over surfaces */
																										 
PEP              { return InputCParser::PEP_TOKEN;  }    /* Projection over points */
													
PEC              { return InputCParser::PEC_TOKEN;  }    /* Perfect electric conductor */
													
PBC              { return InputCParser::PBC_TOKEN;  }    /* Periodic bounday conditions */
													
PMC              { return InputCParser::PMC_TOKEN;  }    /* Perfect magnetic conductor */
																										
TEC              { return InputCParser::TEC_TOKEN;  }    /* Transversal electric condition */
													
PVIE             { return InputCParser::PVIE_TOKEN; }    /* Field integration over a volume */
													
PSIE             { return InputCParser::PSIE_TOKEN; }    /* Field integration over a surface */

CONVECTION_COEFFICIENT {
  	          yylval.double_variable = &CONVECTION_COEFFICIENT;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

FREQUENCY   {
  	          yylval.double_variable = &FREQUENCY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

PHASE       {
  	          yylval.double_variable = &PHASE;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

SINUSOIDAL_VOLUME_CURRENT {
  	          yylval.vector_double_variable = &SINUSOIDAL_VOLUME_CURRENT;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

SINUSOIDAL_SURFACE_CURRENT {
  	          yylval.vector_double_variable = &SINUSOIDAL_SURFACE_CURRENT;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

SINUSOIDAL_POINT_CURRENT {
  	          yylval.vector_double_variable = &SINUSOIDAL_POINT_CURRENT;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

COMPLEX_NEUMANN_FLOW {
  	          yylval.vector_double_variable = &COMPLEX_NEUMANN_FLOW;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

COMPLEX_IBC {
  	          yylval.vector_double_variable = &COMPLEX_IBC;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

COMPLEX_IBC_2o {
  	          yylval.vector_double_variable = &COMPLEX_IBC_2o;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

REAL_ELECTRIC_CONDUCTIVITY {
  	          yylval.double_variable = &REAL_ELECTRIC_CONDUCTIVITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_ELECTRIC_CONDUCTIVITY {
  	          yylval.double_variable = &IMAG_ELECTRIC_CONDUCTIVITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

REAL_ELECTRIC_PERMITTIVITY {
  	          yylval.double_variable = &REAL_ELECTRIC_PERMITTIVITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_ELECTRIC_PERMITTIVITY {
  	          yylval.double_variable = &IMAG_ELECTRIC_PERMITTIVITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

REAL_MAGNETIC_PERMEABILITY {
  	          yylval.double_variable = &REAL_MAGNETIC_PERMEABILITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_MAGNETIC_PERMEABILITY {
  	          yylval.double_variable = &IMAG_MAGNETIC_PERMEABILITY;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

REAL_IMPEDANCE_BOUNDARY_CONDITION {
  	          yylval.double_variable = &REAL_IMPEDANCE_BOUNDARY_CONDITION;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_IMPEDANCE_BOUNDARY_CONDITION {
  	          yylval.double_variable = &IMAG_IMPEDANCE_BOUNDARY_CONDITION;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

CONTACT_NORMALS {
  	          yylval.vector_double_variable = &CONTACT_NORMALS;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

BOUNDARY_NORMALS {
  	          yylval.vector_double_variable = &BOUNDARY_NORMALS;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

E           {
  	          yylval.vector_double_variable = &E;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_E      {
  	          yylval.vector_double_variable = &IMAG_E;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

REAL_E      {
  	          yylval.vector_double_variable = &REAL_E;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

MOD_E       {
  	          yylval.double_variable = &MOD_E;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }
			
B           {
  	          yylval.vector_double_variable = &B;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_B      {
  	          yylval.vector_double_variable = &IMAG_B;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

REAL_B      {
  	          yylval.vector_double_variable = &REAL_B;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

MOD_B       {
  	          yylval.double_variable = &MOD_B;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

H           {
  	          yylval.vector_double_variable = &H;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_H      {
  	          yylval.vector_double_variable = &IMAG_H;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

REAL_H      {
  	          yylval.vector_double_variable = &REAL_H;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

MOD_H       {
  	          yylval.double_variable = &MOD_H;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

J           {
  	          yylval.vector_double_variable = &J;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

IMAG_J      {
  	          yylval.vector_double_variable = &IMAG_J;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

REAL_J      {
  	          yylval.vector_double_variable = &REAL_J;
	          return InputCParser::VECTOR_DOUBLE_VARIABLE_TOKEN;
	        }

MOD_J       {
  	          yylval.double_variable = &MOD_J;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

TIME        {
  	          yylval.double_variable = &TIME;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

Sg          {
  	          yylval.double_variable = &Sg;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }
			
V           {
  	          yylval.double_variable = &VOLTAGE;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

NEUMANN_FLOW {
  	          yylval.double_variable = &NEUMANN_FLOW;
	          return InputCParser::DOUBLE_VARIABLE_TOKEN;
	        }

LinearSolver       { return InputCParser::LINEAR_SOLVER_TOKEN;        }
																      
SolvingStrategy    { return InputCParser::SOLVING_STRATEGY_TOKEN;     }
																      
TransientLoop      { return InputCParser::TRANSIENT_TOKEN;            }
																      
AlphaLoop          { return InputCParser::ALPHA_TOKEN;                }
																      															      
Calculate          { return InputCParser::CALCULATE_TOKEN;            }
																      
CalculateNodal     { return InputCParser::CALCULATE_NODAL_TOKEN;      }
																      
for                { return InputCParser::FOR_TOKEN;                  }
																      
Solve              { return InputCParser::SOLVE_TOKEN;                }
																      
Print              { return InputCParser::PRINT_TOKEN;                }
																      
PrintOnNodes       { return InputCParser::PRINT_ON_NODES_TOKEN;       }

Execute            { return InputCParser::EXECUTE_TOKEN;              }

CreateSolutionStep { return InputCParser::CREATE_SOLUTION_STEP_TOKEN; }

CreateTimeStep     { return InputCParser::CREATE_TIME_STEP_TOKEN;     }

[[:alnum:]_]+ {
  	          sprintf((yylval.string_value),"%s", yytext);
	          return InputCParser::WORD_TOKEN;
	        }
