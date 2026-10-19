************************************************************************************
* ERMES 20.0 - Flex++ Bison++ input_c_scanner and input_c_parser generation
************************************************************************************
* Ruben Otin
*
* United Kingdom Atomic Energy Authority (UKAEA)
*
* E-mail: ruben.otin@ukaea.uk
*
* Oxford (UK) - July 2024
************************************************************************************

To generate the files:

- "input_c_scanner.cpp" 
- "input_c_parser.cpp" 
- "input_c_parser.h"

from the template scripts: 

- "input_c_scanner.l" 
- "input_c_parser.y" 

and include them in ERMES, follow these steps:

1) Execute in Windows the batch file "flex++bison++.bat" located in this folder with a double click.

   The programs flex++ (version 2.3.8-7) and bison++ (version 1.12-8) can be downloaded from: 
   
   http://www.kohsuke.org/flex++bison++/flex++bison++.src.zip

   More info about these programs in: http://www.kohsuke.org/flex++bison++/

2) In the generated file "input_c_scanner.cpp", change the following header:

   class istream;
   #include <unistd.h>
    
   to:
   
   #include <iostream>
   using namespace std;

3) Copy&Paste "input_c_scanner.cpp" from this folder to the folder "/sources".

4) Copy&Paste "input_c_parser.cpp" from this folder to the folder "/sources".

5) Copy&Paste "input_c_parser.h" from this folder to the folder "/includes".