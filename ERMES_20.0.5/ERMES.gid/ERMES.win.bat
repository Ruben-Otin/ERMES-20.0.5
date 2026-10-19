@ECHO OFF

REM set ProblemName     = %1 
REM set ProblemPath     = %2
REM set ProblemTypePath = %3

REM OutputFile: %2\%1.info

REM GiD executes this file with its own interpreter ( not cmd.exe ), already in the problem folder.
REM Keep the commands simple: "cd /d", "2>NUL" and "2>&1" are not supported. Quotes are supported
REM ( tested ) and allow folders and problem names with spaces.
REM ERMES writes its error messages to the .info file as well.

del %2\%1.info
del %2\%1.flavia.res
del %2\%1.post.res
del %2\%1.post.msh

"%3\execs\win\ERMES_20.0.5.exe" "%1" > "%2\%1.info"
