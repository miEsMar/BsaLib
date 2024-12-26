@echo off

setlocal enabledelayedexpansion
set "SCRIPT_DIR=%~d0"
set "basedir_=%~1"
if "!basedir_!"=="" ( set "basedir_=!SCRIPT_DIR!\..\..\" )
set "outdir_=!basedir_!\fortran\build"
set "indir_=!basedir_!\BSACL\build"
set "libname_=BSACLlib_Fortran_cl.lib"
for /f "tokens=* delims=" %%a in (' dir /b "!outdir_!" ') do (
   set "starto_=%%a"
   set "starto_=!starto_:~0,5!"
   if "!starto_!"=="Release" (
      for /f "tokens=* delims=" %%d in (' dir /b "!indir_!" ') do (
         set "starti_=%%d"
         set "starti_=!starti_:~0,5!"
         if "!starti_!"=="!starto_!" (
            xcopy /Y "!indir_!\%%d\lib\!libname_!" "!outdir_!\%%a\lib\"
            xcopy /Y "!indir_!\%%d\mod\bsacl.mod" "!outdir_!\%%a\mods\"
         )
      )
   )
   if "!starto_!"=="Debug" (
      for /f "tokens=* delims=" %%d in (' dir /b "!indir_!" ') do (
         set "starti_=%%d"
         set "starti_=!starti_:~0,5!"
         if "!starti_!"=="!starto_!" (
            xcopy /Y "!indir_!\%%d\lib\!libname_!" "!outdir_!\%%a\lib\"
            xcopy /Y "!indir_!\%%d\mod\bsacl.mod" "!outdir_!\%%a\mods\"
         )
      )
   )
)
endlocal


@REM setlocal enabledelayedexpansion
@REM set "basedir_=%~1"
@REM set "outdir_=!basedir_!\fortran\build"
@REM set "indir_=!basedir_!\ocl\BSACL\build"
@REM for /f "tokens=* delims=" %%a in (' dir /b "!outdir_!" ') do (
@REM    set "starto_=%%a"
@REM    set "starto_=!starto_:~0,5!"
@REM    if "!starto_!"=="Release" goto :ok_
@REM    if "!starto_!"=="Debug" goto :ok_
@REM    goto :loop_
@REM    :ok_
@REM    for /f "tokens=* delims=" %%d in (' dir /b "!indir_!" ') do (
@REM       set "starti_=%%d"
@REM       set "starti_=!starti_:~0,5!"
@REM       if "!starti_!"=="!starto_!" (
@REM          xcopy /Y "!indir_!\%%d\lib\BSACLlib_Fortran.lib" "!outdir_!\%%a\lib\"
@REM          @REM xcopy /Y "!indir_!\%%d\mod\bsacl.mod" "!outdir_!\%%a\mods\"
@REM       )
@REM    )
@REM    :loop_
@REM )
@REM endlocal