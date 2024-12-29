@echo off

setlocal
set "builddir=build_cmake"
if not exist "%builddir%" mkdir "%builddir%"

rem -D CMAKE_BUILD_TYPE=Debug ^

cmake ^
   -D BUILD_SHARED_LIBS=OFF ^
   -D enable-openmp=ON ^
   -D enable-sym-ev-routine=OFF ^
   -D enable-single=OFF ^
   -D enable-gpu-code=OFF ^
   -D enable-cuda=ON ^
   -D enable-gpu-double=OFF ^
   -S . -B %builddir%
if not "%errorlevel%"=="0" exit /b 1

cmake --build %builddir% --config Debug %*
if not "%errorlevel%"=="0" exit /b 1

cmake --build %builddir% --config Release %*
if not "%errorlevel%"=="0" exit /b 1
endlocal
