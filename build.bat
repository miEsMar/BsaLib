@echo off

setlocal
if not exist "build_cmake" mkdir "build_cmake"

set "cleanfirst=--clean-first"
set "cleanfirst="

rem -D CMAKE_BUILD_TYPE=Debug ^

cmake ^
   -D BUILD_SHARED_LIBS=OFF ^
   -D enable-openmp=ON ^
   -D enable-sym-ev-routine=OFF ^
   -D enable-single=OFF ^
   -D enable-gpu-code=OFF ^
   -D enable-cuda=ON ^
   -D enable-gpu-double=OFF ^
   -S . -B build_cmake %~1
if not "%errorlevel%"=="0" exit /b 1

cmake --build build_cmake --config Debug %cleanfirst%
if not "%errorlevel%"=="0" exit /b 1

cmake --build build_cmake --config Release %cleanfirst%
if not "%errorlevel%"=="0" exit /b 1
endlocal
