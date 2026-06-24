@echo off

:: ============================================================
::  test.bat  —  编译 + 测试
::  在"大作业"目录下运行: test.bat [case] [topo]
::  默认: case1 + MFS2
:: ============================================================

set CASE=%1
set TOPO=%2
if "%CASE%"=="" set CASE=case1
if "%TOPO%"=="" set TOPO=MFS2

cd /d "%~dp0"

echo ========================================
echo    Multi-FPGA Topology-Driven Partitioning
echo ========================================
echo Case    : %CASE%
echo Topo    : %TOPO%
echo CurDir  : %CD%
echo ========================================

:: ---- 1. 编译 ----
echo.
echo [BUILD] Compiling...
cd code
mingw32-make clean
mingw32-make all
if errorlevel 1 (
    echo BUILD FAILED!
    cd ..
    pause
    exit /b 1
)
cd ..

echo [BUILD] Success.

:: ---- 2. 运行 ----
echo.
echo [TEST] Running...
echo ========================================

build\main.exe "Generated Benchmarks\%CASE%" "FPGA Graph\%TOPO%"

echo.
echo ========================================
echo   Test finished.
echo ========================================
pause


