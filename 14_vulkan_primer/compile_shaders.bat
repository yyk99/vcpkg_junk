@echo off
REM Compile Vulkan shaders to SPIR-V
REM Requires glslc.exe from Vulkan SDK

echo Compiling shaders...

REM Find glslc in common Vulkan SDK locations
set GLSLC=
if exist "C:\VulkanSDK\*\Bin\glslc.exe" (
    for /f "delims=" %%i in ('dir /b /s "C:\VulkanSDK\*\Bin\glslc.exe" 2^>nul') do set GLSLC=%%i
)

if "%GLSLC%"=="" (
    echo ERROR: glslc.exe not found!
    echo Please install Vulkan SDK from https://vulkan.lunarg.com/sdk/home
    pause
    exit /b 1
)

echo Found: %GLSLC%

"%GLSLC%" shader.vert -o vert.spv
if errorlevel 1 (
    echo ERROR: Failed to compile vertex shader
    pause
    exit /b 1
)

"%GLSLC%" shader.frag -o frag.spv
if errorlevel 1 (
    echo ERROR: Failed to compile fragment shader
    pause
    exit /b 1
)

echo Shaders compiled successfully!
echo Created: vert.spv, frag.spv
pause
