:: To build extensions for 64 bit Python 3.5 or later no special environment
:: needs to be configured.
::
:: To build extensions for 64 bit Python 2, we need to configure environment
:: variables to use the MSVC 2008 C++ compilers from GRMSDKX_EN_DVD.iso of:
:: MS Windows SDK for Windows 7 and .NET Framework 3.5 (SDK v7.0)
::
:: 32 bit builds do not require specific environment configurations.
::
:: Note: this script needs to be run with the /E:ON and /V:ON flags for the
:: cmd interpreter, at least for (SDK v7.0)
::
:: More details at:
:: https://github.com/cython/cython/wiki/CythonExtensionsOnWindows
:: https://stackoverflow.com/a/13751649/163740
::
:: Original Author: Olivier Grisel
:: License: CC0 1.0 Universal: https://creativecommons.org/publicdomain/zero/1.0/
:: This version based on updates for python 3.5 by Phil Elson at:
::     https://github.com/pelson/Obvious-CI/tree/master/scripts
::
:: Adapted from:
:: https://raw.githubusercontent.com/matthew-brett/multibuild/37040e31b1044468027bd86991c805006a92bf09/ci/appveyor/windows_sdk.cmd

@ECHO OFF

SET COMMAND_TO_RUN=%*
SET WIN_SDK_ROOT=C:\Program Files\Microsoft SDKs\Windows
SET WINDOWS_SDK_VERSION="v7.0"

(
    ECHO Configuring Windows SDK %WINDOWS_SDK_VERSION% for Python 2 on a 64 bit architecture
    SET DISTUTILS_USE_SDK=1
    SET MSSdk=1
    "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Setup\WindowsSdkVer.exe" -q -version:%WINDOWS_SDK_VERSION%
    "%WIN_SDK_ROOT%\%WINDOWS_SDK_VERSION%\Bin\SetEnv.cmd" /x64 /release
    ECHO Executing: %COMMAND_TO_RUN%
    call %COMMAND_TO_RUN% || EXIT 1
)
