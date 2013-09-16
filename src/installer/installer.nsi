;Include Modern UI

  !include "MUI2.nsh"

;--------------------------------
;General

  ;Name and file
  Name "TARQUIN"
  OutFile "TARQUINInstaller.exe"

  ;Default installation folder
  InstallDir "$LOCALAPPDATA\TARQUIN"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\TARQUIN" ""

 ;Request application privileges for Windows Vista
  RequestExecutionLevel user

  ;Set the font
  SetFont "Arial" "10"

  ;Reset the subcaptions where required
  SubCaption 0 ": Licence Agreement"

  ;Clear the branding text
;  BrandingText " "

  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "d:\softwares\tarquin\src\LICENCE.txt"
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH

  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
  !insertmacro MUI_UNPAGE_FINISH

  !insertmacro MUI_LANGUAGE "English"
;--------------------------------
;Installer Sections

  Section "TARQUIN"
   
    ;Set the specified install location, creating if required
    SetOutPath "$INSTDIR"
  
    ; Program Executables
    File "d:\softwares\tarquin\src\bin\tarquin.exe"
    File "d:\softwares\tarquin\src\bin\tarquinGUI.exe"
    File "d:\softwares\tarquin\src\bin\tarquin.exe.manifest"
    File "d:\softwares\tarquin\src\bin\tarquinGUI.exe.manifest"
;    File "d:\softwares\tarquin\src\bin\gnuplot.exe"

    ; misc
    File "d:\softwares\tarquin\src\bin\tarquin.bmp"
    File "d:\softwares\tarquin\src\bin\font.ttf"
    File "d:\softwares\tarquin\src\LICENCE.txt"
    File "d:\softwares\tarquin\src\README.txt"
    File "d:\softwares\tarquin\src\CREDITS.txt"

    File "d:\libs\GnuWin32\bin\freetype6.dll"
    File "d:\libs\GnuWin32\bin\zlib1.dll"

    ; cvmlib and Intel MKL dependencies
    File "d:\libs\cvmlib-5.7\cvm_ia32.dll"
    File "d:\libs\cvmlib-5.7\intel32\libguide40.dll"
    File "d:\libs\cvmlib-5.7\intel32\libifcoremd.dll"
    File "d:\libs\cvmlib-5.7\intel32\libifcoremdd.dll"
    File "d:\libs\cvmlib-5.7\intel32\Libimalloc.dll"
    File "d:\libs\cvmlib-5.7\intel32\libiomp5md.dll"
    File "d:\libs\cvmlib-5.7\intel32\libmmd.dll"
    File "d:\libs\cvmlib-5.7\intel32\libmmdd.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_cdft_core.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_def.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_intel_thread.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_core.dll"
;    File "d:\libs\cvmlib-5.7\intel32\mkl_p3.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_p4.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_p4m.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_p4p.dll"
    File "d:\libs\cvmlib-5.7\intel32\mkl_sequential.dll"
    File "d:\libs\cvmlib-5.7\intel32\msvcr71.dll"

    ; cvmlib needs VC runtime
    File "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\redist\x86\Microsoft.VC90.CRT\msvcm90.dll"
    File "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\redist\x86\Microsoft.VC90.CRT\msvcr90.dll"
    File "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\redist\x86\Microsoft.VC90.CRT\msvcp90.dll"
    File "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\redist\x86\Microsoft.VC90.CRT\Microsoft.VC90.CRT.manifest"
	
    ; fftw
    File "d:\libs\fftw\libfftw3-3.dll"
    File "d:\libs\fftw\libfftw3f-3.dll"
    File "d:\libs\fftw\libfftw3l-3.dll"

    ; example data
    SetOutPath "$INSTDIR\example_data"
    File "d:\softwares\tarquin\src\installer\bundle\w.dpt"
    File "d:\softwares\tarquin\src\installer\bundle\ws.dpt"

    SetOutPath "$INSTDIR\basis"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Ala.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Asp.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Cr.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\GABA.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Glc.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Gln.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Glu.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\GPC.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Gua.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Ins.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Lac.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Lip09.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Lip13a.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Lip13b.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Lip20.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\MM09.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\MM12.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\MM14.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\MM17.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\MM20.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\NAA.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\NAAG.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\PCh.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Scyllo.csv"
    File "d:\softwares\tarquin\src\installer\bundle\default_basis\Tau.csv"

   
    ;Store installation folder
    WriteRegStr HKCU "Software\TARQUIN" "" $INSTDIR
  
    ;Create uninstaller
    WriteUninstaller "$INSTDIR\Uninstall.exe"

    ;Create Start Menu Group
    CreateDirectory "$SMPROGRAMS\TARQUIN"

    ;Create Start Menu Shortcuts
    CreateShortCut "$SMPROGRAMS\TARQUIN\TARQUIN.lnk"   "$INSTDIR\tarquinGUI.exe"    
    CreateShortCut "$SMPROGRAMS\TARQUIN\Uninstall.lnk" "$INSTDIR\Uninstall.exe"    
    CreateShortCut "$SMPROGRAMS\TARQUIN\LICENCE.lnk"   "$INSTDIR\LICENCE.txt"    
    CreateShortCut "$SMPROGRAMS\TARQUIN\README.lnk"    "$INSTDIR\README.txt"    
    CreateShortCut "$SMPROGRAMS\TARQUIN\CREDITS.lnk"   "$INSTDIR\CREDITS.txt"    


  SectionEnd

;--------------------------------
;Uninstaller Section

  Section "Uninstall"

    ; uninstaller
    Delete "$INSTDIR\Uninstall.exe"

    ; Program Executable
    Delete "$INSTDIR\tarquin.exe"
    Delete "$INSTDIR\tarquinGUI.exe"

    ; cvmlib and Intel MKL dependencies
    Delete "$INSTDIR\cvm_ia32.dll"
    Delete "$INSTDIR\ftn_mkl_ia32.dll"
    Delete "$INSTDIR\libguide40.dll"
    Delete "$INSTDIR\libifcoremd.dll"
    Delete "$INSTDIR\libifcoremdd.dll"
    Delete "$INSTDIR\Libimalloc.dll"
    Delete "$INSTDIR\libiomp5md.dll"
    Delete "$INSTDIR\libmmd.dll"
    Delete "$INSTDIR\libmmdd.dll"
    Delete "$INSTDIR\mkl_cdft_core.dll"
    Delete "$INSTDIR\mkl_def.dll"
    Delete "$INSTDIR\mkl_intel_thread.dll"
    Delete "$INSTDIR\mkl_lapack.dll"
    Delete "$INSTDIR\mkl_p3.dll"
    Delete "$INSTDIR\mkl_p4.dll"
    Delete "$INSTDIR\mkl_p4m.dll"
    Delete "$INSTDIR\mkl_p4p.dll"
    Delete "$INSTDIR\mkl_sequential.dll"
    Delete "$INSTDIR\msvcr71.dll"

    Delete "$INSTDIR\msvcm90.dll"
    Delete "$INSTDIR\msvcr90.dll"
    Delete "$INSTDIR\msvcp90.dll"
    Delete "$INSTDIR\Microsoft.VC90.CRT.manifest"

    ; fftw
    Delete "$INSTDIR\libfftw3-3.dll"
    Delete "$INSTDIR\libfftw3f-3.dll"
    Delete "$INSTDIR\libfftw3l-3.dll"

    ; misc
    Delete "$INSTDIR\tarquin.bmp"
    Delete "$INSTDIR\font.ttf" 
    Delete "$INSTDIR\LICENCE.txt"   
    Delete "$INSTDIR\README.txt"   
    Delete "$INSTDIR\CREDITS.txt"   

    ; links
    Delete "$SMPROGRAMS\TARQUIN\TARQUIN.lnk"   
    Delete "$SMPROGRAMS\TARQUIN\Uninstall.lnk" 
    Delete "$SMPROGRAMS\TARQUIN\LICENCE.lnk"   
    Delete "$SMPROGRAMS\TARQUIN\README.lnk"    
    Delete "$SMPROGRAMS\TARQUIN\CREDITS.lnk"   

    RMDIR "$SMPROGRAMS\TARQUIN"
    RMDIR "$INSTDIR"

    DeleteRegKey /ifempty HKCU "Software\TARQUIN"

  SectionEnd
