# Build instructions tested for MacOSX Sequoia 15.7.3

#  Build Dependencies                                                                                                                                                                                          
  Essential:                                                                                            
  1. C Compiler: gcc or clang                                                                           
  2. HDF5 Library: Version 1.10.0 or newer                                                              
    - Required headers: hdf5.h                                                                          
    - Required libraries: libhdf5, libz (zlib for compression)                                          
  3. VMD Plugin Headers: From VMD installation                                                          
    - Specifically needs: molfile_plugin.h                                                              
    - Typically located at: $VMD_DIR/plugins/include/                                                   
                                                                                                        ## Platform-Specific Requirements                                                                        
  macOS:                                                                                                
  brew install hdf5                                                                                     
                                                                                                        
  Linux (Ubuntu/Debian):                                                                                
  sudo apt-get install libhdf5-dev                                                                      
                                                                                                        
  Linux (RedHat/CentOS):                                                                                
  sudo yum install hdf5-devel                                                                           
                                                                                                        
  Build Configuration                                                                                   
                                                                                                        
  The Makefile (Makefile.h5md_gromacs) handles:                                                         
  - Compiler flags: -O2 -fPIC -Wall -Wextra                                                             
  - Include paths: VMD plugin headers and HDF5 headers                                                  
  - Linker flags:                                                                                       
    - macOS: -bundle -Wl,-undefined,dynamic_lookup                                                      
    - Linux: -shared                                                                                    
  - Libraries: -lhdf5 -lz                                                                               
                                                                                                        
# Quick Start                                                                                           

```
# Export VMD_DIR to point to VMD source code (1.9.4 or 2.0a9)

export VMD_DIR=/Users/deb0054/gitlab/vmd2prototype/vmd/
                                                                                                        
# Check configuration                                                                                 
make -f Makefile.h5md_gromacs config                                                                  
                                                                                                        
# Build                                                                                               
make -f Makefile.h5md_gromacs                                                                         

# Open VMD
vmd

# Load the plugin (Terminal or TkConsole)
vmd_plugin_scandirectory REPLACE_BY_PLUGIN_FOLDER h5mdplugin_gromacs.so

# Load a molecule
mol new equil.h5md type h5md_gromacs
```
