# GROMACS H5MD Plugin for VMD 2.x (1.9.4)

University of Illinois Open Source License

Reads H5MD trajectory files produced by GROMACS 2026.0+ with the
gromacs_topology module extension.

Author: Diego E. B. Gomes - Auburn University, Alabama, AL - USA 
E-mail: dgomes@auburn.edu 
Version: 0.1

# WARNING: 
This is a temporary location, plugin will be moved to VMD 2.0 plugin source tree when validated.

# WARNING 2:
The modules gromacs.mdwn and gromacs_topology.mdwn are not official GROMACS files, and certainly incomplete specifications. They're my interpretation of the GROMCAS 2026.0 source code and the h5md I had on my hands. 


# Plugin build instructions tested for MacOSX Sequoia 15.7.3

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
