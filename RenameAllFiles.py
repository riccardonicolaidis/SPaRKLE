import os
import sys
import shutil as sh

def FindReplaceFiles(StartingPath, sFind, sReplace):
    for root, dirs, files in os.walk(StartingPath):
        # First rename the file
        # Then do find and replace in the files
        
        for file in files:
            # Rename the file
            if sFind in file:
                sNewFile = file.replace(sFind, sReplace)
                sOldFile = os.path.join(root, file)
                sNewFile = os.path.join(root, sNewFile)
                print(sOldFile, ' -> ', sNewFile)
                os.rename(sOldFile, sNewFile)
                
                # Find and replace in the file
                sNewFile = os.path.join(root, sNewFile)
                with open(sNewFile, 'r') as f:
                    sFile = f.read()
                sFile = sFile.replace(sFind, sReplace)
                with open(sNewFile, 'w') as f:
                    f.write(sFile)


if __name__ == '__main__':
    
    StartPath = '/home/riccardo/Documenti/GeantProjects/SPaRKLE/'
    sFind = 'LEM_'
    sReplace = 'SPaRKLE_'
    FindReplaceFiles(StartPath, sFind, sReplace)