#!/usr/bin/env python3

"""
Create files CMakeLists.txt for the given folder and all subfolders,
including the add_subdirectory(...) and install(...) commands.
Defaults to the folder `dumux` that contains the header files,
if no folder was specified.
"""

import os
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=str, nargs='?', help='the folder to create CMakeLists.txt\'s for', default=None)
    args = vars(parser.parse_args())

    # default to the dumux folder (relative path to the location of this script)
    if args['folder'] is None:
        rootDir = os.path.dirname(os.path.abspath(__file__)) + "/../../dumux"
    else:
        rootDir = args['folder']


    ignore_folders = ["", "io/format/fmt", "io/xml"]
    extensions = [".hh", ".inc"]
    for fullFolderName, subFolders, files in os.walk(rootDir):
        # alphabetically sort
        subFolders = sorted(subFolders)
        files = sorted(files)
        # get folder name relative to dumux
        folderName = fullFolderName.replace(rootDir + '/', '').replace(rootDir, '')
        if folderName not in ignore_folders:
            with open(fullFolderName + "/CMakeLists.txt", "w") as cmakelists:
                # add subfolders
                for subFolder in subFolders:
                    cmakelists.write("add_subdirectory({})\n".format(subFolder))

                headersExist = False
                for fileName in files:
                    ext = os.path.splitext(fileName)[1]
                    if ext in extensions:
                        headersExist = True
                        break

                if headersExist:
                    if subFolders: cmakelists.write("\n")
                    # collect all files to be installed in a CMake variable
                    headers_variable = "DUMUX_" + folderName.upper().replace("/", "_") + "_HEADERS"
                    cmakelists.write("file(GLOB {}{})\n".format(headers_variable, " *".join([''] + extensions)))
                    cmakelists.write("install(FILES ${{{}}}\n".format(headers_variable))
                    cmakelists.write("        DESTINATION ${{CMAKE_INSTALL_INCLUDEDIR}}/dumux/{})\n".format(folderName))
