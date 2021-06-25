#!/usr/bin/env python3

# script to generate an install script for a dune-module,
# accounting for non-published commits and local changes

import os
import sys
import argparse
import subprocess

from util import getPersistentVersions
from util import printVersionTable
from util import getPatches

try:
    path = os.path.split(os.path.abspath(__file__))[0]
    sys.path.append(os.path.join(path, '../bin/util'))

    from getmoduleinfo import getModuleInfo, getDependencies
except Exception:
    sys.exit('Could not import getModuleInfo')

if sys.version_info[0] < 3:
    sys.exit("\nError': Python3 required")


def makeInstallScript(path,
                      fileName=None,
                      ignoreUntracked=False,
                      topFolderName='DUMUX',
                      optsFile=None,
                      skipFolders=None,
                      suppressHints=False):

    cwd = os.getcwd()
    modPath = os.path.abspath(os.path.join(cwd, path))
    modParentPath = os.path.abspath(os.path.join(modPath, '../'))
    modFolder = os.path.relpath(modPath, modParentPath)

    try:
        modName = getModuleInfo(modPath, 'Module')
    except Exception:
        sys.exit("\nError: Could not determine module name. Make sure the\n"
                 "         module path is correct and contains 'dune.module'.")

    if not fileName:
        instFileName = 'install_' + modName + '.sh'
    else:
        instFileName = fileName
    print("\n-- Creating install script '{}' for module '{}' in folder '{}'"
          .format(instFileName, modName, modPath))

    print("\n-- Determining the dependencies")
    deps = getDependencies(modPath)
    if len(deps) > 0:
        print("-> Found the following dependencies")
        print("\t| {:^50} | {:^50} |".format('module name', 'module folder'))
        print("\t" + 107*'-')
        for dep in deps:
            print('\t| {:^50} | {:^50} |'.format(dep['name'], dep['folder']))
    else:
        sys.exit("Error: Could not determine module dependencies.")

    if not skipFolders:
        depNames = [dep['name'] for dep in deps]
        depFolders = [dep['folder'] for dep in deps]
        depFolderPaths = [os.path.abspath(os.path.join(modParentPath, d)) for d in depFolders]
    else:
        depNames = [dep['name'] for dep in deps if dep['folder'] != skipFolders]
        depFolders = [dep['folder'] for dep in deps if dep['folder'] != skipFolders]
        depFolderPaths = [os.path.abspath(os.path.join(modParentPath, d)) for d in depFolders if d != skipFolders]

    print("\n-- Determining the module versions")
    try:
        versions = getPersistentVersions(depFolderPaths, ignoreUntracked)
    except Exception as err:
        print('\nCaught exception: ' + str(err))
        if 'untracked' in str(err):
            print('If you are certain that the untracked files are not needed for '
                  'the module installation, run this script with the -i flag.')
        sys.exit(1)

    if len(versions) != len(depNames):
        sys.exit("Error': Could not determine versions of all modules")

    print("-> The following (remotely available) versions are used as basis")
    print("   on top of which we will generate the required patches")
    printVersionTable(versions)

    print("\n-- Creating patches for unpublished commits and uncommitted changes")
    patches = getPatches(versions)

    # write installation shell script (switch to relative paths)
    versions = {os.path.relpath(p, modParentPath): v for p, v in versions.items()}
    patches = {os.path.relpath(p, modParentPath): v for p, v in patches.items()}

    optsRelPath = 'dumux/cmake.opts'
    if optsFile:
        optsPath = os.path.abspath(os.path.join(os.getcwd(), optsFile))
        optsRelPath = os.path.relpath(optsPath, modParentPath)

    # TODO: add support for different install script languages (e.g. Python).
    with open(instFileName, 'w') as installFile:
        installFile.write("#!/bin/bash\n\n")
        installFile.write("#"*80 + "\n")
        installFile.write("# This script installs the module '{}'"
                        " together with all dependencies.\n"
                        "\n".format(modName))

        exitFunc = 'exitWith'
        installFile.write(
            '# defines a function to exit with error message\n'
            '{} ()\n'.format(exitFunc)
            +
            '{\n'
            '    echo "\\n$1"\n'
            '    exit 1\n'
            '}\n\n'
        )

        # write patch information into install script
        if len(patches) > 0:
            for depModPath in patches.keys():
                if 'unpublished' in patches[depModPath]:
                    installFile.write("cat >> unpublished.patch <<'EOF'\n"
                                      + patches[depModPath]['unpublished']
                                      + "\nEOF\n")
                if 'uncommitted' in patches[depModPath]:
                    installFile.write("cat >> uncommitted.patch <<'EOF'\n"
                                      + patches[depModPath]['uncommitted']
                                      + "\nEOF\n")

        unitIndentation = ' '*4

        def writeCommandWithErrorCheck(command, errorMessage, indentationLevel=0):
            indent = unitIndentation*indentationLevel
            nextIndent = indent + unitIndentation
            installFile.write(
                f"{indent}if ! {command}; then\n"
                f"{nextIndent}echo \"{errorMessage}\"\n"
                "fi\n"
            )

        if topFolderName:
            installFile.write(
                '# Everything will be installed into a new folder "{}".\n'
                .format(topFolderName)
            )
            installFile.write(
                'echo "Creating the folder {} to install the modules in"\n'
                .format(topFolderName)

            )

            writeCommandWithErrorCheck(
                'mkdir -p {}'.format(topFolderName),
                '--Error: could not create top folder {}'.format(topFolderName)
            )

            writeCommandWithErrorCheck(
                'cd {}'.format(topFolderName),
                '--Error: could not enter top folder {}'.format(topFolderName))

            installFile.write('\n')
        else:
            installFile.write(
                '# Everything will be installed inside the folder from which the'
                ' script is executed.\n\n'
            )

        def writeCloneModule(depModName, depModFolder):
            installFile.write('# ' + depModName + '\n')
            installFile.write('# ' + versions[depModFolder]['branch'] + ' # '
                                + versions[depModFolder]['revision'] + ' # '
                                + versions[depModFolder]['date'] + ' # '
                                + versions[depModFolder]['author'] + '\n')

            writeCommandWithErrorCheck(
                'git clone {}'.format(versions[depModFolder]['remote']),
                '-- Error: failed to clone {}.'.format(depModName)
            )
            writeCommandWithErrorCheck(
                'cd {}'.format(depModFolder),
                '-- Error: could not enter folder {}.'.format(depModFolder)
            )
            writeCommandWithErrorCheck(
                'git checkout {}'.format(versions[depModFolder]['branch']),
                '-- Error: failed to check out branch {} in module {}.'
                .format(versions[depModFolder]['branch'], depModName)
            )
            writeCommandWithErrorCheck(
                'git reset --hard {}'.format(versions[depModFolder]['revision']),
                '-- Error: failed to check out commit {} in module {}.'
                .format(versions[depModFolder]['revision'], depModName)
            )

            installFile.write('echo "-- Successfully set up the module {} \\n"\n'
                            .format(depModName) + '\n')
            installFile.write('cd ..\n\n')

        # write the module clone first in order for the patches to be present
        if modName in depNames:
            writeCloneModule(modName, modFolder)
        for depModName, depModFolder in zip(depNames, depFolders):
            if depModName != modName:
                writeCloneModule(depModName, depModFolder)

        for patchtype in ['unpublished', 'uncommited']:
            installFile.write('if [ -f {}.patch ]; then\n'.format(patchtype))
            writeCommandWithErrorCheck(
                'git apply {}.patch'.format(patchtype),
                '--Error: failed to apply patch {}.patch'.format(patchtype),
                indentationLevel=1)
            installFile.write("else\n"
                              f"{unitIndentation} "
                              + 'echo "No patch {}.patch is found".\n'.format(patchtype))
            installFile.write('fi\n')

        # write configure command
        installFile.write('echo "-- All modules haven been cloned successfully. '
                        'Configuring project..."\n')
        writeCommandWithErrorCheck(
            './dune-common/bin/dunecontrol --opts={} all'.format(optsRelPath),
            '--Error: could not configure project'
        )

        # build the tests of the module
        installFile.write(
            '\n'
            'echo "-- Configuring successful. Compiling applications..."\n'
        )
        writeCommandWithErrorCheck(
            'cd {}/build-cmake'.format(modFolder),
            '--Error: could not enter build directory at {}/build-cmake'
            .format(modFolder)
        )
        writeCommandWithErrorCheck(
            'make build_tests',
            '--Error: applications could not be compiled. '
            'Please try to compile them manually.'
        )

    print("\n-- Successfully created install script file " + instFileName)

    subprocess.call(['chmod', 'u+x', instFileName])  # make script executable

    if len(patches) > 0:
        if not suppressHints:
            print("-> You should now commit and publish the 'patches' folder and this install script in your module such that others can use it.\n"
                  "   IMPORTANT: After you committed the patches, you have to adjust the line of the install script in which your module is checked out to a specific commit.\n"
                  "              That is, in the line 'git reset --hard COMMIT_SHA' for your module, replace COMMIT_SHA by the commit in which you added the patches.\n"
                  "              If patches had to be created for your own module, please think about comitting and pushing your local changes and rerunning this script again.")

    if not suppressHints:
        print(f"\n-- You might want to put installation instructions into the README.md file of your module, for instance:\n"
              f"     ## Installation\n"
              f"     The easiest way of installation is to use the script `{instFileName}` provided in this repository.\n"
              f"     Using `wget`, you can simply install all dependent modules by typing:\n"
              f"\n"
              f"     ```sh\n"
              f"     wget {versions[modFolder]['remote']}/{instFileName}\n"
              f"     bash {instFileName}\n"
              f"     ```\n")

    if topFolderName:
        print(f"   This will create a sub-folder `{topFolderName}`, clone all modules into it, configure the entire project and build the applications contained in this module.")
    else:
        print("   This will clone all modules into the folder from which the script is called, configure the entire project and build the applications contained in this module.")


if __name__ == '__main__':

    ###################
    # parse arguments
    parser = argparse.ArgumentParser(
        description='This script generates an install script for your dune module,'
                    'taking into account non-published commits & local changes.\n'
                    'This expects that your module is a git repository and that a '
                    'remote origin exists and has been set already.'
    )
    parser.add_argument('-p', '--path',
                        required=True,
                        help='The path to your dune module')
    parser.add_argument('-f', '--filename',
                        required=False,
                        help='File in which to write the install script')
    parser.add_argument('-i', '--ignoreuntracked',
                        required=False, action='store_true',
                        help='Use this to ignore untracked files present')
    parser.add_argument('-t', '--topfoldername',
                        required=False, default='DUMUX',
                        help='Name of the folder that the install script creates '
                            'upon execution to install the modules in. If you '
                            'pass an empty string, no folder will be created '
                            'and installation happens in place.')
    parser.add_argument('-o', '--optsfile',
                        required=False,
                        help='Provide custom opts file to be used for the call to '
                            'dunecontrol. Note that this file is required to be '
                            'contained and committed within the module or its '
                            'dependencies.')
    parser.add_argument('-s', '--skipfolders',
                        required=False, default=None,
                        help='a list of module folders to be skipped')

    parser.add_argument('-d', '--suppresshints',
                        required=False, default=False,
                        help='if output needs to be suppressed')

    cmdArgs = vars(parser.parse_args())

    makeInstallScript(
        path=cmdArgs['path'],
        fileName=cmdArgs.get('filename', None),
        ignoreUntracked=cmdArgs.get('ignoreuntracked', False),
        topFolderName=cmdArgs.get('topfoldername', None),
        optsFile=cmdArgs.get('optsFile', None),
        skipFolders=cmdArgs.get('skipfolders', None),
        suppressHints=cmdArgs.get('suppresshints', False),
    )
