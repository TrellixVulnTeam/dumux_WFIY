""""
Helper functions to generate an install script for a dune-module,
accounting for non-published commits and local changes
"""

import os
import sys
import textwrap

from util.common import getPersistentVersions, versionTable, getPatches
from util.moduleinfo import getModuleInfo
from util.installscript_writer import InstallScriptWriterBash
from util.installscript_writer import InstallScriptWriterPython

if sys.version_info[0] < 3:
    sys.exit("\nError': Python3 required")


def supportedLanguages():
    return ['python', 'bash']


def getScriptExtension(language):
    assert language in supportedLanguages()
    ext = {
        'python': '.py',
        'bash': '.sh'
    }
    return ext[language]


def makeScriptWriter(language):
    if language == 'bash':
        return InstallScriptWriterBash()
    elif language == 'python':
        return InstallScriptWriterPython()
    raise ValueError(f'Could not create writer for language {language}')


def getDefaultScriptName(modName, language):
    return 'install_{}{}'.format(
        modName,
        getScriptExtension(language)
    )


def printProgressInfo(infoLines, indLevel=0):
    firstPrefix = '\n' + '--'*(indLevel+1)
    emptyPrefix = firstPrefix.replace('-', ' ').strip('\n')
    print(f"{firstPrefix} {infoLines[0]}")
    for line in infoLines[1:]:
        print(f"{emptyPrefix} {line}")


def filterDependencies(dependencies, skipFolders=[]):
    if not skipFolders:
        return dependencies
    else:
        def skipFolder(folderName):
            return any(folderName == os.path.basename(path) for path in skipFolders)
        return [
            dep for dep in dependencies if not skipFolder(dep['folder'])
        ]


def addDependencyVersions(dependencies, ignoreUntracked=False):
    def getKey(dependency):
        return dependency['path']

    versions = getPersistentVersions(
        [getKey(d) for d in dependencies], ignoreUntracked
    )
    if len(versions) != len(dependencies):
        raise Exception("Not all versions of all modules could be found.")

    mergedResult = []
    for depInfo in dependencies:
        versionInfo = versions[getKey(depInfo)]
        mergedResult.append({**depInfo, **versionInfo})
    return mergedResult


def addDependencyPatches(dependenciesWithVersions):
    def getKey(dependency):
        return dependency['path']

    patches = getPatches({
        getKey(d): d for d in dependenciesWithVersions
    })

    mergedResult = []
    for depInfo in dependenciesWithVersions:
        patch = patches[getKey(depInfo)]
        mergedResult.append({**depInfo, **patch})
    return mergedResult


def makeInstallScript(modPath,
                      dependencies,
                      scriptName,
                      writer,
                      topFolderName='DUMUX',
                      optsFile=None):

    modPath = os.path.abspath(modPath)
    modName = getModuleInfo(modPath, 'Module')

    modOptsFile = '{}/cmake.opts'.format(modPath)
    if not optsFile:
        if os.path.isfile(modOptsFile):
            optsFile = '{}/cmake.opts'.format(os.path.relpath(modPath))
        else:
            optsFile = 'dumux/cmake.opts'
    if os.path.isabs(optsFile):
        raise ValueError("Opts file must be given as relative path")
    if not any(optsFile.startswith(d['folder']) for d in dependencies):
        print("Warning: opts file is not contained in any of the dependencies")

    with open(scriptName, 'w') as script:

        writer.setOutputStream(script)
        writer.writeSheBang()

        script.write('\n')
        writer.writeComment(textwrap.dedent(f"""\

            This installs the module {modName} and its dependencies.
            The exact revisions used are listed in the table below.
            However, note that this script may also apply further patches.
            If so, all patches are required to be the current folder, or,
            in the one that you specified as argument to this script.

        """))

        script.write('\n')
        writer.writeComment(versionTable(dependencies))

        script.write('\n')
        writer.writePreamble(topFolderName)

        for dep in dependencies:
            script.write('\n')
            writer.writeMessageOutput('Installing {}'.format(dep['name']))
            writer.writeInstallation(dep)

        for dep in dependencies:
            def writePatch(patch, moduleName, description):
                script.write('\n')
                writer.writeMessageOutput(
                    f'Applying patch for {description} in {moduleName}'
                )
                writer.writePatchApplication(dep['folder'], patch)

            if dep['untracked'] is not None:
                description = 'untracked files'
                writePatch(dep['untracked'], description, dep['name'])
            if dep['unpublished'] is not None:
                description = 'unpublished commits'
                writePatch(dep['unpublished'], description, dep['name'])
            if dep['uncommitted'] is not None:
                description = 'uncommitted changes'
                writePatch(dep['uncommitted'], description, dep['name'])

        script.write('\n')
        writer.writeMessageOutput('Configuring project')
        writer.writeConfiguration(optsFile)


def printFoundDependencies(deps):
    if len(deps) > 0:
        infoText = ["Found the following dependencies"]
        infoText.extend(
            versionTable(
                deps, {'name': 'module name', 'path': 'folder'}
            ).split('\n')
        )
        printProgressInfo(infoText)


def printFoundVersionInfo(dependenciesWithVersions):
    table = versionTable(dependenciesWithVersions)
    printProgressInfo(
        ["The following (remotely available) versions are used as a basis",
         "on top of which the required patches will be automatically created:",
         "\n{}".format(table)]
    )


def printFinalMessage(scriptName,
                      topFolderName=None):

    if topFolderName:
        description = textwrap.dedent(f"""\
            Running this script will create a folder `{topFolderName}`, clone all modules
            into it, configure the entire project and build the contained applications
        """)
    else:
        description = textwrap.dedent(f"""\
            Running this script will clone all modules into the folder from which it is
            called, configure the entire project and build the contained applications
        """)

    printProgressInfo(['Info:', description])
