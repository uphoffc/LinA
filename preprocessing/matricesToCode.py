#!/usr/bin/env python

import os
import numpy
import math

def readMatrixMarket(pathToMatrix):
  matrixFile = open(pathToMatrix)
  if (not matrixFile.readline().startswith('%%MatrixMarket matrix array real general')):
    print('Wrong matrix market format in file {}.'.format(pathToMatrix))
    exit(1)

  dimensions = matrixFile.readline().split()
  numberOfRows = int(dimensions[0])
  numberOfColumns = int(dimensions[1])
    
  matrix = numberOfRows * numberOfColumns * [(0, 0, '')]
  entry = 0
  for line in matrixFile:
    # format: row, column, value
    matrix[entry] = (entry % numberOfRows + 1, entry / numberOfRows + 1, line.strip());
    entry = entry + 1
  
  return { '#rows':     numberOfRows,
           '#columns':  numberOfColumns,
           'matrix':    matrix
         }

matrixFiles = os.listdir('matrices')
maxDegree = 11
with open('GlobalMatrices.h', 'w') as header:
  with open('GlobalMatrices.cpp', 'w') as cpp:
    header.write('/** This file is generated. Do not edit.\n')
    header.write('    All matrices are NxN and stored column-major.\n')
    header.write('    N = CONVERGENCE_ORDER*(CONVERGENCE_ORDER+1)/2. */\n\n')
    header.write('#ifndef GLOBALMATRICES_H_\n')
    header.write('#define GLOBALMATRICES_H_\n')
    header.write('namespace GlobalMatrices {\n')
    header.write('#if !defined(CONVERGENCE_ORDER)\n')
    header.write('#error CONVERGENCE_ORDER must be set.\n')
    cpp.write('/** This file is generated. Do not edit. */\n')
    cpp.write('#include "GlobalMatrices.h"\n')
    cpp.write('namespace GlobalMatrices {\n')
    cpp.write('#if !defined(CONVERGENCE_ORDER)\n')
    cpp.write('#error CONVERGENCE_ORDER must be set.\n')
    for degree in range(1,maxDegree+1):
      header.write('#elif CONVERGENCE_ORDER == {}\n'.format(degree+1))
      cpp.write('#elif CONVERGENCE_ORDER == {}\n'.format(degree+1))
      matrixFilesForDegree = filter(lambda x: int(x.split('_')[-1].split('.')[0]) == degree, matrixFiles)
      for matrixFile in matrixFilesForDegree:
        matrixMarket = readMatrixMarket('matrices/' + matrixFile)
        name = matrixFile.split('_')[0]
        header.write('  extern double const {}[];\n'.format(name))
        cpp.write('  double const {}[] = {{ {} }};\n'.format(name, ','.join([entry[2] for entry in matrixMarket['matrix']]) ))
    header.write('#endif\n')
    header.write('}\n')
    header.write('#endif // GLOBALMATRICES_H_\n')
    cpp.write('#endif\n')
    cpp.write('}\n')
