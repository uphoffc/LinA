#!/usr/bin/env python

import os
import numpy
import math

matrixFiles = os.listdir('matrices')
maxDegree = 12
with open('basisfunctions.h', 'w') as out:
  with open('basisfunctions.c', 'r') as bf:
    out.write('/** This file is generated. Do not edit. */\n')
    out.write('#ifndef BASISFUNCTIONS_H_\n')
    out.write('#define BASISFUNCTIONS_H_\n')
    out.write('#include <cmath>\n')
    nbf = 0;
    names = list()
    for line in bf:
      name = 'basisFunction' + str(nbf)
      out.write('static double {}(double xi, double eta) {{\n'.format(name))
      out.write('  double ' + line)
      out.write('  return phi;\n')
      out.write('}\n')
      nbf = nbf + 1
      names.append(name)
    out.write('static double (* const basisFunctions[])(double, double) = {{ {} }};\n'.format(','.join(names)))
    out.write('#endif // BASISFUNCTIONS_H_\n')
