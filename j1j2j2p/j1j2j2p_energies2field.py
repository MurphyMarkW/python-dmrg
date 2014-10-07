#!/usr/bin/env python2

import logging
import argparse

def EnergySpinFieldIntersections(pairs):
  zeroes = {}

  # Find the field strengths that cause intersections.
  for i,(energy1,spin1) in enumerate(pairs[:-1]):
    for energy2,spin2 in pairs[i+1:]:
      if energy1-energy2:
        try: field = (energy1 - energy2)/(spin2-spin1)
        except (ZeroDivisionError):
          raise ZeroDivisionError('encountered differing energies for duplicate spins')
      else: continue
      try: zeroes[field].add((energy1,spin1))
      except(KeyError): zeroes[field] = set([(energy1,spin1)])
      zeroes[field].add((energy2,spin2))

  return zeroes

def EnergySpinFieldMinima(pairs):
  # Find the discontinuities in the energy spin pairs for different field strengths.
  intersections = EnergySpinFieldIntersections(pairs)
  
  # Check each intersection to see if the members are at a local minima.
  for i,field in enumerate(intersections.keys()):
    if not (min(pairs,key=lambda x:x[0]+x[1]*field) in intersections[field]):
      del intersections[field]

  return intersections

def main(args):
  # Set up logging.
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.debug(args)

  fs = open(args.file)
  data = {}
  for l in fs:
    j2,j2p,sz,e0 = map(float,l.strip('\n').split('\t'))
    try: data[j2].append((e0,sz))
    except (KeyError): data[j2] = [(e0,sz)]

  ofs = open('.'.join(map(str,[args.prefix,'dat'])),'w')

  for j2 in sorted(data):
    disc = EnergySpinFieldMinima(data[j2])
    for field in sorted(disc):
      for i in disc[field]:
        ofs.write('\t'.join(map(str,[j2,field,i[1]]))+'\n')

  ofs.close()

  # Done.
  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='Example template multithreaded Python application.')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default='.log',
                      help='File in which to store log output. Default: EnergySpinFieldMinima.log')
  parser.add_argument('--nthreads',type=int,default=1,help='Number of threads to use.')

  parser.add_argument('file',type=str,help='Input file.')
  parser.add_argument('prefix',type=str,help='Output file prefix.')

  # Parse arguments and begin execution.
  main(parser.parse_args())
