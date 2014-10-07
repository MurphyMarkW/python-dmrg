#!/usr/bin/env python2

import logging
import pickle
import argparse

def J1J2J2PExtractEnergy(data):
  streams = {}
  for i in sorted(data):
    for j in sorted(data[i]):
      try: streams[i]
      except (KeyError): streams[i] = open(str(i)+'_j2energy.tsv','w')
      for (sz,energy) in data[i][j]:
          streams[i].write('\t'.join(map(str,[j,sz,energy]))+'\n')
          streams[j].write('\t'.join(map(str,[i,sz,energy]))+'\n')
  return

def main(args):
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.debug(args)

  data = pickle.load(open(args.file))
  J1J2J2PExtractEnergy(data)

  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='J1 J2 J2` Energy Calculater')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default='j1j2j2p_extract.log',
                      help='Logfile to use. Default: j1j2j2p_extract.log')

  parser.add_argument('file',type=str,help='Pickled file to extract energies from.')

  main(parser.parse_args())
