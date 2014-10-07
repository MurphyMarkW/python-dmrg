#!/usr/bin/env python2

import logging

try:
  import argparse
  import subprocess
  pass
except (ImportError) as err:
  logging.error(err)

def J1J2J2P(nsites,sz,jz,jx,j2z,j2x,j2pz,j2px,nstates,exe='j1j2j2p_heis',dir='/tmp/'):
  '''
  Calculate the ground state energy of the J1 J2 J2` model specified.
  Uses DMRG for energy minimization.

  Input:
    nsites - number of sites to use
    sz - polarized spin of the system
    jz - nearest neighbor spin coupling
    jx - nearest neighbor spin flip coupling
    j2z - odd next nearest neighbor spin coupling
    j2x - odd next nearest neighbor spin flip coupling
    j2pz - even next nearest neighbor spin coupling
    j2px - even next nearest neighbor spin flip coupling
    nstates - number of states to keep for each iteration of DMRG

  Optional:
    exe - path to J1 J2 J2` DMRG executable
    dir - directory in which to execute

  Output: (str,float)
    str - string printout of the calculation
    float - floating point value representing the minimum energy of the system
  '''
  logging.debug('Calculating ground state energy of J1 J2 J2` chain.')

  # Generate the argument string and pass it to the executable.
  sub = subprocess.Popen(exe,stdin=subprocess.PIPE,stderr=subprocess.STDOUT,stdout=subprocess.PIPE,cwd=dir)
  out,err = sub.communicate('\n'.join(map(str,[nsites,sz,jz,jx,j2z,j2x,j2pz,j2px,nstates,0,0])))
  #return out,float(out.split('\n')[-2].split(' ')[-1]) # Doesn't always seem to be placing result in the same spot?
  return out,out[out.find('RESULTS ')+len('RESULTS '):].split('\n')[0]

def main(args):
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.info('Starting run for arguments: ')
  logging.info(args)

  # Build the specified system and calculate the ground state energy.
  output,energy = J1J2J2P(
    args.nsites,
    args.sz,
    args.jz,
    args.jx,
    args.j2z,
    args.j2x,
    args.j2pz,
    args.j2px,
    args.nstates,
  )

  logging.info('Run finished. Printing output file.')

  # Open output file and print output to file.
  try: ofs = open(args.out,'w')
  except (IOError) as err:
    logging.error(err)
    raise
  ofs.write(output)
  ofs.close()

  logging.info('Printing finished. Exiting.')

  # Finished.
  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='J1 J2 J2` Energy Calculater')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default='j1j2j2p.log',
                      help='Logfile to use. Default: j1j2j2p.log')
  parser.add_argument('--out',type=str,default='j1j2j2p.out',
                      help='File to output results. Default: j1j2j2p.out')

  parser.add_argument('nsites',type=int,help='Number of sites to use. [2+]')
  parser.add_argument('sz',type=int,help='Net spin of the system. [-nsites,nsites]')
  parser.add_argument('jz',type=float,help='Spin coupling for nearest-neighbors. [0.0,1.0]')
  parser.add_argument('jx',type=float,help='Spin flip coupling for nearest-neighbors. [0.0,1.0]')
  parser.add_argument('j2z',type=float,help='Spin coupling for odd sites\' next-nearest-neighbors. [0.0,1.0]')
  parser.add_argument('j2x',type=float,help='Spin flip coupling for odd sites\' next-nearest-neighbors. [0.0,1.0]')
  parser.add_argument('j2pz',type=float,help='Spin coupling for even sites\' next-nearest-neighbors. [0.0,1.0]')
  parser.add_argument('j2px',type=float,help='Spin flip coupling for even sites\' next-nearest-neighbors. [0.0,1.0]')
  parser.add_argument('nstates',type=int,help='Number of states to keep for every truncation. [1+]')

  main(parser.parse_args())
