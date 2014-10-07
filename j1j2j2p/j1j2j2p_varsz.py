#!/usr/bin/env python2

import logging

try:
  import argparse
  import threading
  import time
  import pickle
  import os
  import shutil
  import tempfile
  from j1j2j2p import J1J2J2P
  pass
except (ImportError) as err:
  logging.error(err)

class Foreman:
  '''
  Foreman - Worker handler.
  Calculate the ground state energy for each possible Sz of the J1 J2 J2` model specified.
  Uses DMRG for energy minimization.
  '''

  def __init__(self,nsites,nstates,exe='j1j2j2p_heis',nthreads=1,tmpdir='/tmp/'):
    '''
    Initialization method.

    Input:
      nsites - number of sites to use
      nstates - number of states to keep for each iteration of DMRG

    Optional:
      exe = 'j1j2j2p_heis' - path to J1 J2 J2` DMRG executable
      nthreads = 1 - maximum number of worker threads allowable

    Output:
      f - Foreman object
    '''
    self._nsites = nsites
    self._nstates = nstates
    self._exe = exe
    self._nthreads = nthreads
    self._tmpdir = tmpdir
    return

  def start(self):
    '''
    Start the Foreman and Workers.

    Input:
      None

    Output:
      d - dict of dicts of 2-tuples :  [j2z][j2pz] -> (sz, energy_sz)
    '''
    # Worker queue
    wlist = []

    # Consider only jz = jx = 1.0.
    jz = 1.0
    jx = jz

    # Calculate all magnetizations of the lower triangle of j2, j2p pairs.
    for i in xrange(5,6):
      for j in xrange(i+1):
        # Lower triangle of z pairs
        j2z  = i*0.2
        j2pz = j*0.2

        # Only consider the z == x set.
        j2x  = j2z
        j2px = j2pz

        # Start a worker thread, if any are available. Otherwise, wait.
        while threading.active_count() > self._nthreads: time.sleep(2)
        logging.debug('Starting thread for:\n\tj2 = '+str(j2z)+'\n\tj2p = '+str(j2pz))
        wlist.append(Worker(self._nsites,jz,jx,j2z,j2x,j2pz,j2px,self._nstates,exe=self._exe,tmpdir=self._tmpdir))
        wlist[-1].start()

    # Wait for all threads to finish, then dump results as a dict of dicts.
    while threading.active_count()-1: time.sleep(2) # May be volatile - technically, should use join and wait...
    output = {}
    for w in wlist:
      try: output[w._j2z][w._j2pz] = w.result
      except(KeyError) as err: output[w._j2z] = {w._j2pz : w.result}

    # Return the dict of dicts of 2-tuples.
    return output

class Worker(threading.Thread):
  '''
  Worker thread.
  Calculate the ground state energy for each possible Sz of the J1 J2 J2` model specified.
  Uses DMRG for energy minimization.
  '''

  def __init__(self,nsites,jz,jx,j2z,j2x,j2pz,j2px,nstates,exe='j1j2j2p_heis',tmpdir='/tmp/'):
    '''
    Input:
      nsites - number of sites to use
      jz - nearest neighbor spin coupling
      jx - nearest neighbor spin flip coupling
      j2z - odd next nearest neighbor spin coupling
      j2x - odd next nearest neighbor spin flip coupling
      j2pz - even next nearest neighbor spin coupling
      j2px - even next nearest neighbor spin flip coupling
      nstates - number of states to keep for each iteration of DMRG

    Optional:
      exe = 'j1j2j2p_heis' - path to J1 J2 J2` DMRG executable
      tmpdir = directory to use for temporary directories

    Output:
      w - Worker object.
    '''
    threading.Thread.__init__(self)
    self._nsites = nsites
    self._jz = jz
    self._jx = jx
    self._j2z = j2z
    self._j2x = j2x
    self._j2pz = j2pz
    self._j2px = j2px
    self._nstates = nstates
    self._exe = exe
    self._tmpdir = tmpdir
    self.result = []
    return

  def run(self):
    '''
    Run the Worker.
    Not to be called directly - used by Thread superclass.

    Input:
      None

    Output:
      None

    Internal: result = [(sz,str,float),...]
      sz - particular spin polarization
      str - string output from the DMRG minimization
      float - floating point value representing the minimum energy of the system for the associated sz
    '''
    # Create a temporary directory in which to execute.
    self._tmpdir = tempfile.mkdtemp(prefix=self._tmpdir)

    # Perform energy minimization for each possible value of sz (positive half).
    for sz in xrange(self._nsites/2+1):
      logging.debug('Calculating for sz = '+str(sz*2))
      sz_out,sz_energy = J1J2J2P(self._nsites,
                                 sz*2,
                                 self._jz,
                                 self._jx,
                                 self._j2z,
                                 self._j2x,
                                 self._j2pz,
                                 self._j2px,
                                 self._nstates,
                                 exe = self._exe,
                                 dir = self._tmpdir,)
      self.result.append((sz*2,sz_energy))
      fs = open('_'.join(map(str,[self._j2z,self._j2pz,sz*2]))+'.out','w')
      fs.write(sz_out)
      fs.close()

    # Remove the temporary directory tree.
    shutil.rmtree(self._tmpdir)

    # Finished.
    logging.debug('Thread complete.')
    return

def main(args):
  # Set up logging.
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.info('Starting run for arguments: ')
  logging.info(args)

  # Instantiate the Foreman and start the workers.
  foreman = Foreman(args.nsites,args.nstates,exe=args.exe,nthreads=args.nthreads,tmpdir=args.tmpdir)
  output = foreman.start()

  logging.info('Run finished. Generating pickle.')

  # Dump the dict of dict output to file.
  pickle.dump(output,open(args.out,'w'))

  logging.info('Pickle generation complete. Exiting.')

  # Done.
  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='J1 J2 J2` Energy Calculater')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default='j1j2j2p_varsz.log',
                      help='File in which to store log output. Default: j1j2j2p_varsz.log')

  parser.add_argument('nsites',type=int,help='Number of sites to use. [2+]')
  parser.add_argument('nstates',type=int,help='Number of states to keep for every truncation. [1+]')

  parser.add_argument('--exe',type=str,default='j1j2j2p_heis',help='J1 J2 J2` DMRG executable to use.')
  parser.add_argument('--nthreads',type=int,default=1,help='Number of threads to use.')
  parser.add_argument('--out',type=str,default='j1j2j2p_varsz.out',help='Output file to use.')
  parser.add_argument('--tmpdir',type=str,default='/scratch_global/',help='Directory to use for threads.')

  # Parse arguments and begin execution.
  main(parser.parse_args())
