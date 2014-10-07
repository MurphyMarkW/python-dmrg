#!/usr/bin/env python2

import logging
import argparse
import sqlite3

def main(args):
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.debug(args)

  logging.info('Connecting to SQLite3 database %s.',args.db)
  db = sqlite3.connect(args.db)

  for corr in db.execute('SELECT DISTINCT corr FROM j1j2j2p_correlations'):
    for j2,j2p,sz in db.execute('SELECT DISTINCT j2,j2p,sz FROM j1j2j2p_correlations WHERE corr=?',corr):
      ofs = open('_'.join([corr[0],str(j2),str(j2p),str(sz)])+'.corr','w')
      ofs.write('\t'.join(['Site1','Site2','Value'])+'\n')
      for site1,site2,val in db.execute('SELECT DISTINCT site1,site2,val FROM j1j2j2p_correlations WHERE corr=? AND j2=? AND j2p=? AND sz=?',corr+(j2,j2p,sz)):
        ofs.write('\t'.join(map(str,[site1,site2,val]))+'\n')
      ofs.close()

  db.close()

  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='DMRG output to sqlite3 db.')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default=__file__.split('/')[-1].rsplit('.py',1)[0]+'.log',
                      help='Logfile to use. Default: '+__file__.split('/')[-1].rsplit('.py',1)[0]+'.log')

  parser.add_argument('db',type=str,help='Database file.')

  main(parser.parse_args())
