#!/usr/bin/env python2

import logging
import argparse
import sqlite3

import j1j2j2p_parser

def parseDMRG2DB(db,itt):
  # Verify that the table j1j2j2p_correlations exists.
  try: db.execute('CREATE TABLE j1j2j2p_correlations (nsites INTEGER, \
                                                       sz INTEGER, \
                                                       j1 FLOAT, \
                                                       j2 FLOAT, \
                                                       j2p FLOAT, \
                                                       corr TEXT, \
                                                       site1 INTEGER, \
                                                       site2 INTEGER, \
                                                       val FLOAT)')
  except (sqlite3.OperationalError) as err: pass
  finally: db.commit()

  # Verify that the table j1j2j2p_energies exists.
  try: db.execute('CREATE TABLE j1j2j2p_energies (nsites INTEGER, \
                                                  sz INTEGER, \
                                                  j1 FLOAT, \
                                                  j2 FLOAT, \
                                                  j2p FLOAT, \
                                                  energy FLOAT)')
  except (sqlite3.OperationalError) as err: pass
  finally: db.commit()

  # Parse the itterable.
  nsites,sz,j1,j2,j2p,energy,data = j1j2j2p_parser.parseDMRG(itt)

  # Insert the energy measurement into table j1j2j2p_energies.
  db.execute('INSERT INTO j1j2j2p_energies VALUES (?,?,?,?,?,?)',
             (nsites,sz,j1,j2,j2p,energy,))
  db.commit()

  # Insert the correlations into table j1j2j2p_correlations.
  for corr in sorted(data):
    print corr,len(data[corr])
    for sites in sorted(data[corr]):
      db.execute('INSERT INTO j1j2j2p_correlations VALUES (?,?,?,?,?,?,?,?,?)',
                 (nsites,sz,j1,j2,j2p,corr,sites[0],sites[1],data[corr][sites],))
  db.commit()

  # Finished.
  return


def main(args):
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.debug(args)

  logging.info('Connecting to SQLite3 database %s.',args.dbprefix+'.db')
  db = sqlite3.connect(args.dbprefix+'.db')

  for filename in args.filename:
    logging.info('Inserting %s into database.',filename);
    try: ifs = open(filename)
    except (IOError) as err:
      logging.error(err)
      continue

    parseDMRG2DB(db,ifs)
    ifs.close()

  db.close()

  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='DMRG output to sqlite3 db.')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default=__file__.split('/')[-1].rsplit('.py',1)[0]+'.log',
                      help='Logfile to use. Default: '+__file__.split('/')[-1].rsplit('.py',1)[0]+'.log')

  parser.add_argument('filename',nargs='+',type=str,help='DMRG output files to parse.')
  parser.add_argument('dbprefix',type=str,help='Prefix for database file.')

  main(parser.parse_args())
