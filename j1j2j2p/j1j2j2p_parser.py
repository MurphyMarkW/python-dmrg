#!/usr/bin/env python2

import logging
import argparse
import re


def isCorrSZ(line):
  try: return bool(isCorrSZ.patt.match(line.strip('\n')))
  except: isCorrSZ.patt = re.compile('^Sz\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrSZ.patt.match(line.strip('\n')))

def parseCorrSZ(line):
  try: parseCorrSZ.patt
  except: parseCorrSZ.patt = re.compile('^Sz\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrSZ.patt.findall(line): return []
  groups = parseCorrSZ.patt.search(line.strip('\n'))
  id = 'Sz'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrSP(line):
  try: return bool(isCorrSP.patt.match(line.strip('\n')))
  except: isCorrSP.patt = re.compile('^S\+\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrSP.patt.match(line.strip('\n')))

def parseCorrSP(line):
  try: parseCorrSP.patt
  except: parseCorrSP.patt = re.compile('^S\+\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrSP.patt.findall(line): return []
  groups = parseCorrSP.patt.search(line.strip('\n'))
  id = 'S+'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrSM(line):
  try: return bool(isCorrSM.patt.match(line.strip('\n')))
  except: isCorrSM.patt = re.compile('^S\-\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrSM.patt.match(line.strip('\n')))

def parseCorrSM(line):
  try: parseCorrSM.patt
  except: parseCorrSM.patt = re.compile('^S\-\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrSM.patt.findall(line): return []
  groups = parseCorrSM.patt.search(line.strip('\n'))
  id = 'S+'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrHDIMER(line):
  try: return bool(isCorrHDIMER.patt.match(line.strip('\n')))
  except: isCorrHDIMER.patt = re.compile('^Hdimer\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrHDIMER.patt.match(line.strip('\n')))

def parseCorrHDIMER(line):
  try: parseCorrHDIMER.patt
  except: parseCorrHDIMER.patt = re.compile('^Hdimer\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrHDIMER.patt.findall(line): return []
  groups = parseCorrHDIMER.patt.search(line.strip('\n'))
  id = 'HDimer'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrSZSZ(line):
  try: return bool(isCorrSZSZ.patt.match(line.strip('\n')))
  except: isCorrSZSZ.patt = re.compile('^\(\d+,\d+\)\*Sz\(\d+,\d+\)Sz\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrSZSZ.patt.match(line.strip('\n')))

def parseCorrSZSZ(line):
  try: parseCorrSZSZ.patt
  except: parseCorrSZSZ.patt = re.compile('^\(\d+,\d+\)\*Sz\((\d+),\d+\)Sz\((\d+),\d+\) = \((.+),.+\)')
  if not parseCorrSZSZ.patt.findall(line): return []
  groups = parseCorrSZSZ.patt.search(line.strip('\n'))
  id = 'SzSz'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrSPSM(line):
  try: return bool(isCorrSPSM.patt.match(line.strip('\n')))
  except: isCorrSPSM.patt = re.compile('^\(\d+,\d+\)\*S\+\(\d+,\d+\)S\-\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrSPSM.patt.match(line.strip('\n')))

def parseCorrSPSM(line):
  try: parseCorrSPSM.patt
  except: parseCorrSPSM.patt = re.compile('^\(\d+,\d+\)\*S\+\((\d+),\d+\)S\-\((\d+),\d+\) = \((.+),.+\)')
  if not parseCorrSPSM.patt.findall(line): return []
  groups = parseCorrSPSM.patt.search(line.strip('\n'))
  id = 'S+S-'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrHSCALAR(line):
  try: return bool(isCorrHSCALAR.patt.match(line.strip('\n')))
  except: isCorrHSCALAR.patt = re.compile('^Hscalar_scalar\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrHSCALAR.patt.match(line.strip('\n')))

def parseCorrHSCALAR(line):
  try: parseCorrHSCALAR.patt
  except: parseCorrHSCALAR.patt = re.compile('^Hscalar_scalar\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrHSCALAR.patt.findall(line): return []
  groups = parseCorrHSCALAR.patt.search(line.strip('\n'))
  id = 'HScalar'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrHVECTOR(line):
  try: return bool(isCorrHVECTOR.patt.match(line.strip('\n')))
  except: isCorrHVECTOR.patt = re.compile('^Hvector_vector\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrHVECTOR.patt.match(line.strip('\n')))

def parseCorrHVECTOR(line):
  try: parseCorrHVECTOR.patt
  except: parseCorrHVECTOR.patt = re.compile('^Hvector_vector\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrHVECTOR.patt.findall(line): return []
  groups = parseCorrHVECTOR.patt.search(line.strip('\n'))
  id = 'HVector'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isCorrHPOLAR(line):
  try: return bool(isCorrHPOLAR.patt.match(line.strip('\n')))
  except: isCorrHPOLAR.patt = re.compile('^Hpolar\(\d+,\d+\) = \(.+,.+\)')
  finally: return bool(isCorrHPOLAR.patt.match(line.strip('\n')))

def parseCorrHPOLAR(line):
  try: parseCorrHPOLAR.patt
  except: parseCorrHPOLAR.patt = re.compile('^Hpolar\((\d+),(\d+)\) = \((.+),.+\)')
  if not parseCorrHPOLAR.patt.findall(line): return []
  groups = parseCorrHPOLAR.patt.search(line.strip('\n'))
  id = 'HPolar'
  site1 = int(groups.group(1))
  site2 = int(groups.group(2))
  val = groups.group(3)
  return [id,site1,site2,val]

def isENERGY(line):
  try: return bool(isENERGY.patt.match(line.strip('\n')))
  except: isENERGY.patt = re.compile('^RESULTS .+')
  finally: return bool(isENERGY.patt.match(line.strip('\n')))

def parseENERGY(line):
  try: parseENERGY.patt
  except: parseENERGY.patt = re.compile('^RESULTS (.+)')
  if not parseENERGY.patt.findall(line): return []
  groups = parseENERGY.patt.search(line.strip('\n'))
  id = 'Energy'
  val = float(groups.group(1))
  return [id,val]


_DMRG_PARSERS_PAIRS = [
  (isCorrSZ,parseCorrSZ),
  (isCorrSP,parseCorrSP),
  (isCorrSM,parseCorrSM),
  (isCorrHDIMER,parseCorrHDIMER),
  (isCorrSZSZ,parseCorrSZSZ),
  (isCorrSPSM,parseCorrSPSM),
  (isCorrHSCALAR,parseCorrHSCALAR),
  (isCorrHVECTOR,parseCorrHVECTOR),
  (isCorrHPOLAR,parseCorrHPOLAR),
]

def parseDMRG(itt):
  # Variables to populate and return.
  sites   = 0
  sz      = 0
  j1      = 0
  j2      = 0
  j2p     = 0
  energy  = 0
  data    = {}

  # Grab the header lines for the particular run.
  last = ''
  for line in itt:
    line = line.strip('\n')
    if 'Number of sites:' in last: sites = int(line)
    if 'Sz:' in last: sz = int(line)
    if 'Jz =' in last: j1 = float(line)
    if 'J2z =' in last: j2 = float(line)
    if 'J2pz =' in last: j2p = float(line)
    if 'J2pz =' in last: break # Nothing more of interest to grab. Change as needed.
    last = line

  # Parse the remaining lines.
  for line in itt:
    # Check if the line is an energy line.
    if isENERGY(line): # We found the energy - last item to be parsed from the itterable.
      id,energy = parseENERGY(line)
      break # Stop parsing.

    # Try each parser.
    for check,parser in _DMRG_PARSERS_PAIRS:
      if check(line):
        id,site1,site2,val = parser(line)
        try: data[id][(site1,site2)] = val
        except (KeyError): data[id] = {(site1,site2) : val}

  # Finished. Return the data dict.
  return sites,sz,j1,j2,j2p,energy,data


def main(args):
  logging.basicConfig(level=args.loglevel,
                      filename=args.logfile,
                      format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',)
  logging.debug(args)

  logging.info('Opening file %s for read only.',args.filename)
  try: ifs = open(args.filename)
  except (IOError) as err:
    logging.error(err)
    raise

  logging.info('Parsing file %s as DMRG output.',args.filename)
  sites,sz,j1,j2,j2p,energy,data = parseDMRG(ifs)

  print
  print 'Sites: ',sites
  print 'Sz:    ',sz
  print 'J1:    ',j1
  print 'J2:    ',j2
  print 'J2P:   ',j2p
  print 'Energy:',energy
  print

  for corr in sorted(data):
    print 'Correlation:',corr
    print '\t'.join(['','Sites','Value'])
    for sites in sorted(data[corr]):
      print '\t'.join(map(str,['',sites,data[corr][sites]]))
    print

  return

if __name__ == '__main__':
  # Set up command line arguments.
  parser = argparse.ArgumentParser(description='DMRG output to sqlite3 db.')
  parser.add_argument('--debug',dest='loglevel',action='store_const',const=logging.DEBUG,default=logging.INFO,
                      help='Set logging level to debug.')
  parser.add_argument('--logfile',type=str,default=__file__.split('/')[-1].rsplit('.py',1)[0]+'.log',
                      help='Logfile to use. Default: '+__file__.split('/')[-1].rsplit('.py',1)[0]+'.log')

  parser.add_argument('filename',type=str,help='DMRG output files to parse.')

  main(parser.parse_args())
