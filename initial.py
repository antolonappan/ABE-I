import os
import glob as g
import datetime as dt
from numpy import arange
from configparser import ConfigParser as cnf

now = dt.datetime.now()
day = now.day
month = now.month
year = now.year
hour = now.hour
minute = now.minute

inpu = cnf()
inpu.read('abe.ini')

snapshot = inpu.get('inputs', 'snapshot')
root = inpu.get('outputs', 'root_dir')
no_of_cores = int(inpu.get('misc', 'no_of_cores'))


if not os.path.isdir(snapshot):
    raise Exception('Snapshot not found!')
    
files = g.glob(snapshot+'/snap*')

if len(files) == 0:
    raise Exception('No files in Snpashot directory!')

if len(files)%no_of_cores != 0:
    for pro in reversed(arange(1,no_of_cores+1,1)):
        if len(files)%pro == 0:
            cm_pro = pro
            break

    print("""
          I found No.of Cores you given in the abe.ini is  not
          a total divisor of No.of files in snapshot directory.
          So for running CatalogMaker ABE uses only {} cores.
          No of cores for CatalogAnalyser and DataAnalysis 
          remain unchange.
          """.format(cm_pro))
else:
    cm_pro = no_of_cores

if not os.path.isdir(root):
    os.system('mkdir -p %s' % root)

current = root + '/' +'%d-%d-%d' % (day,month,year)
if not os.path.isdir(current):
    os.system('mkdir %s' % current)
else:
    current = current + '_%d:%d' % (hour, minute)
    os.system('mkdir %s' % current)

temp_out = current+'/out'
os.system('mkdir %s' % temp_out)
os.system('mkdir %s' % temp_out+'/CatMak')
os.system('mkdir %s' % temp_out+'/CatAnl')
os.system('mkdir %s' % temp_out+'/DatAnl')

temp_out = current+'/log'
os.system('mkdir %s' % temp_out)
os.system('mkdir %s' % temp_out+'/CatMak')
os.system('mkdir %s' % temp_out+'/CatAnl')
os.system('mkdir %s' % temp_out+'/DatAnl')

inpu.set('live', 'cat_mak_out', current+'/out/'+'CatMak')
inpu.set('live', 'cat_anl_out', current+'/out/'+'CatAnl')
inpu.set('live', 'dat_anl_out', current+'/out/'+'DatAnl')

inpu.set('live', 'cat_mak_log', current+'/log/'+'CatMak')
inpu.set('live', 'cat_anl_log', current+'/log/'+'CatAnl')
inpu.set('live', 'dat_anl_log', current+'/log/'+'DatAnl')
inpu.set('live', 'no_of_cores_cm', cm_pro)

inpu.set('live', 'last_program', 'initial.py')


with open('abe.ini', 'w') as fw:
    inpu.write(fw)
