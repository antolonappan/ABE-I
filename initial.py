import os
import glob as g
import datetime as dt
import multiprocessing as m
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
bin_size = float(inpu.get('misc', 'bin'))
box = float(inpu.get('misc', 'box'))
mass_cutoff = inpu.get('misc', 'mass_cutoff')
cutoff_value = float(inpu.get('misc', 'cutoff_value'))

if not os.path.isdir(snapshot):
    raise Exception('Snapshot not found!')
    
files = g.glob(snapshot+'/snap*')

if len(files) == 0:
    raise Exception('No files in Snpashot directory!')

if not len(files)%m.cpu_count() == 0:
    print("""
          I found No.of Cores available in this system is {}.
          If you are using full cores, Cataloge_Maker cannot 
          be parallelized properly. Please find No of Cores
          which can be a proper divisor of {}(No of snapshot files)
          Note: This prob;em is only for cataloge_Maker. you
          use full cores for Cataloge_Analyser and DataAnalyser
          """.format(m.cpu_count(),len(files)))
    

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

inpu.set('live', 'last_program', 'initial.py')


with open('abe.ini', 'w') as fw:
    inpu.write(fw)

writer = open('%s/summary.txt' % current, 'w')
writer.write("Time at which execution started was %d hrs: %d mins. " % (hour, minute))
writer.write("Snapshot used for analysis is '%s'. " % snapshot)
writer.write("%d Kpc was the box size for selecting gas particles. " % box)  
if mass_cutoff == 'T':
   writer.write("Only BHs with mass greater than %.e solar mass are considered for analysis." % cutoff_value)
writer.write(" Bin size chosed was %.2f" % bin_size)
writer.close() 
