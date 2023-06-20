import logging
from logging.handlers import RotatingFileHandler

OLIGO_FILE = r'/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/Salmonella/'\
                    'cgMLST_pilot/M347-21-026_fastq/M347-21-026.oligos'

LOG_FILE = "confusion_matrix.log"

# threshold for BLAST
PIDENT = 100
PCOV = 100
P_ALIGN = 0.99

# list of 'negative' controls
# CONTROL_LIST = ['2013K_0676',
#                 '2013K_1246',
#                 '2014K_0979',
#                 '2016K_0878',
#                 'Blank_IRP1',
#                 'Blank_IRP2',
#                 'Water1']

CONTROL_LIST = ['Blank_IRP1',
                'Blank_IRP2',
                '2014K_0979',
                'Water1']

# CONTROL_LIST = ['2013K_0676',
#                 '2013K_1246',
#                 '2014K_0979',
#                 '2016K_0878',
#                 'Blank_IRP1',
#                 'Blank_IRP2','2010K_0968','2011K_0052','2013K_0463','2013K_1067','2014K_0421','2015K_0074',
#                 'Water1']


log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(filename)s(%(lineno)d) - %(message)s')
log_handler = RotatingFileHandler(LOG_FILE, mode='a', maxBytes=5*1024*1024,
                                 backupCount=5, encoding=None, delay=0)
log_handler.setFormatter(log_formatter)
log_handler.setLevel(logging.INFO)
