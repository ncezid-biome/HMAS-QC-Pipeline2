import re
import argparse
import time

'''
This script is to remove space between seq_id and =adapter
turn this 
>M03235:53:000000000-K394R:1:1101:16683:11154  adapter=OG0000294primerGroup8=AR_0409 
into
>M03235:53:000000000-K394R:1:1101:16683:11154=OG0000294primerGroup8=AR_0409
This is necessary because packages like vsearch demarcate the seq_id by space
'''
def parse_argument():

    parser = argparse.ArgumentParser(prog = 'remove_space.py')
    parser.add_argument('-f', '--input_file', metavar = '', required = True, help = 'Specify input file')
    return parser.parse_args()

def parse(file_to_parse):
    
    with open(file_to_parse, 'r', errors='ignore') as f:
        file = f.read()
        
    pattern = re.compile(r'\s+adapter=', re.I)
    file = re.sub(pattern, '=', file)
    
    with open(file_to_parse, 'w') as f:
        f.write(file)
    
if __name__ == "__main__":

    args = parse_argument()
    # tic = time.perf_counter()
    parse(args.input_file)
    # toc = time.perf_counter()
    # print (f"took {toc-tic:.2f} seconds")

