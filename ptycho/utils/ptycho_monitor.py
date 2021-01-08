

import curses
import time
import glob
import os
import numpy as np
from datetime import datetime



def update_max_len(max_len_list, val_list):
    for it, val in enumerate(val_list):
        if max_len_list[it]< len(val):
	    max_len_list[it] = len(val)

    return max_len_list


def main(sc):
    sc.nodelay(1)
    header = 'Reconstruction monitor - cSAXS beamline'

    curses.init_pair(2, curses.COLOR_GREEN, curses.COLOR_BLACK)
    offset_left = 1
    offset_min = 5
    file_path = os.path.dirname(os.path.realpath(__file__))
    while True:
        sc.clear()
        stdscr = curses.initscr()
        ln = curses.LINES
        sc.addstr(0,curses.COLS//2-len(header)//2, header, curses.color_pair(2))
        #sc.addstr(1, 1, time.strftime("%H:%M:%S"))


        fn_list = glob.glob(file_path + '/.tmp_procID/proc_*.dat')

	data = []
	max_len_data = np.zeros([4])
	title_lst = ['Node', 'Scan number', 'Elapsed time', 'Caller']
	max_len_data = update_max_len(max_len_data, title_lst)

        fn_list.sort(key=lambda x: os.stat(os.path.join('./', x)).st_mtime)
        for ii,fn in enumerate(fn_list):
            with open(fn, 'r') as f:
                for line in f:
                    tmp_data = line.split()

            time_diff = datetime.strptime(time.strftime("%H:%M:%S"), '%H:%M:%S') - datetime.strptime(tmp_data[3], '%H:%M:%S')
	    tmp_lst = [tmp_data[0], tmp_data[1], str(time_diff), tmp_data[4]]
	    data.append(tmp_lst)
	    max_len_data = update_max_len(max_len_data, tmp_lst)

	

        sum_offset = np.zeros([5])
        sum_offset[0] = offset_left


	for ii,val in enumerate(title_lst):
	    sc.addstr(3, int(sum_offset[ii] + max_len_data[ii]//2-len(val)//2), val)
	    sum_offset[ii+1] = offset_min + sum_offset[ii] + max_len_data[ii]
	
	for ii,lst in enumerate(data):
	    for jj, lst_val in enumerate(lst):
                sc.addstr(ii+4,int(sum_offset[jj] + max_len_data[jj]//2-len(lst_val)//2), lst_val)

	
	
        sc.refresh()
        
        key = sc.getch()
        if key == ord('q'):
            break
        elif key < 0:
            time.sleep(1)

if __name__=='__main__': curses.wrapper(main)

