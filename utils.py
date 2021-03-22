#!/usr/bin/env python3
'''
Some utilities that can be easily adapted by other scripts
'''

import sys, datetime

def formated_output(func):
    '''
    This is a decorator to make standard outputs split into chunks
    '''
    def inner(*args, **kwargs):
        print('*'*30)
        print(' '* round((30-len(func.__name__))/2) + func.__name__ + ' '* round((30-len(func.__name__))/2))
        print('*'*30)
        func(*args, **kwargs)
        print('\n')
    return inner


def print_time(string_to_print):
    now=datetime.datetime.now()
    print(f'{string_to_print}: {str(now)[:16]} \n')



if __name__ == "__main__":
    print("Nothing to do, Please call from other py files")
    sys.exit()
