import os
import re
PARTTEN_RULENUMBER = re.compile(r'int numRules = \d+;')
PARTTEN_NOISETYPE = re.compile(r'int noisetype = \d+;') #
PARTTEN_FUNCTION = re.compile(r'fxn1\(i,2\) = \S+ \+ n\[i\]')
PARTTEN_INPUTRANGE = re.compile(r'double min_x = [\d\.]+, max_x = [\d\.]+;') #

MY_FUNCTION = '(sin(fxn1(i,0))*cos(fxn1(i,1)))'
NOISE_TYPE = '1'; ## 1 -- No noise, 2 -- Uniform Noise, 3 -- Gaussian, 4-- Cauchy Noise
INPUT_RANGE_START = '0.00'
INPUT_RANGE_END = '6.414'
RULE_NUMBER = '20'
file_main_bak = open('ASAM-Main.cpp', 'r')
bak_data = file_main_bak.readlines()
for l in bak_data:
    if(PARTTEN_FUNCTION.search(l)):
        print('FUNCTION: \n', l)
        m = re.compile(r'= \S+').search(l)
        l_new = l[:m.start()+2] + MY_FUNCTION + l[m.end():]
        print('Change to: \n', l_new)
    elif (PARTTEN_NOISETYPE.search(l)):
        print('NOISETYPE: \n', l)
        m = re.compile(r'\d+').search(l)
        l_new = l[:m.start()] + NOISE_TYPE + l[m.end():]
        print('Change to: \n', l_new)
    elif (PARTTEN_RULENUMBER.search(l)):
        print('RULENUMBER: \n', l)
        m = re.compile(r'\d+').search(l)
        l_new = l[:m.start()] + RULE_NUMBER + l[m.end():]
        print('Change to: \n', l_new)
    elif (PARTTEN_INPUTRANGE.search(l)):
        print('INPUTRANGE: \n', l)
        m = re.compile(r'[\d\.]+').search(l)
        l_new = l[:m.start()] + INPUT_RANGE_START + l[m.end():]
        print('Change to1: \n', l_new)
        m = re.compile(r'[\d\.]+').search(l_new)
        second_pos = m.end();
        m = re.compile(r'[\d+\.]+').search(l_new[m.end():])
        l_new = l_new[:second_pos+m.start()] + INPUT_RANGE_END + l[second_pos+m.end():]
        print('Change to2: \n', l_new)
