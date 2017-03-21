import os
import shutil
import re

def modify_para(cpp,bak, MY_FUNCTION, NOISE_TYPE, INPUT_RANGE_START, INPUT_RANGE_END, RULE_NUMBER):
    bak_data = file_main_bak.readlines()
    for l in bak_data:
        if(PARTTEN_FUNCTION.search(l)):
            #print('FUNCTION: \n', l)
            m = re.compile(r'= \S+').search(l)
            l_new = l[:m.start()+2] + MY_FUNCTION + l[m.end():]
            #print('Change to: \n', l_new)
        elif (PARTTEN_NOISETYPE.search(l)):
            #print('NOISETYPE: \n', l)
            m = re.compile(r'\d+').search(l)
            l_new = l[:m.start()] + NOISE_TYPE + l[m.end():]
            #print('Change to: \n', l_new)
        elif (PARTTEN_RULENUMBER.search(l)):
            #print('RULENUMBER: \n', l)
            m = re.compile(r'\d+').search(l)
            l_new = l[:m.start()] + RULE_NUMBER + l[m.end():]
            #print('Change to: \n', l_new)
        elif (PARTTEN_INPUTRANGE.search(l)):
            #print('INPUTRANGE: \n', l)
            m = re.compile(r'[\d\.]+').search(l)
            l_new = l[:m.start()] + INPUT_RANGE_START + l[m.end():]
            #print('Change to1: \n', l_new)
            m = re.compile(r'[\d\.]+').search(l_new)
            second_pos = m.end();
            m = re.compile(r'[\d+\.]+').search(l_new[m.end():])
            l_new = l_new[:second_pos+m.start()] + INPUT_RANGE_END + l[second_pos+m.end():]
            #print('Change to2: \n', l_new)
        else:
            file_main_cpp.write(l)
            continue
        file_main_cpp.write(l_new)

if not os.path.exists('ASAM-Main.bak'):
    print('Create bak')
    shutil.copyfile('ASAM-Main.cpp', 'ASAM-Main.bak')

PARTTEN_RULENUMBER = re.compile(r'int numRules = \d+;')
PARTTEN_NOISETYPE = re.compile(r'int noisetype = \d+;') #
PARTTEN_FUNCTION = re.compile(r'fxn1\(i,2\) = \S+ \+ n\[i\]')
PARTTEN_INPUTRANGE = re.compile(r'double min_x = [\d\.]+, max_x = [\d\.]+;') #

MY_FUNCTION = '(sin(fxn1(i,0))*cos(fxn1(i,1)))'
NOISE_TYPE = '1'; ## 1 -- No noise, 2 -- Uniform Noise, 3 -- Gaussian, 4-- Cauchy Noise
INPUT_RANGE_START = '0.00'
INPUT_RANGE_END = '6.414'
RULE_NUMBER = '25'

rule_number_array = ['5', '10', '15', '20', '25', '30', '50', '100', '200']
noise_type_array = ['1' , '2', '3','4']
my_function_array = ['gggg', 'fxn1(i,0)*fxn1(i,0) + fxn1(i,1)*fxn1(i,1)']

for rule_number in rule_number_array:
    for fun_index in [1]:
        # clean the file
        os.system('make clean')

        # open the source file
        file_main_bak = open('ASAM-Main.bak', 'r')
        file_main_cpp = open('ASAM-Main.cpp', 'w')

        # modify the source file
        modify_para(file_main_cpp, file_main_bak, my_function_array[fun_index], NOISE_TYPE, INPUT_RANGE_START, INPUT_RANGE_END, rule_number)

        # close the modified source file
        file_main_bak.close()
        file_main_cpp.close()

        # prepare to write result
        dirName = 'N' + NOISE_TYPE + 'F' + str(fun_index) + 'R' + rule_number
        os.system('make ASAM2D | tee log.txt')
        os.system('mkdir ' + dirName)
        os.system('mv Gauss ' + dirName + '/')     # {"Gauss", "Sinc", "Cauchy"}
        os.system('mv Sinc ' + dirName + '/')
        os.system('mv Cauchy ' + dirName + '/')
        with open('log.txt', 'a') as log:
            log.write('MY_FUNCTION = ' + my_function_array[fun_index] + '\n')
            log.write('NOISE_TYPE = ' + NOISE_TYPE + '\n')
            log.write('INPUT_RANGE_START = ' + INPUT_RANGE_START + '\n')
            log.write('INPUT_RANGE_END = ' + INPUT_RANGE_END + '\n')
            log.write('RULE_NUMBER = ' + rule_number + '\n')
        os.system('mv log.txt ' + dirName + '/')
