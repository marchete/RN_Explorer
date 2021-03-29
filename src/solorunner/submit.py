#!/usr/bin/env python3
import sys
import time
import string
import random
import glob
import re
import traceback
import zlib
import base64
program_name=sys.argv[1]

try:
	threads=int(sys.argv[2])
except:
	threads=2
try:
	reset_seconds=int(sys.argv[3])
except:
	reset_seconds=39
try:
	lFa=int(sys.argv[4])
except:
	lFa=3200

try:
	K_A=int(sys.argv[5])
except:
	K_A=50000
try:
	K_B=int(sys.argv[6])
except:
	K_B=0
try:
	K_C=int(sys.argv[7])
except:
	K_C=2000
try:
	K_D=int(sys.argv[8])
except:
	K_D=40000
try:
	INC_TIME=int(sys.argv[9])
except:
	INC_TIME=100
try:
	INC_LFA=int(sys.argv[10])
except:
	INC_LFA=40
try:
	RECOMB=int(sys.argv[11])
except:
	RECOMB=0

program_execute = './'+program_name+' '+str(threads)+' '+str(reset_seconds)+' '+str(lFa)+' '+str(K_A)+' '+str(K_B)+' '+str(K_C)+' '+str(K_D)+' '+str(INC_TIME)+' '+str(INC_LFA)+' '+str(RECOMB)

import os
import io
import requests
import subprocess

salida = "solution.txt"

#put your email on a text file named cg_email.txt
with open('cg_email.txt', 'r') as f:
	email = f.read().strip()
#put your password on a text file named cg_pass.txt
with open('cg_pass.txt', 'r') as f:
	password = f.read().strip()
with open('downloader.py', 'r') as f:
	codeDL = f.read().strip()
number_level=0

userId=''
handle=''
session = requests.Session()

parameters=''
runningprocess=''

def killProcess():
	global runningprocess
	global program_name
	subprocess.run("killall RN_Explo* >/dev/null 2>&1", shell=True)
	if (runningprocess!=""):
		subprocess.run("killall "+runningprocess+" >/dev/null 2>&1", shell=True)
	if (program_name!=""):
		subprocess.run("killall "+program_name+" >/dev/null 2>&1", shell=True)

def recompileCode():
	global runningprocess
	global program_name
	
	if (runningprocess!=""):
		cpp_filename=runningprocess
	else:
		cpp_filename=program_name
	with open(cpp_filename+'.cpp', 'r') as f:
		original_content = f.read().strip()
	try:
		if os.path.exists('level.txt'):
			with open('level.txt', 'r') as f:
				level = f.read().strip()
			subprocess.run('cat level.txt |  tr " " "\\n" |tail -n+2| grep -v ^0$ | wc -l > cuentaNUM.txt', shell=True)
			ele=level.splitlines()
			#******************************* CODE RELOAD *************************************
			if (len(ele)>=4):
				DIM=ele[0].split()
				if (len(DIM)>=2):
					W=int(DIM[0])
					H=int(DIM[1])
					with open('cuentaNUM.txt', 'r') as f:
						MAX_NUMBERS = f.read().strip()
					if (MAX_NUMBERS!="" and MAX_NUMBERS!="0"):
						MAX_NUMBERS=int(MAX_NUMBERS)
						print(' -->Recompiling '+cpp_filename+'.cpp with W='+str(W)+' H='+str(H)+' MAX_NUMBERS='+str(MAX_NUMBERS)+" ");
						if (MAX_NUMBERS<256):
							MAX_NUMBERS=256
						new_content=re.sub(r"const\s+int\s+MAX_W\s*=.*;", "const int MAX_W = "+str(W)+";", original_content)
						new_content=re.sub(r"const\s+int\s+MAX_H\s*=.*;", "const int MAX_H = "+str(H)+";", new_content)
						new_content=re.sub(r"const\s+int\s+MAX_NUMBERS\s*=.*;", "const int MAX_NUMBERS = "+str(MAX_NUMBERS)+";", new_content)
					else:
						print("Wrong numbers:"+MAX_NUMBERS)
			subprocess.run("rm cuentaNUM.txt >/dev/null 2>&1", shell=True)
		else:
			print('Recompiling '+cpp_filename+'.cpp without changes.... ');

	except Exception as error:
		print('Error compiling: %s' % error)
	if (len(new_content)>50000):
		with open(cpp_filename+'.cpp', 'w') as f:
			f.write(new_content)
	if (new_content!=original_content or not os.path.exists(cpp_filename)):
		subprocess.run("./CLANG17.sh "+cpp_filename, shell=True)
	else:
		print('No need to recompile '+cpp_filename+'');
	
# for each level of the game
if os.path.exists('runningprocess.txt'):
	with open('runningprocess.txt', 'r') as f:
		runningprocess = f.read().strip()
while True:
	try:
		recompileCode()
		if os.path.exists('parameters.txt'):
			with open('parameters.txt', 'r') as f:
				parameters = f.read().strip()
		if os.path.exists('runningprocess.txt'):
			with open('runningprocess.txt', 'r') as f:
				runningprocess = f.read().strip()
		if ((runningprocess!='') and (parameters!='')):
			program_execute='./'+runningprocess+' '+str(threads)+' '+parameters
		if ((runningprocess=='') and (parameters!='')):
			program_execute='./'+program_name+' '+str(threads)+' '+parameters
		if ((runningprocess!='') and (parameters=='')):
			program_execute = './'+runningprocess+' '+str(threads)+' '+str(reset_seconds)+' '+str(lFa)+' '+str(K_A)+' '+str(K_B)+' '+str(K_C)+' '+str(K_D)+' '+str(INC_TIME)+' '+str(INC_LFA)+' '+str(RECOMB)
		print("Program execute:"+program_execute)
		# run the solver on level.txt and save output to solution.txt
		subprocess.run(program_execute + " > "+salida, shell=True)
		#with open('solution.txt', "w") as outfile:
		#	subprocess.run(program_execute, stdout=outfile)
		with open('level_password.txt', 'r') as f:
			level_pass = f.read().strip()
		try:
			with open(salida, 'r') as f:
				solution = f.read().strip()
		except:
			pass
		if solution == '':
			listaSoluciones = glob.glob('SOLUTION_*_'+level_pass+'.txt')
			for archivo_sol in listaSoluciones:
				print('Found solution file '+archivo_sol)
				with open(archivo_sol, 'r') as f:
					solution=''
					lines=f.readlines()
					saltaLinea=True
					for line in (lines): 
						if (not saltaLinea):
							solution=solution+line
						if (level_pass in line):
							saltaLinea=False
				break
			if solution == '':
				print('Empty solution, crashed? Retrying...')
				time.sleep(5)
				continue
		solution = level_pass + '\n' + solution
		with open('log.txt', 'a') as f:
			f.write('\nsolution:\n')
			f.write(solution)
		
		# submit the solution to CodinGame
		# login to CodinGame and get submit ID
		session = requests.Session()
		r = session.post('https://www.codingame.com/services/Codingamer/loginSiteV2', json=[email, password, True])
		userId = r.json()['codinGamer']['userId']
		r = session.post('https://www.codingame.com/services/Puzzle/generateSessionFromPuzzlePrettyId', json=[userId, "number-shifting", False])
		handle = r.json()['handle']	
		r = session.post('https://www.codingame.com/services/TestSession/play', json=[handle, {'code':solution, 'programmingLanguageId':'PHP', 'multipleLanguages':{'testIndex':1}}])
		print('replay: https://www.codingame.com/replay/' + str(r.json()['gameId']))
		next_level = ''
		if 'gameInformation' in r.json()['frames'][-2]:
			next_level = r.json()['frames'][-2]['gameInformation']
		if os.path.exists('number_level.txt'):
			with open('number_level.txt', 'r') as f:
				clls=f.read().strip()
				current_number_level=int(clls) if clls.isdigit() else 0
		if current_number_level==931:
			print('Level 931 solved. WARNING! Level 932+ needs a different level.txt downloader, writting on stderr a compressed text in zip->base64')
			break
		if not 'Code for next level' in next_level:
			next_level=''
			try:
				next_level=re.findall(r'Code for next level .level [0-9]+.: ([a-z]+)', r.text)[0]
				level_password=next_level
			except:
				pass
			if next_level == '':
				print('The solution was wrong, watch the replay for details')
				time.sleep(30)
				continue
		else:
			next_level = next_level[next_level.find(':')+2:]
			level_password = next_level.split('\n')[0]
		number_level = int(1 + r.json()['metadata']['Level'])
		if (number_level < current_number_level):
			print('Error: lower than current levelr?')
			number_level=current_number_level+1
		with open('level_password.txt', 'w') as f:
			f.write(level_password)
		with open('number_level.txt', 'w') as f:
			f.write(str(number_level))
		# get the full level
		level_input='\n'.join(next_level.split('\n')[1:])
		if (number_level > 258): #fix for CG stderr limitations. TODO: On level 932+ you need a different approach. I went to zip->base64 print to stderr.
			r = session.post('https://www.codingame.com/services/TestSession/play', json=[handle, {'code':'#Level:'+str(number_level)+'\r\necho "'+level_password+'";cat >&2', 'programmingLanguageId':'Bash', 'multipleLanguages':{'testIndex':1}}])
			level_input = r.json()['frames'][2]['stderr']
		if (level_input!='bajaNivel'):
			with open('level.txt', 'w') as f:
				f.write(level_input + '\n')
		# save input for next level
		with open('log.txt', 'a') as f:
			f.write('\nreplay: https://www.codingame.com/replay/' + str(r.json()['gameId']))
			f.write('\n\nLevel ' + str(number_level) + ':\n')
			f.write(level_password)
			f.write(level_input)
		subprocess.run("rm SAFE_*.txt >/dev/null 2>&1", shell=True)
		subprocess.run("rm APROX_*.txt >/dev/null 2>&1", shell=True)
		subprocess.run("rm EXTERN_*.txt >/dev/null 2>&1", shell=True)
		subprocess.run("rm solution.txt >/dev/null 2>&1", shell=True)
		killProcess()
	except Exception as e:
		with open('log.txt', 'a') as f:
			f.write("Exception {0}\n".format(str(e))+" "+traceback.format_exc())
		session = requests.Session()
		time.sleep(10)
